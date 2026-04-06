# --- 1. BIBLIOTHÈQUE STANDARD PYTHON ---
import os
import io
import re
import logging
import warnings
import json
import traceback
import webbrowser
from datetime import datetime, timedelta
from pathlib import Path
import pickle

# --- 2. BIBLIOTHÈQUES TIERCES GÉNÉRALES ---
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.spatial import KDTree

# --- 3. INTERFACE GRAPHIQUE (TKINTER) ---
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, Toplevel, simpledialog
from tkinter import scrolledtext

# --- 4. VISUALISATION (MATPLOTLIB) ---
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# --- 5. ASTRONOMIE (ASTROPY) ---
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points, proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, match_coordinates_sky, EarthLocation, AltAz
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.utils.exceptions import AstropyWarning

# --- 6. PHOTOMÉTRIE (PHOTUTILS) ---
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from photutils.centroids import centroid_2dg, centroid_quadratic, centroid_sources
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats

# --- 7. REQUÊTES EXTERNES (ASTROQUERY) ---
from astroquery.vizier import Vizier
from astroquery.jplhorizons import Horizons
# --- 8. MODULES LOCAUX ---
import config
# --- NOUVEAU : Import du module d'astrométrie optimisé ---
from core.fast_astrometry import FastAstrometrySolver, GaiaCache, smart_aperture_selection, optimized_query_gaia, calculate_astrometric_statistics, open_fits_with_fixed_wcs
# Tentative d'import CuPy pour GPU
try:
    import cupy as cp
    HAS_GPU = True
except ImportError:
    HAS_GPU = False
# Détection KBMOD : exécutée via WSL (pas d'import kbmod sous Windows)
from utils.wsl_utils import windows_path_to_wsl
from core.periodogram_tools import run_lomb_scargle
from core.asteroid_lightcurve_model import light_curve_model, fit_light_curve
from core.asteroid_shape_model import load_shape
from utils.logging_handler import TextHandler as UILogTextHandler
# --- CONFIGURATION ---
warnings.filterwarnings("ignore", category=AstropyWarning)
logger = logging.getLogger("AsteroidPhotometryTab")
Vizier.server = 'http://vizier.cfa.harvard.edu/'

# =============================================================================
# FONCTIONS UTILITAIRES (MATHS & ASTROMÉTRIE)
# =============================================================================
def _compute_fov_center_and_radius_for_astrometry(wcs, shape):
    ny, nx = shape
    px = np.array([0, nx, 0, nx])
    py = np.array([0, 0, ny, ny])
    corners_sky = wcs.pixel_to_world(px, py)
    center_ra = np.mean(corners_sky.ra.deg)
    center_dec = np.mean(corners_sky.dec.deg)
    center = SkyCoord(center_ra * u.deg, center_dec * u.deg, frame="icrs")
    radius = center.separation(corners_sky).max()
    return center, radius * 1.1

def gauss(x, *p):
    A, mu, sigma, C = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2)) + C

def refine_centroid(data, x_init, y_init, box_size=11):
    h, w = data.shape
    x_init = float(x_init) if isinstance(x_init, (np.ndarray, np.generic)) else float(x_init)
    y_init = float(y_init) if isinstance(y_init, (np.ndarray, np.generic)) else float(y_init)
    x_int, y_int = int(round(x_init)), int(round(y_init))
    r = box_size // 2
    if x_int < r or x_int >= w - r or y_int < r or y_int >= h - r:
        return x_init, y_init
    cutout = data[y_int - r : y_int + r + 1, x_int - r : x_int + r + 1]
    try:
        xc, yc = centroid_2dg(cutout)
        xc = float(xc) if isinstance(xc, (np.ndarray, np.generic)) else float(xc)
        yc = float(yc) if isinstance(yc, (np.ndarray, np.generic)) else float(yc)
        return x_int - r + xc, y_int - r + yc
    except Exception as e:
        logger.debug(f"Échec centroïde : {e}")
        return x_init, y_init

def estimate_fwhm_marginal(data, x, y, box_size=25):
    """
    Estime le FWHM d'une source en pixels (version complète comme pour les exoplanètes).
    Retourne (fwhm_val, (r, row, pr), (c, col, pc)) ou (None, None, None) si l'estimation échoue.
    """
    try:
        x, y = int(x), int(y)
        half = box_size // 2
        if y-half < 0 or y+half+1 > data.shape[0] or x-half < 0 or x+half+1 > data.shape[1]: 
            return None, None, None
        
        sub = data[y-half:y+half+1, x-half:x+half+1]
        if sub.size == 0 or not np.isfinite(sub).any():
            return None, None, None
        
        row = np.sum(sub, axis=1)
        col = np.sum(sub, axis=0)
        
        # Vérifier qu'il y a un signal significatif
        if row.max() <= row.min() or col.max() <= col.min():
            return None, None, None
        
        r = np.arange(len(row))
        c = np.arange(len(col))
        gauss_func = lambda x, a, x0, s, off: a * np.exp(-(x-x0)**2/(2*s**2)) + off
        
        # Fit gaussien avec contraintes raisonnables
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pr, _ = curve_fit(gauss_func, r, row, p0=[row.max(), half, 2.0, np.median(row)], 
                                  maxfev=1000, bounds=([0, 0, 0.5, -np.inf], [np.inf, box_size, 10.0, np.inf]))
                pc, _ = curve_fit(gauss_func, c, col, p0=[col.max(), half, 2.0, np.median(col)], 
                                  maxfev=1000, bounds=([0, 0, 0.5, -np.inf], [np.inf, box_size, 10.0, np.inf]))
            
            # Calcul FWHM (2.355 * sigma)
            sigma_r = abs(pr[2])
            sigma_c = abs(pc[2])
            
            # Validation : sigma doit être raisonnable (entre 0.5 et 8 pixels)
            if sigma_r < 0.5 or sigma_r > 8.0 or sigma_c < 0.5 or sigma_c > 8.0:
                return None, None, None
            
            fwhm_val = (2.355*sigma_r + 2.355*sigma_c) / 2
            
            # Validation finale : FWHM doit être entre 2 et 8 pixels
            if fwhm_val < 2.0 or fwhm_val > 8.0:
                return None, None, None
            
            return fwhm_val, (r, row, pr), (c, col, pc)
        except (RuntimeError, ValueError, TypeError):
            # Échec du fit gaussien
            return None, None, None
    except Exception: 
        return None, None, None

def calculate_aperture_radii_from_fwhm(fwhm, default_fwhm=4.0):
    """
    Calcule les valeurs de départ des ouvertures photométriques basées sur le FWHM.
    
    Paramètres:
    -----------
    fwhm : float ou None
        Valeur du FWHM en pixels. Si None ou invalide, utilise default_fwhm.
    default_fwhm : float
        Valeur par défaut du FWHM si l'estimation échoue (défaut: 4.0 pixels).
    
    Retourne:
    --------
    tuple : (r_ap, r_in, r_out)
        Rayons d'ouverture : aperture, inner annulus, outer annulus (en pixels)
    """
    # Utiliser la valeur par défaut si FWHM invalide
    if fwhm is None or not np.isfinite(fwhm) or fwhm < 2.0 or fwhm > 8.0:
        logger.debug(f"FWHM invalide ({fwhm}), utilisation valeur par défaut {default_fwhm}")
        fwhm = default_fwhm
    
    # Calcul des apertures avec limites raisonnables (même logique que pour exoplanètes)
    r_ap = max(2.0, min(round(2.0*fwhm, 1), 20.0))  # Ouverture = FWHM × 2
    r_in = max(4.0, min(round(2.5*fwhm, 1),40.0))  # Entre 4 et 25 pixels
    r_out = max(6.0, min(round(3.5*fwhm, 1), 50.0))  # Entre 6 et 35 pixels
    
    # Vérification de cohérence : r_in > r_ap et r_out > r_in
    if r_in <= r_ap:
        r_in = r_ap + 2.0
    if r_out <= r_in:
        r_out = r_in + 2.0
    
    return r_ap, r_in, r_out

def _query_gaia_for_astrometry(center_coord, search_radius, mag_limit):
    logger.debug(f"Requête Gaia via Vizier avec box search (workaround timeout CIRCLE)")
    try:
        # Calcul du box autour du centre
        ra_center = center_coord.ra.deg
        dec_center = center_coord.dec.deg
        radius_deg = search_radius.to(u.deg).value
        ra_min = ra_center - radius_deg
        ra_max = ra_center + radius_deg
        dec_min = dec_center - radius_deg
        dec_max = dec_center + radius_deg

        v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag"], row_limit=5000)
        res = v.query_region(
            center_coord,
            width=2 * radius_deg * u.deg,   # largeur/hauteur du box
            height=2 * radius_deg * u.deg,
            catalog="I/355/gaiadr3"
        )
        if res and len(res[0]) > 0:
            tab = res[0]
            # Filtre manuel box + magnitude (plus sûr)
            mask = (
                (tab['RA_ICRS'] >= ra_min) & (tab['RA_ICRS'] <= ra_max) &
                (tab['DE_ICRS'] >= dec_min) & (tab['DE_ICRS'] <= dec_max)
            )
            if 'Gmag' in tab.colnames:
                valid_mask = ~tab['Gmag'].mask if hasattr(tab['Gmag'], 'mask') else ~np.isnan(tab['Gmag'])
                tab = tab[mask & valid_mask & (tab['Gmag'] <= mag_limit)]
            else:
                tab = tab[mask]
            logger.info(f"Gaia (box search) : {len(tab)} étoiles récupérées.")
            return tab
        logger.warning("Gaia : Aucune étoile trouvée.")
        return Table()
    except Exception as e:
        logger.error(f"Erreur requête Gaia (Vizier box) : {e}")
        return Table()

def _read_header_time(path):
    with open_fits_with_fixed_wcs(path) as hdul:
        header = hdul[0].header
        date_str = header.get('DATE-OBS', header.get('DATE'))
        exptime = float(header.get('EXPTIME', header.get('EXPOSURE', 0.0)))
        if not date_str:
            raise ValueError(f"Pas de date dans le header de {os.path.basename(path)}")
        date_str = date_str.replace(' ', 'T').replace('Z', '')
        if '.' in date_str:
            parts = date_str.split('.')
            date_str = f"{parts[0]}.{parts[1][:6]}"
        try:
            start_dt = datetime.fromisoformat(date_str)
        except ValueError:
            start_dt = datetime.strptime(date_str, '%Y-%m-%d')
            logger.warning(f"Heure manquante dans {path}")
        mid_dt = start_dt + timedelta(seconds=exptime / 2)
        return mid_dt, exptime

def _interpolate_ephemeris(eph_table, fits_date_time):
        eph_times_jd = eph_table['datetime_jd'].data
        obs_time = Time(fits_date_time, scale='utc')
        target_jd = obs_time.jd
        coords_eph = SkyCoord(ra=eph_table['RA'], dec=eph_table['DEC'], unit=(u.deg, u.deg))
        cart_eph = coords_eph.cartesian.xyz.value
        x_i = np.interp(target_jd, eph_times_jd, cart_eph[0])
        y_i = np.interp(target_jd, eph_times_jd, cart_eph[1])
        z_i = np.interp(target_jd, eph_times_jd, cart_eph[2])
        vec = np.array([x_i, y_i, z_i])
        norm = np.linalg.norm(vec)
        if norm == 0:
            return None
        unit_vec = vec / norm
        return SkyCoord(x=unit_vec[0], y=unit_vec[1], z=unit_vec[2], 
            representation_type='cartesian', frame='icrs')

def _interpolate_magnitude_V(eph_table, fits_date_time):
    # Vérifier si la magnitude V est disponible
    if 'V' not in eph_table.colnames:
        logger.warning("Magnitude V non disponible dans les éphémérides")
        return 99.0  # Valeur par défaut

    eph_times_jd = eph_table['datetime_jd'].data
    obs_jd = Time(fits_date_time, scale='utc').jd
    return float(np.interp(obs_jd, eph_times_jd, eph_table['V']))

def _match_sources_gaia_for_astrometry(wcs_init, sources, gaia_table, dist_limit=5.0*u.arcsec):
    if len(sources) == 0 or len(gaia_table) == 0:
        return []
    detected_sky = wcs_init.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
    gaia_sky = SkyCoord(ra=gaia_table['RA_ICRS'], dec=gaia_table['DE_ICRS'], unit=(u.deg, u.deg))
    idx, d2d, _ = match_coordinates_sky(detected_sky, gaia_sky)
    mask = d2d < dist_limit
    matches = Table()
    matches['xcentroid'] = sources['xcentroid'][mask]
    matches['ycentroid'] = sources['ycentroid'][mask]
    matches['ra_cat'] = gaia_table['RA_ICRS'][idx][mask]
    matches['dec_cat'] = gaia_table['DE_ICRS'][idx][mask]
    if 'Gmag' in gaia_table.colnames:
        matches['mag'] = gaia_table['Gmag'][idx][mask]
    return matches

def _fit_wcs_from_matches_for_astrometry(matches, wcs_init):
    if len(matches) < 3:
        return wcs_init, np.nan
    try:
        xy_pixels = np.array([matches['xcentroid'], matches['ycentroid']]).T
        world_coords = SkyCoord(ra=matches['ra_cat'], dec=matches['dec_cat'], unit=(u.deg, u.deg))
        new_wcs = fit_wcs_from_points(xy_pixels.T, world_coords, proj_point='center', sip_degree=None)
        calc_sky = new_wcs.pixel_to_world(matches['xcentroid'], matches['ycentroid'])
        sep = world_coords.separation(calc_sky)
        rms = np.sqrt(np.mean(sep.arcsec**2))
        return new_wcs, rms
    except Exception as e:
        logger.debug(f"Erreur fit WCS: {e}")
        return wcs_init, np.nan

def _detect_sources_for_astrometry(data, fwhm, threshold_sigma, max_sources):
    _, _, std = sigma_clipped_stats(data)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma*std)
    sources = daofind(data)
    if sources: 
        sources.sort('flux', reverse=True)
        return sources[:max_sources]
    return []

# =============================================================================
# MODULES D'OPTIMISATION & UTILITAIRES
# =============================================================================
# Note: GaiaCache et FastStarMatcher sont maintenant importés de core.fast_astrometry

def _interpolate_target_position(eph_table, target_jd):
    """Interpolation précise RA/Dec pour le JD cible."""
    eph_jds = eph_table['datetime_jd']
    ra_values = eph_table['RA']
    dec_values = eph_table['DEC']
    
    # Trouver les indices des points d'éphémérides autour du JD cible pour debug
    idx_before = np.searchsorted(eph_jds, target_jd) - 1
    idx_after = idx_before + 1
    
    # Vérifier les limites
    if idx_before < 0:
        idx_before = 0
    if idx_after >= len(eph_jds):
        idx_after = len(eph_jds) - 1
    
    # Log pour debug (seulement si les valeurs sont très proches)
    if idx_before < idx_after:
        jd_before = eph_jds[idx_before]
        jd_after = eph_jds[idx_after]
        ra_before = ra_values[idx_before]
        ra_after = ra_values[idx_after]
        dec_before = dec_values[idx_before]
        dec_after = dec_values[idx_after]
        
        # Si les deux points d'éphémérides ont les mêmes coordonnées, l'interpolation donnera le même résultat
        if abs(ra_before - ra_after) < 1e-10 and abs(dec_before - dec_after) < 1e-10:
            logger.debug(f"Points éphémérides identiques autour de JD={target_jd:.12f}: RA={ra_before:.12f}°, Dec={dec_before:.12f}°")
    
    # Vérifier s'il y a une discontinuité RA (passage de 359° à 1°)
    # Si oui, utiliser les coordonnées cartésiennes pour l'interpolation
    ra_diff = np.diff(ra_values)
    has_ra_wrap = np.any(np.abs(ra_diff) > 180.0)
    
    if has_ra_wrap:
        # Utiliser les coordonnées cartésiennes pour éviter les discontinuités
        coords = SkyCoord(ra=ra_values, dec=dec_values, unit=u.deg)
        xyz = coords.cartesian.xyz.value
        
        # Interpolation sur les 3 composantes cartésiennes
        x_i = np.interp(target_jd, eph_jds, xyz[0])
        y_i = np.interp(target_jd, eph_jds, xyz[1])
        z_i = np.interp(target_jd, eph_jds, xyz[2])
        
        # Convertir en RA/Dec (Astropy normalise automatiquement, mais préserve la direction)
        coord_cart = SkyCoord(x=x_i, y=y_i, z=z_i, 
                             representation_type='cartesian', frame='icrs')
        coord_sph = coord_cart.represent_as('spherical')
        return SkyCoord(ra=coord_sph.lon, dec=coord_sph.lat, frame='icrs')
    else:
        # Interpolation directe des RA/Dec (plus précise et plus rapide)
        # Utiliser une interpolation linéaire avec plus de précision
        ra_interp = float(np.interp(target_jd, eph_jds, ra_values))
        dec_interp = float(np.interp(target_jd, eph_jds, dec_values))
        
        # Log pour vérifier si les points d'éphémérides sont identiques
        if idx_before < idx_after:
            jd_before = eph_jds[idx_before]
            jd_after = eph_jds[idx_after]
            ra_before = ra_values[idx_before]
            ra_after = ra_values[idx_after]
            dec_before = dec_values[idx_before]
            dec_after = dec_values[idx_after]
            
            # Si les deux points d'éphémérides ont les mêmes coordonnées, l'interpolation donnera le même résultat
            if abs(ra_before - ra_after) < 1e-10 and abs(dec_before - dec_after) < 1e-10:
                logger.debug(f"Points éphémérides identiques autour de JD={target_jd:.12f}: RA={ra_before:.12f}°, Dec={dec_before:.12f}°")
            else:
                # Calculer le facteur d'interpolation pour vérifier
                if jd_after != jd_before:
                    alpha = (target_jd - jd_before) / (jd_after - jd_before)
                    logger.debug(f"Interpolation JD={target_jd:.12f} (alpha={alpha:.6f}): RA={ra_before:.12f}->{ra_after:.12f}°, Dec={dec_before:.12f}->{dec_after:.12f}°")
        
        return SkyCoord(ra=ra_interp * u.deg, dec=dec_interp * u.deg, frame='icrs')

class AsteroidPhotometryTab:
    def __init__(self, parent):
        self.frame = ttk.Frame(parent)
        self.frame.pack(fill=tk.BOTH, expand=True)

        # Données
        self.directory = None
        self.image_files = []
        self.current_image_index = 0  # Index de l'image courante
        self.current_image_path = None
        self.current_header = None
        self.current_data = None
        self.ephemeris_data = None
        self.wcs = None
        self.current_selections = []
        self.image_times = []
        self.last_gaia_table = None
        self.last_wcs = None

        # Position T1 sélectionnée par l'utilisateur (clic)
        self.target_t1_px = None  # (x, y) en pixels
        self.target_t1_sky = None  # SkyCoord
        self.last_t1_px_for_fwhm = None  # Dernière position T1 utilisée pour FWHM
        self.manual_t1_anchor_first = None  # {"jd": float, "ra_deg": float, "dec_deg": float}
        self.manual_t1_anchor_last = None   # {"jd": float, "ra_deg": float, "dec_deg": float}

        # Cache Gaia persistant
        self.gaia_cache = GaiaCache()

        # Variables Tkinter (Synchronisées avec votre demande)
        self.asteroid_id_var = tk.StringVar(value="")
        self.observatory_code_var = tk.StringVar(value="")
        self.ephemeris_step_var = tk.StringVar(value="2m")  # Pas d'éphémérides par défaut: 2 minutes
        self.fwhm_var = tk.DoubleVar(value=5.0)
        self.threshold_sigma_var = tk.DoubleVar(value=2.0)
        self.max_sources_var = tk.IntVar(value=1000)
        self.gaia_mag_limit_var = tk.DoubleVar(value=18.0)
        self.match_radius_var = tk.DoubleVar(value=10.0)
        self.centroid_box_size_var = tk.IntVar(value=15)

        self.aperture_radius_var = tk.DoubleVar(value=6.0)
        self.annulus_inner_var = tk.DoubleVar(value=12.0)
        self.annulus_outer_var = tk.DoubleVar(value=20.0)

        # Option GPU
        self.use_gpu_var = tk.BooleanVar(value=False)
        
        # Méthode d'astrométrie
        self.astrometry_method_var = tk.StringVar(value="zero_aperture")  # "zero_aperture" ou "classical"
        
        # Défilement automatique
        self.auto_play_active = False
        self.auto_play_delay = 25  # Délai en millisecondes entre les images (25 ms par défaut)
        self.auto_play_job = None

        # Série zéro-ouverture (6 lignes par image : une par ouverture)
        self.zero_aperture_series_by_image = {}
        self.zero_aperture_rows = []
        self.zero_aperture_window = None
        self.zero_aperture_tree = None
        self.photometry_manual_apertures = {}
        self.photometry_filter_used = None
        self.photometry_delta_to_gaia_g = 0.0
        self.photometry_settings_confirmed = False
        self.photometry_gaia_mag_cache = {}
        self.photometry_filter_var = tk.StringVar(value="G")
        self.photometry_delta_var = tk.DoubleVar(value=0.0)
        
        # Console de logs live dans l'onglet
        self.log_text_widget = None
        self.log_text_handler = None
        self.log_process_depth = 0

        self.setup_gui()

    def setup_gui(self):
        ttk.Label(self.frame, text="Photométrie d'astéroïdes", font=("Helvetica", 12, "bold")).pack(pady=5)

        # Conteneur gauche : Canvas + scrollbar pour zone de contrôle défilante
        left_container = ttk.Frame(self.frame)
        left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(10, 0), pady=10)

        control_canvas = tk.Canvas(left_container, highlightthickness=0)
        control_scrollbar = ttk.Scrollbar(left_container)

        control_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        control_canvas.pack(side=tk.LEFT, fill=tk.Y)

        control_canvas.config(yscrollcommand=control_scrollbar.set)
        control_scrollbar.config(command=control_canvas.yview)

        # Cadre intérieur (contenu défilant) ; parent = canvas pour create_window
        control_frame = ttk.Frame(control_canvas)
        canvas_window = control_canvas.create_window(0, 0, window=control_frame, anchor="nw")

        def _on_control_frame_configure(event):
            control_canvas.configure(scrollregion=control_canvas.bbox("all"))

        def _on_canvas_configure(event):
            control_canvas.itemconfig(canvas_window, width=event.width)

        control_frame.bind("<Configure>", _on_control_frame_configure)
        control_canvas.bind("<Configure>", _on_canvas_configure)

        def _on_mousewheel(event):
            control_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _bind_mousewheel(event):
            control_canvas.bind_all("<MouseWheel>", _on_mousewheel)

        def _unbind_mousewheel(event):
            control_canvas.unbind_all("<MouseWheel>")

        control_canvas.bind("<Enter>", _bind_mousewheel)
        control_canvas.bind("<Leave>", _unbind_mousewheel)
        control_frame.bind("<Enter>", _bind_mousewheel)
        control_frame.bind("<Leave>", _unbind_mousewheel)

        display_frame = ttk.Frame(self.frame)
        display_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Zone de droite en 2 panneaux (50/50) : visualisation + log live
        right_paned = tk.PanedWindow(display_frame, orient=tk.VERTICAL, sashrelief=tk.RAISED, sashwidth=6)
        right_paned.pack(fill=tk.BOTH, expand=True)
        viz_frame = ttk.Frame(right_paned)
        log_frame = ttk.LabelFrame(right_paned, text="Journal du processus")
        right_paned.add(viz_frame, minsize=220)
        right_paned.add(log_frame, minsize=180)

       # ttk.Button(control_frame, text="📁 Charger images", command=self.load_folder).pack(fill=tk.X, pady=2)

        # Navigation entre images
        nav_frame = ttk.LabelFrame(control_frame, text="Navigation Images")
        nav_frame.pack(fill=tk.X, pady=2, ipady=2)
        ttk.Button(nav_frame, text="📁 Charger images", command=self.load_folder).pack(fill=tk.X, pady=1)
        
        # Slider pour naviguer entre les images
        slider_frame = ttk.Frame(nav_frame)
        slider_frame.pack(fill=tk.X, padx=5, pady=2)
        
        # Slider en base 1 côté UI (première image = 1)
        self.image_slider_var = tk.IntVar(value=1)
        self.image_slider = tk.Scale(slider_frame, 
                                     from_=1, 
                                     to=1, 
                                     orient=tk.HORIZONTAL,
                                     variable=self.image_slider_var,
                                     command=self.on_slider_change,
                                     length=200)
        self.image_slider.pack(fill=tk.X, padx=5)
        
        # Contrôles de navigation
        nav_buttons_frame = ttk.Frame(nav_frame)
        nav_buttons_frame.pack(fill=tk.X, padx=5, pady=1)
        
        ttk.Button(nav_buttons_frame, text="⏮", command=self.load_first_image, width=3).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_buttons_frame, text="◀", command=self.load_previous_image, width=3).pack(side=tk.LEFT, padx=2)
        
        # Bouton play/pause pour défilement automatique
        self.play_pause_button = ttk.Button(nav_buttons_frame, text="▶", command=self.toggle_auto_play, width=3)
        self.play_pause_button.pack(side=tk.LEFT, padx=2)
        
        ttk.Button(nav_buttons_frame, text="▶", command=self.load_next_image, width=3).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_buttons_frame, text="⏭", command=self.load_last_image, width=3).pack(side=tk.LEFT, padx=2)
        
        # Contrôle de la vitesse de défilement
        speed_frame = ttk.Frame(nav_frame)
        speed_frame.pack(fill=tk.X, padx=5, pady=1)
        ttk.Label(speed_frame, text="Vitesse (ms):", font=("Helvetica", 8)).pack(side=tk.LEFT)
        self.auto_play_delay_var = tk.IntVar(value=self.auto_play_delay)
        delay_entry = ttk.Entry(speed_frame, textvariable=self.auto_play_delay_var, width=6)
        delay_entry.pack(side=tk.LEFT, padx=2)
        
        self.image_info_label = ttk.Label(nav_frame, text="Aucune image chargée", foreground="gray")
        self.image_info_label.pack(pady=1)

        id_obs_frame = ttk.Frame(control_frame)
        id_obs_frame.pack(pady=2)

        f1 = ttk.Frame(id_obs_frame); f1.pack(side=tk.LEFT, padx=2)
        ttk.Label(f1, text="ID Cible").pack()
        ttk.Entry(f1, width=10, textvariable=self.asteroid_id_var).pack()

        f2 = ttk.Frame(id_obs_frame); f2.pack(side=tk.LEFT, padx=2)
        ttk.Label(f2, text="Code Obs.").pack()
        ttk.Entry(f2, width=6, textvariable=self.observatory_code_var).pack()

        # Frame pour le pas d'éphémérides avec explication
        eph_step_frame = ttk.LabelFrame(control_frame, text="Pas Éphémérides")
        eph_step_frame.pack(fill=tk.X, pady=2, ipady=2)
        
        step_entry_frame = ttk.Frame(eph_step_frame)
        step_entry_frame.pack(fill=tk.X, padx=5, pady=1)
        ttk.Label(step_entry_frame, text="Pas:").pack(side=tk.LEFT, padx=(0, 5))
        step_entry = ttk.Entry(step_entry_frame, textvariable=self.ephemeris_step_var, width=8)
        step_entry.pack(side=tk.LEFT)
        ttk.Label(step_entry_frame, text="(ex: 1m, 2m, 5m, 10m)", font=("Helvetica", 8), foreground="gray").pack(side=tk.LEFT, padx=(5, 0))
        
        # Explication du pas
        explanation = "Plus le pas est petit, plus la résolution est fine:\n" \
                     "• 1m = ~180 points pour 3h (très précis)\n" \
                     "• 2m = ~90 points pour 3h (recommandé)\n" \
                     "• 5m = ~36 points pour 3h (rapide)\n" \
                     "• 10m = ~18 points pour 3h (moins précis)"
        help_label = ttk.Label(eph_step_frame, text=explanation, font=("Helvetica", 8), 
                               foreground="blue", justify=tk.LEFT, wraplength=250)
        help_label.pack(padx=5, pady=1, anchor="w")

        ttk.Button(control_frame, text="🔭 Récupérer Éphémérides", command=self.fetch_ephemeris).pack(fill=tk.X, pady=1)

        ast_frame = ttk.LabelFrame(control_frame, text="Paramètres Astrométrie")
        ast_frame.pack(fill=tk.X, pady=2, ipady=2)

        params = [("FWHM [px]", self.fwhm_var), ("Seuil [σ]", self.threshold_sigma_var),
                  ("Gaia Gmax", self.gaia_mag_limit_var), ("Match [\"]", self.match_radius_var)]
        
        for i, (txt, var) in enumerate(params):
            ttk.Label(ast_frame, text=txt).grid(row=i, column=0, sticky="w", padx=5, pady=0)
            ttk.Entry(ast_frame, textvariable=var, width=6).grid(row=i, column=1, padx=5, pady=0)
        
        # Option GPU
        gpu_status = "✅ Disponible" if HAS_GPU else "❌ Non disponible"
        ttk.Label(ast_frame, text=f"GPU: {gpu_status}", foreground="green" if HAS_GPU else "gray", font=("Helvetica", 8)).grid(row=4, column=0, sticky="w", padx=5, pady=0)
        ttk.Checkbutton(ast_frame, text="Utiliser GPU", variable=self.use_gpu_var, state="normal" if HAS_GPU else "disabled").grid(row=4, column=1, padx=5, pady=0)

        # Sélection de la méthode d'astrométrie
        ttk.Label(ast_frame, text="Méthode:").grid(row=5, column=0, sticky="w", padx=5, pady=1)
        
        # Création d'un combobox avec valeurs lisibles et mapping interne
        method_display_map = {
            "Zero-Aperture (extrapolation)": "zero_aperture",
            "Classique (FWHM)": "classical"
        }
        method_display_values = list(method_display_map.keys())
        
        method_combo = ttk.Combobox(ast_frame, 
                                   values=method_display_values,
                                   state="readonly", width=25)
        method_combo.current(0)  # Sélectionner "Zero-Aperture" par défaut
        method_combo.grid(row=5, column=1, padx=5, pady=1, sticky="w")
        
        # Synchroniser la sélection avec la variable interne
        def on_method_change(event=None):
            display_value = method_combo.get()
            if display_value in method_display_map:
                self.astrometry_method_var.set(method_display_map[display_value])
        
        method_combo.bind("<<ComboboxSelected>>", on_method_change)
        on_method_change()  # Initialisation pour synchroniser la valeur par défaut
        
        # Labels explicatifs pour chaque méthode
        method_labels_frame = ttk.Frame(ast_frame)
        method_labels_frame.grid(row=6, column=0, columnspan=2, padx=5, pady=0, sticky="w")
        
        zero_ap_label = ttk.Label(method_labels_frame, 
                                  text="• Zero-Aperture: Extrapolation (plusieurs ouvertures, plus précis)",
                                  font=("Helvetica", 7), foreground="blue")
        zero_ap_label.pack(anchor="w")
        
        classical_label = ttk.Label(method_labels_frame, 
                                   text="• Classique: Une ouverture optimisée FWHM (plus rapide)",
                                   font=("Helvetica", 7), foreground="green")
        classical_label.pack(anchor="w")

        # Bouton unique pour lancer l'astrométrie
        ttk.Button(ast_frame, text="🧭 Lancer l'Astrométrie", command=self.run_single_astrometry_with_method).grid(row=7, column=0, columnspan=2, pady=2)
        ttk.Button(ast_frame, text="🧭 Batch Astrométrie (Thread)", command=self._batch_astrometry_worker).grid(row=8, column=0, columnspan=2, pady=2)
        ttk.Button(ast_frame, text="📋 Tableau ZA (6 ouvertures)", command=self.open_zero_aperture_table).grid(row=9, column=0, columnspan=2, pady=2)

        photo_frame = ttk.LabelFrame(control_frame, text="Photométrie d'Ouverture")
        photo_frame.pack(fill=tk.X, pady=2, ipady=2)

        # Réglages photométrie Gaia intégrés en tête du cadre Photométrie d'Ouverture
        ttk.Label(photo_frame, text="Filtre utilisé:").grid(row=0, column=0, sticky="w", padx=5, pady=(2, 1))
        ttk.Entry(photo_frame, textvariable=self.photometry_filter_var, width=8).grid(row=0, column=1, sticky="w", padx=4, pady=(2, 1))
        ttk.Label(photo_frame, text="Δ(G-filtre) [mag]:").grid(row=1, column=0, sticky="w", padx=5, pady=1)
        ttk.Entry(photo_frame, textvariable=self.photometry_delta_var, width=8).grid(row=1, column=1, sticky="w", padx=4, pady=1)
        ttk.Label(photo_frame, text="mag_G = mag_mesurée + Δ", font=("Helvetica", 8), foreground="gray").grid(
            row=1, column=2, columnspan=2, sticky="w", padx=6, pady=1
        )
        ttk.Button(photo_frame, text="Appliquer", command=self._apply_photometry_gaia_settings_from_ui).grid(
            row=2, column=0, sticky="w", padx=5, pady=(1, 2)
        )
        ttk.Button(photo_frame, text="Réinitialiser", command=self._reset_photometry_gaia_settings).grid(
            row=2, column=1, sticky="w", padx=4, pady=(1, 2)
        )

        self.photo_button = ttk.Button(photo_frame, text="⭕ SET-UP PHOTOMÉTRIE ", command=self.run_photometry_setup)
        self.photo_button.grid(row=3, column=0, columnspan=4, pady=2)

        self.photo_single_button = ttk.Button(
            photo_frame,
            text="📸 PHOTOMÉTRIE IMAGE COURANTE (Comètes)",
            command=self.run_photometry_single_image
        )
        self.photo_single_button.grid(row=4, column=0, columnspan=4, pady=2)

        self.photo_batch_button = ttk.Button(
            photo_frame,
            text="📊 PHOTOMÉTRIE BATCH (Astéroïdes)",
            command=self.run_photometry_batch
        )
        self.photo_batch_button.grid(row=5, column=0, columnspan=4, pady=2)

        # Bouton KBMOD sorti du sous-cadre Photométrie d'Ouverture, placé juste en dessous
        self.kbmod_button = ttk.Button(
            control_frame, text="🔍 Détection KBMOD (via WSL)",
            command=self._open_kbmod_wsl_detection_dialog
        )
        self.kbmod_button.pack(fill=tk.X, pady=1)

        # Modélisation & Inversion (courbe de lumière → période, amplitude, etc.)
        model_frame = ttk.LabelFrame(control_frame, text="Modélisation & Inversion")
        model_frame.pack(fill=tk.X, pady=2, ipady=2)
        ttk.Label(
            model_frame, text="Courbe de lumière → modèle rotation + ajustement",
            font=("Helvetica", 8), foreground="gray"
        ).pack(anchor="w", padx=5, pady=0)
        ttk.Button(
            model_frame, text="📈 Modélisation & Inversion",
            command=self._open_lightcurve_modeling_window
        ).pack(fill=tk.X, padx=5, pady=2)

        # Forme 3D (DAMIT : inversion de courbes de lumière → modèle polyédrique)
        damit_frame = ttk.LabelFrame(control_frame, text="Forme (DAMIT)")
        damit_frame.pack(fill=tk.X, pady=2, ipady=2)
        ttk.Label(
            damit_frame,
            text="DAMIT : forme 3D, période, axe de spin. Voir doc pour soumettre.",
            font=("Helvetica", 8),
            foreground="gray",
            justify=tk.LEFT,
        ).pack(anchor="w", padx=5, pady=0)
        btn_doc = ttk.Button(
            damit_frame, text="Documentation DAMIT",
            command=lambda: webbrowser.open("https://damit.cuni.cz/projects/damit/pages/documentation")
        )
        btn_doc.pack(fill=tk.X, padx=5, pady=1)
        btn_browse = ttk.Button(
            damit_frame, text="Parcourir les modèles DAMIT",
            command=lambda: webbrowser.open("https://astro.troja.mff.cuni.cz/projects/damit/asteroids/browse")
        )
        btn_browse.pack(fill=tk.X, padx=5, pady=1)
        ttk.Button(
            damit_frame, text="Charger modèle 3D (fichier .obj / shape.txt)",
            command=self._open_damit_shape_window
        ).pack(fill=tk.X, padx=5, pady=2)

        # Visualisation (fenêtre agrandie)
        self.fig, self.ax = plt.subplots(figsize=(12, 9))  # taille plus grande pour la visualisation
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, viz_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Console de log live (zone libérée sous la visualisation)
        self.log_text_widget = scrolledtext.ScrolledText(
            log_frame,
            wrap=tk.WORD,
            height=12,
            state="disabled",
            font=("Consolas", 9),
        )
        self.log_text_widget.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self._append_log_panel_message("Journal inactif. Il s'active uniquement pendant astrométrie/photométrie.")
        
        # Positionner le séparateur vers 50% de la hauteur
        def _place_half_sash():
            try:
                total_h = max(display_frame.winfo_height(), 1)
                right_paned.sash_place(0, 0, int(total_h * 0.5))
            except Exception:
                pass
        display_frame.after(150, _place_half_sash)
        
        # Connexion du gestionnaire de clic pour sélectionner T1
        self.canvas.mpl_connect("button_press_event", self.on_image_click)
    
    def _attach_live_log_view(self):
        """Branche le logger racine sur la console de logs de l'onglet."""
        if not self.log_text_widget or self.log_text_handler is not None:
            return
        try:
            handler = UILogTextHandler(self.log_text_widget)
            handler.setLevel(logging.INFO)
            handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
            logging.getLogger().addHandler(handler)
            self.log_text_handler = handler
            self.frame.bind("<Destroy>", self._on_tab_destroy, add="+")
        except Exception as e:
            logger.debug(f"Impossible d'attacher la vue log live: {e}")
    
    def _append_log_panel_message(self, msg: str):
        """Ajoute un message local dans le panneau de log (sans logger global)."""
        if not self.log_text_widget:
            return
        try:
            self.log_text_widget.configure(state="normal")
            self.log_text_widget.insert("end", msg + "\n")
            self.log_text_widget.see("end")
            self.log_text_widget.configure(state="disabled")
        except Exception:
            pass
    
    def _start_process_log_capture(self, process_label: str):
        """Active la capture live uniquement pendant les traitements demandés."""
        if self.log_process_depth == 0:
            self._attach_live_log_view()
            self._append_log_panel_message(f"=== Début {process_label} ===")
        self.log_process_depth += 1
    
    def _stop_process_log_capture(self, process_label: str):
        """Désactive la capture live quand le dernier traitement se termine."""
        if self.log_process_depth > 0:
            self.log_process_depth -= 1
        if self.log_process_depth == 0:
            self._append_log_panel_message(f"=== Fin {process_label} ===")
            self._detach_live_log_view()
    
    def _on_tab_destroy(self, event):
        """Nettoie le handler de log quand l'onglet est détruit."""
        if event.widget is self.frame:
            self._detach_live_log_view()
    
    def _detach_live_log_view(self):
        if self.log_text_handler is None:
            return
        try:
            logging.getLogger().removeHandler(self.log_text_handler)
        except Exception:
            pass
        self.log_text_handler = None

    def load_folder(self):
        from tkinter import filedialog
        directory = filedialog.askdirectory()
        if not directory: return
        
        self.directory = directory
        self.image_files = sorted([str(f) for f in Path(directory).glob("*.fits")])
        logger.info(f"Dossier : {len(self.image_files)} fichiers trouvés.")
        if self.image_files: 
            # Mettre à jour le slider
            self.update_image_slider()
            self.current_image_index = 0
            self.load_image(self.image_files[0])
            self.update_image_info()

    def load_image(self, path):
        logger.info(f"Chargement image : {Path(path).name}")
        self.current_image_path = path
        
        # Mettre à jour l'index et le slider (si le slider existe)
        if path in self.image_files:
            self.current_image_index = self.image_files.index(path)
            # Mettre à jour le slider sans déclencher l'événement
            if hasattr(self, 'image_slider_var'):
                self.image_slider_var.set(self.current_image_index + 1)
        
        with open_fits_with_fixed_wcs(path) as hdul:
            self.current_header = hdul[0].header
            self.current_data = hdul[0].data.astype(float)
            self.wcs = WCS(self.current_header)
        
        # Essayer de positionner T1 automatiquement depuis les éphémérides si disponibles
        # Sinon, réinitialiser T1 pour la nouvelle image (sera défini par clic gauche)
        self.target_t1_px = None
        self.target_t1_sky = None
        
        # Mode manuel : ne pas positionner T1 automatiquement depuis éphémérides
        
        # Réinitialiser le mode de la toolbar pour s'assurer que les clics fonctionnent
        if hasattr(self.toolbar, 'mode'):
            self.toolbar.mode = ''
        if hasattr(self.toolbar, '_active'):
            self.toolbar._active = None
        
        self.refresh_display()  # refresh_display() appelle déjà canvas.draw()
        self.update_image_info()
    
    def update_image_slider(self):
        """Met à jour le slider pour refléter le nombre d'images disponibles."""
        num_images = len(self.image_files)
        if num_images > 0:
            self.image_slider.config(from_=1, to=num_images)
            self.image_slider.config(state="normal")
        else:
            self.image_slider.config(from_=1, to=1)
            self.image_slider.config(state="disabled")
    
    def on_slider_change(self, value):
        """Gère le changement de position du slider pour charger l'image correspondante."""
        if not self.image_files:
            return
        
        try:
            # Slider affiché en base 1, index interne en base 0
            display_idx = int(float(value))
            index = display_idx - 1
            if 0 <= index < len(self.image_files):
                # Ne charger que si l'index a vraiment changé (évite les boucles infinies)
                if index != self.current_image_index:
                    self.current_image_index = index
                    self.load_image(self.image_files[index])
                    self.update_image_info()
        except (ValueError, IndexError) as e:
            logger.debug(f"Erreur lors du changement de slider : {e}")
    
    def load_first_image(self):
        """Charge la première image du dossier."""
        if not self.image_files:
            messagebox.showwarning("Avertissement", "Aucune image chargée. Utilisez 'Charger images' d'abord.")
            return
        self.stop_auto_play()
        self.current_image_index = 0
        self.image_slider_var.set(1)
        self.load_image(self.image_files[0])

    def load_previous_image(self):
        """Charge l'image précédente."""
        if not self.image_files:
            return
        self.stop_auto_play()
        if self.current_image_index > 0:
            self.current_image_index -= 1
            self.image_slider_var.set(self.current_image_index + 1)
            self.load_image(self.image_files[self.current_image_index])
        else:
            # Aller à la dernière image si on est à la première
            self.current_image_index = len(self.image_files) - 1
            self.image_slider_var.set(self.current_image_index + 1)
            self.load_image(self.image_files[self.current_image_index])
    
    def load_next_image(self):
        """Charge l'image suivante."""
        if not self.image_files:
            return
        self.stop_auto_play()
        if self.current_image_index < len(self.image_files) - 1:
            self.current_image_index += 1
            self.image_slider_var.set(self.current_image_index + 1)
            self.load_image(self.image_files[self.current_image_index])
        else:
            # Revenir à la première image si on est à la dernière
            self.current_image_index = 0
            self.image_slider_var.set(1)
            self.load_image(self.image_files[0])
    
    def toggle_auto_play(self):
        """Démarre ou arrête le défilement automatique des images."""
        if not self.image_files:
            messagebox.showwarning("Avertissement", "Aucune image chargée. Utilisez 'Charger images' d'abord.")
            return
        
        if self.auto_play_active:
            self.stop_auto_play()
        else:
            self.start_auto_play()
    
    def start_auto_play(self):
        """Démarre le défilement automatique."""
        if not self.image_files:
            return

        self.auto_play_active = True
        self.play_pause_button.config(text="⏸")
        self.auto_play_delay = self.auto_play_delay_var.get()
        if self.auto_play_delay < 25:  # Minimum 25ms
            self.auto_play_delay = 25
            self.auto_play_delay_var.set(25)
        
        # Démarrer la boucle de défilement
        self._auto_play_step()
    
    def stop_auto_play(self):
        """Arrête le défilement automatique."""
        self.auto_play_active = False
        if hasattr(self, 'play_pause_button'):
            self.play_pause_button.config(text="▶")
        if self.auto_play_job:
            self.frame.after_cancel(self.auto_play_job)
            self.auto_play_job = None
    
    def _auto_play_step(self):
        """Étape du défilement automatique (appelée récursivement)."""
        if not self.auto_play_active or not self.image_files:
            self.stop_auto_play()
            return
        
        # Passer à l'image suivante
        if self.current_image_index < len(self.image_files) - 1:
            self.current_image_index += 1
        else:
            # Revenir au début si on est à la fin
            self.current_image_index = 0
        
        # Mettre à jour le slider et charger l'image
        self.image_slider_var.set(self.current_image_index + 1)
        self.load_image(self.image_files[self.current_image_index])
        
        # Programmer la prochaine étape
        self.auto_play_job = self.frame.after(self.auto_play_delay, self._auto_play_step)
    
    def load_last_image(self):
        """Charge la dernière image du dossier."""
        if not self.image_files:
            messagebox.showwarning("Avertissement", "Aucune image chargée. Utilisez 'Charger images' d'abord.")
            return
        self.stop_auto_play()
        try:
            self.current_image_index = len(self.image_files) - 1
            self.image_slider_var.set(self.current_image_index + 1)
            self.load_image(self.image_files[-1])
            # Forcer la mise à jour de l'affichage et s'assurer que le canvas est actualisé
            self.refresh_display()
            self.update_image_info()
            logger.info(f"Dernière image chargée : {Path(self.image_files[-1]).name}")
        except Exception as e:
            logger.error(f"Erreur lors du chargement de la dernière image : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Impossible de charger la dernière image : {e}")
    
    def update_image_info(self):
        """Met à jour le label d'information sur l'image courante."""
        if not self.image_files or not self.current_image_path:
            if hasattr(self, 'image_info_label'):
                self.image_info_label.config(text="Aucune image chargée", foreground="gray")
            return

        try:
            current_idx = self.image_files.index(self.current_image_path)
            total = len(self.image_files)
            image_name = Path(self.current_image_path).name
            if hasattr(self, 'image_info_label'):
                self.image_info_label.config(
                    text=f"Image {current_idx + 1}/{total}\n{image_name}", 
                    foreground="black"
                )
        except ValueError:
            if hasattr(self, 'image_info_label'):
                self.image_info_label.config(text="Image inconnue", foreground="orange")
    
    def refresh_display(self, target_px=None):
        """Affiche l'image avec ZScale et optionnellement un cercle sur T1 (sans zoom)."""
        if self.current_data is None:
            logger.warning("Aucune donnée image à afficher")
            return

        self.ax.clear()
        
        # Contraste ZScale automatique
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(self.current_data)
        
        self.ax.imshow(self.current_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
        
        # Pas de zoom - affichage de l'image complète
        ny, nx = self.current_data.shape
        self.ax.set_xlim(0, nx)
        self.ax.set_ylim(0, ny)
        
        # Afficher T1 si sélectionné (priorité à target_px passé en paramètre, sinon self.target_t1_px)
        display_t1_px = target_px if target_px is not None else self.target_t1_px
        
        if display_t1_px:
            tx, ty = display_t1_px
            # Vérification que les coordonnées sont valides
            if np.isfinite(tx) and np.isfinite(ty) and 0 <= tx < nx and 0 <= ty < ny:
                # Dessin du cercle rouge (tx, ty) correspond à (x, y) matplotlib avec origin='lower'
                circ = plt.Circle((tx, ty), 15, color='red', fill=False, lw=2, label="T1")
                self.ax.add_patch(circ)
                logger.info(f"Cercle tracé sur T1 à [{tx:.1f}, {ty:.1f}]")
            else:
                logger.warning(f"Position T1 hors limites : [{tx:.1f}, {ty:.1f}] (image: {nx}x{ny})")
        
        # Toujours mettre à jour le canvas pour afficher l'image
        self.canvas.draw()
        self.canvas.flush_events()
        self.frame.update_idletasks()

    def _update_t1_from_ephemeris(self):
        """Met à jour T1 depuis les éphémérides pour l'image courante."""
        if not self.ephemeris_data or not self.wcs or not self.current_header:
            return

        try:
            # Calcul du temps d'observation (milieu d'exposition)
            date_obs_str = self.current_header.get('DATE-OBS', self.current_header.get('DATE'))
            if not date_obs_str:
                return
            
            exptime = float(self.current_header.get('EXPTIME', self.current_header.get('EXPOSURE', 0.0)))
            t_start = Time(date_obs_str, scale='utc')
            t_img = t_start + exptime / 2.0 * u.second
            
            logger.info(f"Mise à jour T1 depuis éphémérides pour {Path(self.current_image_path).name} : {t_img.isot}")
            
            # Interpolation de la position T1 depuis les éphémérides
            target_sky_cart = _interpolate_target_position(self.ephemeris_data, t_img.jd)
            if target_sky_cart is None:
                logger.warning("Impossible d'interpoler T1 depuis les éphémérides")
                return
            
            # Conversion en représentation sphérique pour obtenir RA/Dec
            target_sky_spherical = target_sky_cart.represent_as('spherical')
            target_sky = SkyCoord(ra=target_sky_spherical.lon, dec=target_sky_spherical.lat, frame='icrs')
            
            # Conversion en pixels
            tx, ty = self.wcs.world_to_pixel(target_sky)

            # Vérifier que T1 est dans l'image
            ny, nx = self.current_data.shape
            if 0 <= tx < nx and 0 <= ty < ny:
                self.target_t1_px = (tx, ty)
                self.target_t1_sky = target_sky
                logger.info(f"T1 mis à jour : RA={target_sky.ra.deg:.6f}°, Dec={target_sky.dec.deg:.6f}°, pixels=[{tx:.1f}, {ty:.1f}]")
            else:
                logger.warning(f"T1 hors champ pour cette image : [{tx:.1f}, {ty:.1f}] (image: {nx}x{ny})")
                self.target_t1_px = None
                self.target_t1_sky = None

        except Exception as e:
            logger.error(f"Erreur lors de la mise à jour de T1 depuis éphémérides : {e}", exc_info=True)
            self.target_t1_px = None
            self.target_t1_sky = None

    def on_image_click(self, event):
        """Gestionnaire de clic sur l'image pour sélectionner T1."""
        logger.debug(f"Événement clic reçu : inaxes={event.inaxes}, button={event.button}")
        
        if event.inaxes != self.ax:
            return

        if self.current_data is None or self.wcs is None:
            logger.warning("Clic ignoré : aucune image chargée ou WCS manquant")
            messagebox.showwarning("Avertissement", "Aucune image chargée ou WCS manquant")
            return

        toolbar_mode = getattr(self.toolbar, 'mode', None) or getattr(self.toolbar, '_active', None)
        if toolbar_mode and str(toolbar_mode).upper() in ['ZOOM', 'PAN', 'ZOOM RECT', 'PAN/ZOOM']:
            return

        if event.xdata is None or event.ydata is None:
            return

        px_click, py_click = event.xdata, event.ydata
        ny, nx = self.current_data.shape
        
        if not (0 <= px_click < nx and 0 <= py_click < ny):
            return

        # Clic droit: diagnostic centroïde + profils gaussiens, puis validation optionnelle de T1
        if event.button == 3:
            self._open_t1_centroid_diagnostic(px_click, py_click)
            return

        if event.button != 1:
            return
        
        logger.info(f"Clic détecté sur {Path(self.current_image_path).name} à [{px_click:.1f}, {py_click:.1f}]")
        
        try:
            # 1. Affiner le centroïde
            box_size = self.centroid_box_size_var.get()
            px_refined, py_refined = refine_centroid(self.current_data, px_click, py_click, box_size=box_size)
            self._finalize_t1_selection(px_refined, py_refined, show_message=True, context_label="clic gauche")

        except Exception as e:
            logger.error(f"Erreur lors de la sélection de T1 : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Impossible de sélectionner T1 : {e}")

    def _finalize_t1_selection(self, px_refined, py_refined, show_message=True, context_label="sélection"):
        """Finalise la sélection T1 à partir d'un centroïde déjà raffiné."""
        # 1) Estimation FWHM et apertures
        fwhm_result = estimate_fwhm_marginal(self.current_data, px_refined, py_refined, box_size=25)
        fwhm_val = None
        if fwhm_result:
            try:
                fwhm_candidate = float(fwhm_result[0])
                if np.isfinite(fwhm_candidate):
                    fwhm_val = fwhm_candidate
            except (TypeError, ValueError):
                fwhm_val = None

        r_ap, r_in, r_out = calculate_aperture_radii_from_fwhm(fwhm_val, default_fwhm=4.0)
        self.aperture_radius_var.set(r_ap)
        self.annulus_inner_var.set(r_in)
        self.annulus_outer_var.set(r_out)

        # 2) Conversion en coordonnées célestes
        px_scalar = float(px_refined)
        py_scalar = float(py_refined)
        target_sky_raw = self.wcs.pixel_to_world(px_scalar, py_scalar)
        if isinstance(target_sky_raw, SkyCoord):
            target_sky = target_sky_raw
        else:
            try:
                ra_deg, dec_deg = self.wcs.pixel_to_world_values(px_scalar, py_scalar)
                if not np.isfinite(ra_deg) or not np.isfinite(dec_deg):
                    raise ValueError("Coordonnées WCS non finies")
                if dec_deg < -90 or dec_deg > 90:
                    raise ValueError("Déc hors limites (-90..90)")
                target_sky = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
            except Exception:
                raise ValueError("WCS invalide. Lancez l'astrométrie avant de sélectionner T1.")

        # 3) Enregistrer ancres manuelles (première / dernière image)
        if self.image_files and self.current_image_path in self.image_files and self.current_header:
            try:
                idx = self.image_files.index(self.current_image_path)
                date_obs = self.current_header.get('DATE-OBS', self.current_header.get('DATE'))
                if date_obs:
                    exptime = float(self.current_header.get('EXPTIME', self.current_header.get('EXPOSURE', 0.0)))
                    t_start = Time(date_obs, scale='utc')
                    t_img = t_start + exptime / 2.0 * u.second
                    anchor = {
                        "jd": float(t_img.jd),
                        "ra_deg": float(target_sky.ra.deg),
                        "dec_deg": float(target_sky.dec.deg)
                    }
                    if idx == 0:
                        self.manual_t1_anchor_first = anchor
                    if idx == len(self.image_files) - 1:
                        self.manual_t1_anchor_last = anchor
            except Exception as e:
                logger.debug(f"Ancre manuelle T1 échouée : {e}")

        # 4) Stockage + affichage
        ra_val = target_sky.ra
        dec_val = target_sky.dec
        ra_str = ra_val.to_string(u.hour, precision=3)
        dec_str = dec_val.to_string(unit=u.deg, precision=2)

        self.target_t1_px = (px_scalar, py_scalar)
        self.target_t1_sky = target_sky
        self.refresh_display()

        logger.info(
            f"T1 sélectionné ({context_label}) : px=[{px_scalar:.2f},{py_scalar:.2f}], "
            f"RA={ra_val.deg:.6f}°, Dec={dec_val.deg:.6f}°"
        )
        if show_message:
            fwhm_msg = (
                f"\nFWHM estimé : {fwhm_val:.2f} px\nOuvertures : {r_ap:.1f} / {r_in:.1f} / {r_out:.1f} px"
                if fwhm_val is not None else ""
            )
            messagebox.showinfo(
                "Succès",
                f"T1 sélectionné à [{px_scalar:.1f}, {py_scalar:.1f}]\nRA: {ra_str}\nDec: {dec_str}{fwhm_msg}"
            )

    def _open_t1_centroid_diagnostic(self, px_click, py_click):
        """Clic droit: ouvre un diagnostic centroïde (cutout + profils gaussiens)."""
        try:
            box_size = self.centroid_box_size_var.get()
            px_refined, py_refined = refine_centroid(self.current_data, px_click, py_click, box_size=box_size)
            fwhm_val, row_fit, col_fit = estimate_fwhm_marginal(self.current_data, px_refined, py_refined, box_size=25)

            # Cutout centré sur la position raffinée
            half = max(7, box_size // 2)
            x0, y0 = int(round(px_refined)), int(round(py_refined))
            y1, y2 = max(0, y0 - half), min(self.current_data.shape[0], y0 + half + 1)
            x1, x2 = max(0, x0 - half), min(self.current_data.shape[1], x0 + half + 1)
            cutout = self.current_data[y1:y2, x1:x2]

            win = Toplevel(self.frame.winfo_toplevel())
            win.title("Diagnostic centroïde T1 (clic droit)")
            win.geometry("980x460")
            win.transient(self.frame.winfo_toplevel())

            fig = Figure(figsize=(11, 4.2))
            ax0 = fig.add_subplot(131)
            ax1 = fig.add_subplot(132)
            ax2 = fig.add_subplot(133)

            # Vue cutout + marqueur pixel central
            ax0.imshow(cutout, origin='lower', cmap='gray')
            ax0.set_title("Cutout local")
            cx = px_refined - x1
            cy = py_refined - y1
            ax0.plot(cx, cy, marker='x', color='red', markersize=10, mew=2)
            ax0.plot(round(cx), round(cy), marker='+', color='cyan', markersize=10, mew=2)
            ax0.text(0.02, 0.98, f"clic=({px_click:.2f},{py_click:.2f})\nrefine=({px_refined:.2f},{py_refined:.2f})",
                     transform=ax0.transAxes, va='top', ha='left', color='yellow', fontsize=8)

            # Profils gaussiens (si fit dispo)
            if row_fit and col_fit:
                r, row, pr = row_fit
                c, col, pc = col_fit
                ax1.plot(r, row, 'k.-', label='profil Y')
                ax1.plot(r, gauss(r, *pr), 'r-', label='fit')
                ax1.set_title("Profil Y (gaussien)")
                ax1.legend(fontsize=8)

                ax2.plot(c, col, 'k.-', label='profil X')
                ax2.plot(c, gauss(c, *pc), 'r-', label='fit')
                ax2.set_title("Profil X (gaussien)")
                ax2.legend(fontsize=8)
            else:
                ax1.text(0.5, 0.5, "Fit gaussien indisponible", ha='center', va='center')
                ax2.text(0.5, 0.5, "Fit gaussien indisponible", ha='center', va='center')
                ax1.set_axis_off()
                ax2.set_axis_off()

            fig.suptitle(f"FWHM estimé: {fwhm_val:.2f} px" if fwhm_val is not None else "FWHM indisponible")
            fig.tight_layout()

            canvas = FigureCanvasTkAgg(fig, master=win)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

            btns = ttk.Frame(win)
            btns.pack(fill=tk.X, padx=6, pady=(0, 6))
            ttk.Button(
                btns,
                text="Utiliser ce centre pour T1",
                command=lambda: (self._finalize_t1_selection(px_refined, py_refined, show_message=True, context_label="diagnostic clic droit"), win.destroy())
            ).pack(side=tk.LEFT, padx=4)
            ttk.Button(btns, text="Fermer", command=win.destroy).pack(side=tk.RIGHT, padx=4)
        except Exception as e:
            logger.error(f"Diagnostic centroïde T1 impossible: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Impossible d'ouvrir le diagnostic T1:\n{e}")
            
    def _get_observation_time_range(self):
        """Retourne (t_start, t_stop) couvrant toute la série d'images."""
        if not self.image_files:
            return None

        t_min = None
        t_max = None
        for img_path in self.image_files:
            try:
                with open_fits_with_fixed_wcs(img_path) as hdul:
                    header = hdul[0].header
                date_obs = header.get('DATE-OBS', header.get('DATE'))
                if not date_obs:
                    continue
                exptime = float(header.get('EXPTIME', header.get('EXPOSURE', 0.0)))
                t_start = Time(date_obs, scale='utc')
                t_mid = t_start + exptime / 2.0 * u.second
                if t_min is None or t_mid < t_min:
                    t_min = t_mid
                if t_max is None or t_mid > t_max:
                    t_max = t_mid
            except Exception:
                continue

        if t_min is None or t_max is None:
            return None

        margin = 5 * u.minute
        return t_min - margin, t_max + margin

    def fetch_ephemeris(self):
        """Récupère les éphémérides depuis JPL Horizons avec coordonnées depuis config.py."""
        target_id = self.asteroid_id_var.get().strip()
        if not target_id:
            messagebox.showerror("Erreur", "ID cible vide")
            return
        if not self.current_image_path or not self.current_header:
            messagebox.showerror("Erreur", "Aucune image chargée")
            return

        try:
            logger.info(f"Requête JPL Horizons pour {target_id}...")

            # Lecture de la date d'observation
            date_obs = self.current_header.get('DATE-OBS', self.current_header.get('DATE'))
            if not date_obs:
                raise ValueError("Pas de DATE-OBS dans le header")

            t_obs = Time(date_obs, scale='utc')

            # Détermination de la localisation
            obs_code = self.observatory_code_var.get().strip().upper()
            if obs_code:
                location = obs_code
                logger.info(f"Utilisation du code observatoire: {obs_code}")
            else:
                # Récupération depuis config.py - JPL Horizons attend élévation en km
                from astropy.coordinates import EarthLocation
                lon = float(config.OBSERVATORY['lon'])
                lat = float(config.OBSERVATORY['lat'])
                elev_m = float(config.OBSERVATORY['elev'])  # En mètres dans config
                elev_km = elev_m / 1000.0  # Conversion en km pour JPL Horizons
                location = {'lon': lon, 'lat': lat, 'elevation': elev_km}
                logger.info(f"Utilisation des coordonnées config.py: Lon={lon}°, Lat={lat}°, Alt={elev_km}km")

            # Requête JPL Horizons - utiliser start/stop/step pour éviter les URI trop longues
            time_range = self._get_observation_time_range()
            if time_range:
                t_start, t_stop = time_range
                logger.info(
                    f"Début récupération éphémérides pour {target_id} "
                    f"(plage images: {t_start.isot} → {t_stop.isot})"
                )
            else:
                t_start = t_obs - 1.5 * u.hour
                t_stop = t_obs + 1.5 * u.hour
                logger.info(f"Début récupération éphémérides pour {target_id} autour de {t_obs.isot}")

            # Utiliser start/stop/step dès la première tentative (évite les URI trop longues)
            step_value = self.ephemeris_step_var.get().strip() or "2m"
            if step_value and step_value[-1].isdigit():
                step_value = f"{step_value}m"
                logger.info(f"Pas éphémérides normalisé: {step_value}")
            eph_table = None
            quantities_full = "1,9,20,23,24,25"  # RA/Dec, V, delta/deldot, elongation, phase, S-brt
            try:
                # Essai 1: start/stop/step avec valeurs par défaut (comme le webtool)
                logger.debug(f"Tentative 1: start/stop/step, step={step_value}")
                obj = Horizons(
                    id=target_id,
                    location=location,
                    epochs={
                        'start': t_start.isot,
                        'stop': t_stop.isot,
                        'step': step_value
                    }
                )
                eph_table = obj.ephemerides(quantities=quantities_full)
                logger.info(f"Succès avec tentative 1: {len(eph_table)} points éphémérides récupérés (pas={step_value})")
            except Exception as e:
                logger.warning(f"Échec tentative 1 (start/stop/step): {str(e)[:100]}...")
                try:
                    # Essai 2: avec quantities minimales si les valeurs par défaut échouent
                    logger.debug("Tentative 2: start/stop/step avec quantities complètes (fallback)")
                    eph_table = obj.ephemerides(quantities=quantities_full)
                    logger.info(f"Succès avec tentative 2: {len(eph_table)} points éphémérides récupérés")
                except Exception as e2:
                    logger.warning(f"Échec tentative 2: {str(e2)[:100]}...")
                    try:
                        # Essai 3: période plus courte avec step configurable par l'utilisateur
                        step_value_fallback = self.ephemeris_step_var.get().strip() or "5m"
                        if step_value_fallback and step_value_fallback[-1].isdigit():
                            step_value_fallback = f"{step_value_fallback}m"
                            logger.info(f"Pas éphémérides (fallback) normalisé: {step_value_fallback}")
                        logger.debug(f"Tentative 3: période 30min, step={step_value_fallback}")
                        obj_short = Horizons(
                            id=target_id,
                            location=location,
                            epochs={
                                'start': (t_obs - 15*u.minute).isot,
                                'stop': (t_obs + 15*u.minute).isot,
                                'step': step_value_fallback
                            }
                        )
                        eph_table = obj_short.ephemerides(quantities=quantities_full)
                        logger.info(f"Succès avec tentative 3: {len(eph_table)} points éphémérides récupérés (pas={step_value_fallback})")

                    except Exception as e3:
                        logger.warning(f"Échec tentative 3 (période courte): {str(e3)[:100]}...")
                        try:
                            # Essai 4: approche minimaliste - une seule date
                            logger.debug(f"Tentative 4: une seule date {t_obs.jd}")
                            obj_single = Horizons(id=target_id, location=location, epochs=[t_obs.jd])
                            eph_table = obj_single.ephemerides(quantities=quantities_full)
                            logger.info(f"Succès avec tentative 4 (une date): {len(eph_table)} points éphémérides récupérés")

                        except Exception as e4:
                            logger.error(f"Échec tentative 4 (une date): {str(e4)[:100]}...")
                            raise ValueError(f"Toutes les tentatives d'éphémérides ont échoué. Dernière erreur: {str(e4)[:200]}")

            if len(eph_table) == 0:
                raise ValueError("Aucune éphéméride retournée par JPL Horizons")

            # Conversion en format interne
            eph_data = Table()

            if 'datetime_jd' in eph_table.colnames:
                eph_data['datetime_jd'] = eph_table['datetime_jd']
            else:
                datetime_jd_list = []
                try:
                    for dt_str in eph_table['datetime_str']:
                        dt_str_clean = str(dt_str).strip()
                        try:
                            dt = datetime.strptime(dt_str_clean, '%Y-%b-%d %H:%M:%S.%f')
                        except ValueError:
                            dt = datetime.strptime(dt_str_clean, '%Y-%b-%d %H:%M:%S')
                        
                        t = Time(dt, scale='utc')
                        datetime_jd_list.append(t.jd)
                    
                    eph_data['datetime_jd'] = datetime_jd_list

                except Exception as e:
                    logger.error(f"Erreur parsing date: {e}")
                    raise ValueError(f"Impossible de parser la date : {e}")

            # Extraire RA et DEC avec précision maximale (en degrés)
            # Les valeurs de JPL Horizons sont des objets Angle, il faut les convertir explicitement en degrés
            if hasattr(eph_table['RA'], 'deg'):
                eph_data['RA'] = eph_table['RA'].deg  # Conversion explicite en degrés
            else:
                eph_data['RA'] = np.array([float(ra) for ra in eph_table['RA']])

            if hasattr(eph_table['DEC'], 'deg'):
                eph_data['DEC'] = eph_table['DEC'].deg  # Conversion explicite en degrés
            else:
                eph_data['DEC'] = np.array([float(dec) for dec in eph_table['DEC']])

            # Extraire la magnitude V si disponible (seulement si quantities incluait '9')
            if 'V' in eph_table.colnames:
                eph_data['V'] = eph_table['V']
            else:
                # Si pas de magnitude disponible, utiliser une valeur par défaut
                eph_data['V'] = np.full(len(eph_data), 99.0)  # Valeur par défaut
                logger.warning("Magnitude V non disponible dans les éphémérides, utilisation de valeur par défaut")
            
            # Colonnes complémentaires Horizons (si disponibles)
            def _get_col(table, *names):
                cols_lower = {c.lower(): c for c in table.colnames}
                for name in names:
                    key = name.lower()
                    if key in cols_lower:
                        return table[cols_lower[key]]
                return None
            
            col_delta = _get_col(eph_table, 'delta')
            if col_delta is not None:
                eph_data['delta'] = col_delta
            
            col_deldot = _get_col(eph_table, 'deldot')
            if col_deldot is not None:
                eph_data['deldot'] = col_deldot
            
            col_sbrt = _get_col(eph_table, 'S-brt', 'SBR', 'S_BRT', 'surf_brt')
            if col_sbrt is not None:
                eph_data['S-brt'] = col_sbrt
            
            col_elong = _get_col(eph_table, 'S-O-T', 'elong', 'ELONG')
            if col_elong is not None:
                eph_data['S-O-T'] = col_elong
            
            col_phase = _get_col(eph_table, 'S-T-O', 'alpha', 'PHASE', 'phase')
            if col_phase is not None:
                eph_data['S-T-O'] = col_phase

            self.ephemeris_data = eph_data
            logger.info(f"Éphémérides récupérées : {len(eph_data)} points")

            messagebox.showinfo("Succès", f"Éphémérides récupérées ({len(eph_data)} points).")

        except Exception as e:
            logger.error(f"Échec récupération éphémérides : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Impossible de récupérer les éphémérides :\n{e}")

    def _open_kbmod_wsl_detection_dialog(self):
        """Ouvre la boîte de dialogue de détection KBMOD via WSL (script scripts/kbmod_wsl_detect.py)."""
        if not self.directory or not self.image_files:
            logger.error("KBMOD WSL: Chargez d'abord un dossier d'images FITS.")
            return
        first_path = self.image_files[0]
        scale = 1.0
        # Priorité 1 : échelle pixel depuis config (Configuration du matériel, onglet Accueil)
        try:
            cfg_scale = config.EQUIPMENT_OBSERVATION.get("pixel_scale_arcsec")
            if cfg_scale is not None and float(cfg_scale) > 0:
                scale = float(cfg_scale)
        except (TypeError, ValueError):
            pass
        # Priorité 2 : si pas en config, déduire du WCS de la première image
        if scale <= 0 or scale == 1.0:
            try:
                with open_fits_with_fixed_wcs(first_path) as hdul:
                    from astropy.wcs import WCS
                    wcs_first = WCS(hdul[0].header)
                    try:
                        scales = proj_plane_pixel_scales(wcs_first)
                        scale = float(np.mean(scales).to(u.arcsec).value)
                        if not np.isfinite(scale) or scale <= 0:
                            scale = 1.0
                    except Exception:
                        if scale <= 0:
                            scale = 1.0
            except Exception as e:
                logger.error("KBMOD WSL: Impossible de lire le WCS de la premiere image : %s", e)
                if scale <= 0:
                    scale = 1.0

        win = tk.Toplevel(self.frame.winfo_toplevel())
        win.title("Détection KBMOD (via WSL)")
        win.geometry("720x480")
        win.transient(self.frame.winfo_toplevel())

        ttk.Label(win, text=f"Dossier : {self.directory}\n{len(self.image_files)} FITS — exécution sous WSL", font=("Helvetica", 9)).pack(pady=5)
        params_f = ttk.Frame(win)
        params_f.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(params_f, text="Vitesse min (\"/min):").grid(row=0, column=0, sticky="w", padx=2)
        v_min_var = tk.DoubleVar(value=0.02)  # bas pour inclure mouvements lents (ex. 17 Téthys ~0.03 "/min)
        ttk.Entry(params_f, textvariable=v_min_var, width=8).grid(row=0, column=1, padx=2)
        ttk.Label(params_f, text="Vitesse max (\"/min):").grid(row=0, column=2, sticky="w", padx=2)
        v_max_var = tk.DoubleVar(value=60.0)  # 60 pour astéroïdes ceinture (ex. 17 Téthys) ; NEO = plus élevé
        ttk.Entry(params_f, textvariable=v_max_var, width=8).grid(row=0, column=3, padx=2)
        ttk.Label(params_f, text="Échelle (\"/px):").grid(row=1, column=0, sticky="w", padx=2)
        scale_var = tk.DoubleVar(value=round(scale, 4))
        ttk.Entry(params_f, textvariable=scale_var, width=8).grid(row=1, column=1, padx=2)
        ttk.Label(params_f, text="Max candidats:").grid(row=1, column=2, sticky="w", padx=2)
        max_res_var = tk.IntVar(value=200)
        ttk.Entry(params_f, textvariable=max_res_var, width=8).grid(row=1, column=3, padx=2)
        ttk.Label(params_f, text="Seuil LH min:").grid(row=2, column=2, sticky="w", padx=2)
        min_lh_var = tk.DoubleVar(value=10.0)
        ttk.Entry(params_f, textvariable=min_lh_var, width=8).grid(row=2, column=3, padx=2)
        kbmod_use_gpu_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(params_f, text="Utiliser GPU (CUDA sous WSL)", variable=kbmod_use_gpu_var).grid(row=2, column=0, columnspan=2, sticky="w", padx=2, pady=2)
        kbmod_static_mask_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            params_f,
            text="Masque statique robuste",
            variable=kbmod_static_mask_var
        ).grid(row=3, column=0, columnspan=2, sticky="w", padx=2, pady=2)
        ttk.Label(params_f, text="Profil cible:").grid(row=4, column=0, sticky="w", padx=2, pady=(4, 2))
        kbmod_profiles = {
            "Astéroïde lent (ceinture)": {"vmin": 0.02, "vmax": 80.0, "min_lh": 6.0, "static_mask": False},
            "NEO rapide": {"vmin": 1.0, "vmax": 300.0, "min_lh": 12.0, "static_mask": True},
        }
        profile_var = tk.StringVar(value="Astéroïde lent (ceinture)")
        profile_combo = ttk.Combobox(
            params_f,
            textvariable=profile_var,
            values=list(kbmod_profiles.keys()),
            width=27,
            state="readonly",
        )
        profile_combo.grid(row=4, column=1, columnspan=2, sticky="w", padx=2, pady=(4, 2))

        def _apply_kbmod_profile():
            selected = profile_var.get()
            preset = kbmod_profiles.get(selected)
            if not preset:
                return
            v_min_var.set(float(preset["vmin"]))
            v_max_var.set(float(preset["vmax"]))
            min_lh_var.set(float(preset["min_lh"]))
            kbmod_static_mask_var.set(bool(preset.get("static_mask", True)))
            logger.info(
                "KBMOD WSL: profil '%s' appliqué (vmin=%.3f\"/min, vmax=%.1f\"/min, min_lh=%.1f, static_mask=%s)",
                selected,
                float(preset["vmin"]),
                float(preset["vmax"]),
                float(preset["min_lh"]),
                "on" if bool(preset.get("static_mask", True)) else "off",
            )

        ttk.Button(params_f, text="Appliquer profil", command=_apply_kbmod_profile).grid(
            row=4, column=3, sticky="w", padx=2, pady=(4, 2)
        )
        profile_combo.bind("<<ComboboxSelected>>", lambda _event: _apply_kbmod_profile())
        _apply_kbmod_profile()

        progress_var = tk.DoubleVar(value=0.0)
        progress = ttk.Progressbar(win, variable=progress_var, maximum=100)
        progress.pack(fill=tk.X, padx=10, pady=5)
        status_var = tk.StringVar(value="")
        ttk.Label(win, textvariable=status_var, font=("Helvetica", 9)).pack(pady=2)

        tree_f = ttk.Frame(win)
        tree_f.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        cols = ("RA", "Dec", "v_px/j", "likelihood")
        tree = ttk.Treeview(tree_f, columns=cols, show="headings", height=10, selectmode="browse")
        for c in cols:
            tree.heading(c, text=c)
            tree.column(c, width=120)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll = ttk.Scrollbar(tree_f, command=tree.yview)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        tree.configure(yscrollcommand=scroll.set)

        results_store = []

        def _fill_tree(res_list):
            for i in tree.get_children():
                tree.delete(i)
            for r in res_list:
                v_px_day = (float(r.get("vx_px_per_day", 0)) ** 2 + float(r.get("vy_px_per_day", 0)) ** 2) ** 0.5
                tree.insert("", tk.END, values=(
                    f"{float(r.get('ra_deg', 0)):.6f}",
                    f"{float(r.get('dec_deg', 0)):.6f}",
                    f"{v_px_day:.2f}",
                    f"{float(r.get('likelihood', 0)):.4f}",
                ))

        def do_search():
            import subprocess
            import csv as csv_module
            try:
                win.after(0, lambda: status_var.set("Lancement de KBMOD sous WSL..."))
                project_root = Path(__file__).resolve().parent.parent
                script_path = project_root / "scripts" / "kbmod_wsl_detect.py"
                if not script_path.is_file():
                    msg = f"Script introuvable : {script_path}. Voir docs/INSTALL_KBMOD_WSL.md."
                    logger.error("KBMOD WSL: %s", msg)
                    win.after(0, lambda: status_var.set("ERREUR: " + msg))
                    return
                dir_wsl = windows_path_to_wsl(self.directory)
                script_wsl = windows_path_to_wsl(script_path)
                import importlib
                importlib.reload(config)
                python_wsl = getattr(config, "KBMOD_WSL_PYTHON", "python3")
                logger.info("KBMOD WSL: Python WSL utilise = %s", python_wsl)
                cmd = [
                    "wsl", python_wsl, script_wsl, dir_wsl,
                    "--vmin", str(v_min_var.get()),
                    "--vmax", str(v_max_var.get()),
                    "--scale", str(scale_var.get()),
                    "--max-results", str(max_res_var.get()),
                    "--min-lh", str(min_lh_var.get()),
                ]
                if not kbmod_static_mask_var.get():
                    cmd.append("--no-static-mask")
                if kbmod_use_gpu_var.get():
                    cmd.append("--gpu")
                logger.info("KBMOD WSL: commande = %s", " ".join(cmd[:4]) + " ...")
                win.after(0, lambda: progress_var.set(10.0))
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=3600,
                )
                win.after(0, lambda: progress_var.set(90.0))
                if proc.returncode != 0:
                    err = (proc.stderr or "").strip()
                    out = (proc.stdout or "").strip()
                    if err:
                        logger.error("KBMOD WSL stderr: %s", err)
                    if out:
                        logger.error("KBMOD WSL stdout: %s", out)
                    msg = err or out or ""
                    if proc.returncode == 9:
                        msg = "Processus tué (code 9 = SIGKILL). Souvent manque de mémoire (OOM) sous WSL. Réduisez le nombre d'images ou augmentez la RAM WSL (.wslconfig)." + (" " + msg if msg else "")
                    elif proc.returncode == 15:
                        msg = "Processus terminé (code 15 = SIGTERM)." + (" " + msg if msg else "")
                    else:
                        msg = msg or f"Code sortie {proc.returncode}"
                    logger.error("KBMOD WSL: %s", msg)
                    if "No module named 'kbmod'" in msg:
                        hint = "Le Python WSL (%s) n'a pas kbmod. Dans WSL: conda activate astroenv (ou l'env ou vous avez installe KBMOD), puis cd ~/kbmod && pip install -e . Voir config.KBMOD_WSL_PYTHON et docs/INSTALL_KBMOD_WSL.md." % (python_wsl,)
                        logger.error("KBMOD WSL: %s", hint)
                        win.after(0, lambda: status_var.set("ERREUR: Python WSL sans kbmod. Installez KBMOD dans cet env (voir logs)."))
                    else:
                        win.after(0, lambda: status_var.set("ERREUR: " + msg[:280]))
                    return
                # Succès : logger stderr (params KBMOD, ex. plage vitesses px/jour) pour debug
                if (proc.stderr or "").strip():
                    for line in (proc.stderr or "").strip().splitlines():
                        logger.info("KBMOD WSL: %s", line)
                csv_path = Path(self.directory) / "kbmod_candidates.csv"
                if not csv_path.is_file():
                    logger.warning("KBMOD WSL: Aucun fichier kbmod_candidates.csv produit par le script WSL.")
                    win.after(0, lambda: status_var.set("ATTENTION: Aucun kbmod_candidates.csv produit."))
                    return
                res = []
                with open(csv_path, newline="", encoding="utf-8") as f:
                    reader = csv_module.DictReader(f)
                    for row in reader:
                        res.append({
                            "ra_deg": float(row.get("ra_deg", 0)),
                            "dec_deg": float(row.get("dec_deg", 0)),
                            "x0": float(row.get("x0", 0)),
                            "y0": float(row.get("y0", 0)),
                            "vx_px_per_day": float(row.get("vx_px_per_day", 0)),
                            "vy_px_per_day": float(row.get("vy_px_per_day", 0)),
                            "likelihood": float(row.get("likelihood", 0)),
                            "jd_ref": float(row.get("jd_ref", 0)),
                        })
                results_store.clear()
                results_store.extend(res)
                win.after(0, lambda: _fill_tree(res))
                win.after(0, lambda: progress_var.set(100.0))
                win.after(0, lambda: status_var.set(f"Terminé — {len(res)} candidats."))
                if len(res) == 0:
                    logger.warning("KBMOD WSL: 0 candidats. Vérifiez Vitesse max (\"/min) (ex. 60 pour ceinture), Échelle (\"/px), et voir les lignes KBMOD WSL ci-dessus dans ce log.")
            except subprocess.TimeoutExpired:
                logger.error("KBMOD WSL: Délai dépassé.")
                win.after(0, lambda: status_var.set("ERREUR: Délai dépassé."))
            except Exception as e:
                logger.error("KBMOD WSL: %s", e, exc_info=True)
                win.after(0, lambda: status_var.set("ERREUR: " + str(e)[:200]))

        def run_in_thread():
            progress_var.set(0.0)
            status_var.set("")
            import threading
            t = threading.Thread(target=do_search, daemon=True)
            t.start()

        def select_as_t1():
            sel = tree.selection()
            if not sel or not results_store:
                logger.warning("KBMOD WSL: Sélectionnez un candidat dans la liste.")
                status_var.set("ATTENTION: Sélectionnez un candidat dans la liste.")
                return
            idx = tree.index(sel[0])
            if idx < 0 or idx >= len(results_store):
                return
            r = results_store[idx]
            ra_deg = r.get("ra_deg")
            dec_deg = r.get("dec_deg")
            x0 = r.get("x0", 0)
            y0 = r.get("y0", 0)
            vx = r.get("vx_px_per_day", 0)
            vy = r.get("vy_px_per_day", 0)
            jd_ref = r.get("jd_ref", 0)
            if ra_deg is None or dec_deg is None or not np.isfinite(ra_deg) or not np.isfinite(dec_deg):
                logger.error("KBMOD WSL: Coordonnées invalides pour ce candidat.")
                status_var.set("ERREUR: Coordonnées invalides pour ce candidat.")
                return
            target_sky = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
            self.target_t1_sky = target_sky
            self.target_t1_px = (float(x0), float(y0))
            # Ancres manuelles T1 (première et dernière image) pour le batch
            try:
                with open_fits_with_fixed_wcs(first_path) as h1:
                    jd_first = float(h1[0].header.get("JD-UTC", 0) or h1[0].header.get("JD", 0))
                    if jd_first == 0 and h1[0].header.get("DATE-OBS"):
                        t = Time(h1[0].header["DATE-OBS"], scale="utc")
                        jd_first = t.jd
                    wcs_first_img = WCS(h1[0].header)
                last_path = self.image_files[-1]
                with open_fits_with_fixed_wcs(last_path) as h2:
                    jd_last = float(h2[0].header.get("JD-UTC", 0) or h2[0].header.get("JD", 0))
                    if jd_last == 0 and h2[0].header.get("DATE-OBS"):
                        t = Time(h2[0].header["DATE-OBS"], scale="utc")
                        jd_last = t.jd
                    wcs_last_img = WCS(h2[0].header)
                dt_days = jd_last - jd_ref
                x_last = x0 + vx * dt_days
                y_last = y0 + vy * dt_days
                sky_first = wcs_first_img.pixel_to_world(x0, y0)
                sky_last = wcs_last_img.pixel_to_world(x_last, y_last)
                self.manual_t1_anchor_first = {"jd": jd_first, "ra_deg": float(sky_first.ra.deg), "dec_deg": float(sky_first.dec.deg)}
                self.manual_t1_anchor_last = {"jd": jd_last, "ra_deg": float(sky_last.ra.deg), "dec_deg": float(sky_last.dec.deg)}
            except Exception as e:
                logger.warning(f"Ancres T1 KBMOD: {e}")
            # Sélection photométrique T1 avec ouvertures par défaut
            r_ap = getattr(self, "aperture_radius_var", None)
            r_ap = r_ap.get() if r_ap else 8.0
            r_in = getattr(self, "annulus_inner_var", None)
            r_in = r_in.get() if r_in else 12.0
            r_out = getattr(self, "annulus_outer_var", None)
            r_out = r_out.get() if r_out else 20.0
            self.current_selections = [
                {"label": "T1", "coord": target_sky, "r_ap": r_ap, "r_in": r_in, "r_out": r_out}
            ]
            self.refresh_display()
            logger.info("KBMOD WSL: Candidat KBMOD défini comme T1. Ajoutez des comparateurs (SET-UP PHOTOMÉTRIE) puis lancez le batch.")
            win.destroy()

        btn_f = ttk.Frame(win)
        btn_f.pack(fill=tk.X, padx=10, pady=10)
        ttk.Button(btn_f, text="Lancer la détection", command=run_in_thread).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_f, text="Sélectionner comme T1", command=select_as_t1).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_f, text="Fermer", command=win.destroy).pack(side=tk.LEFT, padx=5)

    def run_photometry_setup(self):
        """Ouverture de la fenêtre de sélection photométrique avec T1 et étoiles de comparaison."""
        if not self.current_image_path or self.current_data is None:
            messagebox.showerror("Erreur", "Aucune image chargée")
            return

        if not self.wcs:
            messagebox.showerror("Erreur", "WCS manquant. Lancez l'astrométrie d'abord.")
            return

        self._launch_manual_comps_selection()

    def _launch_manual_comps_selection(self):
        """Fenêtre simple de sélection manuelle des étoiles de comparaison (astéroïdes)."""
        from core.photometry_pipeline_asteroids import launch_photometry_aperture, refine_centroid
        from astropy.visualization import ZScaleInterval

        data = self.current_data
        wcs = self.wcs
        target_coord = self.target_t1_sky
        t1_patch = None
        t1_text = None

        if target_coord is None:
            messagebox.showerror("Erreur", "T1 non défini. Cliquez d'abord sur T1 dans l'image principale.")
            return

        root = tk.Toplevel()
        root.title("Sélection manuelle des comparateurs")
        root.geometry("1350x800")

        # ZScale
        try:
            interval = ZScaleInterval()
            vmin, vmax = interval.get_limits(data)
        except Exception:
            vmin, vmax = np.percentile(data, [10, 90])

        # --- Zone Gauche (Liste) ---
        left_frame = tk.Frame(root, width=420, bg="#f0f0f0")
        left_frame.pack(side=tk.LEFT, fill=tk.Y)
        left_frame.pack_propagate(False)

        tk.Label(left_frame, text=" Étoiles de comparaison ", bg="#ddd", font=("Arial", 10, "bold"), pady=5).pack(fill=tk.X)

        list_container = tk.Frame(left_frame, bg="#f0f0f0")
        list_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=5)

        btn_frame = tk.Frame(left_frame, bg="#e0e0e0", pady=10)
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # --- Zone Droite (Graphique) ---
        right_frame = tk.Frame(root, bg="black")
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        status_label = tk.Label(
            right_frame,
            text="Cliquez sur les étoiles de comparaison (gauche).",
            bg="black",
            fg="#00ff00",
            font=("Consolas", 11, "bold"),
            pady=5
        )
        status_label.pack(side=tk.TOP, fill=tk.X)

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': wcs})
        ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)

        # Afficher T1 si déjà défini
        if target_coord is not None:
            try:
                tx, ty = wcs.world_to_pixel(target_coord)
                if np.isfinite(tx) and np.isfinite(ty):
                    t1_patch = Circle((tx, ty), 18, edgecolor='red', lw=2, fill=False)
                    t1_text = ax.annotate("T1", (tx+20, ty+20), color='red', fontweight='bold')
                    ax.add_patch(t1_patch)
            except Exception as e:
                logger.debug(f"Affichage T1 échoué: {e}")

        canvas = FigureCanvasTkAgg(fig, master=right_frame)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, right_frame)
        toolbar.update()
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Changer le curseur pour le rendre plus visible en mode zoom
        def update_cursor(event):
            mode = toolbar.mode if hasattr(toolbar, 'mode') else ""
            if mode and mode != "" and 'zoom' in mode.lower():
                canvas_widget.config(cursor="pencil")
            else:
                canvas_widget.config(cursor="")

        canvas.mpl_connect('motion_notify_event', update_cursor)

        comps = []
        comp_patches = []
        comp_texts = []

        def update_list():
            for w in list_container.winfo_children():
                w.destroy()
            for i, comp in enumerate(comps):
                row = tk.Frame(list_container, bg="#f0f0f0")
                row.pack(fill="x", pady=2)
                tk.Label(row, text=f"C{i+1}", width=4, fg="blue", bg="#f0f0f0").pack(side="left")
                tk.Label(row, text=f"{comp.ra.deg:.5f}", width=12, bg="#fff").pack(side="left", padx=2)
                tk.Label(row, text=f"{comp.dec.deg:.5f}", width=12, bg="#fff").pack(side="left", padx=2)

        def add_comp(x, y):
            try:
                x_c, y_c = refine_centroid(data, x, y, box_size=15)
                if not (np.isfinite(x_c) and np.isfinite(y_c)):
                    return

                coord = wcs.pixel_to_world(x_c, y_c)

                # Éviter T1 et doublons
                if target_coord is not None and coord.separation(target_coord).arcsec < 5.0:
                    status_label.config(text="⚠️ Cliquez sur une étoile différente de T1", fg="white", bg="red")
                    return
                for c in comps:
                    if coord.separation(c).arcsec < 5.0:
                        status_label.config(text="⚠️ Comparateur déjà sélectionné", fg="white", bg="red")
                        return

                comps.append(coord)
                idx = len(comps)

                patch = Rectangle((x_c-15, y_c-15), 30, 30, edgecolor='cyan', fill=False)
                text = ax.annotate(f"C{idx}", (x_c+18, y_c), color='cyan', fontsize=9)
                ax.add_patch(patch)
                comp_patches.append(patch)
                comp_texts.append(text)
                canvas.draw()
                update_list()
                status_label.config(text=f"Comparateur ajouté: C{idx}", fg="#00ff00", bg="black")
            except Exception as e:
                logger.error(f"Erreur ajout comparateur: {e}", exc_info=True)
                status_label.config(text=f"Erreur ajout comparateur: {e}", fg="white", bg="red")

        def on_click(event):
            if event.button != 1:
                return
            if toolbar.mode != "" and toolbar.mode is not None:
                status_label.config(text=f"⚠️ MODE {str(toolbar.mode).upper()} ACTIF ! Désactivez la loupe.", fg="white", bg="red")
                return
            if not event.inaxes:
                return
            x_raw, y_raw = event.xdata, event.ydata
            if not np.isfinite(x_raw) or not np.isfinite(y_raw):
                return
            add_comp(x_raw, y_raw)

        fig.canvas.mpl_connect("button_press_event", on_click)

        def remove_last():
            if not comps:
                return
            comps.pop()
            if comp_patches:
                try:
                    comp_patches.pop().remove()
                except Exception:
                    pass
            if comp_texts:
                try:
                    comp_texts.pop().remove()
                except Exception:
                    pass
            canvas.draw()
            update_list()

        def clear_all():
            comps.clear()
            for p in comp_patches:
                try:
                    p.remove()
                except Exception:
                    pass
            for t in comp_texts:
                try:
                    t.remove()
                except Exception:
                    pass
            comp_patches.clear()
            comp_texts.clear()
            canvas.draw()
            update_list()

        def validate():
            if not comps:
                messagebox.showwarning("Stop", "Sélectionnez au moins un comparateur.")
                return

            target_data = {'coord': target_coord}
            launch_photometry_aperture(
                fits_path=self.current_image_path,
                target_data=target_data,
                comp_coords_data=[{'coord': c} for c in comps],
                on_finish=self._on_aperture_finished
            )
            root.destroy()

        tk.Button(btn_frame, text="↩ Supprimer dernier", command=remove_last).pack(fill=tk.X, padx=5, pady=2)
        tk.Button(btn_frame, text="🧹 Effacer tout", command=clear_all).pack(fill=tk.X, padx=5, pady=2)
        tk.Button(btn_frame, text="✅ VALIDER LA SÉLECTION", command=validate).pack(fill=tk.X, padx=5, pady=5)

    def _on_aperture_finished(self, selections):
        """Callback appelé quand les apertures sont validées."""
        self.current_selections = selections
        logger.info(f"Photométrie configurée : {len(selections)} étoiles avec apertures")
        messagebox.showinfo("Succès", f"Photométrie configurée : {len(selections)} étoiles sélectionnées")
    
    def launch_selection_window(self):
        """Fenêtre de sélection photométrique avec T1 et étoiles de comparaison."""
        from core.photometry_pipeline_asteroids import launch_photometry_aperture
        
        # Vérifier que T1 a été sélectionné par clic
        if self.target_t1_sky is None:
            messagebox.showerror("Erreur", "T1 non sélectionné. Cliquez sur l'image pour désigner T1.")
            return
        
        # Utiliser les coordonnées T1 sélectionnées par l'utilisateur
        target_coord = self.target_t1_sky
        logger.info(f"Utilisation de T1 sélectionné : RA={target_coord.ra.deg:.6f}°, Dec={target_coord.dec.deg:.6f}°")
        
        # Initialisation de comp_coords
        comp_coords = []
        nx, ny = self.current_data.shape[1], self.current_data.shape[0]

        # Marge de sécurité pour éviter les bords
        margin = 50  # Marge importante pour s'assurer que les étoiles sont bien visibles

        # Recherche d'étoiles de comparaison : autour de T1 avec un rayon de recherche configurable
        try:
            # Recherche autour de T1 plutôt que sur tout le FOV
            # Rayon de recherche en arcmin (par défaut 10 arcmin autour de T1)
            search_radius_arcmin = 10.0  # Rayon de recherche autour de T1 (en arcmin)
            search_radius = search_radius_arcmin * u.arcmin
            
            logger.info(f"Recherche de comparateurs autour de T1 (rayon: {search_radius_arcmin} arcmin)")

            # Requête Gaia avec cache autour de T1
            gaia_table = optimized_query_gaia(
                target_coord, search_radius,
                mag_limit=self.gaia_mag_limit_var.get(),
                gaia_cache=self.gaia_cache
            )

            if len(gaia_table) == 0:
                messagebox.showwarning("Avertissement", "Aucune étoile Gaia trouvée pour comparaison")
            else:
                # Détection directe des étoiles dans l'image pour s'assurer qu'elles sont visibles
                user_fwhm = self.fwhm_var.get()
                threshold_sigma = self.threshold_sigma_var.get()
                detected_sources = _detect_sources_for_astrometry(
                    self.current_data, user_fwhm, threshold_sigma, max_sources=500
                )
                
                if len(detected_sources) == 0:
                    logger.warning("Aucune source détectée dans l'image pour comparaison")
                else:
                    # Conversion des sources détectées en coordonnées célestes
                    detected_sky = self.wcs.pixel_to_world(detected_sources['xcentroid'], detected_sources['ycentroid'])
                    
                    # Matching avec Gaia pour obtenir les magnitudes
                    gaia_sky = SkyCoord(ra=gaia_table['RA_ICRS'], dec=gaia_table['DE_ICRS'], unit=(u.deg, u.deg))
                    idx, d2d, _ = match_coordinates_sky(detected_sky, gaia_sky)
                    match_mask = d2d < 3.0 * u.arcsec  # Tolérance de 3" pour le matching
                    
                    # Magnitude cible
                    t_img = Time(self.current_header.get('DATE-OBS', self.current_header.get('DATE')), scale='utc')
                    target_mag = _interpolate_magnitude_V(self.ephemeris_data, t_img)
                    
                    # Filtrage des comparateurs candidats
                    for i in range(len(detected_sky)):
                        if not match_mask[i]:
                            continue  # Pas de match Gaia
                        
                        # Vérifier que la source est dans l'image (pas sur les bords)
                        px = float(detected_sources['xcentroid'][i])
                        py = float(detected_sources['ycentroid'][i])
                        
                        if not (margin <= px < nx - margin and margin <= py < ny - margin):
                            continue  # Trop proche des bords
                        
                        # Correspondance Gaia
                        gaia_idx = int(idx[i])
                        if gaia_idx < 0 or gaia_idx >= len(gaia_table):
                            continue  # Index invalide
                        
                        comp_mag = gaia_table['Gmag'][gaia_idx] if 'Gmag' in gaia_table.colnames else None
                        
                        if comp_mag is None or not np.isfinite(comp_mag):
                            continue
                        
                        # Exclure les étoiles variables
                        if 'phot_variable_flag' in gaia_table.colnames:
                            var_flag = gaia_table['phot_variable_flag'][gaia_idx]
                            # Le flag peut être 'VARIABLE', True, ou une chaîne indiquant la variabilité
                            # Exclure si variable (différentes représentations possibles)
                            if var_flag is not None and var_flag != '':
                                var_str = str(var_flag).strip().upper()
                                if var_str in ['VARIABLE', 'TRUE', '1', 'Y', 'YES']:
                                    logger.debug(f"Comparateur exclu (variable) : px={px:.1f}, py={py:.1f}, mag={comp_mag:.2f}, flag={var_flag}")
                                    continue
                        
                        # Coordonnées célestes de la source détectée
                        src_sky = detected_sky[i]
                        
                        # Exclure T1 (séparation < 5")
                        sep_from_t1 = target_coord.separation(src_sky).arcsec
                        if sep_from_t1 < 5.0:
                            continue  # C'est probablement T1
                        
                        # Chercher les comparateurs dans une zone raisonnable autour de T1 (5-300")
                        # avec une magnitude proche de celle de T1 (±2.5 mag)
                        if 5.0 < sep_from_t1 < 300.0 and abs(comp_mag - target_mag) <= 2.5:
                            comp_coords.append(src_sky)
                            logger.debug(f"Comparateur ajouté : px={px:.1f}, py={py:.1f}, mag={comp_mag:.2f}, sep={sep_from_t1:.1f}\"")
                            if len(comp_coords) >= 10:  # Limite à 10 comparateurs
                                break

                logger.info(f"{len(comp_coords)} étoiles de comparaison trouvées (détectées dans l'image)")

        except Exception as e:
            logger.error(f"Erreur recherche comparateurs : {e}")
            comp_coords = []
        
        # Lancement de la fenêtre de sélection
        def on_finish(selections):
            """Callback appelé quand la sélection est validée."""
            self.current_selections = selections
            logger.info(f"Sélection validée : {len(selections)} étoiles")
            messagebox.showinfo("Succès", f"Photométrie configurée : {len(selections)} étoiles sélectionnées")
        
        launch_photometry_aperture(
            fits_path=self.current_image_path,
            target_data={'coord': target_coord},
            comp_coords_data=[{'coord': c} for c in comp_coords],
            on_finish=on_finish
        ) 
        

    def solve_astrometry_zero_aperture(self, path=None, skip_zero_aperture=False, extrapolation_params=None):
        """Résout l'astrométrie avec optimisations (cache Gaia, KD-tree, apertures intelligentes)."""
        if path is None:
            path = self.current_image_path
        
        if not path or not os.path.exists(path):
            if path == self.current_image_path:
                messagebox.showerror("Erreur", "Aucune image sélectionnée")
            return
            
        try:
            logger.info(f"=== Astrométrie optimisée pour {Path(path).name} ===")
            
            # Chargement de l'image (header corrigé SIP/CTYPE à la source)
            with open_fits_with_fixed_wcs(path) as hdul:
                data = hdul[0].data.astype(float)
                header = hdul[0].header.copy()
            
            # Paramètres
            user_fwhm = self.fwhm_var.get()
            threshold_sigma = self.threshold_sigma_var.get()
            max_sources = self.max_sources_var.get()
            match_radius = self.match_radius_var.get() * u.arcsec
            mag_limit = self.gaia_mag_limit_var.get()
            
            # WCS initial depuis header
            try:
                wcs_init = WCS(header)
                if not wcs_init.is_celestial:
                    raise ValueError("WCS initial invalide")
            except Exception:
                if path == self.current_image_path:
                    logger.warning("WCS initial requis dans le header (astrométrie FWHM).")
                return

            # Tentative avec le solveur optimisé
            try:
                logger.debug(f"skip_zero_aperture={skip_zero_aperture} pour {Path(path).name}")
                fast_solver = FastAstrometrySolver(
                    gaia_cache=self.gaia_cache,
                    fwhm=user_fwhm,
                    threshold_sigma=threshold_sigma,
                    max_sources=max_sources,
                    match_radius=match_radius,
                    mag_limit=mag_limit,
                    use_gpu=self.use_gpu_var.get() and HAS_GPU  # Option GPU depuis l'interface
                )
                logger.debug(f"Appel de fast_solver.solve avec skip_zero_aperture={skip_zero_aperture}")
                wcs_fitted, zero_rms, stats = fast_solver.solve(data, header, skip_zero_aperture=skip_zero_aperture, extrapolation_params=extrapolation_params)
                logger.debug(f"fast_solver.solve terminé : zero_rms={zero_rms:.4f}, method={stats.get('method', 'unknown')}")
                
                # Mise à jour du header
                wcs_header = wcs_fitted.to_header()
                header.update(wcs_header)
                
                # Stocker toutes les métadonnées astrométriques enrichies dans le header
                header['ASTREF'] = ('Gaia DR3', 'Astrometric reference catalog')
                header['ASTNREF'] = (stats['n_stars'], 'Number of reference stars used')
                header['ASTRRMS'] = (stats['rms_total'], 'Astrometry RMS total (arcsec)')
                header['ASTRRMSR'] = (stats['rms_ra'], 'Astrometry RMS RA (arcsec)')
                header['ASTRRMSD'] = (stats['rms_dec'], 'Astrometry RMS Dec (arcsec)')
                header['ASTRMED'] = (stats['median_residual'], 'Median residual (arcsec)')
                header['ASTNOUT'] = (stats['n_outliers'], 'Number of outliers (>3-sigma)')
                header['ASTMETHOD'] = (stats['method'], 'Method used (classical/zero-aperture)')
                header['ASTRMSC'] = (stats['rms_classical'], 'RMS classical method (arcsec)')
                header['ASTRMSZ'] = (stats['rms_zero_aperture'], 'RMS zero-aperture method (arcsec)')
                header['ASTRCOR'] = (stats['rms_corr'], 'RMS correlation coefficient RA/Dec (for ADES rmsCorr)')
                if stats.get('fit_r2') is not None:
                    header['ASTR2'] = (stats['fit_r2'], 'R2 of extrapolation fit')
                
                # Sauvegarde
                with open_fits_with_fixed_wcs(path, mode='update') as hdul:
                    hdul[0].header = header
                    hdul.flush()
                
                # Log complet avec tous les détails RMS
                logger.info(f"=== Astrométrie zéro-aperture réussie pour {Path(path).name} ===")
                logger.info(f"Méthode utilisée: {stats['method']}")
                logger.info(f"RMS total: {stats['rms_total']:.4f}\"")
                logger.info(f"RMS RA: {stats['rms_ra']:.4f}\"")
                logger.info(f"RMS Dec: {stats['rms_dec']:.4f}\"")
                logger.info(f"RMS zéro-aperture: {stats['rms_zero_aperture']:.4f}\"")
                logger.info(f"RMS classique: {stats['rms_classical']:.4f}\"")
                logger.info(f"RMS final: {zero_rms:.4f}\"")
                logger.info(f"Résidu médian: {stats['median_residual']:.4f}\"")
                logger.info(f"Coefficient de corrélation RA/Dec: {stats['rms_corr']:.4f}")
                if stats.get('fit_r2') is not None:
                    logger.info(f"R² de l'extrapolation: {stats['fit_r2']:.4f}")
                logger.info(f"Nombre d'étoiles de référence: {stats['n_stars']}")
                logger.info(f"Nombre d'outliers (>3σ): {stats['n_outliers']}")
                logger.info(f"Catalogue de référence: Gaia DR3")

                # Mémoriser la série zéro-ouverture (6 ouvertures) pour affichage tableau
                try:
                    image_key = str(Path(path).resolve())
                    za_series = stats.get('zero_aperture_series', []) if isinstance(stats, dict) else []
                    if isinstance(za_series, list) and len(za_series) > 0:
                        self.zero_aperture_series_by_image[image_key] = za_series
                except Exception as e:
                    logger.debug(f"Stockage série zéro-ouverture impossible pour {Path(path).name}: {e}")

                # Mise à jour du WCS dans l'interface si c'est l'image courante
                if path == self.current_image_path:
                    self.wcs = wcs_fitted
                    self.current_header = header
                    # Ne pas repositionner T1 automatiquement, il reste à sa position si déjà défini
                    self.refresh_display()

                # Pour les images batch, mettre à jour T1 temporairement pour le log
                if self.ephemeris_data:
                    try:
                        date_obs_str = header.get('DATE-OBS', header.get('DATE'))
                        if date_obs_str:
                            exptime = float(header.get('EXPTIME', header.get('EXPOSURE', 0.0)))
                            t_start = Time(date_obs_str, scale='utc')
                            t_img = t_start + exptime / 2.0 * u.second
                            target_sky = _interpolate_target_position(self.ephemeris_data, t_img.jd)
                            if target_sky:
                                tx, ty = wcs_fitted.world_to_pixel(target_sky)
                                ny, nx = data.shape
                                if 0 <= tx < nx and 0 <= ty < ny:
                                    logger.debug(f"T1 pour {Path(path).name} : [{tx:.1f}, {ty:.1f}]")
                    except Exception as e:
                        logger.debug(f"Impossible de calculer T1 pour {Path(path).name} : {e}")
                return stats

            except Exception as e:
                logger.error(f"Erreur astrométrie optimisée : {e}", exc_info=True)
                # Message d'erreur uniquement pour l'image courante (pas en batch)
                if path == self.current_image_path:
                    messagebox.showerror("Erreur", f"Erreur lors de l'astrométrie :\n{e}")
                return None
        
        except Exception as e:
            logger.error(f"Erreur astrométrie : {e}", exc_info=True)
            # Message d'erreur uniquement pour l'image courante (pas en batch)
            if path == self.current_image_path:
                messagebox.showerror("Erreur", f"Erreur lors de l'astrométrie :\n{e}")
            return None
    
    def solve_astrometry_classical(self, path=None):
        """
        Résout l'astrométrie classique avec une ouverture unique calculée à partir du FWHM estimé.
        Utilise les ouvertures proposées depuis le FWHM (r_ap = 1.4 * FWHM).
        """
        if path is None:
            path = self.current_image_path
        
        if not path or not os.path.exists(path):
            if path == self.current_image_path:
                messagebox.showerror("Erreur", "Aucune image sélectionnée")
            return

        try:
            logger.info(f"=== Astrométrie classique pour {Path(path).name} ===")
            
            # Chargement de l'image (header corrigé SIP/CTYPE à la source)
            with open_fits_with_fixed_wcs(path) as hdul:
                data = hdul[0].data.astype(float)
                header = hdul[0].header.copy()
            
            # Paramètres
            threshold_sigma = self.threshold_sigma_var.get()
            max_sources = self.max_sources_var.get()
            match_radius = self.match_radius_var.get() * u.arcsec
            mag_limit = self.gaia_mag_limit_var.get()
            
            # WCS initial depuis header
            try:
                wcs_init = WCS(header)
                if not wcs_init.is_celestial:
                    raise ValueError("WCS initial invalide")
            except Exception:
                if path == self.current_image_path:
                    messagebox.showerror("Erreur", "WCS initial requis dans le header")
                return

            # Estimation du FWHM : priorité à T1 si sélectionné, sinon via éphémérides
            fwhm_estimated = None
            t1_px_for_fwhm = None
            if path == self.current_image_path and self.target_t1_px:
                t1_px_for_fwhm = self.target_t1_px
            elif self.ephemeris_data:
                try:
                    date_obs_str = header.get('DATE-OBS', header.get('DATE'))
                    if date_obs_str:
                        exptime = float(header.get('EXPTIME', header.get('EXPOSURE', 0.0)))
                        t_start = Time(date_obs_str, scale='utc')
                        t_img = t_start + exptime / 2.0 * u.second
                        target_sky = _interpolate_target_position(self.ephemeris_data, t_img.jd)
                        if target_sky:
                            tx, ty = wcs_init.world_to_pixel(target_sky)
                            ny, nx = data.shape
                            if 0 <= tx < nx and 0 <= ty < ny:
                                t1_px_for_fwhm = (tx, ty)
                                logger.info(f"T1 éphémérides pour FWHM : [{tx:.1f}, {ty:.1f}]")
                except Exception as e:
                    logger.debug(f"Impossible d'estimer T1 via éphémérides : {e}")

            if t1_px_for_fwhm is None and self.last_t1_px_for_fwhm:
                t1_px_for_fwhm = self.last_t1_px_for_fwhm
                logger.info(
                    f"T1 fallback (image précédente) pour FWHM : "
                    f"[{t1_px_for_fwhm[0]:.1f}, {t1_px_for_fwhm[1]:.1f}]"
                )

            if t1_px_for_fwhm:
                try:
                    tx, ty = t1_px_for_fwhm
                    tx_ref, ty_ref = refine_centroid(data, tx, ty, box_size=11)
                    fwhm_result = estimate_fwhm_marginal(data, tx_ref, ty_ref, box_size=25)
                    if fwhm_result and fwhm_result[0] and np.isfinite(fwhm_result[0]):
                        fwhm_estimated = fwhm_result[0]
                        self.last_t1_px_for_fwhm = (tx_ref, ty_ref)
                        logger.info(f"FWHM estimé depuis T1 : {fwhm_estimated:.2f} pixels")
                    elif self.last_t1_px_for_fwhm:
                        tx, ty = self.last_t1_px_for_fwhm
                        tx_ref, ty_ref = refine_centroid(data, tx, ty, box_size=11)
                        fwhm_result = estimate_fwhm_marginal(data, tx_ref, ty_ref, box_size=25)
                        if fwhm_result and fwhm_result[0] and np.isfinite(fwhm_result[0]):
                            fwhm_estimated = fwhm_result[0]
                            self.last_t1_px_for_fwhm = (tx_ref, ty_ref)
                            logger.info(f"FWHM estimé depuis T1 (fallback) : {fwhm_estimated:.2f} pixels")
                except Exception as e:
                    logger.debug(f"Impossible d'estimer FWHM depuis T1 : {e}")

            # Si pas de FWHM depuis T1, utiliser la valeur par défaut
            if fwhm_estimated is None or not np.isfinite(fwhm_estimated):
                fwhm_estimated = self.fwhm_var.get()
                logger.info(f"Utilisation FWHM par défaut : {fwhm_estimated:.2f} pixels")

            # Calcul de l'ouverture à utiliser (r_ap = 1.4 * FWHM)
            r_ap, r_in, r_out = calculate_aperture_radii_from_fwhm(fwhm_estimated, default_fwhm=4.0)
            logger.info(f"Ouverture classique calculée : {r_ap:.1f} px (FWHM={fwhm_estimated:.2f} px)")
            
            # 1. Calcul centre et rayon FOV
            ny, nx = data.shape
            px = np.array([0, nx, 0, nx])
            py = np.array([0, 0, ny, ny])
            corners_sky = wcs_init.pixel_to_world(px, py)
            center_ra = np.mean(corners_sky.ra.deg)
            center_dec = np.mean(corners_sky.dec.deg)
            center_coord = SkyCoord(ra=center_ra*u.deg, dec=center_dec*u.deg, frame='icrs')
            radius = center_coord.separation(corners_sky).max() * 1.1
            
            # 2. Requête Gaia avec cache
            gaia_table = optimized_query_gaia(
                center_coord, radius,
                mag_limit=mag_limit,
                gaia_cache=self.gaia_cache
            )
            
            if len(gaia_table) == 0:
                if path == self.current_image_path:
                    messagebox.showerror("Erreur", "Aucune étoile Gaia trouvée")
                return
            
            logger.info(f"Gaia : {len(gaia_table)} étoiles récupérées")
            
            # 3. Détection des sources dans l'image
            _, _, std = sigma_clipped_stats(data)
            daofind = DAOStarFinder(fwhm=fwhm_estimated, threshold=threshold_sigma*std)
            sources = daofind(data)
            
            if not sources or len(sources) == 0:
                if path == self.current_image_path:
                    messagebox.showerror("Erreur", "Aucune source détectée dans l'image")
                return
            
            sources.sort('flux', reverse=True)
            sources = sources[:max_sources]
            logger.info(f"Détection : {len(sources)} sources")
            
            # 4. Matching avec Gaia (coordonnées brutes)
            detected_sky = wcs_init.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
            gaia_sky = SkyCoord(ra=gaia_table['RA_ICRS'], dec=gaia_table['DE_ICRS'], unit=u.deg)
            
            idx, d2d, _ = match_coordinates_sky(detected_sky, gaia_sky)
            valid_mask = d2d < match_radius
            
            if np.sum(valid_mask) < 3:
                if path == self.current_image_path:
                    messagebox.showerror("Erreur", f"Pas assez de correspondances (< 3) entre image et Gaia")
                return
            
            matches_count = np.sum(valid_mask)
            logger.info(f"Matching : {matches_count} correspondances trouvées")
            
            # 5. Affinage des centroïdes avec l'ouverture classique (r_ap)
            box_size = int(2 * r_ap) + 1
            if box_size % 2 == 0:
                box_size += 1
            
            refined_x = []
            refined_y = []
            
            for src in sources[valid_mask]:
                x_init = src['xcentroid']
                y_init = src['ycentroid']
                
                try:
                    x_int, y_int = int(round(x_init)), int(round(y_init))
                    r = box_size // 2
                    
                    if r <= x_int < nx - r and r <= y_int < ny - r:
                        cutout = data[y_int - r:y_int + r + 1, x_int - r:x_int + r + 1]
                        if cutout.size >= 6 and np.isfinite(cutout).any():
                            xc, yc = centroid_2dg(cutout)
                            refined_x.append(x_int - r + xc)
                            refined_y.append(y_int - r + yc)
                        else:
                            refined_x.append(x_init)
                            refined_y.append(y_init)
                    else:
                        refined_x.append(x_init)
                        refined_y.append(y_init)
                except Exception as e:
                    logger.debug(f"Échec affinage centroïde : {e}")
                    refined_x.append(x_init)
                    refined_y.append(y_init)
            
            # 6. Fit WCS avec centroïdes affinés
            xy_pixels = np.array([refined_x, refined_y]).T
            world_coords = SkyCoord(
                ra=gaia_table['RA_ICRS'][idx[valid_mask]], 
                dec=gaia_table['DE_ICRS'][idx[valid_mask]], 
                unit=u.deg
            )
            
            wcs_fitted = fit_wcs_from_points(
                xy_pixels.T, world_coords, proj_point='center', sip_degree=None
            )
            
            # 7. Calcul des statistiques astrométriques détaillées
            stats = calculate_astrometric_statistics(
                wcs_fitted,
                xy_pixels,
                world_coords,
                outlier_threshold=3.0
            )
            rms = stats['rms_total']
            
            logger.info(f"Astrométrie classique réussie : RMS={rms:.4f}\" (RA={stats['rms_ra']:.4f}\", Dec={stats['rms_dec']:.4f}\", corr={stats['rms_corr']:.3f}, ouverture {r_ap:.1f} px)")
            
            # Mise à jour du header
            wcs_header = wcs_fitted.to_header()
            header.update(wcs_header)
            
            # Stocker toutes les métadonnées astrométriques enrichies dans le header
            header['ASTREF'] = ('Gaia DR3', 'Astrometric reference catalog')
            header['ASTNREF'] = (stats['n_stars'], 'Number of reference stars used')
            header['ASTRRMS'] = (stats['rms_total'], 'Astrometry RMS total (arcsec) - classical method')
            header['ASTRRMSR'] = (stats['rms_ra'], 'Astrometry RMS RA (arcsec) - classical method')
            header['ASTRRMSD'] = (stats['rms_dec'], 'Astrometry RMS Dec (arcsec) - classical method')
            header['ASTRMED'] = (stats['median_residual'], 'Median residual (arcsec)')
            header['ASTNOUT'] = (stats['n_outliers'], 'Number of outliers (>3-sigma)')
            header['ASTMETHOD'] = ('classical', 'Method used')
            header['ASTRCOR'] = (stats['rms_corr'], 'RMS correlation coefficient RA/Dec (for ADES rmsCorr)')
            
            # Sauvegarde
            with open_fits_with_fixed_wcs(path, mode='update') as hdul:
                hdul[0].header = header
                hdul.flush()
            
            # Mise à jour du WCS dans l'interface si c'est l'image courante
            if path == self.current_image_path:
                self.wcs = wcs_fitted
                self.current_header = header
                self.refresh_display()
            
        except Exception as e:
            logger.error(f"Erreur astrométrie classique : {e}", exc_info=True)
            if path == self.current_image_path:
                messagebox.showerror("Erreur", f"Erreur lors de l'astrométrie classique :\n{e}")
    
    def run_single_astrometry_with_method(self):
        """Lance l'astrométrie sur l'image courante selon la méthode sélectionnée."""
        method = self.astrometry_method_var.get()
        self._start_process_log_capture("astrométrie (image courante)")
        try:
            if method == "classical":
                self.solve_astrometry_classical()
            else:  # "zero_aperture" par défaut
                self.solve_astrometry_zero_aperture()
        finally:
            self._stop_process_log_capture("astrométrie (image courante)")
    
    def run_single_astrometry(self):
        """Lance l'astrométrie sur l'image courante (ancien nom, conservé pour compatibilité)."""
        self.run_single_astrometry_with_method()
    
    def _batch_astrometry_worker(self):
        """Lance l'astrométrie en batch sur toutes les images (threadé)."""
        if not self.image_files:
            messagebox.showerror("Erreur", "Aucune image chargée")
            return
            
        import threading
        
        # Réinitialiser les variables de suivi de date pour chaque batch
        self._first_obs_date = None
        self._first_obs_jd = None
        self.zero_aperture_series_by_image = {}
        self._start_process_log_capture("astrométrie batch")
        
        def worker():
            try:
                total = len(self.image_files)
                success = 0
                failed = 0
                observations = []  # Collecte des observations pour les rapports
                
                for i, img_path in enumerate(self.image_files):
                    try:
                        logger.info(f"Traitement {i+1}/{total}: {Path(img_path).name}")
                        
                        # Résoudre l'astrométrie selon la méthode sélectionnée
                        method = self.astrometry_method_var.get()
                        if method == "classical":
                            # Utiliser l'astrométrie classique pour toutes les images
                            self.solve_astrometry_classical(img_path)
                        else:  # "zero_aperture" par défaut
                            # Toujours faire le test complet des apertures avec extrapolation zero-aperture
                            logger.info(f"Test complet des apertures zero-aperture pour {Path(img_path).name} (image {i+1}/{total})")
                            self.solve_astrometry_zero_aperture(img_path, skip_zero_aperture=False)
                        
                        # Collecte des données d'observation
                        try:
                            obs_data = self._collect_observation_data(img_path)
                            if obs_data:
                                observations.append(obs_data)
                        except Exception as e:
                            logger.warning(f"Impossible de collecter les données pour {Path(img_path).name}: {e}")
                        
                        success += 1
                    except Exception as e:
                        logger.error(f"Échec {Path(img_path).name}: {e}")
                        failed += 1
                
                # Génération des rapports ADES et JSON
                if observations and self.directory:
                    try:
                        self._generate_ades_report(observations)
                        self._generate_json_archive(observations)
                        self._generate_horizons_report()
                        self._generate_photometry_report(observations, method=self.astrometry_method_var.get())
                        # Mettre à jour le tableau zéro-ouverture (6 ouvertures/image)
                        self.frame.after(0, self.refresh_zero_aperture_table)
                        logger.info(f"Rapports générés : {len(observations)} observations")
                    except Exception as e:
                        logger.error(f"Erreur lors de la génération des rapports : {e}", exc_info=True)
                
                # Utilisation de 'self.after' pour éviter des soucis de thread
                report_msg = f"\n\nRapports générés ({len(observations)} observations)" if observations else ""
                self.frame.after(0, lambda: messagebox.showinfo("Batch terminé", 
                    f"Astrométrie batch terminée :\n{success} réussies\n{failed} échouées{report_msg}"))
            finally:
                self.frame.after(0, lambda: self._stop_process_log_capture("astrométrie batch"))
        
        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        messagebox.showinfo("Info", f"Batch astrométrie lancé en arrière-plan sur {len(self.image_files)} images")
    
    def _compute_deviant_threshold(self, values):
        """Calcule un seuil robuste de déviation (median + 2.5*MAD)."""
        if not values:
            return None
        arr = np.asarray(values, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return None
        med = float(np.median(arr))
        mad = float(np.median(np.abs(arr - med)))
        if mad <= 0:
            return float(np.percentile(arr, 85))
        return med + 2.5 * 1.4826 * mad
    
    def _build_zero_aperture_rows(self):
        """Construit les lignes de tableau: une ligne par image et ouverture."""
        rows = []
        if not self.image_files:
            return rows
        
        for img_path in self.image_files:
            obs = self._collect_observation_data(img_path)
            if not obs:
                continue
            image_key = str(Path(img_path).resolve())
            series = self.zero_aperture_series_by_image.get(image_key, [])
            if not series:
                # Fallback : 6 ouvertures standards si la série n'est pas disponible
                fallback_ap = np.linspace(2.0, 8.0, 6)
                for ap in fallback_ap:
                    rows.append({
                        'image': Path(img_path).name,
                        'obsTime': obs.get('date_obs', ''),
                        'photAp': float(ap),
                        'ra': obs.get('ra_deg'),
                        'dec': obs.get('dec_deg'),
                        'rmsRA': float(obs.get('rms_ra', 0.0)),
                        'rmsDec': float(obs.get('rms_dec', 0.0)),
                        'rmsFit': float(obs.get('rms_total', 0.0)),
                        'selected': True,
                    })
                continue
            
            for item in series[:6]:
                rows.append({
                    'image': Path(img_path).name,
                    'obsTime': obs.get('date_obs', ''),
                    'photAp': float(item.get('phot_ap', 0.0)),
                    'ra': obs.get('ra_deg'),
                    'dec': obs.get('dec_deg'),
                    'rmsRA': float(item.get('rms_ra', obs.get('rms_ra', 0.0))),
                    'rmsDec': float(item.get('rms_dec', obs.get('rms_dec', 0.0))),
                    'rmsFit': float(item.get('rms_fit', obs.get('rms_total', 0.0))),
                    'selected': bool(item.get('selected', True)),
                })
        return rows
    
    def _fill_zero_aperture_tree(self):
        """Remplit la table ZA et applique le surlignage des points déviants."""
        if not self.zero_aperture_tree:
            return
        tree = self.zero_aperture_tree
        for iid in tree.get_children():
            tree.delete(iid)
        
        self.zero_aperture_rows = self._build_zero_aperture_rows()
        dev_threshold = self._compute_deviant_threshold([r['rmsFit'] for r in self.zero_aperture_rows])
        
        for idx, row in enumerate(self.zero_aperture_rows):
            is_deviant = False
            if dev_threshold is not None and np.isfinite(row['rmsFit']):
                is_deviant = row['rmsFit'] >= dev_threshold
            tags = ("deviant",) if is_deviant else ()
            tree.insert(
                "",
                tk.END,
                iid=str(idx),
                values=(
                    row['image'],
                    row['obsTime'],
                    f"{row['photAp']:.2f}",
                    f"{row['ra']:.6f}" if row['ra'] is not None else "",
                    f"{row['dec']:.6f}" if row['dec'] is not None else "",
                    f"{row['rmsRA']:.3f}",
                    f"{row['rmsDec']:.3f}",
                    f"{row['rmsFit']:.3f}",
                    "1" if row['selected'] else "0",
                ),
                tags=tags,
            )
    
    def _on_zero_aperture_double_click(self, event):
        """Bascule selected (0/1) au double-clic sur la colonne selected."""
        if not self.zero_aperture_tree:
            return
        tree = self.zero_aperture_tree
        item = tree.identify_row(event.y)
        column = tree.identify_column(event.x)
        if not item or column != "#9":
            return
        try:
            idx = int(item)
        except (TypeError, ValueError):
            return
        if idx < 0 or idx >= len(self.zero_aperture_rows):
            return
        self.zero_aperture_rows[idx]['selected'] = not self.zero_aperture_rows[idx]['selected']
        vals = list(tree.item(item, "values"))
        vals[8] = "1" if self.zero_aperture_rows[idx]['selected'] else "0"
        tree.item(item, values=vals)
        self._sync_zero_aperture_selection_to_series(self.zero_aperture_rows[idx])
    
    def _sync_zero_aperture_selection_to_series(self, row):
        """Synchronise la sélection UI vers la série interne par image/ouverture."""
        try:
            image_name = row.get('image', '')
            phot_ap = float(row.get('photAp', 0.0))
            selected = bool(row.get('selected', True))
            img_key = None
            for p in self.image_files:
                if Path(p).name == image_name:
                    img_key = str(Path(p).resolve())
                    break
            if not img_key:
                return
            series = self.zero_aperture_series_by_image.get(img_key, [])
            if not isinstance(series, list) or len(series) == 0:
                return
            # Met à jour l'entrée d'ouverture la plus proche
            best_idx = None
            best_dist = None
            for i, s in enumerate(series):
                ap = float(s.get('phot_ap', 0.0))
                d = abs(ap - phot_ap)
                if best_dist is None or d < best_dist:
                    best_dist = d
                    best_idx = i
            if best_idx is not None:
                series[best_idx]['selected'] = selected
        except Exception as e:
            logger.debug(f"Sync selected ZA impossible: {e}")
    
    def refresh_zero_aperture_table(self):
        """Rafraîchit le tableau ZA si la fenêtre est ouverte."""
        if self.zero_aperture_window and self.zero_aperture_window.winfo_exists():
            self._fill_zero_aperture_tree()
    
    def _set_all_zero_aperture_selected(self, selected: bool):
        """Applique selected (0/1) à toutes les lignes du tableau ZA."""
        if not self.zero_aperture_rows:
            return
        for idx, row in enumerate(self.zero_aperture_rows):
            row['selected'] = bool(selected)
            self._sync_zero_aperture_selection_to_series(row)
            if self.zero_aperture_tree and self.zero_aperture_tree.exists(str(idx)):
                vals = list(self.zero_aperture_tree.item(str(idx), "values"))
                vals[8] = "1" if selected else "0"
                self.zero_aperture_tree.item(str(idx), values=vals)
    
    def _get_results_dir(self):
        """Retourne le dossier results et le crée si nécessaire."""
        if not self.directory:
            return None
        results_dir = Path(self.directory) / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        return results_dir
    
    def _open_results_folder(self):
        """Ouvre le dossier results dans l'explorateur."""
        results_dir = self._get_results_dir()
        if not results_dir:
            messagebox.showerror("Erreur", "Aucun dossier d'images actif.")
            return
        try:
            os.startfile(str(results_dir))
        except Exception:
            try:
                webbrowser.open(results_dir.as_uri())
            except Exception as e:
                logger.error(f"Impossible d'ouvrir le dossier results: {e}")
                messagebox.showerror("Erreur", f"Impossible d'ouvrir le dossier results:\n{e}")
    
    def _export_zero_aperture_full_table(self):
        """Exporte le tableau complet ZA (toutes ouvertures, selected inclus) dans results."""
        results_dir = self._get_results_dir()
        if not results_dir:
            return None
        rows = self._build_zero_aperture_rows()
        if not rows:
            return None
        df = pd.DataFrame(rows, columns=["image", "obsTime", "photAp", "ra", "dec", "rmsRA", "rmsDec", "rmsFit", "selected"])
        asteroid_id = self.asteroid_id_var.get().strip() or "UNKNOWN"
        output_path = results_dir / f"{asteroid_id}_ZA_table_complete.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Tableau ZA complet exporté : {output_path} ({len(df)} lignes)")
        return output_path
    
    def _compute_zero_aperture_extrapolated_observations(self):
        """
        Construit des observations finales extrapolées à ouverture 0 (une ligne par image),
        en utilisant uniquement les ouvertures selected=1 du tableau ZA.
        """
        if not self.image_files:
            return []
        
        # Récupérer les observations de base (temps, mag, etc.)
        base_obs = {}
        for img_path in self.image_files:
            obs = self._collect_observation_data(img_path)
            if obs:
                base_obs[Path(img_path).name] = obs
        
        if not base_obs:
            return []
        
        def _fit_intercept(x_vals, y_vals, w_vals=None):
            x = np.asarray(x_vals, dtype=float)
            y = np.asarray(y_vals, dtype=float)
            if x.size == 0 or y.size == 0:
                return 0.0
            if x.size == 1:
                return float(y[0])
            try:
                if w_vals is not None:
                    w = np.asarray(w_vals, dtype=float)
                    w = 1.0 / np.maximum(w, 1e-10)
                    coeff = np.polyfit(x, y, 1, w=w)
                else:
                    coeff = np.polyfit(x, y, 1)
                return float(coeff[1])  # ordonnée à l'origine (aperture=0)
            except Exception:
                return float(np.mean(y))
        
        # Magnitudes photométriques (image par image) pour enrichir ADES ZA0
        phot_mag_map = self._load_photometry_magnitude_map()

        extrapolated = []
        for img_path in self.image_files:
            img_name = Path(img_path).name
            obs = base_obs.get(img_name)
            if not obs:
                continue
            img_key = str(Path(img_path).resolve())
            series = self.zero_aperture_series_by_image.get(img_key, [])
            if not isinstance(series, list) or len(series) == 0:
                continue
            
            selected = [s for s in series[:6] if bool(s.get('selected', True))]
            if len(selected) < 1:
                continue
            
            ap = [float(s.get('phot_ap', 0.0)) for s in selected]
            rms_fit = [float(s.get('rms_fit', obs.get('rms_total', 0.2))) for s in selected]
            rms_ra = [float(s.get('rms_ra', obs.get('rms_ra', 0.2))) for s in selected]
            rms_dec = [float(s.get('rms_dec', obs.get('rms_dec', 0.2))) for s in selected]
            # Résidus moyens par ouverture (ajoutés côté solveur)
            res_ra = [float(s.get('mean_residual_ra_arcsec', 0.0)) for s in selected]
            res_dec = [float(s.get('mean_residual_dec_arcsec', 0.0)) for s in selected]
            
            offset_ra0_arcsec = _fit_intercept(ap, res_ra, rms_fit)
            offset_dec0_arcsec = _fit_intercept(ap, res_dec, rms_fit)
            rms_fit0 = max(0.0, _fit_intercept(ap, rms_fit, rms_fit))
            rms_ra0 = max(0.0, _fit_intercept(ap, rms_ra, rms_fit))
            rms_dec0 = max(0.0, _fit_intercept(ap, rms_dec, rms_fit))
            
            ra_base = float(obs.get('ra_deg', 0.0))
            dec_base = float(obs.get('dec_deg', 0.0))
            cos_dec = np.cos(np.radians(dec_base))
            if abs(cos_dec) < 1e-8:
                cos_dec = 1e-8
            
            # Même convention de signe que dans le solveur ZA (on enlève l'offset)
            ra0 = ra_base - (offset_ra0_arcsec / (3600.0 * cos_dec))
            dec0 = dec_base - (offset_dec0_arcsec / 3600.0)
            
            obs_out = dict(obs)
            obs_out['ra_deg'] = float(ra0)
            obs_out['dec_deg'] = float(dec0)
            obs_out['rms_ra'] = float(rms_ra0)
            obs_out['rms_dec'] = float(rms_dec0)
            obs_out['rms_corr'] = float(obs.get('rms_corr', 0.0))
            obs_out['notes'] = "e"
            obs_out['remarks'] = f"ZA0 extrapolation ({len(selected)} apertures)"
            obs_out['za_selected_count'] = len(selected)
            obs_out['za_rms_fit0'] = float(rms_fit0)
            if phot_mag_map and img_name in phot_mag_map:
                m = phot_mag_map[img_name]
                obs_out['mag'] = m.get('mag')
                obs_out['rms_mag'] = m.get('rms_mag')
                obs_out['log_snr'] = m.get('log_snr')
                obs_out['filter'] = 'G'
                obs_out['phot_cat'] = 'Gaia'
            extrapolated.append(obs_out)
        
        return extrapolated
    
    def _apply_photometry_mags_to_observations(self, observations):
        """Injecte mag_T1_G/rmsMag_T1 du CSV photométrie dans une liste d'observations."""
        if not observations:
            return observations
        mag_map = self._load_photometry_magnitude_map()
        if not mag_map:
            return observations
        enriched = []
        for obs in observations:
            obs2 = dict(obs)
            fn = obs2.get("filename")
            if fn in mag_map:
                m = mag_map[fn]
                obs2["mag"] = m.get("mag")
                obs2["rms_mag"] = m.get("rms_mag")
                obs2["log_snr"] = m.get("log_snr")
                obs2["filter"] = "G"
                obs2["phot_cat"] = "Gaia"
            enriched.append(obs2)
        return enriched
    
    def _export_zero_aperture_extrapolated_final(self, show_message=True):
        """Génère le résultat final ZA(0) puis exporte un ADES final dédié."""
        if not self.directory or not self.image_files:
            if show_message:
                messagebox.showerror("Erreur", "Chargez d'abord un dossier d'images FITS.")
            return None, None
        
        try:
            # Export tableau complet ZA pour audit
            self._export_zero_aperture_full_table()
            extrapolated_obs = self._compute_zero_aperture_extrapolated_observations()
            if not extrapolated_obs:
                if show_message:
                    messagebox.showerror(
                        "Erreur",
                        "Aucune observation extrapolée ZA(0) générée.\n"
                        "Vérifiez les sélections selected=1 dans le tableau."
                    )
                return None, None
            
            asteroid_id = self.asteroid_id_var.get().strip() or "UNKNOWN"
            # Archive CSV dédiée ZA(0)
            results_dir = self._get_results_dir()
            csv_path = None
            if results_dir is not None:
                df = pd.DataFrame(extrapolated_obs)
                csv_path = results_dir / f"{asteroid_id}_ZA0_extrapolated.csv"
                df.to_csv(csv_path, index=False)
            
            ades_path = self._generate_ades_report(
                extrapolated_obs,
                apply_selected_filter=False,
                output_filename=f"{asteroid_id}_ADES_final_MPC_ZA0.psv",
            )
            
            if show_message:
                msg = "Résultat final ZA(0) généré :"
                if csv_path:
                    msg += f"\n- Détails extrapolation : {csv_path}"
                if ades_path:
                    msg += f"\n- ADES final ZA(0) : {ades_path}"
                messagebox.showinfo("ZA(0) terminé", msg)
            return ades_path, csv_path
        except Exception as e:
            logger.error(f"Erreur génération résultat final ZA(0): {e}", exc_info=True)
            if show_message:
                messagebox.showerror("Erreur", f"Impossible de générer le résultat final ZA(0) :\n{e}")
            return None, None
    
    def _export_ades_final_mpc(self, show_message=True):
        """Génère le rapport ADES final à partir des images chargées."""
        if not self.directory or not self.image_files:
            if show_message:
                messagebox.showerror("Erreur", "Chargez d'abord un dossier d'images FITS.")
            return None
        
        observations = []
        for img_path in self.image_files:
            obs = self._collect_observation_data(img_path)
            if obs:
                observations.append(obs)
        
        if not observations:
            if show_message:
                messagebox.showerror("Erreur", "Aucune observation exploitable pour générer le rapport ADES.")
            return None
        
        try:
            za_path = self._export_zero_aperture_full_table()
            observations = self._apply_photometry_mags_to_observations(observations)
            ades_path = self._generate_ades_report(observations)
            if show_message:
                msg = "Rapports générés dans le dossier results :"
                if za_path:
                    msg += f"\n- Tableau ZA complet : {za_path}"
                if ades_path:
                    msg += f"\n- ADES final MPC : {ades_path}"
                messagebox.showinfo("Exports results", msg)
            return ades_path
        except Exception as e:
            logger.error(f"Erreur génération ADES final MPC: {e}", exc_info=True)
            if show_message:
                messagebox.showerror("Erreur", f"Impossible de générer le rapport ADES final :\n{e}")
            return None

    def _create_final_ades_reports(self):
        """Génère ADES et ADES_ZA0 finaux avec mag_T1_G intégrée si disponible."""
        ades_path = self._export_ades_final_mpc(show_message=False)
        ades_za0_path, za0_csv = self._export_zero_aperture_extrapolated_final(show_message=False)
        if not ades_path and not ades_za0_path:
            messagebox.showerror("Erreur", "Aucun rapport ADES généré.")
            return
        msg = "Création des rapports ADES terminée :"
        if ades_path:
            msg += f"\n- ADES final : {ades_path}"
        if ades_za0_path:
            msg += f"\n- ADES final ZA0 : {ades_za0_path}"
        if za0_csv:
            msg += f"\n- Détails ZA0 : {za0_csv}"
        messagebox.showinfo("Rapports ADES", msg)
    
    def open_zero_aperture_table(self):
        """Ouvre la fenêtre de tableau zéro-ouverture (6 lignes/image)."""
        if self.zero_aperture_window and self.zero_aperture_window.winfo_exists():
            self.zero_aperture_window.lift()
            self.refresh_zero_aperture_table()
            return
        
        win = tk.Toplevel(self.frame.winfo_toplevel())
        win.title("Tableau Astrométrie Zéro-Ouverture")
        win.geometry("1200x430")
        self.zero_aperture_window = win
        
        ttk.Label(
            win,
            text="Colonnes: image, obsTime, photAp, ra, dec, rmsRA, rmsDec, rmsFit, selected (6 ouvertures/image)",
            font=("Helvetica", 9),
        ).pack(anchor="w", padx=8, pady=(8, 4))
        
        tree_frame = ttk.Frame(win)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)
        
        cols = ("image", "obsTime", "photAp", "ra", "dec", "rmsRA", "rmsDec", "rmsFit", "selected")
        tree = ttk.Treeview(tree_frame, columns=cols, show="headings", height=14, selectmode="browse")
        self.zero_aperture_tree = tree
        
        widths = {
            "image": 170,
            "obsTime": 210,
            "photAp": 75,
            "ra": 110,
            "dec": 110,
            "rmsRA": 75,
            "rmsDec": 75,
            "rmsFit": 75,
            "selected": 70,
        }
        for c in cols:
            tree.heading(c, text=c)
            tree.column(c, width=widths.get(c, 90), anchor="center")
        
        tree.tag_configure("deviant", background="#ffd6d6")
        tree.bind("<Double-1>", self._on_zero_aperture_double_click)
        
        yscroll = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=tree.yview)
        xscroll = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL, command=tree.xview)
        tree.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
        tree.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)
        
        btn_frame = ttk.Frame(win)
        btn_frame.pack(fill=tk.X, padx=8, pady=(0, 8))
        ttk.Button(btn_frame, text="🧾 Création des rapports ADES", command=self._create_final_ades_reports).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="📝 Export ADES final (MPC)", command=self._export_ades_final_mpc).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="🧮 Générer résultat final ZA(0)", command=self._export_zero_aperture_extrapolated_final).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="📤 Ouvrir dossier results", command=self._open_results_folder).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Tout sélectionner", command=lambda: self._set_all_zero_aperture_selected(True)).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Tout désélectionner", command=lambda: self._set_all_zero_aperture_selected(False)).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Rafraîchir", command=self.refresh_zero_aperture_table).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Fermer", command=win.destroy).pack(side=tk.RIGHT, padx=4)
        
        self._fill_zero_aperture_tree()
    
    def _collect_observation_data(self, img_path):
        """Collecte les données d'observation d'une image pour les rapports."""
        try:
            with open_fits_with_fixed_wcs(img_path) as hdul:
                header = hdul[0].header.copy()
                wcs = WCS(hdul[0].header)
            
            # Date d'observation - PRIORITÉ dans l'ordre (Horizons utilise JD-UTC) :
            # 1. JD-UTC au milieu de l'exposition (compatible avec Horizons)
            # 2. BJD-TDB au milieu de l'exposition
            # 3. DATE-OBS (début d'exposition, nécessite calcul du milieu)
            # 4. Extraction depuis le nom de fichier
            
            exptime = float(header.get('EXPTIME', header.get('EXPOSURE', 0.0)))
            time_used = False
            
            # 1. Essayer JD-UTC d'abord (compatible avec éphémérides JPL Horizons)
            # Chercher différentes variations du nom de champ JD-UTC
            jd_utc_raw = (header.get('JD-UTC', None) or 
                         header.get('JD_UTC', None) or
                         header.get('JD-UTC at mid exposure', None) or
                         header.get('JD_UTC_at_mid_exposure', None) or
                         header.get('JDUTC', None))
            if jd_utc_raw is not None:
                try:
                    jd_utc_val = float(jd_utc_raw)
                    t_mid = Time(jd_utc_val, format='jd', scale='utc')
                    jd_utc = jd_utc_val
                    date_obs_str = t_mid.isot
                    time_used = True
                    logger.info(f"JD-UTC lu depuis header pour {Path(img_path).name}: JD={jd_utc_val:.10f}, UTC={date_obs_str}")
                except (ValueError, TypeError) as e:
                    logger.warning(f"Erreur parsing JD-UTC pour {Path(img_path).name}: {e}")
            
            # 2. Si JD-UTC non disponible, essayer BJD-TDB au milieu de l'exposition
            if not time_used:
                bjd_tdb_raw = header.get('BJD-TDB', header.get('BJD_TDB', None))
                if bjd_tdb_raw is not None:
                    try:
                        bjd_tdb_val = float(bjd_tdb_raw)
                        t_mid_tdb = Time(bjd_tdb_val, format='jd', scale='tdb')
                        # Convertir en UTC pour compatibilité avec les éphémérides Horizons
                        t_mid = t_mid_tdb.utc
                        jd_utc = t_mid.jd
                        date_obs_str = t_mid.isot
                        time_used = True
                        logger.info(f"BJD-TDB lu depuis header pour {Path(img_path).name}: BJD={bjd_tdb_val:.10f}, UTC={date_obs_str}, JD={jd_utc:.10f}")
                    except (ValueError, TypeError) as e:
                        logger.warning(f"Erreur parsing BJD-TDB pour {Path(img_path).name}: {e}")
            
            # 3. Si ni JD-UTC ni BJD-TDB disponibles, utiliser DATE-OBS ou extraction depuis nom de fichier
            if not time_used:
                date_obs_str = header.get('DATE-OBS', header.get('DATE'))
                
                # Si DATE-OBS manquant, essayer d'extraire depuis le nom de fichier
                # Format attendu : LIGHT_17-Thetis_2025-12-12-05-34-45_FilterV_Exp45.00_GAIN0_solved.fits
                if not date_obs_str:
                    filename = Path(img_path).name
                    # Pattern : YYYY-MM-DD-HH-MM-SS
                    date_match = re.search(r'(\d{4})-(\d{2})-(\d{2})-(\d{2})-(\d{2})-(\d{2})', filename)
                    if date_match:
                        year, month, day, hour, minute, second = date_match.groups()
                        date_obs_str = f"{year}-{month}-{day}T{hour}:{minute}:{second}"
                        logger.info(f"Date extraite depuis nom fichier pour {filename}: {date_obs_str}")
                    else:
                        logger.warning(f"Pas de DATE-OBS dans {filename} et impossible d'extraire depuis le nom")
                        return None
                
                # Calculer le temps au milieu de l'exposition
                t_start = Time(date_obs_str, scale='utc')
                t_mid = t_start + exptime / 2.0 * u.second
                jd_utc = t_mid.jd
            
            # Log pour vérifier que chaque image a bien une date différente
            logger.info(f"Collecte {Path(img_path).name}: DATE-OBS={date_obs_str}, JD={jd_utc:.10f}, t_mid={t_mid.isot}")
            
            # Position T1 - calcul pour chaque image avec son WCS
            target_ra = None
            target_dec = None
            target_mag = None
            
            # Interpolation depuis éphémérides (RA/Dec observés) même si WCS manquant
            if self.ephemeris_data:
                try:
                    # 1. Interpoler position RA/Dec depuis éphémérides pour le temps d'observation
                    target_sky_eph = _interpolate_target_position(self.ephemeris_data, t_mid.jd)
                    if target_sky_eph:
                        # Utiliser directement les coordonnées interpolées depuis les éphémérides
                        # (les éphémérides JPL Horizons sont précises et tiennent compte du mouvement de l'astéroïde)
                        target_ra = target_sky_eph.ra.deg
                        target_dec = target_sky_eph.dec.deg
                        
                        # Pour le log, convertir en pixels pour vérification
                        try:
                            if wcs and wcs.is_celestial:
                                tx_px, ty_px = wcs.world_to_pixel(target_sky_eph)
                                # Log avec plus de précision pour voir les différences subtiles
                                logger.info(f"Position T1 pour {Path(img_path).name}: JD={t_mid.jd:.15f}, RA={target_ra:.15f}°, Dec={target_dec:.15f}° (px={tx_px:.2f},{ty_px:.2f})")
                            else:
                                logger.info(f"Position T1 pour {Path(img_path).name}: JD={t_mid.jd:.15f}, RA={target_ra:.15f}°, Dec={target_dec:.15f}°")
                        except Exception:
                            logger.info(f"Position T1 pour {Path(img_path).name}: JD={t_mid.jd:.15f}, RA={target_ra:.15f}°, Dec={target_dec:.15f}°")
                except Exception as e:
                    logger.warning(f"Impossible de calculer position T1 pour {Path(img_path).name}: {e}")
            
            # Sans éphémérides : interpolation depuis ancres manuelles T1 (première + dernière image)
            if target_ra is None and target_dec is None and self.manual_t1_anchor_first and self.manual_t1_anchor_last:
                try:
                    jd1 = self.manual_t1_anchor_first.get("jd")
                    jd2 = self.manual_t1_anchor_last.get("jd")
                    if jd1 is not None and jd2 is not None and jd2 != jd1:
                        frac = (jd_utc - jd1) / (jd2 - jd1)
                        c1 = SkyCoord(
                            self.manual_t1_anchor_first["ra_deg"] * u.deg,
                            self.manual_t1_anchor_first["dec_deg"] * u.deg,
                            frame="icrs"
                        )
                        c2 = SkyCoord(
                            self.manual_t1_anchor_last["ra_deg"] * u.deg,
                            self.manual_t1_anchor_last["dec_deg"] * u.deg,
                            frame="icrs"
                        )
                        v1 = c1.cartesian.xyz.value
                        v2 = c2.cartesian.xyz.value
                        vx = v1[0] + frac * (v2[0] - v1[0])
                        vy = v1[1] + frac * (v2[1] - v1[1])
                        vz = v1[2] + frac * (v2[2] - v1[2])
                        coord_cart = SkyCoord(x=vx, y=vy, z=vz, representation_type='cartesian', frame='icrs')
                        coord_sph = coord_cart.represent_as('spherical')
                        target_ra = float(coord_sph.lon.deg)
                        target_dec = float(coord_sph.lat.deg)
                        logger.debug(f"Position T1 interpolée (ancres manuelles) pour {Path(img_path).name}: RA={target_ra:.6f}°, Dec={target_dec:.6f}°")
                except Exception as e:
                    logger.debug(f"Interpolation ancres T1 pour {Path(img_path).name}: {e}")
            
            # Magnitude depuis éphémérides (toujours nécessaire)
            if self.ephemeris_data:
                try:
                    target_mag = _interpolate_magnitude_V(self.ephemeris_data, t_mid)
                except Exception as e:
                    logger.debug(f"Impossible d'interpoler magnitude T1: {e}")
            
            # Airmass
            airmass = header.get('AIRMASS', None)
            if airmass is None and target_ra is not None and target_dec is not None:
                try:
                    from core.photometry_pipeline_asteroids import airmass
                    airmass = airmass(target_ra, target_dec, t_mid)
                except:
                    airmass = 1.0
            if airmass is None:
                airmass = 1.0
            
            # Récupérer toutes les données RMS depuis le header (stocké après astrométrie)
            rms_total = header.get('ASTRRMS', 0.2)  # RMS total
            rms_ra = header.get('ASTRRMSR', header.get('ASTRRMS', 0.2))  # RMS RA, sinon RMS global, sinon 0.2 par défaut
            rms_dec = header.get('ASTRRMSD', header.get('ASTRRMS', 0.2))  # RMS Dec, sinon RMS global, sinon 0.2 par défaut
            rms_corr = header.get('ASTRCOR', 0.0)  # Corrélation RMS (stockée après astrométrie)
            rms_zero_aperture = header.get('ASTRMSZ', None)  # RMS zéro-aperture
            rms_classical = header.get('ASTRMSC', None)  # RMS classique
            method = header.get('ASTMETHOD', 'unknown')  # Méthode utilisée
            fit_r2 = header.get('ASTR2', None)  # R² de l'extrapolation
            median_residual = header.get('ASTRMED', None)  # Résidu médian
            n_stars = header.get('ASTNREF', None)  # Nombre d'étoiles de référence
            n_outliers = header.get('ASTNOUT', None)  # Nombre d'outliers
            
            # S'assurer que les valeurs sont numériques
            try:
                rms_total = float(rms_total)
                rms_ra = float(rms_ra)
                rms_dec = float(rms_dec)
                rms_corr = float(rms_corr)
                if rms_zero_aperture is not None:
                    rms_zero_aperture = float(rms_zero_aperture)
                if rms_classical is not None:
                    rms_classical = float(rms_classical)
                if fit_r2 is not None:
                    fit_r2 = float(fit_r2)
                if median_residual is not None:
                    median_residual = float(median_residual)
                if n_stars is not None:
                    n_stars = int(n_stars)
                if n_outliers is not None:
                    n_outliers = int(n_outliers)
            except (ValueError, TypeError):
                rms_total = 0.2
                rms_ra = 0.2
                rms_dec = 0.2
                rms_corr = 0.0
            
            return {
                'filename': Path(img_path).name,
                'jd_utc': jd_utc,
                'date_obs': t_mid.isot,
                'exptime': exptime,
                'ra_deg': target_ra,
                'dec_deg': target_dec,
                'mag': target_mag,
                'airmass': float(airmass),
                'filter': header.get('FILTER', header.get('FILT', '')),
                'wcs_valid': wcs.is_celestial if wcs else False,
                'rms_total': rms_total,
                'rms_ra': rms_ra,
                'rms_dec': rms_dec,
                'rms_corr': rms_corr,  # Corrélation RMS (coefficient de corrélation RA/Dec)
                'rms_zero_aperture': rms_zero_aperture,
                'rms_classical': rms_classical,
                'astrometry_method': method,
                'fit_r2': fit_r2,
                'median_residual': median_residual,
                'n_stars': n_stars,
                'n_outliers': n_outliers
            }
        except Exception as e:
            logger.error(f"Erreur collecte données observation {Path(img_path).name}: {e}")
            return None
            
    def _generate_ades_report(self, observations, apply_selected_filter=True, output_filename=None):
        """Génère un rapport ADES au format PSV (Pipe-Separated Values) dans le dossier des images."""
        if not self.directory or not observations:
            return None

        import config
        from datetime import datetime
        
        asteroid_id = self.asteroid_id_var.get().strip()
        if not asteroid_id:
            asteroid_id = "UNKNOWN"
        
        obs_code = self.observatory_code_var.get().strip().upper()
        if not obs_code:
            obs_code = "XXX"  # Code par défaut
        
        obs_name = config.OBSERVATORY.get('name', 'Observatoire Personnel')
        
        results_dir = self._get_results_dir()
        if results_dir is None:
            return None
        output_path = results_dir / (output_filename or f"{asteroid_id}_ADES_final_MPC.psv")
        
        # En mode zero-aperture, exporter uniquement les ouvertures marquées selected=1.
        # Chaque ouverture sélectionnée génère une ligne ADES.
        observations_to_export = []
        method = self.astrometry_method_var.get()
        if method == "zero_aperture" and apply_selected_filter:
            image_key_by_name = {}
            for p in self.image_files:
                image_key_by_name[Path(p).name] = str(Path(p).resolve())
            
            for obs in observations:
                img_name = obs.get('filename', '')
                img_key = image_key_by_name.get(img_name)
                series = self.zero_aperture_series_by_image.get(img_key, []) if img_key else []
                if isinstance(series, list) and len(series) > 0:
                    selected_series = [s for s in series[:6] if bool(s.get('selected', True))]
                    for s in selected_series:
                        obs_copy = dict(obs)
                        obs_copy['rms_ra'] = float(s.get('rms_ra', obs_copy.get('rms_ra', 0.2)))
                        obs_copy['rms_dec'] = float(s.get('rms_dec', obs_copy.get('rms_dec', 0.2)))
                        obs_copy['rms_corr'] = float(s.get('rms_fit', obs_copy.get('rms_corr', 0.0)))
                        # Signalement ZA (convention tutoriel/IAWN-MPC)
                        notes_prev = (obs_copy.get('notes') or '').strip()
                        obs_copy['notes'] = notes_prev if notes_prev else "e"
                        observations_to_export.append(obs_copy)
                else:
                    # Fallback : si aucune série ZA disponible, garder l'observation brute
                    obs_copy = dict(obs)
                    notes_prev = (obs_copy.get('notes') or '').strip()
                    obs_copy['notes'] = notes_prev if notes_prev else "e"
                    observations_to_export.append(obs_copy)
        else:
            observations_to_export = observations
        
        # En-têtes ADES selon le format standard
        psv_lines = [
            "# version=2017",
            "# observatory",
            f"! mpcCode {obs_code}",
            f"! name {obs_name}",
            "# submitter",
            "! name NPOAP",
            "! institution NPOAP",
            "# telescope",
            "! aperture 0.0",  # À compléter si disponible
            "! design Reflector",  # À compléter si disponible
            "! detector CCD",
            "# observers",
            "! name NPOAP",
            "# measurers",
            "! name NPOAP",
            "# comment",
            "! line Generated by NPOAP",
            "",  # Ligne vide
            # En-tête des colonnes (format aligné avec espaces)
            "permID |trkSub |mode|stn |obsTime                |ra         |dec        |rmsRA|rmsDec|rmsCorr|astCat|mag  |rmsMag|band|photCat|logSNR|notes|remarks"
        ]
        
        for obs in observations_to_export:
            if obs['ra_deg'] is None or obs['dec_deg'] is None:
                    continue

            # RA en degrés (format ADES standard)
            ra_deg = obs['ra_deg']
            dec_deg = obs['dec_deg']
            
            # Format date ADES (YYYY-MM-DDThh:mm:ss.sssZ)
            date_ades = obs['date_obs'].replace(' ', 'T')
            if '.' not in date_ades.split('T')[1]:
                date_ades += '.000'
            if not date_ades.endswith('Z'):
                date_ades += 'Z'
            
            # RMS par défaut (en arcsec)
            rms_ra = obs.get('rms_ra', 0.2)
            rms_dec = obs.get('rms_dec', 0.2)
            rms_corr = obs.get('rms_corr', 0.0)  # Corrélation RMS (par défaut 0)
            
            # Magnitude (ne pas afficher si valeur par défaut 99.0 = magnitude non disponible)
            mag = obs.get('mag', None)
            if mag is not None and abs(mag - 99.0) < 0.1:  # Valeur par défaut = magnitude non disponible
                mag_str = ""
            else:
                mag_str = f"{mag:.2f}" if mag is not None else ""
            rms_mag = obs.get('rms_mag', 0.0)
            rms_mag_str = f"{rms_mag:.2f}" if rms_mag > 0 else ""
            
            # Bande/filtre
            band = obs.get('filter', '').strip()
            if not band:
                band = "V"  # Par défaut V si non spécifié
            
            # PhotCat (catalogue photométrique)
            phot_cat = "Gaia"  # Par défaut Gaia
            
            # logSNR (log du rapport signal/bruit)
            log_snr = obs.get('log_snr', 0.0)
            log_snr_str = f"{log_snr:.3f}" if log_snr > 0 else ""
            
            # Notes et remarks
            notes = obs.get('notes', '').strip()
            remarks = obs.get('remarks', '').strip()
            
            # Construction de la ligne PSV avec alignement (format du modèle)
            # Format: permID |trkSub |mode|stn |obsTime                |ra         |dec        |rmsRA|rmsDec|rmsCorr|astCat|mag  |rmsMag|band|photCat|logSNR|notes|remarks
            line_parts = [
                f"{asteroid_id:<8}",      # permID (aligné à gauche, 8 caractères)
                "     ",                  # trkSub (vide, 5 espaces)
                "CCD",                    # mode
                f"{obs_code:<4}",         # stn (aligné à gauche, 4 caractères)
                f"{date_ades:<24}",       # obsTime (aligné à gauche, 24 caractères)
                f"{ra_deg:>10.6f}",      # ra (aligné à droite, 10 caractères, 6 décimales)
                f"{dec_deg:>10.6f}",     # dec (aligné à droite, 10 caractères, 6 décimales)
                f"{rms_ra:.2f}",         # rmsRA
                f"{rms_dec:.2f}",        # rmsDec
                f"{rms_corr:>6.2f}",     # rmsCorr (aligné à droite, 6 caractères)
                "Gaia",                   # astCat
                f"{mag_str:<5}" if mag_str else "     ",  # mag (aligné à gauche, 5 caractères)
                f"{rms_mag_str:<6}" if rms_mag_str else "      ",  # rmsMag (aligné à gauche, 6 caractères)
                f"{band:<5}",            # band (aligné à gauche, 5 caractères)
                f"{phot_cat:<8}",        # photCat (aligné à gauche, 8 caractères)
                f"{log_snr_str:<6}" if log_snr_str else "      ",  # logSNR (aligné à gauche, 6 caractères)
                notes[:5] if notes else "",  # notes (max 5 caractères)
                remarks[:50] if remarks else ""  # remarks (max 50 caractères)
            ]
            
            psv_lines.append("|".join(line_parts))
        
        # Créer le répertoire parent si nécessaire (évite FileNotFoundError)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(psv_lines))
        
        num_observations = len(psv_lines) - 18  # 18 lignes d'en-tête
        logger.info(f"Rapport ADES PSV généré : {output_path} ({num_observations} observations)")
        return output_path
    
    def _generate_json_archive(self, observations):
        """Génère un fichier JSON d'archive avec toutes les observations."""
        if not self.directory or not observations:
            return
        
        asteroid_id = self.asteroid_id_var.get().strip()
        if not asteroid_id:
            asteroid_id = "UNKNOWN"
        
        # Informations observatoire
        obs_info = {
            'name': config.OBSERVATORY.get('name', 'Unknown Observatory'),
            'lat': config.OBSERVATORY.get('lat', 0.0),
            'lon': config.OBSERVATORY.get('lon', 0.0),
            'elev_m': config.OBSERVATORY.get('elev', 0.0),
            'code': self.observatory_code_var.get().strip().upper() or "XXX"
        }
        
        # Structure JSON
        archive_data = {
            'asteroid_id': asteroid_id,
            'generated': datetime.now().isoformat(),
            'observatory': obs_info,
            'observation_count': len(observations),
            'observations': observations
        }
        
        output_path = Path(self.directory) / "rapport" / f"{asteroid_id}_observations.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(archive_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Archive JSON générée : {output_path}")
    
    def _generate_photometry_report(self, observations, method=None):
        """Génère un rapport TXT astrométrique (positions observées)."""
        if not self.directory or not observations:
            return None
        
        import config
        
        output_dir = self._get_results_dir()
        if output_dir is None:
            return None
        output_dir.mkdir(parents=True, exist_ok=True)
        asteroid_id = (self.asteroid_id_var.get() or "unknown").strip()
        if not asteroid_id:
            asteroid_id = "unknown"
        method_label = "ZA"
        if method == "classical":
            method_label = "FWHM"
        elif method == "zero_aperture":
            method_label = "ZA"
        else:
            # Fallback: tenter depuis les observations
            obs_method = observations[0].get('astrometry_method', '')
            method_label = "FWHM" if obs_method == "classical" else "ZA"
        
        output_file = output_dir / f"{asteroid_id}_astrometrie_{method_label}.txt"
        
        try:
            # Créer un DataFrame avec les données disponibles
            data = []
            for obs in observations:
                row = {
                    'filename': obs.get('filename', ''),
                    'JD-UTC': obs.get('jd_utc', ''),
                    'date_obs': obs.get('date_obs', ''),
                    'mag_ephemeris': obs.get('mag', None),  # Magnitude des éphémérides
                    'filter': obs.get('filter', ''),
                    'airmass': obs.get('airmass', None),
                    'ra_deg': obs.get('ra_deg', None),
                    'dec_deg': obs.get('dec_deg', None),
                }
                data.append(row)
            
            df = pd.DataFrame(data)
            
            # Réorganiser les colonnes
            cols_order = ['filename', 'JD-UTC', 'date_obs', 'ra_deg', 'dec_deg', 'mag_ephemeris', 'filter', 'airmass']
            df = df[[c for c in cols_order if c in df.columns]]
            
            # Sauvegarder en TXT (CSV lisible)
            df.to_csv(output_file, index=False)
            logger.info(f"Rapport astrométrie généré : {output_file}")
            return output_file
            
        except Exception as e:
            logger.error(f"Erreur génération rapport photométrie : {e}", exc_info=True)
            return None
    
    def _generate_horizons_report(self):
        """Génère un rapport Horizons au format TXT avec les éphémérides."""
        if not self.directory or not self.ephemeris_data:
            return
        
        import config
        
        asteroid_id = self.asteroid_id_var.get().strip()
        if not asteroid_id:
            asteroid_id = "UNKNOWN"
        
        obs_code = self.observatory_code_var.get().strip().upper()
        if not obs_code:
            obs_code = "XXX"
        
        obs_name = config.OBSERVATORY.get('name', 'Observatoire Personnel')
        obs_lat = config.OBSERVATORY.get('lat', 0.0)
        obs_lon = config.OBSERVATORY.get('lon', 0.0)
        obs_elev = config.OBSERVATORY.get('elev', 0.0)
        
        output_path = Path(self.directory) / "rapport" / f"{asteroid_id}_Horizons.txt"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # En-tête du rapport Horizons
        lines = [
            "=" * 80,
            f"EPHEMERIS REPORT - JPL HORIZONS FORMAT",
            f"Generated by NPOAP",
            f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "=" * 80,
            "",
            f"Target Body: Asteroid {asteroid_id}",
            f"Observatory: {obs_name} (Code: {obs_code})",
            f"Location: Lat={obs_lat:.6f}°, Lon={obs_lon:.6f}°, Elev={obs_elev:.1f}m",
            "",
            "-" * 80,
            "Date__(UT)__HR:MN     R.A._____(ICRF)_____DEC  APmag  S-brt    delta      deldot    S-O-T /r    S-T-O",
            "-" * 80
        ]
        
        # Données des éphémérides
        eph_data = self.ephemeris_data
        
        def _format_val(value, width, precision=2):
            if value is None:
                return "  N.A.".rjust(width)
            try:
                return f"{float(value):>{width}.{precision}f}"
            except Exception:
                return f"{str(value):>{width}}"
        
        def _col_value(table, key, idx):
            if key in table.colnames:
                return table[key][idx]
            return None

        for i in range(len(eph_data)):
            # Date au format JD -> UTC
            jd = eph_data['datetime_jd'][i]
            t = Time(jd, format='jd', scale='utc')
            date_str = t.strftime('%Y-%b-%d %H:%M')
            
            # RA en format HH MM SS.sss
            ra_deg = eph_data['RA'][i]
            ra_h = int(ra_deg / 15.0)
            ra_m = int((ra_deg / 15.0 - ra_h) * 60.0)
            ra_s = ((ra_deg / 15.0 - ra_h) * 60.0 - ra_m) * 60.0
            ra_str = f"{ra_h:02d} {ra_m:02d} {ra_s:05.2f}"
            
            # Dec en format +/-DD MM SS.ss
            dec_deg = eph_data['DEC'][i]
            dec_sign = '+' if dec_deg >= 0 else '-'
            dec_deg_abs = abs(dec_deg)
            dec_d = int(dec_deg_abs)
            dec_m = int((dec_deg_abs - dec_d) * 60.0)
            dec_s = ((dec_deg_abs - dec_d) * 60.0 - dec_m) * 60.0
            dec_str = f"{dec_sign}{dec_d:02d} {dec_m:02d} {dec_s:04.1f}"
            
            # Magnitude
            mag = eph_data['V'][i] if 'V' in eph_data.colnames else None
            mag_str = f"{mag:.2f}" if mag is not None else "  N.A."
            
            # Colonnes supplémentaires Horizons si disponibles
            s_brt = _format_val(_col_value(eph_data, 'S-brt', i), width=6, precision=2)
            delta = _format_val(_col_value(eph_data, 'delta', i), width=10, precision=6)
            deldot = _format_val(_col_value(eph_data, 'deldot', i), width=10, precision=6)
            s_o_t = _format_val(_col_value(eph_data, 'S-O-T', i), width=8, precision=2)
            s_t_o = _format_val(_col_value(eph_data, 'S-T-O', i), width=8, precision=2)
            
            # Format de la ligne (format Horizons standard)
            line = f"{date_str:<20} {ra_str:<15} {dec_str:<15} {mag_str:>6} {s_brt} {delta} {deldot} {s_o_t} {s_t_o}"
            lines.append(line)

        lines.extend([
            "-" * 80,
            "",
            f"Total ephemeris points: {len(eph_data)}",
            f"Time span: {Time(eph_data['datetime_jd'][0], format='jd').strftime('%Y-%b-%d %H:%M')} to {Time(eph_data['datetime_jd'][-1], format='jd').strftime('%Y-%b-%d %H:%M')}",
            "",
            "=" * 80
        ])
        
        with open(output_path, 'w', encoding='utf-8') as f:
                f.write('\n'.join(lines))
        
        logger.info(f"Rapport Horizons TXT généré : {output_path} ({len(eph_data)} points)")
        
    def _load_astrometry_positions(self):
        """Charge les positions astrométriques avec priorité ZA0, puis fallback JSON."""
        if not self.directory:
            return None
        
        asteroid_id = self.asteroid_id_var.get().strip()
        if not asteroid_id:
            asteroid_id = "UNKNOWN"

        def _build_map_from_observations(observations):
            positions_map = {}
            for obs in observations:
                filename = obs.get('filename')
                ra_deg = obs.get('ra_deg')
                dec_deg = obs.get('dec_deg')
                jd_utc = obs.get('jd_utc')
                if filename and ra_deg is not None and dec_deg is not None:
                    positions_map[filename] = {
                        'ra_deg': float(ra_deg),
                        'dec_deg': float(dec_deg),
                        'jd_utc': float(jd_utc) if jd_utc is not None else None
                    }
            return positions_map

        # 1) Priorité absolue : extrapolation ZA0 (results/*_ZA0_extrapolated.csv)
        try:
            results_dir = self._get_results_dir()
            if results_dir is not None:
                za0_path = results_dir / f"{asteroid_id}_ZA0_extrapolated.csv"
                if za0_path.exists():
                    df = pd.read_csv(za0_path)
                    obs = df.to_dict(orient='records')
                    positions_map = _build_map_from_observations(obs)
                    if positions_map:
                        logger.info(f"Positions ZA0 chargées : {len(positions_map)} images depuis {za0_path}")
                        return positions_map
                    logger.warning(f"Fichier ZA0 présent mais sans positions valides : {za0_path}")
        except Exception as e:
            logger.error(f"Erreur chargement positions ZA0 : {e}", exc_info=True)

        # 2) Fallback : archive JSON astrométrique classique
        json_path = Path(self.directory) / "rapport" / f"{asteroid_id}_observations.json"
        if not json_path.exists():
            logger.warning(f"Aucun fichier positions trouvé (ZA0/JSON). JSON absent : {json_path}")
            return None
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            observations = data.get('observations', [])
            if not observations:
                logger.warning(f"Aucune observation dans {json_path}")
                return None
            positions_map = _build_map_from_observations(observations)
            logger.info(f"Positions astrométriques (fallback JSON) chargées : {len(positions_map)} images depuis {json_path}")
            return positions_map if positions_map else None
        except Exception as e:
            logger.error(f"Erreur chargement positions astrométriques JSON : {e}", exc_info=True)
            return None
    
    def run_photometry_batch(self):
        """Lance la photométrie en batch sur toutes les images."""
        if not self.directory or not self.image_files:
            messagebox.showerror("Erreur", "Aucun dossier d'images chargé")
            return
        
        if not self.current_selections:
            messagebox.showerror("Erreur", "Aucune sélection photométrique configurée. Utilisez 'LANCER PHOTOMÉTRIE (T1)' d'abord pour configurer les apertures.")
            return
        
        if not self.current_image_path:
            messagebox.showerror("Erreur", "Aucune image de référence sélectionnée")
            return

        # Mode sans éphémérides: autoriser la photométrie avec suivi astrométrique/ancres/fixe
        astrometry_positions = self._load_astrometry_positions()
        has_manual_anchors = bool(self.manual_t1_anchor_first and self.manual_t1_anchor_last)
        has_ephemerides = bool(self.ephemeris_data is not None)
        if not has_ephemerides and not astrometry_positions and not has_manual_anchors:
            proceed = messagebox.askyesno(
                "Mode sans éphémérides",
                "Aucune éphéméride, aucune position astrométrique ADES et aucune ancre T1 (first/last) détectées.\n\n"
                "La photométrie utilisera la position T1 fixe de l'image de référence (peut être insuffisant si l'objet se déplace vite).\n\n"
                "Continuer quand même ?"
            )
            if not proceed:
                return
        
        # Confirmer avec l'utilisateur
        if not messagebox.askyesno("Confirmation", 
            f"Lancer la photométrie batch sur {len(self.image_files)} images ?\n\n"
            f"Référence: {Path(self.current_image_path).name}\n"
            f"Étoiles configurées: {len(self.current_selections)}"):
            return
        
        # Lancer dans un thread pour ne pas bloquer l'interface
        import threading
        import shutil
        from core.photometry_pipeline_asteroids import process_photometry_series
        self._start_process_log_capture("photométrie batch")
        
        def worker():
            try:
                logger.info(f"Début photométrie batch sur {len(self.image_files)} images")
                
                # Charger les positions astrométriques depuis le JSON
                astrometry_positions = self._load_astrometry_positions()
                if astrometry_positions:
                    logger.info(f"Utilisation des positions astrométriques ADES pour suivre T1 ({len(astrometry_positions)} images)")
                elif self.manual_t1_anchor_first and self.manual_t1_anchor_last:
                    logger.info("Utilisation des ancres manuelles T1 (first/last) pour suivre T1")
                elif self.ephemeris_data:
                    logger.info("Utilisation des éphémérides pour suivre T1")
                else:
                    logger.warning("Aucune éphéméride/ADES/ancres: utilisation de la position T1 fixe (image de référence)")
                
                # Trouver T1 dans les sélections
                target_coord = None
                comp_coords = []
                for sel in self.current_selections:
                    if sel['label'] == 'T1':
                        target_coord = sel['coord']
                    elif sel['label'].startswith('C'):
                        comp_coords.append(sel['coord'])
                
                if not target_coord:
                    self.frame.after(0, lambda: messagebox.showerror("Erreur", "T1 non trouvé dans les sélections"))
                    return

                # Les ancres ne sont requises que si on n'a ni positions astrométriques ni éphémérides.
                if (not astrometry_positions) and (self.ephemeris_data is None) and (not self.manual_t1_anchor_first or not self.manual_t1_anchor_last):
                    self.frame.after(
                        0,
                        lambda: messagebox.showerror(
                            "Erreur",
                            "Sans positions astrométriques ZA0/ADES ni éphémérides, cliquez T1 sur la première et la dernière image avant le batch."
                        )
                    )
                    return
                
                # Exécuter la photométrie avec les positions astrométriques (prioritaires) ou éphémérides
                process_photometry_series(
                    folder=self.directory,
                    target_coord=target_coord,
                    comps=comp_coords,
                    ref_image=self.current_image_path,
                    selections=self.current_selections,
                    progress_callback=lambda pct: logger.info(f"Progression: {pct:.1f}%"),
                    ephemeris_data=self.ephemeris_data,  # Éphémérides en fallback
                    astrometry_positions=astrometry_positions,  # Positions astrométriques ADES (prioritaires)
                    manual_t1_anchors={
                        "first": self.manual_t1_anchor_first,
                        "last": self.manual_t1_anchor_last
                    },
                    manual_aperture_overrides=self.photometry_manual_apertures,
                    manual_aperture_callback=self._request_manual_apertures_for_image
                )
                
                # Générer le fichier light_curve.txt
                result_path = Path(self.directory) / "photometrie" / "results.csv"
                light_curve_path = Path(self.directory) / "photometrie" / "light_curve.txt"
                try:
                    if result_path.exists():
                        df = pd.read_csv(result_path)
                        
                        # Vérifier que les colonnes nécessaires existent
                        if 'JD-UTC' in df.columns and 'rel_flux_T1_fn' in df.columns:
                            # Créer le fichier light_curve.txt avec Time, Relative_flux_fn, relative_flux_fn_err
                            with open(light_curve_path, 'w', encoding='utf-8') as f:
                                # En-tête
                                f.write("Time Relative_flux_fn relative_flux_fn_err\n")
                                
                                # Données
                                for idx, row in df.iterrows():
                                    time = row['JD-UTC']
                                    flux = row['rel_flux_T1_fn']
                                    # L'erreur normalisée est dans rel_flux_err_T1 (déjà normalisée par process_photometry_series)
                                    flux_err = row.get('rel_flux_err_T1', 0.0)
                                    
                                    f.write(f"{time:.6f} {flux:.6f} {flux_err:.6f}\n")
                            
                            logger.info(f"Fichier light_curve.txt généré : {light_curve_path}")
                        else:
                            logger.warning(f"Colonnes manquantes dans results.csv pour générer light_curve.txt")
                except Exception as e:
                    logger.error(f"Erreur lors de la génération de light_curve.txt : {e}", exc_info=True)
                
                # Message de succès
                compiled_csv = None
                try:
                    if result_path.exists():
                        results_dir = self._get_results_dir()
                        if results_dir is not None:
                            asteroid_id = (self.asteroid_id_var.get() or "UNKNOWN").strip() or "UNKNOWN"
                            compiled_csv = results_dir / f"{asteroid_id}_photometrie_compilee.csv"
                            shutil.copy2(result_path, compiled_csv)
                            logger.info(f"CSV photométrie compilé exporté : {compiled_csv}")
                except Exception as e:
                    logger.warning(f"Impossible d'exporter le CSV photométrie compilé : {e}")

                self.frame.after(0, lambda: messagebox.showinfo("Succès", 
                    f"Photométrie batch terminée !\n\nRésultats: {result_path}\nLight curve: {light_curve_path}"
                    + (f"\nCSV compilé: {compiled_csv}" if compiled_csv else "")))

            except Exception as e:
                logger.error(f"Erreur photométrie batch: {e}", exc_info=True)
                self.frame.after(0, lambda: messagebox.showerror("Erreur", f"Erreur lors de la photométrie batch:\n{e}"))
            finally:
                self.frame.after(0, lambda: self._stop_process_log_capture("photométrie batch"))
        
        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        messagebox.showinfo("Info", f"Photométrie batch lancée en arrière-plan sur {len(self.image_files)} images")

    def run_photometry_single_image(self):
        """Photométrie sur l'image courante (mode comètes, image par image)."""
        if not self.current_image_path or self.current_data is None or self.wcs is None:
            messagebox.showerror("Erreur", "Chargez d'abord une image avec WCS valide.")
            return
        if not self.current_selections:
            messagebox.showerror("Erreur", "Aucune sélection photométrique. Utilisez d'abord SET-UP PHOTOMÉTRIE.")
            return
        if not self._ensure_photometry_gaia_settings():
            return

        self._start_process_log_capture("photométrie image courante")
        try:
            data = self.current_data
            wcs = self.wcs
            header = self.current_header or {}
            gain = float(header.get("GAIN", 1.0) or 1.0)
            if gain <= 0:
                gain = 1.0

            jd_utc = header.get("JD-UTC", None)
            if jd_utc is None:
                try:
                    date_obs = header.get("DATE-OBS") or header.get("DATE")
                    jd_utc = Time(date_obs, scale="utc").jd if date_obs else np.nan
                except Exception:
                    jd_utc = np.nan

            row = {
                "filename": Path(self.current_image_path).name,
                "JD-UTC": float(jd_utc) if jd_utc is not None and np.isfinite(jd_utc) else np.nan,
                "date_obs": header.get("DATE-OBS", header.get("DATE", "")),
                "AIRMASS": float(header.get("AIRMASS", 1.0) or 1.0),
            }

            flux_t1 = 0.0
            err_t1_sq = 0.0
            flux_comps = 0.0
            var_comps = 0.0
            n_comps_valid = 0
            comp_mag_data = []

            for sel in self.current_selections:
                label = sel.get("label", "")
                coord = sel.get("coord")
                if coord is None:
                    continue

                r_ap = float(sel.get("r_ap", self.aperture_radius_var.get()))
                r_in = float(sel.get("r_in", self.annulus_inner_var.get()))
                r_out = float(sel.get("r_out", self.annulus_outer_var.get()))
                if not (r_ap > 0 and r_in > r_ap and r_out > r_in):
                    logger.warning(f"Apertures invalides pour {label}, ignoré.")
                    continue

                try:
                    px, py = wcs.world_to_pixel(coord)
                except Exception:
                    logger.debug(f"{label}: conversion WCS impossible.")
                    continue

                ny, nx = data.shape
                if not (0 <= px < nx and 0 <= py < ny):
                    logger.debug(f"{label}: hors champ ({px:.1f},{py:.1f}).")
                    continue

                # Centroid local : plus large pour T1 comète
                box = 45 if label == "T1" else 25
                px, py = refine_centroid(data, px, py, box_size=box)
                if not np.isfinite(px) or not np.isfinite(py) or not (0 <= px < nx and 0 <= py < ny):
                    logger.debug(f"{label}: centroid invalide.")
                    continue

                ap = CircularAperture((px, py), r=r_ap)
                an = CircularAnnulus((px, py), r_in=r_in, r_out=r_out)
                stats = ApertureStats(data, an)
                bkg_mean = float(stats.mean)
                bkg_std = float(stats.std) if np.isfinite(stats.std) else 0.0
                bg_area = float(stats.sum_aper_area.value) if np.isfinite(stats.sum_aper_area.value) else 1.0
                if bg_area <= 0:
                    bg_area = 1.0

                phot = aperture_photometry(data, ap)
                raw_flux = float(phot["aperture_sum"][0])
                ap_area = float(ap.area)
                net_flux = raw_flux - (bkg_mean * ap_area)

                term_source = max(net_flux, 0.0) / gain
                term_sky = ap_area * (bkg_std ** 2)
                term_bkg_err = (ap_area ** 2) * (bkg_std ** 2) / bg_area
                flux_err = float(np.sqrt(max(term_source + term_sky + term_bkg_err, 0.0)))

                row[f"X_{label}"] = float(px)
                row[f"Y_{label}"] = float(py)
                row[f"r_ap_{label}"] = float(r_ap)
                row[f"r_in_{label}"] = float(r_in)
                row[f"r_out_{label}"] = float(r_out)
                row[f"Flux_{label}"] = float(net_flux)
                row[f"Err_{label}"] = float(flux_err)
                row[f"Sky_{label}"] = float(bkg_mean)

                if label == "T1":
                    flux_t1 = net_flux
                    err_t1_sq = flux_err ** 2
                elif label.startswith("C") and np.isfinite(net_flux) and net_flux > 0:
                    flux_comps += net_flux
                    var_comps += flux_err ** 2
                    n_comps_valid += 1
                    gaia_gmag = self._query_gaia_gmag_for_coord(coord)
                    if gaia_gmag is not None and np.isfinite(gaia_gmag):
                        comp_mag_data.append({
                            "label": label,
                            "flux": float(net_flux),
                            "flux_err": float(flux_err),
                            "gaia_gmag": float(gaia_gmag),
                        })

            # Flux relatif
            if flux_t1 > 0 and flux_comps > 0 and n_comps_valid > 0:
                rel_flux = flux_t1 / flux_comps
                term_a = (np.sqrt(err_t1_sq) / flux_t1) ** 2
                term_b = (np.sqrt(var_comps) / flux_comps) ** 2
                rel_err = rel_flux * np.sqrt(term_a + term_b)
            else:
                rel_flux = 0.0
                rel_err = 0.0

            row["tot_C_cnts"] = float(flux_comps)
            row["rel_flux_T1"] = float(rel_flux)
            row["rel_flux_err_T1"] = float(rel_err)

            # Magnitude T1 calibrée sur Gaia G via comparateurs
            mag_t1_g = np.nan
            rms_mag_t1 = np.nan
            zp_std = np.nan
            delta_auto = np.nan
            if flux_t1 > 0 and len(comp_mag_data) > 0:
                m_inst_t1 = -2.5 * np.log10(flux_t1)
                zp_vals = []
                for c in comp_mag_data:
                    if c["flux"] > 0:
                        m_inst_c = -2.5 * np.log10(c["flux"])
                        zp_vals.append(c["gaia_gmag"] - m_inst_c)
                if len(zp_vals) > 0:
                    zp_med = float(np.median(zp_vals))
                    zp_std = float(np.std(zp_vals)) if len(zp_vals) > 1 else 0.03
                    # Delta G-filtre auto estimé depuis comparateurs Gaia
                    delta_auto = float(zp_med)
                    self.photometry_delta_to_gaia_g = delta_auto
                    self.photometry_delta_var.set(delta_auto)
                    mag_t1_g = float(m_inst_t1 + delta_auto)
                    t1_mag_err = float(1.0857 * (np.sqrt(err_t1_sq) / flux_t1)) if flux_t1 > 0 else 0.0
                    rms_mag_t1 = float(np.sqrt(max(t1_mag_err**2 + zp_std**2, 0.0)))

            row["filter_used"] = self.photometry_filter_used or ""
            row["delta-to_G"] = float(self.photometry_delta_to_gaia_g)
            row["n_comp_gaia"] = int(len(comp_mag_data))
            row["mag_T1_G"] = float(mag_t1_g) if np.isfinite(mag_t1_g) else np.nan
            row["rmsMag_T1"] = float(rms_mag_t1) if np.isfinite(rms_mag_t1) else np.nan
            row["zp_std"] = float(zp_std) if np.isfinite(zp_std) else np.nan

            # Export CSV cumulatif (compilation progressive)
            results_dir = self._get_results_dir()
            if results_dir is None:
                raise RuntimeError("Dossier results indisponible.")
            asteroid_id = (self.asteroid_id_var.get() or "UNKNOWN").strip() or "UNKNOWN"
            out_csv = results_dir / f"{asteroid_id}_photometrie_image_par_image.csv"

            row_df = pd.DataFrame([row])
            if out_csv.exists():
                prev_df = pd.read_csv(out_csv)
                combined = pd.concat([prev_df, row_df], ignore_index=True)
                if "filename" in combined.columns:
                    combined = combined.drop_duplicates(subset=["filename"], keep="last")
                if "JD-UTC" in combined.columns:
                    combined = combined.sort_values("JD-UTC", na_position="last")
                # Conserver uniquement les colonnes demandées
                cols_keep = ["filename", "JD-UTC", "date_obs", "filter_used", "delta-to_G", "mag_T1_G", "rmsMag_T1"]
                for c in cols_keep:
                    if c not in combined.columns:
                        combined[c] = np.nan
                combined = combined[cols_keep]
                combined.to_csv(out_csv, index=False)
            else:
                cols_keep = ["filename", "JD-UTC", "date_obs", "filter_used", "delta-to_G", "mag_T1_G", "rmsMag_T1"]
                for c in cols_keep:
                    if c not in row_df.columns:
                        row_df[c] = np.nan
                row_df = row_df[cols_keep]
                row_df.to_csv(out_csv, index=False)

            logger.info(f"Photométrie image enregistrée : {out_csv}")
            messagebox.showinfo(
                "Photométrie image",
                (
                    f"Image traitée : {Path(self.current_image_path).name}\n"
                    f"rel_flux_T1 = {rel_flux:.6f}\n"
                    + (
                        f"mag_T1(G) = {mag_t1_g:.3f} ± {rms_mag_t1:.3f}\n"
                        f"Δ(G-filtre) auto = {delta_auto:+.3f} mag\n"
                        if np.isfinite(mag_t1_g)
                        else "mag_T1(G) indisponible (comparateurs Gaia insuffisants)\n"
                    )
                    + f"CSV compilé : {out_csv}"
                )
            )
        except Exception as e:
            logger.error(f"Erreur photométrie image courante : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Impossible de traiter l'image courante :\n{e}")
        finally:
            self._stop_process_log_capture("photométrie image courante")

    def _ensure_photometry_gaia_settings(self):
        """Demande le filtre utilisé. Le delta vers Gaia G est calculé automatiquement."""
        if self.photometry_settings_confirmed and self.photometry_filter_used:
            return True
        try:
            filt = simpledialog.askstring(
                "Filtre photométrie",
                "Filtre utilisé pour cette série (ex: G, V, r', clear) :",
                initialvalue=self.photometry_filter_used or "G",
                parent=self.frame.winfo_toplevel(),
            )
            if filt is None:
                return False
            filt = filt.strip()
            if not filt:
                filt = "G"
            self.photometry_filter_used = filt
            self.photometry_delta_to_gaia_g = 0.0
            self.photometry_filter_var.set(self.photometry_filter_used)
            self.photometry_delta_var.set(self.photometry_delta_to_gaia_g)
            self.photometry_settings_confirmed = True
            logger.info(
                f"Réglages photométrie: filtre={self.photometry_filter_used}, "
                "delta_to_G=auto (calculé depuis comparateurs Gaia)"
            )
            return True
        except Exception as e:
            logger.error(f"Paramétrage filtre photométrie impossible: {e}")
            return False

    def _apply_photometry_gaia_settings_from_ui(self):
        """Applique le filtre UI. Le delta vers Gaia G est recalculé automatiquement."""
        try:
            filt = (self.photometry_filter_var.get() or "").strip()
            if not filt:
                filt = "G"
            self.photometry_filter_used = filt
            self.photometry_settings_confirmed = True
            logger.info(
                f"Réglages Gaia appliqués (UI): filtre={self.photometry_filter_used}, delta_to_G=auto"
            )
            messagebox.showinfo(
                "Réglages Gaia",
                f"Réglages appliqués :\nFiltre = {self.photometry_filter_used}\n"
                "Δ(G-filtre) sera recalculé automatiquement depuis les comparateurs Gaia."
            )
        except Exception as e:
            messagebox.showerror("Erreur", f"Réglages Gaia invalides : {e}")

    def _reset_photometry_gaia_settings(self):
        """Réinitialise les réglages photométrie Gaia (retour à G, delta=0)."""
        self.photometry_filter_used = "G"
        self.photometry_delta_to_gaia_g = 0.0
        self.photometry_settings_confirmed = True
        self.photometry_filter_var.set("G")
        self.photometry_delta_var.set(0.0)
        logger.info("Réglages Gaia réinitialisés : filtre=G, delta_to_G=0.000")

    def _query_gaia_gmag_for_coord(self, coord, radius_arcsec=2.0):
        """Retourne la magnitude Gaia G la plus proche d'une coordonnée."""
        try:
            key = (round(float(coord.ra.deg), 6), round(float(coord.dec.deg), 6), round(float(radius_arcsec), 1))
            if key in self.photometry_gaia_mag_cache:
                return self.photometry_gaia_mag_cache[key]
            v = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag"], row_limit=20)
            res = v.query_region(coord, radius=radius_arcsec * u.arcsec, catalog="I/355/gaiadr3")
            if not res or len(res[0]) == 0:
                self.photometry_gaia_mag_cache[key] = None
                return None
            tab = res[0]
            ctab = SkyCoord(ra=tab["RA_ICRS"], dec=tab["DE_ICRS"], unit=u.deg)
            idx = int(np.argmin(coord.separation(ctab).arcsec))
            gmag = tab["Gmag"][idx] if "Gmag" in tab.colnames else None
            gmag_val = float(gmag) if gmag is not None and np.isfinite(gmag) else None
            self.photometry_gaia_mag_cache[key] = gmag_val
            return gmag_val
        except Exception:
            return None

    def _load_photometry_magnitude_map(self):
        """Charge les magnitudes T1 issues de la photométrie image-par-image."""
        if not self.directory:
            return {}
        asteroid_id = (self.asteroid_id_var.get() or "UNKNOWN").strip() or "UNKNOWN"
        results_dir = self._get_results_dir()
        if results_dir is None:
            return {}
        csv_path = results_dir / f"{asteroid_id}_photometrie_image_par_image.csv"
        if not csv_path.exists():
            return {}
        try:
            df = pd.read_csv(csv_path)
            out = {}
            for _, r in df.iterrows():
                fn = str(r.get("filename", "")).strip()
                if not fn:
                    continue
                mag = r.get("mag_T1_G", np.nan)
                rms_mag = r.get("rmsMag_T1", np.nan)
                if np.isfinite(mag):
                    log_snr = None
                    if np.isfinite(rms_mag) and rms_mag > 0:
                        # sigma_mag ≈ 1.0857 / SNR  =>  SNR ≈ 1.0857 / sigma_mag
                        snr = 1.0857 / float(rms_mag)
                        if np.isfinite(snr) and snr > 0:
                            log_snr = float(np.log10(snr))
                    out[fn] = {
                        "mag": float(mag),
                        "rms_mag": float(rms_mag) if np.isfinite(rms_mag) else None,
                        "log_snr": log_snr,
                    }
            if out:
                logger.info(f"Magnitudes photométriques chargées: {len(out)} images depuis {csv_path}")
            return out
        except Exception as e:
            logger.warning(f"Impossible de charger les magnitudes photométriques: {e}")
            return {}

    def _request_manual_apertures_for_image(self, filename, suggested):
        """Demande des apertures manuelles pour une image quand le FWHM est indisponible."""
        existing = self.photometry_manual_apertures.get(filename)
        if existing:
            return existing

        import threading
        done = threading.Event()
        out = {"value": None}

        def _prompt():
            try:
                proceed = messagebox.askyesno(
                    "FWHM indisponible",
                    f"Image: {filename}\n\nFWHM non estimable.\n"
                    "Voulez-vous saisir des rayons manuels (comète) pour cette image ?"
                )
                if not proceed:
                    done.set()
                    return

                init_rap = float((suggested or {}).get('r_ap', self.aperture_radius_var.get()))
                init_rin = float((suggested or {}).get('r_in', self.annulus_inner_var.get()))
                init_rout = float((suggested or {}).get('r_out', self.annulus_outer_var.get()))

                r_ap = simpledialog.askfloat("Aperture manuelle", f"{filename}\nRayon aperture (r_ap) [px]:", initialvalue=init_rap, parent=self.frame.winfo_toplevel(), minvalue=0.5)
                if r_ap is None:
                    done.set()
                    return
                r_in = simpledialog.askfloat("Aperture manuelle", f"{filename}\nRayon annulus interne (r_in) [px]:", initialvalue=max(init_rin, r_ap + 1.0), parent=self.frame.winfo_toplevel(), minvalue=r_ap + 0.1)
                if r_in is None:
                    done.set()
                    return
                r_out = simpledialog.askfloat("Aperture manuelle", f"{filename}\nRayon annulus externe (r_out) [px]:", initialvalue=max(init_rout, r_in + 1.0), parent=self.frame.winfo_toplevel(), minvalue=r_in + 0.1)
                if r_out is None:
                    done.set()
                    return

                values = {"r_ap": float(r_ap), "r_in": float(r_in), "r_out": float(r_out)}
                self.photometry_manual_apertures[filename] = values
                out["value"] = values
                logger.info(f"Apertures manuelles mémorisées pour {filename}: {values}")
            except Exception as e:
                logger.error(f"Saisie manuelle apertures impossible pour {filename}: {e}")
            finally:
                done.set()

        self.frame.after(0, _prompt)
        done.wait()
        return out["value"]

    def _open_lightcurve_modeling_window(self):
        """Ouvre la fenêtre Modélisation & Inversion : charger courbe, périodogramme, modèle, ajustement."""
        win = Toplevel(self.frame.winfo_toplevel())
        win.title("Modélisation & Inversion - Courbe de lumière astéroïde")
        win.geometry("900x720")
        win.transient(self.frame.winfo_toplevel())

        # Données chargées (dans la fenêtre)
        lc_data = {"time": None, "flux": None, "flux_err": None}
        period_best = {"value": None}
        fit_result = {"P": None, "t0": None, "A": None, "F0": None}

        # --- Fichier ---
        file_frame = ttk.LabelFrame(win, text="1. Fichier courbe de lumière", padding=8)
        file_frame.pack(fill=tk.X, padx=10, pady=5)
        default_path = ""
        if self.directory:
            for name in ("light_curve.txt", "results.csv"):
                p = Path(self.directory) / "photometrie" / name
                if p.exists():
                    default_path = str(p)
                    break
        path_var = tk.StringVar(value=default_path)
        ttk.Entry(file_frame, textvariable=path_var, width=70).pack(side=tk.LEFT, padx=2, fill=tk.X, expand=True)

        def load_file():
            path = path_var.get().strip()
            if not path or not Path(path).exists():
                messagebox.showwarning("Attention", "Fichier invalide ou absent.")
                return
            path = Path(path)
            try:
                if path.suffix.lower() == ".csv":
                    df = pd.read_csv(path)
                    time_col = next((c for c in ["JD-UTC", "JD_UTC", "time", "Time"] if c in df.columns), None)
                    flux_col = next((c for c in ["rel_flux_T1_fn", "flux", "rel_flux"] if c in df.columns), None)
                    err_col = next((c for c in ["rel_flux_err_T1", "flux_err", "err"] if c in df.columns), None)
                    if time_col is None or flux_col is None:
                        messagebox.showerror("Erreur", "CSV doit contenir JD-UTC (ou time) et rel_flux_T1_fn (ou flux).")
                        return
                    lc_data["time"] = np.array(df[time_col], dtype=float)
                    lc_data["flux"] = np.array(df[flux_col], dtype=float)
                    lc_data["flux_err"] = np.array(df[err_col], dtype=float) if err_col else None
                else:
                    # light_curve.txt : Time Relative_flux_fn relative_flux_fn_err
                    data = np.loadtxt(path)
                    if data.ndim == 1:
                        data = data.reshape(-1, 1)
                    lc_data["time"] = data[:, 0]
                    lc_data["flux"] = data[:, 1]
                    lc_data["flux_err"] = data[:, 2] if data.shape[1] >= 3 else None
                n = len(lc_data["time"])
                lc_data["time"] = np.asarray(lc_data["time"])
                lc_data["flux"] = np.asarray(lc_data["flux"])
                if lc_data["flux_err"] is not None:
                    lc_data["flux_err"] = np.asarray(lc_data["flux_err"])
                status_var.set(f"Chargé : {n} points")
                _refresh_plots()
            except Exception as e:
                logger.exception("Chargement courbe: %s", e)
                messagebox.showerror("Erreur", f"Impossible de charger : {e}")

        ttk.Button(file_frame, text="Charger", command=load_file).pack(side=tk.LEFT, padx=2)
        status_var = tk.StringVar(value="")
        ttk.Label(file_frame, textvariable=status_var, font=("Helvetica", 9)).pack(side=tk.LEFT, padx=5)

        # --- Périodogramme ---
        per_frame = ttk.LabelFrame(win, text="2. Périodogramme (Lomb-Scargle)", padding=8)
        per_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(per_frame, text="P min (j):").pack(side=tk.LEFT, padx=2)
        min_per_var = tk.StringVar(value="0.05")
        ttk.Entry(per_frame, textvariable=min_per_var, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Label(per_frame, text="P max (j):").pack(side=tk.LEFT, padx=2)
        max_per_var = tk.StringVar(value="20.0")
        ttk.Entry(per_frame, textvariable=max_per_var, width=8).pack(side=tk.LEFT, padx=2)

        def run_periodogram():
            if lc_data["time"] is None or lc_data["flux"] is None:
                messagebox.showwarning("Attention", "Chargez d'abord un fichier.")
                return
            try:
                min_per = float(min_per_var.get())
                max_per = float(max_per_var.get())
            except ValueError:
                messagebox.showwarning("Attention", "P min et P max doivent être des nombres.")
                return
            period, power, best = run_lomb_scargle(lc_data["time"], lc_data["flux"], min_period=min_per, max_period=max_per)
            period_best["value"] = best
            per_result_var.set(f"Période retenue : {best:.6f} j")
            ax_per.clear()
            ax_per.plot(period, power, "b-", lw=0.8)
            ax_per.axvline(best, color="red", ls="--", label=f"P = {best:.4f} j")
            ax_per.set_xlabel("Période (j)")
            ax_per.set_ylabel("Puissance")
            ax_per.set_title("Lomb-Scargle")
            ax_per.legend(loc="upper right", fontsize=8)
            ax_per.grid(True, alpha=0.3)
            P_var.set(f"{best:.6f}")
            fig_per.canvas.draw_idle()

        ttk.Button(per_frame, text="Lancer Lomb-Scargle", command=run_periodogram).pack(side=tk.LEFT, padx=5)
        per_result_var = tk.StringVar(value="")
        ttk.Label(per_frame, textvariable=per_result_var, font=("Helvetica", 9)).pack(side=tk.LEFT, padx=5)

        # --- Modèle ---
        param_frame = ttk.LabelFrame(win, text="3. Modèle F = F0·(1 + A·cos(2π(t-t0)/P))", padding=8)
        param_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(param_frame, text="P (j):").grid(row=0, column=0, sticky="w", padx=2)
        P_var = tk.StringVar(value="0.5")
        ttk.Entry(param_frame, textvariable=P_var, width=12).grid(row=0, column=1, padx=2)
        ttk.Label(param_frame, text="t0 (JD):").grid(row=0, column=2, sticky="w", padx=2)
        t0_var = tk.StringVar(value="")
        ttk.Entry(param_frame, textvariable=t0_var, width=14).grid(row=0, column=3, padx=2)
        ttk.Label(param_frame, text="A:").grid(row=1, column=0, sticky="w", padx=2)
        A_var = tk.StringVar(value="0.1")
        ttk.Entry(param_frame, textvariable=A_var, width=12).grid(row=1, column=1, padx=2)
        ttk.Label(param_frame, text="F0:").grid(row=1, column=2, sticky="w", padx=2)
        F0_var = tk.StringVar(value="1.0")
        ttk.Entry(param_frame, textvariable=F0_var, width=12).grid(row=1, column=3, padx=2)

        def compute_model():
            if lc_data["time"] is None:
                messagebox.showwarning("Attention", "Chargez d'abord un fichier.")
                return
            try:
                P = float(P_var.get())
                t0 = float(t0_var.get()) if t0_var.get().strip() else float(np.median(lc_data["time"]))
                A = float(A_var.get())
                F0 = float(F0_var.get())
            except ValueError:
                messagebox.showwarning("Attention", "Paramètres invalides.")
                return
            t0_var.set(f"{t0:.6f}")
            t_obs = lc_data["time"]
            f_mod = light_curve_model(t_obs, P, t0, A, F0)
            ax_lc.clear()
            ax_lc.errorbar(
                t_obs, lc_data["flux"],
                yerr=lc_data["flux_err"] if lc_data["flux_err"] is not None else None,
                fmt="o", markersize=4, alpha=0.7, label="Observations", capsize=2
            )
            ax_lc.plot(t_obs, f_mod, "r-", lw=2, label="Modèle")
            ax_lc.set_xlabel("Temps (JD)")
            ax_lc.set_ylabel("Flux")
            ax_lc.legend(loc="upper right", fontsize=8)
            ax_lc.grid(True, alpha=0.3)
            ax_lc.set_title("Courbe de lumière (temps)")
            # Phase
            phase = (t_obs - t0) % P / P
            ax_ph.clear()
            ax_ph.errorbar(
                phase, lc_data["flux"],
                yerr=lc_data["flux_err"] if lc_data["flux_err"] is not None else None,
                fmt="o", markersize=4, alpha=0.7, label="Observations", capsize=2
            )
            ax_ph.plot(phase, f_mod, "r-", lw=2, label="Modèle")
            ax_ph.set_xlabel("Phase")
            ax_ph.set_ylabel("Flux")
            ax_ph.legend(loc="upper right", fontsize=8)
            ax_ph.grid(True, alpha=0.3)
            ax_ph.set_title("Courbe de phase")
            fig_lc.canvas.draw_idle()

        def run_fit():
            if lc_data["time"] is None or lc_data["flux"] is None:
                messagebox.showwarning("Attention", "Chargez d'abord un fichier.")
                return
            try:
                P_init = float(P_var.get())
                t0_init = float(t0_var.get()) if t0_var.get().strip() else float(np.median(lc_data["time"]))
                A_init = float(A_var.get())
                F0_init = float(F0_var.get())
            except ValueError:
                P_init = period_best["value"] or 0.5
                t0_init = float(np.median(lc_data["time"]))
                A_init = 0.1
                F0_init = float(np.median(lc_data["flux"]))
            res = fit_light_curve(
                lc_data["time"], lc_data["flux"], lc_data["flux_err"],
                P_init=P_init, t0_init=t0_init, A_init=A_init, F0_init=F0_init,
            )
            fit_result.update(res)
            P_var.set(f"{res['P']:.6f}")
            t0_var.set(f"{res['t0']:.6f}")
            A_var.set(f"{res['A']:.4f}")
            F0_var.set(f"{res['F0']:.6f}")
            chi2_str = f"χ² = {res['chi2']:.2f}" + (f" ({res['n_dof']} ddl)" if res['n_dof'] else "")
            fit_status_var.set(f"Ajustement: {chi2_str} — " + ("OK" if res["success"] else res.get("message", "")))
            compute_model()

        ttk.Button(param_frame, text="Calculer le modèle", command=compute_model).grid(row=2, column=0, columnspan=2, padx=2, pady=4)
        ttk.Button(param_frame, text="4. Ajuster (inversion)", command=run_fit).grid(row=2, column=2, columnspan=2, padx=2, pady=4)
        fit_status_var = tk.StringVar(value="")
        ttk.Label(param_frame, textvariable=fit_status_var, font=("Helvetica", 8), foreground="blue").grid(row=3, column=0, columnspan=4, sticky="w", padx=2)

        # --- Figures ---
        fig_per = Figure(figsize=(4, 2))
        ax_per = fig_per.add_subplot(111)
        canvas_per = FigureCanvasTkAgg(fig_per, master=win)
        canvas_per.get_tk_widget().pack(fill=tk.BOTH, expand=False, padx=10, pady=5)
        fig_per.tight_layout()

        fig_lc = Figure(figsize=(8, 5))
        ax_lc = fig_lc.add_subplot(211)
        ax_ph = fig_lc.add_subplot(212)
        canvas_lc = FigureCanvasTkAgg(fig_lc, master=win)
        canvas_lc.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        fig_lc.tight_layout()

        def _refresh_plots():
            if lc_data["time"] is None:
                return
            t_med = float(np.median(lc_data["time"]))
            if not t0_var.get().strip():
                t0_var.set(f"{t_med:.6f}")
            compute_model()

        if default_path and Path(default_path).exists():
            win.after(100, load_file)

    def _open_damit_shape_window(self):
        """Ouvre un dialogue pour charger un fichier modèle 3D (.obj ou shape.txt DAMIT) et l'afficher."""
        path = filedialog.askopenfilename(
            title="Modèle de forme 3D (DAMIT)",
            initialdir=self.directory or "",
            filetypes=[
                ("Fichiers forme 3D", "*.obj;*.txt"),
                ("OBJ", "*.obj"),
                ("DAMIT shape.txt", "*.txt"),
                ("Tous", "*.*"),
            ],
        )
        if not path:
            return
        try:
            vertices, faces = load_shape(path)
        except Exception as e:
            logger.exception("Chargement modèle 3D: %s", e)
            messagebox.showerror("Erreur", f"Impossible de charger le modèle :\n{e}")
            return
        # Fenêtre 3D
        win = Toplevel(self.frame.winfo_toplevel())
        win.title(f"Modèle 3D — {Path(path).name}")
        win.geometry("700x600")
        win.transient(self.frame.winfo_toplevel())
        try:
            from mpl_toolkits.mplot3d import Axes3D
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        except ImportError:
            messagebox.showerror("Erreur", "matplotlib 3D (mplot3d) requis pour la visualisation.")
            return
        fig = Figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
        # Centrer et normaliser pour affichage
        v = vertices - vertices.mean(axis=0)
        r = np.sqrt((v ** 2).sum(axis=1)).max()
        if r > 0:
            v = v / r
        triangles = v[faces]  # (n_faces, 3, 3)
        poly = Poly3DCollection(
            [triangles[i] for i in range(len(triangles))],
            alpha=0.85,
            facecolor="lightsalmon",
            edgecolor="gray",
            linewidths=0.3,
        )
        ax.add_collection3d(poly)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_zlim(-1.2, 1.2)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(f"{Path(path).name} — {len(vertices)} sommets, {len(faces)} facettes")
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()