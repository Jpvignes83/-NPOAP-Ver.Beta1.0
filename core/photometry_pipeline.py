import config
import logging
import tkinter as tk
from tkinter import ttk, Toplevel, messagebox
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from astropy.visualization import ZScaleInterval, ImageNormalize, MinMaxInterval, HistEqStretch
from astroquery.gaia import Gaia
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats
from photutils.centroids import centroid_quadratic
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import Circle
import warnings

warnings.filterwarnings("ignore", category=AstropyWarning)
logger = logging.getLogger(__name__)

# -------------------------------------------------------------------
# 1. UTILITAIRES
# -------------------------------------------------------------------

def airmass(ra, dec, obstime):
    """Calcule l'airmass à partir des coordonnées RA/Dec et du temps d'observation."""
    try:
        # 1. Récupération Config
        obs = getattr(config, "OBSERVATORY", {})
        lat = obs.get("lat", -30.52)
        lon = obs.get("lon", -70.82)
        elev = obs.get("elev", 1710.0)

        # 2. Création du lieu
        location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=elev*u.m)

        # 3. Validation du temps
        if not isinstance(obstime, Time):
            obstime = Time(obstime)

        # 4. Calcul AltAz
        frame = AltAz(obstime=obstime, location=location)
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        altaz = coord.transform_to(frame)
        alt_deg = altaz.alt.deg

        # 5. Calcul Airmass
        if alt_deg <= 0:
            return 99.0
        
        zenith_angle_rad = np.radians(90.0 - alt_deg)
        am = 1.0 / np.cos(zenith_angle_rad)
        return float(am)
    except Exception as e:
        logger.warning(f"Erreur calcul airmass: {e}")
        return 1.0

def compute_zero_point(matched_table, airmass):
    try:
        zp_vals = matched_table['phot_g_mean_mag'] - matched_table['instrumental_mag']
        zp_clipped = sigma_clip(zp_vals, sigma=3, cenfunc='median')
        zp_median = np.median(zp_clipped)
        zp_corrected = zp_median + 0.1 * airmass
        zp_std = np.std(zp_clipped)
        return zp_corrected, zp_std
    except: return None, None

def match_sources_with_gaia(phot_df, search_radius_arcsec=2.0, gaia_mag_column="phot_g_mean_mag"):
    if not {"ra", "dec"}.issubset(phot_df.columns): return phot_df
    coords = SkyCoord(ra=phot_df["ra"].values * u.deg, dec=phot_df["dec"].values * u.deg)
    results = []
    for star in coords:
        try:
            r = Gaia.cone_search_async(star, radius=search_radius_arcsec * u.arcsec)
            res = r.get_results()
            if len(res) > 0:
                best = res[0]
                results.append({"ra_gaia": best["ra"], "dec_gaia": best["dec"], 
                                gaia_mag_column: best[gaia_mag_column], "source_id": best["source_id"]})
            else: results.append({"ra_gaia": np.nan, "dec_gaia": np.nan, gaia_mag_column: np.nan, "source_id": None})
        except: results.append({"ra_gaia": np.nan, "dec_gaia": np.nan, gaia_mag_column: np.nan, "source_id": None})
    return pd.concat([phot_df.reset_index(drop=True), pd.DataFrame(results)], axis=1)

def two_d_gaussian(x_array, y_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta):
    """Fonction gaussienne 2D pour le fit"""
    xt_array = x_array - model_x_mean
    yt_array = y_array - model_y_mean
    coss = np.cos(model_theta)
    sinn = np.sin(model_theta)
    
    return model_floor + model_norm * np.exp(-0.5 * (((-xt_array * sinn + yt_array * coss) / model_y_sigma) ** 2 + ((xt_array * coss + yt_array * sinn) / model_x_sigma) ** 2))

def estimate_fwhm_marginal(data, x, y, box_size=25):
    """
    Estime le FWHM d'une source en pixels en utilisant un fit gaussien 2D.
    Retourne None si l'estimation échoue ou est invalide.
    """
    try:
        x, y = int(x), int(y)
        half = box_size // 2
        if y-half < 0 or y+half+1 > data.shape[0] or x-half < 0 or x+half+1 > data.shape[1]: 
            return None, None, None
        
        sub = data[y-half:y+half+1, x-half:x+half+1]
        if sub.size == 0 or not np.isfinite(sub).any():
            return None, None, None
        
        # Créer les grilles de coordonnées pour le fit 2D
        x_indices = np.arange(sub.shape[1], dtype=float)
        y_indices = np.arange(sub.shape[0], dtype=float)
        x_grid, y_grid = np.meshgrid(x_indices, y_indices, indexing='xy')
        
        # Valeurs initiales pour le fit
        sub_flat = sub.flatten()
        bg_estimate = np.percentile(sub_flat, 10)  # Estimation du fond
        peak_value = np.max(sub_flat) - bg_estimate
        x_peak, y_peak = np.unravel_index(np.argmax(sub_flat), sub.shape)
        
        # Fonction wrapper pour curve_fit (prend coords comme premier argument)
        def gaussian_wrapper(coords, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta):
            x_array, y_array = coords
            result = two_d_gaussian(x_array, y_array, model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta)
            return result.flatten()
        
        # Paramètres initiaux : norm, floor, x_mean, y_mean, x_sigma, y_sigma, theta
        p0 = [
            peak_value,           # model_norm
            bg_estimate,          # model_floor
            float(x_peak),        # model_x_mean
            float(y_peak),        # model_y_mean
            3.0,                  # model_x_sigma (valeur initiale raisonnable)
            3.0,                  # model_y_sigma
            0.0                   # model_theta (angle de rotation, 0 = pas de rotation)
        ]
        
        # Bornes pour le fit
        bounds = (
            [0, -np.inf, 0, 0, 0.5, 0.5, -np.pi/4],  # Min
            [np.inf, np.inf, sub.shape[1], sub.shape[0], 15.0, 15.0, np.pi/4]  # Max
        )
        
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # Fit 2D gaussien
                popt, _ = curve_fit(
                    gaussian_wrapper,
                    (x_grid, y_grid),
                    sub.flatten(),
                    p0=p0,
                    bounds=bounds,
                    maxfev=5000
                )
                
                model_norm, model_floor, model_x_mean, model_y_mean, model_x_sigma, model_y_sigma, model_theta = popt
                
                # Validation : sigma doit être raisonnable
                if model_x_sigma < 0.5 or model_x_sigma > 15.0 or model_y_sigma < 0.5 or model_y_sigma > 15.0:
                    return None, None, None
                
                # Calcul du FWHM depuis les sigmas 2D
                # Pour une gaussienne 2D, on utilise la moyenne géométrique des sigmas
                # FWHM = 2.355 * sqrt(sigma_x * sigma_y) pour l'aire efficace
                # Ou moyenne : FWHM = 2.355 * (sigma_x + sigma_y) / 2
                # Ici on utilise la moyenne pour être cohérent avec l'ancien code
                sigma_effective = (model_x_sigma + model_y_sigma) / 2.0
                fwhm_val = 2.355 * sigma_effective
                
                # Validation finale : FWHM doit être raisonnable (entre 1 et 30 pixels)
                if fwhm_val < 1.0 or fwhm_val > 30.0:
                    return None, None, None
                
                # Calculer les projections marginales pour show_diagnostics_windows
                # À partir du fit 2D, on peut calculer les projections en intégrant sur une direction
                # Pour simplifier, on utilise les sigmas du fit 2D pour créer des projections gaussiennes
                row_sum = np.sum(sub, axis=1)
                col_sum = np.sum(sub, axis=0)
                
                # Créer des arrays d'indices pour les projections
                r = np.arange(len(row_sum))
                c = np.arange(len(col_sum))
                
                # Paramètres des gaussiennes 1D pour les projections (approximation)
                # On utilise les sigmas du fit 2D
                # Pour une gaussienne 2D, les projections sont aussi des gaussiennes
                row_peak = np.argmax(row_sum)
                col_peak = np.argmax(col_sum)
                row_max = np.max(row_sum)
                col_max = np.max(col_sum)
                row_bg = np.median(row_sum)
                col_bg = np.median(col_sum)
                
                # Paramètres pour les projections (approximation depuis le fit 2D)
                # sigma_y pour la projection Y (row), sigma_x pour la projection X (col)
                p_opt_row = [row_max - row_bg, float(row_peak), model_y_sigma, row_bg]
                p_opt_col = [col_max - col_bg, float(col_peak), model_x_sigma, col_bg]
                
                return fwhm_val, (r, row_sum, p_opt_row), (c, col_sum, p_opt_col)
                
        except (RuntimeError, ValueError, TypeError) as e:
            # Échec du fit gaussien
            logger.warning(f"Fit gaussien échoué pour ({x}, {y}): {e}")
            return None, None, None
    except Exception as e:
        logger.warning(f"Erreur dans estimate_fwhm_marginal pour ({x}, {y}): {e}")
        return None, None, None

def refine_centroid(data, x0, y0, box_size=25):
    xi, yi = int(x0), int(y0)
    half = box_size // 2
    if yi-half < 0 or xi-half < 0: return x0, y0
    cut = data[yi-half:yi+half+1, xi-half:xi+half+1]
    if cut.size == 0 or not np.isfinite(cut).all(): return x0, y0
    try:
        dx, dy = centroid_quadratic(cut)
        if not np.isfinite(dx): return x0, y0
        return xi - half + dx, yi - half + dy
    except Exception: 
        return x0, y0

def show_diagnostics_windows(data, x_t1, y_t1, fwhm_info, radii=None, parent=None):
    """Affiche les profils FWHM et Radial. Retourne la liste des fenêtres créées."""
    created_windows = []
    fwhm_val, row_data, col_data = fwhm_info
    if not (fwhm_val and row_data and col_data): return []
    if radii: r_ap, r_in, r_out = radii
    else: r_ap, r_in, r_out = 1.4*fwhm_val, 2.6*fwhm_val, 3.4*fwhm_val 

    gauss = lambda x, a, x0, s, off: a * np.exp(-(x-x0)**2/(2*s**2)) + off
    def draw_limits(ax, center, fw, rap, rin, rout):
        ax.axvline(center, color='k', ls='--', alpha=0.3)
        ax.axvline(center - fw/2, color='b', ls=':', label=f'FWHM={fw:.2f}')
        ax.axvline(center + fw/2, color='b', ls=':')
        ax.axvline(center + rap, color='g', ls='-', label=f'Ap={rap:.1f}')
        ax.axvline(center - rap, color='g', ls='-')
        ax.axvline(center + rin, color='orange', ls='--', label=f'In={rin:.1f}')
        ax.axvline(center + rout, color='r', ls='--', label=f'Out={rout:.1f}')
        ax.legend(fontsize=8)

    try:
        win_f = Toplevel(parent) if parent else Toplevel()
        win_f.title(f"Profils FWHM (T1) ~ {fwhm_val:.2f} px")
        created_windows.append(win_f)
       
        fig = Figure(figsize=(11, 5.5))
        ax1, ax2 = fig.subplots(1, 2)

        # === Coupe Y (verticale) ===
        r, y_dat, p_opt = row_data
        r_shifted = r - p_opt[1]
        p_opt_shifted = list(p_opt)
        p_opt_shifted[1] = 0

        ax1.plot(r_shifted, y_dat, 'o-', markersize=4, alpha=0.7, label='Données Y')
        
        # Domaine élargi pour prolonger la gaussienne (±6 sigma par défaut)
        pad = 6.0
        x_fit_y = np.linspace(r_shifted.min() - pad * abs(p_opt[2]),
                              r_shifted.max() + pad * abs(p_opt[2]), 800)
        ax1.plot(x_fit_y, gauss(x_fit_y, *p_opt_shifted), 'r-', lw=1.5,
                 label=f'Fit Gaussien (FWHM = {2.355*abs(p_opt[2]):.3f} px)')

        draw_limits(ax1, 0, 2.355*abs(p_opt[2]), r_ap, r_in, r_out)
        ax1.set_title("Coupe Y (Centrée)")
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)

        # === Coupe X (horizontale) ===
        c, x_dat, p_opt2 = col_data
        c_shifted = c - p_opt2[1]
        p_opt2_shifted = list(p_opt2)
        p_opt2_shifted[1] = 0

        ax2.plot(c_shifted, x_dat, 'o-', markersize=4, alpha=0.7, label='Données X')
        
        x_fit_x = np.linspace(c_shifted.min() - pad * abs(p_opt2[2]),
                              c_shifted.max() + pad * abs(p_opt2[2]), 800)
        ax2.plot(x_fit_x, gauss(x_fit_x, *p_opt2_shifted), 'r-', lw=1.5,
                 label=f'Fit Gaussien (FWHM = {2.355*abs(p_opt2[2]):.3f} px)')

        draw_limits(ax2, 0, 2.355*abs(p_opt2[2]), r_ap, r_in, r_out)
        ax2.set_title("Coupe X (Centrée)")
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

        FigureCanvasTkAgg(fig, win_f).get_tk_widget().pack(fill="both", expand=True)

    except Exception as e:
        logger.error(f"Erreur Profils FWHM: {e}")


    # =============================================================================
    # Deuxième fenêtre : Profil Radial PSF
    # =============================================================================
    try:
        win_p = Toplevel(parent) if parent else Toplevel()
        win_p.title("Profil Radial PSF (T1)")
        created_windows.append(win_p)
       
        fig = Figure(figsize=(7, 5.5))
        ax = fig.add_subplot(111)

        r_vals, profile, popt = col_data
        r_shifted = r_vals - popt[1]
        popt_shifted = list(popt)
        popt_shifted[1] = 0

        # Données expérimentales
        ax.plot(r_shifted, profile, 'bo-', markersize=5, alpha=0.8, label='Données radiales')

        # Gaussienne prolongée très loin (±8 sigma pour bien voir la ligne de base)
        pad_radial = 8.0
        r_fit = np.linspace(-pad_radial * abs(popt[2]), 
                            +pad_radial * abs(popt[2]), 1000)
        ax.plot(r_fit, gauss(r_fit, *popt_shifted), 'r-', lw=2.8,
                label=f'Fit Gaussien\nFWHM = {2.355*abs(popt[2]):.3f} px')

        # Lignes de référence
        ax.axvline(0, color='k', ls='--', alpha=0.4, lw=1.2)
        ax.axvline(r_ap, color='g', ls='-', lw=2, label=f'Apperture = {r_ap:.2f}')
        ax.axvline(r_in, color='orange', ls='--', lw=1.2, label=f'Inner = {r_in:.2f}')
        ax.axvline(r_out, color='red', ls='--', lw=1.2, label=f'Outer = {r_out:.2f}')

        # Optionnel : ligne horizontale de l'offset
        ax.axhline(popt_shifted[3], color='gray', ls=':', alpha=0.7, lw=1)

        ax.set_xlim(-r_out*1.3, r_out*1.3)
        ax.set_title("Profil Radial Centré (PSF)", fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        FigureCanvasTkAgg(fig, win_p).get_tk_widget().pack(fill="both", expand=True)

    except Exception as e:
        logger.error(f"Erreur Profil Radial: {e}")

    return created_windows

# -------------------------------------------------------------------
# 3. INTERFACE INTERACTIVE (Choix Ouvertures)
# -------------------------------------------------------------------
def launch_photometry_aperture(fits_path, target_data, comp_coords_data, variable_flags=None, on_finish=None):
    
    # --- 1. Extraction des Données ---
    fwhm_t1_provided = None
    if isinstance(target_data, dict):
        target_coord_in = target_data['coord']
        fwhm_t1_provided = target_data.get('fwhm')  # FWHM mesuré par clic droit (optionnel)
    else:
        target_coord_in = target_data

    comp_coords_in = []
    if comp_coords_data:
        for c in comp_coords_data:
            if isinstance(c, dict):
                comp_coords_in.append(c['coord'])
            else:
                comp_coords_in.append(c)

    variable_flags = variable_flags or {}

    try:
        with fits.open(fits_path) as hdul:
            data = hdul[0].data.astype(float)
            header = hdul[0].header
            wcs = WCS(header)
    except Exception as e:
        return messagebox.showerror("Erreur", f"Lecture FITS: {e}")

    root = tk.Toplevel()
    root.title(f"Sélection Photométrique - {Path(fits_path).name}")
    root.geometry("1250x800")
    
    paned = tk.PanedWindow(root, orient=tk.HORIZONTAL)
    paned.pack(fill=tk.BOTH, expand=True)
    
    frame_left = tk.Frame(paned, width=400, bg="#f0f0f0")
    frame_right = tk.Frame(paned, bg="black")
    paned.add(frame_left)
    paned.add(frame_right)

    tk.Label(frame_left, text="Cibles & Comparateurs", font=("Arial", 12, "bold"), bg="#f0f0f0").pack(pady=(10,5))
    tk.Label(frame_left, text="Coordonnées = Centroïde sur l'image", font=("Arial", 9), bg="#f0f0f0", fg="blue").pack(pady=2)
    
    frame_bot = tk.Frame(frame_left, bg="#e0e0e0", bd=1, relief="sunken")
    frame_bot.pack(side="bottom", fill="x")

    canvas_list = tk.Canvas(frame_left, bg="#f0f0f0", highlightthickness=0)
    scroll = tk.Scrollbar(frame_left, command=canvas_list.yview)
    frame_inner = tk.Frame(canvas_list, bg="#f0f0f0")
    
    canvas_list.configure(yscrollcommand=scroll.set)
    scroll.pack(side="right", fill="y")
    canvas_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    
    win_can = canvas_list.create_window((0,0), window=frame_inner, anchor="nw")
    
    def update_scroll(e):
        canvas_list.configure(scrollregion=canvas_list.bbox("all"))
        canvas_list.itemconfigure(win_can, width=canvas_list.winfo_width())
        
    frame_inner.bind("<Configure>", update_scroll)
    canvas_list.bind("<Configure>", update_scroll)

    h = tk.Frame(frame_inner, bg="#d9d9d9")
    h.pack(fill="x", pady=2)
    headers = [("Act",3), ("Nom",5), ("RA (deg)", 9), ("Dec (deg)", 9), ("Ap",4), ("In",4), ("Out",4)]
    for t, w in headers:
        tk.Label(h, text=t, width=w, bg="#d9d9d9", font=("Arial",8,"bold")).pack(side="left", padx=1)

    fig = Figure(figsize=(6,6), dpi=100)
    ax = fig.add_subplot(111, projection=wcs)
    ax.format_coord = lambda x, y: ""

    valid = data[np.isfinite(data)]
    vm, vM = ZScaleInterval().get_limits(valid) if valid.size > 0 else (0,1)
    im = ax.imshow(data, origin='lower', cmap='gray', vmin=vm, vmax=vM)
    fig.colorbar(im, ax=ax, shrink=0.8)
    
    canvas_mpl = FigureCanvasTkAgg(fig, master=frame_right)
    canvas_mpl.draw()
    canvas_widget = canvas_mpl.get_tk_widget()
    canvas_widget.pack(fill="both", expand=True)    
    
    toolbar = NavigationToolbar2Tk(canvas_mpl, frame_right)
    toolbar.update()
    toolbar.pack(fill="x", side="bottom")
    
    # Changer le curseur pour le rendre plus visible en mode zoom
    # Utiliser "pencil" qui est plus visible
    def update_cursor(event):
        mode = toolbar.mode if hasattr(toolbar, 'mode') else ""
        if mode and mode != "" and 'zoom' in mode.lower():
            canvas_widget.config(cursor="pencil")
        else:
            canvas_widget.config(cursor="")
    
    canvas_mpl.mpl_connect('motion_notify_event', update_cursor)

    stars_data = []

    def _format_aperture_field(val):
        if val is None or (isinstance(val, float) and not np.isfinite(val)):
            return "8.0"
        x = float(val)
        s = f"{x:.3f}".rstrip("0").rstrip(".")
        return s if s else "0"

    def _parse_aperture_field(svar):
        try:
            return float(str(svar.get()).strip().replace(",", "."))
        except (TypeError, ValueError):
            return None

    def update_aperture_visibility(aperture_data):
        state = aperture_data['active'].get()
        for p in aperture_data['patches']:
            p.set_visible(state)
        canvas_mpl.draw_idle()

    def add_ui(label, x, y, ra, ri, ro, color, coord):
        c_ap = Circle((x, y), ra, ec=color, ls="-", fill=False)
        c_in = Circle((x, y), ri, ec=color, ls=":", fill=False)
        c_out = Circle((x, y), ro, ec=color, ls=":", fill=False)
        t_label = ax.text(x, y + 15, label, color=color, fontweight="bold", ha="center")
        ax.add_patch(c_ap)
        ax.add_patch(c_in)
        ax.add_patch(c_out)

        row = tk.Frame(frame_inner, bg="#f0f0f0")
        row.pack(fill="x", pady=1)
        act = tk.BooleanVar(value=True)
        vap = tk.StringVar(value=_format_aperture_field(ra))
        vin = tk.StringVar(value=_format_aperture_field(ri))
        vout = tk.StringVar(value=_format_aperture_field(ro))

        star_entry = {
            "label": label,
            "coord": coord,
            "x": x,
            "y": y,
            "vars": (vap, vin, vout),
            "active": act,
            "patches": (c_ap, c_in, c_out, t_label),
        }
        tk.Checkbutton(row, variable=act, bg="#f0f0f0").pack(side="left", padx=1)
        act.trace_add("write", lambda *a: update_aperture_visibility(star_entry))
        lbl = tk.Label(
            row,
            text=label,
            fg=color,
            font=("Consolas", 9, "bold"),
            width=5,
            bg="#f0f0f0",
            cursor="hand2",
        )
        lbl.pack(side="left")
        tk.Label(row, text=f"{coord.ra.deg:.5f}", width=9, bg="#eef").pack(side="left", padx=1)
        tk.Label(row, text=f"{coord.dec.deg:.5f}", width=9, bg="#eef").pack(side="left", padx=1)

        def refresh_circles(*_):
            try:
                a = _parse_aperture_field(vap)
                b = _parse_aperture_field(vin)
                c = _parse_aperture_field(vout)
                if a is None or b is None or c is None:
                    return
                if not (a > 0 and b > a and c > b):
                    return
                c_ap.set_radius(a)
                c_in.set_radius(b)
                c_out.set_radius(c)
                canvas_mpl.draw_idle()
            except Exception:
                pass

        def on_entry_commit(event=None):
            refresh_circles()

        for sv in (vap, vin, vout):
            ent = tk.Entry(row, textvariable=sv, width=7, justify="right", takefocus=True)
            ent.pack(side="left", padx=1)
            ent.bind("<FocusOut>", on_entry_commit)
            ent.bind("<Return>", on_entry_commit)

        stars_data.append(star_entry)
        lbl.bind(
            "<Button-1>",
            lambda e: (
                ax.set_xlim(x - 100, x + 100),
                ax.set_ylim(y - 100, y + 100),
                canvas_mpl.draw_idle(),
            ),
        )

        if label == "T1":

            def on_right_click(e):
                fwhm_info = estimate_fwhm_marginal(data, x, y)
                if fwhm_info and fwhm_info[0] is not None:
                    r_ap_curr = _parse_aperture_field(vap)
                    r_in_curr = _parse_aperture_field(vin)
                    r_out_curr = _parse_aperture_field(vout)
                    if r_ap_curr is None or r_in_curr is None or r_out_curr is None:
                        r_ap_curr = float(ra)
                        r_in_curr = float(ri)
                        r_out_curr = float(ro)
                    show_diagnostics_windows(
                        data, x, y, fwhm_info, radii=(r_ap_curr, r_in_curr, r_out_curr)
                    )

            lbl.bind("<Button-3>", on_right_click)

    def process_coord(wx, wy, label, fwhm_override=None):
        x, y = refine_centroid(data, wx, wy, box_size=25)
        if not np.isfinite(x): return None
        
        # Si FWHM fourni en override (pour T1 depuis clic droit), l'utiliser
        if fwhm_override is not None and np.isfinite(fwhm_override) and 1.0 <= fwhm_override <= 30.0:
            fwhm = fwhm_override
            logger.info(f"[{label}] Utilisation FWHM fourni : {fwhm:.2f} px")
        else:
            # Sinon, estimer le FWHM
            fwhm_res = estimate_fwhm_marginal(data, x, y)
            
            # Validation du FWHM : doit être dans une plage raisonnable (1-30 pixels)
            if fwhm_res and fwhm_res[0] and np.isfinite(fwhm_res[0]):
                fwhm = fwhm_res[0]
                # Limiter le FWHM à une plage raisonnable
                if fwhm < 1.0 or fwhm > 30.0:
                    logger.warning(f"FWHM invalide pour {label}: {fwhm:.2f} px, utilisation valeur par défaut")
                    fwhm = 4.0
            else:
                fwhm = 4.0  # Valeur par défaut
        
        # Calcul des apertures basées sur le FWHM
        r_ap = max(2.0, min(round(2.0*fwhm, 1), 50.0))  # Ouverture = FWHM × 2
        r_in = max(4.0, min(round(3.6*fwhm, 1), 60.0))  # Coefficient 3.6
        r_out = max(6.0, min(round(5.0*fwhm, 1), 70.0))  # Coefficient 5.0
        
        logger.info(f"[{label}] FWHM={fwhm:.2f} px → Apertures: r_ap={r_ap:.1f}, r_in={r_in:.1f}, r_out={r_out:.1f}")
        
        # Vérification de cohérence : r_in > r_ap et r_out > r_in
        if r_in <= r_ap:
            r_in = r_ap + 2.0
            logger.warning(f"[{label}] Ajustement r_in: {r_in:.1f} (r_in <= r_ap)")
        if r_out <= r_in:
            r_out = r_in + 2.0
            logger.warning(f"[{label}] Ajustement r_out: {r_out:.1f} (r_out <= r_in)")
        
        col = "red" if label == "T1" else ("orange" if variable_flags.get(label) else "cyan")
        real_coord = wcs.pixel_to_world(x, y)
        add_ui(label, x, y, r_ap, r_in, r_out, col, real_coord)

    def on_click(e):
        if toolbar.mode != "" and toolbar.mode is not None: return
        if e.button != 1 or e.xdata is None: return
        t1 = any(s['label'] == 'T1' for s in stars_data)
        if not t1: process_coord(e.xdata, e.ydata, "T1")
        else: 
            idx = len([s for s in stars_data if s['label'].startswith('C')]) + 1
            process_coord(e.xdata, e.ydata, f"C{idx}")
        canvas_mpl.draw_idle()

    canvas_mpl.mpl_connect("button_press_event", on_click)

    def validate():
        if not any(s["label"] == "T1" for s in stars_data):
            return messagebox.showwarning("Erreur", "Pas de T1")
        sels = []
        for s in stars_data:
            if not s["active"].get():
                continue
            ra = _parse_aperture_field(s["vars"][0])
            ri = _parse_aperture_field(s["vars"][1])
            ro = _parse_aperture_field(s["vars"][2])
            if ra is None or ri is None or ro is None:
                return messagebox.showwarning(
                    "Rayons invalides",
                    f"{s['label']} : saisissez des nombres valides pour Ap / In / Out.",
                )
            if not (ra > 0 and ri > ra and ro > ri):
                return messagebox.showwarning(
                    "Rayons incohérents",
                    f"{s['label']} : il faut Ap > 0, In > Ap, Out > In (en pixels).",
                )
            sels.append(
                {
                    "label": s["label"],
                    "coord": s["coord"],
                    "r_ap": ra,
                    "r_in": ri,
                    "r_out": ro,
                }
            )
        root.destroy()
        if on_finish:
            on_finish(sels)

    tk.Button(frame_bot, text="VALIDER & ANALYSER", bg="#ccffcc", font=("Arial",11,"bold"), command=validate).pack(side="left", padx=5, pady=5)
    
    if target_coord_in:
        tx, ty = wcs.world_to_pixel(target_coord_in)
        if 0 <= tx < data.shape[1] and 0 <= ty < data.shape[0]: 
            process_coord(tx, ty, "T1", fwhm_override=fwhm_t1_provided)

    if comp_coords_in:
        for i, cc in enumerate(comp_coords_in, 1):
            try:
                cx, cy = wcs.world_to_pixel(cc)
                if 0 <= cx < data.shape[1] and 0 <= cy < data.shape[0]: process_coord(cx, cy, f"C{i}")
            except: pass
    canvas_mpl.draw()

# -------------------------------------------------------------------
# 3. PIPELINE BATCH
# -------------------------------------------------------------------

def process_photometry_series(folder, target_coord, comps, ref_image, selections=None, min_snr=5.0, variable_flags=None, progress_callback=None, ephemeris_data=None, astrometry_positions=None, manual_t1_anchors=None):
    """
    Exécute la photométrie en série.
    
    Parameters
    ----------
    ephemeris_data : Table, optional
        Table d'éphémérides avec colonnes 'datetime_jd', 'RA', 'DEC' pour interpoler T1 à chaque image.
        Si fourni, la position de T1 sera interpolée depuis les éphémérides pour chaque image.
    astrometry_positions : dict, optional
        Dictionnaire {nom_fichier: {'ra_deg': float, 'dec_deg': float, 'jd_utc': float}}.
        Si fourni, les positions astrométriques mesurées (ADES) seront utilisées en priorité pour suivre T1.
    """
    
    folder = Path(folder)
    fits_files = sorted(folder.glob("*.fits"))
    total = len(fits_files)

    if not fits_files:
        raise FileNotFoundError("Aucun fichier FITS trouvé.")

    logger.info(f"Début Batch sur {total} fichiers.")

    results = []

    def _has_signal(data, x, y, box_size=25, snr_thresh=5.0):
        try:
            xi, yi = int(x), int(y)
            half = box_size // 2
            if yi-half < 0 or xi-half < 0 or yi+half+1 > data.shape[0] or xi+half+1 > data.shape[1]:
                return False
            sub = data[yi-half:yi+half+1, xi-half:xi+half+1]
            if sub.size == 0 or not np.isfinite(sub).any():
                return False
            med = np.nanmedian(sub)
            std = np.nanstd(sub)
            if not np.isfinite(std) or std <= 0:
                return False
            peak = np.nanmax(sub)
            return (peak - med) >= (snr_thresh * std)
        except Exception:
            return False

    def _find_peak(data, x, y, box_size=65):
        try:
            xi, yi = int(x), int(y)
            half = box_size // 2
            if yi-half < 0 or xi-half < 0 or yi+half+1 > data.shape[0] or xi-half < 0 or xi+half+1 > data.shape[1]:
                return None
            sub = data[yi-half:yi+half+1, xi-half:xi+half+1]
            if sub.size == 0 or not np.isfinite(sub).any():
                return None
            max_idx = np.nanargmax(sub)
            dy, dx = np.unravel_index(max_idx, sub.shape)
            return xi - half + dx, yi - half + dy
        except Exception:
            return None

    # Définition des étoiles à mesurer
    stars_to_process = selections or []
    if not stars_to_process:
        if target_coord:
            stars_to_process.append({'label': 'T1', 'coord': target_coord, 'r_ap': 8, 'r_in': 12, 'r_out': 18})
        if comps:
            for i, c in enumerate(comps):
                stars_to_process.append({'label': f'C{i+1}', 'coord': c, 'r_ap': 8, 'r_in': 12, 'r_out': 18})

    comp_labels = [s['label'] for s in stars_to_process if s.get('label', '').startswith('C')]
    comp_valid_flux_counts = {lab: 0 for lab in comp_labels}
    comp_invalid_flux_counts = {lab: 0 for lab in comp_labels}
    bad_t1_images = 0
    bad_comp_images = 0
    comp_images_total = 0

    # Référence WCS
    try:
        with fits.open(ref_image) as hdul:
            wcs_ref = WCS(hdul[0].header)
    except Exception:
        wcs_ref = None

    def _get_jd_from_header(header):
        jd_val = header.get("JD-UTC", 0.0)
        if jd_val == 0.0:
            try:
                date_obs = header.get("DATE-OBS") or header.get("DATE-LOC")
                if date_obs:
                    jd_val = Time(date_obs).jd
            except Exception:
                pass
        return jd_val

    manual_anchors_mode = manual_t1_anchors is not None
    t1_anchor_first = manual_t1_anchors.get("first") if manual_anchors_mode else None
    t1_anchor_last = manual_t1_anchors.get("last") if manual_anchors_mode else None
    if manual_anchors_mode:
        logger.info(
            "Ancres T1 manuelles utilisées. "
            f"first={t1_anchor_first}, last={t1_anchor_last}"
        )

    def _interp_manual_t1_coord(jd_val):
        if not (t1_anchor_first and t1_anchor_last and jd_val):
            return None
        jd1 = t1_anchor_first.get("jd")
        jd2 = t1_anchor_last.get("jd")
        if jd1 is None or jd2 is None or jd2 == jd1:
            return None
        frac = (jd_val - jd1) / (jd2 - jd1)
        try:
            c1 = SkyCoord(t1_anchor_first["ra_deg"] * u.deg, t1_anchor_first["dec_deg"] * u.deg, frame="icrs")
            c2 = SkyCoord(t1_anchor_last["ra_deg"] * u.deg, t1_anchor_last["dec_deg"] * u.deg, frame="icrs")
            v1 = c1.cartesian.xyz.value
            v2 = c2.cartesian.xyz.value
            vx = v1[0] + frac * (v2[0] - v1[0])
            vy = v1[1] + frac * (v2[1] - v1[1])
            vz = v1[2] + frac * (v2[2] - v1[2])
            coord_cart = SkyCoord(x=vx, y=vy, z=vz, representation_type='cartesian', frame='icrs')
            coord_sph = coord_cart.represent_as('spherical')
            return SkyCoord(ra=coord_sph.lon, dec=coord_sph.lat, frame='icrs')
        except Exception:
            return None

    def _get_t1_coord_for_file(fpath, jd_val, base_coord):
        if base_coord is None:
            return None

        if manual_anchors_mode:
            manual_coord = _interp_manual_t1_coord(jd_val)
            return manual_coord or base_coord
        coord_updated = None

        if astrometry_positions is not None:
            filename = Path(fpath).name
            if filename in astrometry_positions:
                pos_data = astrometry_positions[filename]
                ra_deg = pos_data.get('ra_deg')
                dec_deg = pos_data.get('dec_deg')
                if ra_deg is not None and dec_deg is not None:
                    try:
                        coord_updated = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
                    except Exception:
                        coord_updated = None

        if coord_updated is None and ephemeris_data is not None and jd_val > 0:
            try:
                eph_jds = ephemeris_data['datetime_jd']
                ra_values = ephemeris_data['RA']
                dec_values = ephemeris_data['DEC']
                ra_diff = np.diff(ra_values)
                has_ra_wrap = np.any(np.abs(ra_diff) > 180.0)
                if has_ra_wrap:
                    coords = SkyCoord(ra=ra_values, dec=dec_values, unit=u.deg)
                    xyz = coords.cartesian.xyz.value
                    x_i = np.interp(jd_val, eph_jds, xyz[0])
                    y_i = np.interp(jd_val, eph_jds, xyz[1])
                    z_i = np.interp(jd_val, eph_jds, xyz[2])
                    coord_cart = SkyCoord(x=x_i, y=y_i, z=z_i, representation_type='cartesian', frame='icrs')
                    coord_sph = coord_cart.represent_as('spherical')
                    coord_updated = SkyCoord(ra=coord_sph.lon, dec=coord_sph.lat, frame='icrs')
                else:
                    ra_interp = float(np.interp(jd_val, eph_jds, ra_values))
                    dec_interp = float(np.interp(jd_val, eph_jds, dec_values))
                    coord_updated = SkyCoord(ra=ra_interp * u.deg, dec=dec_interp * u.deg, frame='icrs')
            except Exception:
                coord_updated = None

        return coord_updated or base_coord

    t1_def = next((s for s in stars_to_process if s['label'] == 'T1'), None)

    # Traitement image par image
    for i, fpath in enumerate(fits_files):
        try:
            with fits.open(fpath) as hdul:
                header = hdul[0].header
                data = hdul[0].data.astype(float)
                wcs = WCS(header)

                jd = header.get("JD-UTC", 0.0)
                if jd == 0.0:
                    try:
                        date_obs = header.get("DATE-OBS") or header.get("DATE-LOC")
                        if date_obs:
                            jd = Time(date_obs).jd
                    except Exception:
                        pass

                # Calcul de l'airmass : d'abord essayer depuis le header, sinon calculer
                am_header = header.get("AIRMASS", None)
                am = None
                
                # Si la valeur du header est négative ou invalide, calculer l'airmass
                if am_header is not None:
                    try:
                        am_val = float(am_header)
                        if am_val > 0 and am_val < 100:  # Valeur raisonnable
                            am = am_val
                    except (ValueError, TypeError):
                        pass
                
                # Si pas de valeur valide du header, calculer depuis les coordonnées
                if am is None and jd > 0:
                    try:
                        # Récupérer les coordonnées de la cible
                        t1_def = next((s for s in stars_to_process if s['label'] == 'T1'), None)
                        if t1_def:
                            t1_coord = t1_def['coord']
                            # Pour T1, utiliser positions astrométriques ADES en priorité, sinon éphémérides
                            t1_coord_updated = None
                            
                            if manual_anchors_mode:
                                t1_coord_updated = _interp_manual_t1_coord(jd)
                            else:
                                # 1. Priorité aux positions astrométriques ADES
                                if astrometry_positions is not None:
                                    filename = Path(fpath).name
                                    if filename in astrometry_positions:
                                        pos_data = astrometry_positions[filename]
                                        ra_deg = pos_data.get('ra_deg')
                                        dec_deg = pos_data.get('dec_deg')
                                        if ra_deg is not None and dec_deg is not None:
                                            try:
                                                t1_coord_updated = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
                                            except Exception:
                                                pass
                                
                                # 2. Fallback sur éphémérides si positions astrométriques non disponibles
                                if t1_coord_updated is None and ephemeris_data is not None and jd > 0:
                                    try:
                                        eph_jds = ephemeris_data['datetime_jd']
                                        ra_values = ephemeris_data['RA']
                                        dec_values = ephemeris_data['DEC']
                                        ra_diff = np.diff(ra_values)
                                        has_ra_wrap = np.any(np.abs(ra_diff) > 180.0)
                                        if has_ra_wrap:
                                            coords = SkyCoord(ra=ra_values, dec=dec_values, unit=u.deg)
                                            xyz = coords.cartesian.xyz.value
                                            x_i = np.interp(jd, eph_jds, xyz[0])
                                            y_i = np.interp(jd, eph_jds, xyz[1])
                                            z_i = np.interp(jd, eph_jds, xyz[2])
                                            coord_cart = SkyCoord(x=x_i, y=y_i, z=z_i, 
                                                                 representation_type='cartesian', frame='icrs')
                                            coord_sph = coord_cart.represent_as('spherical')
                                            t1_coord_updated = SkyCoord(ra=coord_sph.lon, dec=coord_sph.lat, frame='icrs')
                                        else:
                                            ra_interp = float(np.interp(jd, eph_jds, ra_values))
                                            dec_interp = float(np.interp(jd, eph_jds, dec_values))
                                            t1_coord_updated = SkyCoord(ra=ra_interp * u.deg, dec=dec_interp * u.deg, frame='icrs')
                                    except Exception:
                                        pass  # Utiliser coordonnées fixes
                            
                            # Utiliser la coordonnée mise à jour si disponible
                            if t1_coord_updated is not None:
                                t1_coord = t1_coord_updated
                            
                            # Calculer l'airmass depuis les coordonnées
                            if hasattr(t1_coord, 'ra') and hasattr(t1_coord, 'dec'):
                                am = airmass(t1_coord.ra.deg, t1_coord.dec.deg, Time(jd, format='jd'))
                    except Exception as e:
                        logger.debug(f"Erreur calcul airmass pour {fpath.name}: {e}")
                
                # Fallback si aucun calcul n'a fonctionné
                if am is None:
                    am = 1.0

                gain = header.get("GAIN", 1.0)
                if gain <= 0:
                    gain = 1.0

                row = {
                    "slice": i+1,
                    "JD-UTC": jd,
                    "AIRMASS": float(am)
                }

                flux_t1_net = err_t1_sq = flux_comps_total = var_comps = 0
                comps_total_this_image = 0
                comps_valid_this_image = 0
                fwhm_list = []
                
                # Estimer le FWHM sur T1 pour cette image (avant de traiter les étoiles)
                # Ce FWHM sera utilisé pour calculer les apertures dynamiquement
                image_fwhm = None
                t1_px_refined = None
                t1_pred_px = None
                t1_pred_py = None
                centroid_max_shift_px = 15.0
                t1_def = next((s for s in stars_to_process if s['label'] == 'T1'), None)
                if t1_def:
                    try:
                        t1_coord = t1_def['coord']
                        # Pour T1, utiliser ancres manuelles si disponibles, sinon ADES/éphémérides
                        t1_coord_updated = None
                        if manual_anchors_mode:
                            t1_coord_updated = _interp_manual_t1_coord(jd)
                        else:
                            # 1. Priorité aux positions astrométriques ADES
                            if astrometry_positions is not None:
                                filename = Path(fpath).name
                                if filename in astrometry_positions:
                                    pos_data = astrometry_positions[filename]
                                    ra_deg = pos_data.get('ra_deg')
                                    dec_deg = pos_data.get('dec_deg')
                                    if ra_deg is not None and dec_deg is not None:
                                        try:
                                            t1_coord_updated = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
                                        except Exception:
                                            pass
                            
                            # 2. Fallback sur éphémérides si positions astrométriques non disponibles
                            if t1_coord_updated is None and ephemeris_data is not None and jd > 0:
                                try:
                                    eph_jds = ephemeris_data['datetime_jd']
                                    ra_values = ephemeris_data['RA']
                                    dec_values = ephemeris_data['DEC']
                                    ra_diff = np.diff(ra_values)
                                    has_ra_wrap = np.any(np.abs(ra_diff) > 180.0)
                                    if has_ra_wrap:
                                        coords = SkyCoord(ra=ra_values, dec=dec_values, unit=u.deg)
                                        xyz = coords.cartesian.xyz.value
                                        x_i = np.interp(jd, eph_jds, xyz[0])
                                        y_i = np.interp(jd, eph_jds, xyz[1])
                                        z_i = np.interp(jd, eph_jds, xyz[2])
                                        coord_cart = SkyCoord(x=x_i, y=y_i, z=z_i, 
                                                             representation_type='cartesian', frame='icrs')
                                        coord_sph = coord_cart.represent_as('spherical')
                                        t1_coord_updated = SkyCoord(ra=coord_sph.lon, dec=coord_sph.lat, frame='icrs')
                                    else:
                                        ra_interp = float(np.interp(jd, eph_jds, ra_values))
                                        dec_interp = float(np.interp(jd, eph_jds, dec_values))
                                        t1_coord_updated = SkyCoord(ra=ra_interp * u.deg, dec=dec_interp * u.deg, frame='icrs')
                                except Exception:
                                    pass  # Utiliser coordonnées fixes
                        
                        # Utiliser la coordonnée mise à jour si disponible
                        if t1_coord_updated is not None:
                            t1_coord = t1_coord_updated
                        
                        px_t1, py_t1 = wcs.world_to_pixel(t1_coord)
                        ny, nx = data.shape
                        if 0 <= px_t1 < nx and 0 <= py_t1 < ny:
                            t1_pred_px, t1_pred_py = px_t1, py_t1
                            px_t1, py_t1 = refine_centroid(data, px_t1, py_t1, box_size=45)
                            if np.isfinite(px_t1) and np.isfinite(py_t1) and 0 <= px_t1 < nx and 0 <= py_t1 < ny:
                                if np.isfinite(px_t1) and np.isfinite(py_t1) and 0 <= px_t1 < nx and 0 <= py_t1 < ny:
                                    t1_px_refined = (px_t1, py_t1)
                                elif t1_pred_px is not None and t1_pred_py is not None:
                                    t1_px_refined = (t1_pred_px, t1_pred_py)
                                fwhm_x, fwhm_y = t1_px_refined if t1_px_refined else (px_t1, py_t1)
                                fwhm_result = estimate_fwhm_marginal(data, fwhm_x, fwhm_y)
                                if fwhm_result and fwhm_result[0] and np.isfinite(fwhm_result[0]):
                                    image_fwhm = fwhm_result[0]
                    except Exception:
                        pass
                
                # Si FWHM non estimé, utiliser valeurs par défaut
                if image_fwhm is None or image_fwhm <= 0:
                    image_fwhm = 4.0  # FWHM par défaut
                
                # Calcul des apertures basées sur le FWHM
                r_ap_default = 1.4 * image_fwhm
                r_in_default = 3.6 * image_fwhm
                r_out_default = 5.0 * image_fwhm
                logger.info(f"[BATCH] Image {Path(fpath).name}: FWHM={image_fwhm:.2f} px → Apertures: r_ap={r_ap_default:.1f}, r_in={r_in_default:.1f}, r_out={r_out_default:.1f}")

                for star_def in stars_to_process:
                    label = star_def['label']
                    coord = star_def['coord']
                    # Utiliser les apertures calculées depuis le FWHM (pas les valeurs par défaut du star_def)
                    r_ap = r_ap_default
                    r_in = r_in_default
                    r_out = r_out_default

                    # Pour T1, utiliser ancres manuelles si disponibles, sinon ADES/éphémérides
                    if label == 'T1':
                        coord_updated = None
                        if manual_anchors_mode:
                            coord_updated = _interp_manual_t1_coord(jd)
                        else:
                            # 1. Priorité aux positions astrométriques ADES
                            if astrometry_positions is not None:
                                filename = Path(fpath).name
                                if filename in astrometry_positions:
                                    pos_data = astrometry_positions[filename]
                                    ra_deg = pos_data.get('ra_deg')
                                    dec_deg = pos_data.get('dec_deg')
                                    if ra_deg is not None and dec_deg is not None:
                                        try:
                                            coord_updated = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
                                            logger.debug(f"T1 position depuis astrométrie ADES pour {filename}: RA={ra_deg:.6f}°, Dec={dec_deg:.6f}°")
                                        except Exception as e:
                                            logger.debug(f"Erreur création coord depuis astrométrie ADES: {e}")
                            
                            # 2. Fallback sur éphémérides si positions astrométriques non disponibles
                            if coord_updated is None and ephemeris_data is not None and jd > 0:
                                try:
                                    # Interpolation de la position T1 depuis les éphémérides
                                    eph_jds = ephemeris_data['datetime_jd']
                                    ra_values = ephemeris_data['RA']
                                    dec_values = ephemeris_data['DEC']
                                    
                                    # Vérifier s'il y a une discontinuité RA
                                    ra_diff = np.diff(ra_values)
                                    has_ra_wrap = np.any(np.abs(ra_diff) > 180.0)
                                    
                                    if has_ra_wrap:
                                        # Interpolation cartésienne
                                        coords = SkyCoord(ra=ra_values, dec=dec_values, unit=u.deg)
                                        xyz = coords.cartesian.xyz.value
                                        x_i = np.interp(jd, eph_jds, xyz[0])
                                        y_i = np.interp(jd, eph_jds, xyz[1])
                                        z_i = np.interp(jd, eph_jds, xyz[2])
                                        coord_cart = SkyCoord(x=x_i, y=y_i, z=z_i, 
                                                             representation_type='cartesian', frame='icrs')
                                        coord_sph = coord_cart.represent_as('spherical')
                                        coord_updated = SkyCoord(ra=coord_sph.lon, dec=coord_sph.lat, frame='icrs')
                                    else:
                                        # Interpolation directe RA/Dec
                                        ra_interp = float(np.interp(jd, eph_jds, ra_values))
                                        dec_interp = float(np.interp(jd, eph_jds, dec_values))
                                        coord_updated = SkyCoord(ra=ra_interp * u.deg, dec=dec_interp * u.deg, frame='icrs')
                                    
                                    logger.debug(f"T1 interpolé depuis éphémérides pour {Path(fpath).name}: JD={jd:.10f}, RA={coord_updated.ra.deg:.6f}°, Dec={coord_updated.dec.deg:.6f}°")
                                except Exception as e:
                                    logger.warning(f"Impossible d'interpoler T1 depuis éphémérides pour {Path(fpath).name}: {e}, utilisation coordonnées fixes")
                        
                        # Utiliser la coordonnée mise à jour si disponible
                        if coord_updated is not None:
                            coord = coord_updated

                    # Vérifier que la source est dans l'image avant la photométrie
                    try:
                        if label == 'T1' and t1_px_refined is not None:
                            px, py = t1_px_refined
                            logger.debug(f"T1 recentré par centroid pour {Path(fpath).name}: px={px:.2f}, py={py:.2f}")
                        else:
                            px, py = wcs.world_to_pixel(coord)
                        # Vérifier que les coordonnées pixel sont dans l'image
                        ny, nx = data.shape
                        if not (0 <= px < nx and 0 <= py < ny):
                            logger.debug(f"{label} hors champ sur {Path(fpath).name}: px={px:.1f}, py={py:.1f}, image={nx}x{ny}")
                            continue
                        
                        # Pour T1, utiliser un centroïde élargi (box_size=35) pour mieux capturer le PSF
                        # Pour les comparateurs, utiliser la taille par défaut (box_size=25)
                        centroid_box_size = 45 if label == 'T1' else 25
                        if not (label == 'T1' and t1_px_refined is not None):
                            if label == 'T1' and t1_pred_px is not None and t1_pred_py is not None:
                                px, py = refine_centroid(data, px, py, box_size=centroid_box_size)
                            else:
                                px, py = refine_centroid(data, px, py, box_size=centroid_box_size)
                            if label == 'T1' and t1_pred_px is not None and t1_pred_py is not None:
                                shift = np.hypot(px - t1_pred_px, py - t1_pred_py)
                                if shift > centroid_max_shift_px:
                                    logger.warning(
                                        f"T1 centroid rejeté (shift={shift:.1f}px) pour {Path(fpath).name}. "
                                        "Utilisation position prédite."
                                    )
                                    px, py = t1_pred_px, t1_pred_py
                        if not np.isfinite(px) or not np.isfinite(py): 
                            raise ValueError("NaN centroid")
                        # Vérifier à nouveau après centroidage
                        if not (0 <= px < nx and 0 <= py < ny):
                            logger.debug(f"{label} hors champ après centroidage sur {Path(fpath).name}")
                            continue
                    except Exception:
                        if wcs_ref:
                            try:
                                px, py = wcs_ref.world_to_pixel(coord)
                                ny, nx = data.shape
                                if not (0 <= px < nx and 0 <= py < ny):
                                    logger.debug(f"{label} hors champ (WCS ref) sur {Path(fpath).name}")
                                    continue
                                centroid_box_size = 35 if label == 'T1' else 25
                                px, py = refine_centroid(data, px, py, box_size=centroid_box_size)
                                if not np.isfinite(px) or not np.isfinite(py):
                                    continue
                                if not (0 <= px < nx and 0 <= py < ny):
                                    continue
                            except Exception:
                                logger.debug(f"{label} non trouvé sur {Path(fpath).name}")
                                continue
                        else:
                            logger.debug(f"{label} non trouvé (pas de WCS ref) sur {Path(fpath).name}")
                            continue

                    ap = CircularAperture((px, py), r=r_ap)
                    an = CircularAnnulus((px, py), r_in=r_in, r_out=r_out)
                    stats = ApertureStats(data, an)
                    bkg_mean = stats.mean
                    bkg_std = stats.std
                    bg_area = stats.sum_aper_area.value

                    phot = aperture_photometry(data, ap)
                    raw_flux = phot['aperture_sum'][0]
                    ap_area = ap.area
                    net_flux = raw_flux - (bkg_mean * ap_area)

                    # Erreur photométrique
                    term_source = max(net_flux, 0) / gain
                    term_sky = ap_area * (bkg_std ** 2)
                    term_bkg_err = (ap_area ** 2) * (bkg_std ** 2) / bg_area
                    flux_err = np.sqrt(term_source + term_sky + term_bkg_err)

                    # Enregistrement des résultats
                    row[f"Flux_{label}"] = net_flux
                    row[f"Err_{label}"] = flux_err
                    row[f"Sky_{label}"] = bkg_mean
                    row[f"X_{label}"] = px
                    row[f"Y_{label}"] = py

                    # Estimation FWHM (si possible)
                    try:
                        fwhm_val, *_ = estimate_fwhm_marginal(data, px, py)
                        if np.isfinite(fwhm_val) and fwhm_val > 0:
                            fwhm_list.append(fwhm_val)
                    except Exception:
                        pass

                    if label == "T1":
                        flux_t1_net = net_flux
                        err_t1_sq = flux_err**2
                        row["Sky/Pixel_T1"] = bkg_mean
                        row["Source-Sky_T1"] = net_flux
                        t1_invalid = not (np.isfinite(net_flux) and net_flux > 0)
                        logger.info(
                            f"[T1] {Path(fpath).name} px=({px:.2f},{py:.2f}) "
                            f"net_flux={net_flux:.3f} invalid={t1_invalid}"
                        )
                    elif label.startswith("C"):
                        flux_comps_total += net_flux
                        var_comps += flux_err**2
                        comps_total_this_image += 1
                        if np.isfinite(net_flux) and net_flux > 0:
                            comps_valid_this_image += 1
                            if label in comp_valid_flux_counts:
                                comp_valid_flux_counts[label] += 1
                        else:
                            if label in comp_invalid_flux_counts:
                                comp_invalid_flux_counts[label] += 1

                # Moyenne des FWHM
                row["FWHM_Mean"] = np.mean(fwhm_list) if fwhm_list else np.nan

                # Calcul flux relatif
                row["tot_C_cnts"] = flux_comps_total
                if flux_comps_total > 0 and flux_t1_net > 0:
                    rel_flux = flux_t1_net / flux_comps_total
                    term_A = (np.sqrt(err_t1_sq) / flux_t1_net)**2
                    term_B = (np.sqrt(var_comps) / flux_comps_total)**2
                    err_final = rel_flux * np.sqrt(term_A + term_B)
                else:
                    rel_flux = 0
                    err_final = 0

                row["rel_flux_T1"] = rel_flux
                row["rel_flux_err_T1"] = err_final

                results.append(row)

                if comps_total_this_image > 0:
                    comp_images_total += 1
                    if comps_valid_this_image == 0:
                        bad_comp_images += 1
                if not (np.isfinite(flux_t1_net) and flux_t1_net > 0):
                    bad_t1_images += 1

        except Exception as e:
            logger.warning(f"Skip {fpath.name} → {e}")
            continue

        if progress_callback:
            progress_callback(int((i + 1) / total * 100))

    # --- EXPORT CSV FINAL ---
    if results:
        df = pd.DataFrame(results)

        # Nettoyage colonnes comparateurs inutiles
        cols_to_remove = [col for col in df.columns if any(col.startswith(prefix) for prefix in ["Flux_C", "Err_C", "Sky_C", "X_C", "Y_C"])]
        df.drop(columns=cols_to_remove, inplace=True, errors="ignore")
        total_images = len(results)
        if comp_labels:
            comps_without_valid_flux = [
                f"{lab} ({comp_valid_flux_counts.get(lab, 0)}/{total_images})"
                for lab in comp_labels
                if comp_valid_flux_counts.get(lab, 0) == 0
            ]
            logger.info(
                "Suppression colonnes comparateurs (nettoyage CSV). "
                f"Images T1 invalides: {bad_t1_images}/{total_images}; "
                f"Images comps sans flux valide: {bad_comp_images}/{comp_images_total if comp_images_total else total_images}; "
                f"Comps sans flux valide: {', '.join(comps_without_valid_flux) if comps_without_valid_flux else 'aucun'}"
            )
        logger.info(f"Colonnes supprimées (comparateurs) : {cols_to_remove}")

        # Normalisation flux relatif
        median_rel = df["rel_flux_T1"].median()
        if median_rel > 0:
            df["rel_flux_T1_fn"] = df["rel_flux_T1"] / median_rel
            df["rel_flux_err_T1"] = df["rel_flux_err_T1"] / median_rel
            # Calcul des résidus : différence avec la médiane normalisée (1.0)
            # Les résidus indiquent l'écart par rapport à la valeur médiane
            df["rel_flux_T1_fn_residual"] = df["rel_flux_T1_fn"] - 1.0
        else:
            df["rel_flux_T1_fn"] = 0
            df["rel_flux_err_T1"] = 0
            df["rel_flux_T1_fn_residual"] = 0

        # Ordre des colonnes préférentiel
        cols_priority = [
            'slice', 'JD-UTC', 'AIRMASS', 'FWHM_Mean',
            'rel_flux_T1_fn', 'rel_flux_T1_fn_residual',
            'rel_flux_T1', 'rel_flux_err_T1'
        ]
        df = df[cols_priority + [c for c in df.columns if c not in cols_priority]]

        # Export final
        phot_dir = folder / "photometrie"
        phot_dir.mkdir(exist_ok=True)
        results_path = phot_dir / "results.csv"
        try:
            df.to_csv(results_path, index=False)
            logger.info(f"Résultats : {results_path}")
        except PermissionError:
            fallback_path = phot_dir / f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            df.to_csv(fallback_path, index=False)
            logger.warning(f"Permission refusée pour {results_path}. Résultats écrits dans {fallback_path}")

    return results



class PhotometryPipeline:
    def __init__(self, logger=None): self.logger = logger or logging.getLogger(__name__)
    @staticmethod
    def airmass(ra, dec, obstime): return airmass(ra, dec, obstime)
    @staticmethod
    def compute_zero_point(mt, am): return compute_zero_point(mt, am)
    @staticmethod
    def match_sources_with_gaia(df, sr=2.0, col="phot_g_mean_mag"): return match_sources_with_gaia(df, sr, col)
    def launch_photometry_aperture(self, fits_path, target_coord, comp_coords, variable_flags=None, on_finish=None):
        return launch_photometry_aperture(fits_path, target_coord, comp_coords, variable_flags, on_finish)
    def process_photometry_series(self, folder, target_coord, comps, 
                                  ref_image, selections=None, 
                                  min_snr=5.0, variable_flags=None, progress_callback=None, ephemeris_data=None):
        return process_photometry_series(
            folder=folder,
            target_coord=target_coord,
            comps=comps,
            ref_image=ref_image,
            selections=selections,
            min_snr=min_snr,
            variable_flags=variable_flags,
            progress_callback=progress_callback,
            ephemeris_data=ephemeris_data
        )