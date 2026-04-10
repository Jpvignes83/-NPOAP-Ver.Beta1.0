"""
Module d'Astrométrie Optimisé (Fast Astrometry)

Ce module optimise le processus d'astrométrie dans l'onglet Photométrie Astéroïdes,
inspiré des techniques utilisées par Tycho Tracker.

La méthodologie zero-aperture (extrapolation des positions/incertitudes à aperture 0)
est alignée sur le projet Zero-Aperture-Astrometry :
https://github.com/bensharkey/Zero-Aperture-Astrometry

Optimisations:
- Cache local persistant des catalogues Gaia
- Matching rapide avec KD-Tree (10-100x plus rapide)
- Sélection intelligente d'apertures (6 valeurs par image)
- Support GPU optionnel avec CuPy
- Fit linéaire pondéré (RA/Dec vs aperture) pour extrapolation zero-aperture
"""

import logging
import pickle
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Optional, Tuple, Dict

import numpy as np
from scipy.spatial import KDTree
from scipy.optimize import curve_fit

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import sigma_clipped_stats
from astropy.utils.exceptions import AstropyWarning

from photutils.detection import DAOStarFinder
from photutils.centroids import centroid_2dg
from astroquery.vizier import Vizier

logger = logging.getLogger(__name__)

# Tentative d'import CuPy pour GPU
try:
    # Supprimer l'avertissement CuPy concernant CUDA_PATH si CUDA n'est pas configuré
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*CUDA path could not be detected.*", category=UserWarning)
        import cupy as cp
    HAS_GPU = True
    # Vérifier si CUDA est réellement disponible
    try:
        cp.cuda.is_available()
        logger.debug("CuPy installé et CUDA disponible")
    except Exception:
        HAS_GPU = False
        logger.debug("CuPy installé mais CUDA non disponible - utilisation CPU")
except Exception as e:
    # ImportError si absent ; AttributeError/RuntimeError possibles si site-packages
    # incohérent (ex. metadata None dans _detect_duplicate_installation).
    HAS_GPU = False
    cp = None
    logger.debug("CuPy indisponible (%s) — astrométrie en CPU uniquement", type(e).__name__)
warnings.filterwarnings("ignore", category=AstropyWarning)
Vizier.server = 'http://vizier.cfa.harvard.edu/'


# =============================================================================
# 0. COHÉRENCE SIP / CTYPE (Éviter l'avertissement Astropy)
# =============================================================================
# Réf. https://docs.astropy.org/en/stable/wcs/note_sip.html
# Si le header contient des coefficients SIP mais que CTYPE n'a pas le suffixe "-SIP",
# Astropy émet un INFO. On corrige le header avant de construire le WCS.


def fix_header_sip_ctype(header: fits.Header) -> fits.Header:
    """
    Retourne une copie du header avec CTYPE1/CTYPE2 cohérents avec la présence de SIP.

    Si des coefficients SIP sont présents (A_ORDER, B_ORDER, etc.) et que CTYPE
    n'a pas le suffixe "-SIP", ajoute "-SIP" pour supprimer l'avertissement
    Astropy et garantir des coordonnées cohérentes.

    Parameters
    ----------
    header : fits.Header
        Header FITS d'entrée (non modifié).

    Returns
    -------
    fits.Header
        Nouveau header (copie) avec CTYPE mis à jour si nécessaire.
    """
    out = header.copy()
    # Présence de coefficients SIP (convention FITS)
    has_sip = 'A_ORDER' in out or 'B_ORDER' in out
    if not has_sip:
        return out

    def add_sip_suffix(ctype: str) -> str:
        if not ctype or str(ctype).strip().endswith('-SIP'):
            return ctype
        c = str(ctype).strip()
        if '-TAN' in c:
            return c.replace('-TAN', '-TAN-SIP')
        if '-SIN' in c:
            return c.replace('-SIN', '-SIN-SIP')
        return c + '-SIP'

    if 'CTYPE1' in out:
        out['CTYPE1'] = add_sip_suffix(out['CTYPE1'])
    if 'CTYPE2' in out:
        out['CTYPE2'] = add_sip_suffix(out['CTYPE2'])
    logger.debug("Header: CTYPE1/CTYPE2 alignés avec les coefficients SIP (-SIP ajouté)")
    return out


@contextmanager
def open_fits_with_fixed_wcs(path, mode='readonly', **kwargs):
    """
    Ouvre un fichier FITS et corrige le header (SIP/CTYPE) à la source.

    À utiliser partout où on charge un FITS dont le header servira à créer un WCS,
    pour éviter l'avertissement Astropy "Inconsistent SIP distortion information".
    Le header de la première extension est remplacé par une copie avec CTYPE
    cohérent (suffixe -SIP si coefficients SIP présents).

    Parameters
    ----------
    path : str ou Path
        Chemin du fichier FITS.
    mode : str
        Mode d'ouverture ('readonly', 'update', etc.), passé à fits.open.
    **kwargs
        Arguments additionnels pour fits.open.

    Yields
    ------
    astropy.io.fits.HDUList
        HDU list ; hdul[0].header est déjà corrigé pour SIP/CTYPE.
    """
    with fits.open(path, mode=mode, **kwargs) as hdul:
        if len(hdul) > 0 and hasattr(hdul[0], 'header'):
            hdul[0].header = fix_header_sip_ctype(hdul[0].header.copy())
        yield hdul


# =============================================================================
# 1. CACHE GAIA PERSISTANT
# =============================================================================

class GaiaCache:
    """Gestion du cache local persistant pour les données Gaia."""
    
    CACHE_DIR = Path.home() / ".npoap" / "gaia_cache"
    
    def __init__(self, cache_dir=None):
        self.cache_dir = cache_dir or self.CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Cache Gaia initialisé : {self.cache_dir}")
    
    def _cache_key(self, coord: SkyCoord, radius: u.Quantity) -> str:
        """Génère une clé de cache basée sur la région du ciel (grille 0.5°)."""
        ra_idx = int(coord.ra.deg / 0.5)
        dec_idx = int(coord.dec.deg / 0.5)
        radius_deg = radius.to(u.deg).value
        # Arrondir le rayon à 0.1° près pour éviter les variations mineures
        radius_rounded = round(radius_deg, 1)
        return f"gaia_tile_{ra_idx}_{dec_idx}_r{radius_rounded:.1f}.pkl"
    
    def _find_cached_with_larger_radius(self, coord: SkyCoord, radius: u.Quantity) -> Optional[Table]:
        """Cherche dans le cache une entrée avec un rayon plus grand qui pourrait contenir les données."""
        ra_idx = int(coord.ra.deg / 0.5)
        dec_idx = int(coord.dec.deg / 0.5)
        radius_deg = radius.to(u.deg).value
        radius_rounded = round(radius_deg, 1)
        
        # Chercher des fichiers avec le même RA/Dec mais rayon >= radius
        pattern = f"gaia_tile_{ra_idx}_{dec_idx}_r*.pkl"
        for cache_file in self.cache_dir.glob(pattern):
            try:
                # Extraire le rayon du nom de fichier
                radius_str = cache_file.stem.split('_r')[1]
                cached_radius = float(radius_str)
                if cached_radius >= radius_rounded:
                    # Le rayon du cache est suffisant, utiliser ces données
                    with open(cache_file, 'rb') as f:
                        table = pickle.load(f)
                    logger.debug(f"Cache hit (rayon plus grand) : {cache_file.name} ({len(table)} étoiles, rayon={cached_radius:.1f}° >= {radius_rounded:.1f}°)")
                    return table
            except (ValueError, IndexError, Exception) as e:
                logger.debug(f"Erreur parsing cache file {cache_file.name}: {e}")
                continue
        return None
    
    def _cache_path(self, coord: SkyCoord, radius: u.Quantity) -> Path:
        """Retourne le chemin du fichier de cache."""
        return self.cache_dir / self._cache_key(coord, radius)
    
    def get_cached(self, coord: SkyCoord, radius: u.Quantity) -> Optional[Table]:
        """Récupère les données Gaia depuis le cache."""
        cache_path = self._cache_path(coord, radius)
        if cache_path.exists():
            try:
                with open(cache_path, 'rb') as f:
                    table = pickle.load(f)
                logger.info(f"Cache hit : {cache_path.name} ({len(table)} étoiles)")
                return table
            except Exception as e:
                logger.warning(f"Erreur lecture cache {cache_path}: {e}")
        
        # Si pas trouvé exactement, chercher avec un rayon plus grand
        larger_cache = self._find_cached_with_larger_radius(coord, radius)
        if larger_cache is not None:
            return larger_cache
        
        logger.debug(f"Cache miss pour {coord} avec rayon {radius}")
        return None
    
    def save_cache(self, coord: SkyCoord, radius: u.Quantity, table: Table):
        """Sauvegarde les données Gaia dans le cache."""
        cache_path = self._cache_path(coord, radius)
        try:
            with open(cache_path, 'wb') as f:
                pickle.dump(table, f)
            logger.debug(f"Cache sauvegardé : {cache_path.name} ({len(table)} étoiles)")
        except Exception as e:
            logger.error(f"Erreur sauvegarde cache {cache_path}: {e}")
    
    def build_kdtree(self, gaia_table: Table) -> Optional[KDTree]:
        """Construit un KD-tree à partir d'une table Gaia."""
        if len(gaia_table) == 0:
            return None
        coords = SkyCoord(ra=gaia_table['RA_ICRS'], dec=gaia_table['DE_ICRS'], unit=u.deg)
        xyz = coords.cartesian.xyz.value.T
        return KDTree(xyz)
    
    def get_cache_info(self) -> dict:
        """Retourne des informations sur le cache."""
        cache_files = list(self.cache_dir.glob("*.pkl"))
        total_size = sum(f.stat().st_size for f in cache_files)
        return {
            "count": len(cache_files),
            "total_size_mb": total_size / (1024 * 1024),
            "cache_dir": str(self.cache_dir)
        }
    
    def clear_cache(self):
        """Supprime tous les fichiers de cache."""
        cache_files = list(self.cache_dir.glob("*.pkl"))
        for f in cache_files:
            f.unlink()
        logger.info(f"Cache nettoyé : {len(cache_files)} fichiers supprimés")


# =============================================================================
# 2. MATCHING RAPIDE AVEC KD-TREE
# =============================================================================

class FastStarMatcher:
    """Matching utilisant un KD-Tree spatial pour accélération 10-100x."""
    
    def __init__(self, gaia_table: Table, gaia_cache: Optional[GaiaCache] = None):
        """
        Initialise le matcher avec une table Gaia.
        
        Parameters
        ----------
        gaia_table : Table
            Table Gaia avec colonnes RA_ICRS, DE_ICRS
        gaia_cache : GaiaCache, optional
            Cache pour réutiliser le KD-tree si disponible
        """
        self.gaia_table = gaia_table
        self.gaia_sky = SkyCoord(
            ra=gaia_table['RA_ICRS'], 
            dec=gaia_table['DE_ICRS'], 
            unit=u.deg
        )
        # Construction du KD-tree
        xyz = self.gaia_sky.cartesian.xyz.value.T
        self.kdtree = KDTree(xyz)
        logger.debug(f"KD-tree construit pour {len(gaia_table)} étoiles")
    
    def match(self, detected_coords: SkyCoord, dist_limit: u.Quantity = 5.0*u.arcsec) -> Tuple[np.ndarray, np.ndarray]:
        """
        Match les sources détectées avec le catalogue Gaia.
        
        Parameters
        ----------
        detected_coords : SkyCoord
            Coordonnées des sources détectées
        dist_limit : Quantity
            Distance maximale de matching en arcsec
        
        Returns
        -------
        indices : np.ndarray
            Indices des correspondances dans gaia_table (-1 si pas de match)
        distances : np.ndarray
            Distances angulaires en arcsec
        """
        if len(detected_coords) == 0 or len(self.gaia_table) == 0:
            return np.array([]), np.array([])
        
        # Conversion en coordonnées cartésiennes
        detected_xyz = detected_coords.cartesian.xyz.value.T
        
        # Rayon de recherche en distance cartésienne (approximation pour petits angles)
        # Pour un angle θ en radians, distance cartésienne ≈ 2*sin(θ/2)
        radius_rad = dist_limit.to(u.rad).value
        search_radius = 2 * np.sin(radius_rad / 2.0)
        
        # Recherche avec KD-tree
        matches = self.kdtree.query_ball_point(detected_xyz, search_radius)
        
        # Sélection du meilleur match pour chaque source (le plus proche)
        indices = np.full(len(detected_coords), -1, dtype=int)
        distances = np.full(len(detected_coords), np.inf)
        
        for i, match_list in enumerate(matches):
            if len(match_list) > 0:
                # Calcul des distances angulaires réelles
                match_coords = self.gaia_sky[match_list]
                sep = detected_coords[i].separation(match_coords)
                sep_arcsec = sep.arcsec
                
                # Sélection du meilleur match
                best_idx = np.argmin(sep_arcsec)
                if sep_arcsec[best_idx] < dist_limit.to(u.arcsec).value:
                    indices[i] = match_list[best_idx]
                    distances[i] = sep_arcsec[best_idx]
        
        valid_mask = indices >= 0
        logger.debug(f"Matching : {np.sum(valid_mask)}/{len(detected_coords)} sources matchées")
        
        return indices, distances


# =============================================================================
# 3. REQUÊTE GAIA OPTIMISÉE AVEC CACHE
# =============================================================================

def gaia_box_radius_for_imaging_fov(
    target_coord: SkyCoord,
    wcs: WCS,
    shape: Tuple[int, int],
    margin: float = 1.25,
    min_arcmin: float = 20.0,
    max_arcmin: float = 72.0,
) -> u.Quantity:
    """
    Rayon (demi-côté de la boîte Vizier) centré sur T1 pour couvrir tout le CCD.

    Les comparateurs près des bords peuvent être à > 12′ de T1 ; un rayon fixe trop
    petit exclut ces étoiles du catalogue téléchargé (G_Gaia vide dans le SET-UP).
    """
    ny, nx = int(shape[0]), int(shape[1])
    sep_max_deg = 0.0
    for px, py in ((0, 0), (nx - 1, 0), (0, ny - 1), (nx - 1, ny - 1)):
        try:
            cc = wcs.pixel_to_world(px, py)
            sep_max_deg = max(sep_max_deg, float(target_coord.separation(cc).to(u.deg).value))
        except Exception:
            continue
    r_deg = sep_max_deg * float(margin)
    lo = float((min_arcmin * u.arcmin).to(u.deg).value)
    hi = float((max_arcmin * u.arcmin).to(u.deg).value)
    r_deg = float(np.clip(r_deg, lo, hi))
    return r_deg * u.deg


def optimized_query_gaia(
    center_coord: SkyCoord,
    search_radius: u.Quantity,
    mag_limit: float = 18.0,
    gaia_cache: Optional[GaiaCache] = None,
    bypass_cache: bool = False,
    row_limit: int = 5000,
) -> Table:
    """
    Requête Gaia avec cache local persistant.
    
    Parameters
    ----------
    center_coord : SkyCoord
        Centre de la recherche
    search_radius : Quantity
        Rayon de recherche
    mag_limit : float
        Magnitude limite G
    gaia_cache : GaiaCache, optional
        Cache pour éviter les requêtes répétées
    bypass_cache : bool
        Si True, ignore le cache (utile quand une requête antérieure, ex. astrométrie,
        a enregistré seulement les étoiles G ≤ limite plus stricte : les comparateurs
        plus faibles seraient absents du cache).
    row_limit : int
        Limite Vizier (régions denses : au-delà de quelques milliers, des étoiles du bord
        du champ peuvent être absentes de la réponse tronquée).
    
    Returns
    -------
    Table
        Table Gaia avec colonnes Source, RA_ICRS, DE_ICRS, Gmag
    """
    # Vérification du cache
    if gaia_cache and not bypass_cache:
        cached = gaia_cache.get_cached(center_coord, search_radius)
        if cached is not None:
            # Filtre par magnitude si nécessaire
            if 'Gmag' in cached.colnames:
                mask = cached['Gmag'] <= mag_limit
                return cached[mask]
            return cached
    
    # Requête Vizier si pas en cache
    logger.warning(f"⚠️ REQUÊTE RÉSEAU Gaia via Vizier (box search) pour {center_coord} - rayon {search_radius:.3f}°")
    try:
        ra_center = center_coord.ra.deg
        dec_center = center_coord.dec.deg
        radius_deg = search_radius.to(u.deg).value
        
        v = Vizier(
            columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "phot_variable_flag"],
            row_limit=int(row_limit),
        )
        res = v.query_region(
            center_coord,
            width=2 * radius_deg * u.deg,
            height=2 * radius_deg * u.deg,
            catalog="I/355/gaiadr3"
        )
        
        if res and len(res[0]) > 0:
            tab = res[0]
            # Filtre par magnitude
            if 'Gmag' in tab.colnames:
                valid_mask = ~tab['Gmag'].mask if hasattr(tab['Gmag'], 'mask') else ~np.isnan(tab['Gmag'])
                tab = tab[valid_mask & (tab['Gmag'] <= mag_limit)]
            
            # Sauvegarde dans le cache
            if gaia_cache:
                gaia_cache.save_cache(center_coord, search_radius, tab)
            
            logger.info(f"Gaia : {len(tab)} étoiles récupérées")
            return tab
        
        logger.warning("Gaia : Aucune étoile trouvée")
        return Table()
    
    except Exception as e:
        logger.error(f"Erreur requête Gaia : {e}")
        return Table()


def query_gaia_dr2_vizier_box(
    center_coord: SkyCoord,
    search_radius: u.Quantity,
    mag_limit: float = 21.0,
    row_limit: int = 5000,
) -> Table:
    """
    Gaia DR2 via Vizier (I/345/gaia2) — complément lorsque DR3 n’a pas la source ou G vide.

    Retourne une table minimale RA_ICRS, DE_ICRS, Gmag (+ Source si présent).
    """
    try:
        radius_deg = search_radius.to(u.deg).value
        logger.info(
            "Gaia DR2 (Vizier I/345/gaia2), rayon %.4f°, G ≤ %.1f",
            radius_deg,
            mag_limit,
        )
        v = Vizier(
            columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag"],
            row_limit=int(row_limit),
        )
        res = v.query_region(
            center_coord,
            width=2 * radius_deg * u.deg,
            height=2 * radius_deg * u.deg,
            catalog="I/345/gaia2",
        )
        if not res or len(res[0]) == 0:
            return Table(names=("RA_ICRS", "DE_ICRS", "Gmag"), dtype=("f8", "f8", "f8"))
        tab = res[0]
        if "Gmag" not in tab.colnames or "RA_ICRS" not in tab.colnames or "DE_ICRS" not in tab.colnames:
            logger.warning("Gaia DR2 Vizier : colonnes inattendues %s", list(tab.colnames)[:12])
            return Table(names=("RA_ICRS", "DE_ICRS", "Gmag"), dtype=("f8", "f8", "f8"))
        gm = np.asarray(tab["Gmag"], dtype=float)
        ok = np.isfinite(gm) & (gm <= float(mag_limit))
        tab = tab[ok]
        if len(tab) == 0:
            return Table(names=("RA_ICRS", "DE_ICRS", "Gmag"), dtype=("f8", "f8", "f8"))
        out = Table()
        out["RA_ICRS"] = np.asarray(tab["RA_ICRS"], dtype=float)
        out["DE_ICRS"] = np.asarray(tab["DE_ICRS"], dtype=float)
        out["Gmag"] = np.asarray(tab["Gmag"], dtype=float)
        if "Source" in tab.colnames:
            out["Source"] = tab["Source"]
        logger.info("Gaia DR2 : %d étoiles retenues (G ≤ %.1f)", len(out), mag_limit)
        return out
    except Exception as e:
        logger.error("Erreur requête Gaia DR2 (Vizier) : %s", e)
        return Table(names=("RA_ICRS", "DE_ICRS", "Gmag"), dtype=("f8", "f8", "f8"))


def merge_gaia_dr3_dr2_photometry_field(
    tab_dr3: Table,
    tab_dr2: Table,
    dedup_arcsec: float = 1.0,
) -> Table:
    """
    Fusion pour magnitudes comparateurs : DR3 (prioritaire) + sources DR2 sans voisin DR3.

    Sortie : RA_ICRS, DE_ICRS, Gmag uniquement.
    """
    cols = ("RA_ICRS", "DE_ICRS", "Gmag")
    empty = Table(names=cols, dtype=("f8", "f8", "f8"))

    def _slim_dr3(t: Table) -> Table:
        if t is None or len(t) == 0:
            return empty.copy()
        if not all(c in t.colnames for c in cols):
            return empty.copy()
        gm = np.asarray(t["Gmag"], dtype=float)
        ok = np.isfinite(gm)
        if not np.any(ok):
            return empty.copy()
        return Table(
            {
                "RA_ICRS": np.asarray(t["RA_ICRS"], dtype=float)[ok],
                "DE_ICRS": np.asarray(t["DE_ICRS"], dtype=float)[ok],
                "Gmag": gm[ok],
            }
        )

    rows3 = _slim_dr3(tab_dr3)
    if tab_dr2 is None or len(tab_dr2) == 0:
        return rows3
    if not all(c in tab_dr2.colnames for c in cols):
        return rows3
    if len(rows3) == 0:
        return Table(
            {
                "RA_ICRS": np.asarray(tab_dr2["RA_ICRS"], dtype=float),
                "DE_ICRS": np.asarray(tab_dr2["DE_ICRS"], dtype=float),
                "Gmag": np.asarray(tab_dr2["Gmag"], dtype=float),
            }
        )

    c3 = SkyCoord(ra=rows3["RA_ICRS"] * u.deg, dec=rows3["DE_ICRS"] * u.deg)
    extra_ra, extra_de, extra_g = [], [], []
    lim = float(dedup_arcsec)
    for i in range(len(tab_dr2)):
        c2 = SkyCoord(
            ra=float(tab_dr2["RA_ICRS"][i]) * u.deg,
            dec=float(tab_dr2["DE_ICRS"][i]) * u.deg,
        )
        sep_min = float(np.min(c2.separation(c3).arcsec))
        if sep_min > lim:
            extra_ra.append(float(tab_dr2["RA_ICRS"][i]))
            extra_de.append(float(tab_dr2["DE_ICRS"][i]))
            extra_g.append(float(tab_dr2["Gmag"][i]))
    if not extra_ra:
        return rows3
    t_add = Table(
        {
            "RA_ICRS": np.array(extra_ra, dtype=float),
            "DE_ICRS": np.array(extra_de, dtype=float),
            "Gmag": np.array(extra_g, dtype=float),
        }
    )
    return vstack([rows3, t_add])


# =============================================================================
# 4. SÉLECTION INTELLIGENTE D'APERTURES
# =============================================================================

def smart_aperture_selection(fwhm: float, min_ap: float = 2.0, max_ap: float = 8.0) -> np.ndarray:
    """
    Sélection intelligente de 6 apertures pour extrapolation zero-aperture.
    
    Garantit toujours exactement 6 apertures distinctes :
    - 2 apertures petites
    - 2 apertures autour du FWHM
    - 2 apertures plus grandes
    
    Parameters
    ----------
    fwhm : float
        FWHM estimé en pixels
    min_ap : float
        Aperture minimale
    max_ap : float
        Aperture maximale
    
    Returns
    -------
    np.ndarray
        Array de 6 apertures sélectionnées (toujours 6 éléments)
    """
    # Bornes robustes
    min_ap = float(min_ap)
    max_ap = float(max_ap)
    if max_ap <= min_ap:
        max_ap = min_ap + 1.0
    fwhm = float(np.clip(fwhm, min_ap, max_ap))

    # Candidats centrés sur le FWHM + bornes
    apertures_list = [
        min_ap,
        min(max_ap, max(min_ap, fwhm * 0.70)),
        min(max_ap, max(min_ap, fwhm * 0.90)),
        min(max_ap, max(min_ap, fwhm * 1.10)),
        min(max_ap, max(min_ap, fwhm * 1.35)),
        max_ap,
    ]
    
    # Supprimer les doublons en préservant l'ordre, puis compléter si nécessaire
    unique_apertures = []
    seen = set()
    for ap in apertures_list:
        ap_rounded = round(ap, 2)  # Arrondir pour éviter les problèmes de précision
        if ap_rounded not in seen:
            unique_apertures.append(ap)
            seen.add(ap_rounded)
    
    # Si on a moins de 6 apertures uniques, compléter uniformément
    if len(unique_apertures) < 6:
        uniform_apertures = np.linspace(min_ap, max_ap, 6).tolist()
        all_candidates = sorted(set([round(v, 4) for v in (unique_apertures + uniform_apertures)]))
        apertures = np.array(all_candidates[:6], dtype=float)
    else:
        apertures = np.array(sorted(unique_apertures)[:6], dtype=float)
    
    # Garantir exactement 6 apertures (fallback robuste)
    if len(apertures) != 6:
        apertures = np.linspace(min_ap, max_ap, 6)
    apertures = np.sort(apertures)
    
    logger.debug(f"Apertures sélectionnées (FWHM={fwhm:.1f}): {apertures}")
    return apertures


# =============================================================================
# 4.5 ZERO-APERTURE (MÉTHODE ZERO-APERTURE-ASTROMETRY)
# =============================================================================
# Aligné sur https://github.com/bensharkey/Zero-Aperture-Astrometry :
# fit linéaire pondéré RA/Dec vs photAp, extrapolation à aperture=0.

def zero_aperture_fit_radec(
    phot_ap: np.ndarray,
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    rms_ra_arcsec: Optional[np.ndarray] = None,
    rms_dec_arcsec: Optional[np.ndarray] = None,
) -> Tuple[float, float, float, float]:
    """
    Extrapole les positions astrométriques à aperture zéro (méthode Zero-Aperture-Astrometry).

    Fit linéaire pondéré (poids = 1/rms) de RA et Dec en fonction de l'aperture
    photométrique, puis évaluation à aperture=0. Si les RMS ne sont pas fournis,
    le fit est non pondéré.

    Référence: https://github.com/bensharkey/Zero-Aperture-Astrometry

    Parameters
    ----------
    phot_ap : np.ndarray
        Apertures photométriques (même unité que les données d'entrée, ex. pixels ou arcsec).
    ra_deg : np.ndarray
        Ascension droite en degrés pour chaque aperture.
    dec_deg : np.ndarray
        Déclinaison en degrés pour chaque aperture.
    rms_ra_arcsec : np.ndarray, optional
        Incertitude RA en arcsec (pour pondération). Doit être > 0.
    rms_dec_arcsec : np.ndarray, optional
        Incertitude Dec en arcsec (pour pondération). Doit être > 0.

    Returns
    -------
    ra0 : float
        RA extrapolée à aperture=0 (deg).
    dec0 : float
        Dec extrapolée à aperture=0 (deg).
    rms_ra0 : float
        Incertitude RA à aperture=0 (arcsec), estimée (2× ordre de grandeur typique si non fourni).
    rms_dec0 : float
        Incertitude Dec à aperture=0 (arcsec), idem.
    """
    x = np.asarray(phot_ap, dtype=float)
    ra = np.asarray(ra_deg, dtype=float)
    dec = np.asarray(dec_deg, dtype=float)
    if len(x) < 2 or len(x) != len(ra) or len(ra) != len(dec):
        raise ValueError("phot_ap, ra_deg, dec_deg doivent avoir la même longueur (>= 2)")

    # Poids pour fit pondéré (comme Zero-Aperture-Astrometry: w = 1 / rms)
    use_weights_ra = rms_ra_arcsec is not None and len(rms_ra_arcsec) == len(x)
    use_weights_dec = rms_dec_arcsec is not None and len(rms_dec_arcsec) == len(x)
    if use_weights_ra:
        rms_ra = np.asarray(rms_ra_arcsec, dtype=float)
        w_ra = 1.0 / np.maximum(rms_ra, 1e-10)
    else:
        w_ra = None
    if use_weights_dec:
        rms_dec = np.asarray(rms_dec_arcsec, dtype=float)
        w_dec = 1.0 / np.maximum(rms_dec, 1e-10)
    else:
        w_dec = None

    try:
        if w_ra is not None:
            ra_fit, _ = np.polyfit(x, ra, 1, w=w_ra, cov="unscaled")
        else:
            ra_fit, _ = np.polyfit(x, ra, 1, cov="unscaled")
        if w_dec is not None:
            dec_fit, _ = np.polyfit(x, dec, 1, w=w_dec, cov="unscaled")
        else:
            dec_fit, _ = np.polyfit(x, dec, 1, cov="unscaled")
    except Exception as e:
        logger.warning(f"zero_aperture_fit_radec: polyfit échoué, fallback non pondéré: {e}")
        ra_fit, _ = np.polyfit(x, ra, 1, cov="unscaled")
        dec_fit, _ = np.polyfit(x, dec, 1, cov="unscaled")

    ra0 = float(np.polyval(ra_fit, 0.0))
    dec0 = float(np.polyval(dec_fit, 0.0))

    # Incertitudes à aperture 0 : si on a des rms d'entrée, prendre ~2× le min ou la médiane
    if use_weights_ra and np.all(np.isfinite(rms_ra)):
        rms_ra0 = float(2.0 * np.min(rms_ra))  # convention du repo (ordre de grandeur)
    else:
        rms_ra0 = 0.1
    if use_weights_dec and np.all(np.isfinite(rms_dec)):
        rms_dec0 = float(2.0 * np.min(rms_dec))
    else:
        rms_dec0 = 0.1

    return ra0, dec0, rms_ra0, rms_dec0


# =============================================================================
# 5. DÉTECTION GPU OPTIONNELLE
# =============================================================================

def gpu_star_finder(data: np.ndarray, fwhm: float, threshold_sigma: float, use_gpu: bool = False):
    """
    Détection d'étoiles avec support GPU optionnel via CuPy.
    
    Parameters
    ----------
    data : np.ndarray
        Image 2D
    fwhm : float
        FWHM en pixels
    threshold_sigma : float
        Seuil en sigma
    use_gpu : bool
        Utiliser CuPy si disponible
    
    Returns
    -------
    Table ou None
        Table des sources détectées
    """
    if use_gpu and HAS_GPU:
        try:
            logger.debug("Utilisation GPU pour détection d'étoiles")
            
            # Transfert vers GPU
            data_gpu = cp.asarray(data.astype(np.float32))
            
            # Calcul statistiques sur GPU (sigma-clipped)
            # Approximation : on utilise percentile au lieu de sigma-clipping complet
            p25 = cp.percentile(data_gpu, 25)
            p75 = cp.percentile(data_gpu, 75)
            median_gpu = cp.median(data_gpu)
            mad = cp.median(cp.abs(data_gpu - median_gpu))
            std_gpu = 1.4826 * mad  # Conversion MAD vers std pour distribution normale
            
            threshold = float(median_gpu + threshold_sigma * std_gpu)
            
            # Création d'un kernel gaussien pour la convolution (sur GPU)
            kernel_size = int(2 * fwhm) + 1
            if kernel_size % 2 == 0:
                kernel_size += 1
            
            # Kernel gaussien 2D sur GPU
            y, x = cp.ogrid[:kernel_size, :kernel_size]
            center = kernel_size // 2
            sigma = fwhm / 2.355  # Conversion FWHM -> sigma
            kernel = cp.exp(-((x - center)**2 + (y - center)**2) / (2 * sigma**2))
            kernel = kernel / cp.sum(kernel)
            
            # Convolution sur GPU avec cupyx.scipy.ndimage
            from cupyx.scipy import ndimage as ndi_gpu
            convolved = ndi_gpu.convolve(data_gpu, kernel, mode='constant')
            
            # Détection de pics locaux sur GPU
            # On cherche les maxima locaux dans une fenêtre de taille ~FWHM
            window_size = int(fwhm) + 1
            if window_size % 2 == 0:
                window_size += 1
            
            # Maxima locaux : chaque pixel doit être le maximum dans sa fenêtre
            local_max = ndi_gpu.maximum_filter(convolved, size=window_size)
            peaks = (convolved == local_max) & (convolved > threshold)
            
            # Extraction des positions et flux des pics
            peak_coords = cp.argwhere(peaks)
            peak_values = convolved[peaks]
            
            if len(peak_coords) == 0:
                logger.warning("Aucune source détectée sur GPU")
                # Fallback CPU
                _, _, std = sigma_clipped_stats(data)
                daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
                return daofind(data)
            
            # Conversion en Table Astropy (retour CPU)
            peak_coords_cpu = cp.asnumpy(peak_coords)
            peak_values_cpu = cp.asnumpy(peak_values)
            
            # Création d'une table avec les résultats
            sources_table = Table()
            sources_table['xcentroid'] = peak_coords_cpu[:, 1].astype(float)
            sources_table['ycentroid'] = peak_coords_cpu[:, 0].astype(float)
            sources_table['flux'] = peak_values_cpu
            
            # Tri par flux décroissant
            sources_table.sort('flux', reverse=True)
            
            logger.debug(f"GPU : {len(sources_table)} sources détectées")
            return sources_table
            
        except Exception as e:
            logger.warning(f"Erreur GPU, fallback CPU : {e}")
            # Fallback vers CPU
            use_gpu = False
    
    # Détection CPU standard (fallback ou si GPU désactivé)
    _, _, std = sigma_clipped_stats(data)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
    sources = daofind(data)
    
    return sources


# =============================================================================
# 5.5 CALCUL DES STATISTIQUES ASTRONOMIQUES
# =============================================================================

def calculate_astrometric_statistics(
    wcs: WCS,
    pixel_coords: np.ndarray,
    catalog_coords: SkyCoord,
    outlier_threshold: float = 3.0
) -> dict:
    """
    Calcule les statistiques astrométriques détaillées.
    
    Parameters
    ----------
    wcs : WCS
        WCS ajusté
    pixel_coords : np.ndarray
        Coordonnées pixels (N, 2) ou (2, N)
    catalog_coords : SkyCoord
        Coordonnées catalogues (Gaia)
    outlier_threshold : float
        Seuil pour identification des outliers (en sigma)
    
    Returns
    -------
    dict
        Dictionnaire contenant toutes les statistiques :
        - rms_ra, rms_dec, rms_total : RMS séparés et total
        - median_residual : Médiane des résidus totaux
        - n_stars : Nombre d'étoiles
        - n_outliers : Nombre d'outliers
        - residuals_ra, residuals_dec : Résidus individuels (arcsec)
        - outlier_mask : Masque booléen pour les outliers
    """
    # Conversion pixel -> sky
    if pixel_coords.shape[0] == 2 and pixel_coords.shape[1] > 2:
        # Format (2, N) -> transposer
        pixel_coords = pixel_coords.T
    
    calc_coords = wcs.pixel_to_world(pixel_coords[:, 0], pixel_coords[:, 1])
    
    # Calcul des résidus en arcsec
    # RA : multiplier par cos(dec) pour la distance angulaire
    ra_calc = calc_coords.ra.deg
    dec_calc = calc_coords.dec.deg
    ra_cat = catalog_coords.ra.deg
    dec_cat = catalog_coords.dec.deg
    
    # Résidus RA (corrigés par cos(dec))
    residuals_ra = (ra_calc - ra_cat) * np.cos(np.radians(dec_cat)) * 3600.0  # arcsec
    residuals_dec = (dec_calc - dec_cat) * 3600.0  # arcsec
    
    # Résidus totaux (séparation angulaire)
    sep = catalog_coords.separation(calc_coords)
    residuals_total = sep.arcsec
    
    # Statistiques RMS
    rms_ra = np.sqrt(np.mean(residuals_ra**2))
    rms_dec = np.sqrt(np.mean(residuals_dec**2))
    rms_total = np.sqrt(np.mean(residuals_total**2))
    
    # Médiane
    median_residual = np.median(residuals_total)
    
    # Identification des outliers (> threshold * sigma)
    std_total = np.std(residuals_total)
    outlier_mask = residuals_total > (outlier_threshold * std_total)
    n_outliers = np.sum(outlier_mask)
    
    # Calcul du coefficient de corrélation RMS (rmsCorr pour format ADES)
    # rmsCorr = corrélation entre les résidus RA et Dec
    # Formule : r = cov(RA, Dec) / (σ_RA * σ_Dec)
    if len(residuals_ra) > 1 and len(residuals_dec) > 1:
        std_ra = np.std(residuals_ra)
        std_dec = np.std(residuals_dec)
        if std_ra > 0 and std_dec > 0:
            # Covariance normalisée (coefficient de corrélation de Pearson)
            correlation = np.corrcoef(residuals_ra, residuals_dec)[0, 1]
            # Si NaN (peut arriver si variance nulle), mettre à 0
            if not np.isfinite(correlation):
                correlation = 0.0
        else:
            correlation = 0.0
    else:
        correlation = 0.0
    
    return {
        'rms_ra': float(rms_ra),
        'rms_dec': float(rms_dec),
        'rms_total': float(rms_total),
        'median_residual': float(median_residual),
        'n_stars': len(residuals_total),
        'n_outliers': int(n_outliers),
        'rms_corr': float(correlation),  # Coefficient de corrélation RA/Dec
        'residuals_ra': residuals_ra,
        'residuals_dec': residuals_dec,
        'residuals_total': residuals_total,
        'outlier_mask': outlier_mask
    }


# =============================================================================
# 6. PIPELINE D'ASTROMÉTRIE OPTIMISÉ COMPLET
# =============================================================================

def fast_astrometry_solve(
    data: np.ndarray,
    header: fits.Header,
    gaia_cache: Optional[GaiaCache] = None,
    fwhm: float = 5.0,
    threshold_sigma: float = 3.0,
    max_sources: int = 500,
    match_radius: u.Quantity = 5.0*u.arcsec,
    mag_limit: float = 18.0,
    use_gpu: bool = False,
    skip_zero_aperture: bool = False,
    extrapolation_params: Optional[Tuple[float, float]] = None
) -> Tuple[WCS, float, Dict]:
    """
    Pipeline d'astrométrie optimisé avec toutes les améliorations.
    
    Parameters
    ----------
    data : np.ndarray
        Image 2D
    header : fits.Header
        Header FITS avec WCS initial (optionnel)
    gaia_cache : GaiaCache, optional
        Cache Gaia pour éviter requêtes répétées
    fwhm : float
        FWHM estimé en pixels
    threshold_sigma : float
        Seuil de détection en sigma
    max_sources : int
        Nombre maximum de sources à détecter
    match_radius : Quantity
        Rayon de matching en arcsec
    mag_limit : float
        Magnitude limite Gaia
    use_gpu : bool
        Utiliser GPU si disponible
    
    Returns
    -------
    wcs : WCS
        WCS ajusté
    zero_rms : float
        RMS extrapolé à zero-aperture en arcsec
    """
    logger.info(f"=== Démarrage astrométrie optimisée (skip_zero_aperture={skip_zero_aperture}) ===")
    
    # 1. WCS initial depuis header (correction SIP/CTYPE pour éviter l'avertissement Astropy)
    try:
        header_fixed = fix_header_sip_ctype(header)
        wcs_init = WCS(header_fixed)
        if not wcs_init.is_celestial:
            raise ValueError("WCS initial invalide")
    except Exception:
        logger.warning("Pas de WCS initial, création d'un WCS par défaut")
        # Création d'un WCS minimal (nécessiterait estimation FOV)
        raise ValueError("WCS initial requis pour l'instant")
    
    # 2. Calcul centre et rayon FOV
    ny, nx = data.shape
    px = np.array([0, nx, 0, nx])
    py = np.array([0, 0, ny, ny])
    corners_sky = wcs_init.pixel_to_world(px, py)
    center_ra = np.mean(corners_sky.ra.deg)
    center_dec = np.mean(corners_sky.dec.deg)
    center_coord = SkyCoord(center_ra * u.deg, center_dec * u.deg, frame="icrs")
    radius = center_coord.separation(corners_sky).max() * 1.1
    
    # 3. Requête Gaia avec cache
    gaia_table = optimized_query_gaia(
        center_coord, radius, mag_limit=mag_limit, gaia_cache=gaia_cache
    )
    
    if len(gaia_table) == 0:
        raise ValueError("Aucune étoile Gaia trouvée")
    
    # 4. Détection des sources
    sources = gpu_star_finder(data, fwhm, threshold_sigma, use_gpu=use_gpu)
    if sources is None or len(sources) == 0:
        raise ValueError("Aucune source détectée")
    
    sources.sort('flux', reverse=True)
    sources = sources[:max_sources]
    
    # Affinage des centroïdes
    refined_x = []
    refined_y = []
    for src in sources:
        try:
            y_int = int(src['ycentroid'])
            x_int = int(src['xcentroid'])
            y_slice = slice(max(0, y_int-5), min(data.shape[0], y_int+6))
            x_slice = slice(max(0, x_int-5), min(data.shape[1], x_int+6))
            cutout = data[y_slice, x_slice]
            
            # Vérifier que le cutout a suffisamment de pixels (> 6)
            if cutout.size < 6:
                refined_x.append(src['xcentroid'])
                refined_y.append(src['ycentroid'])
                continue
                
            xc, yc = centroid_2dg(cutout)
            refined_x.append(src['xcentroid'] - 5 + xc)
            refined_y.append(src['ycentroid'] - 5 + yc)
        except Exception as e:
            logger.debug(f"Échec affinage centroïde source {len(refined_x)}: {e}")
            refined_x.append(src['xcentroid'])
            refined_y.append(src['ycentroid'])
    
    sources['xcentroid'] = refined_x
    sources['ycentroid'] = refined_y
    
    # 5. Matching rapide avec KD-tree
    detected_sky = wcs_init.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
    matcher = FastStarMatcher(gaia_table, gaia_cache)
    indices, distances = matcher.match(detected_sky, match_radius)
    
    valid_mask = indices >= 0
    n_matches = np.sum(valid_mask)
    
    # Diagnostic détaillé si pas assez de matches
    if n_matches < 3:
        logger.warning(f"⚠️ Seulement {n_matches} matches trouvés (requis: 3 minimum)")
        logger.warning(f"   - Sources détectées : {len(sources)}")
        logger.warning(f"   - Étoiles Gaia dans région : {len(gaia_table)}")
        logger.warning(f"   - Rayon de matching : {match_radius:.1f}")
        logger.warning(f"   - Magnitude limite Gaia : {mag_limit:.1f}")
        if n_matches > 0:
            mean_dist = np.mean(distances[valid_mask])
            logger.warning(f"   - Distance moyenne des matches : {mean_dist:.2f} arcsec")
        
        # Tentative avec un rayon de matching plus large
        if n_matches < 3 and match_radius.value < 30.0:
            larger_radius = match_radius * 2.0
            logger.info(f"   → Tentative avec rayon augmenté à {larger_radius:.1f}")
            indices, distances = matcher.match(detected_sky, larger_radius)
            valid_mask = indices >= 0
            n_matches = np.sum(valid_mask)
            if n_matches >= 3:
                logger.info(f"   ✓ Succès avec rayon élargi : {n_matches} matches")
            else:
                logger.warning(f"   ✗ Échec même avec rayon élargi : {n_matches} matches")
        
        if n_matches < 3:
            raise ValueError(f"Pas assez de matches ({n_matches} < 3). Suggestions : augmenter le rayon de matching, réduire le seuil de détection, ou vérifier la qualité du WCS initial.")
    
    # 6. Construction table de matches
    matches = Table()
    matches['xcentroid'] = sources['xcentroid'][valid_mask]
    matches['ycentroid'] = sources['ycentroid'][valid_mask]
    matches['ra_cat'] = gaia_table['RA_ICRS'][indices[valid_mask]]
    matches['dec_cat'] = gaia_table['DE_ICRS'][indices[valid_mask]]
    
    # 7. Fit WCS
    xy_pixels = np.array([matches['xcentroid'], matches['ycentroid']]).T
    world_coords = SkyCoord(
        ra=matches['ra_cat'], 
        dec=matches['dec_cat'], 
        unit=u.deg
    )
    wcs_fitted = fit_wcs_from_points(
        xy_pixels.T, world_coords, proj_point='center', sip_degree=None
    )
    
    # 8. Extrapolation zero-aperture (sélection intelligente de 6 apertures)
    # Si skip_zero_aperture est True, on saute l'extrapolation et on utilise le RMS de base
    fit_r2 = None  # Initialiser fit_r2
    offset_ra0_arcsec = None
    offset_dec0_arcsec = None
    wcs_zero_aperture = None
    if skip_zero_aperture:
        logger.debug("Extrapolation zero-aperture ignorée (optimisation)")
        calc_sky = wcs_fitted.pixel_to_world(matches['xcentroid'], matches['ycentroid'])
        sep = world_coords.separation(calc_sky)
        rms_base = np.sqrt(np.mean(sep.arcsec**2))
        zero_rms = rms_base
    else:
        apertures = smart_aperture_selection(fwhm)
        rms_list = []
        # Résidus moyens RA/Dec par aperture (approche Zero-Aperture-Astrometry : fit puis extrapolation à 0)
        mean_residual_ra_list = []
        mean_residual_dec_list = []
        # Série exploitable côté UI (une ligne par ouverture)
        aperture_metrics = []
        
        logger.info(f"Extrapolation zero-aperture avec {len(apertures)} apertures : {apertures}")
        
        for ap_radius in apertures:
            try:
                # Affinage des centroïdes avec box_size = aperture
                box_size = int(2 * ap_radius) + 1
                if box_size % 2 == 0:
                    box_size += 1
                
                refined_x_ap = []
                refined_y_ap = []
                
                for src in sources[valid_mask]:
                    x_init = src['xcentroid']
                    y_init = src['ycentroid']
                    
                    # Affinage avec box_size adapté à l'aperture
                    try:
                        h, w = data.shape
                        x_int, y_int = int(round(x_init)), int(round(y_init))
                        r = box_size // 2
                        
                        # Vérifier les limites avec marges de sécurité
                        if r <= x_int < w - r and r <= y_int < h - r:
                            cutout = data[y_int - r:y_int + r + 1, x_int - r:x_int + r + 1]
                            
                            # Vérifier que le cutout a suffisamment de pixels non masqués (> 6)
                            if cutout.size < 6:
                                refined_x_ap.append(x_init)
                                refined_y_ap.append(y_init)
                                continue
                            
                            # Vérifier les valeurs NaN/inf
                            if np.any(~np.isfinite(cutout)):
                                valid_mask_cutout = np.isfinite(cutout)
                                if np.sum(valid_mask_cutout) < 6:
                                    refined_x_ap.append(x_init)
                                    refined_y_ap.append(y_init)
                                    continue
                            
                            xc, yc = centroid_2dg(cutout)
                            refined_x_ap.append(x_int - r + xc)
                            refined_y_ap.append(y_int - r + yc)
                        else:
                            refined_x_ap.append(x_init)
                            refined_y_ap.append(y_init)
                    except Exception as e:
                        logger.debug(f"Échec affinage centroïde aperture {ap_radius:.1f}px: {e}")
                        refined_x_ap.append(x_init)
                        refined_y_ap.append(y_init)
                
                # Fit WCS avec centroïdes affinés
                xy_pixels_ap = np.array([refined_x_ap, refined_y_ap]).T
                wcs_ap = fit_wcs_from_points(
                    xy_pixels_ap.T, world_coords, proj_point='center', sip_degree=None
                )
                
                # Calcul RMS pour cette aperture
                calc_sky_ap = wcs_ap.pixel_to_world(refined_x_ap, refined_y_ap)
                sep_ap = world_coords.separation(calc_sky_ap)
                rms_ap = np.sqrt(np.mean(sep_ap.arcsec**2))
                rms_list.append(rms_ap)
                
                # Résidus RA/Dec moyens (arcsec) pour cette aperture — approche type Zero-Aperture-Astrometry
                ra_calc_ap = calc_sky_ap.ra.deg
                dec_calc_ap = calc_sky_ap.dec.deg
                ra_cat = world_coords.ra.deg
                dec_cat = world_coords.dec.deg
                res_ra_ap = (ra_calc_ap - ra_cat) * np.cos(np.radians(dec_cat)) * 3600.0  # arcsec
                res_dec_ap = (dec_calc_ap - dec_cat) * 3600.0  # arcsec
                mean_residual_ra_list.append(float(np.mean(res_ra_ap)))
                mean_residual_dec_list.append(float(np.mean(res_dec_ap)))
                rms_ra_ap = float(np.sqrt(np.mean(res_ra_ap**2)))
                rms_dec_ap = float(np.sqrt(np.mean(res_dec_ap**2)))
                aperture_metrics.append({
                    'phot_ap': float(ap_radius),
                    'rms_ra': rms_ra_ap,
                    'rms_dec': rms_dec_ap,
                    'rms_fit': float(rms_ap),
                    'mean_residual_ra_arcsec': float(mean_residual_ra_list[-1]),
                    'mean_residual_dec_arcsec': float(mean_residual_dec_list[-1]),
                    'selected': True,
                })
                
                logger.info(f"  Aperture {ap_radius:.1f}px : RMS = {rms_ap:.4f}\"")
            
            except Exception as e:
                logger.warning(f"Erreur aperture {ap_radius:.1f}px : {e}")
                # Utilise le RMS de base si échec
                if len(rms_list) == 0:
                    calc_sky = wcs_fitted.pixel_to_world(matches['xcentroid'], matches['ycentroid'])
                    sep = world_coords.separation(calc_sky)
                    rms_base = np.sqrt(np.mean(sep.arcsec**2))
                    rms_list.append(rms_base)
                    mean_residual_ra_list.append(0.0)
                    mean_residual_dec_list.append(0.0)
                else:
                    rms_list.append(rms_list[-1])  # Répète la dernière valeur
                    mean_residual_ra_list.append(mean_residual_ra_list[-1])
                    mean_residual_dec_list.append(mean_residual_dec_list[-1])
                # Conserver une ligne même en fallback
                if len(aperture_metrics) > 0:
                    prev = aperture_metrics[-1]
                    aperture_metrics.append({
                        'phot_ap': float(ap_radius),
                        'rms_ra': float(prev['rms_ra']),
                        'rms_dec': float(prev['rms_dec']),
                        'rms_fit': float(prev['rms_fit']),
                        'mean_residual_ra_arcsec': float(prev.get('mean_residual_ra_arcsec', 0.0)),
                        'mean_residual_dec_arcsec': float(prev.get('mean_residual_dec_arcsec', 0.0)),
                        'selected': True,
                    })
                else:
                    aperture_metrics.append({
                        'phot_ap': float(ap_radius),
                        'rms_ra': 0.0,
                        'rms_dec': 0.0,
                        'rms_fit': float(rms_list[-1]),
                        'mean_residual_ra_arcsec': 0.0,
                        'mean_residual_dec_arcsec': 0.0,
                        'selected': True,
                    })
        
        # Filtrage robuste des ouvertures aberrantes avant extrapolation.
        # Exemple typique: une ou deux apertures avec RMS catastrophique (fit WCS instable).
        rms_arr_all = np.asarray(rms_list, dtype=float)
        valid_rms_mask = np.isfinite(rms_arr_all) & (rms_arr_all > 0)
        if np.any(valid_rms_mask):
            med_rms = float(np.median(rms_arr_all[valid_rms_mask]))
            mad_rms = float(np.median(np.abs(rms_arr_all[valid_rms_mask] - med_rms)))
            # Seuils souples: borne relative + borne robuste (MAD), et borne absolue anti-divergence.
            rel_cap = max(5.0 * med_rms, med_rms + 10.0)
            robust_cap = med_rms + 6.0 * max(mad_rms, 1e-6)
            hard_cap = min(1e4, max(rel_cap, robust_cap))
            valid_rms_mask &= (rms_arr_all <= hard_cap)
        else:
            valid_rms_mask = np.zeros_like(rms_arr_all, dtype=bool)

        # Mettre à jour le marquage selected dans la série UI.
        for i_m, m in enumerate(aperture_metrics):
            if i_m < len(valid_rms_mask):
                m["selected"] = bool(valid_rms_mask[i_m])

        if np.sum(valid_rms_mask) < 2:
            logger.warning(
                "Zero-aperture: pas assez d'ouvertures valides après filtrage robuste; "
                "fallback sur toutes les ouvertures valides numériques."
            )
            valid_rms_mask = np.isfinite(rms_arr_all) & (rms_arr_all > 0)

        apertures_fit = np.asarray(apertures, dtype=float)[valid_rms_mask]
        rms_fit = rms_arr_all[valid_rms_mask]
        ra_fit = np.asarray(mean_residual_ra_list, dtype=float)[valid_rms_mask]
        dec_fit = np.asarray(mean_residual_dec_list, dtype=float)[valid_rms_mask]

        rejected_idx = np.where(~valid_rms_mask)[0]
        for ridx in rejected_idx:
            logger.warning(
                f"  Aperture {float(apertures[ridx]):.2f}px rejetée (RMS={float(rms_arr_all[ridx]):.4f}\")"
            )

        # Extrapolation vers zero-aperture (linéaire ou quadratique)
        if len(rms_fit) >= 2 and len(apertures_fit) == len(rms_fit):
            # Essayer d'abord extrapolation quadratique si on a au moins 3 points
            use_quadratic = len(rms_fit) >= 3
            zero_rms = None
            fit_r2 = None
            
            quad_r2_min = 0.9
            lin_r2_min = 0.5
            if use_quadratic:
                try:
                    # Fit quadratique pondéré (w = 1/rms, méthode Zero-Aperture-Astrometry)
                    rms_arr = np.array(rms_fit, dtype=float)
                    w_rms = 1.0 / np.maximum(rms_arr, 1e-10)
                    coeffs_quad = np.polyfit(apertures_fit, rms_fit, 2, w=w_rms)
                    zero_rms_quad = coeffs_quad[2]  # Ordonnée à l'origine (aperture = 0)
                    
                    # Calcul R² pour validation
                    rms_pred_quad = np.polyval(coeffs_quad, apertures_fit)
                    ss_res = np.sum((rms_fit - rms_pred_quad)**2)
                    ss_tot = np.sum((rms_fit - np.mean(rms_fit))**2)
                    fit_r2_quad = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    
                    # Validation : zero_rms doit être positif et raisonnable, R² > seuil
                    if zero_rms_quad >= 0 and zero_rms_quad <= max(rms_fit) * 2 and fit_r2_quad > quad_r2_min:
                        zero_rms = zero_rms_quad
                        fit_r2 = fit_r2_quad
                        logger.info(f"Extrapolation quadratique : RMS(0) = {zero_rms_quad:.4f}\" (R²={fit_r2_quad:.3f})")
                    else:
                        logger.info(
                            f"Extrapolation quadratique rejetée (R²={fit_r2_quad:.3f} < {quad_r2_min})"
                        )
                except Exception as e:
                    logger.debug(f"Extrapolation quadratique échouée : {e}")
            
            # Si quadratique n'a pas fonctionné ou pas assez de points, utiliser linéaire
            if zero_rms is None:
                # Fit linéaire pondéré (w = 1/rms, méthode Zero-Aperture-Astrometry)
                rms_arr = np.array(rms_fit, dtype=float)
                w_rms = 1.0 / np.maximum(rms_arr, 1e-10)
                coeffs = np.polyfit(apertures_fit, rms_fit, 1, w=w_rms)
                zero_rms = coeffs[1]  # Ordonnée à l'origine (aperture = 0)
                
                # Calcul R²
                rms_pred = np.polyval(coeffs, apertures_fit)
                ss_res = np.sum((rms_fit - rms_pred)**2)
                ss_tot = np.sum((rms_fit - np.mean(rms_fit))**2)
                fit_r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                
                logger.info(f"Extrapolation linéaire : RMS(0) = {zero_rms:.4f}\" (pente={coeffs[0]:.4f}, R²={fit_r2:.3f})")
                if fit_r2 < lin_r2_min:
                    zero_rms = min(rms_fit)
                    logger.warning(
                        f"Extrapolation linéaire instable (R²={fit_r2:.3f} < {lin_r2_min}), "
                        "utilisation du minimum observé"
                    )
            
            # Validation finale : zero_rms doit être positif et raisonnable
            if zero_rms < 0:
                zero_rms = min(rms_fit)  # Prend le minimum si extrapolation négative
                logger.warning("Extrapolation zero-aperture négative, utilisation du minimum observé")
            if zero_rms > max(rms_fit) * 2:
                zero_rms = min(rms_fit)  # Limite si extrapolation trop grande
                logger.warning("Extrapolation zero-aperture excessive, utilisation du minimum observé")
        else:
            # Fallback : utilise le RMS minimum
            calc_sky = wcs_fitted.pixel_to_world(matches['xcentroid'], matches['ycentroid'])
            sep = world_coords.separation(calc_sky)
            rms_base = np.sqrt(np.mean(sep.arcsec**2))
            zero_rms = rms_base
            fit_r2 = None
            logger.warning("Extrapolation zero-aperture échouée, utilisation RMS de base")
        
        # Approche Zero-Aperture-Astrometry : fit résidus moyens RA/Dec vs aperture, extrapolation à 0
        offset_ra0_arcsec = None
        offset_dec0_arcsec = None
        wcs_zero_aperture = None
        if len(ra_fit) >= 2 and len(ra_fit) == len(apertures_fit) and len(dec_fit) == len(apertures_fit):
            try:
                rms_arr = np.array(rms_fit, dtype=float)
                w_rms = 1.0 / np.maximum(rms_arr, 1e-10)
                coeffs_ra = np.polyfit(apertures_fit, ra_fit, 1, w=w_rms)
                coeffs_dec = np.polyfit(apertures_fit, dec_fit, 1, w=w_rms)
                offset_ra0_arcsec = float(np.polyval(coeffs_ra, 0.0))
                offset_dec0_arcsec = float(np.polyval(coeffs_dec, 0.0))
                # WCS corrigé : décalage crval pour annuler l'offset systématique à aperture 0
                wcs_zero_aperture = wcs_fitted.deepcopy()
                crval_ra, crval_dec = wcs_zero_aperture.wcs.crval[0], wcs_zero_aperture.wcs.crval[1]
                wcs_zero_aperture.wcs.crval[0] -= offset_ra0_arcsec / (3600.0 * np.cos(np.radians(crval_dec)))
                wcs_zero_aperture.wcs.crval[1] -= offset_dec0_arcsec / 3600.0
                logger.info(
                    f"Résidus RA/Dec vs aperture (extrap. 0) : offset_ra0={offset_ra0_arcsec:.4f}\", "
                    f"offset_dec0={offset_dec0_arcsec:.4f}\""
                )
            except Exception as e:
                logger.debug(f"Fit résidus RA/Dec vs aperture échoué : {e}")
    
    # RMS classique avec centroïdes affinés (comme dans solve_astrometry_classical)
    # Calcul de l'ouverture classique (r_ap = 1.4 * FWHM)
    r_ap_classical = max(2.0, min(round(1.4 * fwhm, 1), 15.0))
    box_size_classical = int(2 * r_ap_classical) + 1
    if box_size_classical % 2 == 0:
        box_size_classical += 1
    
    refined_x_classical = []
    refined_y_classical = []
    
    for src in sources[valid_mask]:
        x_init = src['xcentroid']
        y_init = src['ycentroid']
        
        try:
            h, w = data.shape
            x_int, y_int = int(round(x_init)), int(round(y_init))
            r = box_size_classical // 2
            
            if r <= x_int < w - r and r <= y_int < h - r:
                cutout = data[y_int - r:y_int + r + 1, x_int - r:x_int + r + 1]
                if cutout.size >= 6 and np.isfinite(cutout).any():
                    xc, yc = centroid_2dg(cutout)
                    refined_x_classical.append(x_int - r + xc)
                    refined_y_classical.append(y_int - r + yc)
                else:
                    refined_x_classical.append(x_init)
                    refined_y_classical.append(y_init)
            else:
                refined_x_classical.append(x_init)
                refined_y_classical.append(y_init)
        except Exception as e:
            logger.debug(f"Échec affinage centroïde classique : {e}")
            refined_x_classical.append(x_init)
            refined_y_classical.append(y_init)
    
    # Fit WCS avec centroïdes affinés classiques
    xy_pixels_classical = np.array([refined_x_classical, refined_y_classical]).T
    wcs_classical = fit_wcs_from_points(
        xy_pixels_classical.T, world_coords, proj_point='center', sip_degree=None
    )
    
    # Calcul des statistiques classiques avec centroïdes affinés
    stats_classical = calculate_astrometric_statistics(
        wcs_classical,
        xy_pixels_classical,
        world_coords,
        outlier_threshold=3.0
    )
    rms = stats_classical['rms_total']
    
    # Choisir le WCS final : le RMS « zero-aperture » est une extrapolation à r=0 sur la courbe
    # RMS(aperture) ; il peut être plus pessimiste que le RMS du fit WCS à ouverture classique.
    # Ancien critère (zero > 1.5 × classique) gardait souvent le zéro alors que le classique était déjà meilleur.
    _prefer_classical = float(zero_rms) > float(rms) + 1e-9
    if _prefer_classical:
        logger.info(
            f"Zero-aperture RMS ({zero_rms:.4f}\") > RMS classique ({rms:.4f}\"); "
            "on conserve le WCS et le RMS classiques (meilleure cohérence métrique)."
        )
        final_wcs = wcs_classical
        final_stats = stats_classical
        final_rms = rms
        method_used = "classical"
    else:
        # Méthode zero-aperture : utiliser le WCS corrigé (résidus RA/Dec extrapolés à 0) si disponible
        if wcs_zero_aperture is not None:
            final_wcs = wcs_zero_aperture
        else:
            final_wcs = wcs_fitted
        final_stats = stats_classical
        final_rms = zero_rms
        method_used = "zero-aperture"
    
    # Ajouter les informations supplémentaires au dictionnaire de statistiques
    final_stats['method'] = method_used
    final_stats['rms_zero_aperture'] = float(zero_rms)
    final_stats['rms_classical'] = rms
    final_stats['fit_r2'] = fit_r2
    # Offsets résidus RA/Dec extrapolés à aperture 0 (approche Zero-Aperture-Astrometry)
    final_stats['residual_offset_ra0_arcsec'] = offset_ra0_arcsec
    final_stats['residual_offset_dec0_arcsec'] = offset_dec0_arcsec
    final_stats['wcs_zero_aperture'] = wcs_zero_aperture
    final_stats['zero_aperture_series'] = aperture_metrics if not skip_zero_aperture else []
    
    logger.info(f"Astrométrie optimisée réussie : RMS classique={rms:.4f}\" (RA={stats_classical['rms_ra']:.4f}\", Dec={stats_classical['rms_dec']:.4f}\"), "
                f"Zero-aperture RMS={zero_rms:.4f}\", RMS final={final_rms:.4f}\" ({method_used}), "
                f"{stats_classical['n_outliers']}/{stats_classical['n_stars']} outliers")
    
    return final_wcs, final_rms, final_stats


# =============================================================================
# 7. CLASSE PRINCIPALE POUR INTÉGRATION
# =============================================================================

class FastAstrometrySolver:
    """Classe principale pour résolution astrométrie optimisée."""
    
    def __init__(
        self,
        gaia_cache: Optional[GaiaCache] = None,
        fwhm: float = 5.0,
        threshold_sigma: float = 3.0,
        max_sources: int = 500,
        match_radius: u.Quantity = 5.0*u.arcsec,
        mag_limit: float = 18.0,
        use_gpu: bool = False
    ):
        self.gaia_cache = gaia_cache or GaiaCache()
        self.fwhm = fwhm
        self.threshold_sigma = threshold_sigma
        self.max_sources = max_sources
        self.match_radius = match_radius
        self.mag_limit = mag_limit
        self.use_gpu = use_gpu and HAS_GPU
    
    def solve(self, data: np.ndarray, header: fits.Header, skip_zero_aperture: bool = False, extrapolation_params: Optional[Tuple[float, float]] = None) -> Tuple[WCS, float, Dict]:
        """Résout l'astrométrie avec toutes les optimisations."""
        return fast_astrometry_solve(
            data=data,
            header=header,
            gaia_cache=self.gaia_cache,
            fwhm=self.fwhm,
            threshold_sigma=self.threshold_sigma,
            max_sources=self.max_sources,
            match_radius=self.match_radius,
            mag_limit=self.mag_limit,
            use_gpu=self.use_gpu,
            skip_zero_aperture=skip_zero_aperture,
            extrapolation_params=extrapolation_params
        )

