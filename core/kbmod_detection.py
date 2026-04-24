# core/kbmod_detection.py
"""
Module optionnel pour la d?tection d'ast?ro?des par Synthetic Tracking (KBMOD).
N?cessite: pip install (depuis https://github.com/dirac-institute/kbmod), GPU NVIDIA, CUDA.
Si KBMOD n'est pas install?, HAS_KBMOD est False et run_kbmod_detection l?ve ImportError.
"""
import logging
from pathlib import Path
from typing import List, Dict, Optional, Callable, Any

import numpy as np

logger = logging.getLogger(__name__)

HAS_KBMOD = False
_kbmod_search = None
_kbmod_psf = None
_kbmod_util = None

try:
    import kbmod.search as _kbmod_search
    from kbmod.core.psf import PSF as _PSF
    _kbmod_psf = _PSF
    try:
        from kbmod.util_functions import load_layered_image
        _kbmod_util = {"load_layered_image": load_layered_image}
    except ImportError:
        try:
            from kbmod.util_functions import load_deccam_layered_image
            _kbmod_util = {"load_deccam_layered_image": load_deccam_layered_image}
        except ImportError:
            _kbmod_util = {}
    HAS_KBMOD = True
except ImportError as e:
    logger.debug(f"KBMOD non disponible: {e}")


def run_kbmod_detection(
    fits_paths: List[str],
    wcs_ref: Any,
    scale_arcsec_per_px: float = 1.0,
    velocity_arcsec_per_min_min: float = 0.1,
    velocity_arcsec_per_min_max: float = 10.0,
    max_results: int = 500,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> List[Dict[str, Any]]:
    """
    Lance une d?tection KBMOD (Synthetic Tracking) sur une liste de FITS.

    Args:
        fits_paths: Chemins vers les fichiers FITS (m?me champ, temps diff?rents).
        wcs_ref: WCS de la premi?re image (astropy.wcs.WCS) pour convertir (x0, y0) en RA/Dec.
        scale_arcsec_per_px: ?chelle en "/px pour convertir vitesses (arcsec/min -> px/day).
        velocity_arcsec_per_min_min: Vitesse min en "/min (recherche).
        velocity_arcsec_per_min_max: Vitesse max en "/min (recherche).
        max_results: Nombre max de candidats ? retourner.
        progress_callback: Optional callback(step, total_steps) pendant le calcul.

    Returns:
        Liste de dicts: ra_deg, dec_deg, x0, y0, vx_px_per_day, vy_px_per_day, likelihood, jd_ref.
    """
    if not HAS_KBMOD or _kbmod_search is None or _kbmod_psf is None:
        raise ImportError(
            "KBMOD n'est pas install? ou pas importable. "
            "Installation: voir docs/INTEGRATION_KBMOD_NPOAP.md. "
            "Pr?requis: GPU NVIDIA, CUDA, puis git clone https://github.com/dirac-institute/kbmod && pip install ."
        )

    from astropy.io import fits
    from astropy.time import Time

    fits_paths = [Path(p) for p in fits_paths]
    if not fits_paths:
        return []

    # PSF simple (sigma ~1 px)
    psf = _kbmod_psf(1.0)
    if hasattr(psf, "kernel"):
        kernel = psf.kernel
    else:
        kernel = np.ones((3, 3)) / 9.0
        if hasattr(psf, "set_kernel"):
            psf.set_kernel(kernel)

    # Construire ImageStack ? partir des FITS
    imstack = _kbmod_search.ImageStack()
    times_used = []

    for i, fpath in enumerate(fits_paths):
        if progress_callback:
            progress_callback(i, len(fits_paths) + 2)
        try:
            with fits.open(str(fpath)) as hdul:
                data = hdul[0].data.astype(np.float32)
                if data is None or not data.size:
                    continue
                header = hdul[0].header
                # JD ou DATE-OBS
                jd = header.get("JD-UTC", 0.0) or header.get("JD", 0.0)
                if jd == 0.0:
                    date_obs = header.get("DATE-OBS", header.get("DATE"))
                    if date_obs:
                        exptime = float(header.get("EXPTIME", header.get("EXPOSURE", 0.0)))
                        t = Time(date_obs, scale="utc")
                        jd = float(t.jd) + exptime / (2.0 * 86400.0)
                if jd == 0.0:
                    continue
                # Variance: approximation (read_noise^2 ou Poisson)
                gain = float(header.get("GAIN", header.get("EGAIN", 1.0)))
                rn = float(header.get("RDNOISE", header.get("READNOIS", 5.0)))
                var = np.abs(data) / np.clip(gain, 1e-6, None) + rn * rn
                var = np.clip(var, 1e-6, None).astype(np.float32)
                mask = np.zeros_like(data, dtype=np.float32)
        except Exception as e:
            logger.warning(f"Skip {fpath}: {e}")
            continue

        # Cr?er LayeredImage: l'API KBMOD peut exiger RawImage ou buffers
        try:
            if hasattr(_kbmod_search, "LayeredImage"):
                # Certaines versions: LayeredImage(science, mask, variance, time)
                if hasattr(_kbmod_search, "RawImage"):
                    sci = _kbmod_search.RawImage(data)
                    msk = _kbmod_search.RawImage(mask)
                    var = _kbmod_search.RawImage(var)
                    limg = _kbmod_search.LayeredImage(sci, msk, var, jd, psf)
                else:
                    limg = _kbmod_search.LayeredImage(data, mask, var, jd, psf)
            else:
                # Fallback: ?crire FITS temporaire 3 extensions et charger
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as tmp:
                    tmp_path = tmp.name
                try:
                    hdu_sci = fits.PrimaryHDU(data)
                    hdu_var = fits.ImageHDU(var, name="VARIANCE")
                    hdu_mask = fits.ImageHDU(mask, name="MASK")
                    hdul = fits.HDUList([hdu_sci, hdu_var, hdu_mask])
                    hdul.writeto(tmp_path, overwrite=True)
                    hdul.close()
                    if _kbmod_util.get("load_layered_image"):
                        limg = _kbmod_util["load_layered_image"](tmp_path, kernel)
                    elif _kbmod_util.get("load_deccam_layered_image"):
                        limg = _kbmod_util["load_deccam_layered_image"](tmp_path, kernel)
                    else:
                        raise ImportError("Aucune fonction de chargement LayeredImage trouv?e dans kbmod")
                    # D?finir le temps (MJD souvent)
                    if hasattr(limg, "set_obstime"):
                        limg.set_obstime(Time(jd, format="jd").mjd)
                finally:
                    Path(tmp_path).unlink(missing_ok=True)
            imstack.add_image(limg)
            times_used.append(jd)
        except Exception as e:
            logger.warning(f"LayeredImage pour {fpath.name}: {e}")
            continue

    if imstack.img_count() == 0:
        logger.warning("Aucune image valide pour KBMOD")
        return []

    # Vitesse en px/jour depuis arcsec/min: v_px_day = v_arcmin_min / scale * 60 * 24
    scale = max(1e-6, scale_arcsec_per_px)
    v_min_px_day = velocity_arcsec_per_min_min / scale * 60.0 * 24.0
    v_max_px_day = velocity_arcsec_per_min_max / scale * 60.0 * 24.0

    if progress_callback:
        progress_callback(len(fits_paths), len(fits_paths) + 2)

    # Recherche
    try:
        search = _kbmod_search.StackSearch(imstack)
        # Bornes de vitesse (vx, vy en px/day). Recherche sur tout le champ.
        if hasattr(search, "search"):
            search.search(
                v_min_px_day, v_max_px_day,
                -v_max_px_day, v_max_px_day,  # vy
                max_results,
            )
        elif hasattr(search, "run_search"):
            # API alternative: param?tres via config
            cfg = getattr(_kbmod_search, "SearchConfiguration", None)
            if cfg:
                config = cfg()
                if hasattr(config, "min_velocity"):
                    config.min_velocity = v_min_px_day
                if hasattr(config, "max_velocity"):
                    config.max_velocity = v_max_px_day
                search.run_search(config)
            else:
                search.run_search(v_min_px_day, v_max_px_day, max_results)
        else:
            logger.warning("StackSearch sans m?thode search/run_search connue")
            return []

        if progress_callback:
            progress_callback(len(fits_paths) + 1, len(fits_paths) + 2)

        # R?cup?rer les r?sultats (trajectoires)
        if hasattr(search, "get_results"):
            raw_results = search.get_results(max_results)
        else:
            raw_results = []

        jd_ref = times_used[0] if times_used else 0.0
        results = []
        for r in raw_results[:max_results]:
            try:
                if hasattr(r, "x") and hasattr(r, "y"):
                    x0, y0 = float(r.x), float(r.y)
                elif hasattr(r, "trajectory"):
                    tr = r.trajectory
                    x0 = float(getattr(tr, "x", 0))
                    y0 = float(getattr(tr, "y", 0))
                else:
                    x0, y0 = 0.0, 0.0
                if hasattr(r, "vx") and hasattr(r, "vy"):
                    vx = float(r.vx)
                    vy = float(r.vy)
                elif hasattr(r, "trajectory"):
                    tr = r.trajectory
                    vx = float(getattr(tr, "vx", 0))
                    vy = float(getattr(tr, "vy", 0))
                else:
                    vx = vy = 0.0
                lh = float(getattr(r, "lh", getattr(r, "likelihood", 0.0)))
                try:
                    sky = wcs_ref.pixel_to_world(x0, y0)
                    ra_deg = float(sky.ra.deg)
                    dec_deg = float(sky.dec.deg)
                except Exception:
                    ra_deg = dec_deg = np.nan
                results.append({
                    "ra_deg": ra_deg,
                    "dec_deg": dec_deg,
                    "x0": x0,
                    "y0": y0,
                    "vx_px_per_day": vx,
                    "vy_px_per_day": vy,
                    "likelihood": lh,
                    "jd_ref": jd_ref,
                })
            except Exception as e:
                logger.debug(f"Skip result: {e}")
                continue

        return results

    except Exception as e:
        logger.error(f"Erreur KBMOD run_search: {e}", exc_info=True)
        raise RuntimeError(f"KBMOD search failed: {e}") from e
