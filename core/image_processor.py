# core/image_processor.py
import dataclasses
import os
import logging
from pathlib import Path
from astropy.io import fits
import numpy as np
from astropy.coordinates import (
    EarthLocation,
    SkyCoord,
    solar_system_ephemeris,
)
from astropy.time import Time
import astropy.units as u
from config import OBSERVATORY
from core.calculate_master import CalculateMaster
from core.astrometry import (
    AstrometrySolverNova,
    AstrometrySolverLocal,
    SolverConfig,
    run_plate_solve_parallel,
)
from utils.progress_manager import ProgressManager
import shutil
import threading

logger = logging.getLogger(__name__)


class PipelineControl:
    """Contrôle pause / reprise / arrêt entre étapes (thread-safe)."""

    def __init__(self):
        self._lock = threading.Lock()
        self._pause_event = threading.Event()
        self._pause_event.set()
        self._stop_requested = False

    def pause(self) -> None:
        with self._lock:
            self._pause_event.clear()

    def resume(self) -> None:
        with self._lock:
            self._pause_event.set()

    def stop(self) -> None:
        with self._lock:
            self._stop_requested = True

    def reset(self) -> None:
        with self._lock:
            self._stop_requested = False
            self._pause_event.set()

    def wait_if_paused(self) -> None:
        self._pause_event.wait()

    def should_stop(self) -> bool:
        with self._lock:
            return self._stop_requested


def _wcs_pixels_aligned(ref_wcs, img_wcs, shape, tol_arcsec: float = 1.0) -> bool:
    """True si les 4 coins ont la même position ciel (~tol) entre les deux WCS."""
    from astropy.coordinates import SkyCoord

    h, w = int(shape[0]), int(shape[1])
    if w < 2 or h < 2:
        return False
    corners = [
        (1.0, 1.0),
        (float(w - 1), 1.0),
        (1.0, float(h - 1)),
        (float(w - 1), float(h - 1)),
    ]
    for x, y in corners:
        try:
            s0 = ref_wcs.pixel_to_world(x, y)
            s1 = img_wcs.pixel_to_world(x, y)
            if s0.separation(s1).arcsec > tol_arcsec:
                return False
        except Exception:
            return False
    return True


def _get_binning(header) -> tuple[int, int]:
    """
    Extrait le binning (x, y) depuis l'en-tête FITS.
    Cherche XBINNING/YBINNING, BINAXIS1/BINAXIS2, ou CCDSUM (ex: "2x2").
    Retourne (1, 1) si non trouvé.
    """
    h = header
    # XBINNING / YBINNING (MaxIm DL, etc.)
    x = h.get("XBINNING")
    y = h.get("YBINNING")
    if x is not None and y is not None:
        try:
            return (int(float(x)), int(float(y)))
        except (TypeError, ValueError):
            pass
    # BINAXIS1 / BINAXIS2
    x = h.get("BINAXIS1")
    y = h.get("BINAXIS2")
    if x is not None and y is not None:
        try:
            return (int(float(x)), int(float(y)))
        except (TypeError, ValueError):
            pass
    # CCDSUM "1x1" ou "2x2"
    ccdsum = h.get("CCDSUM", "")
    if isinstance(ccdsum, str) and "x" in ccdsum.lower():
        try:
            parts = ccdsum.strip().lower().split("x")
            if len(parts) >= 2:
                return (int(float(parts[0].strip())), int(float(parts[1].strip())))
        except (TypeError, ValueError):
            pass
    return (1, 1)


def _get_filter(header) -> str:
    """
    Extrait le nom du filtre depuis l'en-tête FITS (FILTER, FILT, etc.).
    Retourne une chaîne normalisée (strip, vide si absent).
    """
    for key in ("FILTER", "FILT", "FILTER1", "INSFLTNM"):
        val = header.get(key)
        if val is not None and str(val).strip():
            return str(val).strip()
    return ""


def _get_exptime(header) -> float:
    """
    Extrait le temps d'exposition en secondes (EXPTIME, EXPOSURE).
    Retourne 0.0 si absent ou invalide.
    """
    for key in ("EXPTIME", "EXPOSURE"):
        val = header.get(key)
        if val is not None:
            try:
                return float(val)
            except (TypeError, ValueError):
                pass
    return 0.0


def _build_master_header(first_file, data_shape, exptime_override=None, obstype=None):
    """
    Construit un en-tête FITS minimal pour un master (bias/dark/flat)
    afin qu'il soit réutilisable (binning, date, temps d'exposition).
    """
    hdr_src = fits.getheader(first_file, 0)
    binx, biny = _get_binning(hdr_src)
    exptime = _get_exptime(hdr_src) if exptime_override is None else exptime_override
    date_obs = hdr_src.get("DATE-OBS") or hdr_src.get("DATE")
    if not date_obs:
        date_obs = Time.now().to_value("fits", subfmt="date")
    ny, nx = data_shape[0], data_shape[1]
    new_hdr = fits.Header()
    new_hdr["SIMPLE"] = True
    new_hdr["BITPIX"] = -32
    new_hdr["NAXIS"] = 2
    new_hdr["NAXIS1"] = nx
    new_hdr["NAXIS2"] = ny
    new_hdr["XBINNING"] = binx
    new_hdr["YBINNING"] = biny
    new_hdr["DATE-OBS"] = str(date_obs).strip()
    new_hdr["EXPTIME"] = exptime
    if obstype:
        new_hdr["OBSTYPE"] = obstype
    return new_hdr


def _check_binning_consistency(science_files, bias_files, dark_files, flat_files) -> None:
    """
    Vérifie que toutes les images de calibration (bias, darks, flats) ont le même
    binning que les images science. Lève ValueError avec message explicite si non.
    """
    if not science_files:
        return
    ref_header = fits.getheader(science_files[0], 0)
    ref_bin = _get_binning(ref_header)

    def check_set(files, label):
        for path in (files or []):
            try:
                h = fits.getheader(path, 0)
                bin_here = _get_binning(h)
                # Anciens masters sans en-tête sont lus en 1x1 : les accepter si réutilisés (on suppose même binning que les lights)
                is_master = (
                    path.name.lower().startswith("master_")
                    or str(h.get("OBSTYPE", "")).upper().startswith("MASTER")
                )
                if is_master and bin_here == (1, 1) and ref_bin != (1, 1):
                    bin_here = ref_bin
                if bin_here != ref_bin:
                    raise ValueError(
                        f"Binning incohérent : les images {label} doivent avoir le même binning que les lights.\n"
                        f"  Lights (référence) : {ref_bin[0]}x{ref_bin[1]} | "
                        f"{path.name} : {bin_here[0]}x{bin_here[1]}."
                    )
            except Exception as e:
                if isinstance(e, ValueError):
                    raise
                logger.warning(f"Impossible de lire le binning de {path}: {e}")

    check_set(bias_files, "bias")
    check_set(dark_files, "darks")
    check_set(flat_files, "flats")
    logger.info(f"Vérification binning OK (référence lights : {ref_bin[0]}x{ref_bin[1]}).")


def _check_filter_consistency(science_files, flat_files) -> str | None:
    """
    Vérifie la cohérence filtres lights / flats.
    En cas d'incohérence : retourne un message d'avertissement (pas d'exception).
    """
    if not science_files or not flat_files:
        return None
    # Filtres présents dans les lights (ensemble)
    light_filters = set()
    for path in science_files:
        try:
            h = fits.getheader(path, 0)
            f = _get_filter(h)
            if f:
                light_filters.add(f)
        except Exception as e:
            logger.warning(f"Impossible de lire FILTER de {path}: {e}")
    # Filtre(s) des flats
    flat_filters = set()
    for path in flat_files:
        try:
            h = fits.getheader(path, 0)
            f = _get_filter(h)
            if f:
                flat_filters.add(f)
        except Exception as e:
            logger.warning(f"Impossible de lire FILTER de {path}: {e}")

    if not flat_filters:
        logger.warning("Aucune clé FILTER trouvée dans les flats ; vérification filtre ignorée.")
        return None
    if not light_filters:
        logger.warning("Aucune clé FILTER trouvée dans les lights ; vérification filtre ignorée.")
        return None

    if not flat_filters.issubset(light_filters):
        msg = (
            "Filtre incohérent entre lights et flats : le master flat sera quand même appliqué. "
            f"Filtres lights : {sorted(light_filters) or '(aucun)'} | "
            f"filtres flats : {sorted(flat_filters)}."
        )
        logger.warning(msg)
        return f"⚠️ {msg}"

    if len(light_filters) > 1 and len(flat_filters) == 1:
        logger.warning(
            "Les lights contiennent plusieurs filtres (%s) alors que les flats n'en ont qu'un (%s). "
            "Le même master flat sera appliqué à toutes les images.",
            sorted(light_filters), next(iter(flat_filters)),
        )
    logger.info("Vérification filtre lights/flats OK.")
    return None


def _gain_readnoise_electrons(header) -> tuple[float, float]:
    """Gain (e⁻/ADU) et bruit de lecture (e⁻/pixel) depuis l'en-tête FITS."""
    g = header.get("GAIN", header.get("EGAIN", 1.0))
    try:
        gain = float(g) if g is not None else 1.0
    except (TypeError, ValueError):
        gain = 1.0
    gain = max(gain, 1e-6)
    rn = header.get("RDNOISE", header.get("READNOIS", header.get("READNOISE", 0.0)))
    try:
        readnoise = float(rn) if rn is not None else 0.0
    except (TypeError, ValueError):
        readnoise = 0.0
    readnoise = max(readnoise, 0.0)
    return gain, readnoise


def _variance_poisson_readnoise_adu(
    data_adu: np.ndarray, gain: float, readnoise_e: float
) -> np.ndarray:
    """
    Variance (ADU**2) modèle Poisson + bruit de lecture : max(I,0)/G + (RN_e/G)**2.
    ``readnoise_e`` en électrons par pixel (convention FITS courante).
    """
    g = max(float(gain), 1e-6)
    rn_adu = max(float(readnoise_e), 0.0) / g
    d = np.asarray(data_adu, dtype=float)
    return np.maximum(d, 0.0) / g + rn_adu**2


def _stack_median_and_variance(stack: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """
    Médiane et variance d'échantillon (ddof=1) par pixel sur une pile 2D.
    Une seule image → variance nulle.
    """
    if not stack:
        raise ValueError("stack vide")
    arr = np.stack(stack, axis=0).astype(float)
    med = np.median(arr, axis=0)
    if arr.shape[0] < 2:
        var = np.zeros_like(med, dtype=float)
    else:
        var = np.var(arr, axis=0, ddof=1)
    return med, var


def _write_calibrated_with_variance(
    out_path: Path,
    data: np.ndarray,
    variance: np.ndarray,
    header: fits.Header,
) -> None:
    """Écrit image calibrée (HDU0) + carte de variance par pixel (HDU1 ``VARIANCE``)."""
    data = np.asarray(data, dtype=np.float32)
    variance = np.asarray(variance, dtype=np.float32)
    variance = np.clip(variance, 0.0, np.finfo(np.float32).max)
    phdu = fits.PrimaryHDU(data=data, header=header)
    vhdr = fits.Header()
    vhdr["EXTNAME"] = "VARIANCE"
    vhdr["BUNIT"] = ("ADU**2", "Variance (science units squared)")
    vhdr["COMMENT"] = "Variance ADU**2 (Poisson+RN + masters + division flat, numpy)"
    ihdu = fits.ImageHDU(data=variance, header=vhdr)
    fits.HDUList([phdu, ihdu]).writeto(str(out_path), overwrite=True)


class ImageProcessor:
    def __init__(self, base_dir: str | Path, bash_path: str | None = None):
        self.base_dir = Path(base_dir)
        self.science_dir = self.base_dir / "science"
        self.science_dir.mkdir(parents=True, exist_ok=True)
        self.observatory = OBSERVATORY
        self.output_dir = self.base_dir / "output"
        self.calibrated_dir = self.output_dir / "calibrated"
        self.astrometry_dir = self.output_dir / "astrometry"
        self.aligned_dir = self.science_dir / "aligned"
        for d in (self.calibrated_dir, self.astrometry_dir, self.aligned_dir):
            d.mkdir(parents=True, exist_ok=True)

        # --- CORRECTION ICI : ON FORCE WSL ---
        # On ne cherche plus cygwin. On assume que WSL est installé.
        self.bash_path = "wsl" 
        logger.info(f"Environnement configuré pour WSL (Ubuntu)")

        # Position par défaut de l'observatoire (config.py)
        # ... (le reste de l'init ne change pas) ...
        lat = self.observatory["lat"]
        lon = self.observatory["lon"]
        elev = self.observatory["elev"]
        self.location = EarthLocation(lat=float(lat) * u.deg,
                              lon=float(lon) * u.deg,
                              height=float(elev) * u.m)

    def process_calibration(
        self,
        science_files,
        bias_files=None,
        dark_files=None,
        flat_files=None,
        progress_callback=None,
        scale_darks=False,
    ):
        """
        Calibration bias/dark/flat avec contrôles de cohérence (binning, filtre)
        et option de scaling des darks au temps d'exposition des lights (type AstroImageJ).

        Réduction **numpy** (médiane des masters, soustractions / division comme avant).
        Carte de variance par pixel (HDU1 ``VARIANCE``, ADU**2) : Poisson + bruit de
        lecture (en-tête GAIN / RDNOISE), variance d'empilement des masters, puis
        propagation pour soustraction du dark scalé et division par le flat.

        Returns
        -------
        str | None
            Message d'avertissement filtres lights/flats si incohérence (calibration poursuivie), sinon None.
        """
        logging.info("Début de la calibration (ImageProcessor.process_calibration)")

        science_files = [Path(f) for f in (science_files or [])]
        bias_files = [Path(f) for f in (bias_files or [])]
        dark_files = [Path(f) for f in (dark_files or [])]
        flat_files = [Path(f) for f in (flat_files or [])]

        if not science_files:
            logging.warning("Aucun fichier science fourni.")
            if progress_callback:
                progress_callback(0)
            return None

        # ─── Sécurités : binning et filtre (filtre : avertissement seulement) ─
        _check_binning_consistency(science_files, bias_files, dark_files, flat_files)
        filter_warn = _check_filter_consistency(science_files, flat_files)

        if progress_callback:
            progress_callback(0)

        def _parse_date_obs(date_str):
            date_str = str(date_str).strip()
            if not date_str:
                return None
            for fmt in ("isot", "fits"):
                try:
                    return Time(date_str, format=fmt, scale="utc")
                except Exception:
                    continue
            logging.error(f"DATE-OBS illisible : {date_str}")
            return None

        def _dark_exptime_paths(files):
            if not files:
                return 0.0
            times = []
            for fp in files:
                try:
                    t = _get_exptime(fits.getheader(fp, 0))
                    if t > 0:
                        times.append(t)
                except Exception:
                    pass
            return float(np.median(times)) if times else 0.0

        def build_master_stack(
            files, bias=None, dark=None, dark_exptime_ref=0.0, scale_dark_by_exptime=False
        ):
            """Retourne (médiane, variance d'échantillon par pixel) ou (None, None)."""
            stack = []
            for fp in files:
                try:
                    data = fits.getdata(fp, 0).astype(float)
                    h = fits.getheader(fp, 0)
                    if bias is not None:
                        data = data - bias
                    if dark is not None:
                        dark_to_subtract = dark
                        if scale_dark_by_exptime and dark_exptime_ref > 0:
                            frame_exptime = _get_exptime(h)
                            if frame_exptime > 0:
                                dark_to_subtract = dark * (frame_exptime / dark_exptime_ref)
                            else:
                                logger.warning(
                                    "%s : EXPTIME absent/nul, dark non scalé pour construction master.",
                                    Path(fp).name,
                                )
                        data = data - dark_to_subtract
                    stack.append(data)
                except Exception as e:
                    logging.error("Erreur lecture %s : %s", fp, e)
            if not stack:
                return None, None
            return _stack_median_and_variance(stack)

        dark_exptime = _dark_exptime_paths(dark_files) if dark_files else 0.0
        if scale_darks and dark_files and dark_exptime <= 0:
            logging.warning(
                "Scale darks activé mais EXPTIME des darks introuvable ; scaling désactivé pour les darks."
            )
            scale_darks = False
        if scale_darks and dark_exptime > 0:
            logger.info("Scaling des darks activé (temps de référence dark : %.2f s).", dark_exptime)

        master_bias, master_bias_var = build_master_stack(bias_files)
        master_dark, master_dark_var = build_master_stack(
            dark_files, bias=master_bias, dark=None
        )
        master_flat, master_flat_var = build_master_stack(
            flat_files,
            bias=master_bias,
            dark=master_dark,
            dark_exptime_ref=dark_exptime,
            scale_dark_by_exptime=bool(scale_darks),
        )

        if master_flat is not None:
            mean_flat = np.mean(master_flat[master_flat > 0])
            if mean_flat > 0:
                master_flat = master_flat / mean_flat
                if master_flat_var is not None:
                    master_flat_var = master_flat_var / (mean_flat**2)

        try:
            if master_bias is not None and bias_files:
                hdr_bias = _build_master_header(
                    bias_files[0], master_bias.shape, exptime_override=0.0, obstype="MASTER_BIAS"
                )
                fits.writeto(self.output_dir / "master_bias.fits", master_bias, hdr_bias, overwrite=True)
            if master_dark is not None and dark_files:
                hdr_dark = _build_master_header(
                    dark_files[0], master_dark.shape, exptime_override=dark_exptime, obstype="MASTER_DARK"
                )
                fits.writeto(self.output_dir / "master_dark.fits", master_dark, hdr_dark, overwrite=True)
            if master_flat is not None and flat_files:
                hdr_flat = _build_master_header(
                    flat_files[0], master_flat.shape, exptime_override=None, obstype="MASTER_FLAT"
                )
                fits.writeto(self.output_dir / "master_flat.fits", master_flat, hdr_flat, overwrite=True)
        except Exception as e:
            logging.warning(f"Erreur sauvegarde masters : {e}")

        total = len(science_files)
        failed_files = []
        eps = 1e-12
        for i, f in enumerate(science_files, 1):
            try:
                data = fits.getdata(f, 0).astype(float)
                hdr = fits.getheader(f, 0)
                gain, rn_e = _gain_readnoise_electrons(hdr)
                var = _variance_poisson_readnoise_adu(data, gain, rn_e)

                if master_bias is not None:
                    data = data - master_bias
                    if master_bias_var is not None:
                        var = var + master_bias_var
                if master_dark is not None:
                    if scale_darks and dark_exptime > 0:
                        science_exptime = _get_exptime(hdr)
                        if science_exptime > 0:
                            scale = science_exptime / dark_exptime
                            data = data - master_dark * scale
                            if master_dark_var is not None:
                                var = var + (scale**2) * master_dark_var
                        else:
                            data = data - master_dark
                            logging.warning(f"{f.name} : EXPTIME absent ou nul, dark non scalé.")
                            if master_dark_var is not None:
                                var = var + master_dark_var
                    else:
                        data = data - master_dark
                        if master_dark_var is not None:
                            var = var + master_dark_var
                if master_flat is not None:
                    T = data.astype(float)
                    var_T = var.astype(float)
                    F = master_flat.astype(float)
                    Vf = (
                        master_flat_var.astype(float)
                        if master_flat_var is not None
                        else np.zeros_like(F)
                    )
                    mask = F != 0
                    data = data.copy()
                    data[mask] = T[mask] / F[mask]
                    rel = var_T / (T * T + eps) + Vf / (F * F + eps)
                    var_out = var_T.copy()
                    var_out[mask] = (data[mask] ** 2) * rel[mask]
                    var = var_out

                hdr["HISTORY"] = (
                    "Calibrated (bias/dark/flat + VARIANCE HDU, numpy) by ImageProcessor"
                )

                date_obs_str = hdr.get("DATE-OBS")
                if date_obs_str:
                    date_obs = _parse_date_obs(date_obs_str)
                    if date_obs:
                        exptime = float(hdr.get("EXPTIME", 0))
                        midpoint = Time(
                            date_obs.utc + (exptime / 2.0) * u.s,
                            scale="utc",
                            location=self.location,
                        )
                        hdr["JD-UTC"] = (midpoint.jd, "Julian Date (UTC) at mid-exposure")
                        target = None
                        try:
                            ra_str = hdr.get("OBJCTRA")
                            dec_str = hdr.get("OBJCTDEC")
                            if ra_str and dec_str:
                                target = SkyCoord(
                                    ra_str,
                                    dec_str,
                                    unit=(u.hourangle, u.deg),
                                    frame="icrs",
                                )
                        except Exception as e:
                            logging.debug(f"{f.name} : erreur parsing RA/DEC : {e}")
                        if target:
                            try:
                                midpoint_tdb = midpoint.tdb
                                with solar_system_ephemeris.set("builtin"):
                                    ltt = midpoint_tdb.light_travel_time(target)
                                bjd_tdb = midpoint_tdb + ltt
                                hdr["BJD-TDB"] = (bjd_tdb.jd, "Barycentric Julian Date (TDB)")
                            except Exception as e_bjd:
                                logging.warning(f"Échec calcul BJD-TDB pour {f.name}: {e_bjd}")

                out_path = self.calibrated_dir / f.name
                _write_calibrated_with_variance(out_path, data, var, hdr)
                logging.info("Image calibrée + VARIANCE : %s", out_path)

            except Exception as e:
                logging.error(f"Erreur calibration {f.name} : {e}")
                failed_files.append(f.name)

            if progress_callback:
                progress_callback(min(100, int(i / total * 100)))

        if failed_files:
            failed_count = len(failed_files)
            fail_ratio = failed_count / float(total) if total else 1.0
            preview = ", ".join(failed_files[:5])
            if failed_count > 5:
                preview += ", ..."
            summary = (
                f"Calibration terminée avec {failed_count} échec(s) sur {total} image(s). "
                f"Fichiers en échec: {preview}"
            )
            # Seuil bloquant: au-delà de 20% d'échecs, on force un arrêt explicite.
            if fail_ratio > 0.20:
                logging.error(summary)
                raise RuntimeError(summary)
            logging.warning(summary)
        else:
            logging.info("Calibration terminée avec succès !")

        return filter_warn

    def process_astrometry(
        self,
        method: str = "NOVA",
        progress_callback=None,
        pipeline_control: PipelineControl | None = None,
    ):
        """
        Astrométrie des images calibrées.

        - Les images d'entrée sont prises dans self.calibrated_dir
        - Pour la méthode LOCAL :
            * solve-field travaille sur les images calibrées
            * les images astrométrées finales sont enregistrées dans self.science_dir
        - La progression est renvoyée via progress_callback (0–100).
        - Plusieurs images sont résolues en parallèle (``SolverConfig.parallel_workers``),
          avec OMP_NUM_THREADS réduit par solve-field pour limiter la surcharge CPU.
        """
        calibrated_dir = self.calibrated_dir
        calibrated_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"📂 Astrométrie sur le dossier calibré : {calibrated_dir}")

        calibrated_files = sorted(calibrated_dir.glob("*.fits"))
        if not calibrated_files:
            raise ValueError("Aucune image calibrée. Impossible de lancer l'astrométrie.")

        solve_cfg = SolverConfig()
        total = len(calibrated_files)
        pw = max(1, min(solve_cfg.parallel_workers, total))

        # ----------------------- choix du solveur -----------------------
        method = method.upper()

        if method == "NOVA":
            solver = AstrometrySolverNova(
                output_dir=self.astrometry_dir,
                parallel_workers=pw,
            )

        elif method == "LOCAL":
            omp_each = max(1, solve_cfg.omp_threads // pw)
            cfg_adj = dataclasses.replace(solve_cfg, omp_threads=omp_each)
            solver = AstrometrySolverLocal(science_dir=self.science_dir, config=cfg_adj)
            logger.info(
                "Astrométrie locale (WSL) : jusqu'à %d images en parallèle, "
                "OMP_NUM_THREADS=%d par solve-field.",
                pw,
                omp_each,
            )
        else:
            raise ValueError(
                "Méthode d'astrométrie non supportée : utilisez uniquement 'NOVA' ou 'LOCAL' (WSL / solve-field)."
            )

        if progress_callback:
            progress_callback(0)

        if method == "NOVA":
            solver.solve_directory(
                calibrated_dir,
                progress_callback=progress_callback,
                pipeline_control=pipeline_control,
            )
            if pipeline_control and pipeline_control.should_stop():
                logger.info("Astrométrie NOVA interrompue (Stop).")
            else:
                logger.info("Astrométrie NOVA terminée.")
            return

        run_plate_solve_parallel(
            calibrated_files,
            lambda p: solver.solve_file(p, None),
            max_workers=pw,
            progress_callback=progress_callback,
            pipeline_control=pipeline_control,
            log_prefix="LOCAL solve-field",
        )
        if pipeline_control and pipeline_control.should_stop():
            logger.info("Astrométrie locale interrompue (Stop).")
        else:
            logger.info("Astrométrie LOCAL terminée.")
        
    def process_alignment_wcs(
        self,
        input_dir=None,
        output_dir=None,
        progress_callback=None,
        pipeline_control: PipelineControl | None = None,
    ):
        """
        Aligne les images FITS dans input_dir en utilisant les informations WCS.
        Les images alignées/reprojetées sont sauvegardées dans output_dir.
        """
        input_dir = input_dir if input_dir else self.science_dir
        output_dir = output_dir if output_dir else self.aligned_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        from astropy.wcs import WCS
        try:
            from reproject import reproject_interp
        except ImportError:
            logger.error("Module 'reproject' manquant. Installez-le avec 'pip install reproject'")
            raise

        input_files = sorted(input_dir.glob("*.fits"))
        if not input_files:
            raise ValueError(f"Aucune image FITS résolue trouvée dans {input_dir}.")

        ref_path = input_files[0]
        with fits.open(ref_path, memmap=True) as hdul:
            ref_wcs = WCS(hdul[0].header)
            ref_shape = hdul[0].data.shape

        logger.info(f"Début de l'alignement WCS sur {len(input_files)} fichiers.")

        total = len(input_files)

        for i, f_path in enumerate(input_files):
            if pipeline_control:
                pipeline_control.wait_if_paused()
                if pipeline_control.should_stop():
                    logger.info("Alignement WCS interrompu (Stop).")
                    return

            try:
                with fits.open(f_path, memmap=True) as hdul:
                    raw = hdul[0].data
                    hdr = hdul[0].header
                    current_wcs = WCS(hdr)

                if i == 0:
                    new_data = np.asarray(raw, dtype=float)
                    new_hdr = hdr.copy()
                    new_hdr["HISTORY"] = "Aligned (reference) by ImageProcessor"
                else:
                    if _wcs_pixels_aligned(ref_wcs, current_wcs, raw.shape):
                        new_data = np.asarray(raw, dtype=float)
                        new_hdr = hdr.copy()
                        ref_wcs_hdr = ref_wcs.to_header()
                        for key in ref_wcs_hdr:
                            try:
                                new_hdr[key] = ref_wcs_hdr[key]
                            except Exception:
                                pass
                        new_hdr["NAXIS"] = len(ref_shape)
                        if len(ref_shape) >= 1:
                            new_hdr["NAXIS1"] = ref_shape[1]
                        if len(ref_shape) >= 2:
                            new_hdr["NAXIS2"] = ref_shape[0]
                        new_hdr["HISTORY"] = "Aligned (WCS corners match ref, no reproject) by ImageProcessor"
                    else:
                        arr_in = np.asarray(raw, dtype=float)
                        try:
                            new_data, _footprint = reproject_interp(
                                (arr_in, current_wcs),
                                ref_wcs,
                                shape_out=ref_shape,
                                order="bilinear",
                                parallel=True,
                            )
                        except TypeError:
                            new_data, _footprint = reproject_interp(
                                (arr_in, current_wcs),
                                ref_wcs,
                                shape_out=ref_shape,
                                order="bilinear",
                            )
                        new_hdr = hdr.copy()
                        ref_wcs_hdr = ref_wcs.to_header()
                        for key in ref_wcs_hdr:
                            try:
                                new_hdr[key] = ref_wcs_hdr[key]
                            except Exception:
                                pass
                        new_hdr["NAXIS"] = len(ref_shape)
                        if len(ref_shape) >= 1:
                            new_hdr["NAXIS1"] = ref_shape[1]
                        if len(ref_shape) >= 2:
                            new_hdr["NAXIS2"] = ref_shape[0]
                        new_hdr["HISTORY"] = "Aligned (WCS reprojected) by ImageProcessor"

                out_path = output_dir / f_path.name
                fits.writeto(out_path, new_data, new_hdr, overwrite=True)
                logger.info(f"Image alignée (WCS): {out_path.name}")

            except Exception as e:
                logger.error(f"Erreur alignment {f_path.name}: {e}")

            if progress_callback:
                progress_callback(min(100, int((i + 1) / total * 100)))

        logger.info("Alignement WCS terminé.")

    def process_stacking(
        self,
        input_files: list[Path],
        output_path: Path,
        progress_callback=None,
        pipeline_control: PipelineControl | None = None,
    ):
        """
        Empile (median stacking) les images de la liste fournie.
        """
        if not input_files:
            raise ValueError("Aucun fichier d'entrée fourni pour l'empilement.")

        logger.info(f"Début de l'empilement (median stacking) sur {len(input_files)} fichiers.")

        stack = []
        total = len(input_files)

        for i, f_path in enumerate(input_files):
            if pipeline_control:
                pipeline_control.wait_if_paused()
                if pipeline_control.should_stop():
                    raise RuntimeError("Empilement interrompu par l'utilisateur.")
            try:
                try:
                    data = fits.getdata(f_path, 0, memmap=True)
                except TypeError:
                    data = fits.getdata(f_path, 0)
                stack.append(np.asarray(data, dtype=float))
            except Exception as e:
                logger.error(f"Erreur lecture pour stacking {f_path.name}: {e}")

            if progress_callback:
                progress_callback(min(90, int((i + 1) / total * 90)))

        if not stack:
            raise ValueError("Aucune image valide n'a pu être lue pour l'empilement.")

        median_stacked_data = np.median(np.array(stack, dtype=float), axis=0)
        
        # Copie de l'en-tête de la première image
        first_hdr = fits.getheader(input_files[0], 0)
        first_hdr['HISTORY'] = 'Median Stacked by ImageProcessor'
        first_hdr['N_STACK'] = (len(stack), 'Number of frames stacked')

        # Sauvegarde
        fits.writeto(output_path, median_stacked_data, first_hdr, overwrite=True)
        logger.info(f"Image empilée sauvegardée : {output_path}")

        if progress_callback:
            progress_callback(100)

        logger.info("Empilement terminé avec succès.")



