# core/exotethys_ldc.py
"""
Coefficients de limb darkening via **ExoTETHyS** pour pylightcurve (``add_filter``).

Référence : https://github.com/ucl-exoplanets/ExoTETHyS — Morello et al., AJ 160, 112 (2020).

``ldc_calculate`` peut appeler ``exit()`` : le calcul est exécuté dans un sous-processus
(``exotethys_ldc_worker.py``).
"""

from __future__ import annotations

import importlib.resources
import logging
import pickle
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

try:
    import h5py  # noqa: F401
    import exotethys  # noqa: F401
    import click  # noqa: F401

    EXOTETHYS_AVAILABLE = True
except ImportError:
    EXOTETHYS_AVAILABLE = False

EXOTETHYS_STELLAR_MODELS: tuple[str, ...] = (
    "Phoenix_2018",
    "Phoenix_2012_13",
    "Phoenix_drift_2012",
    "Atlas_2000",
    "Stagger_2015",
    "Stagger_2018",
    "MPS_Atlas_set1_2023",
    "MPS_Atlas_set2_2023",
)

_EXOTETHYS_BUILTIN_PASSBAND: Dict[str, str] = {
    "TESS": "TESS",
    "COUSINS_V": "Johnson_V",
    "COUSINS_R": "Johnson_R",
    "COUSINS_I": "Johnson_R",
}

_PYLC_PASS_ALIASES: Dict[str, str] = {
    "SLOAN_g": "sdss_g",
    "SLOAN_r": "sdss_r",
    "SLOAN_i": "sdss_i",
}

_NPOAP_EXO_TARGET = "npoap_ld"
_DEFAULT_TIMEOUT_SEC = 400


def _exotethys_passbands_dir() -> Path:
    return Path(str(importlib.resources.files("exotethys") / "Passbands"))


def _worker_script_path() -> Path:
    return Path(__file__).resolve().parent / "exotethys_ldc_worker.py"


def ensure_passband_for_pylightcurve(filter_ui: str) -> None:
    """Garantit ``filter_ui.pass`` dans le répertoire photometry pylightcurve."""
    from pylightcurve.__databases__ import plc_data

    dest = Path(plc_data.photometry()) / f"{filter_ui}.pass"
    if dest.is_file():
        return

    dest.parent.mkdir(parents=True, exist_ok=True)
    builtin = _EXOTETHYS_BUILTIN_PASSBAND.get(filter_ui)
    if builtin is not None:
        src = _exotethys_passbands_dir() / f"{builtin}.pass"
        if src.is_file():
            shutil.copy2(src, dest)
            logger.info("Passband %s : copie depuis ExoTETHyS (%s.pass).", filter_ui, builtin)
            plc_data.all_filters_data = None
            return

    plc_base = _PYLC_PASS_ALIASES.get(filter_ui, filter_ui)
    src_plc = Path(plc_data.photometry()) / f"{plc_base}.pass"
    if src_plc.is_file():
        shutil.copy2(src_plc, dest)
        logger.info("Passband %s : copie depuis pylightcurve (%s.pass).", filter_ui, plc_base)
        plc_data.all_filters_data = None
        return

    if filter_ui in ("G", "G_BP", "G_RP"):
        from core.gaia_pylightcurve_support import ensure_gaia_passbands_installed

        ensure_gaia_passbands_installed()
        return

    raise FileNotFoundError(
        f"Aucun fichier passband pour « {filter_ui} ». Utilisez un filtre reconnu ou ajoutez un .pass."
    )


def _copy_passband_for_exotethys_run(filter_ui: str, work_dir: Path) -> None:
    """
    Copie la courbe pour ExoTETHyS dans ``work_dir/{filter_ui}`` (sans extension).

    Si ``passbands_path`` est fourni dans le cfg, ExoTETHyS force ``passbands_ext`` à
    une chaîne vide : le fichier lu est ``join(path, passband)``, pas ``… .pass``.
    Le mot-clé ``passbands_ext`` n'est pas autorisé dans le fichier de config.
    """
    dest = work_dir / filter_ui
    if dest.is_file():
        return

    if filter_ui in ("TESS",) or filter_ui in _EXOTETHYS_BUILTIN_PASSBAND:
        builtin = _EXOTETHYS_BUILTIN_PASSBAND.get(filter_ui, filter_ui)
        src = _exotethys_passbands_dir() / f"{builtin}.pass"
        if not src.is_file():
            raise FileNotFoundError(f"Passband ExoTETHyS introuvable : {builtin}.pass")
        shutil.copy2(src, dest)
        return

    ensure_passband_for_pylightcurve(filter_ui)
    from pylightcurve.__databases__ import plc_data

    src = Path(plc_data.photometry()) / f"{filter_ui}.pass"
    if not src.is_file():
        raise FileNotFoundError(src)
    shutil.copy2(src, dest)


def run_exotethys_ldc_claret4(
    filter_ui: str,
    teff: float,
    logg: float,
    mh: float,
    stellar_model: str = "Phoenix_2018",
    timeout_sec: int = _DEFAULT_TIMEOUT_SEC,
) -> np.ndarray:
    """
    Lance ExoTETHyS (calculation_type individual, loi claret4) et retourne les 4 coefficients.
    """
    if not EXOTETHYS_AVAILABLE:
        raise RuntimeError("ExoTETHyS n'est pas installé (pip install exotethys h5py click).")
    if stellar_model not in EXOTETHYS_STELLAR_MODELS:
        raise ValueError(f"Modèle stellaire invalide : {stellar_model!r}")

    work_dir = Path(tempfile.mkdtemp(prefix="npoap_exotethys_"))
    out_dir = work_dir / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        _copy_passband_for_exotethys_run(filter_ui, work_dir)
        wp = str(work_dir.resolve()).replace("\\", "/")
        op = str(out_dir.resolve()).replace("\\", "/")
        cfg_path = work_dir / "cfg.txt"
        cfg_path.write_text(
            textwrap.dedent(
                f"""
                calculation_type individual
                stellar_models_grid {stellar_model}
                limb_darkening_laws claret4
                passbands {filter_ui}
                passbands_path {wp}
                star_effective_temperature {float(teff)}
                star_log_gravity {float(logg)}
                star_metallicity {float(mh)}
                target_names {_NPOAP_EXO_TARGET}
                output_path {op}
                user_output basic
                """
            ).strip()
            + "\n",
            encoding="utf-8",
        )

        worker = _worker_script_path()
        if not worker.is_file():
            raise FileNotFoundError(str(worker))

        proc = subprocess.run(
            [sys.executable, str(worker), str(cfg_path)],
            capture_output=True,
            text=True,
            timeout=timeout_sec,
        )
        if proc.returncode != 0:
            msg = (proc.stderr or proc.stdout or "")[:4000]
            raise RuntimeError(f"ExoTETHyS a échoué (code {proc.returncode}) :\n{msg}")

        pkl = out_dir / f"{_NPOAP_EXO_TARGET}_ldc.pickle"
        if not pkl.is_file():
            msg = (proc.stderr or proc.stdout or "")[:4000]
            raise FileNotFoundError(
                f"Sortie ExoTETHyS absente : {pkl}. Sortie du worker :\n{msg}"
            )

        with open(pkl, "rb") as f:
            data = pickle.load(f)
        coeffs = data["passbands"][filter_ui]["claret4"]["coefficients"]
        return np.asarray(coeffs, dtype=float).ravel()
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def inject_exotethys_claret4_into_planet(planet, filter_ui: str, coeffs4: np.ndarray) -> None:
    """Enregistre le filtre sur ``Planet`` avec les 4 coefficients Claret et Fp/Fs pylightcurve."""
    from pylightcurve.models.exoplanet_lc import fp_over_fs

    rp = float(planet.rp_over_rs)
    fp = fp_over_fs(
        rp,
        planet.sma_over_rs,
        planet.albedo,
        planet.emissivity,
        planet.stellar_temperature,
        filter_ui,
    )
    c = np.asarray(coeffs4, dtype=float).ravel()
    out = np.zeros(4, dtype=float)
    out[: min(4, len(c))] = c[:4]
    planet.filters.pop(filter_ui, None)
    planet.add_filter(filter_ui, rp, float(out[0]), float(out[1]), float(out[2]), float(out[3]), fp)
    planet.ldc_method = "claret"
    logger.info(
        "LDC ExoTETHyS (claret4) pour %s : %.5f %.5f %.5f %.5f",
        filter_ui,
        out[0],
        out[1],
        out[2],
        out[3],
    )


def prepare_planet_exotethys_claret4(
    planet,
    filter_ui: str,
    teff: float,
    logg: float,
    mh: float,
    stellar_model: str = "Phoenix_2018",
) -> None:
    """Passeband pylightcurve + coefficients ExoTETHyS + ``ldc_method='claret'``."""
    ensure_passband_for_pylightcurve(filter_ui)
    coeffs = run_exotethys_ldc_claret4(
        filter_ui, teff, logg, mh, stellar_model=stellar_model
    )
    inject_exotethys_claret4_into_planet(planet, filter_ui, coeffs)


def run_exotethys_ldc_power2(
    filter_ui: str,
    teff: float,
    logg: float,
    mh: float,
    stellar_model: str = "Phoenix_2018",
    timeout_sec: int = _DEFAULT_TIMEOUT_SEC,
) -> np.ndarray:
    """Retourne ``[c, alpha]`` (loi power-2 ExoTETHyS, format original)."""
    if not EXOTETHYS_AVAILABLE:
        raise RuntimeError("ExoTETHyS n'est pas installé (pip install exotethys h5py click).")
    if stellar_model not in EXOTETHYS_STELLAR_MODELS:
        raise ValueError(f"Modèle stellaire invalide : {stellar_model!r}")

    work_dir = Path(tempfile.mkdtemp(prefix="npoap_exotethys_p2_"))
    out_dir = work_dir / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        _copy_passband_for_exotethys_run(filter_ui, work_dir)
        wp = str(work_dir.resolve()).replace("\\", "/")
        op = str(out_dir.resolve()).replace("\\", "/")
        cfg_path = work_dir / "cfg.txt"
        cfg_path.write_text(
            textwrap.dedent(
                f"""
                calculation_type individual
                stellar_models_grid {stellar_model}
                limb_darkening_laws power2
                passbands {filter_ui}
                passbands_path {wp}
                star_effective_temperature {float(teff)}
                star_log_gravity {float(logg)}
                star_metallicity {float(mh)}
                target_names {_NPOAP_EXO_TARGET}
                output_path {op}
                user_output basic
                """
            ).strip()
            + "\n",
            encoding="utf-8",
        )
        worker = _worker_script_path()
        proc = subprocess.run(
            [sys.executable, str(worker), str(cfg_path)],
            capture_output=True,
            text=True,
            timeout=timeout_sec,
        )
        if proc.returncode != 0:
            msg = (proc.stderr or proc.stdout or "")[:4000]
            raise RuntimeError(f"ExoTETHyS (power2) a échoué (code {proc.returncode}) :\n{msg}")
        pkl = out_dir / f"{_NPOAP_EXO_TARGET}_ldc.pickle"
        if not pkl.is_file():
            msg = (proc.stderr or proc.stdout or "")[:4000]
            raise FileNotFoundError(
                f"Sortie ExoTETHyS (power2) absente : {pkl}. Sortie du worker :\n{msg}"
            )
        with open(pkl, "rb") as f:
            data = pickle.load(f)
        coeffs = data["passbands"][filter_ui]["power2"]["coefficients"]
        return np.asarray(coeffs, dtype=float).ravel()
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
