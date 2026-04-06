# core/gaia_pylightcurve_support.py
"""
Support des filtres Gaia (G, G_BP, G_RP) pour pylightcurve.

pylightcurve ne fournit ni fichiers ``*.pass`` ni grilles exotethys pour ces bandes.
Ce module :

1. **Passbands** — installe des courbes de transmission (Angstrom, réponse normalisée)
   dans le répertoire photometry de pylightcurve (``~/.pylightcurve4/photometry/``).
   Par défaut : copies depuis ``data/gaia_pylightcurve/*.pass`` (profils analytiques
   calés sur les plages effectives type EDR3). Vous pouvez remplacer ces fichiers
   par les tables ASCII officielles ESA (page « Gaia EDR3 passbands » sur cosmos.esa.int)
   puis relancer l’application.

2. **Coefficients de limb darkening** — coefficients **power-2** tabulés (Claret & Southworth
   2022, VizieR ``J/A+A/664/A128``), identiques à ceux utilisés pour la loi power-2 dans NPOAP.
   pylightcurve n’implémentant pas la loi power-2 dans ``transit_integrated``, on projette
   cette loi sur la **base quadratique** utilisée par pylightcurve
   (``I = 1 - a1(1-μ) - a2(1-μ)²`` avec ``μ = √(1-r²)``), par moindres carrés sur μ∈(0,1].
   La méthode d’intégration pylightcurve est alors forcée à ``quad`` pour ces filtres.

Références : Claret & Southworth 2022, A&A, 664, A128 ; Gaia EDR3 passbands (ESA).
"""

from __future__ import annotations

import logging
import os
import shutil
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

GAIA_PYLIGHTCURVE_FILTERS = frozenset({"G", "G_BP", "G_RP"})

# SVO (Spanish Virtual Observatory) — profils Gaia EDR3 / DR3 ; remplace le bundle si le téléchargement réussit.
_GAIA_PASSBAND_SOURCES: dict[str, tuple[str, ...]] = {
    "G": (
        "https://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=GAIA/GAIA3.G",
    ),
    "G_BP": (
        "https://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=GAIA/GAIA3.Gbp",
    ),
    "G_RP": (
        "https://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=GAIA/GAIA3.Grp",
    ),
}

_BUNDLE_DIR = Path(__file__).resolve().parent.parent / "data" / "gaia_pylightcurve"


def is_gaia_pylightcurve_filter(filter_name: str) -> bool:
    return bool(filter_name) and filter_name.strip() in GAIA_PYLIGHTCURVE_FILTERS


def _ui_to_claret_band(ui: str) -> str:
    u = ui.strip()
    if u == "G":
        return "G"
    if u == "G_BP":
        return "BP"
    if u == "G_RP":
        return "RP"
    raise ValueError(f"Filtre Gaia UI inconnu : {ui!r}")


def power2_to_pylightcurve_quadratic(c: float, alpha: float) -> Tuple[float, float]:
    """
    Projette la loi power-2 de Claret, I(μ)/I(1) = 1 - c(1 - μ^α), sur la forme quadratique
    pylightcurve : I(μ) = 1 - a1(1-μ) - a2(1-μ)² (avec μ le cosinus de l’angle au bord).
    """
    c = float(c)
    alpha = float(alpha)
    mu = np.linspace(0.04, 1.0, 200)
    x = 1.0 - mu
    y = c * (1.0 - np.power(mu, alpha))
    a = np.column_stack([x, x * x])
    u1, u2 = np.linalg.lstsq(a, y, rcond=None)[0]
    u1 = float(u1)
    u2 = float(u2)
    # Bornes douces pour rester dans un régime physique courant
    u1 = float(np.clip(u1, 0.0, 1.0))
    u2 = float(np.clip(u2, 0.0, 1.0))
    if u1 + u2 > 1.05:
        s = u1 + u2
        u1, u2 = u1 / s * 1.0, u2 / s * 1.0
    return u1, u2


def _parse_two_column_spectrum(text: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Interprète un fichier type synphot / ASCII : colonnes longueur d’onde (nm ou Å), transmission."""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            w = float(parts[0])
            t = float(parts[1])
        except ValueError:
            continue
        rows.append((w, t))
    if len(rows) < 5:
        return None
    w = np.array([r[0] for r in rows], dtype=float)
    t = np.array([r[1] for r in rows], dtype=float)
    t = np.clip(t, 0.0, None)
    if np.nanmax(t) <= 0:
        return None
    t = t / np.nanmax(t)
    wmax = float(np.nanmax(w))
    # pylightcurve attend des longueurs d’onde en Å (cf. JOHNSON_V.pass ~ 4700–7000).
    if wmax < 2.0:
        w = w * 1.0e4
    elif wmax < 4000.0:
        w = w * 10.0
    return w, t


def _try_download_passband(ui_name: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    for url in _GAIA_PASSBAND_SOURCES.get(ui_name, ()):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "NPOAP/1.0"})
            with urllib.request.urlopen(req, timeout=25) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
            parsed = _parse_two_column_spectrum(raw)
            if parsed is not None:
                logger.info("Passband Gaia %s téléchargée depuis %s", ui_name, url)
                return parsed
        except (urllib.error.URLError, OSError, ValueError, UnicodeDecodeError) as e:
            logger.debug("Téléchargement passband %s (%s) : %s", ui_name, url, e)
    return None


def _build_parametric_passband(ui_name: str) -> Tuple[np.ndarray, np.ndarray]:
    """Profil analytique de secours (λ en Å), normalisé à 1."""
    lam = np.linspace(2800.0, 11800.0, 120)

    def L(z: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-z))

    if ui_name == "G":
        center, w, lo, hi = 5830.0, 2100.0, 3300.0, 10500.0
    elif ui_name == "G_BP":
        center, w, lo, hi = 5030.0, 450.0, 3300.0, 7000.0
    elif ui_name == "G_RP":
        center, w, lo, hi = 7970.0, 380.0, 6100.0, 10800.0
    else:
        raise ValueError(ui_name)
    t = np.exp(-0.5 * ((lam - center) / w) ** 2) * L((lam - lo) / 80.0) * L((hi - lam) / 80.0)
    t = t / np.max(t)
    return lam, t


def _write_pass_file(path: str, wave_angstrom: np.ndarray, trans: np.ndarray) -> None:
    arr = np.column_stack([wave_angstrom, trans])
    np.savetxt(path, arr, fmt="%.8e")


def ensure_gaia_passbands_installed(force: bool = False) -> None:
    """
    Garantit la présence de G.pass, G_BP.pass, G_RP.pass dans le répertoire photometry pylightcurve.
    Invalide le cache interne ``all_filters`` de pylightcurve pour prise en compte immédiate.
    """
    from pylightcurve.__databases__ import plc_data

    dest_dir = plc_data.photometry()
    os.makedirs(dest_dir, exist_ok=True)

    for name in ("G", "G_BP", "G_RP"):
        dest = os.path.join(dest_dir, f"{name}.pass")
        if os.path.isfile(dest) and not force:
            continue
        data = _try_download_passband(name)
        if data is None:
            bundle = _BUNDLE_DIR / f"{name}.pass"
            if bundle.is_file():
                shutil.copy2(bundle, dest)
                logger.info("Passband Gaia %s installée depuis le bundle NPOAP.", name)
            else:
                w, t = _build_parametric_passband(name)
                _write_pass_file(dest, w, t)
                logger.warning(
                    "Passband Gaia %s générée analytiquement (bundle %s absent).",
                    name,
                    bundle,
                )
        else:
            w, t = data
            _write_pass_file(dest, w, t)

    plc_data.all_filters_data = None


def refresh_gaia_passbands_from_network() -> None:
    """Réécrit G.pass / G_BP.pass / G_RP.pass (téléchargement SVO si disponible, sinon bundle / analytique)."""
    ensure_gaia_passbands_installed(force=True)


def inject_gaia_filter_from_claret2022(
    planet,
    filter_ui: str,
    vturb_kms: float = 2.0,
) -> None:
    """
    Enregistre le filtre Gaia sur l’objet Planet : LDC quadratiques (projection power-2 Claret 2022)
    et Fp/Fs via la réponse passband installée.
    """
    from pylightcurve.models.exoplanet_lc import fp_over_fs

    from core.claret_ld_vizier import nearest_power2_gaia

    if not is_gaia_pylightcurve_filter(filter_ui):
        raise ValueError(filter_ui)

    ensure_gaia_passbands_installed()

    teff = float(planet.stellar_temperature)
    logg = float(planet.stellar_logg)
    met = float(planet.stellar_metallicity)
    band = _ui_to_claret_band(filter_ui)

    info = nearest_power2_gaia(teff, logg, met, vturb_kms, band=band)
    u1, u2 = power2_to_pylightcurve_quadratic(info["c"], info["alpha"])

    rp = float(planet.rp_over_rs)
    fp = fp_over_fs(
        rp,
        planet.sma_over_rs,
        planet.albedo,
        planet.emissivity,
        planet.stellar_temperature,
        filter_ui,
    )

    planet.filters.pop(filter_ui, None)
    planet.add_filter(filter_ui, rp, u1, u2, 0.0, 0.0, fp)
    logger.debug(
        "Filtre Gaia %s : c=%.5f α=%.5f → quad pylightcurve u1=%.5f u2=%.5f (Claret+2022)",
        filter_ui,
        info["c"],
        info["alpha"],
        u1,
        u2,
    )


def prepare_gaia_pylightcurve_transit(planet, filter_ui: str, vturb_kms: float = 2.0) -> None:
    """
    Avant ``transit_integrated`` : restaure le ``ldc_method`` d’origine hors Gaia ;
    pour G / G_BP / G_RP, installe passbands + filtre (LDC Claret 2022 → quad) et force ``ldc_method='quad'``.
    """
    if planet is None or not filter_ui:
        return

    if not hasattr(planet, "_npoap_ldc_original"):
        planet._npoap_ldc_original = planet.ldc_method

    if not is_gaia_pylightcurve_filter(filter_ui):
        planet.ldc_method = planet._npoap_ldc_original
        return

    inject_gaia_filter_from_claret2022(planet, filter_ui, vturb_kms=vturb_kms)
    planet.ldc_method = "quad"
