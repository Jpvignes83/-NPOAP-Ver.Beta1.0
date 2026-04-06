# core/claret_ld_vizier.py
"""
Coefficients de limb-darkening, loi power-2, depuis Claret & Southworth (2022),
catalogue VizieR **J/A+A/664/A128** (table1 : Gaia, Kepler, TESS ; modèles ATLAS).

Forme adoptée dans l’article (et dans NPOAP pour power-2) :
    I(μ) / I(1) = 1 − g · (1 − μ^h)

Les colonnes VizieR **gG**, **hG** (etc.) correspondent donc à **c** et **α** dans
`limb_darkening_power2` (c ≡ g, α ≡ h).

Référence : Claret A., Southworth J., 2022, A&A, 664, A128
https://doi.org/10.1051/0004-6361/202243827
VizieR : https://vizier.cds.unistra.fr/viz-bin/cat/J/A+A/664/A128
"""

from __future__ import annotations

import logging
from typing import Any, Dict

import numpy as np

logger = logging.getLogger(__name__)

_TABLE1_ID = "J/A+A/664/A128/table1"
_TABLE1_CACHE = None

# Passbands Gaia DR3 (réponse intégrée dans Claret+ 2022, table1)
GAIA_BAND_COLUMNS: Dict[str, tuple] = {
    "G": ("gG", "hG", "chi2G"),
    "BP": ("gBP", "hBP", "chi2BP"),
    "RP": ("gRP", "hRP", "chi2RP"),
}


def clear_claret_2022_cache() -> None:
    """Vide le cache mémoire de la table (tests / rechargement)."""
    global _TABLE1_CACHE
    _TABLE1_CACHE = None


def load_claret_2022_table1():
    """
    Télécharge (une fois) la table1 du catalogue J/A+A/664/A128 via astroquery.
    """
    global _TABLE1_CACHE
    if _TABLE1_CACHE is not None:
        return _TABLE1_CACHE
    try:
        from astroquery.vizier import Vizier
    except ImportError as e:
        raise ImportError(
            "Le paquet astroquery est requis pour interroger VizieR "
            "(pip install astroquery)."
        ) from e

    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("J/A+A/664/A128")
    if _TABLE1_ID not in result.keys():
        raise KeyError(f"Table attendue {_TABLE1_ID}, obtenu : {list(result.keys())}")
    _TABLE1_CACHE = result[_TABLE1_ID]
    logger.info("Claret & Southworth (2022) table1 chargée : %s lignes", len(_TABLE1_CACHE))
    return _TABLE1_CACHE


def nearest_power2_gaia(
    teff: float,
    logg: float,
    feh_dex: float,
    vturb_kms: float = 2.0,
    band: str = "G",
) -> Dict[str, Any]:
    """
    Retourne les coefficients power-2 pour un passage Gaia à partir de la grille ATLAS
    (ligne la plus proche en Teff, log g, [M/H], v_turb).

    Parameters
    ----------
    teff : float
        Température effective (K). Grille tabulée ~3500–50 000 K.
    logg : float
        log g (cgs).
    feh_dex : float
        [M/H] ou [Fe/H] en dex (colonne **Z** du catalogue).
    vturb_kms : float
        Vitesse de microturbulence (km/s). Grille : 0, 1, 2, 4, 8.
    band : str
        ``'G'``, ``'BP'`` (GBP), ``'RP'`` (GRP).

    Returns
    -------
    dict
        ``c``, ``alpha`` (pour NPOAP), métadonnées de la ligne utilisée, source biblio.
    """
    b = band.strip().upper()
    if b in ("GBP", "G_BP"):
        b = "BP"
    if b in ("GRP", "G_RP"):
        b = "RP"
    if b not in GAIA_BAND_COLUMNS:
        raise ValueError(f"band doit être parmi {list(GAIA_BAND_COLUMNS)} (reçu : {band!r})")

    tbl = load_claret_2022_table1()
    gcol, hcol, xcol = GAIA_BAND_COLUMNS[b]

    Teff = np.asarray(tbl["Teff"], dtype=float)
    Logg = np.asarray(tbl["logg"], dtype=float)
    Ztab = np.asarray(tbl["Z"], dtype=float)
    Vel = np.asarray(tbl["Vel"], dtype=float)

    dist = (
        ((Teff - float(teff)) / 500.0) ** 2
        + ((Logg - float(logg)) / 0.35) ** 2
        + ((Ztab - float(feh_dex)) / 0.2) ** 2
        + ((Vel - float(vturb_kms)) / 2.0) ** 2
    )
    i = int(np.argmin(dist))
    row = tbl[i]

    return {
        "c": float(row[gcol]),
        "alpha": float(row[hcol]),
        "chi2_fit": float(row[xcol]) if xcol in row.colnames else None,
        "band_gaia": b,
        "row_teff": float(row["Teff"]),
        "row_logg": float(row["logg"]),
        "row_z": float(row["Z"]),
        "row_vel": float(row["Vel"]),
        "bibcode": "2022A&A...664A.128C",
        "vizier_table": _TABLE1_ID,
        "note": "g,h Claret ≡ c,α NPOAP ; I = 1 − c(1 − μ^α)",
    }


def format_lookup_summary(info: Dict[str, Any]) -> str:
    """Texte court pour boîte de dialogue."""
    return (
        f"Passband Gaia {info['band_gaia']}\n"
        f"c = g = {info['c']:.6f}\n"
        f"α = h = {info['alpha']:.6f}\n"
        f"Ligne grille : Teff={info['row_teff']:.0f} K, log g={info['row_logg']:.2f}, "
        f"Z={info['row_z']:.2f}, vt={info['row_vel']:.1f} km/s\n"
        f"Source : Claret & Southworth (2022), {_TABLE1_ID}"
    )
