# core/alcdef_gaia_comparators.py
"""
Résolution des étoiles de comparaison via Gaia DR3 pour les métadonnées ALCDEF.

Les champs COMP* suivent l’usage courant ALCDEF ; magnitude Gaia G → COMPMAGBAND=GG ;
couleur BP−RP → COMPCI / COMPCIBAND=GBpGrp (indice Gaia, pas Johnson).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

logger = logging.getLogger(__name__)

# ALCDEF 2024+ : distinguer Gaia G du Sloan g
ALCDEF_GAIA_G_MAGBAND = "GG"
# Couleur BP−RP Gaia (libellé fréquent côté ALCDEF / catalogues Gaia)
ALCDEF_GAIA_BP_RP_CIBAND = "GBpGrp"


@dataclass(frozen=True)
class GaiaComparatorAlcdef:
    """Une étoile de comparaison alignée sur une entrée Gaia DR3."""

    source_id: int
    ra_deg: float
    dec_deg: float
    phot_g_mean_mag: float
    bp_rp: Optional[float]

    def comp_name(self) -> str:
        return f"Gaia DR3 {self.source_id}"

    def format_ra_hms(self) -> str:
        c = SkyCoord(ra=self.ra_deg * u.deg, dec=self.dec_deg * u.deg, frame="icrs")
        return c.ra.to_string(unit=u.hourangle, sep=":", pad=True, precision=3)

    def format_dec_dms(self) -> str:
        c = SkyCoord(ra=self.ra_deg * u.deg, dec=self.dec_deg * u.deg, frame="icrs")
        return c.dec.to_string(unit=u.deg, sep=":", pad=True, alwayssign=True, precision=2)


def _nearest_gaia_row_index(coord: SkyCoord, tab, max_sep_deg: float) -> Optional[int]:
    if tab is None or len(tab) == 0:
        return None
    gmag = np.asarray(tab["phot_g_mean_mag"], dtype=float)
    ok = np.isfinite(gmag)
    if not np.any(ok):
        return None
    ra = np.asarray(tab["ra"], dtype=float)
    dec = np.asarray(tab["dec"], dtype=float)
    cat = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    sep_deg = coord.separation(cat).deg
    sep_masked = np.where(ok, sep_deg, np.inf)
    idx = int(np.argmin(sep_masked))
    if not np.isfinite(sep_masked[idx]) or float(sep_masked[idx]) > max_sep_deg:
        return None
    return idx


def resolve_comparator_gaia_dr3(
    coord: SkyCoord,
    *,
    radius_arcsec: float = 3.0,
) -> Optional[GaiaComparatorAlcdef]:
    """
    Trouve la source Gaia DR3 la plus proche (parmi les entrées avec G fini) dans un cône.
    """
    try:
        from astroquery.gaia import Gaia
    except ImportError:
        logger.error("astroquery.gaia requis pour résoudre les comparateurs Gaia.")
        return None

    c = coord.icrs if hasattr(coord, "icrs") else coord
    max_sep_deg = max(float(radius_arcsec) / 3600.0, 1e-12)
    try:
        job = Gaia.cone_search_async(c, radius=radius_arcsec * u.arcsec)
        tab = job.get_results()
    except Exception as e:
        logger.warning("Gaia cone_search (comparateur) : %s", e)
        return None

    idx = _nearest_gaia_row_index(c, tab, max_sep_deg)
    if idx is None:
        return None

    row = tab[idx]
    try:
        sid = int(row["source_id"])
        ra = float(row["ra"])
        de = float(row["dec"])
        gmag = float(row["phot_g_mean_mag"])
    except (TypeError, ValueError, KeyError):
        return None

    bp_rp: Optional[float] = None
    if "bp_rp" in tab.colnames:
        v = row["bp_rp"]
        try:
            fv = float(v)
            if np.isfinite(fv):
                bp_rp = fv
        except (TypeError, ValueError):
            pass

    return GaiaComparatorAlcdef(
        source_id=sid,
        ra_deg=ra,
        dec_deg=de,
        phot_g_mean_mag=gmag,
        bp_rp=bp_rp,
    )


def build_alcdef_comparator_metadata_lines(records: Sequence[GaiaComparatorAlcdef]) -> List[str]:
    """Lignes KEY=value à placer après STARTBLOCK et avant STARTDATA."""
    lines: List[str] = ["EQUINOX=J2000.0"]
    for i, r in enumerate(records, start=1):
        lines.append(f"COMPNAME{i}={r.comp_name()}")
        lines.append(f"COMPRA{i}={r.format_ra_hms()}")
        lines.append(f"COMPDEC{i}={r.format_dec_dms()}")
        lines.append(f"COMPMAG{i}={r.phot_g_mean_mag:.3f}")
        lines.append(f"COMPMAGBAND{i}={ALCDEF_GAIA_G_MAGBAND}")
        if r.bp_rp is not None:
            lines.append(f"COMPCI{i}={r.bp_rp:.3f}")
            lines.append(f"COMPCIBAND{i}={ALCDEF_GAIA_BP_RP_CIBAND}")
    return lines


def serializable_dict_list(records: Sequence[GaiaComparatorAlcdef]) -> List[Dict[str, object]]:
    return [
        {
            "source_id": r.source_id,
            "ra_deg": r.ra_deg,
            "dec_deg": r.dec_deg,
            "phot_g_mean_mag": r.phot_g_mean_mag,
            "bp_rp": r.bp_rp,
        }
        for r in records
    ]


def records_from_serializable(rows: Sequence[Dict[str, object]]) -> List[GaiaComparatorAlcdef]:
    out: List[GaiaComparatorAlcdef] = []
    for o in rows:
        br = o.get("bp_rp")
        out.append(
            GaiaComparatorAlcdef(
                source_id=int(o["source_id"]),
                ra_deg=float(o["ra_deg"]),
                dec_deg=float(o["dec_deg"]),
                phot_g_mean_mag=float(o["phot_g_mean_mag"]),
                bp_rp=None if br is None else float(br),
            )
        )
    return out
