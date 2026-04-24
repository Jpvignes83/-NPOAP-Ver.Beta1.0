"""
Fusion d’épémérides transit (P, T₀ BTJD, durée) pour le flux Catalogues → TESS → LcTools.

Sources : exoplanet.eu (TAP Paris / exoplanet.epn_core) et NASA Exoplanet Archive (astroquery).
Priorité par champ : CLI > exoplanet.eu > NASA (sauf durée : CLI > NASA ; l’EU n’a pas T₁₄ standard).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from io import BytesIO
from typing import Any, Dict, List, Optional
from urllib.error import URLError, HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np

logger = logging.getLogger(__name__)

TAP_OBSPM_SYNC = "http://voparis-tap-planeto.obspm.fr/tap/sync"


@dataclass
class MergedEphemeris:
    period_days: Optional[float]
    t0_btjd: Optional[float]
    duration_days: Optional[float]
    sources: Dict[str, str] = field(default_factory=dict)


def _finite(x: Any) -> bool:
    try:
        return bool(np.isfinite(float(x)))
    except (TypeError, ValueError):
        return False


def _as_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        if hasattr(x, "mask") and np.ma.is_masked(x):
            return None
        v = float(np.asarray(x).reshape(-1)[0])
        return v if np.isfinite(v) else None
    except Exception:
        return None


def _planet_name_variants(name: str) -> List[str]:
    n = " ".join(str(name).strip().split())
    out: List[str] = []
    for v in (
        n,
        n.upper(),
        n.lower(),
        n.replace(" ", "-"),
        n.replace("-", " "),
        re.sub(r"\s+", "", n),
    ):
        if v and v not in out:
            out.append(v)
    return out


def _adql_escape(s: str) -> str:
    return s.replace("'", "''")


def _tap_query_epn_core(planet_name: str, timeout_s: float) -> Optional[dict]:
    """Interroge exoplanet.epn_core (TAP Obspm). Retourne dict period, tzero_tr (JD) ou None."""
    try:
        from astropy.io.votable import parse
    except ImportError:
        logger.warning("astropy requis pour le TAP exoplanet.eu")
        return None

    safe = _adql_escape(planet_name.strip())
    queries = [
        f"SELECT TOP 1 period, tzero_tr, target_name FROM exoplanet.epn_core WHERE target_name = '{safe}'",
        f"SELECT TOP 1 period, tzero_tr, target_name FROM exoplanet.epn_core "
        f"WHERE alt_target_name LIKE '%{_adql_escape(planet_name.split()[0] if planet_name.split() else safe)}%'",
    ]

    for adql in queries:
        try:
            data = urlencode(
                {
                    "REQUEST": "doQuery",
                    "LANG": "ADQL",
                    "FORMAT": "votable",
                    "QUERY": adql,
                }
            ).encode("ascii")
            req = Request(
                TAP_OBSPM_SYNC,
                data=data,
                method="POST",
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            with urlopen(req, timeout=timeout_s) as resp:
                raw = resp.read()
            vo = parse(BytesIO(raw))
            tab = vo.get_first_table().to_table()
            if len(tab) == 0:
                continue
            row = tab[0]
            period = _as_float(row["period"]) if "period" in row.colnames else None
            t0 = _as_float(row["tzero_tr"]) if "tzero_tr" in row.colnames else None
            if period and t0:
                return {"period": period, "tzero_tr_jd": t0, "target": str(row["target_name"]) if "target_name" in row.colnames else ""}
        except (URLError, HTTPError, TimeoutError, ValueError, KeyError) as e:
            logger.debug("TAP EU échec (%s): %s", adql[:80], e)
            continue
    return None


def _nasa_ephemeris_row(planet_name: str) -> Optional[Any]:
    try:
        from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    except ImportError:
        return None
    for variant in _planet_name_variants(planet_name):
        for table in ("pscomppars", "exoplanets"):
            try:
                t = NasaExoplanetArchive.query_object(variant, table=table)
                if t is not None and len(t) > 0:
                    return t[0]
            except Exception as e:
                logger.debug("NASA %s %s: %s", table, variant, e)
                continue
    return None


def _jd_to_btjd(t_jd: float, btjd_offset: float) -> float:
    return float(t_jd) - float(btjd_offset)


def merge_transit_ephemeris(
    planet_name: str,
    *,
    cli_period_days: Optional[float] = None,
    cli_t0_btjd: Optional[float] = None,
    cli_duration_days: Optional[float] = None,
    btjd_offset: float = 2457000.0,
    eu_timeout_s: float = 120.0,
    query_exoplanet_eu: bool = True,
    query_nasa: bool = True,
) -> MergedEphemeris:
    """
    Fusionne P, T₀ (échelle BTJD comme la colonne TIME TESS SPOC), T₁₄ en jours.
    """
    sources: Dict[str, str] = {}

    eu_period = eu_t0_btjd = None
    if query_exoplanet_eu and planet_name.strip():
        eu = _tap_query_epn_core(planet_name.strip(), float(eu_timeout_s))
        if eu:
            eu_period = eu["period"]
            eu_t0_btjd = _jd_to_btjd(eu["tzero_tr_jd"], btjd_offset)
            sources["eu_note"] = f"EU target={eu.get('target', '')}"

    nasa_row = None
    if query_nasa and planet_name.strip():
        nasa_row = _nasa_ephemeris_row(planet_name.strip())

    nasa_period = nasa_t0_btjd = nasa_dur_days = None
    if nasa_row is not None:
        row = nasa_row
        nasa_period = _as_float(row["pl_orbper"]) if "pl_orbper" in row.colnames else None
        if "pl_tranmid" in row.colnames:
            tm = _as_float(row["pl_tranmid"])
            if tm is not None:
                nasa_t0_btjd = _jd_to_btjd(tm, btjd_offset)
        if "pl_trandur" in row.colnames:
            hours = _as_float(row["pl_trandur"])
            if hours is not None and hours > 0:
                nasa_dur_days = hours / 24.0

    # Priorité CLI > EU > NASA pour P et T0
    period = cli_period_days if _finite(cli_period_days) else None
    src_p = "cli" if period is not None else None
    if period is None and eu_period is not None:
        period = eu_period
        src_p = "exoplanet.eu (TAP)"
    if period is None and nasa_period is not None:
        period = nasa_period
        src_p = "NASA Exoplanet Archive"
    if src_p:
        sources["period"] = src_p

    t0_btjd = cli_t0_btjd if _finite(cli_t0_btjd) else None
    src_t0 = "cli" if t0_btjd is not None else None
    if t0_btjd is None and eu_t0_btjd is not None:
        t0_btjd = eu_t0_btjd
        src_t0 = "exoplanet.eu (TAP)"
    if t0_btjd is None and nasa_t0_btjd is not None:
        t0_btjd = nasa_t0_btjd
        src_t0 = "NASA Exoplanet Archive"
    if src_t0:
        sources["t0"] = src_t0

    # Durée : CLI > NASA (EU : pas de colonne T14 fiable dans epn_core)
    duration = cli_duration_days if _finite(cli_duration_days) else None
    src_d = "cli" if duration is not None else None
    if duration is None and nasa_dur_days is not None:
        duration = nasa_dur_days
        src_d = "NASA (pl_trandur h → j)"
    if src_d:
        sources["duration"] = src_d

    return MergedEphemeris(
        period_days=period,
        t0_btjd=t0_btjd,
        duration_days=duration,
        sources=sources,
    )
