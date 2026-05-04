# core/glade_alert_targets.py
"""
Associe une alerte (RA/Dec + incertitude gaussienne circulaire en degrés) au catalogue
GLADE+ (HDF5 ``Gladeplus.h5``) : galaxies avec d_L ≤ coupure (500 Mpc typique au build),
classement par vraisemblance sous la cloche d'erreur et qualité de distance, export NINA.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

import config
from core.nina_sequence_export import (
    build_nina_deep_sky_container_dict,
    sanitize_nina_filename_component,
)

logger = logging.getLogger(__name__)

DEFAULT_SIGMA_DEG = 0.05
DEFAULT_MAX_DISTANCE_MPC = 500.0
DEFAULT_CONE_DEG = 5.0
DEFAULT_TOP_N = 10
CHUNK_ROWS = 250_000


def _glade_h5_table_dataset(f: h5py.File) -> h5py.Dataset:
    """
    Fichier ``Gladeplus.h5`` produit par le convertisseur PyTables (``/catalog`` + ``table``).
    La lecture côté NPOAP se fait en HDF5 pur via h5py (pas de dépendance PyTables à l'exécution).
    """
    try:
        node = f["catalog"]
    except KeyError as e:
        raise FileNotFoundError("Gladeplus.h5 : groupe /catalog absent.") from e
    if isinstance(node, h5py.Dataset):
        return node
    if "table" in node:
        return node["table"]  # type: ignore[return-value]
    raise FileNotFoundError(
        "Gladeplus.h5 : format inattendu (attendu : PyTables, groupe /catalog avec dataset table)."
    )


def _decode_str(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="replace").strip()
    return str(x).strip()


def best_glade_display_name(row: Dict[str, Any]) -> str:
    if _decode_str(row.get("name_SDSS_DR16Q")) and _decode_str(row.get("name_SDSS_DR16Q")) != "null":
        return "SDSS " + _decode_str(row["name_SDSS_DR16Q"])
    if _decode_str(row.get("name_GWGC")) and _decode_str(row.get("name_GWGC")) != "null":
        return _decode_str(row["name_GWGC"])
    if _decode_str(row.get("name_HyperLEDA")) and _decode_str(row.get("name_HyperLEDA")) != "null":
        return _decode_str(row["name_HyperLEDA"])
    pg = _decode_str(row.get("no_PGC"))
    if pg and pg != "null":
        return "PGC " + pg
    if _decode_str(row.get("name_2MASS")) and _decode_str(row.get("name_2MASS")) != "null":
        return "2MASS " + pg if pg else "2MASS"
    if _decode_str(row.get("name_WISExSCOS")) and _decode_str(row.get("name_WISExSCOS")) != "null":
        return "WISE " + _decode_str(row["name_WISExSCOS"])
    ng = row.get("no_GLADE", 0)
    try:
        return "GLADE+ %s" % int(ng)
    except (TypeError, ValueError):
        return "GLADE+"


def _row_dict(block: np.ndarray, i: int) -> Dict[str, Any]:
    return {name: block[name][i] for name in block.dtype.names}


def ranked_glade_targets_for_alert(
    ra_deg: float,
    dec_deg: float,
    *,
    sigma_deg: float = DEFAULT_SIGMA_DEG,
    h5_path: Optional[Path] = None,
    max_distance_mpc: float = DEFAULT_MAX_DISTANCE_MPC,
    cone_deg: float = DEFAULT_CONE_DEG,
    top_n: int = DEFAULT_TOP_N,
) -> List[Dict[str, Any]]:
    path = Path(h5_path) if h5_path is not None else Path(getattr(config, "GLADE_CATALOG_H5_PATH"))
    if not path.is_file():
        raise FileNotFoundError(
            "Catalogue GLADE+ introuvable : %s\n"
            "Copiez Gladeplus.h5 dans ce dossier (fourni avec la distribution NPOAP : "
            "%%USERPROFILE%%\\.npoap\\catalogues\\)." % path
        )

    sigma_deg = max(float(sigma_deg), 1e-6)
    cone_deg = max(float(cone_deg), sigma_deg * 3)

    alert = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)

    scores: List[float] = []
    metas: List[Dict[str, Any]] = []

    with h5py.File(str(path), mode="r") as h5:
        dset = _glade_h5_table_dataset(h5)
        nrows = int(dset.shape[0])
        start = 0
        while start < nrows:
            stop = min(start + CHUNK_ROWS, nrows)
            block = dset[start:stop]

            ra_g = np.asarray(block["RA"], dtype=np.float64)
            dec_g = np.asarray(block["Dec"], dtype=np.float64)
            d_l = np.asarray(block["d_L"], dtype=np.float64)
            d_l_err = np.asarray(block["d_L_err"], dtype=np.float64)

            good = np.isfinite(ra_g) & np.isfinite(dec_g) & np.isfinite(d_l)
            good &= d_l <= max_distance_mpc

            if not np.any(good):
                start = stop
                continue

            ix = np.flatnonzero(good)
            coords = SkyCoord(ra_g[good] * u.deg, dec_g[good] * u.deg)
            sep = alert.separation(coords).deg
            in_cone = sep <= cone_deg
            if not np.any(in_cone):
                start = stop
                continue

            sel = ix[in_cone]
            sep_c = sep[in_cone]

            for k in range(len(sel)):
                bi = sel[k]
                theta = float(sep_c[k])
                err = float(d_l_err[bi]) if np.isfinite(d_l_err[bi]) else 1e3
                lik = float(np.exp(-0.5 * (theta / sigma_deg) ** 2))
                qual = 1.0 / (1.0 + max(err, 0.0))
                score = lik * qual

                row = _row_dict(block, int(bi))
                scores.append(score)
                metas.append(
                    {
                        "ra_deg": float(ra_g[bi]),
                        "dec_deg": float(dec_g[bi]),
                        "d_L": float(d_l[bi]),
                        "d_L_err": err,
                        "no_GLADE": int(row["no_GLADE"]),
                        "name": best_glade_display_name(row),
                        "score": score,
                        "sep_deg": theta,
                    }
                )

            start = stop

    if not scores:
        return []

    order = np.argsort(-np.asarray(scores))
    out: List[Dict[str, Any]] = []
    for k in order[:top_n]:
        out.append(metas[int(k)])
    return out


def export_nina_for_glade_targets(
    targets: List[Dict[str, Any]],
    dest_dir: Path,
    *,
    filename_prefix: str = "",
) -> List[Path]:
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    prefix = sanitize_nina_filename_component(filename_prefix) + "_" if filename_prefix else ""
    written: List[Path] = []
    for rank, t in enumerate(targets, start=1):
        name = sanitize_nina_filename_component(t["name"], max_len=60)
        fname = "%s%02d_%s.json" % (prefix, rank, name)
        path = dest_dir / fname
        payload = build_nina_deep_sky_container_dict(t["name"], t["ra_deg"], t["dec_deg"])
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        written.append(path)
    return written


def transient_details_to_glade_export(
    details: Dict[str, Any],
    dest_dir: Path,
    *,
    sigma_deg: Optional[float] = None,
    top_n: int = DEFAULT_TOP_N,
) -> Tuple[List[Dict[str, Any]], List[Path]]:
    ra = float(details["ra"])
    dec = float(details["dec"])
    if sigma_deg is not None:
        sig = float(sigma_deg)
    elif details.get("localization_sigma_deg") is not None:
        sig = float(details["localization_sigma_deg"])
    else:
        sig = DEFAULT_SIGMA_DEG
    evt = (
        details.get("source_name")
        or details.get("trigger_id")
        or details.get("astro_colibri_id")
        or "alert"
    )
    evt_safe = sanitize_nina_filename_component(str(evt), max_len=40)
    targets = ranked_glade_targets_for_alert(ra, dec, sigma_deg=sig, top_n=top_n)
    paths = export_nina_for_glade_targets(targets, dest_dir, filename_prefix=evt_safe) if targets else []
    return targets, paths
