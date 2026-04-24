# core/sncosmo_sn_ia.py
"""Ajustement SN Ia via sncosmo (dépendance optionnelle, profil cosmologie SNe).

Mapping des noms de filtres Gaia (comme dans l’onglet transitoires : G, G_Bp, G_Rp)
vers les bandes intégrées sncosmo eDR3 : ``gaia::g``, ``gaia::gbp``, ``gaia::grp``.
"""

from __future__ import annotations

import logging
import re
import warnings
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

try:
    import sncosmo

    SNCOSMO_AVAILABLE = True
    SNCOSMO_IMPORT_ERROR: Optional[str] = None
except ImportError as e:
    sncosmo = None  # type: ignore
    SNCOSMO_AVAILABLE = False
    SNCOSMO_IMPORT_ERROR = str(e)


def gaia_filter_label_to_sncosmo_band(label: str) -> str:
    """
    Convertit un libellé filtre Gaia (ex. G, G_Bp, G_Rp) en nom de bande sncosmo eDR3.
    """
    raw = str(label).strip()
    u = raw.upper().replace(" ", "_").replace("-", "_")
    if u in ("G", "GAIA_G"):
        return "gaia::g"
    if u in ("G_BP", "GBP", "BP", "G_BP_MAG", "PHOT_BP_MEAN_MAG", "BP"):
        return "gaia::gbp"
    if u in ("G_RP", "GRP", "RP", "G_RP_MAG", "PHOT_RP_MEAN_MAG"):
        return "gaia::grp"
    raise ValueError(
        f"Filtre Gaia non reconnu : {label!r}. Utilisez G, G_Bp ou G_Rp (ou GBP, GRP)."
    )


def mag_to_flux_ab(mag: np.ndarray, mag_err: np.ndarray, zp: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Conversion magnitudes AB → flux et erreur pour sncosmo (relation mag = zp - 2.5 log10(flux)).
    """
    mag = np.asarray(mag, dtype=float)
    mag_err = np.asarray(mag_err, dtype=float)
    flux = np.power(10.0, 0.4 * (zp - mag))
    coeff = np.log(10.0) / 2.5
    fluxerr = flux * coeff * np.maximum(mag_err, 1e-12)
    return flux, fluxerr


def get_mw_ebv_sfd(ra_deg: float, dec_deg: float) -> float:
    """
    Récupère E(B-V) MW (SFD) sur la ligne de visée via IRSA dust service.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.irsa_dust import IrsaDust

    coord = SkyCoord(float(ra_deg) * u.deg, float(dec_deg) * u.deg, frame="icrs")
    table = IrsaDust.get_extinction_table(coord)
    if table is None or len(table) == 0:
        raise RuntimeError("Aucune donnée d'extinction retournée par IRSA.")

    # Plusieurs variantes selon versions/services: on cherche d'abord une colonne E(B-V).
    ebv_cols = [c for c in table.colnames if "ebv" in c.lower()]
    for col in ebv_cols:
        try:
            v = float(table[col][0])
            if np.isfinite(v):
                return v
        except Exception:
            continue

    # Fallback IRSA courant: colonnes par filtre avec A_SFD et A_over_E_B_V_SFD.
    # On reconstruit E(B-V) = A_SFD / (A_over_E_B_V_SFD), puis on prend la médiane.
    a_col = None
    r_col = None
    for c in table.colnames:
        cl = c.lower()
        if a_col is None and ("a_sfd" == cl or cl.endswith("a_sfd")):
            a_col = c
        if r_col is None and ("a_over_e_b_v_sfd" == cl or cl.endswith("a_over_e_b_v_sfd")):
            r_col = c
    if a_col is not None and r_col is not None:
        vals = []
        for a, r in zip(table[a_col], table[r_col]):
            try:
                a_f = float(a)
                r_f = float(r)
                if np.isfinite(a_f) and np.isfinite(r_f) and abs(r_f) > 1e-12:
                    vals.append(a_f / r_f)
            except Exception:
                continue
        if vals:
            return float(np.median(np.asarray(vals, dtype=float)))

    raise RuntimeError(
        "Impossible de déduire E(B-V) depuis la table IRSA. "
        f"Colonnes disponibles: {table.colnames}"
    )


def apply_mw_extinction_correction_gaia(
    mag: np.ndarray,
    band_labels: np.ndarray,
    ebv_mw: float,
    *,
    r_g: float = 2.74,
    r_bp: float = 3.10,
    r_rp: float = 1.97,
) -> np.ndarray:
    """
    Corrige les magnitudes observées de l'extinction MW:
      mag_corr = mag_obs - A_band, avec A_band = R_band * E(B-V).
    Coefficients Gaia (approx.) ajustables.
    """
    mag = np.asarray(mag, dtype=float)
    out = mag.copy()
    for i, b in enumerate(band_labels):
        u = str(b).strip().upper().replace(" ", "_").replace("-", "_")
        if u in ("G", "GAIA_G"):
            out[i] = mag[i] - r_g * ebv_mw
        elif u in ("G_BP", "GBP", "BP", "G_BP_MAG", "PHOT_BP_MEAN_MAG"):
            out[i] = mag[i] - r_bp * ebv_mw
        elif u in ("G_RP", "GRP", "RP", "G_RP_MAG", "PHOT_RP_MEAN_MAG"):
            out[i] = mag[i] - r_rp * ebv_mw
        else:
            # Bande inconnue côté Gaia: on laisse inchangé.
            out[i] = mag[i]
    return out


def build_photometric_table(
    time_mjd: np.ndarray,
    band_labels: np.ndarray,
    mag: np.ndarray,
    mag_err: np.ndarray,
    *,
    zp: float = 25.0,
    zpsys: str = "ab",
) -> "Any":
    """Construit une table Astropy utilisable par ``sncosmo.PhotometricData``."""
    from astropy.table import Table

    if len(time_mjd) != len(band_labels) or len(time_mjd) != len(mag):
        raise ValueError("Les colonnes time, band et mag doivent avoir la même longueur.")
    bands_out: List[str] = []
    for b in band_labels:
        bands_out.append(gaia_filter_label_to_sncosmo_band(b))
    flux, fluxerr = mag_to_flux_ab(mag, mag_err, zp)
    t = Table(
        {
            "time": np.asarray(time_mjd, dtype=float),
            "band": np.array(bands_out, dtype=str),
            "flux": flux,
            "fluxerr": fluxerr,
            "zp": np.full(len(time_mjd), float(zp)),
            "zpsys": np.full(len(time_mjd), str(zpsys).lower()),
        }
    )
    return t


def run_salt_fit(
    phot_table: Any,
    *,
    source: str = "salt2",
    z_fixed: Optional[float] = None,
    fit_z: bool = False,
    z_bounds: Tuple[float, float] = (0.001, 1.2),
) -> Dict[str, Any]:
    """
    Lance ``sncosmo.fit_lc`` sur une table photométrique (colonnes time, band, flux, …).

    Retourne un dictionnaire sérialisable pour l’interface (paramètres, erreurs, chi², message).
    """
    if not SNCOSMO_AVAILABLE or sncosmo is None:
        raise RuntimeError(
            "sncosmo n’est pas installé. Installez le profil optionnel : "
            "pip install -r requirements-cosmology-sne.txt"
        )

    # sncosmo>=2.12 n'expose plus PhotometricData en attribut public.
    # fit_lc accepte directement une table Astropy avec colonnes:
    # time, band, flux, fluxerr, zp, zpsys.
    data = phot_table
    model = sncosmo.Model(source=source)

    if fit_z:
        params = ["z", "t0", "x0", "x1", "c"]
        if z_fixed is not None:
            model.set(z=float(z_fixed))
        bounds_kw: Optional[Dict[str, Tuple[float, float]]] = {"z": z_bounds}
    else:
        if z_fixed is None:
            raise ValueError("Redshift fixe requis lorsque fit_z=False.")
        model.set(z=float(z_fixed))
        params = ["t0", "x0", "x1", "c"]
        bounds_kw = None

    try:
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            if bounds_kw is not None:
                result, fitted_model = sncosmo.fit_lc(data, model, params, bounds=bounds_kw)
            else:
                result, fitted_model = sncosmo.fit_lc(data, model, params)
    except Exception as e:
        logger.exception("sncosmo.fit_lc a échoué")
        return {"ok": False, "error": str(e), "source": source}

    chisq_raw = getattr(result, "chisq", None)
    if chisq_raw is None:
        chisq_raw = getattr(result, "chi2", None)
    try:
        chisq_f = float(chisq_raw) if chisq_raw is not None else float("nan")
    except (TypeError, ValueError):
        chisq_f = float("nan")
    ndof_raw = getattr(result, "ndof", None)
    try:
        ndof_i = int(ndof_raw) if ndof_raw is not None else -1
    except (TypeError, ValueError):
        ndof_i = -1

    out: Dict[str, Any] = {
        "ok": bool(getattr(result, "success", True)),
        "source": source,
        "message": getattr(result, "message", ""),
        "chisq": chisq_f,
        "ndof": ndof_i,
        "parameters": {},
        "errors": {},
        "covariance": None,
        "dropped_bands": [],
        "warnings": [],
    }
    for w in caught_warnings:
        wtxt = str(getattr(w, "message", "") or "").strip()
        if not wtxt:
            continue
        out["warnings"].append(wtxt)
        m = re.search(r"Dropping following bands from data:\s*(.+)$", wtxt)
        if m:
            for b in m.group(1).split(","):
                bn = b.strip()
                if bn and bn not in out["dropped_bands"]:
                    out["dropped_bands"].append(bn)
    names = list(getattr(result, "param_names", []) or [])
    par = np.asarray(getattr(result, "parameters", []), dtype=float)
    err = getattr(result, "errors", None)
    for i, name in enumerate(names):
        if i < len(par):
            out["parameters"][name] = float(par[i])
        if err is not None:
            # sncosmo peut retourner:
            # - un tableau/list indexé
            # - un dict indexé par nom de paramètre
            try:
                if isinstance(err, dict):
                    if name in err:
                        out["errors"][name] = float(err[name])
                else:
                    if i < len(err):
                        out["errors"][name] = float(err[i])
            except (TypeError, ValueError, KeyError, IndexError):
                pass
    cov = getattr(result, "covariance", None)
    if cov is not None:
        try:
            out["covariance"] = np.asarray(cov, dtype=float).tolist()
        except Exception:
            out["covariance"] = None
    out["fitted_model_repr"] = repr(fitted_model)
    return out


def load_lc_csv(
    path: str,
    *,
    col_time: str = "mjd",
    col_band: str = "band",
    col_mag: str = "mag",
    col_mag_err: str = "mag_err",
    delimiter: str = ",",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Charge un CSV minimal pour SN Ia (temps MJD, libellé filtre Gaia, mag, erreur)."""
    import csv

    times: List[float] = []
    bands: List[str] = []
    mags: List[float] = []
    merrs: List[float] = []

    def norm_key(s: str) -> str:
        return s.strip().lower().replace(" ", "_")

    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("CSV sans en-tête ou fichier vide.")
        fmap = {norm_key(h): h for h in reader.fieldnames if h}
        ct = fmap.get(norm_key(col_time))
        cb = fmap.get(norm_key(col_band))
        cm = fmap.get(norm_key(col_mag))
        ce = fmap.get(norm_key(col_mag_err))
        if not all([ct, cb, cm, ce]):
            raise ValueError(
                f"Colonnes requises introuvables. Attendu : {col_time}, {col_band}, {col_mag}, {col_mag_err}. "
                f"Trouvé : {reader.fieldnames}"
            )
        for row in reader:
            try:
                t_val = float(str(row[ct]).replace(",", "."))
                m_val = float(str(row[cm]).replace(",", "."))
                e_val = float(str(row[ce]).replace(",", "."))
            except (KeyError, TypeError, ValueError):
                continue
            b_val = row.get(cb, "").strip()
            if b_val == "":
                continue
            times.append(t_val)
            bands.append(b_val)
            mags.append(m_val)
            merrs.append(max(e_val, 1e-6))

    if len(times) < 5:
        raise ValueError(f"Trop peu de points valides ({len(times)}). Vérifiez le séparateur et les colonnes.")

    return (
        np.asarray(times, dtype=float),
        np.asarray(bands, dtype=object),
        np.asarray(mags, dtype=float),
        np.asarray(merrs, dtype=float),
    )
