"""
Coefficients de transformation photométrique (réduction NPOAP).

À partir d'une image FITS astrométrisée : photométrie d'ouverture sur étoiles détectées,
magnitudes de référence **Gaia DR3** (G, BP, RP) ou **Johnson B, V** issues du catalogue
**APASS DR9** (Vizier II/336, photométrie Johnson / Cousins côté AAVSO),
correction d'extinction grossière via E(B-V), ajustement linéaire M_ref = a * m_inst + b.

**E(B-V) automatique** (côté interface réduction) : au lancement du calcul des magnitudes,
requête **IRSA Galactic Dust Reddening & Extinction** (NASA/IPAC) — cartes **Schlegel et al. (1998)**
et recalibration **Schlafly & Finkbeiner (2011)** au **centre géométrique** du FITS (WCS) ;
en cas d’échec IRSA : le champ E(B-V) n’est pris en compte que s’il contient un nombre ≥ 0 ;
sinon 0 ; une ligne est ajoutée au journal et un commentaire peut apparaître dans le tableau.

**Johnson Rc Cousins** : APASS ne publie pas Rc directement ; on utilise les magnitudes **g′** et **r′**
(proches du système SDSS) et la transformation **ugriz → Rc** de **Lupton (2005)** (SDSS DR5),
documentée avec les magnitudes SDSS ; l’écart APASS/SDSS introduit une incertitude systémique
supplémentaire (typ. ~0,01–0,03 mag selon les champs).
"""

from __future__ import annotations

import csv
import logging
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
import astropy.units as u

logger = logging.getLogger(__name__)

# Seuil pour le calcul des coefficients : rejette les points trop « brillants » (saturés / aberrants).
MAGNITUDE_MIN_FOR_COEFF = -15.0

# SNR minimal sur le flux net d'ouverture (ciel soustrait) pour entrer dans l'ajustement transformation.
MIN_SNR_FOR_TRANS_COEFF = 20.0

# Ligne de synthèse (médiane après clipping 2σ) écrite dans le CSV paire ; exclue des recalculs suivants.
AGGREGATE_ROW_FITS_MARKER = "__AGG_MEDIAN_2SIG__"


def safe_coeff_csv_token(s: str) -> str:
    """Fragment de nom de fichier pour ``{filtre}__{bande_ref}.csv`` (comme ``append_coefficient_record``)."""
    t = re.sub(r"[^\w.\-]+", "_", str(s), flags=re.UNICODE).strip("_")
    return t[:80] if t else "unknown"


def transformation_storage_dir() -> Path:
    d = Path.home() / ".npoap" / "Coeftransformation"
    d.mkdir(parents=True, exist_ok=True)
    return d


# Bandes : (libellé, identifiant interne)
REFERENCE_BAND_CHOICES: list[tuple[str, str]] = [
    ("Gaia G", "gaia_g"),
    ("Gaia BP", "gaia_bp"),
    ("Gaia RP", "gaia_rp"),
    ("Johnson B (APASS DR9)", "johnson_b"),
    ("Johnson V (APASS DR9)", "johnson_v"),
    (
        "Johnson Rc (APASS g′r′ → Lupton 2005)",
        "johnson_rc",
    ),
]


def field_center_icrs_from_fits(fits_path: Path) -> SkyCoord:
    """
    Centre du champ (moyenne des quatre coins pixel → ciel) en ICRS, pour requêtes catalogue.
    """
    fits_path = Path(fits_path)
    with fits.open(fits_path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    if data is None:
        raise ValueError("HDU primaire sans données.")
    arr = np.asarray(data)
    if arr.ndim != 2:
        raise ValueError("Image 2D requise.")
    wcs = WCS(header)
    if hasattr(wcs, "has_celestial") and not wcs.has_celestial:
        raise ValueError("WCS céleste absent.")
    ny, nx = arr.shape
    corners = np.array([[0, 0], [nx, 0], [0, ny], [nx, ny]], dtype=float)
    cc = pixel_to_skycoord(corners[:, 0], corners[:, 1], wcs, origin=0)
    return SkyCoord(ra=cc.ra.mean(), dec=cc.dec.mean(), frame="icrs")


def query_ebv_irsa_dust_at_position(center: SkyCoord) -> tuple[float, str]:
    """
    E(B-V) à une position céleste via **IRSA Dust Extinction Service** (astroquery).

    Priorité aux valeurs **Schlafly & Finkbeiner (2011)** (recalibration des cartes SFD),
    puis secours **Schlegel et al. (1998) SFD** brut.

    Returns
    -------
    ebv : float
        E(B-V) (mag).
    provenance : str
        Courte description pour journal / interface.
    """
    try:
        from astroquery.ipac.irsa.irsa_dust import IrsaDust
    except ImportError as e:
        raise RuntimeError(
            "astroquery avec le module ``ipac.irsa.irsa_dust`` est requis "
            "(requête HTTP vers https://irsa.ipac.caltech.edu/applications/DUST/)."
        ) from e

    c = center.icrs
    tab = IrsaDust.get_query_table(c, section="ebv")
    if tab is None or len(tab) == 0:
        raise RuntimeError("IRSA Dust : réponse vide (section E(B-V)).")

    row = tab[0]

    def _try_float(col: str) -> float | None:
        if col not in tab.colnames:
            return None
        try:
            v = float(row[col])
        except (TypeError, ValueError):
            return None
        return v if np.isfinite(v) else None

    for col, prov in (
        ("ext SandF ref", "Schlafly & Finkbeiner (2011), ref. IRSA Dust"),
        ("ext SandF mean", "Schlafly & Finkbeiner (2011), moyenne IRSA Dust"),
        ("ext SFD ref", "Schlegel et al. (1998) SFD, IRSA Dust"),
    ):
        v = _try_float(col)
        if v is not None and v >= 0.0:
            return v, f"{prov} ({col})"

    raise RuntimeError(
        "IRSA Dust : impossible de lire E(B-V) (colonnes ext SandF / ext SFD absentes ou invalides)."
    )


def query_ebv_irsa_for_fits_center(fits_path: Path) -> tuple[float, str]:
    """E(B-V) IRSA au centre géométrique du champ défini par ``fits_path``."""
    center = field_center_icrs_from_fits(fits_path)
    ebv, prov = query_ebv_irsa_dust_at_position(center)
    return ebv, f"{prov} — centre ICRS {center.ra.deg:.5f}°, {center.dec.deg:.5f}°"


def read_obs_filter_from_header(header: fits.Header) -> str:
    """Lit un nom de filtre courant dans les en-têtes FITS (première valeur non vide)."""
    keys = ("FILTER", "FILTNAM", "FILTER1", "INSFLNAM", "FILTER2", "FW_NAME")
    for k in keys:
        v = header.get(k)
        if v is None or v == "":
            continue
        s = str(v).strip()
        if s and s.upper() not in ("CLEAR", "NONE", "OPEN"):
            return s
    return ""


def _ebv_to_extinction_ebv(band_id: str) -> float:
    """R_λ = A_λ / E(B-V) (ordres de grandeur Fitzpatrick+1999 / Schlafly pour Gaia / Johnson)."""
    r = {
        "gaia_g": 0.789,
        "gaia_bp": 1.072,
        "gaia_rp": 0.662,
        "johnson_b": 4.1,
        "johnson_v": 3.1,
        "johnson_rc": 2.31,
    }
    return r.get(band_id, 3.1)


def _as_float_col(tab: Table, *candidates: str) -> np.ndarray:
    for name in candidates:
        if name in tab.colnames:
            return np.asarray(tab[name], dtype=float)
    raise KeyError(f"Aucune colonne parmi {candidates} dans {list(tab.colnames)[:20]}…")


def _as_float_col_optional(tab: Table, n: int, *candidates: str) -> np.ndarray:
    """Comme _as_float_col mais retourne NaN si aucune colonne trouvée."""
    for name in candidates:
        if name in tab.colnames:
            return np.asarray(tab[name], dtype=float)
    return np.full(n, np.nan, dtype=float)


def johnson_rc_from_apass_gr_lupton2005(g: np.ndarray, r: np.ndarray) -> np.ndarray:
    """
    R magnitude Johnson **Cousins** à partir de magnitudes type SDSS **g**, **r**.

    Formule **Lupton (2005)** (équations SDSS DR5, étoiles) :
    ``R = r - 0.1837 * (g - r) - 0.0971`` (RMS ~0,011 mag sur l’échantillon Stetson).

    Ici **g** et **r** sont les **g′**, **r′** APASS (proches SDSS, pas identiques).
    """
    g = np.asarray(g, dtype=float)
    r = np.asarray(r, dtype=float)
    gr = g - r
    out = np.full_like(g, np.nan, dtype=float)
    ok = np.isfinite(g) & np.isfinite(r)
    out[ok] = r[ok] - 0.1837 * gr[ok] - 0.0971
    return out


def reference_magnitude_gaia(
    G: np.ndarray, BP: np.ndarray, RP: np.ndarray, band_id: str
) -> np.ndarray:
    """Magnitudes Gaia natives uniquement (G, BP, RP)."""
    G = np.asarray(G, dtype=float)
    BP = np.asarray(BP, dtype=float)
    RP = np.asarray(RP, dtype=float)
    out = np.full_like(G, np.nan, dtype=float)
    if band_id == "gaia_g":
        ok = np.isfinite(G)
        out[ok] = G[ok]
    elif band_id == "gaia_bp":
        ok = np.isfinite(BP)
        out[ok] = BP[ok]
    elif band_id == "gaia_rp":
        ok = np.isfinite(RP)
        out[ok] = RP[ok]
    return out


def query_gaia_dr3_field(center: SkyCoord, radius_deg: float, mag_limit: float = 18.5) -> Table:
    """Requête Vizier Gaia DR3 (I/355/gaiadr3) sur un champ."""
    try:
        from astroquery.vizier import Vizier
    except ImportError as e:
        raise RuntimeError("astroquery est requis pour le catalogue Gaia.") from e

    v = Vizier(
        columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "BPmag", "RPmag"],
        row_limit=8000,
    )
    res = v.query_region(
        center,
        width=2 * radius_deg * u.deg,
        height=2 * radius_deg * u.deg,
        catalog="I/355/gaiadr3",
    )
    if not res or len(res[0]) == 0:
        return Table(
            names=("RA_ICRS", "DE_ICRS", "Gmag", "BPmag", "RPmag"),
            dtype=("f8", "f8", "f8", "f8", "f8"),
        )

    tab = res[0]
    ra = _as_float_col(tab, "RA_ICRS", "_RA.icrs", "RAICRS")
    dec = _as_float_col(tab, "DE_ICRS", "_DE.icrs", "DEICRS")
    gm = _as_float_col(tab, "Gmag", "GMAG", "phot_g_mean_mag")
    try:
        bp = _as_float_col(tab, "BPmag", "BPMAG", "phot_bp_mean_mag")
    except KeyError:
        bp = np.full_like(gm, np.nan)
    try:
        rp = _as_float_col(tab, "RPmag", "RPMAG", "phot_rp_mean_mag")
    except KeyError:
        rp = np.full_like(gm, np.nan)

    ok = np.isfinite(gm) & (gm <= mag_limit)
    out = Table()
    out["RA_ICRS"] = ra[ok]
    out["DE_ICRS"] = dec[ok]
    out["Gmag"] = gm[ok]
    out["BPmag"] = bp[ok]
    out["RPmag"] = rp[ok]
    return out


def query_apass9_field(
    center: SkyCoord,
    radius_deg: float,
    *,
    mag_limit_v: float = 16.5,
) -> Table:
    """
    Requête Vizier APASS DR9 (II/336/apass9).

    Colonnes : positions, **B**, **V** (Johnson / Cousins AAVSO), **g′**, **r′** (système proche SDSS)
    lorsque le serveur Vizier les renvoie (nécessaire pour **Johnson Rc** via Lupton 2005).

    Référence : AAVSO Photometric All Sky Survey, DR9 (Henden+ 2016, Vizier II/336).
    """
    try:
        from astroquery.vizier import Vizier
    except ImportError as e:
        raise RuntimeError("astroquery est requis pour le catalogue APASS.") from e

    def _do_query(cols: list[str]):
        v = Vizier(columns=cols, row_limit=8000)
        return v.query_region(
            center,
            width=2 * radius_deg * u.deg,
            height=2 * radius_deg * u.deg,
            catalog="II/336/apass9",
        )

    res = None
    for cols in (
        ["RAJ2000", "DEJ2000", "Bmag", "Vmag", "g'mag", "r'mag"],
        ["RAJ2000", "DEJ2000", "Bmag", "Vmag"],
    ):
        try:
            res = _do_query(cols)
            if res and len(res[0]) > 0:
                break
        except Exception as e:
            logger.debug("Requête APASS colonnes %s : %s", cols, e)
            res = None

    if not res or len(res[0]) == 0:
        return Table(
            names=("RAJ2000", "DEJ2000", "Bmag", "Vmag", "gpmag", "rpmag"),
            dtype=("f8", "f8", "f8", "f8", "f8", "f8"),
        )

    tab = res[0]
    n = len(tab)
    try:
        ra = _as_float_col(tab, "RAJ2000", "_RAJ2000", "RA_ICRS")
        dec = _as_float_col(tab, "DEJ2000", "_DEJ2000", "DE_ICRS")
    except KeyError as e:
        raise RuntimeError(
            "Colonnes de position APASS inattendues. Vérifiez la table II/336/apass9 sur Vizier."
        ) from e
    bmag = _as_float_col_optional(tab, n, "Bmag", "BMAG")
    vmag = _as_float_col_optional(tab, n, "Vmag", "VMAG")
    gpmag = _as_float_col_optional(tab, n, "g'mag", "gpmag", "gmag", "_g'mag")
    rpmag = _as_float_col_optional(tab, n, "r'mag", "rpmag", "rmag", "_r'mag")
    if (not np.any(np.isfinite(gpmag))) or (not np.any(np.isfinite(rpmag))):
        for c in tab.colnames:
            cl = c.lower()
            if "mag" not in cl:
                continue
            if not np.any(np.isfinite(gpmag)):
                if ("g'" in c or "g′" in c or "gpm" in cl) and "r'" not in c and "r′" not in c:
                    try:
                        t = np.asarray(tab[c], dtype=float)
                        if np.any(np.isfinite(t)):
                            gpmag = t
                    except Exception:
                        pass
            if not np.any(np.isfinite(rpmag)):
                if "r'" in c or "r′" in c or "rpm" in cl:
                    try:
                        t = np.asarray(tab[c], dtype=float)
                        if np.any(np.isfinite(t)):
                            rpmag = t
                    except Exception:
                        pass

    bad_v = (vmag > 90.0) | ~np.isfinite(vmag)
    bad_b = (bmag > 90.0) | ~np.isfinite(bmag)
    bad_g = (gpmag > 90.0) | ~np.isfinite(gpmag)
    bad_r = (rpmag > 90.0) | ~np.isfinite(rpmag)

    ok_bv = (
        np.isfinite(bmag)
        & np.isfinite(vmag)
        & (vmag <= mag_limit_v)
        & ~bad_b
        & ~bad_v
    )
    ok_gr = (
        np.isfinite(gpmag)
        & np.isfinite(rpmag)
        & np.isfinite(vmag)
        & (vmag <= mag_limit_v)
        & ~bad_g
        & ~bad_r
        & ~bad_v
    )
    ok = ok_bv | ok_gr
    if not np.any(ok):
        return Table(
            names=("RAJ2000", "DEJ2000", "Bmag", "Vmag", "gpmag", "rpmag"),
            dtype=("f8", "f8", "f8", "f8", "f8", "f8"),
        )

    out = Table()
    out["RAJ2000"] = ra[ok]
    out["DEJ2000"] = dec[ok]
    out["Bmag"] = bmag[ok]
    out["Vmag"] = vmag[ok]
    out["gpmag"] = gpmag[ok]
    out["rpmag"] = rpmag[ok]
    return out


def compute_instrumental_and_catalog(
    fits_path: Path,
    ref_band_id: str,
    ebv: float,
    *,
    max_sep_arcsec: float = 2.0,
    dao_threshold_sigma: float = 5.0,
    mag_limit_gaia: float = 18.0,
    mag_limit_apass_v: float = 16.5,
    min_snr: float = MIN_SNR_FOR_TRANS_COEFF,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Retourne (x, y, mag_inst, mag_cat_obs, mag_cat_deredden, flux) pour chaque étoile appariée.

    On ne conserve que les étoiles pour lesquelles la magnitude instrumentale et la magnitude
    catalogue observée sont strictement supérieures à ``MAGNITUDE_MIN_FOR_COEFF`` (-15).

    Si ``min_snr`` > 0, on écarte les sources dont le SNR du flux net (ouverture moins ciel
    annulaire) est strictement inférieur à ce seuil. Le bruit est estimé à partir de l'écart-type
    par pixel ``std`` (``sigma_clipped_stats`` sur l'image), en supposant un bruit blanc :
    ``σ_net ≈ std × √(N_ap + N_ap² / N_an)`` avec ``N_ap``, ``N_an`` les aires en pixels.
    """
    try:
        from photutils.detection import DAOStarFinder
        from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
    except ImportError as e:
        raise RuntimeError("photutils est requis pour la détection et la photométrie.") from e

    from astropy.stats import sigma_clipped_stats

    fits_path = Path(fits_path)
    with fits.open(fits_path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header

    if data.ndim != 2:
        raise ValueError("Image 2D attendue.")

    try:
        wcs = WCS(header)
        if hasattr(wcs, "has_celestial") and not wcs.has_celestial:
            raise ValueError("WCS céleste absent : astrométrie requise.")
    except Exception as e:
        raise ValueError("WCS invalide ou absent : astrométrie requise.") from e

    _, median, std = sigma_clipped_stats(data, sigma_lower=3.0, sigma_upper=3.0)
    if not np.isfinite(std) or std <= 0:
        std = max(1e-6, float(np.nanstd(data)))
    fwhm = 3.0
    daofind = DAOStarFinder(fwhm=fwhm, threshold=dao_threshold_sigma * std)
    sources = daofind(data - median)
    if sources is None or len(sources) == 0:
        raise RuntimeError("Aucune source détectée (essayez une image calibrée / moins bruitée).")

    ny, nx = data.shape
    half = int(max(3.0, 2.5 * fwhm))
    r_ap = max(2.5, 1.8 * fwhm)
    r_in = r_ap + 2.0
    r_out = r_ap + 4.0

    xs = np.asarray(sources["xcentroid"], dtype=float)
    ys = np.asarray(sources["ycentroid"], dtype=float)
    fluxes = []
    valid_idx = []
    for i in range(len(xs)):
        xi, yi = int(round(xs[i])), int(round(ys[i]))
        if xi < half or xi >= nx - half or yi < half or yi >= ny - half:
            continue
        ap = CircularAperture((xs[i], ys[i]), r=r_ap)
        an = CircularAnnulus((xs[i], ys[i]), r_in=r_in, r_out=r_out)
        phot = aperture_photometry(data, ap, error=None)
        sky = aperture_photometry(data, an)
        bkg = float(sky["aperture_sum"][0] / an.area)
        fl = float(phot["aperture_sum"][0] - bkg * ap.area)
        if not (fl > 0 and np.isfinite(fl)):
            continue
        if min_snr > 0.0:
            n_ap = float(ap.area)
            n_an = float(an.area)
            if n_ap <= 0.0 or n_an <= 0.0:
                continue
            sigma_net = float(std) * float(np.sqrt(n_ap + (n_ap**2) / n_an))
            if not (sigma_net > 0.0 and np.isfinite(sigma_net)):
                continue
            snr = fl / sigma_net
            if not (np.isfinite(snr) and snr >= float(min_snr)):
                continue
        fluxes.append(fl)
        valid_idx.append(i)

    if not fluxes:
        raise RuntimeError("Photométrie d'ouverture : aucun flux positif.")

    xs = xs[valid_idx]
    ys = ys[valid_idx]
    fluxes = np.asarray(fluxes, dtype=float)
    mag_inst = -2.5 * np.log10(fluxes)

    coords = pixel_to_skycoord(xs, ys, wcs, origin=0)
    corners = np.array([[0, 0], [nx, 0], [0, ny], [nx, ny]])
    cc = pixel_to_skycoord(corners[:, 0], corners[:, 1], wcs, origin=0)
    center = SkyCoord(ra=cc.ra.mean(), dec=cc.dec.mean(), frame="icrs")
    sep = center.separation(cc).deg
    radius_deg = float(np.max(sep) * 1.05 + 0.02)

    max_sep = max_sep_arcsec * u.arcsec

    if ref_band_id in ("gaia_g", "gaia_bp", "gaia_rp"):
        cat = query_gaia_dr3_field(center, radius_deg, mag_limit=mag_limit_gaia)
        if len(cat) == 0:
            raise RuntimeError("Catalogue Gaia vide pour ce champ.")
        cat_coords = SkyCoord(ra=cat["RA_ICRS"] * u.deg, dec=cat["DE_ICRS"] * u.deg, frame="icrs")
        idx, d2d, _ = match_coordinates_sky(coords, cat_coords)
        m = d2d < max_sep
        if not np.any(m):
            raise RuntimeError("Aucun appariement Gaia (rayon trop strict ou champ hors catalogue).")
        idx_m = idx[m]
        G = np.asarray(cat["Gmag"][idx_m], dtype=float)
        BP = np.asarray(cat["BPmag"][idx_m], dtype=float)
        RP = np.asarray(cat["RPmag"][idx_m], dtype=float)
        mag_cat = reference_magnitude_gaia(G, BP, RP, ref_band_id)
    elif ref_band_id in ("johnson_b", "johnson_v", "johnson_rc"):
        cat = query_apass9_field(center, radius_deg, mag_limit_v=mag_limit_apass_v)
        if len(cat) == 0:
            raise RuntimeError(
                "Catalogue APASS vide pour ce champ (ou magnitude V limite trop stricte). "
                "Réduisez la profondeur demandée ou élargissez le rayon si possible."
            )
        cat_coords = SkyCoord(ra=cat["RAJ2000"] * u.deg, dec=cat["DEJ2000"] * u.deg, frame="icrs")
        idx, d2d, _ = match_coordinates_sky(coords, cat_coords)
        m = d2d < max_sep
        if not np.any(m):
            raise RuntimeError(
                "Aucun appariement APASS (positions ou rayon de recherche). "
                "Vérifiez l'astrométrie du champ et la couverture APASS."
            )
        idx_m = idx[m]
        if ref_band_id == "johnson_b":
            mag_cat = np.asarray(cat["Bmag"][idx_m], dtype=float)
        elif ref_band_id == "johnson_v":
            mag_cat = np.asarray(cat["Vmag"][idx_m], dtype=float)
        else:
            g = np.asarray(cat["gpmag"][idx_m], dtype=float)
            rr = np.asarray(cat["rpmag"][idx_m], dtype=float)
            mag_cat = johnson_rc_from_apass_gr_lupton2005(g, rr)
            if not np.any(np.isfinite(mag_cat)):
                raise RuntimeError(
                    "Impossible de calculer Johnson Rc : magnitudes APASS g′/r′ absentes ou invalides "
                    "pour les étoiles appariées (catalogue Vizier II/336)."
                )
    else:
        raise ValueError(f"Bande de référence inconnue : {ref_band_id}")

    Rl = _ebv_to_extinction_ebv(ref_band_id)
    mag_der = mag_cat - Rl * float(ebv)

    xs_out = xs[m]
    ys_out = ys[m]
    mi = mag_inst[m]
    fl = fluxes[m]

    floor = MAGNITUDE_MIN_FOR_COEFF
    keep = (
        np.isfinite(mi)
        & np.isfinite(mag_cat)
        & (mi > floor)
        & (mag_cat > floor)
    )
    if not np.any(keep):
        raise RuntimeError(
            f"Aucune étoile ne satisfait la condition magnitude instrumentale et catalogue > {floor} "
            "(seuil pour exclure les mesures trop brillantes / aberrantes)."
        )

    return (
        xs_out[keep],
        ys_out[keep],
        mi[keep],
        mag_cat[keep],
        mag_der[keep],
        fl[keep],
    )


def fit_transformation(mag_inst: np.ndarray, mag_cat_der: np.ndarray) -> tuple[float, float, float, int]:
    """Régression linéaire mag_cat = a * mag_inst + b ; retourne (a, b, rms, n)."""
    x = np.asarray(mag_inst, dtype=float)
    y = np.asarray(mag_cat_der, dtype=float)
    floor = MAGNITUDE_MIN_FOR_COEFF
    ok = np.isfinite(x) & np.isfinite(y) & (x > floor) & (y > floor)
    x, y = x[ok], y[ok]
    if len(x) < 5:
        raise RuntimeError("Au moins 5 étoiles appariées valides sont nécessaires pour l'ajustement.")
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    rms = float(np.std(y - pred))
    return float(a), float(b), rms, int(len(x))


def _norm_obs_filter(s: str) -> str:
    return re.sub(r"\s+", "", (s or "").strip().lower())


def filters_match_for_trans_coeff(obs_ui: str, obs_journal: str) -> bool:
    """Correspondance souple entre le filtre saisi (astéroïdes) et celui du journal réduction."""
    a = _norm_obs_filter(obs_ui)
    b = _norm_obs_filter(obs_journal)
    if not a or not b:
        return False
    return a == b or a in b or b in a


def ref_band_label_for_id(ref_band_id: str) -> str:
    for lab, bid in REFERENCE_BAND_CHOICES:
        if bid == ref_band_id:
            return lab
    return str(ref_band_id)


def ref_band_id_for_observer_filter(obs_filter: str) -> str:
    """
    Bande catalogue cible pour les coefficients de transformation, déduite du **filtre observateur**
    (champ « Filtre utilisé » / en-tête FITS).

    Correspondances : **G** → Gaia G ; **BP**, **g_bp**, etc. → Gaia BP ; **RP**, **g_rp** → Gaia RP ;
    **B** → Johnson B (APASS) ; **V** → Johnson V ; **R**, **Rc** (seuls ou explicites) → Johnson Rc.

    Filtre non reconnu : retourne ``gaia_g`` (comportement historique) avec trace en log.
    """
    n = _norm_obs_filter(obs_filter)
    if not n:
        return "gaia_g"
    # Sous-chaînes Gaia RP/BP avant tout (évite qu'une lettre « b » absorbe « bp »).
    if "rp" in n or n in ("grp", "g_rp", "gaia_rp"):
        return "gaia_rp"
    if "bp" in n or n in ("gbp", "g_bp", "gaia_bp"):
        return "gaia_bp"
    if n in (
        "johnson_rc",
        "johnsonr",
        "rc",
        "r_c",
        "cousinsr",
        "cousins_rc",
        "r",
        "rcousins",
    ):
        return "johnson_rc"
    if n in ("johnson_v", "johnsonv", "v"):
        return "johnson_v"
    if n in ("johnson_b", "johnsonb", "b"):
        return "johnson_b"
    if n in ("gaia_g", "gaiag", "gg", "g", "gaia"):
        return "gaia_g"
    logger.info(
        "Filtre observateur %r : aucune règle dédiée ; bande catalogue par défaut gaia_g pour les coefficients.",
        obs_filter,
    )
    return "gaia_g"


def pair_coefficient_csv_path(obs_filter: str, ref_band_id: str) -> Path:
    """Chemin du CSV paire ex. ``G__Gaia_G.csv`` pour filtre observateur ``G`` et bande ``gaia_g``."""
    root = transformation_storage_dir()
    obs_t = safe_coeff_csv_token(obs_filter)
    ref_t = safe_coeff_csv_token(ref_band_label_for_id(ref_band_id))
    return root / f"{obs_t}__{ref_t}.csv"


def _read_coefficient_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return [dict(r) for r in reader]


def _mask_within_2sigma_1d(vals: np.ndarray) -> np.ndarray:
    """Masque booléen : valeurs finies et dans [médiane − 2σ, médiane + 2σ] (σ = écart-type échantillon)."""
    v = np.asarray(vals, dtype=float)
    ok = np.isfinite(v)
    if int(np.count_nonzero(ok)) < 3:
        return ok
    vv = v[ok]
    med = float(np.median(vv))
    std = float(np.std(vv))
    if not np.isfinite(std) or std <= 0:
        return ok
    return ok & (np.abs(v - med) <= 2.0 * std)


def median_slope_intercept_after_2sigma_clip(
    rows: list[dict[str, str]],
    *,
    exclude_aggregate_marker: str = AGGREGATE_ROW_FITS_MARKER,
) -> tuple[float, float, int, float]:
    """
    Médiane de la pente et de l'ordonnée à l'origine après exclusion des valeurs > 2σ
    (indépendamment sur les pentes et sur les intercepts, intersection des masques).

    Returns
    -------
    slope_med, intercept_med, n_rows_used, mag_std_med
        ``mag_std_med`` = médiane des colonnes ``mag_std`` ou ``rms`` des lignes retenues
        (incertitude-type pour les magnitudes issues du fit), ou NaN si aucune valeur.
    """
    slopes: list[float] = []
    inter: list[float] = []
    rms_or_ms: list[float] = []
    for r in rows:
        if (r.get("fits") or "").strip() == exclude_aggregate_marker:
            continue
        try:
            a = float(r["slope"])
            b = float(r["intercept"])
        except (KeyError, ValueError, TypeError):
            continue
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        slopes.append(a)
        inter.append(b)
        ms = None
        for key in ("mag_std", "rms"):
            if key not in r or r[key] is None or str(r[key]).strip() == "":
                continue
            try:
                ms = float(r[key])
            except (ValueError, TypeError):
                continue
            if np.isfinite(ms):
                break
            ms = None
        rms_or_ms.append(float(ms) if ms is not None else np.nan)

    n = len(slopes)
    if n == 0:
        raise ValueError("Aucune ligne exploitable (slope/intercept numériques).")
    sa = np.asarray(slopes, dtype=float)
    sb = np.asarray(inter, dtype=float)
    ma = _mask_within_2sigma_1d(sa)
    mb = _mask_within_2sigma_1d(sb)
    m = ma & mb
    if not np.any(m):
        m = np.ones(n, dtype=bool)

    a_med = float(np.median(sa[m]))
    b_med = float(np.median(sb[m]))
    n_used = int(np.count_nonzero(m))
    rms_arr = np.asarray(rms_or_ms, dtype=float)[m]
    rms_fin = rms_arr[np.isfinite(rms_arr)]
    mag_std_med = float(np.median(rms_fin)) if rms_fin.size > 0 else float("nan")
    return a_med, b_med, n_used, mag_std_med


def append_median_2sigma_aggregate_row_to_pair_csv(
    obs_filter: str,
    ref_band_id: str,
    *,
    ref_band_label: str | None = None,
) -> Path:
    """
    Réécrit le CSV paire ``{filtre}__{bande}.csv`` : retire les anciennes lignes ``__AGG_MEDIAN_2SIG__``,
    recalcule la médiane (pente, intercept) après clipping 2σ sur les lignes restantes,
    puis ajoute **une** ligne de synthèse (dernière ligne = « la plus récente » pour le chargeur).

    La ligne de synthèse inclut ``mag_std`` (médiane des ``mag_std``/``rms`` des lignes retenues) lorsque
    l'en-tête du fichier est migré pour contenir cette colonne.
    """
    obs_filter = (obs_filter or "").strip() or "UNKNOWN"
    ref_lab = ref_band_label or ref_band_label_for_id(ref_band_id)
    path = pair_coefficient_csv_path(obs_filter, ref_band_id)
    rows = _read_coefficient_csv_rows(path)
    if not rows:
        raise FileNotFoundError(
            f"Aucune donnée dans {path.name}. Enregistrez au moins un coefficient "
            f"(Réduction → Sauvegarder le coefficient) pour ce filtre et cette bande."
        )

    a_med, b_med, n_used, mag_std_med = median_slope_intercept_after_2sigma_clip(rows)
    base_fieldnames = [
        "timestamp_utc",
        "obs_filter",
        "ref_band",
        "ref_band_id",
        "slope",
        "intercept",
        "rms",
        "n_stars",
        "ebv",
        "fits",
    ]
    fieldnames = list(base_fieldnames) + ["mag_std"]

    data_rows = [r for r in rows if (r.get("fits") or "").strip() != AGGREGATE_ROW_FITS_MARKER]
    for r in data_rows:
        if "mag_std" in fieldnames and "mag_std" not in r:
            r["mag_std"] = r.get("mag_std", "") or ""

    ts = datetime.now(timezone.utc).isoformat()
    rms_str = f"{mag_std_med:.6f}" if np.isfinite(mag_std_med) else ""
    mag_std_str = rms_str
    new_row: dict[str, str] = {
        "timestamp_utc": ts,
        "obs_filter": obs_filter,
        "ref_band": ref_lab,
        "ref_band_id": ref_band_id,
        "slope": f"{a_med:.6f}",
        "intercept": f"{b_med:.6f}",
        "rms": rms_str,
        "n_stars": str(n_used),
        "ebv": "",
        "fits": AGGREGATE_ROW_FITS_MARKER,
        "mag_std": mag_std_str,
    }
    for k in fieldnames:
        if k not in new_row:
            new_row[k] = ""

    out_rows: list[dict[str, str]] = []
    for r in data_rows:
        out = {fn: (r.get(fn) or "") for fn in fieldnames}
        out_rows.append(out)
    out_rows.append({fn: new_row.get(fn, "") for fn in fieldnames})

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    logger.info(
        "CSV paire agrégé médiane 2σ : %s (n_retenu=%s, a=%.6f, b=%.6f, mag_std=%s)",
        path,
        n_used,
        a_med,
        b_med,
        mag_std_str or "—",
    )
    return path


def _parse_mag_std_from_row(row: dict[str, str]) -> float | None:
    """Lit ``mag_std`` ou, à défaut, ``rms`` (RMS du fit réduction) comme incertitude magnitude."""
    raw_ms = (row.get("mag_std") or "").strip()
    if raw_ms:
        try:
            v = float(raw_ms)
            return v if np.isfinite(v) else None
        except (ValueError, TypeError):
            pass
    raw_r = (row.get("rms") or "").strip()
    if raw_r:
        try:
            v = float(raw_r)
            return v if np.isfinite(v) else None
        except (ValueError, TypeError):
            pass
    return None


def load_latest_transformation_coefficient(
    obs_filter: str,
    *,
    ref_band_id: str = "gaia_g",
    journal_path: Path | None = None,
) -> tuple[float, float, str, float | None] | None:
    """
    Dernière ligne pertinente : d'abord le CSV paire ``{filtre}__{bande}.csv`` s'il existe,
    sinon le journal ``coefficients_journal.csv``.

    Retourne ``(a, b, description_courte, mag_std)`` avec ``mag_cat ≈ a × m_inst + b`` (réduction).
    ``mag_std`` provient de la colonne du même nom si présente, sinon ``None`` (utiliser seulement
    l'erreur instrumentale côté flux).
    """
    obs_filter = (obs_filter or "").strip()
    if not obs_filter:
        return None

    pp = pair_coefficient_csv_path(obs_filter, ref_band_id)
    if pp.is_file():
        prow = _read_coefficient_csv_rows(pp)
        for row in reversed(prow):
            rid = (row.get("ref_band_id") or "").strip()
            if rid and rid != ref_band_id:
                continue
            jf = row.get("obs_filter") or ""
            if jf and not filters_match_for_trans_coeff(obs_filter, jf):
                continue
            try:
                a = float(row["slope"])
                b = float(row["intercept"])
            except (KeyError, ValueError, TypeError):
                continue
            if not (np.isfinite(a) and np.isfinite(b)):
                continue
            ts = (row.get("timestamp_utc") or "")[:19]
            rms = row.get("rms", "")
            desc = f"{pp.name} @ {ts}, RMS={rms}"
            mag_std = _parse_mag_std_from_row(row)
            return (a, b, desc, mag_std)

    jp = journal_path or (transformation_storage_dir() / "coefficients_journal.csv")
    if not jp.is_file():
        return None
    rows: list[dict[str, str]] = []
    with jp.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(dict(row))
    for row in reversed(rows):
        rid = (row.get("ref_band_id") or "").strip()
        if rid != ref_band_id:
            continue
        jf = row.get("obs_filter") or ""
        if not filters_match_for_trans_coeff(obs_filter, jf):
            continue
        try:
            a = float(row["slope"])
            b = float(row["intercept"])
        except (KeyError, ValueError, TypeError):
            continue
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        ts = (row.get("timestamp_utc") or "")[:19]
        rms = row.get("rms", "")
        desc = f"{jp.name} @ {ts}, RMS={rms}"
        return (a, b, desc, _parse_mag_std_from_row(row))
    return None


def append_coefficient_record(
    obs_filter: str,
    ref_band_label: str,
    ref_band_id: str,
    slope: float,
    intercept: float,
    rms: float,
    n_stars: int,
    ebv: float,
    fits_path: str,
) -> Path:
    root = transformation_storage_dir()
    ts = datetime.now(timezone.utc).isoformat()

    row = {
        "timestamp_utc": ts,
        "obs_filter": obs_filter,
        "ref_band": ref_band_label,
        "ref_band_id": ref_band_id,
        "slope": f"{slope:.6f}",
        "intercept": f"{intercept:.6f}",
        "rms": f"{rms:.6f}",
        "n_stars": str(n_stars),
        "ebv": f"{ebv:.5f}",
        "fits": fits_path,
    }
    journal = root / "coefficients_journal.csv"
    fieldnames = list(row.keys())
    write_header = not journal.exists()
    with journal.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            w.writeheader()
        w.writerow(row)

    pair_path = root / f"{safe_coeff_csv_token(obs_filter)}__{safe_coeff_csv_token(ref_band_label)}.csv"
    pair_new = not pair_path.exists() or pair_path.stat().st_size == 0
    with pair_path.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if pair_new:
            w.writeheader()
        w.writerow(row)

    return pair_path
