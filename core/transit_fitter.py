"""
Ajustement local d'un transit (modèle trapézoïdal + ligne de base affine) pour affiner Tc.

Utilise scipy.optimize.least_squares (numpy/scipy uniquement).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares

logger = logging.getLogger(__name__)


def trapezoid_shape(u: np.ndarray, half_width: float, ingress_half: float) -> np.ndarray:
    """
    Forme d'occultation : 0 hors transit, rampe d'ingress/egress, plateau à 1.
    u = |t - t0|, half_width = demi-largeur totale (contact extérieur), ingress_half = demi-ingress.
    """
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u, dtype=float)
    hw = float(max(half_width, 1e-9))
    ing = float(np.clip(ingress_half, 1e-9, hw * 0.49))
    inner = hw - ing
    mask_full = u <= inner
    out[mask_full] = 1.0
    mask_ing = (u > inner) & (u < hw)
    out[mask_ing] = (hw - u[mask_ing]) / ing
    return out


def model_flux(
    t: np.ndarray,
    t0: float,
    depth: float,
    a0: float,
    a1: float,
    width: float,
    ingress_ratio: float,
    t_ref: float,
) -> np.ndarray:
    """Flux modèle : baseline affine moins une bosse trapézoïdale (profondeur * shape)."""
    t = np.asarray(t, dtype=float)
    half_w = max(width / 2.0, 1e-9)
    ing_half = max(ingress_ratio * half_w, 1e-9)
    shape = trapezoid_shape(np.abs(t - t0), half_w, ing_half)
    return a0 + a1 * (t - t_ref) - depth * shape


def fit_trapezoid_transit(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: Optional[np.ndarray],
    t0_init: float,
    half_window: float,
    width_fixed: float,
    ingress_ratio: float = 0.25,
    fit_slope: bool = True,
) -> Dict[str, Any]:
    """
    Ajuste t0, profondeur, offset (et optionnellement pente) avec largeur d'ingress fixée.

    Parameters
    ----------
    time, flux : courbe découpée autour du transit (strictement croissant en temps recommandé).
    flux_err : incertitudes (si None, poids uniformes).
    t0_init : premier guess du milieu de transit.
    half_window : demi-largeur de la fenêtre (pour bornes sur t0, t0 ± half_window/2 utilisé ci-dessous non — bornes t0 ± 0.45*half_window).
    width_fixed : largeur totale du modèle trapèze (T14-like), fixe.
    ingress_ratio : fraction de la *demi-largeur* occupée par l'ingress (0–1).
    fit_slope : si False, pente baseline fixée à 0.

    Returns
    -------
    dict avec t0, t0_err, depth, width, ingress_ratio, a0, a1, rms, success, message
    """
    time = np.asarray(time, dtype=float).ravel()
    flux = np.asarray(flux, dtype=float).ravel()
    if flux_err is not None:
        ferr = np.asarray(flux_err, dtype=float).ravel()
    else:
        ferr = None

    ok = np.isfinite(time) & np.isfinite(flux)
    if ferr is not None and len(ferr) == len(time):
        ok &= np.isfinite(ferr) & (ferr > 0)
    time = time[ok]
    flux = flux[ok]
    if ferr is not None and len(ferr) == len(time):
        ferr = ferr[ok]
    else:
        ferr = None

    n = len(time)
    if n < 8:
        return {
            "success": False,
            "message": f"Pas assez de points valides ({n}).",
            "t0": float(t0_init),
            "t0_err": np.nan,
        }

    t_ref = float(np.mean(time))
    hw = float(max(half_window, 1e-6))
    t_lo = float(t0_init - 0.45 * hw)
    t_hi = float(t0_init + 0.45 * hw)

    med = float(np.nanmedian(flux))
    fmin = float(np.nanmin(flux))
    depth0 = float(np.clip(med - fmin, 1e-6, 0.99))
    a0_0 = float(med)
    a1_0 = 0.0

    w_fix = float(max(width_fixed, 1e-6))
    ing_r = float(np.clip(ingress_ratio, 0.05, 0.45))

    if fit_slope:

        def resid(x):
            t0, depth, a0, a1 = x
            m = model_flux(time, t0, depth, a0, a1, w_fix, ing_r, t_ref)
            r = flux - m
            if ferr is not None:
                return r / ferr
            return r

        x0 = np.array([t0_init, depth0, a0_0, a1_0], dtype=float)
        lb = np.array([t_lo, 1e-8, -np.inf, -np.inf], dtype=float)
        ub = np.array([t_hi, 1.5, np.inf, np.inf], dtype=float)
    else:

        def resid(x):
            t0, depth, a0 = x
            m = model_flux(time, t0, depth, a0, 0.0, w_fix, ing_r, t_ref)
            r = flux - m
            if ferr is not None:
                return r / ferr
            return r

        x0 = np.array([t0_init, depth0, a0_0], dtype=float)
        lb = np.array([t_lo, 1e-8, -np.inf], dtype=float)
        ub = np.array([t_hi, 1.5, np.inf], dtype=float)

    try:
        res = least_squares(resid, x0, bounds=(lb, ub), method="trf", max_nfev=400)
    except Exception as e:
        logger.exception("least_squares transit")
        return {"success": False, "message": str(e), "t0": float(t0_init), "t0_err": np.nan}

    if not res.success:
        return {
            "success": False,
            "message": res.message or "Échec du fit",
            "t0": float(x0[0]),
            "t0_err": np.nan,
        }

    if fit_slope:
        t0_f, depth_f, a0_f, a1_f = res.x
    else:
        t0_f, depth_f, a0_f = res.x
        a1_f = 0.0

    npar = len(res.x)
    rms = float(np.sqrt(np.mean(res.fun**2)))
    t0_err = np.nan
    try:
        J = res.jac
        if J is not None and J.shape[0] >= J.shape[1]:
            JTJ = J.T @ J
            cov = np.linalg.inv(JTJ)
            chi2 = float(np.sum(res.fun**2))
            dof = max(len(res.fun) - npar, 1)
            scale = chi2 / dof
            t0_err = float(np.sqrt(max(cov[0, 0] * scale, 0.0)))
    except Exception:
        pass

    if not np.isfinite(t0_err) or t0_err <= 0:
        t0_err = float(hw / max(n, 1) * rms / max(depth_f, 1e-6) * 0.1)

    return {
        "success": True,
        "message": res.message,
        "t0": float(t0_f),
        "t0_err": float(t0_err),
        "depth": float(depth_f),
        "width": w_fix,
        "ingress_ratio": ing_r,
        "a0": float(a0_f),
        "a1": float(a1_f),
        "rms": rms,
        "t_ref": t_ref,
    }


def slice_window(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: Optional[np.ndarray],
    center: float,
    half_width: float,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Extrait les points dans [center - half_width, center + half_width]."""
    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)
    m = np.abs(time - center) <= half_width
    t = time[m]
    f = flux[m]
    e = None
    if flux_err is not None:
        fe = np.asarray(flux_err, dtype=float)
        if len(fe) == len(time):
            e = fe[m]
    return t, f, e
