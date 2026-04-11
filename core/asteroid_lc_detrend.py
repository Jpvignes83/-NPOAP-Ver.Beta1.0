# core/asteroid_lc_detrend.py
"""Détrendage lent pour courbes de lumière d'astéroïdes (séries longues, extinction, etc.)."""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

try:
    from scipy.signal import savgol_filter

    _HAS_SCIPY_SIGNAL = True
except ImportError:
    _HAS_SCIPY_SIGNAL = False


def _additive_poly_detrend(
    time: np.ndarray,
    flux: np.ndarray,
    degree: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Retire un polynôme en temps (centré) ; ramène le niveau médian du flux initial."""
    t = np.asarray(time, dtype=float)
    f = np.asarray(flux, dtype=float)
    n = len(t)
    if n <= degree + 1:
        return f.copy(), np.ones_like(f) * float(np.nanmedian(f))
    t0 = float(np.nanmedian(t))
    tc = t - t0
    X = np.vander(tc, N=degree + 1, increasing=True)
    mask = np.isfinite(t) & np.isfinite(f)
    if np.sum(mask) <= degree + 1:
        return f.copy(), np.ones_like(f) * float(np.nanmedian(f))
    coef, *_ = np.linalg.lstsq(X[mask], f[mask], rcond=None)
    trend = X @ coef
    med_f = float(np.nanmedian(f[mask]))
    f_out = f - trend + med_f
    return f_out, trend


def _savgol_additive_detrend(
    flux: np.ndarray,
    window_points: Optional[int],
    polyorder: int,
) -> Tuple[np.ndarray, np.ndarray, str]:
    if not _HAS_SCIPY_SIGNAL:
        raise RuntimeError("scipy.signal indisponible (savgol_filter).")
    f = np.asarray(flux, dtype=float)
    n = len(f)
    if n < 5:
        return f.copy(), np.ones_like(f) * float(np.nanmedian(f)), "Savgol ignoré (< 5 points)."
    if window_points is None:
        wl = max(5, (n // 6) | 1)
        wl = min(wl, n - (1 - n % 2))
        if wl % 2 == 0:
            wl -= 1
        wl = max(3, wl)
    else:
        wl = int(window_points)
        if wl % 2 == 0:
            wl += 1
        wl = int(np.clip(wl, 3, n - (1 - n % 2)))
    po = int(min(max(1, polyorder), wl - 1))
    smooth = savgol_filter(
        np.where(np.isfinite(f), f, np.nanmedian(f)),
        window_length=wl,
        polyorder=po,
        mode="interp",
    )
    smooth = np.asarray(smooth, dtype=float)
    med_s = float(np.nanmedian(smooth[np.isfinite(smooth)]))
    f_out = f - smooth + med_s
    return f_out, smooth, f"Savgol (fenêtre {wl} pts, ordre {po})."


def _wotan_detrend(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: Optional[np.ndarray],
    method: str,
    window_days: Optional[float],
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], np.ndarray, str]:
    try:
        import wotan
    except ImportError as e:
        raise RuntimeError("Le paquet « wotan » n'est pas installé.") from e

    from core.wotan_detrend import default_window_length

    t = np.asarray(time, dtype=float)
    f = np.asarray(flux, dtype=float)
    span = float(np.ptp(t[np.isfinite(t)])) if np.any(np.isfinite(t)) else 0.0
    if window_days is not None and np.isfinite(window_days) and float(window_days) > 0:
        wl = float(window_days)
        if span > 0:
            wl = float(np.clip(wl, 1e-4, 0.49 * span))
    else:
        wl = default_window_length(t)

    flatten_lc, trend_lc = wotan.flatten(
        t,
        f,
        window_length=wl,
        method=method,
        return_trend=True,
    )
    flatten_lc = np.asarray(flatten_lc, dtype=float)
    trend_lc = np.asarray(trend_lc, dtype=float)
    bad = ~np.isfinite(flatten_lc)
    if np.any(bad):
        flatten_lc = flatten_lc.copy()
        flatten_lc[bad] = f[bad]

    med_f = float(np.nanmedian(f[np.isfinite(f)]))
    med_fl = float(np.nanmedian(flatten_lc[np.isfinite(flatten_lc)]))
    if np.isfinite(med_fl) and abs(med_fl) > 1e-15:
        scale = med_f / med_fl
        f_out = flatten_lc * scale
    else:
        f_out = flatten_lc.copy()

    if flux_err is not None:
        e = np.asarray(flux_err, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(np.abs(f) > 1e-15, f_out / f, 1.0)
        ratio = np.where(np.isfinite(ratio), ratio, 1.0)
        e_out = e * ratio
    else:
        e_out = None

    msg = f"Wōtan ({method}, fenêtre ≈ {wl:.5g} j)."
    return f_out, e_out, trend_lc, msg


def detrend_asteroid_lc(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: Optional[np.ndarray],
    method: str,
    *,
    savgol_window_points: Optional[int] = None,
    savgol_polyorder: int = 2,
    wotan_window_days: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], str]:
    """
    Applique un détrendage lent sur une courbe (temps inchangé).

    Returns
    -------
    time_out, flux_out, flux_err_out, message
    """
    t = np.asarray(time, dtype=float)
    f = np.asarray(flux, dtype=float)
    fe = np.asarray(flux_err, dtype=float) if flux_err is not None else None
    method = (method or "none").strip().lower()

    if method in ("none", "", "aucun"):
        return t.copy(), f.copy(), fe.copy() if fe is not None else None, "Pas de détrendage."

    if method == "linear":
        f_out, _trend = _additive_poly_detrend(t, f, 1)
        return t.copy(), f_out, fe.copy() if fe is not None else None, "Tendance polynomiale degré 1 retirée (niveau médian conservé)."

    if method == "poly2":
        f_out, _trend = _additive_poly_detrend(t, f, 2)
        return t.copy(), f_out, fe.copy() if fe is not None else None, "Tendance polynomiale degré 2 retirée."

    if method == "poly3":
        f_out, _trend = _additive_poly_detrend(t, f, 3)
        return t.copy(), f_out, fe.copy() if fe is not None else None, "Tendance polynomiale degré 3 retirée."

    if method == "savgol":
        f_out, _sm, extra = _savgol_additive_detrend(f, savgol_window_points, savgol_polyorder)
        return t.copy(), f_out, fe.copy() if fe is not None else None, extra

    if method == "wotan_biweight":
        f_out, e_out, _tr, msg = _wotan_detrend(t, f, fe, "biweight", wotan_window_days)
        return t.copy(), f_out, e_out, msg

    if method == "wotan_median":
        f_out, e_out, _tr, msg = _wotan_detrend(t, f, fe, "median", wotan_window_days)
        return t.copy(), f_out, e_out, msg

    return t.copy(), f.copy(), fe.copy() if fe is not None else None, f"Méthode inconnue « {method} » — données inchangées."


def list_detrend_methods() -> Tuple[Tuple[str, str], ...]:
    """(clé interne, libellé UI)."""
    base = (
        ("none", "Aucun"),
        ("linear", "Polynôme degré 1 (linéaire)"),
        ("poly2", "Polynôme degré 2"),
        ("poly3", "Polynôme degré 3"),
        ("savgol", "Savitzky–Golay (lissage)"),
    )
    try:
        import wotan  # noqa: F401

        extra = (
            ("wotan_biweight", "Wōtan biweight"),
            ("wotan_median", "Wōtan median"),
        )
        return base + extra
    except ImportError:
        return base
