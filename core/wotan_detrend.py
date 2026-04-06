"""
Détrendage optionnel avec Wotan (Hippke et al. 2019) sur un segment local de LC.

Référence : https://github.com/hippke/wotan
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def is_wotan_available() -> bool:
    try:
        import wotan  # noqa: F401

        return True
    except ImportError:
        return False


def default_window_length(time: np.ndarray) -> float:
    """Longueur de fenêtre Wotan (jours) adaptée à l’étendue temporelle du segment."""
    t = np.asarray(time, dtype=float)
    ok = np.isfinite(t)
    if not np.any(ok):
        return 0.1
    span = float(np.ptp(t[ok]))
    if span <= 0:
        return 0.05
    # ~22 % de la plage, bornée pour garder plusieurs points par fenêtre
    wl = 0.22 * span
    wl = max(wl, min(0.02, 0.5 * span))
    wl = min(wl, 0.48 * span)
    return float(max(wl, 1e-4))


def transit_boolean_mask(
    time: np.ndarray,
    tc: float,
    t14_full_width: float,
    edge_factor: float = 1.15,
) -> np.ndarray:
    """
    True pour les échantillons à l’intérieur du transit.

    Parameters
    ----------
    t14_full_width : largeur T14 (contact à contact), en unités de *time* (souvent jours).
    edge_factor : marge (>1) pour couvrir un léger décalage de Tc (TTV).
    """
    half = 0.5 * float(max(t14_full_width, 1e-12)) * float(edge_factor)
    return np.abs(np.asarray(time, dtype=float) - float(tc)) <= half


def flatten_segment_with_transit_mask(
    time: np.ndarray,
    flux: np.ndarray,
    tc: float,
    t14_width: float,
    window_length: Optional[float] = None,
    method: str = "biweight",
    break_tolerance: float = 0.0,
    edge_cutoff: float = 0.0,
    mask_edge_factor: float = 1.15,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Applique ``wotan.flatten`` sur le segment avec ``mask`` = transit (True = dans le transit).

    Returns
    -------
    flux_flat, trend, in_transit_mask, scale
        *flux_flat* est divisé par la médiane hors transit pour viser ~1 hors transit.
        *scale* est ce facteur (médiane avant division), pour rescaling des erreurs.
    """
    import wotan

    t = np.asarray(time, dtype=float)
    f = np.asarray(flux, dtype=float)
    wl = float(window_length) if window_length is not None else default_window_length(t)
    span = float(np.ptp(t[np.isfinite(t)])) if np.any(np.isfinite(t)) else 0.0
    if span > 0:
        wl = float(np.clip(wl, 1e-4, 0.49 * span))

    in_transit = transit_boolean_mask(t, tc, t14_width, edge_factor=mask_edge_factor)

    flatten_lc, trend_lc = wotan.flatten(
        t,
        f,
        window_length=wl,
        method=method,
        mask=in_transit,
        return_trend=True,
        break_tolerance=break_tolerance,
        edge_cutoff=edge_cutoff,
    )

    flatten_lc = np.asarray(flatten_lc, dtype=float)
    trend_lc = np.asarray(trend_lc, dtype=float)

    bad = ~np.isfinite(flatten_lc)
    if np.any(bad):
        flatten_lc = flatten_lc.copy()
        flatten_lc[bad] = f[bad]

    oot = (~in_transit) & np.isfinite(flatten_lc)
    med = np.nanmedian(flatten_lc[oot]) if np.any(oot) else np.nanmedian(flatten_lc)
    if med is None or not np.isfinite(med) or abs(med) < 1e-15:
        med = 1.0
    flux_norm = flatten_lc / med

    return flux_norm, trend_lc, in_transit, float(med)
