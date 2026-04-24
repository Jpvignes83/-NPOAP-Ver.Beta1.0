from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


@dataclass
class PTCResult:
    means_adu: np.ndarray
    vars_adu2: np.ndarray
    slope: float
    intercept: float
    gain_e_per_adu: float | None
    readnoise_e: float | None
    used_points: int
    total_points: int


def _roi_slices(shape: tuple[int, int], frac: float = 0.5) -> tuple[slice, slice]:
    ny, nx = int(shape[0]), int(shape[1])
    frac = float(np.clip(frac, 0.1, 1.0))
    hy = max(8, int(0.5 * ny * frac))
    hx = max(8, int(0.5 * nx * frac))
    cy, cx = ny // 2, nx // 2
    return slice(max(0, cy - hy), min(ny, cy + hy)), slice(max(0, cx - hx), min(nx, cx + hx))


def _read_data(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=False) as hdul:
        return np.asarray(hdul[0].data, dtype=float), hdul[0].header


def _master_bias(bias_files: list[Path]) -> np.ndarray | None:
    if not bias_files:
        return None
    stack = []
    for p in bias_files:
        data, _ = _read_data(p)
        stack.append(data)
    if not stack:
        return None
    return np.median(np.stack(stack, axis=0), axis=0)


def _group_by_exptime(flat_files: list[Path]) -> dict[float, list[Path]]:
    groups: dict[float, list[Path]] = {}
    for p in flat_files:
        try:
            hdr = fits.getheader(p, 0)
            exp = float(hdr.get("EXPTIME", hdr.get("EXPOSURE", 0.0)) or 0.0)
        except Exception:
            exp = 0.0
        groups.setdefault(exp, []).append(p)
    return dict(sorted(groups.items(), key=lambda kv: kv[0]))


def analyze_ptc(
    flat_files: list[str | Path],
    bias_files: list[str | Path] | None = None,
    *,
    roi_fraction: float = 0.5,
) -> PTCResult:
    flats = [Path(f) for f in flat_files]
    biases = [Path(f) for f in (bias_files or [])]
    if len(flats) < 2:
        raise ValueError("PTC: au moins 2 flats sont nécessaires.")

    master_bias = _master_bias(biases)
    groups = _group_by_exptime(flats)
    means: list[float] = []
    variances: list[float] = []

    for _, files in groups.items():
        n_pairs = len(files) // 2
        if n_pairs < 1:
            continue
        for i in range(n_pairs):
            p1, p2 = files[2 * i], files[2 * i + 1]
            d1, _ = _read_data(p1)
            d2, _ = _read_data(p2)
            if master_bias is not None and master_bias.shape == d1.shape:
                d1 = d1 - master_bias
                d2 = d2 - master_bias
            ys, xs = _roi_slices(d1.shape, frac=roi_fraction)
            a = d1[ys, xs]
            b = d2[ys, xs]
            mean = float(np.nanmedian(0.5 * (a + b)))
            diff = (a - b) / np.sqrt(2.0)
            var = float(np.nanvar(diff, ddof=1))
            if np.isfinite(mean) and np.isfinite(var) and var > 0:
                means.append(mean)
                variances.append(var)

    if len(means) < 3:
        raise ValueError("PTC: pas assez de points valides (>=3 requis).")

    x = np.asarray(means, dtype=float)
    y = np.asarray(variances, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    # Fit sur la zone centrale (évite bas signal et haut signal potentiellement non linéaire)
    qlo, qhi = np.quantile(x, [0.15, 0.85])
    fit_mask = (x >= qlo) & (x <= qhi)
    if np.sum(fit_mask) < 3:
        fit_mask = np.ones_like(x, dtype=bool)

    slope, intercept = np.polyfit(x[fit_mask], y[fit_mask], 1)
    gain = (1.0 / slope) if slope > 0 else None
    readnoise = (gain * np.sqrt(intercept)) if (gain is not None and intercept > 0) else None

    return PTCResult(
        means_adu=x,
        vars_adu2=y,
        slope=float(slope),
        intercept=float(intercept),
        gain_e_per_adu=(float(gain) if gain is not None and np.isfinite(gain) else None),
        readnoise_e=(float(readnoise) if readnoise is not None and np.isfinite(readnoise) else None),
        used_points=int(np.sum(fit_mask)),
        total_points=len(x),
    )

