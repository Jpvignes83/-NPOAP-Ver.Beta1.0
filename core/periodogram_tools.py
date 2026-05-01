import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle, BoxLeastSquares
from scipy.signal import find_peaks
import warnings
import logging

logger = logging.getLogger(__name__)


def _prepare_lc_arrays(time, flux, sort_time=True):
    """Convertit en float64, retire NaN / Inf, optionnellement trie par temps."""
    t = np.asarray(time, dtype=float)
    y = np.asarray(flux, dtype=float)
    mask = np.isfinite(t) & np.isfinite(y)
    t, y = t[mask], y[mask]
    if len(t) == 0:
        return t, y
    if sort_time:
        order = np.argsort(t)
        t, y = t[order], y[order]
    return t, y


def _poly_detrend(time, flux, degree):
    """
    Retire un trend polynomial en temps (temps recentré et mis à l’échelle pour la stabilité numérique).
    À appliquer sur le flux brut ou normalisé avant normalisation médiane pour le BLS.
    """
    t = np.asarray(time, dtype=float)
    y = np.asarray(flux, dtype=float)
    degree = int(degree)
    if degree < 1:
        return y
    n = len(t)
    if n < degree + 1:
        logger.warning(
            "[Périodogramme] Détrend ignoré: %s points < degré+1 (%s)",
            n,
            degree + 1,
        )
        return y
    span = float(np.ptp(t))
    if not np.isfinite(span) or span <= 0:
        return y
    tn = (t - np.mean(t)) / (0.5 * span)
    try:
        coeffs = np.polyfit(tn, y, degree)
    except (np.linalg.LinAlgError, ValueError, TypeError) as e:
        logger.warning("[Périodogramme] polyfit détrend échoué: %s", e)
        return y
    trend = np.polyval(coeffs, tn)
    return y - trend


def min_duration_frac_estimate(time, min_period, max_period):
    """
    Borne inférieure réaliste de la fraction d'orbite occupée par un transit bref,
    dérivée de la cadence médiane (pour calibrer phase_box_size dans Plavchan).
    """
    del max_period  # réservé pour extensions (p.ex. durée max vs P_max)
    span = float(np.max(time) - np.min(time))
    dt = float(np.median(np.diff(time))) if len(time) > 1 else span
    if not np.isfinite(dt) or dt <= 0:
        dt = span / max(len(time) - 1, 1)
    dur_guess = max(3.0 * dt, 2.0 / 24.0)
    return float(min(0.2, dur_guess / min_period))


def run_lomb_scargle(time, flux, min_period=0.1, max_period=20.0, samples_per_peak=15):
    """
    Périodogramme de Lomb-Scargle (Adaptation NASA).
    Utilise la normalisation 'standard' (inverse de la variance).
    Flux recentré (soustraction de la médiane) pour renforcer les signatures
    quasi sinusoïdales / harmoniques sans changer l’échelle des pics relatifs.
    """
    logger.debug("[Périodogramme] run_lomb_scargle n=%s min_period=%s max_period=%s", len(time), min_period, max_period)
    time, flux = _prepare_lc_arrays(time, flux)
    if len(time) < 3:
        return np.array([]), np.array([]), 0.0
    flux_work = flux - np.median(flux)
    ls = LombScargle(time, flux_work, normalization='standard')
    
    # Grille de fréquence uniforme
    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    
    frequency, power = ls.autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
        samples_per_peak=samples_per_peak,
    )
    
    period = 1.0 / frequency
    
    best_idx = np.argmax(power)
    best_period = period[best_idx]
    logger.info("[Périodogramme] Lomb-Scargle terminé, meilleure période=%.6f j", best_period)
    return period, power, best_period


def run_bls(
    time,
    flux,
    min_period=0.1,
    max_period=20.0,
    n_durations=24,
    n_periods=50000,
    detrend=False,
    detrend_degree=2,
):
    """
    Box-fitting Least Squares (BLS), adapté aux transits (Kovács et al. 2002).

    Prétraitement léger : médiane unitaire (photométrie relative), temps trié,
    comme recommandé dans les exercices comparatifs CoRoT / analyses avec variabilité
    stellaire (voir aussi Defaÿ et al. 2001 sur la séparation signal / bruit).

    Si ``detrend`` est True, un polynôme de temps de degré ``detrend_degree`` (1–5)
    est soustrait avant la normalisation médiane (réduit dérives / variabilité lente).

    Grille de périodes logarithmique (échantillonnage plus fin aux courtes périodes),
    conforme aux usages du type « log-frequency » dans les blind tests transit.

    La puissance renvoyée pour le tracé est, pour chaque période testée, le maximum
    sur la grille des durées de transit — obligatoire car ``power`` Astropy est 2D.
    """
    logger.debug(
        "[Périodogramme] run_bls n=%s min_period=%s max_period=%s detrend=%s",
        len(time),
        min_period,
        max_period,
        detrend,
    )
    time, flux = _prepare_lc_arrays(time, flux)
    if len(time) < 3:
        return np.array([]), np.array([]), 0.0

    if detrend:
        deg = int(detrend_degree)
        deg = max(1, min(5, deg))
        flux = _poly_detrend(time, flux, deg)
        logger.info("[Périodogramme] BLS : détrend polynomial actif (degré %s)", deg)

    med = np.median(flux)
    if not np.isfinite(med) or med == 0.0:
        flux_norm = flux - np.mean(flux)
    else:
        flux_norm = flux / med

    bls = BoxLeastSquares(time, flux_norm)

    if min_period <= 0 or max_period <= 0 or min_period >= max_period:
        raise ValueError("run_bls: attendu 0 < min_period < max_period")
    periods = np.geomspace(min_period, max_period, int(n_periods))

    span = float(np.max(time) - np.min(time))
    cadence_floor = span / max(len(time) - 1, 1)
    min_duration = max(min_period * 0.005, cadence_floor * 2.0)
    limit_physical = min_period * 0.5
    target_max = max_period * 0.12
    max_duration = min(target_max, limit_physical)

    if min_duration >= max_duration:
        max_duration = min_duration * 1.25

    durations = np.geomspace(min_duration, max_duration, int(n_durations))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = bls.power(periods, durations)

    power = np.asarray(results.power, dtype=float)
    if power.ndim != 2:
        power = power.reshape(len(periods), len(durations))

    finite = np.isfinite(power)
    if not np.any(finite):
        return periods, np.zeros_like(periods), 0.0

    power_masked = np.where(finite, power, np.nan)
    power_1d = np.nanmax(power_masked, axis=1)
    power_1d = np.nan_to_num(power_1d, nan=0.0, posinf=0.0, neginf=0.0)

    best_ij = np.unravel_index(np.nanargmax(np.where(finite, power, -np.inf)), power.shape)
    best_period = float(periods[best_ij[0]])
    best_dur_h = float(durations[best_ij[1]]) * 24.0
    logger.info(
        "[Périodogramme] BLS terminé, meilleure période=%.6f j (~durée %.3f h sur grille)",
        best_period,
        best_dur_h,
    )
    return periods, power_1d, best_period


def run_plavchan(time, flux, min_period=0.1, max_period=20.0, n_periods=5000, phase_box_size=None):
    """
    Périodogramme de Plavchan (2008), grille de périodes log-uniforme.
    Si phase_box_size est None, une largeur de phase ~max(0.015, min_frac) est utilisée,
    min_frac étant une borne inférieure réaliste pour la fraction orbitale occupée par un transit court.
    """
    logger.debug("[Périodogramme] run_plavchan n=%s min_period=%s max_period=%s", len(time), min_period, max_period)
    time, flux = _prepare_lc_arrays(time, flux)
    if len(time) < 3:
        return np.array([]), np.array([]), 0.0

    if min_period <= 0 or max_period <= 0 or min_period >= max_period:
        raise ValueError("run_plavchan: attendu 0 < min_period < max_period")
    periods = np.geomspace(min_period, max_period, int(n_periods))

    if phase_box_size is None:
        min_frac_orbit = min_duration_frac_estimate(time, min_period, max_period)
        phase_box_size = float(max(0.015, min(0.12, min_frac_orbit * 1.5)))
    powers = np.zeros(len(periods))
    
    flux_mean = np.mean(flux)
    normalization = np.sum((flux - flux_mean)**2)
    
    for i, P in enumerate(periods):
        phases = (time % P) / P
        sort_order = np.argsort(phases)
        ph_sorted = phases[sort_order]
        fl_sorted = flux[sort_order]
        
        ph_ext = np.concatenate([ph_sorted - 1.0, ph_sorted, ph_sorted + 1.0])
        fl_ext = np.concatenate([fl_sorted, fl_sorted, fl_sorted])
        
        w = phase_box_size
        left_edges = ph_sorted - w / 2.0
        right_edges = ph_sorted + w / 2.0
        
        idx_L = np.searchsorted(ph_ext, left_edges, side='left')
        idx_R = np.searchsorted(ph_ext, right_edges, side='right')
        
        cum_flux = np.insert(np.cumsum(fl_ext), 0, 0.0)
        
        sum_in_window = cum_flux[idx_R] - cum_flux[idx_L]
        count_in_window = idx_R - idx_L
        
        valid = count_in_window > 0
        smoothed = np.zeros_like(fl_sorted)
        smoothed[valid] = sum_in_window[valid] / count_in_window[valid]
        smoothed[~valid] = fl_sorted[~valid] 
        
        residuals_sq = (fl_sorted - smoothed)**2
        ssr = np.sum(residuals_sq)
        
        if ssr > 0:
            powers[i] = normalization / ssr
        else:
            powers[i] = 0 
            
    best_idx = np.argmax(powers)
    best_period = periods[best_idx]
    logger.info("[Périodogramme] Plavchan terminé, meilleure période=%.6f j", best_period)
    return periods, powers, best_period