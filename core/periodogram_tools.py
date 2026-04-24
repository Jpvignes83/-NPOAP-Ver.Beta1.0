import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle, BoxLeastSquares
from scipy.signal import find_peaks
import warnings
import logging

logger = logging.getLogger(__name__)


def run_lomb_scargle(time, flux, min_period=0.1, max_period=20.0):
    """
    Périodogramme de Lomb-Scargle (Adaptation NASA).
    Utilise la normalisation 'standard' (inverse de la variance).
    """
    logger.debug("[Périodogramme] run_lomb_scargle n=%s min_period=%s max_period=%s", len(time), min_period, max_period)
    ls = LombScargle(time, flux, normalization='standard')
    
    # Grille de fréquence uniforme
    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    
    frequency, power = ls.autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
        samples_per_peak=10
    )
    
    period = 1.0 / frequency
    
    best_idx = np.argmax(power)
    best_period = period[best_idx]
    logger.info("[Périodogramme] Lomb-Scargle terminé, meilleure période=%.6f j", best_period)
    return period, power, best_period


def run_bls(time, flux, min_period=0.1, max_period=20.0, n_durations=20):
    """
    Box-fitting Least Squares (BLS) - Adapté NASA.
    Optimisé pour les transits.
    """
    logger.debug("[Périodogramme] run_bls n=%s min_period=%s max_period=%s", len(time), min_period, max_period)
    bls = BoxLeastSquares(time, flux)
    
    # Grille de Périodes
    periods = np.linspace(min_period, max_period, 50000)
    
    # --- CORRECTION DU BUG ---
    # La durée max ne doit pas dépasser la plus petite période testée.
    # Sinon Astropy lève "max transit duration must be shorter than min period"
    
    # Durée min : 1% de la petite période
    min_duration = min_period * 0.01
    
    # Durée max : On prend le MINIMUM entre :
    # 1. 10% de la période MAX (ce qu'on voudrait idéalement)
    # 2. 50% de la période MIN (la limite physique pour ne pas crasher)
    limit_physical = min_period * 0.5 
    target_max = max_period * 0.1
    
    max_duration = min(target_max, limit_physical)
    
    # Sécurité si min_duration >= max_duration (cas limites)
    if min_duration >= max_duration:
        max_duration = min_duration * 1.1

    # Création de la grille de durées (log-uniforme)
    durations = np.geomspace(min_duration, max_duration, n_durations)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = bls.power(periods, durations)
    
    mask = np.isfinite(results.power)
    clean_power = results.power[mask]
    clean_period = results.period[mask]
    
    if len(clean_power) == 0:
        return periods, np.zeros_like(periods), 0.0
        
    best_idx = np.argmax(clean_power)
    best_period = clean_period[best_idx]
    logger.info("[Périodogramme] BLS terminé, meilleure période=%.6f j", best_period)
    return clean_period, clean_power, best_period


def run_plavchan(time, flux, min_period=0.1, max_period=20.0, n_periods=5000, phase_box_size=0.05):
    """
    Périodogramme de Plavchan (2008).
    """
    logger.debug("[Périodogramme] run_plavchan n=%s min_period=%s max_period=%s", len(time), min_period, max_period)
    periods = np.linspace(min_period, max_period, n_periods)
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