"""
Module de binning temporel pour les courbes de lumière d'exoplanètes.

Basé sur les études montrant qu'un temps d'intégration de ~1 min (comme Kepler)
n'affecte pas la précision des paramètres de transit comparé à des temps plus courts.

Le binning optimal dépend du temps d'exposition et permet de :
- Réduire le bruit tout en préservant la forme du transit
- Accélérer les calculs de modélisation
- Améliorer le SNR sans perte significative d'information
"""

import numpy as np
import logging
from typing import Tuple, Optional
from scipy import stats

logger = logging.getLogger(__name__)


def optimal_bin_time(exposure_time: float, 
                     cadence: Optional[float] = None,
                     target_bin_time: float = 60.0) -> float:
    """
    Calcule le temps de binning optimal basé sur le temps d'exposition.
    
    Parameters
    ----------
    exposure_time : float
        Temps d'exposition en secondes
    cadence : float, optional
        Temps entre deux expositions en secondes. Si None, assume cadence = exposure_time
    target_bin_time : float, optional
        Temps de binning cible en secondes (défaut: 60.0 comme Kepler)
        
    Returns
    -------
    float
        Temps de binning optimal en secondes
    """
    if cadence is None:
        cadence = exposure_time
    
    # Le binning optimal devrait être proche du temps d'exposition
    # mais pas trop court (minimum 3x l'exposition pour réduire efficacement le bruit)
    min_bin_time = max(3 * exposure_time, 10.0)  # Minimum 10 secondes
    
    # Si le temps d'exposition est déjà long, utiliser un binning plus petit
    if exposure_time >= target_bin_time / 2:
        optimal_bin = max(cadence, exposure_time)
    else:
        # Calculer le nombre de points à moyenner pour atteindre le temps cible
        n_bins = int(np.ceil(target_bin_time / cadence))
        optimal_bin = n_bins * cadence
    
    # S'assurer que le binning est au moins égal au minimum
    optimal_bin = max(optimal_bin, min_bin_time)
    
    return optimal_bin


def bin_lightcurve(time: np.ndarray, 
                   flux: np.ndarray,
                   flux_err: Optional[np.ndarray] = None,
                   bin_time: float = 60.0,
                   method: str = 'mean',
                   preserve_transit: bool = True,
                   transit_duration: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Effectue un binning temporel d'une courbe de lumière.
    
    Parameters
    ----------
    time : np.ndarray
        Tableau des temps (JD)
    flux : np.ndarray
        Tableau des flux
    flux_err : np.ndarray, optional
        Tableau des erreurs sur le flux
    bin_time : float
        Temps de binning en secondes
    method : str
        Méthode de binning : 'mean' (moyenne), 'median' (médiane), 'weighted' (pondéré)
    preserve_transit : bool
        Si True, utilise un binning plus fin pendant le transit pour préserver la forme
    transit_duration : float, optional
        Durée du transit en jours. Utilisé si preserve_transit=True
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        (time_binned, flux_binned, flux_err_binned)
    """
    if len(time) == 0:
        return time, flux, flux_err if flux_err is not None else np.array([])
    
    # Convertir bin_time de secondes en jours
    bin_time_days = bin_time / 86400.0
    
    # Si preserve_transit et transit_duration fourni, utiliser un binning adaptatif
    if preserve_transit and transit_duration is not None:
        return _adaptive_binning(time, flux, flux_err, bin_time_days, 
                                transit_duration, method)
    
    # Binning uniforme
    return _uniform_binning(time, flux, flux_err, bin_time_days, method)


def _uniform_binning(time: np.ndarray,
                     flux: np.ndarray,
                     flux_err: Optional[np.ndarray],
                     bin_time_days: float,
                     method: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Binning uniforme de la courbe de lumière.
    """
    if len(time) == 0:
        return time, flux, flux_err if flux_err is not None else np.array([])
    
    # Créer les bords des bins
    time_min = np.min(time)
    time_max = np.max(time)
    n_bins = int(np.ceil((time_max - time_min) / bin_time_days)) + 1
    
    bin_edges = np.linspace(time_min, time_max, n_bins)
    bin_indices = np.digitize(time, bin_edges) - 1
    bin_indices = np.clip(bin_indices, 0, len(bin_edges) - 2)
    
    # Calculer les valeurs binnées pour chaque bin
    time_binned = []
    flux_binned = []
    flux_err_binned = []
    
    for i in range(len(bin_edges) - 1):
        mask = bin_indices == i
        if not np.any(mask):
            continue
        
        time_bin = time[mask]
        flux_bin = flux[mask]
        
        # Temps binné = moyenne pondérée par le flux (ou simplement moyenne)
        time_binned.append(np.mean(time_bin))
        
        # Flux binné selon la méthode
        if method == 'median':
            flux_binned.append(np.median(flux_bin))
            if flux_err is not None:
                # Erreur de la médiane approximée comme 1.253 * σ / sqrt(N)
                flux_err_binned.append(1.253 * np.median(np.abs(flux_bin - np.median(flux_bin))) / np.sqrt(len(flux_bin)))
            else:
                flux_err_binned.append(np.nan)
        elif method == 'weighted' and flux_err is not None:
            # Binning pondéré par l'inverse de la variance
            weights = 1.0 / (flux_err[mask]**2)
            weights_sum = np.sum(weights)
            flux_binned.append(np.sum(flux_bin * weights) / weights_sum)
            flux_err_binned.append(1.0 / np.sqrt(weights_sum))
        else:  # method == 'mean' or default
            flux_binned.append(np.mean(flux_bin))
            if flux_err is not None:
                # Erreur de la moyenne = sqrt(Σσ²) / N
                flux_err_binned.append(np.sqrt(np.sum(flux_err[mask]**2)) / len(flux_bin))
            else:
                # Estimation de l'erreur depuis l'écart-type du bin
                flux_err_binned.append(np.std(flux_bin) / np.sqrt(len(flux_bin)))
    
    time_binned = np.array(time_binned)
    flux_binned = np.array(flux_binned)
    flux_err_binned = np.array(flux_err_binned) if flux_err_binned else None
    
    return time_binned, flux_binned, flux_err_binned


def _adaptive_binning(time: np.ndarray,
                      flux: np.ndarray,
                      flux_err: Optional[np.ndarray],
                      bin_time_days: float,
                      transit_duration_days: float,
                      method: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Binning adaptatif : binning fin pendant le transit, plus grossier ailleurs.
    """
    # Estimer la position du transit (minimum de flux)
    transit_idx = np.argmin(flux)
    transit_time = time[transit_idx]
    
    # Définir la région du transit (centre ± 1.5 × durée)
    transit_region = 1.5 * transit_duration_days
    
    # Binning fin pendant le transit (bin_time / 3)
    # Binning normal ailleurs
    fine_bin_time = bin_time_days / 3.0
    
    # Séparer les données
    in_transit = np.abs(time - transit_time) <= transit_region
    out_transit = ~in_transit
    
    time_binned = []
    flux_binned = []
    flux_err_binned = []
    
    # Binner la région hors transit
    if np.any(out_transit):
        t_out, f_out, e_out = _uniform_binning(
            time[out_transit], flux[out_transit], 
            flux_err[out_transit] if flux_err is not None else None,
            bin_time_days, method
        )
        time_binned.append(t_out)
        flux_binned.append(f_out)
        if e_out is not None:
            flux_err_binned.append(e_out)
    
    # Binner la région du transit avec un binning plus fin
    if np.any(in_transit):
        t_in, f_in, e_in = _uniform_binning(
            time[in_transit], flux[in_transit],
            flux_err[in_transit] if flux_err is not None else None,
            fine_bin_time, method
        )
        time_binned.append(t_in)
        flux_binned.append(f_in)
        if e_in is not None:
            flux_err_binned.append(e_in)
    
    # Concaténer et trier par temps
    if len(time_binned) > 0:
        time_binned = np.concatenate(time_binned)
        flux_binned = np.concatenate(flux_binned)
        flux_err_binned = np.concatenate(flux_err_binned) if flux_err_binned else None
        
        sort_idx = np.argsort(time_binned)
        time_binned = time_binned[sort_idx]
        flux_binned = flux_binned[sort_idx]
        if flux_err_binned is not None:
            flux_err_binned = flux_err_binned[sort_idx]
    else:
        time_binned = np.array([])
        flux_binned = np.array([])
        flux_err_binned = None
    
    return time_binned, flux_binned, flux_err_binned


def calculate_binning_statistics(time: np.ndarray,
                                 flux: np.ndarray,
                                 time_binned: np.ndarray,
                                 flux_binned: np.ndarray) -> dict:
    """
    Calcule des statistiques sur l'effet du binning.
    
    Parameters
    ----------
    time : np.ndarray
        Temps original
    flux : np.ndarray
        Flux original
    time_binned : np.ndarray
        Temps binné
    flux_binned : np.ndarray
        Flux binné
        
    Returns
    -------
    dict
        Dictionnaire avec les statistiques
    """
    stats_dict = {
        'n_points_original': len(time),
        'n_points_binned': len(time_binned),
        'compression_ratio': len(time) / len(time_binned) if len(time_binned) > 0 else 1.0,
        'time_span_original': np.max(time) - np.min(time),
        'time_span_binned': np.max(time_binned) - np.min(time_binned) if len(time_binned) > 0 else 0.0,
    }
    
    # SNR amélioration approximative
    if len(flux) > 0 and len(flux_binned) > 0:
        snr_original = np.mean(flux) / np.std(flux) if np.std(flux) > 0 else np.nan
        snr_binned = np.mean(flux_binned) / np.std(flux_binned) if np.std(flux_binned) > 0 else np.nan
        stats_dict['snr_original'] = snr_original
        stats_dict['snr_binned'] = snr_binned
        stats_dict['snr_improvement'] = snr_binned / snr_original if snr_original > 0 else np.nan
    
    return stats_dict
