"""
Module pour l'estimation des masses et rayons stellaires basé sur
Enoch et al. (2010), A&A, 510, A21.

Ce module implémente les équations de calibration qui permettent de calculer
la masse et le rayon stellaires à partir de T_eff, log ρ (densité stellaire)
et [Fe/H], évitant ainsi la dépendance aux modèles évolutifs.

Référence:
Enoch, B., Collier Cameron, A., Parley, N. R., & Hebb, L. (2010)
"An Improved Method for Estimating the Masses of Stars with Transiting Planets"
Astronomy & Astrophysics, 510, A21
DOI: 10.1051/0004-6361/200912675
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)

# Coefficients de calibration pour la masse (Table 1 de l'article)
# log M = a₁ + a₂X + a₃X² + a₄log ρ + a₅(log ρ)² + a₆(log ρ)³ + a₇[Fe/H]
MASS_COEFFICIENTS = {
    'a1': 0.458,  # Constante
    'a2': 1.430,  # X
    'a3': 0.329,  # X²
    'a4': 0.042,  # log ρ
    'a5': 0.067,  # (log ρ)²
    'a6': 0.010,  # (log ρ)³
    'a7': 0.044,  # [Fe/H]
}

MASS_COEFFICIENT_ERRORS = {
    'a1': 0.017,
    'a2': 0.019,
    'a3': 0.128,
    'a4': 0.021,
    'a5': 0.019,
    'a6': 0.004,
    'a7': 0.019,
}

# Coefficients de calibration pour le rayon (Table 1 de l'article)
# log R = b₁ + b₂X + b₃log ρ + b₄[Fe/H]
RADIUS_COEFFICIENTS = {
    'b1': 0.150,  # Constante
    'b2': 0.434,  # X
    'b3': 0.381,  # log ρ
    'b4': 0.012,  # [Fe/H]
}

RADIUS_COEFFICIENT_ERRORS = {
    'b1': 0.002,
    'b2': 0.005,
    'b3': 0.002,
    'b4': 0.004,
}

# Précision intrinsèque de la calibration (scatter observé)
SIGMA_LOG_M = 0.023  # dex
SIGMA_LOG_R = 0.009  # dex


def calculate_stellar_mass(teff, log_rho, feh, teff_err=None, log_rho_err=None, feh_err=None):
    """
    Calcule la masse stellaire en utilisant la calibration d'Enoch et al. (2010).
    
    Équation: log M = a₁ + a₂X + a₃X² + a₄log ρ + a₅(log ρ)² + a₆(log ρ)³ + a₇[Fe/H]
    où X = log(T_eff) - 4.1
    
    Parameters
    ----------
    teff : float
        Température effective (en K)
    log_rho : float
        Logarithme de la densité stellaire (log₁₀(ρ) en g/cm³)
    feh : float
        Métallicité [Fe/H] (décimal)
    teff_err : float, optional
        Erreur sur T_eff (en K). Si fourni, calcule l'erreur sur M
    log_rho_err : float, optional
        Erreur sur log ρ. Si fourni, calcule l'erreur sur M
    feh_err : float, optional
        Erreur sur [Fe/H]. Si fourni, calcule l'erreur sur M
    
    Returns
    -------
    mass : float
        Masse stellaire (en M☉)
    mass_err : float, optional
        Erreur sur la masse (en M☉). Retourné seulement si toutes les erreurs sont fournies
    """
    # Vérifications
    if teff <= 0:
        raise ValueError(f"T_eff doit être positif, reçu: {teff}")
    if not np.isfinite(log_rho):
        raise ValueError(f"log ρ doit être fini, reçu: {log_rho}")
    
    # Calculer X = log(T_eff) - 4.1
    X = np.log10(teff) - 4.1
    X2 = X**2
    
    # Calculer les termes de log ρ
    log_rho2 = log_rho**2
    log_rho3 = log_rho**3
    
    # Calculer log M
    log_M = (MASS_COEFFICIENTS['a1'] +
             MASS_COEFFICIENTS['a2'] * X +
             MASS_COEFFICIENTS['a3'] * X2 +
             MASS_COEFFICIENTS['a4'] * log_rho +
             MASS_COEFFICIENTS['a5'] * log_rho2 +
             MASS_COEFFICIENTS['a6'] * log_rho3 +
             MASS_COEFFICIENTS['a7'] * feh)
    
    # Convertir en masse (M☉)
    mass = 10.0**log_M
    
    # Calcul de l'erreur si toutes les erreurs sont fournies
    if teff_err is not None and log_rho_err is not None and feh_err is not None:
        # Propagation d'erreurs (approximation linéaire)
        # d(log M) = (∂log M/∂X) * dX + (∂log M/∂log ρ) * d(log ρ) + (∂log M/∂[Fe/H]) * d[Fe/H]
        # + erreurs sur les coefficients
        
        # Dérivées partielles
        dX_dTeff = 1.0 / (teff * np.log(10))  # dX/dT_eff
        dlogM_dX = MASS_COEFFICIENTS['a2'] + 2 * MASS_COEFFICIENTS['a3'] * X
        dlogM_dlogrho = (MASS_COEFFICIENTS['a4'] +
                        2 * MASS_COEFFICIENTS['a5'] * log_rho +
                        3 * MASS_COEFFICIENTS['a6'] * log_rho2)
        dlogM_dfeh = MASS_COEFFICIENTS['a7']
        
        # Erreur sur X
        dX = teff_err * dX_dTeff
        
        # Erreur sur log M (propagation)
        log_M_err_sq = (
            (dlogM_dX * dX)**2 +
            (dlogM_dlogrho * log_rho_err)**2 +
            (dlogM_dfeh * feh_err)**2 +
            # Ajouter l'erreur intrinsèque de la calibration
            SIGMA_LOG_M**2
        )
        
        log_M_err = np.sqrt(log_M_err_sq)
        
        # Erreur sur M (en M☉)
        # dM/M = ln(10) * d(log M)
        mass_err = mass * np.log(10) * log_M_err
        
        return mass, mass_err
    else:
        return mass


def calculate_stellar_radius(teff, log_rho, feh, teff_err=None, log_rho_err=None, feh_err=None):
    """
    Calcule le rayon stellaire en utilisant la calibration d'Enoch et al. (2010).
    
    Équation: log R = b₁ + b₂X + b₃log ρ + b₄[Fe/H]
    où X = log(T_eff) - 4.1
    
    Parameters
    ----------
    teff : float
        Température effective (en K)
    log_rho : float
        Logarithme de la densité stellaire (log₁₀(ρ) en g/cm³)
    feh : float
        Métallicité [Fe/H] (décimal)
    teff_err : float, optional
        Erreur sur T_eff (en K). Si fourni, calcule l'erreur sur R
    log_rho_err : float, optional
        Erreur sur log ρ. Si fourni, calcule l'erreur sur R
    feh_err : float, optional
        Erreur sur [Fe/H]. Si fourni, calcule l'erreur sur R
    
    Returns
    -------
    radius : float
        Rayon stellaire (en R☉)
    radius_err : float, optional
        Erreur sur le rayon (en R☉). Retourné seulement si toutes les erreurs sont fournies
    """
    # Vérifications
    if teff <= 0:
        raise ValueError(f"T_eff doit être positif, reçu: {teff}")
    if not np.isfinite(log_rho):
        raise ValueError(f"log ρ doit être fini, reçu: {log_rho}")
    
    # Calculer X = log(T_eff) - 4.1
    X = np.log10(teff) - 4.1
    
    # Calculer log R
    log_R = (RADIUS_COEFFICIENTS['b1'] +
             RADIUS_COEFFICIENTS['b2'] * X +
             RADIUS_COEFFICIENTS['b3'] * log_rho +
             RADIUS_COEFFICIENTS['b4'] * feh)
    
    # Convertir en rayon (R☉)
    radius = 10.0**log_R
    
    # Calcul de l'erreur si toutes les erreurs sont fournies
    if teff_err is not None and log_rho_err is not None and feh_err is not None:
        # Propagation d'erreurs
        dX_dTeff = 1.0 / (teff * np.log(10))
        dlogR_dX = RADIUS_COEFFICIENTS['b2']
        dlogR_dlogrho = RADIUS_COEFFICIENTS['b3']
        dlogR_dfeh = RADIUS_COEFFICIENTS['b4']
        
        dX = teff_err * dX_dTeff
        
        # Erreur sur log R
        log_R_err_sq = (
            (dlogR_dX * dX)**2 +
            (dlogR_dlogrho * log_rho_err)**2 +
            (dlogR_dfeh * feh_err)**2 +
            # Ajouter l'erreur intrinsèque de la calibration
            SIGMA_LOG_R**2
        )
        
        log_R_err = np.sqrt(log_R_err_sq)
        
        # Erreur sur R (en R☉)
        radius_err = radius * np.log(10) * log_R_err
        
        return radius, radius_err
    else:
        return radius


def calculate_stellar_mass_and_radius(teff, log_rho, feh, teff_err=None, log_rho_err=None, feh_err=None):
    """
    Calcule simultanément la masse et le rayon stellaires.
    
    Parameters
    ----------
    teff : float
        Température effective (en K)
    log_rho : float
        Logarithme de la densité stellaire (log₁₀(ρ) en g/cm³)
    feh : float
        Métallicité [Fe/H] (décimal)
    teff_err : float, optional
        Erreur sur T_eff (en K)
    log_rho_err : float, optional
        Erreur sur log ρ
    feh_err : float, optional
        Erreur sur [Fe/H]
    
    Returns
    -------
    result : dict
        Dictionnaire contenant:
        - 'mass' : Masse stellaire (M☉)
        - 'radius' : Rayon stellaire (R☉)
        - 'mass_err' : Erreur sur la masse (M☉), si erreurs fournies
        - 'radius_err' : Erreur sur le rayon (R☉), si erreurs fournies
    """
    if teff_err is not None and log_rho_err is not None and feh_err is not None:
        mass, mass_err = calculate_stellar_mass(teff, log_rho, feh, teff_err, log_rho_err, feh_err)
        radius, radius_err = calculate_stellar_radius(teff, log_rho, feh, teff_err, log_rho_err, feh_err)
        return {
            'mass': mass,
            'radius': radius,
            'mass_err': mass_err,
            'radius_err': radius_err
        }
    else:
        mass = calculate_stellar_mass(teff, log_rho, feh)
        radius = calculate_stellar_radius(teff, log_rho, feh)
        return {
            'mass': mass,
            'radius': radius
        }


def rho_to_log_rho(rho):
    """
    Convertit la densité stellaire en logarithme.
    
    Parameters
    ----------
    rho : float
        Densité stellaire (en g/cm³)
    
    Returns
    -------
    log_rho : float
        Logarithme de la densité (log₁₀(ρ))
    """
    if rho <= 0:
        raise ValueError(f"La densité doit être positive, reçu: {rho}")
    return np.log10(rho)


def log_rho_to_rho(log_rho):
    """
    Convertit le logarithme de la densité en densité.
    
    Parameters
    ----------
    log_rho : float
        Logarithme de la densité (log₁₀(ρ))
    
    Returns
    -------
    rho : float
        Densité stellaire (en g/cm³)
    """
    return 10.0**log_rho
