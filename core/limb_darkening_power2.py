"""
Module pour l'implémentation de la loi de limb-darkening "power-2"
basée sur Morello et al. (2017), AJ, 154, 111.

La loi power-2 est définie comme :
I(μ) / I(1) = 1 - c * (1 - μ^α)

où :
- μ = cos(θ) avec θ l'angle entre la normale à la surface et la ligne de visée
- c et α sont les deux coefficients
- I(1) est l'intensité au centre du disque

Cette loi surpasse les autres lois à deux coefficients (quadratique, square-root)
particulièrement pour les étoiles froides.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
import logging

logger = logging.getLogger(__name__)


def power2_intensity(mu, c, alpha):
    """
    Calcule l'intensité normalisée selon la loi power-2.
    
    Parameters
    ----------
    mu : array-like
        Cosinus de l'angle entre la normale et la ligne de visée (0 ≤ μ ≤ 1)
    c : float
        Premier coefficient de la loi power-2
    alpha : float
        Deuxième coefficient (exposant) de la loi power-2
    
    Returns
    -------
    intensity : array-like
        Intensité normalisée I(μ)/I(1)
    """
    mu = np.asarray(mu)
    # Éviter les valeurs négatives ou nulles de mu
    mu = np.clip(mu, 1e-10, 1.0)
    
    # Loi power-2 : I(μ)/I(1) = 1 - c * (1 - μ^α)
    intensity = 1.0 - c * (1.0 - mu**alpha)
    
    # S'assurer que l'intensité reste positive
    intensity = np.maximum(intensity, 0.0)
    
    return intensity


def quadratic_intensity(mu, u1, u2):
    """
    Loi quadratique standard : I(μ)/I(1) = 1 - u₁(1-μ) - u₂(1-μ)².
    μ = cos θ (centre du disque → μ = 1).
    """
    mu = np.asarray(mu)
    mu = np.clip(mu, 1e-10, 1.0)
    one_m = 1.0 - mu
    intensity = 1.0 - u1 * one_m - u2 * (one_m ** 2)
    return np.maximum(intensity, 0.0)


def square_root_intensity(mu, u1, u2):
    """
    Loi « square-root » : I(μ)/I(1) = 1 - u₁(1-μ) - u₂(1-√μ).
    """
    mu = np.asarray(mu)
    mu = np.clip(mu, 1e-10, 1.0)
    intensity = 1.0 - u1 * (1.0 - mu) - u2 * (1.0 - np.sqrt(mu))
    return np.maximum(intensity, 0.0)


def _annulus_grid(n_annuli):
    dr = 1.0 / n_annuli
    r_centers = np.linspace(0.5 * dr, 1.0 - 0.5 * dr, n_annuli)
    mu = np.sqrt(np.clip(1.0 - r_centers ** 2, 1e-10, 1.0))
    return r_centers, dr, mu


def _occultation_fraction_per_annulus(r_centers, d, p):
    """Fraction de l'anneau de rayon r occultée par un disque de rayon p à distance d des centres."""
    r = np.asarray(r_centers, dtype=float)
    d = float(d)
    p = float(p)
    f_occ = np.zeros_like(r)
    m0 = d >= (r + p)
    f_occ[m0] = 0.0
    m1 = (d + r) <= p
    f_occ[m1] = 1.0
    m2 = (d + p) <= r
    f_occ[m2] = 0.0
    m_part = ~(m0 | m1 | m2)
    if np.any(m_part):
        rr = r[m_part]
        cos_theta = (rr ** 2 + d ** 2 - p ** 2) / (2.0 * rr * d)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        f_occ[m_part] = np.arccos(cos_theta) / np.pi
    return f_occ


def compute_occulted_flux_ld(rp_rs, z, intensity_fn, n_annuli=10000):
    """
    Fraction de flux stellaire occultée pour un profil d'intensité I(μ) donné (loi de limbe).

    Parameters
    ----------
    rp_rs : float
    z : array-like — distance projetée centre étoile – planète (rayons stellaires)
    intensity_fn : callable
        mu -> I(μ)/I(1), accepte un ndarray μ
    """
    z = np.asarray(z, dtype=float)
    p = float(rp_rs)
    r_centers, dr, mu = _annulus_grid(n_annuli)
    I = np.maximum(np.asarray(intensity_fn(mu), dtype=float), 0.0)
    flux_per_annulus = 2.0 * np.pi * I * r_centers * dr
    total_flux = float(np.sum(flux_per_annulus))
    if total_flux <= 0:
        return np.zeros_like(z)

    flux_occulted = np.zeros_like(z)
    for i, z_val in enumerate(z):
        if z_val >= 1.0 + p:
            flux_occulted[i] = 0.0
        elif z_val <= p - 1.0:
            flux_occulted[i] = 1.0
        else:
            f_ann = _occultation_fraction_per_annulus(r_centers, z_val, p)
            occ = float(np.sum(flux_per_annulus * f_ann))
            flux_occulted[i] = occ / total_flux
    return flux_occulted


def _z_from_orbit(time, period, t0, a_rs, inclination_deg, e=0.0, w=0.0):
    """Séparation projetée z(t) en rayons stellaires (orbite circulaire)."""
    time = np.asarray(time, dtype=float)
    i_rad = np.radians(inclination_deg)
    phase = 2.0 * np.pi * (time - t0) / period
    z_sq = (a_rs ** 2) * (
        np.sin(phase) ** 2 + np.cos(i_rad) ** 2 * np.cos(phase) ** 2
    )
    return np.sqrt(z_sq)


def transit_lightcurve_quadratic(
    time,
    period,
    t0,
    rp_rs,
    a_rs,
    inclination,
    u1=0.42,
    u2=0.28,
    n_annuli=10000,
):
    """Transit avec loi quadratique (coefficients u₁, u₂ ; défauts type solaire approximatif, bande V)."""
    z = _z_from_orbit(time, period, t0, a_rs, inclination)
    occ = compute_occulted_flux_ld(
        rp_rs, z, lambda mu: quadratic_intensity(mu, u1, u2), n_annuli=n_annuli
    )
    return 1.0 - occ


def transit_lightcurve_square_root(
    time,
    period,
    t0,
    rp_rs,
    a_rs,
    inclination,
    u1=0.44,
    u2=0.26,
    n_annuli=10000,
):
    """Transit avec loi square-root (défauts approximatifs ; même ordre que le quadratique pour comparaison)."""
    z = _z_from_orbit(time, period, t0, a_rs, inclination)
    occ = compute_occulted_flux_ld(
        rp_rs, z, lambda mu: square_root_intensity(mu, u1, u2), n_annuli=n_annuli
    )
    return 1.0 - occ


def calculate_occulted_flux_power2(rp_rs, z, c, alpha, n_annuli=10000):
    """
    Fraction de flux occulté (loi power-2), via intégration sur anneaux (Morello et al. 2017).
    """
    return compute_occulted_flux_ld(
        rp_rs, z, lambda mu: power2_intensity(mu, c, alpha), n_annuli=n_annuli
    )


def transit_lightcurve_power2(time, period, t0, rp_rs, a_rs, inclination, e=0.0, w=0.0, 
                               c=0.5, alpha=0.5, n_annuli=10000):
    """
    Génère une courbe de lumière de transit en utilisant la loi power-2 pour le limb-darkening.
    
    Parameters
    ----------
    time : array-like
        Temps en jours (JD)
    period : float
        Période orbitale en jours
    t0 : float
        Temps de conjonction (mid-transit) en JD
    rp_rs : float
        Rapport des rayons planète/étoile
    a_rs : float
        Demi-grand axe en unités de rayon stellaire
    inclination : float
        Inclinaison orbitale en degrés
    e : float, optional
        Excentricité (défaut: 0.0)
    w : float, optional
        Argument du périastre en degrés (défaut: 0.0)
    c : float, optional
        Premier coefficient de la loi power-2 (défaut: 0.5)
    alpha : float, optional
        Deuxième coefficient de la loi power-2 (défaut: 0.5)
    n_annuli : int, optional
        Nombre d'anneaux pour l'intégration (défaut: 10000)
    
    Returns
    -------
    flux : array-like
        Flux normalisé (1.0 = hors transit)
    """
    time = np.asarray(time)
    z = _z_from_orbit(time, period, t0, a_rs, inclination, e=e, w=w)
    flux_occulted = calculate_occulted_flux_power2(rp_rs, z, c, alpha, n_annuli)
    
    # Flux normalisé : F = 1 - (flux occulté)
    flux = 1.0 - flux_occulted
    
    return flux


def fit_power2_coefficients(time, flux, flux_err, period, t0, rp_rs, a_rs, inclination,
                            c_initial=0.5, alpha_initial=0.5, n_annuli=5000):
    """
    Ajuste les coefficients c et α de la loi power-2 aux données observées.
    
    Parameters
    ----------
    time : array-like
        Temps en jours
    flux : array-like
        Flux observé normalisé
    flux_err : array-like
        Erreurs sur le flux
    period : float
        Période orbitale
    t0 : float
        Temps de conjonction
    rp_rs : float
        Rapport des rayons planète/étoile
    a_rs : float
        Demi-grand axe en unités de rayon stellaire
    inclination : float
        Inclinaison en degrés
    c_initial : float, optional
        Valeur initiale pour c (défaut: 0.5)
    alpha_initial : float, optional
        Valeur initiale pour α (défaut: 0.5)
    n_annuli : int, optional
        Nombre d'anneaux pour l'intégration (défaut: 5000)
    
    Returns
    -------
    result : dict
        Dictionnaire contenant :
        - 'c': coefficient c ajusté
        - 'alpha': coefficient α ajusté
        - 'chi2': chi2 final
        - 'success': booléen indiquant le succès de l'ajustement
    """
    def chi2_function(params):
        c, alpha = params
        
        # Limites physiques
        if c < 0 or c > 1 or alpha < 0.1 or alpha > 10:
            return 1e10
        
        try:
            # Générer le modèle
            model = transit_lightcurve_power2(
                time, period, t0, rp_rs, a_rs, inclination,
                c=c, alpha=alpha, n_annuli=n_annuli
            )
            
            # Calculer le chi2
            chi2 = np.sum(((flux - model) / flux_err)**2)
            return chi2
        except:
            return 1e10
    
    # Ajustement
    initial_params = [c_initial, alpha_initial]
    bounds = [(0.0, 1.0), (0.1, 10.0)]
    
    result_opt = minimize(chi2_function, initial_params, method='L-BFGS-B', bounds=bounds,
                         options={'maxiter': 100, 'ftol': 1e-6})
    
    if result_opt.success:
        return {
            'c': float(result_opt.x[0]),
            'alpha': float(result_opt.x[1]),
            'chi2': float(result_opt.fun),
            'success': True
        }
    else:
        logger.warning(f"Ajustement power-2 échoué: {result_opt.message}")
        return {
            'c': c_initial,
            'alpha': alpha_initial,
            'chi2': 1e10,
            'success': False
        }
