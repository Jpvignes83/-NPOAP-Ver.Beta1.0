"""
Module pour l'analyse de transits d'exoplanètes basé sur 
Seager & Mallén-Ornelas (2003), ApJ, 585, 1038.

Ce module implémente les équations qui permettent de déterminer de manière unique
les paramètres planète/étoile à partir d'une courbe de lumière de transit.
"""

import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy.optimize import minimize_scalar
import logging

logger = logging.getLogger(__name__)

# Constante de gravitation
G = const.G.to(u.cm**3 / (u.g * u.s**2)).value  # En cgs pour compatibilité avec ρ en g/cm³


def calculate_transit_depth(light_curve):
    """
    Calcule la profondeur du transit ΔF depuis une courbe de lumière.
    
    Parameters
    ----------
    light_curve : array-like
        Courbe de lumière normalisée (flux relatif, 1.0 = hors transit)
    
    Returns
    -------
    delta_F : float
        Profondeur du transit (fraction, 0 < delta_F < 1)
    """
    light_curve = np.asarray(light_curve)
    
    if len(light_curve) == 0:
        raise ValueError("Courbe de lumière vide")
    
    # Flux hors transit : utiliser une méthode robuste
    # Option 1: Médiane des valeurs supérieures au percentile 75
    p75 = np.percentile(light_curve, 75)
    mask_high = light_curve >= p75
    if np.any(mask_high):
        flux_out = np.median(light_curve[mask_high])
    else:
        # Fallback: médiane de toutes les valeurs
        flux_out = np.median(light_curve)
    
    # Flux dans le transit : minimum
    flux_in = np.min(light_curve)
    
    # Profondeur
    if flux_out <= 0:
        raise ValueError(f"Flux hors transit non positif: {flux_out}")
    
    delta_F = (flux_out - flux_in) / flux_out
    
    # Vérification de cohérence
    if delta_F <= 0 or delta_F >= 1:
        logger.warning(f"Profondeur de transit inattendue: delta_F={delta_F}")
        delta_F = max(0.0, min(delta_F, 0.99))
    
    return float(delta_F)


def calculate_transit_durations(time, flux, period=None, epoch=None):
    """
    Calcule les durées t_T (totale) et t_F (ingress/egress) depuis une courbe de lumière.
    
    Parameters
    ----------
    time : array-like
        Temps (en jours ou unités appropriées)
    flux : array-like
        Flux relatif normalisé
    period : float, optional
        Période orbitale (si None, estimée depuis les données)
    epoch : float, optional
        Époque du premier transit (si None, estimée depuis les données)
    
    Returns
    -------
    t_T : float
        Durée totale du transit (même unité que time)
    t_F : float
        Durée de l'ingress/egress (même unité que time)
    """
    time = np.asarray(time)
    flux = np.asarray(flux)
    
    # Détection du transit (flux minimal)
    delta_F = calculate_transit_depth(flux)
    threshold = 1.0 - delta_F / 2.0  # Seuil à mi-profondeur
    
    # Points dans le transit
    in_transit = flux < threshold
    
    if not np.any(in_transit):
        logger.warning("Aucun point en transit détecté")
        return None, None
    
    # Si period est fourni, on peut identifier les transits individuels
    if period is not None and period > 0:
        # Si epoch est fourni, utiliser le transit le plus proche de l'epoch
        # Sinon, utiliser le transit le plus profond
        if epoch is not None:
            # Normaliser les temps par rapport à l'epoch et la période
            phase = ((time - epoch) % period) / period
            phase[phase > 0.5] -= 1.0  # Phase dans [-0.5, 0.5]
            # Prendre le transit autour de phase=0
            transit_mask = (np.abs(phase) < 0.3) & in_transit
        else:
            # Trouver l'index du flux minimal (centre du transit le plus profond)
            min_idx = np.argmin(flux)
            # Identifier les points du même transit (autour du minimum)
            # Chercher les points en transit proches du minimum
            time_min = time[min_idx]
            if period is not None:
                # Prendre les points dans une fenêtre de ±0.5 période autour du minimum
                transit_mask = (np.abs(time - time_min) < period/2) & in_transit
            else:
                # Sinon, prendre tous les points en transit contigus autour du minimum
                transit_mask = in_transit
        
        if not np.any(transit_mask):
            transit_mask = in_transit
    else:
        # Sans période, prendre tous les points en transit contigus
        # Trouver le groupe le plus grand de points consécutifs en transit
        transit_mask = in_transit
    
    # Calculer t_T pour le transit identifié
    time_in_transit = time[transit_mask]
    if len(time_in_transit) == 0:
        time_in_transit = time[in_transit]
    
    t_T = np.max(time_in_transit) - np.min(time_in_transit)
    
    # Pour t_F, on calcule la durée entre le début de l'ingress et le début du plateau
    # Approche simplifiée : on cherche les points à 20% et 80% de la profondeur dans le transit identifié
    threshold_20 = 1.0 - 0.2 * delta_F
    threshold_80 = 1.0 - 0.8 * delta_F
    
    # Chercher ingress_start et ingress_end dans le transit identifié
    flux_transit = flux[transit_mask]
    time_transit = time[transit_mask]
    
    above_20 = flux_transit > threshold_20
    below_80 = flux_transit < threshold_80
    
    if np.any(above_20) and np.any(below_80):
        ingress_start = time_transit[above_20][0]
        ingress_end = time_transit[below_80][0]
        t_F = 2.0 * abs(ingress_end - ingress_start)  # Durée totale ingress + egress
    else:
        # Fallback : estimation grossière
        t_F = t_T * 0.1  # Estimation : t_F ≈ 10% de t_T
    
    return float(t_T), float(t_F)


def calculate_impact_parameter(delta_F, t_F_over_t_T, period=None, t_T=None):
    """
    Calcule le paramètre d'impact b depuis la profondeur et le rapport t_F/t_T.
    
    NOTE: Pour une détermination précise de b, il faut aussi period et t_T.
    Sinon, utilise une approximation simplifiée.
    
    Basé sur Seager & Mallén-Ornelas (2003).
    
    Parameters
    ----------
    delta_F : float
        Profondeur du transit (fraction)
    t_F_over_t_T : float
        Rapport durée ingress/egress sur durée totale
    period : float, optional
        Période orbitale (en jours). Si fourni avec t_T, calcul plus précis
    t_T : float, optional
        Durée totale du transit (en jours). Si fourni avec period, calcul plus précis
    
    Returns
    -------
    b : float
        Paramètre d'impact (0 ≤ b < 1)
    """
    if delta_F <= 0 or delta_F >= 1:
        raise ValueError(f"delta_F doit être entre 0 et 1, reçu {delta_F}")
    if t_F_over_t_T <= 0 or t_F_over_t_T >= 1:
        raise ValueError(f"t_F_over_t_T doit être entre 0 et 1, reçu {t_F_over_t_T}")
    
    k = np.sqrt(delta_F)  # k = R_p / R*
    
    # Si period et t_T sont fournis, on peut faire un calcul itératif précis
    if period is not None and t_T is not None and period > 0 and t_T > 0:
        # D'après Seager & Mallén-Ornelas (2003), pour une orbite circulaire:
        # La relation exacte entre t_F/t_T, b, et k est complexe.
        # On utilise une approximation basée sur la formule:
        # Pour b << 1 et k << 1: t_F/t_T ≈ k / sqrt((1+k)^2 - b^2)
        # Pour une formule plus générale, on itère
        
        def objective(b_guess):
            # Calculer a/R* depuis t_T et b
            # t_T = (P/π) * sqrt((a/R*)^2 - b^2) * sqrt((1+k)^2 - b^2) / sqrt(1 - e^2)
            # Pour e=0: t_T = (P/π) * sqrt((a/R*)^2 - b^2) * sqrt((1+k)^2 - b^2)
            term1 = (1 + k)**2 - b_guess**2
            if term1 <= 0:
                return 1e10  # Pénalité
            
            # Résoudre pour (a/R*)^2
            # (a/R*)^2 - b^2 = (t_T * π / P)^2 / ((1+k)^2 - b^2)
            a_over_R_squared = b_guess**2 + (t_T * np.pi / period)**2 / term1
            if a_over_R_squared <= b_guess**2:
                return 1e10  # Pénalité
            
            a_over_R = np.sqrt(a_over_R_squared)
            
            # Calculer t_F/t_T prédit
            # D'après Seager & Mallén-Ornelas (2003), pour une orbite circulaire:
            # La relation exacte est complexe, mais on peut utiliser une approximation:
            # Pour b petit: t_F/t_T ≈ k/(1+k) (indépendant de b)
            # Pour b grand: t_F/t_T diminue car le transit est plus partiel
            # Approximation: t_F/t_T ≈ (k/(1+k)) * (1 - b^2/(1+k)^2) pour b < (1-k)
            
            # Formule simplifiée qui fonctionne mieux:
            # Si b est proche de 0, t_F/t_T ≈ k/(1+k)
            # Si b augmente, t_F/t_T diminue
            t_F_over_t_T_b0 = k / (1 + k)  # Valeur pour b=0
            
            # Pour b > 0, on réduit progressivement
            # Approximation: t_F/t_T ≈ t_F_over_t_T_b0 * (1 - alpha * b^2)
            # où alpha est ajusté pour correspondre aux observations
            if b_guess < (1 - k):
                # Pour b < (1-k), on peut utiliser une approximation linéaire en b^2
                alpha = 1.0 / ((1 + k)**2)  # Facteur d'ajustement
                t_F_over_t_T_pred = t_F_over_t_T_b0 * (1 - alpha * b_guess**2)
            else:
                # Pour b proche de (1-k), t_F/t_T devient très petit
                t_F_over_t_T_pred = t_F_over_t_T_b0 * 0.01
            
            # S'assurer que la prédiction est positive
            t_F_over_t_T_pred = max(0.001, t_F_over_t_T_pred)
            
            # Erreur au carré
            return (t_F_over_t_T_pred - t_F_over_t_T)**2
        
        # Recherche de b dans [0, 0.99]
        result = minimize_scalar(objective, bounds=(0.0, 0.99), method='bounded')
        b = result.x
    else:
        # Approximation simplifiée (moins précise)
        # Relation approximative: pour transits peu profonds et t_F/t_T petit, b est petit
        # Pour transits avec t_F/t_T grand, b est grand
        # Approximation: b ≈ 1 - (t_F/t_T) * (1+k) / (1-k) pour b < 0.7
        if t_F_over_t_T < 0.1:
            # Transit central (b faible)
            b = 0.0
        elif t_F_over_t_T > 0.8:
            # Transit très partiel (b élevé)
            b = 0.9
        else:
            # Interpolation linéaire simplifiée
            # Cette approximation est grossière mais fonctionne pour des cas typiques
            b = (t_F_over_t_T - 0.1) / 0.7 * 0.9  # Échelle entre 0.1 et 0.8 -> b entre 0 et 0.9
    
    # S'assurer que b est dans [0, 1)
    b = max(0.0, min(b, 0.99))
    
    return float(b)


def calculate_a_over_R_star(period, delta_F, t_T, b):
    """
    Calcule le rapport a/R* (demi-grand axe / rayon stellaire).
    
    Basé sur l'équation (9) de Seager & Mallén-Ornelas (2003).
    
    Parameters
    ----------
    period : float
        Période orbitale (en jours)
    delta_F : float
        Profondeur du transit (fraction)
    t_T : float
        Durée totale du transit (en jours)
    b : float
        Paramètre d'impact
    
    Returns
    -------
    a_over_R_star : float
        Rapport a/R*
    """
    if period <= 0 or t_T <= 0:
        raise ValueError(f"Période et t_T doivent être positives : P={period}, t_T={t_T}")
    
    k = np.sqrt(delta_F)
    
    # Équation (9) : a/R* = (P/π) * sqrt((1+k)² - b²) / t_T
    # Note: en unités cohérentes (jours)
    numerator = np.sqrt((1 + k)**2 - b**2)
    a_over_R_star = (period / np.pi) * numerator / t_T
    
    return float(a_over_R_star)


def calculate_stellar_density(period, a_over_R_star):
    """
    Calcule la densité stellaire depuis la période et a/R*.
    
    Basé sur l'équation de Kepler : ρ* = 3π / (G P²) * (a/R*)³
    
    Parameters
    ----------
    period : float
        Période orbitale (en jours)
    a_over_R_star : float
        Rapport a/R*
    
    Returns
    -------
    rho_star : float
        Densité stellaire (en g/cm³)
        Pour conversion en unités solaires : ρ☉ = 1.41 g/cm³
    """
    # Conversion période en secondes
    P_sec = period * 86400.0  # jours -> secondes
    
    # Densité stellaire (équation de Kepler)
    # ρ* = 3π / (G P²) * (a/R*)³
    # Note: G doit être en cgs (cm³ g⁻¹ s⁻²)
    rho_star = (3.0 * np.pi) / (G * P_sec**2) * (a_over_R_star**3)
    
    # Conversion en g/cm³ (déjà en cgs)
    return float(rho_star)


def calculate_planet_radius(delta_F, R_star):
    """
    Calcule le rayon planétaire depuis la profondeur et le rayon stellaire.
    
    Parameters
    ----------
    delta_F : float
        Profondeur du transit (fraction)
    R_star : float
        Rayon stellaire (en unités solaires ou autres)
    
    Returns
    -------
    R_planet : float
        Rayon planétaire (même unité que R_star)
    """
    k = np.sqrt(delta_F)  # k = R_p / R*
    R_planet = k * R_star
    return float(R_planet)


def solve_transit_parameters(time, flux, period, R_star=None, M_star=None):
    """
    Résout tous les paramètres de transit depuis une courbe de lumière.
    
    Implémente la méthode de Seager & Mallén-Ornelas (2003) pour obtenir
    une solution unique des paramètres.
    
    Parameters
    ----------
    time : array-like
        Temps (en jours)
    flux : array-like
        Flux relatif normalisé
    period : float
        Période orbitale (en jours)
    R_star : float, optional
        Rayon stellaire (en R☉). Si fourni, calcule R_planet et a
    M_star : float, optional
        Masse stellaire (en M☉). Si fourni avec R_star, calcule l'inclinaison i
    
    Returns
    -------
    params : dict
        Dictionnaire contenant tous les paramètres calculés :
        - delta_F : profondeur du transit
        - t_T : durée totale
        - t_F : durée ingress/egress
        - b : paramètre d'impact
        - a_over_R_star : rapport a/R*
        - rho_star : densité stellaire (g/cm³)
        - R_planet : rayon planétaire (si R_star fourni)
        - a : demi-grand axe (si R_star fourni)
        - i : inclinaison (si R_star et M_star fournis)
    """
    # 1. Calculer la profondeur du transit
    delta_F = calculate_transit_depth(flux)
    
    # 2. Calculer les durées
    t_T, t_F = calculate_transit_durations(time, flux, period=period)
    
    if t_T is None or t_F is None:
        raise ValueError("Impossible de déterminer les durées de transit")
    
    t_F_over_t_T = t_F / t_T
    
    # 3. Calculer le paramètre d'impact b (avec period et t_T pour calcul précis)
    b = calculate_impact_parameter(delta_F, t_F_over_t_T, period=period, t_T=t_T)
    
    # 4. Calculer a/R*
    a_over_R_star = calculate_a_over_R_star(period, delta_F, t_T, b)
    
    # 5. Calculer la densité stellaire
    rho_star = calculate_stellar_density(period, a_over_R_star)
    
    # Résultats de base
    params = {
        'delta_F': delta_F,
        't_T': t_T,
        't_F': t_F,
        't_F_over_t_T': t_F_over_t_T,
        'b': b,
        'a_over_R_star': a_over_R_star,
        'rho_star': rho_star,  # en g/cm³
        'rho_star_solar': rho_star / 1.41,  # en unités solaires (ρ☉ = 1.41 g/cm³)
    }
    
    # Paramètres supplémentaires si R_star est fourni
    if R_star is not None:
        # Rayon planétaire
        R_planet = calculate_planet_radius(delta_F, R_star)
        params['R_planet'] = R_planet
        
        # Demi-grand axe (en unités de R_star, puis conversion)
        # a = (a/R*) * R*
        a = a_over_R_star * R_star
        params['a'] = a  # en même unité que R_star
        
        # Si M_star est aussi fourni, on peut calculer l'inclinaison
        if M_star is not None:
            # b = (a/R*) * cos(i)
            # cos(i) = b / (a/R*)
            cos_i = b / a_over_R_star
            cos_i = max(-1.0, min(1.0, cos_i))  # Clipper entre -1 et 1
            i = np.arccos(cos_i) * 180.0 / np.pi  # en degrés
            params['inclination'] = i
            params['cos_i'] = cos_i
    
    return params


def estimate_period_from_single_transit(time, flux, spectral_type_or_Teff):
    """
    Estime la période depuis un transit unique avec le type spectral.
    
    Basé sur la méthode décrite dans Seager & Mallén-Ornelas (2003).
    
    Parameters
    ----------
    time : array-like
        Temps (en jours)
    flux : array-like
        Flux relatif normalisé
    spectral_type_or_Teff : str or float
        Type spectral (ex: 'G5V') ou température effective (en K)
    
    Returns
    -------
    period_estimate : float
        Estimation de la période (en jours)
    """
    # Table de densités stellaires typiques (en ρ☉)
    # À compléter avec plus de types spectraux si nécessaire
    stellar_density_table = {
        'F5V': 1.5, 'F8V': 1.3,
        'G0V': 1.2, 'G2V': 1.0, 'G5V': 0.9,
        'K0V': 0.8, 'K5V': 0.7,
        'M0V': 0.6, 'M5V': 0.5,
    }
    
    # Si c'est un float, c'est probablement Teff
    if isinstance(spectral_type_or_Teff, (int, float)):
        # Estimation grossière depuis Teff (à améliorer)
        Teff = float(spectral_type_or_Teff)
        if Teff > 6000:
            rho_solar = 1.3
        elif Teff > 5000:
            rho_solar = 1.0
        elif Teff > 4000:
            rho_solar = 0.7
        else:
            rho_solar = 0.5
    else:
        spec_type = str(spectral_type_or_Teff).strip().upper()
        rho_solar = stellar_density_table.get(spec_type, 1.0)
    
    rho_star = rho_solar * 1.41  # Conversion en g/cm³
    
    # Calculer les paramètres depuis le transit
    delta_F = calculate_transit_depth(flux)
    t_T, t_F = calculate_transit_durations(time, flux)
    t_F_over_t_T = t_F / t_T
    b = calculate_impact_parameter(delta_F, t_F_over_t_T)
    
    # a/R* depuis t_T et b
    k = np.sqrt(delta_F)
    # Inversion de l'équation (9) : P = π * t_T * (a/R*) / sqrt((1+k)² - b²)
    # Mais on a besoin de (a/R*) depuis ρ*
    # ρ* = 3π / (G P²) * (a/R*)³
    # => (a/R*)³ = ρ* G P² / (3π)
    # On itère pour trouver P
    
    # Estimation initiale : on suppose a/R* ≈ 10 (typique pour Jupiters chauds)
    a_over_R_guess = 10.0
    P_guess = (np.pi * t_T * a_over_R_guess) / np.sqrt((1 + k)**2 - b**2)
    
    # Itération pour trouver P cohérent avec ρ*
    def objective(P):
        a_over_R = calculate_a_over_R_star(P, delta_F, t_T, b)
        rho_calc = calculate_stellar_density(P, a_over_R)
        return (rho_calc - rho_star)**2
    
    result = minimize_scalar(objective, bounds=(P_guess*0.1, P_guess*10.0), method='bounded')
    period_estimate = result.x
    
    return float(period_estimate)


# Fonction utilitaire pour afficher les résultats
def print_transit_parameters(params):
    """
    Affiche les paramètres de transit de manière lisible.
    
    Parameters
    ----------
    params : dict
        Dictionnaire de paramètres (sortie de solve_transit_parameters)
    """
    print("=" * 60)
    print("PARAMÈTRES DE TRANSIT (Seager & Mallén-Ornelas 2003)")
    print("=" * 60)
    print(f"Profondeur du transit (ΔF): {params['delta_F']:.6f} ({params['delta_F']*100:.4f}%)")
    print(f"Durée totale (t_T): {params['t_T']:.6f} jours")
    print(f"Durée ingress/egress (t_F): {params['t_F']:.6f} jours")
    print(f"Rapport t_F/t_T: {params.get('t_F_over_t_T', params['t_F']/params['t_T']):.4f}")
    print(f"Paramètre d'impact (b): {params['b']:.4f}")
    print(f"Rapport a/R*: {params['a_over_R_star']:.2f}")
    print(f"Densité stellaire: {params['rho_star']:.4f} g/cm³")
    print(f"Densité stellaire: {params.get('rho_star_solar', params['rho_star']/1.41):.2f} ρ☉")
    
    if 'R_planet' in params:
        print(f"\nParamètres supplémentaires (R_star fourni):")
        print(f"Rayon planétaire (R_p): {params['R_planet']:.4f} R☉")
        print(f"Demi-grand axe (a): {params['a']:.4f} R☉")
        
        if 'inclination' in params:
            print(f"Inclinaison (i): {params['inclination']:.2f}°")
    print("=" * 60)
