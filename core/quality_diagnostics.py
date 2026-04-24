"""
Module pour les diagnostics de qualité automatiques des courbes de lumière de transit.

Ce module fournit des fonctions pour :
1. Valider les paramètres récupérés depuis des données synthétiques
2. Détecter automatiquement les problèmes de qualité
3. Générer des rapports de diagnostic
"""

import numpy as np
import logging
from typing import Dict, List, Tuple, Optional
from scipy.stats import shapiro, normaltest
from statsmodels.tsa.stattools import acf

logger = logging.getLogger(__name__)


class QualityDiagnostics:
    """
    Classe pour effectuer des diagnostics de qualité automatiques
    sur les courbes de lumière de transit.
    """
    
    def __init__(self):
        self.warnings = []
        self.errors = []
        self.info = []
    
    def validate_transit_depth(self, depth: float, expected_range: Tuple[float, float] = (0.0001, 0.5)) -> bool:
        """
        Valide que la profondeur de transit est dans une plage raisonnable.
        
        Parameters
        ----------
        depth : float
            Profondeur de transit (fraction)
        expected_range : tuple
            Plage attendue (min, max)
        
        Returns
        -------
        bool
            True si valide, False sinon
        """
        if depth < expected_range[0]:
            self.warnings.append(f"Profondeur de transit très faible: {depth:.6f} (< {expected_range[0]})")
            return False
        if depth > expected_range[1]:
            self.errors.append(f"Profondeur de transit invalide: {depth:.6f} (> {expected_range[1]})")
            return False
        return True
    
    def validate_orbital_parameters(self, period: float, a_rs: float, inclination: float) -> bool:
        """
        Valide la cohérence des paramètres orbitaux.
        
        Parameters
        ----------
        period : float
            Période orbitale (jours)
        a_rs : float
            Demi-grand axe en unités de rayon stellaire
        inclination : float
            Inclinaison orbitale (degrés)
        
        Returns
        -------
        bool
            True si valide, False sinon
        """
        valid = True
        
        if period <= 0:
            self.errors.append(f"Période orbitale invalide: {period}")
            valid = False
        
        if a_rs < 1.0:
            self.errors.append(f"a/R* invalide: {a_rs:.2f} (< 1.0, planète à l'intérieur de l'étoile)")
            valid = False
        
        if not (0 < inclination <= 90):
            self.errors.append(f"Inclinaison invalide: {inclination:.2f}° (doit être entre 0° et 90°)")
            valid = False
        
        # Vérifier la cohérence physique
        if a_rs > 1.0 and inclination < 90:
            # Calculer le paramètre d'impact attendu
            b_expected = a_rs * np.cos(np.radians(inclination))
            if b_expected < 0 or b_expected >= 1.0 + 0.1:  # Tolérance pour les erreurs numériques
                self.warnings.append(f"Paramètre d'impact calculé suspect: b = {b_expected:.3f}")
        
        return valid
    
    def check_residuals_quality(self, residuals: np.ndarray, 
                                sigma_oot: float,
                                check_normality: bool = True,
                                check_autocorr: bool = True) -> Dict[str, float]:
        """
        Vérifie la qualité des résidus.
        
        Parameters
        ----------
        residuals : array-like
            Résidus (données - modèle)
        sigma_oot : float
            Écart-type hors transit
        check_normality : bool
            Vérifier la normalité des résidus
        check_autocorr : bool
            Vérifier l'autocorrélation
        
        Returns
        -------
        dict
            Dictionnaire contenant les métriques de qualité
        """
        residuals = np.asarray(residuals)
        metrics = {}
        
        # Statistiques de base
        metrics['mean'] = np.mean(residuals)
        metrics['std'] = np.std(residuals)
        metrics['rms'] = np.sqrt(np.mean(residuals**2))
        
        # Comparaison avec sigma_oot
        if sigma_oot > 0:
            metrics['rms_over_sigma'] = metrics['rms'] / sigma_oot
            if metrics['rms_over_sigma'] > 1.5:
                self.warnings.append(f"RMS des résidus élevé: {metrics['rms']:.6f} ({metrics['rms_over_sigma']:.2f}× sigma_OOT)")
        
        # Test de normalité
        if check_normality and len(residuals) > 3:
            try:
                # Test de Shapiro-Wilk (pour petits échantillons)
                if len(residuals) <= 5000:
                    _, shapiro_p = shapiro(residuals)
                    metrics['shapiro_p'] = shapiro_p
                    if shapiro_p < 0.05:
                        self.warnings.append(
                            f"Résidus non-normaux (Shapiro-Wilk : p={shapiro_p:.4f} ; "
                            "p < 0.05 → le test rejette la normalité, pas une valeur des résidus)"
                        )
                else:
                    # Test de D'Agostino pour grands échantillons
                    _, dagostino_p = normaltest(residuals)
                    metrics['dagostino_p'] = dagostino_p
                    if dagostino_p < 0.05:
                        self.warnings.append(
                            f"Résidus non-normaux (D'Agostino : p={dagostino_p:.4f} ; "
                            "p < 0.05 → le test rejette la normalité)"
                        )
            except Exception as e:
                logger.debug(f"Erreur lors du test de normalité: {e}")
        
        # Autocorrélation
        if check_autocorr and len(residuals) > 10:
            try:
                acf_vals = acf(residuals, nlags=3)
                metrics['acf_1'] = acf_vals[1] if len(acf_vals) > 1 else np.nan
                metrics['acf_2'] = acf_vals[2] if len(acf_vals) > 2 else np.nan
                
                if abs(metrics['acf_1']) > 0.2:
                    self.warnings.append(f"Autocorrélation détectée (ACF lag-1 = {metrics['acf_1']:.3f}, > 0.2)")
            except Exception as e:
                logger.debug(f"Erreur lors du calcul de l'ACF: {e}")
        
        return metrics
    
    def check_chi2_quality(self, chi2: float, n_data: int, n_params: int = 4) -> Dict[str, any]:
        """
        Vérifie la qualité du chi2.
        
        Parameters
        ----------
        chi2 : float
            Chi2 réduit
        n_data : int
            Nombre de points de données
        n_params : int
            Nombre de paramètres ajustés
        
        Returns
        -------
        dict
            Dictionnaire contenant les métriques
        """
        metrics = {'chi2': chi2, 'n_data': n_data, 'n_params': n_params}
        metrics['dof'] = n_data - n_params  # Degrés de liberté
        
        if chi2 < 0.8:
            self.warnings.append(f"Chi2 réduit faible: {chi2:.2f} (< 0.8, erreurs peut-être sur-estimées)")
        elif chi2 > 1.2:
            self.warnings.append(f"Chi2 réduit élevé: {chi2:.2f} (> 1.2, bruit non modélisé ou mauvais ajustement)")
        elif 0.8 <= chi2 <= 1.2:
            self.info.append(f"Chi2 réduit excellent: {chi2:.2f}")
        
        return metrics
    
    def check_limb_darkening_bias(self, depth_blue: float, depth_red: float,
                                   depth_blue_err: Optional[float] = None,
                                   depth_red_err: Optional[float] = None,
                                   threshold_ppm: float = 30.0) -> Dict[str, any]:
        """
        Vérifie les biais potentiels dus au limb-darkening en comparant
        les profondeurs entre filtres bleu et rouge.
        
        Parameters
        ----------
        depth_blue : float
            Profondeur de transit dans le bleu
        depth_red : float
            Profondeur de transit dans le rouge
        depth_blue_err : float, optional
            Erreur sur la profondeur bleue
        depth_red_err : float, optional
            Erreur sur la profondeur rouge
        threshold_ppm : float
            Seuil en ppm pour détecter un biais (défaut: 30 ppm)
        
        Returns
        -------
        dict
            Dictionnaire contenant les métriques de comparaison
        """
        delta_depth = depth_blue - depth_red
        delta_depth_ppm = delta_depth * 1e6  # Convertir en ppm
        
        metrics = {
            'depth_blue': depth_blue,
            'depth_red': depth_red,
            'delta_depth': delta_depth,
            'delta_depth_ppm': delta_depth_ppm
        }
        
        # Calculer l'erreur sur la différence si les erreurs sont fournies
        if depth_blue_err is not None and depth_red_err is not None:
            delta_err = np.sqrt(depth_blue_err**2 + depth_red_err**2)
            delta_err_ppm = delta_err * 1e6
            metrics['delta_err'] = delta_err
            metrics['delta_err_ppm'] = delta_err_ppm
            metrics['significance'] = abs(delta_depth) / delta_err if delta_err > 0 else np.inf
        else:
            metrics['delta_err'] = None
            metrics['delta_err_ppm'] = None
            metrics['significance'] = None
        
        # Détecter les biais
        if abs(delta_depth_ppm) > threshold_ppm:
            if metrics['significance'] is not None and metrics['significance'] > 3:
                self.warnings.append(
                    f"Biais potentiel détecté entre filtres: Δδ = {delta_depth_ppm:.1f} ppm "
                    f"(seuil: {threshold_ppm} ppm, significativité: {metrics['significance']:.1f}σ)"
                )
            else:
                self.warnings.append(
                    f"Différence entre filtres: Δδ = {delta_depth_ppm:.1f} ppm "
                    f"(seuil: {threshold_ppm} ppm, peut être dû au limb-darkening)"
                )
        else:
            self.info.append(
                f"Profondeurs cohérentes entre filtres: Δδ = {delta_depth_ppm:.1f} ppm "
                f"(< {threshold_ppm} ppm)"
            )
        
        return metrics
    
    def validate_parameter_recovery(self, recovered: Dict[str, float],
                                    expected: Dict[str, float],
                                    tolerances: Optional[Dict[str, float]] = None) -> Dict[str, bool]:
        """
        Valide la récupération de paramètres depuis des données synthétiques.
        
        Parameters
        ----------
        recovered : dict
            Paramètres récupérés
        expected : dict
            Paramètres attendus (vrais)
        tolerances : dict, optional
            Tolérances pour chaque paramètre (fraction relative)
        
        Returns
        -------
        dict
            Dictionnaire indiquant si chaque paramètre est dans la tolérance
        """
        if tolerances is None:
            # Tolérances par défaut (en fraction relative)
            tolerances = {
                'delta_F': 0.20,  # 20% pour la profondeur (bruit dans les données)
                'b': 0.50,        # 50% pour b (approximation)
                'a_over_R_star': 0.10,  # 10% pour a/R*
                'rho_star': 0.15,  # 15% pour la densité
                'R_planet': 0.20,  # 20% pour le rayon planétaire
                'inclination': 0.10,  # 10% pour l'inclinaison
            }
        
        results = {}
        
        for param_name in expected:
            if param_name not in recovered:
                self.warnings.append(f"Paramètre {param_name} non récupéré")
                results[param_name] = False
                continue
            
            expected_val = expected[param_name]
            recovered_val = recovered[param_name]
            tolerance = tolerances.get(param_name, 0.20)  # 20% par défaut
            
            if expected_val == 0:
                # Pour les valeurs nulles, vérifier que la valeur récupérée est proche de zéro
                if abs(recovered_val) < tolerance:
                    results[param_name] = True
                else:
                    self.warnings.append(
                        f"Paramètre {param_name}: récupéré={recovered_val:.6f}, "
                        f"attendu≈0, écart={abs(recovered_val):.6f}"
                    )
                    results[param_name] = False
            else:
                # Erreur relative
                rel_error = abs(recovered_val - expected_val) / abs(expected_val)
                
                if rel_error <= tolerance:
                    results[param_name] = True
                    self.info.append(
                        f"Paramètre {param_name}: récupéré={recovered_val:.6f}, "
                        f"attendu={expected_val:.6f}, erreur={rel_error*100:.1f}%"
                    )
                else:
                    results[param_name] = False
                    self.warnings.append(
                        f"Paramètre {param_name}: récupéré={recovered_val:.6f}, "
                        f"attendu={expected_val:.6f}, erreur={rel_error*100:.1f}% "
                        f"(tolérance: {tolerance*100:.0f}%)"
                    )
        
        return results
    
    def generate_report(self) -> str:
        """
        Génère un rapport de diagnostic.
        
        Returns
        -------
        str
            Rapport formaté
        """
        lines = ["=" * 60]
        lines.append("RAPPORT DE DIAGNOSTIC DE QUALITÉ")
        lines.append("=" * 60)
        lines.append("(Pour les tests de normalité : p < 0.05 = le test rejette la normalité des résidus.)")
        lines.append("")
        
        if self.info:
            lines.append("\n✓ INFORMATIONS:")
            for msg in self.info:
                lines.append(f"  • {msg}")
        
        if self.warnings:
            lines.append("\n⚠ AVERTISSEMENTS:")
            for msg in self.warnings:
                lines.append(f"  • {msg}")
        
        if self.errors:
            lines.append("\n❌ ERREURS:")
            for msg in self.errors:
                lines.append(f"  • {msg}")
        
        if not self.info and not self.warnings and not self.errors:
            lines.append("\n✓ Aucun problème détecté")
        
        lines.append("\n" + "=" * 60)
        
        return "\n".join(lines)
    
    def clear(self):
        """Réinitialise les diagnostics."""
        self.warnings.clear()
        self.errors.clear()
        self.info.clear()


def test_parameter_recovery_with_synthetic_data(
    solver_function,
    time: np.ndarray,
    true_params: Dict[str, float],
    noise_level: float = 0.001,
    n_trials: int = 10
) -> Dict[str, any]:
    """
    Teste la récupération de paramètres depuis des données synthétiques.
    
    Parameters
    ----------
    solver_function : callable
        Fonction qui résout les paramètres (ex: solve_transit_parameters)
    time : array-like
        Temps (jours)
    true_params : dict
        Paramètres vrais utilisés pour générer les données
    noise_level : float
        Niveau de bruit à ajouter
    n_trials : int
        Nombre d'essais pour statistiques
    
    Returns
    -------
    dict
        Dictionnaire contenant les statistiques de récupération
    """
    recovered_params_list = []
    diagnostics = QualityDiagnostics()
    
    for trial in range(n_trials):
        try:
            # Générer des données synthétiques avec bruit
            # (Cette fonction devrait être fournie par l'utilisateur)
            # flux = generate_synthetic_lightcurve(time, true_params, noise_level)
            # recovered = solver_function(time, flux, ...)
            # recovered_params_list.append(recovered)
            pass
        except Exception as e:
            logger.warning(f"Essai {trial+1} échoué: {e}")
    
    if not recovered_params_list:
        return {'success': False, 'message': 'Aucun essai réussi'}
    
    # Calculer les statistiques
    stats = {}
    for param_name in true_params:
        recovered_vals = [p.get(param_name, np.nan) for p in recovered_params_list]
        recovered_vals = [v for v in recovered_vals if not np.isnan(v)]
        
        if recovered_vals:
            stats[param_name] = {
                'mean': np.mean(recovered_vals),
                'std': np.std(recovered_vals),
                'median': np.median(recovered_vals),
                'true_value': true_params[param_name],
                'bias': np.mean(recovered_vals) - true_params[param_name],
                'rmse': np.sqrt(np.mean((np.array(recovered_vals) - true_params[param_name])**2))
            }
    
    return {
        'success': True,
        'n_trials': n_trials,
        'n_successful': len(recovered_params_list),
        'statistics': stats,
        'diagnostics': diagnostics
    }
