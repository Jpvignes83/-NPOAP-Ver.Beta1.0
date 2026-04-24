"""
Module LCA (Light Curve Analysis)
Basé sur Dai et al. (2023), Research in Astronomy and Astrophysics, 23:055011

Ce module analyse les courbes de lumière extraites par SPDE pour identifier
les variations potentielles et classifier les étoiles variables.
"""

import numpy as np
import logging
from typing import Dict, List, Tuple, Optional
from astropy.table import Table
from scipy import stats
from scipy.signal import find_peaks
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

logger = logging.getLogger(__name__)


class LCAAnalysis:
    """
    Analyseur de courbes de lumière pour identifier les variations potentielles.
    """
    
    def __init__(self, 
                 st_threshold: float = 0.01,
                 min_amplitude: float = 0.005,
                 min_period: float = 0.1,
                 max_period: float = 100.0):
        """
        Parameters
        ----------
        st_threshold : float
            Seuil de variation (St) pour séparer les variables des constantes
        min_amplitude : float
            Amplitude minimale pour considérer une variation significative
        min_period : float
            Période minimale à rechercher (jours)
        max_period : float
            Période maximale à rechercher (jours)
        """
        self.st_threshold = st_threshold
        self.min_amplitude = min_amplitude
        self.min_period = min_period
        self.max_period = max_period
        
        self.variable_candidates = []
        self.periodic_variables = []
        self.transient_variables = []
        self.peculiar_variables = []
        
    def calculate_variability_index(self, flux: np.ndarray) -> float:
        """
        Calcule l'index de variabilité St (Dai et al. 2023).
        
        St = std(flux) / mean(flux)
        
        Parameters
        ----------
        flux : np.ndarray
            Courbe de lumière (flux)
            
        Returns
        -------
        float
            Index de variabilité St
        """
        if len(flux) == 0 or np.all(flux <= 0):
            return 0.0
        
        flux_positive = flux[flux > 0]
        if len(flux_positive) == 0:
            return 0.0
        
        mean_flux = np.mean(flux_positive)
        std_flux = np.std(flux_positive)
        
        if mean_flux == 0:
            return 0.0
        
        st = std_flux / mean_flux
        return st
    
    def detect_variability(self, light_curves: Dict[str, Dict]) -> List[Dict]:
        """
        Détecte les courbes de lumière variables.
        
        Parameters
        ----------
        light_curves : Dict[str, Dict]
            Dictionnaire des courbes de lumière {star_id: {'flux': ..., 'time': ...}}
            
        Returns
        -------
        List[Dict]
            Liste des étoiles variables détectées
        """
        logger.info(f"Analyse de {len(light_curves)} courbes de lumière")
        
        variable_candidates = []
        
        for star_id, lc_data in light_curves.items():
            flux = lc_data['flux']
            
            # Calculer l'index de variabilité
            st = self.calculate_variability_index(flux)
            
            # Calculer l'amplitude
            if len(flux) > 0:
                amplitude = (np.max(flux) - np.min(flux)) / np.mean(flux)
            else:
                amplitude = 0.0
            
            # Vérifier si variable
            if st >= self.st_threshold and amplitude >= self.min_amplitude:
                variable_candidates.append({
                    'star_id': star_id,
                    'st': st,
                    'amplitude': amplitude,
                    'mean_flux': np.mean(flux),
                    'std_flux': np.std(flux),
                    'n_points': len(flux)
                })
        
        self.variable_candidates = variable_candidates
        logger.info(f"{len(variable_candidates)} étoiles variables détectées")
        
        return variable_candidates
    
    def classify_variables(self, light_curves: Dict[str, Dict],
                          variable_candidates: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Classe les étoiles variables en périodiques, transitoires et particulières.
        
        Parameters
        ----------
        light_curves : Dict[str, Dict]
            Dictionnaire des courbes de lumière
        variable_candidates : List[Dict]
            Liste des candidats variables
            
        Returns
        -------
        Dict[str, List[Dict]]
            Dictionnaire avec les classifications
        """
        logger.info("Classification des étoiles variables")
        
        periodic = []
        transient = []
        peculiar = []
        
        for candidate in variable_candidates:
            star_id = candidate['star_id']
            lc_data = light_curves[star_id]
            flux = lc_data['flux']
            time = lc_data['time']
            
            # Recherche de périodicité
            period = self.find_period(time, flux)
            
            if period is not None:
                # Variable périodique
                candidate['period'] = period
                candidate['type'] = 'periodic'
                periodic.append(candidate)
            else:
                # Vérifier si transitoire (variation brusque)
                is_transient = self.is_transient(time, flux)
                
                if is_transient:
                    candidate['type'] = 'transient'
                    transient.append(candidate)
                else:
                    candidate['type'] = 'peculiar'
                    peculiar.append(candidate)
        
        self.periodic_variables = periodic
        self.transient_variables = transient
        self.peculiar_variables = peculiar
        
        logger.info(f"Classification terminée : {len(periodic)} périodiques, "
                   f"{len(transient)} transitoires, {len(peculiar)} particulières")
        
        return {
            'periodic': periodic,
            'transient': transient,
            'peculiar': peculiar
        }
    
    def find_period(self, time: np.ndarray, flux: np.ndarray,
                   method: str = 'lomb_scargle') -> Optional[float]:
        """
        Trouve la période d'une courbe de lumière.
        
        Parameters
        ----------
        time : np.ndarray
            Temps d'observation
        flux : np.ndarray
            Flux observé
        method : str
            Méthode de recherche ('lomb_scargle' ou 'autocorr')
            
        Returns
        -------
        float or None
            Période trouvée (jours) ou None
        """
        if len(time) < 10:
            return None
        
        # Normaliser le temps
        time_norm = time - time[0]
        time_span = time_norm[-1] - time_norm[0]
        
        if time_span < self.min_period:
            return None
        
        try:
            if method == 'lomb_scargle':
                from astropy.timeseries import LombScargle
                
                # Préparer les données
                flux_norm = (flux - np.mean(flux)) / np.std(flux)
                
                # Périodes à tester
                periods = np.logspace(np.log10(self.min_period), 
                                    np.log10(min(self.max_period, time_span/2)),
                                    1000)
                frequencies = 1.0 / periods
                
                # Calculer le périodogramme de Lomb-Scargle
                ls = LombScargle(time_norm, flux_norm)
                power = ls.power(frequencies)
                
                # Trouver le pic principal
                peaks, properties = find_peaks(power, height=0.3, distance=10)
                
                if len(peaks) > 0:
                    best_peak_idx = peaks[np.argmax(power[peaks])]
                    period = periods[best_peak_idx]
                    
                    # Vérifier la significativité
                    if power[best_peak_idx] > 0.5:  # Seuil arbitraire
                        return period
            
            elif method == 'autocorr':
                # Autocorrélation simple
                flux_norm = (flux - np.mean(flux)) / np.std(flux)
                
                # Calculer l'autocorrélation
                autocorr = np.correlate(flux_norm, flux_norm, mode='full')
                autocorr = autocorr[len(autocorr)//2:]
                autocorr = autocorr / autocorr[0]  # Normaliser
                
                # Trouver les pics
                peaks, _ = find_peaks(autocorr[1:], height=0.3, distance=5)
                
                if len(peaks) > 0:
                    # Prendre le premier pic significatif
                    lag = peaks[0] + 1
                    if lag < len(time_norm):
                        period = time_norm[lag] - time_norm[0]
                        if self.min_period <= period <= self.max_period:
                            return period
        
        except Exception as e:
            logger.debug(f"Erreur lors de la recherche de période : {e}")
        
        return None
    
    def is_transient(self, time: np.ndarray, flux: np.ndarray,
                    threshold_sigma: float = 3.0) -> bool:
        """
        Détecte si une courbe de lumière présente un comportement transitoire.
        
        Parameters
        ----------
        time : np.ndarray
            Temps d'observation
        flux : np.ndarray
            Flux observé
        threshold_sigma : float
            Seuil en multiples de sigma pour détecter une variation brusque
            
        Returns
        -------
        bool
            True si transitoire détecté
        """
        if len(flux) < 5:
            return False
        
        # Normaliser
        flux_norm = (flux - np.mean(flux)) / np.std(flux)
        
        # Détecter les variations brusques
        # Calculer la dérivée (différence entre points consécutifs)
        diff = np.diff(flux_norm)
        
        # Vérifier s'il y a des sauts importants
        if len(diff) > 0:
            max_jump = np.max(np.abs(diff))
            if max_jump > threshold_sigma:
                return True
        
        # Vérifier s'il y a un pic ou un creux isolé
        peaks, _ = find_peaks(flux_norm, height=threshold_sigma)
        troughs, _ = find_peaks(-flux_norm, height=threshold_sigma)
        
        if len(peaks) > 0 or len(troughs) > 0:
            # Vérifier que le pic/creux n'est pas périodique
            # (si c'est périodique, ce n'est pas un transitoire)
            return True
        
        return False
    
    def cross_identify_simbad(self, star_coords: List[SkyCoord],
                              radius: float = 5.0) -> Dict[str, Dict]:
        """
        Identifie les étoiles avec SIMBAD.
        
        Parameters
        ----------
        star_coords : List[SkyCoord]
            Liste des coordonnées des étoiles
        radius : float
            Rayon de recherche en arcsecondes
            
        Returns
        -------
        Dict[str, Dict]
            Dictionnaire {star_id: simbad_data}
        """
        logger.info(f"Identification SIMBAD pour {len(star_coords)} étoiles")
        
        identifications = {}
        custom_simbad = Simbad()
        custom_simbad.add_votable_fields('otype', 'sp', 'flux(V)', 'flux(B)')
        
        for i, coord in enumerate(star_coords):
            try:
                result = custom_simbad.query_region(coord, radius=f"{radius}s")
                
                if result is not None and len(result) > 0:
                    identifications[f"star_{i:04d}"] = {
                        'main_id': result['MAIN_ID'][0],
                        'otype': result['OTYPE'][0] if 'OTYPE' in result.colnames else None,
                        'sp_type': result['SP_TYPE'][0] if 'SP_TYPE' in result.colnames else None,
                        'v_mag': result['FLUX_V'][0] if 'FLUX_V' in result.colnames else None,
                        'b_mag': result['FLUX_B'][0] if 'FLUX_B' in result.colnames else None
                    }
            except Exception as e:
                logger.debug(f"Erreur identification SIMBAD pour étoile {i}: {e}")
                continue
        
        logger.info(f"{len(identifications)} étoiles identifiées dans SIMBAD")
        return identifications
    
    def generate_report(self, light_curves: Dict[str, Dict],
                       classifications: Dict[str, List[Dict]],
                       simbad_ids: Optional[Dict[str, Dict]] = None,
                       star_coordinates: Optional[Dict[str, Dict]] = None) -> str:
        """
        Génère un rapport d'analyse.
        
        Parameters
        ----------
        light_curves : Dict[str, Dict]
            Dictionnaire des courbes de lumière
        classifications : Dict[str, List[Dict]]
            Classifications des variables
        simbad_ids : Dict[str, Dict], optional
            Identifications SIMBAD
        star_coordinates : Dict[str, Dict], optional
            Coordonnées des étoiles {star_id: {'ra': float, 'dec': float}}
            
        Returns
        -------
        str
            Rapport formaté
        """
        report = "="*60 + "\n"
        report += "RAPPORT D'ANALYSE DES COURBES DE LUMIÈRE (LCA)\n"
        report += "="*60 + "\n\n"
        
        report += f"Nombre total de courbes de lumière analysées : {len(light_curves)}\n"
        report += f"Nombre d'étoiles variables détectées : {len(self.variable_candidates)}\n\n"
        
        report += f"Classification :\n"
        report += f"  - Périodiques : {len(classifications['periodic'])}\n"
        report += f"  - Transitoires : {len(classifications['transient'])}\n"
        report += f"  - Particulières : {len(classifications['peculiar'])}\n\n"
        
        # Variables périodiques avec coordonnées
        if len(classifications['periodic']) > 0:
            report += "Variables périodiques :\n"
            for var in classifications['periodic'][:10]:  # Limiter à 10
                star_id = var['star_id']
                report += f"  - {star_id} : St={var['st']:.4f}, "
                report += f"Amplitude={var['amplitude']:.4f}, Période={var.get('period', 'N/A'):.2f} j"
                
                # Ajouter les coordonnées si disponibles
                if star_coordinates and star_id in star_coordinates:
                    coords = star_coordinates[star_id]
                    ra = coords.get('ra', None)
                    dec = coords.get('dec', None)
                    if ra is not None and dec is not None:
                        # Formater en heures:minutes:secondes pour RA
                        ra_h = int(ra / 15)
                        ra_m = int((ra / 15 - ra_h) * 60)
                        ra_s = ((ra / 15 - ra_h) * 60 - ra_m) * 60
                        # Formater en degrés:minutes:secondes pour Dec
                        dec_sign = '+' if dec >= 0 else '-'
                        dec_abs = abs(dec)
                        dec_d = int(dec_abs)
                        dec_m = int((dec_abs - dec_d) * 60)
                        dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
                        report += f"\n    RA: {ra_h:02d}h{ra_m:02d}m{ra_s:05.2f}s ({ra:.6f}°) "
                        report += f"Dec: {dec_sign}{dec_d:02d}°{dec_m:02d}'{dec_s:05.2f}\" ({dec:+.6f}°)"
                report += "\n"
            if len(classifications['periodic']) > 10:
                report += f"  ... et {len(classifications['periodic']) - 10} autres\n"
            report += "\n"
        
        # Variables transitoires avec coordonnées
        if len(classifications['transient']) > 0:
            report += "Variables transitoires :\n"
            for var in classifications['transient'][:10]:
                star_id = var['star_id']
                report += f"  - {star_id} : St={var['st']:.4f}, Amplitude={var['amplitude']:.4f}"
                if star_coordinates and star_id in star_coordinates:
                    coords = star_coordinates[star_id]
                    ra = coords.get('ra', None)
                    dec = coords.get('dec', None)
                    if ra is not None and dec is not None:
                        ra_h = int(ra / 15)
                        ra_m = int((ra / 15 - ra_h) * 60)
                        ra_s = ((ra / 15 - ra_h) * 60 - ra_m) * 60
                        dec_sign = '+' if dec >= 0 else '-'
                        dec_abs = abs(dec)
                        dec_d = int(dec_abs)
                        dec_m = int((dec_abs - dec_d) * 60)
                        dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
                        report += f"\n    RA: {ra_h:02d}h{ra_m:02d}m{ra_s:05.2f}s ({ra:.6f}°) "
                        report += f"Dec: {dec_sign}{dec_d:02d}°{dec_m:02d}'{dec_s:05.2f}\" ({dec:+.6f}°)"
                report += "\n"
            if len(classifications['transient']) > 10:
                report += f"  ... et {len(classifications['transient']) - 10} autres\n"
            report += "\n"
        
        # Variables particulières avec coordonnées
        if len(classifications['peculiar']) > 0:
            report += "Variables particulières :\n"
            for var in classifications['peculiar'][:10]:
                star_id = var['star_id']
                report += f"  - {star_id} : St={var['st']:.4f}, Amplitude={var['amplitude']:.4f}"
                if star_coordinates and star_id in star_coordinates:
                    coords = star_coordinates[star_id]
                    ra = coords.get('ra', None)
                    dec = coords.get('dec', None)
                    if ra is not None and dec is not None:
                        ra_h = int(ra / 15)
                        ra_m = int((ra / 15 - ra_h) * 60)
                        ra_s = ((ra / 15 - ra_h) * 60 - ra_m) * 60
                        dec_sign = '+' if dec >= 0 else '-'
                        dec_abs = abs(dec)
                        dec_d = int(dec_abs)
                        dec_m = int((dec_abs - dec_d) * 60)
                        dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
                        report += f"\n    RA: {ra_h:02d}h{ra_m:02d}m{ra_s:05.2f}s ({ra:.6f}°) "
                        report += f"Dec: {dec_sign}{dec_d:02d}°{dec_m:02d}'{dec_s:05.2f}\" ({dec:+.6f}°)"
                report += "\n"
            if len(classifications['peculiar']) > 10:
                report += f"  ... et {len(classifications['peculiar']) - 10} autres\n"
            report += "\n"
        
        if simbad_ids:
            report += f"Identifications SIMBAD : {len(simbad_ids)} étoiles connues\n\n"
        
        report += "="*60 + "\n"
        
        return report
