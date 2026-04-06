# core/binary_star_reduction.py
"""
Module pour la réduction d'images d'étoiles binaires inspiré des techniques REDUC
(http://www.astrosurf.com/hfosaf/reduc/tutoriel.htm)

Techniques implémentées:
- BestOf: Tri des meilleures images basé sur FWHM, contraste, etc.
- Alignement sub-pixel: Alignement précis pour stacking
- Centroiding: Mesure précise des positions et séparations
- ELI (Easy Lucky Imaging): Amélioration de la résolution
- Mesures astrométriques et LITE (effet temps de parcours) pour ajustement orbital combiné

Toutes les dates sont renvoyées en BJD-TDB (préféré à HJD pour éviter les ambiguïtés).
Référence: Zasche & Wolf (2007), Astron. Nachr. 328, 928.
"""
import logging
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Union, Callable
from dataclasses import dataclass
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from scipy import ndimage
from scipy.optimize import minimize
import warnings

logger = logging.getLogger(__name__)

try:
    from photutils.detection import DAOStarFinder
    from photutils.centroids import centroid_1dg, centroid_2dg, centroid_com
    PHOTUTILS_AVAILABLE = True
except ImportError:
    PHOTUTILS_AVAILABLE = False
    logger.warning("photutils non disponible. Certaines fonctionnalités seront limitées.")

try:
    from astropy.wcs import WCS
    from astropy.stats import sigma_clip
    WCS_AVAILABLE = True
except ImportError:
    WCS_AVAILABLE = False


def _to_bjd_tdb(
    jd_or_hjd: float,
    target_coord: Optional[Union[SkyCoord, Tuple[float, float]]] = None,
    observer: Optional[EarthLocation] = None,
) -> float:
    """
    Convertit une date JD ou HJD en BJD-TDB.

    Si target_coord et observer sont fournis, la correction barycentrique
    est appliquée. Sinon, la date est retournée telle quelle (supposée déjà BJD-TDB).

    Parameters
    ----------
    jd_or_hjd : float
        Julian Date ou HJD (Attention: HJD négatif déconseillé, préférer BJD-TDB).
    target_coord : SkyCoord ou tuple (ra_deg, dec_deg), optional
        Coordonnées célestes de la cible (nécessaires pour conversion correcte).
    observer : EarthLocation, optional
        Position de l'observateur (par défaut Greenwich).

    Returns
    -------
    float
        Date en BJD-TDB.
    """
    if jd_or_hjd < 0:
        warnings.warn(
            "HJD négatif fourni. Préférer BJD-TDB pour éviter les ambiguïtés. "
            "Fournir target_coord et observer pour conversion correcte.",
            UserWarning,
            stacklevel=2,
        )
    t = Time(jd_or_hjd, format="jd", scale="utc")
    if target_coord is not None and observer is not None:
        if isinstance(target_coord, (list, tuple)) and len(target_coord) == 2:
            target_coord = SkyCoord(
                ra=target_coord[0],
                dec=target_coord[1],
                unit="deg",
            )
        ltt_bary = t.light_travel_time(target_coord, location=observer)
        bjd_tdb = t.tdb + ltt_bary
        return float(bjd_tdb.jd)
    return float(t.jd)


@dataclass
class AstrometricMeasurement:
    """Mesure astrométrique (θ, ρ) à un instant donné en BJD-TDB."""
    bjd_tdb: float
    theta_deg: float
    rho_arcsec: float
    sigma_theta: float = 0.0
    sigma_rho: float = 0.0
    source: str = ""


@dataclass
class LiteMeasurement:
    """Temps de minimum (LITE) à un instant donné en BJD-TDB."""
    bjd_tdb: float
    sigma_days: float = 0.0
    primary: bool = True
    source: str = ""


def _eccentric_anomaly(M: float, e: float) -> float:
    """Résout E - e*sin(E) = M (Newton-Raphson)."""
    E = M
    for _ in range(20):
        d = E - e * np.sin(E) - M
        if abs(d) < 1e-12:
            break
        E = E - d / (1 - e * np.cos(E))
    return E


def _true_anomaly(E: float, e: float) -> float:
    """Anomalie vraie ν à partir de E."""
    nu = 2 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E / 2),
        np.sqrt(1 - e) * np.cos(E / 2),
    )
    return nu


class _ThirdBodyOrbitModel:
    """
    Modèle orbital du troisième corps pour astrométrie + LITE.
    Paramètres: P3_yr, e, omega_deg, Omega_deg, i_deg, a_mas, A_days, T0_jd, JD0, P_days, q.
    """

    def __init__(self, params: np.ndarray):
        (
            self.P3_yr,
            self.e,
            self.omega_deg,
            self.Omega_deg,
            self.i_deg,
            self.a_mas,
            self.A_days,
            self.T0_jd,
            self.JD0,
            self.P_days,
            self.q,
        ) = params

    def theta_rho(self, t_jd: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Retourne θ (deg) et ρ (mas) à l'instant t (JD)."""
        P3_days = self.P3_yr * 365.25
        M = 2 * np.pi * (t_jd - self.T0_jd) / P3_days
        E = np.array([_eccentric_anomaly(m, self.e) for m in M])
        nu = _true_anomaly(E, self.e)
        r = self.a_mas * (1 - self.e**2) / (1 + self.e * np.cos(nu))
        o = np.radians(self.omega_deg)
        O = np.radians(self.Omega_deg)
        i = np.radians(self.i_deg)
        X = r * (np.cos(o + nu) * np.cos(O) - np.sin(o + nu) * np.sin(O) * np.cos(i))
        Y = r * (np.cos(o + nu) * np.sin(O) + np.sin(o + nu) * np.cos(O) * np.cos(i))
        theta = np.degrees(np.arctan2(Y, X))
        theta = np.where(theta < 0, theta + 360, theta)
        rho = np.sqrt(X**2 + Y**2)
        return theta, rho

    def oc_lite(self, t_jd: np.ndarray) -> np.ndarray:
        """Retourne O-C (jours) LITE à l'instant t (JD)."""
        P3_days = self.P3_yr * 365.25
        M = 2 * np.pi * (t_jd - self.T0_jd) / P3_days
        E = np.array([_eccentric_anomaly(m, self.e) for m in M])
        nu = _true_anomaly(E, self.e)
        o = np.radians(self.omega_deg)
        denom = 1 + self.e * np.cos(nu)
        oc = self.A_days * np.sqrt(1 - self.e**2) * (
            np.cos(o) * (np.cos(E) - self.e) + np.sin(o) * np.sin(E)
        ) / denom
        return oc

    def ephemeris(self, t_jd: np.ndarray) -> np.ndarray:
        """Temps de minimum calculé (éphéméride linéaire + quadratique)."""
        n = np.round((t_jd - self.JD0) / self.P_days).astype(int)
        return self.JD0 + n * self.P_days + self.q * (n.astype(float) ** 2)


@dataclass
class ThirdBodyFitResult:
    """Résultat d'un ajustement orbital combiné astrométrie + LITE."""
    P3_yr: float
    e: float
    omega_deg: float
    Omega_deg: float
    i_deg: float
    a_mas: float
    A_days: float
    T0_jd: float
    JD0: float
    P_days: float
    q: float = 0.0
    chi2_astr: float = 0.0
    chi2_lite: float = 0.0
    chi2_comb: float = 0.0
    n_astr: int = 0
    n_lite: int = 0
    success: bool = False
    message: str = ""


class BinaryStarReduction:
    """
    Classe principale pour la réduction d'images d'étoiles binaires
    inspirée des techniques REDUC.

    Les mesures astrométriques et LITE sont stockées avec dates en BJD-TDB.
    Pour une conversion correcte, définir target_coord et observer.
    """

    def __init__(
        self,
        target_coord: Optional[Union[SkyCoord, Tuple[float, float]]] = None,
        observer: Optional[EarthLocation] = None,
    ):
        """
        Initialise le processeur de réduction.

        Parameters
        ----------
        target_coord : SkyCoord ou tuple (ra_deg, dec_deg), optional
            Coordonnées célestes de la cible. Requis pour convertir JD/HJD en BJD-TDB.
            Utiliser set_target_and_observer() si non fourni à l'init.
        observer : EarthLocation, optional
            Position de l'observateur. Par défaut Greenwich.
        """
        self.image_quality_metrics = {}
        self.astrometric_measurements: List[AstrometricMeasurement] = []
        self.lite_measurements: List[LiteMeasurement] = []
        self._target_coord = target_coord
        self._observer = observer or EarthLocation(
            lat=51.4769 * u.deg,
            lon=-0.0005 * u.deg,
            height=0 * u.m,
        )  # Greenwich par défaut

    def set_target_and_observer(
        self,
        target_coord: Union[SkyCoord, Tuple[float, float]],
        observer: Optional[EarthLocation] = None,
    ) -> None:
        """
        Définit les coordonnées de la cible et l'observateur pour conversion BJD-TDB.

        Parameters
        ----------
        target_coord : SkyCoord ou tuple (ra_deg, dec_deg)
            Coordonnées célestes de la cible.
        observer : EarthLocation, optional
            Position de l'observateur.
        """
        self._target_coord = target_coord
        if observer is not None:
            self._observer = observer
        logger.info(
            "Coordonnées cible et observateur définies pour conversion BJD-TDB."
        )

    def add_astrometric_measurement(
        self,
        theta_deg: float,
        rho_arcsec: float,
        jd_or_bjd: float,
        sigma_theta: float = 0.0,
        sigma_rho: float = 0.0,
        source: str = "",
    ) -> AstrometricMeasurement:
        """
        Enregistre une mesure astrométrique (θ, ρ). La date est convertie
        en BJD-TDB si target_coord et observer sont définis (préféré à HJD).

        Parameters
        ----------
        theta_deg : float
            Angle de position en degrés.
        rho_arcsec : float
            Séparation en secondes d'arc.
        jd_or_bjd : float
            Date JD/HJD de l'observation (convertie en BJD-TDB si target_coord défini).
        sigma_theta, sigma_rho : float
            Incertitudes sur θ et ρ.
        source : str
            Identifiant de la source (fichier, etc.).

        Returns
        -------
        AstrometricMeasurement
        """
        if self._target_coord is None:
            logger.debug(
                "target_coord non défini: les dates seront utilisées telles quelles. "
                "Utiliser set_target_and_observer(ra_deg, dec_deg) pour conversion BJD-TDB."
            )
        bjd = _to_bjd_tdb(
            jd_or_bjd,
            target_coord=self._target_coord,
            observer=self._observer,
        )
        m = AstrometricMeasurement(
            bjd_tdb=bjd,
            theta_deg=theta_deg,
            rho_arcsec=rho_arcsec,
            sigma_theta=sigma_theta,
            sigma_rho=sigma_rho,
            source=source,
        )
        self.astrometric_measurements.append(m)
        return m

    def add_lite_measurement(
        self,
        jd_or_bjd: float,
        sigma_days: float = 0.0,
        primary: bool = True,
        source: str = "",
    ) -> LiteMeasurement:
        """
        Enregistre un temps de minimum (LITE). La date est convertie en BJD-TDB
        si target_coord et observer sont définis (préféré à HJD).

        Parameters
        ----------
        jd_or_bjd : float
            Date JD/HJD du minimum (convertie en BJD-TDB si target_coord défini).
        sigma_days : float
            Incertitude sur le temps en jours.
        primary : bool
            True si minimum primaire.
        source : str
            Identifiant de la source.

        Returns
        -------
        LiteMeasurement
        """
        if self._target_coord is None:
            logger.debug(
                "target_coord non défini: les dates seront utilisées telles quelles. "
                "Utiliser set_target_and_observer(ra_deg, dec_deg) pour conversion BJD-TDB."
            )
        bjd = _to_bjd_tdb(
            jd_or_bjd,
            target_coord=self._target_coord,
            observer=self._observer,
        )
        m = LiteMeasurement(
            bjd_tdb=bjd,
            sigma_days=sigma_days,
            primary=primary,
            source=source,
        )
        self.lite_measurements.append(m)
        return m

    def export_measurements(self) -> Table:
        """
        Exporte les mesures astrométriques et LITE en astropy.table.Table.

        Les dates sont en BJD-TDB.

        Returns
        -------
        Table
            Table avec colonnes: type, bjd_tdb, theta_deg, rho_arcsec,
            sigma_theta, sigma_rho, sigma_days, primary, source.
        """
        rows = []
        for m in self.astrometric_measurements:
            rows.append({
                "type": "astrometry",
                "bjd_tdb": m.bjd_tdb,
                "theta_deg": m.theta_deg,
                "rho_arcsec": m.rho_arcsec,
                "sigma_theta": m.sigma_theta,
                "sigma_rho": m.sigma_rho,
                "sigma_days": np.nan,
                "primary": True,
                "source": m.source,
            })
        for m in self.lite_measurements:
            rows.append({
                "type": "lite",
                "bjd_tdb": m.bjd_tdb,
                "theta_deg": np.nan,
                "rho_arcsec": np.nan,
                "sigma_theta": np.nan,
                "sigma_rho": np.nan,
                "sigma_days": m.sigma_days,
                "primary": m.primary,
                "source": m.source,
            })
        return Table(rows=rows)

    def calculate_image_quality(self, image_data: np.ndarray,
                               threshold: float = None) -> Dict[str, float]:
        """
        Calcule des métriques de qualité d'image (inspiré de REDUC BestOf)
        
        Parameters
        ----------
        image_data : np.ndarray
            Données de l'image
        threshold : float, optional
            Seuil pour la détection des sources (si None, calculé automatiquement)
        
        Returns
        -------
        Dict[str, float]
            Dictionnaire avec les métriques:
            - fwhm: FWHM moyen des sources
            - contrast: Contraste de l'image
            - snr: Signal-to-noise ratio
            - n_stars: Nombre de sources détectées
            - score: Score de qualité global (plus élevé = meilleur)
        """
        try:
            # Calculer un seuil si non fourni
            if threshold is None:
                median = np.median(image_data)
                std = np.std(image_data)
                threshold = median + 3 * std
            
            # Détection des sources (méthode simple si photutils non disponible)
            if PHOTUTILS_AVAILABLE:
                # Utiliser DAOStarFinder pour une détection précise
                daofind = DAOStarFinder(fwhm=3.0, threshold=threshold)
                sources = daofind(image_data)
                
                if sources is None or len(sources) == 0:
                    return {
                        'fwhm': np.nan,
                        'contrast': 0.0,
                        'snr': 0.0,
                        'n_stars': 0,
                        'score': 0.0
                    }
                
                # FWHM moyen
                if 'fwhm' in sources.colnames:
                    fwhm_mean = np.median(sources['fwhm'])
                else:
                    # Estimer FWHM depuis la largeur
                    if 'x' in sources.colnames and 'y' in sources.colnames:
                        fwhm_mean = 3.0  # Estimation par défaut
                    else:
                        fwhm_mean = 3.0
                
                # Contraste: différence entre max et median
                contrast = (np.max(image_data) - np.median(image_data)) / np.median(image_data)
                
                # SNR: signal moyen / bruit
                signal = np.median(sources['peak'])
                noise = np.std(image_data[image_data < threshold])
                snr = signal / noise if noise > 0 else 0.0
                
                n_stars = len(sources)
                
            else:
                # Méthode simple sans photutils
                # Détection de pics locaux
                from scipy.ndimage import maximum_filter
                local_maxima = maximum_filter(image_data, size=5) == image_data
                peaks = image_data[local_maxima & (image_data > threshold)]
                
                if len(peaks) == 0:
                    return {
                        'fwhm': np.nan,
                        'contrast': 0.0,
                        'snr': 0.0,
                        'n_stars': 0,
                        'score': 0.0
                    }
                
                fwhm_mean = 3.0  # Estimation par défaut
                contrast = (np.max(image_data) - np.median(image_data)) / np.median(image_data)
                signal = np.median(peaks)
                noise = np.std(image_data[image_data < threshold])
                snr = signal / noise if noise > 0 else 0.0
                n_stars = len(peaks)
            
            # Score de qualité: combinaison des métriques
            # Plus le FWHM est petit et le SNR/contraste élevés, mieux c'est
            if fwhm_mean > 0 and not np.isnan(fwhm_mean):
                score = (snr * contrast * n_stars) / (fwhm_mean ** 2)
            else:
                score = 0.0
            
            return {
                'fwhm': float(fwhm_mean),
                'contrast': float(contrast),
                'snr': float(snr),
                'n_stars': int(n_stars),
                'score': float(score)
            }
            
        except Exception as e:
            logger.error(f"Erreur calcul qualité image: {e}")
            return {
                'fwhm': np.nan,
                'contrast': 0.0,
                'snr': 0.0,
                'n_stars': 0,
                'score': 0.0
            }
    
    def bestof_sort(self, image_files: List[Path],
                   threshold: float = None,
                   top_percent: float = 0.5,
                   progress_callback: Optional[Callable[[int, int], None]] = None) -> List[Tuple[Path, Dict[str, float]]]:
        """
        Trie les images par qualité (technique REDUC BestOf)

        Parameters
        ----------
        image_files : List[Path]
            Liste des fichiers images à trier
        threshold : float, optional
            Seuil pour la détection (si None, calculé automatiquement)
        top_percent : float
            Pourcentage des meilleures images à retourner (0.0-1.0)
        progress_callback : callable, optional
            Callback(current_index, total) appelé à chaque image traitée (thread-safe côté GUI).

        Returns
        -------
        List[Tuple[Path, Dict[str, float]]]
            Liste triée (meilleures d'abord) avec les métriques de qualité
        """
        # Dédupliquer les fichiers (en utilisant le chemin absolu normalisé)
        seen = set()
        unique_files = []
        for img_path in image_files:
            # Utiliser resolve() pour normaliser le chemin et gérer les liens symboliques
            resolved = img_path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                unique_files.append(img_path)
        
        if len(unique_files) < len(image_files):
            logger.info(f"Déduplication: {len(image_files)} fichiers → {len(unique_files)} fichiers uniques")
        
        logger.info(f"Tri BestOf de {len(unique_files)} images...")
        
        results = []
        total = len(unique_files)
        for i, img_path in enumerate(unique_files):
            try:
                with fits.open(img_path) as hdul:
                    image_data = hdul[0].data.astype(float)
                
                metrics = self.calculate_image_quality(image_data, threshold)
                results.append((img_path, metrics))
                
                if progress_callback is not None:
                    progress_callback(i + 1, total)
                if (i + 1) % 10 == 0:
                    logger.info(f"  Traité {i+1}/{total} images...")
                    
            except Exception as e:
                logger.warning(f"Erreur traitement {img_path.name}: {e}")
                results.append((img_path, {'score': 0.0}))
        
        # Trier par score décroissant
        results.sort(key=lambda x: x[1].get('score', 0.0), reverse=True)
        
        # Retourner les meilleures (en évitant les doublons finaux)
        n_keep = max(1, int(len(results) * top_percent))
        
        # Dédupliquer une dernière fois les résultats (au cas où)
        final_results = []
        final_seen = set()
        for path, metrics in results:
            resolved = path.resolve()
            if resolved not in final_seen:
                final_seen.add(resolved)
                final_results.append((path, metrics))
                if len(final_results) >= n_keep:
                    break
        
        logger.info(f"BestOf terminé: {len(final_results)} meilleures images sélectionnées")
        
        return final_results
    
    def _estimate_centroid_uncertainty(
        self, cutout: np.ndarray, x_cen: float, y_cen: float
    ) -> Tuple[float, float]:
        """
        Estime l'incertitude sur le centroïde (σ_x, σ_y) en pixels
        à partir des moments du cutout (bruit, SNR).
        """
        h, w = cutout.shape
        total = np.sum(cutout)
        if total <= 0:
            return 1.0, 1.0
        y_idx, x_idx = np.indices(cutout.shape)
        mx = np.sum(x_idx * cutout) / total
        my = np.sum(y_idx * cutout) / total
        vx = np.sum((x_idx - mx) ** 2 * cutout) / total
        vy = np.sum((y_idx - my) ** 2 * cutout) / total
        noise = np.std(cutout[cutout < np.median(cutout)]) or 1e-6
        peak = np.max(cutout)
        snr = peak / noise if noise > 0 else 1.0
        fwhm_est = 2.355 * np.sqrt(max(vx, vy, 0.5))
        sigma_pix = max(0.3, 0.6 * fwhm_est / max(snr, 0.1))
        return float(sigma_pix), float(sigma_pix)

    def find_centroid(self, image_data: np.ndarray, 
                     x_init: float, y_init: float,
                     box_size: int = 15,
                     method: str = '2dg') -> Tuple[float, float]:
        """
        Trouve le centroïde précis d'une source (technique REDUC)
        
        Parameters
        ----------
        image_data : np.ndarray
            Données de l'image
        x_init, y_init : float
            Position initiale approximative
        box_size : int
            Taille de la boîte de recherche
        method : str
            Méthode: '1dg', '2dg', 'com'
        
        Returns
        -------
        Tuple[float, float]
            Position du centroïde (x, y)
        """
        x, y = int(x_init), int(y_init)
        h, w = image_data.shape
        
        # Extraire une région autour de la position initiale
        x_min = max(0, x - box_size // 2)
        x_max = min(w, x + box_size // 2 + 1)
        y_min = max(0, y - box_size // 2)
        y_max = min(h, y + box_size // 2 + 1)
        
        cutout = image_data[y_min:y_max, x_min:x_max]
        
        if cutout.size == 0:
            return float(x_init), float(y_init)
        
        try:
            if PHOTUTILS_AVAILABLE:
                if method == '1dg':
                    x_cen, y_cen = centroid_1dg(cutout)
                elif method == '2dg':
                    x_cen, y_cen = centroid_2dg(cutout)
                else:  # 'com'
                    x_cen, y_cen = centroid_com(cutout)
                
                # Convertir en coordonnées de l'image complète
                x_final = x_min + x_cen
                y_final = y_min + y_cen
                
                return float(x_final), float(y_final)
            else:
                # Méthode simple: centre de masse
                y_indices, x_indices = np.indices(cutout.shape)
                total = np.sum(cutout)
                if total > 0:
                    x_cen = np.sum(x_indices * cutout) / total
                    y_cen = np.sum(y_indices * cutout) / total
                    return float(x_min + x_cen), float(y_min + y_cen)
                else:
                    return float(x_init), float(y_init)
                    
        except Exception as e:
            logger.warning(f"Erreur centroiding: {e}")
            return float(x_init), float(y_init)
    
    def align_subpixel(self, image1: np.ndarray, image2: np.ndarray,
                      reference_pos: Tuple[float, float],
                      search_box: int = 50) -> Tuple[float, float]:
        """
        Alignement sub-pixel de deux images (technique REDUC)
        
        Parameters
        ----------
        image1 : np.ndarray
            Image de référence
        image2 : np.ndarray
            Image à aligner
        reference_pos : Tuple[float, float]
            Position de référence (x, y) dans image1
        search_box : int
            Taille de la boîte de recherche
        
        Returns
        -------
        Tuple[float, float]
            Décalage (dx, dy) pour aligner image2 sur image1
        """
        x_ref, y_ref = reference_pos
        x_ref, y_ref = int(x_ref), int(y_ref)
        
        # Extraire des régions autour de la position de référence
        h, w = image1.shape
        size = search_box
        
        x1_min = max(0, x_ref - size // 2)
        x1_max = min(w, x_ref + size // 2)
        y1_min = max(0, y_ref - size // 2)
        y1_max = min(h, y_ref + size // 2)
        
        ref_cutout = image1[y1_min:y1_max, x1_min:x1_max]
        
        # Fonction de corrélation croisée
        def correlation_loss(shift):
            dx, dy = shift
            x2_min = max(0, x_ref - size // 2 - int(dx))
            x2_max = min(w, x_ref + size // 2 - int(dx))
            y2_min = max(0, y_ref - size // 2 - int(dy))
            y2_max = min(h, y_ref + size // 2 - int(dy))
            
            if (x2_max <= x2_min) or (y2_max <= y2_min):
                return 1e10
            
            # Ajuster les tailles
            actual_w = min(ref_cutout.shape[1], x2_max - x2_min)
            actual_h = min(ref_cutout.shape[0], y2_max - y2_min)
            
            ref_region = ref_cutout[:actual_h, :actual_w]
            test_region = image2[y2_min:y2_min+actual_h, x2_min:x2_min+actual_w]
            
            if ref_region.shape != test_region.shape:
                return 1e10
            
            # Corrélation croisée (négative pour minimisation)
            corr = np.corrcoef(ref_region.flatten(), test_region.flatten())[0, 1]
            return -corr if not np.isnan(corr) else 1e10
        
        # Optimisation
        try:
            result = minimize(correlation_loss, [0.0, 0.0], method='Powell',
                            bounds=((-5, 5), (-5, 5)))
            dx, dy = result.x
            return float(dx), float(dy)
        except:
            return 0.0, 0.0
    
    def eli_lucky_imaging(self, image_files: List[Path],
                         output_path: Path,
                         reference_pos: Tuple[float, float],
                         top_percent: float = 0.1,
                         method: str = 'median') -> bool:
        """
        ELI (Easy Lucky Imaging) - Stacking des meilleures images
        inspiré de REDUC ELI
        
        Parameters
        ----------
        image_files : List[Path]
            Liste des images à traiter
        output_path : Path
            Chemin de sortie pour l'image empilée
        reference_pos : Tuple[float, float]
            Position de référence pour l'alignement (x, y)
        top_percent : float
            Pourcentage des meilleures images à utiliser (0.0-1.0)
        method : str
            Méthode de stacking: 'median', 'mean', 'sigma_clip'
        
        Returns
        -------
        bool
            True si succès
        """
        logger.info(f"ELI Lucky Imaging sur {len(image_files)} images...")
        
        # 1. BestOf: sélectionner les meilleures images
        best_images = self.bestof_sort(image_files, top_percent=top_percent)
        n_images = len(best_images)
        
        if n_images == 0:
            logger.error("Aucune image de qualité suffisante trouvée")
            return False
        
        logger.info(f"  {n_images} meilleures images sélectionnées pour ELI")
        
        # 2. Charger la première image comme référence
        ref_path = best_images[0][0]
        with fits.open(ref_path) as hdul:
            ref_data = hdul[0].data.astype(float)
            ref_header = hdul[0].header.copy()
        
        aligned_images = [ref_data]
        
        # 3. Aligner et empiler les autres images
        for i, (img_path, metrics) in enumerate(best_images[1:], 1):
            try:
                with fits.open(img_path) as hdul:
                    img_data = hdul[0].data.astype(float)
                
                # Alignement sub-pixel
                dx, dy = self.align_subpixel(ref_data, img_data, reference_pos)
                
                # Appliquer le décalage
                if abs(dx) > 0.1 or abs(dy) > 0.1:
                    from scipy.ndimage import shift
                    img_aligned = shift(img_data, (-dy, -dx), order=1)
                else:
                    img_aligned = img_data
                
                aligned_images.append(img_aligned)
                
                if (i + 1) % 10 == 0:
                    logger.info(f"  Aligné {i+1}/{n_images-1} images...")
                    
            except Exception as e:
                logger.warning(f"Erreur alignement {img_path.name}: {e}")
        
        # 4. Stacking
        logger.info(f"Empilement de {len(aligned_images)} images (méthode: {method})...")
        
        stack_array = np.array(aligned_images)
        
        if method == 'median':
            stacked = np.median(stack_array, axis=0)
        elif method == 'mean':
            stacked = np.mean(stack_array, axis=0)
        elif method == 'sigma_clip':
            try:
                clipped = sigma_clip(stack_array, axis=0, sigma=3.0)
                stacked = np.ma.median(clipped, axis=0).filled(np.nanmedian(stack_array, axis=0))
            except:
                # Fallback si sigma_clip échoue
                stacked = np.median(stack_array, axis=0)
        else:
            stacked = np.median(stack_array, axis=0)
        
        # 5. Sauvegarder
        # Utiliser uniquement des caractères ASCII pour les valeurs FITS
        ref_header['HISTORY'] = f'ELI Lucky Imaging - {n_images} images, method: {method}'
        ref_header['N_STACK'] = (n_images, 'Number of stacked images')
        ref_header['ELI_METHOD'] = (method, 'ELI stacking method')
        
        fits.writeto(output_path, stacked.astype(ref_data.dtype), ref_header, overwrite=True)
        
        logger.info(f"Image ELI sauvegardée: {output_path}")
        return True
    
    def measure_binary_separation(self, image_path: Path,
                                 pos1: Tuple[float, float],
                                 pos2: Tuple[float, float],
                                 pixel_scale: float = 1.0) -> Dict[str, float]:
        """
        Mesure la séparation et l'angle de position d'un système binaire,
        avec estimation des incertitudes (σ_θ, σ_ρ).

        Parameters
        ----------
        image_path : Path
            Chemin vers l'image
        pos1, pos2 : Tuple[float, float]
            Positions initiales approximatives (x, y)
        pixel_scale : float
            Échelle pixel en arcsec/pixel

        Returns
        -------
        Dict[str, float]
            Dictionnaire avec: separation_pix, separation_arcsec, position_angle,
            sigma_theta (deg), sigma_rho (arcsec).
        """
        try:
            with fits.open(image_path) as hdul:
                image_data = hdul[0].data.astype(float)

            box = 15
            x1, y1 = int(pos1[0]), int(pos1[1])
            x2, y2 = int(pos2[0]), int(pos2[1])
            h, w = image_data.shape

            def cutout(x, y):
                xm = max(0, x - box // 2)
                xM = min(w, x + box // 2 + 1)
                ym = max(0, y - box // 2)
                yM = min(h, y + box // 2 + 1)
                return image_data[ym:yM, xm:xM], xm, ym

            c1, xm1, ym1 = cutout(x1, y1)
            c2, xm2, ym2 = cutout(x2, y2)

            # Centroïdes précis
            x1_c, y1_c = self.find_centroid(image_data, pos1[0], pos1[1])
            x2_c, y2_c = self.find_centroid(image_data, pos2[0], pos2[1])

            # Incertitudes centroïdes
            sx1, sy1 = self._estimate_centroid_uncertainty(
                c1, x1_c - xm1, y1_c - ym1
            )
            sx2, sy2 = self._estimate_centroid_uncertainty(
                c2, x2_c - xm2, y2_c - ym2
            )

            # Propagation: dx = x2-x1, dy = y2-y1, σ_dx² = σ_x1² + σ_x2²
            sigma_dx = np.sqrt(sx1**2 + sx2**2)
            sigma_dy = np.sqrt(sy1**2 + sy2**2)

            dx = x2_c - x1_c
            dy = y2_c - y1_c
            separation_pix = np.sqrt(dx**2 + dy**2)

            if separation_pix < 1e-9:
                sigma_rho_pix = 1.0
                sigma_theta_rad = 1.0
            else:
                sigma_rho_pix = np.sqrt(
                    (dx / separation_pix * sigma_dx) ** 2
                    + (dy / separation_pix * sigma_dy) ** 2
                )
                rho_sq = separation_pix ** 2
                dtheta_dx = -dy / rho_sq
                dtheta_dy = dx / rho_sq
                sigma_theta_rad = np.sqrt(
                    (dtheta_dx * sigma_dx) ** 2 + (dtheta_dy * sigma_dy) ** 2
                )

            separation_arcsec = separation_pix * pixel_scale
            sigma_rho_arcsec = sigma_rho_pix * pixel_scale
            sigma_theta_deg = np.degrees(sigma_theta_rad)

            position_angle = np.degrees(np.arctan2(dx, -dy))
            if position_angle < 0:
                position_angle += 360.0

            return {
                "x1": float(x1_c),
                "y1": float(y1_c),
                "x2": float(x2_c),
                "y2": float(y2_c),
                "separation_pix": float(separation_pix),
                "separation_arcsec": float(separation_arcsec),
                "position_angle": float(position_angle),
                "sigma_theta": float(sigma_theta_deg),
                "sigma_rho": float(sigma_rho_arcsec),
            }

        except Exception as e:
            logger.error(f"Erreur mesure séparation: {e}")
            return {
                "separation_pix": 0.0,
                "separation_arcsec": 0.0,
                "position_angle": 0.0,
                "sigma_theta": 0.0,
                "sigma_rho": 0.0,
            }

    def fit_third_body_orbit(
        self,
        p3_init_yr: float = 30.0,
        e_init: float = 0.5,
        bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    ) -> ThirdBodyFitResult:
        """
        Ajuste les paramètres orbitaux du troisième corps par moindres carrés
        combinés astrométrie + LITE (Zasche & Wolf 2007).

        Parameters
        ----------
        p3_init_yr : float
            Période initiale du troisième corps (années).
        e_init : float
            Excentricité initiale.
        bounds : dict, optional
            Bornes sur les paramètres (clés: P3_yr, e, omega, Omega, i, a_mas, A_days, T0, JD0, P_days).

        Returns
        -------
        ThirdBodyFitResult
        """
        n_astr = len(self.astrometric_measurements)
        n_lite = len(self.lite_measurements)
        if n_astr == 0 and n_lite == 0:
            return ThirdBodyFitResult(
                P3_yr=0, e=0, omega_deg=0, Omega_deg=0, i_deg=0,
                a_mas=0, A_days=0, T0_jd=0, JD0=0, P_days=0,
                success=False, message="Aucune mesure astrométrique ou LITE.",
            )

        def objective(params: np.ndarray) -> float:
            model = _ThirdBodyOrbitModel(params)
            chi2 = 0.0
            if n_astr > 0:
                scale = np.sqrt(n_astr)
                for m in self.astrometric_measurements:
                    theta0, rho0 = model.theta_rho(np.array([m.bjd_tdb]))
                    st = max(m.sigma_theta * scale, 1e-6)
                    sr = max(m.sigma_rho * scale, 1e-6)
                    chi2 += ((m.theta_deg - theta0[0]) / st) ** 2
                    chi2 += ((m.rho_arcsec - rho0[0]) / sr) ** 2
            if n_lite > 0:
                scale = np.sqrt(n_lite)
                t_lite = np.array([m.bjd_tdb for m in self.lite_measurements])
                oc = model.oc_lite(t_lite)
                ephem = model.ephemeris(t_lite)
                for i, m in enumerate(self.lite_measurements):
                    o_c = m.bjd_tdb - ephem[i] - oc[i]
                    sm = max(m.sigma_days * scale, 1e-6)
                    chi2 += (o_c / sm) ** 2
            return chi2

        t_ref = 2450000.0
        if n_astr > 0:
            t_ref = np.mean([m.bjd_tdb for m in self.astrometric_measurements])
        elif n_lite > 0:
            t_ref = np.mean([m.bjd_tdb for m in self.lite_measurements])

        P_days_init = 0.5
        JD0_init = t_ref
        x0 = np.array([
            p3_init_yr,
            e_init,
            0.0,
            0.0,
            60.0,
            500.0,
            0.1,
            t_ref - 1000,
            JD0_init,
            P_days_init,
            0.0,
        ])

        try:
            res = minimize(
                objective,
                x0,
                method="Nelder-Mead",
                options={"maxfev": 50000, "xatol": 1e-8, "fatol": 1e-8},
            )
            p = res.x
            chi2_comb = res.fun
            chi2_astr = 0.0
            chi2_lite = 0.0
            model = _ThirdBodyOrbitModel(p)
            if n_astr > 0:
                scale = np.sqrt(n_astr)
                for m in self.astrometric_measurements:
                    theta0, rho0 = model.theta_rho(np.array([m.bjd_tdb]))
                    st = max(m.sigma_theta * scale, 1e-6)
                    sr = max(m.sigma_rho * scale, 1e-6)
                    chi2_astr += ((m.theta_deg - theta0[0]) / st) ** 2
                    chi2_astr += ((m.rho_arcsec - rho0[0]) / sr) ** 2
            if n_lite > 0:
                scale = np.sqrt(n_lite)
                t_lite = np.array([m.bjd_tdb for m in self.lite_measurements])
                oc = model.oc_lite(t_lite)
                ephem = model.ephemeris(t_lite)
                for i, m in enumerate(self.lite_measurements):
                    o_c = m.bjd_tdb - ephem[i] - oc[i]
                    sm = max(m.sigma_days * scale, 1e-6)
                    chi2_lite += (o_c / sm) ** 2

            return ThirdBodyFitResult(
                P3_yr=float(p[0]),
                e=float(p[1]),
                omega_deg=float(p[2]),
                Omega_deg=float(p[3]),
                i_deg=float(p[4]),
                a_mas=float(p[5]),
                A_days=float(p[6]),
                T0_jd=float(p[7]),
                JD0=float(p[8]),
                P_days=float(p[9]),
                q=float(p[10]),
                chi2_astr=chi2_astr,
                chi2_lite=chi2_lite,
                chi2_comb=chi2_comb,
                n_astr=n_astr,
                n_lite=n_lite,
                success=res.success,
                message=str(res.message),
            )
        except Exception as e:
            logger.error(f"Erreur fit_third_body_orbit: {e}")
            return ThirdBodyFitResult(
                P3_yr=0, e=0, omega_deg=0, Omega_deg=0, i_deg=0,
                a_mas=0, A_days=0, T0_jd=0, JD0=0, P_days=0,
                success=False, message=str(e),
            )

    def predict_third_body(
        self, fit: ThirdBodyFitResult, t_jd: Union[float, np.ndarray]
    ) -> Dict[str, np.ndarray]:
        """
        Prédit θ, ρ et O-C pour des instants donnés à partir du résultat d'ajustement.

        Parameters
        ----------
        fit : ThirdBodyFitResult
        t_jd : float ou array
            Dates en BJD-TDB (jours).

        Returns
        -------
        dict avec clés: theta_deg, rho_mas, oc_days
        """
        t = np.atleast_1d(t_jd)
        params = np.array([
            fit.P3_yr, fit.e, fit.omega_deg, fit.Omega_deg, fit.i_deg,
            fit.a_mas, fit.A_days, fit.T0_jd, fit.JD0, fit.P_days, fit.q,
        ])
        model = _ThirdBodyOrbitModel(params)
        theta, rho = model.theta_rho(t)
        oc = model.oc_lite(t)
        return {"theta_deg": theta, "rho_mas": rho, "oc_days": oc}

    def residuals_third_body(self, fit: ThirdBodyFitResult) -> Dict[str, np.ndarray]:
        """
        Calcule les résidus (observé - calculé) pour les mesures enregistrées.

        Returns
        -------
        dict avec clés: astrometry_theta, astrometry_rho, lite_oc
        """
        out = {"astrometry_theta": np.array([]), "astrometry_rho": np.array([]), "lite_oc": np.array([])}
        if len(self.astrometric_measurements) > 0:
            t = np.array([m.bjd_tdb for m in self.astrometric_measurements])
            pred = self.predict_third_body(fit, t)
            out["astrometry_theta"] = np.array([m.theta_deg for m in self.astrometric_measurements]) - pred["theta_deg"]
            out["astrometry_rho"] = np.array([m.rho_arcsec for m in self.astrometric_measurements]) - pred["rho_mas"] / 1000.0
        if len(self.lite_measurements) > 0:
            t = np.array([m.bjd_tdb for m in self.lite_measurements])
            pred = self.predict_third_body(fit, t)
            ephem = _ThirdBodyOrbitModel(np.array([
                fit.P3_yr, fit.e, fit.omega_deg, fit.Omega_deg, fit.i_deg,
                fit.a_mas, fit.A_days, fit.T0_jd, fit.JD0, fit.P_days, fit.q,
            ])).ephemeris(t)
            oc_pred = pred["oc_days"]
            out["lite_oc"] = np.array([m.bjd_tdb for m in self.lite_measurements]) - ephem - oc_pred
        return out
