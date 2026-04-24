"""
Module SPDE (Synchronous Photometry Data Extraction)
Basé sur Dai et al. (2023), Research in Astronomy and Astrophysics, 23:055011

Ce module implémente un pipeline de réduction de données photométriques en pleine trame
qui extrait automatiquement des courbes de lumière pour toutes les étoiles détectées
dans une série temporelle d'images astronomiques.

Les 5 étapes du pipeline :
1. Classification automatique des images
2. Pre-processing (calibration)
3. Quality justification (sélection de l'image de référence)
4. Automatic matching (matching basé sur triangle de repères)
5. Annulus aperture photometry (photométrie pour toutes les étoiles)
"""

import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clipped_stats

# Définir le logger avant les imports conditionnels
logger = logging.getLogger(__name__)

try:
    from photutils.detection import DAOStarFinder
    from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
    from photutils.background import Background2D, MedianBackground
    PHOTUTILS_AVAILABLE = True
except ImportError:
    PHOTUTILS_AVAILABLE = False
    logger.warning("photutils non disponible. Le pipeline SPDE nécessite photutils.")
from scipy.spatial.distance import cdist
from scipy.optimize import minimize
try:
    import ccdproc
    from astropy.nddata import CCDData
    CCDPROC_AVAILABLE = True
except ImportError:
    CCDPROC_AVAILABLE = False
    # Le logger est déjà défini maintenant
    logger.warning("ccdproc non disponible. Le pre-processing sera limité.")
    logger.info("Pour installer ccdproc : pip install ccdproc>=2.4.0")


class SPDEPipeline:
    """
    Pipeline SPDE pour l'extraction automatique de courbes de lumière
    à partir d'une série temporelle d'images astronomiques.
    """
    
    def __init__(self, 
                 sigma: float = 3.0,
                 fwhm: float = 3.0,
                 threshold: float = 5.0,
                 n_max_brightest: int = 50):
        """
        Parameters
        ----------
        sigma : float
            Paramètre pour la soustraction du fond (DAOPHOT)
        fwhm : float
            Largeur à mi-hauteur pour la détection d'étoiles
        threshold : float
            Seuil inférieur de détection (multiples de sigma)
        n_max_brightest : int
            Nombre d'étoiles les plus brillantes à utiliser pour le matching
        """
        self.sigma = sigma
        self.fwhm = fwhm
        self.threshold = threshold
        self.n_max_brightest = n_max_brightest
        
        self.reference_image_idx = None
        self.reference_stars = None
        self.reference_coords = None
        self.landmark_triangle = None
        self.all_stars_coords = None
        self.light_curves = {}
        
    def step1_classify_images(self, fits_files: List[Path]) -> Dict[str, any]:
        """
        Étape 1 : Classification automatique des images.
        
        Vérifie la présence de WCS, la qualité des headers, etc.
        
        Parameters
        ----------
        fits_files : List[Path]
            Liste des fichiers FITS à classifier
            
        Returns
        -------
        dict
            Dictionnaire avec les résultats de classification
        """
        logger.info(f"Étape 1 : Classification de {len(fits_files)} images")
        
        classified = {
            'valid': [],
            'invalid': [],
            'has_wcs': [],
            'no_wcs': []
        }
        
        for fpath in fits_files:
            try:
                with fits.open(fpath) as hdul:
                    header = hdul[0].header
                    
                    # Vérifier la présence de WCS
                    try:
                        wcs = WCS(header)
                        if wcs.is_celestial:
                            classified['has_wcs'].append(fpath)
                        else:
                            classified['no_wcs'].append(fpath)
                    except Exception:
                        classified['no_wcs'].append(fpath)
                    
                    # Vérifications basiques
                    if 'NAXIS' in header and header['NAXIS'] == 2:
                        classified['valid'].append(fpath)
                    else:
                        classified['invalid'].append(fpath)
                        
            except Exception as e:
                logger.warning(f"Erreur lors de la classification de {fpath.name}: {e}")
                classified['invalid'].append(fpath)
        
        logger.info(f"Classification terminée : {len(classified['valid'])} valides, "
                   f"{len(classified['has_wcs'])} avec WCS")
        
        return classified
    
    def step2_preprocessing(self, fits_files: List[Path], 
                           apply_dark: bool = False,
                           apply_flat: bool = False,
                           dark_path: Optional[Path] = None,
                           flat_path: Optional[Path] = None) -> List:
        """
        Étape 2 : Pre-processing (calibration des images).
        
        Parameters
        ----------
        fits_files : List[Path]
            Liste des fichiers FITS à traiter
        apply_dark : bool
            Appliquer la soustraction de dark
        apply_flat : bool
            Appliquer la correction de flat-field
        dark_path : Path, optional
            Chemin vers l'image dark
        flat_path : Path, optional
            Chemin vers l'image flat
            
        Returns
        -------
        List[CCDData]
            Liste des images calibrées
        """
        logger.info(f"Étape 2 : Pre-processing de {len(fits_files)} images")
        
        calibrated_images = []
        
        # Charger dark et flat si nécessaire
        dark_data = None
        flat_data = None
        
        if apply_dark and dark_path:
            try:
                dark_ccd = CCDData.read(dark_path)
                dark_data = dark_ccd.data
                logger.info(f"Dark chargé : {dark_path.name}")
            except Exception as e:
                logger.warning(f"Impossible de charger le dark : {e}")
        
        if apply_flat and flat_path:
            try:
                flat_ccd = CCDData.read(flat_path)
                flat_data = flat_ccd.data
                logger.info(f"Flat chargé : {flat_path.name}")
            except Exception as e:
                logger.warning(f"Impossible de charger le flat : {e}")
        
        for fpath in fits_files:
            try:
                if not CCDPROC_AVAILABLE:
                    # Fallback : lire directement avec astropy
                    from astropy.io import fits
                    with fits.open(fpath) as hdul:
                        data = hdul[0].data
                        header = hdul[0].header
                    # Créer un objet similaire à CCDData
                    class SimpleCCD:
                        def __init__(self, data, header):
                            self.data = data
                            self.header = header
                    ccd = SimpleCCD(data, header)
                else:
                    ccd = CCDData.read(fpath)
                
                # S'assurer que ccd a les attributs data et header
                if not hasattr(ccd, 'data'):
                    # Si c'est un tableau numpy, créer un wrapper
                    if isinstance(ccd, np.ndarray):
                        data = ccd
                        header = {}
                        class SimpleCCD:
                            def __init__(self, data, header):
                                self.data = data
                                self.header = header
                        ccd = SimpleCCD(data, header)
                
                # Soustraction de dark
                if dark_data is not None and CCDPROC_AVAILABLE:
                    ccd = ccdproc.subtract_dark(ccd, dark_data, dark_exposure=None, 
                                                data_exposure=None, scale=True)
                elif dark_data is not None:
                    # Fallback simple
                    if hasattr(ccd, 'data'):
                        ccd.data = ccd.data - dark_data
                
                # Correction de flat-field
                if flat_data is not None and CCDPROC_AVAILABLE:
                    ccd = ccdproc.flat_correct(ccd, flat_data)
                elif flat_data is not None:
                    # Fallback simple
                    if hasattr(ccd, 'data'):
                        ccd.data = ccd.data / flat_data
                
                calibrated_images.append(ccd)
                
            except Exception as e:
                logger.warning(f"Erreur lors du pre-processing de {fpath.name}: {e}")
                # Ajouter quand même l'image non calibrée
                try:
                    ccd = CCDData.read(fpath)
                    calibrated_images.append(ccd)
                except Exception:
                    pass
        
        logger.info(f"Pre-processing terminé : {len(calibrated_images)} images calibrées")
        return calibrated_images
    
    def step3_quality_justification(self, calibrated_images: List,
                                    fits_files: List[Path]) -> int:
        """
        Étape 3 : Quality justification (sélection de l'image de référence).
        
        L'image de référence est celle avec le plus grand nombre d'étoiles détectées
        par DAOPHOT.
        
        Parameters
        ----------
        calibrated_images : List[CCDData]
            Liste des images calibrées
        fits_files : List[Path]
            Liste des chemins des fichiers FITS correspondants
            
        Returns
        -------
        int
            Index de l'image de référence
        """
        logger.info("Étape 3 : Quality justification (sélection image de référence)")
        
        max_stars = 0
        ref_idx = 0
        
        for idx, ccd in enumerate(calibrated_images):
            try:
                data = ccd.data
                
                # Soustraction du fond
                if not PHOTUTILS_AVAILABLE:
                    raise ImportError("photutils est requis pour le pipeline SPDE")
                bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                                  bkg_estimator=MedianBackground())
                data_sub = data - bkg.background
                
                # Détection des étoiles avec DAOPHOT
                mean, median, std = sigma_clipped_stats(data_sub)
                daofind = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold * std)
                sources = daofind(data_sub)
                
                n_stars = len(sources) if sources is not None else 0
                
                if n_stars > max_stars:
                    max_stars = n_stars
                    ref_idx = idx
                    
            except Exception as e:
                logger.debug(f"Erreur lors de l'analyse de l'image {idx}: {e}")
                continue
        
        self.reference_image_idx = ref_idx
        logger.info(f"Image de référence sélectionnée : {fits_files[ref_idx].name} "
                   f"(N_max = {max_stars} étoiles)")
        
        return ref_idx
    
    def _detect_stars_reference(self, ccd, header=None) -> Table:
        """
        Détecte les étoiles sur l'image de référence.
        
        Parameters
        ----------
        ccd : CCDData or array-like
            Image de référence
        header : fits.Header, optional
            Header FITS avec WCS pour calculer RA/Dec
            
        Returns
        -------
        Table
            Table des étoiles détectées avec colonnes x, y, flux, ra, dec (si WCS disponible)
        """
        data = ccd.data if hasattr(ccd, 'data') else ccd
        
        # Soustraction du fond
        if not PHOTUTILS_AVAILABLE:
            raise ImportError("photutils est requis pour le pipeline SPDE")
        bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                          bkg_estimator=MedianBackground())
        data_sub = data - bkg.background
        
        # Détection DAOPHOT
        mean, median, std = sigma_clipped_stats(data_sub)
        daofind = DAOStarFinder(fwhm=self.fwhm, threshold=self.threshold * std)
        sources = daofind(data_sub)
        
        if sources is None or len(sources) == 0:
            return Table()
        
        # Trier par flux décroissant
        sources.sort('flux', reverse=True)
        
        # Normaliser les noms de colonnes (DAOStarFinder utilise xcentroid/ycentroid)
        if 'xcentroid' in sources.colnames and 'ycentroid' in sources.colnames:
            sources['x'] = sources['xcentroid']
            sources['y'] = sources['ycentroid']
        elif 'x' not in sources.colnames or 'y' not in sources.colnames:
            # Si ni x/y ni xcentroid/ycentroid, créer des colonnes par défaut
            if 'x' not in sources.colnames:
                sources['x'] = sources['xcentroid'] if 'xcentroid' in sources.colnames else 0.0
            if 'y' not in sources.colnames:
                sources['y'] = sources['ycentroid'] if 'ycentroid' in sources.colnames else 0.0
        
        # Calculer RA/Dec si WCS disponible
        if header is not None:
            try:
                wcs = WCS(header)
                if wcs.is_celestial:
                    coords = wcs.pixel_to_world(sources['x'], sources['y'])
                    sources['ra'] = coords.ra.deg
                    sources['dec'] = coords.dec.deg
            except Exception as e:
                logger.debug(f"Impossible de calculer RA/Dec : {e}")
        elif hasattr(ccd, 'header'):
            try:
                wcs = WCS(ccd.header)
                if wcs.is_celestial:
                    coords = wcs.pixel_to_world(sources['x'], sources['y'])
                    sources['ra'] = coords.ra.deg
                    sources['dec'] = coords.dec.deg
            except Exception as e:
                logger.debug(f"Impossible de calculer RA/Dec depuis ccd.header : {e}")
        
        return sources
    
    def _find_landmark_triangle(self, stars: Table, n_candidates: int = 10) -> Tuple[int, int, int]:
        """
        Trouve un triangle de repères parmi les étoiles les plus brillantes.
        
        Le triangle doit être non dégénéré et avoir des côtés de longueurs raisonnables.
        
        Parameters
        ----------
        stars : Table
            Table des étoiles détectées
        n_candidates : int
            Nombre d'étoiles candidates à considérer
            
        Returns
        -------
        Tuple[int, int, int]
            Indices des trois étoiles formant le triangle
        """
        if len(stars) < 3:
            raise ValueError("Pas assez d'étoiles pour former un triangle")
        
        n = min(n_candidates, len(stars))
        candidates = stars[:n]
        
        coords = np.array([candidates['x'], candidates['y']]).T
        
        # Essayer différentes combinaisons de triangles
        best_triangle = None
        best_score = -np.inf
        
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    # Calculer les côtés du triangle
                    d12 = np.linalg.norm(coords[i] - coords[j])
                    d13 = np.linalg.norm(coords[i] - coords[k])
                    d23 = np.linalg.norm(coords[j] - coords[k])
                    
                    # Vérifier que le triangle n'est pas dégénéré
                    if d12 < 5 or d13 < 5 or d23 < 5:
                        continue
                    
                    # Score basé sur la régularité du triangle
                    # (on préfère des triangles équilatéraux)
                    sides = np.array([d12, d13, d23])
                    score = -np.std(sides) / np.mean(sides)  # Plus petit écart relatif = meilleur
                    
                    if score > best_score:
                        best_score = score
                        best_triangle = (i, j, k)
        
        if best_triangle is None:
            # Fallback : prendre les 3 premières étoiles
            best_triangle = (0, 1, 2)
        
        return best_triangle
    
    def step4_automatic_matching(self, calibrated_images: List,
                                 fits_files: List[Path]) -> Dict[int, Dict]:
        """
        Étape 4 : Automatic matching basé sur un triangle de repères.
        
        Parameters
        ----------
        calibrated_images : List[CCDData]
            Liste des images calibrées
        fits_files : List[Path]
            Liste des chemins des fichiers FITS
            
        Returns
        -------
        Dict[int, Dict]
            Dictionnaire {index_image: {'transform': transform_matrix, 'matched_stars': ...}}
        """
        logger.info("Étape 4 : Automatic matching")
        
        # Détecter les étoiles sur l'image de référence
        ref_ccd = calibrated_images[self.reference_image_idx]
        ref_header = ref_ccd.header if hasattr(ref_ccd, 'header') else None
        ref_stars = self._detect_stars_reference(ref_ccd, header=ref_header)
        
        if len(ref_stars) == 0:
            raise ValueError("Aucune étoile détectée sur l'image de référence")
        
        # Sélectionner les N étoiles les plus brillantes
        n_select = min(self.n_max_brightest, len(ref_stars))
        self.reference_stars = ref_stars[:n_select]
        
        # Trouver le triangle de repères
        triangle_indices = self._find_landmark_triangle(self.reference_stars)
        self.landmark_triangle = triangle_indices
        
        ref_triangle_coords = np.array([
            [self.reference_stars[triangle_indices[0]]['x'], 
             self.reference_stars[triangle_indices[0]]['y']],
            [self.reference_stars[triangle_indices[1]]['x'], 
             self.reference_stars[triangle_indices[1]]['y']],
            [self.reference_stars[triangle_indices[2]]['x'], 
             self.reference_stars[triangle_indices[2]]['y']]
        ])
        
        logger.info(f"Triangle de repères sélectionné : indices {triangle_indices}")
        
        # Coordonnées de référence pour toutes les étoiles
        self.reference_coords = np.array([self.reference_stars['x'], 
                                         self.reference_stars['y']]).T
        
        # Matching pour chaque image
        matching_results = {}
        
        for idx, ccd in enumerate(calibrated_images):
            if idx == self.reference_image_idx:
                # Image de référence : transformation identité
                matching_results[idx] = {
                    'transform': np.eye(3),
                    'matched_stars': self.reference_stars,
                    'success': True
                }
                continue
            
            try:
                # Détecter les étoiles sur cette image
                ccd_header = ccd.header if hasattr(ccd, 'header') else None
                stars = self._detect_stars_reference(ccd, header=ccd_header)
                
                if len(stars) < 3:
                    matching_results[idx] = {'success': False, 'reason': 'Pas assez d\'étoiles'}
                    continue
                
                # Trouver le triangle correspondant
                triangle_found = self._match_triangle(ref_triangle_coords, stars)
                
                if triangle_found is None:
                    matching_results[idx] = {'success': False, 'reason': 'Triangle non trouvé'}
                    continue
                
                # Calculer la transformation affine
                transform = self._compute_affine_transform(ref_triangle_coords, triangle_found)
                
                # Appliquer la transformation à toutes les étoiles de référence
                ref_coords_homogeneous = np.column_stack([self.reference_coords, 
                                                         np.ones(len(self.reference_coords))])
                transformed_coords = (transform @ ref_coords_homogeneous.T).T[:, :2]
                
                # Matcher les étoiles transformées avec les étoiles détectées
                matched_stars = self._match_stars(transformed_coords, stars, max_distance=2.0)
                
                matching_results[idx] = {
                    'transform': transform,
                    'matched_stars': matched_stars,
                    'success': True
                }
                
            except Exception as e:
                logger.warning(f"Erreur lors du matching de l'image {idx}: {e}")
                matching_results[idx] = {'success': False, 'reason': str(e)}
        
        n_success = sum(1 for r in matching_results.values() if r.get('success', False))
        logger.info(f"Matching terminé : {n_success}/{len(calibrated_images)} images matchées")
        
        return matching_results
    
    def _match_triangle(self, ref_triangle: np.ndarray, stars: Table, 
                       tolerance: float = 0.1) -> Optional[np.ndarray]:
        """
        Trouve le triangle correspondant dans une image.
        
        Parameters
        ----------
        ref_triangle : np.ndarray
            Coordonnées du triangle de référence (3x2)
        stars : Table
            Table des étoiles détectées
        tolerance : float
            Tolérance relative pour la correspondance des distances
            
        Returns
        -------
        np.ndarray or None
            Coordonnées du triangle trouvé (3x2) ou None
        """
        # Calculer les distances du triangle de référence
        ref_distances = np.array([
            np.linalg.norm(ref_triangle[0] - ref_triangle[1]),
            np.linalg.norm(ref_triangle[0] - ref_triangle[2]),
            np.linalg.norm(ref_triangle[1] - ref_triangle[2])
        ])
        ref_distances = np.sort(ref_distances)
        
        n = min(20, len(stars))  # Limiter la recherche
        # Normaliser les colonnes si nécessaire
        if 'x' not in stars.colnames and 'xcentroid' in stars.colnames:
            stars['x'] = stars['xcentroid']
        if 'y' not in stars.colnames and 'ycentroid' in stars.colnames:
            stars['y'] = stars['ycentroid']
        coords = np.array([stars[:n]['x'], stars[:n]['y']]).T
        
        # Essayer toutes les combinaisons de triangles
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    triangle = np.array([coords[i], coords[j], coords[k]])
                    distances = np.array([
                        np.linalg.norm(triangle[0] - triangle[1]),
                        np.linalg.norm(triangle[0] - triangle[2]),
                        np.linalg.norm(triangle[1] - triangle[2])
                    ])
                    distances = np.sort(distances)
                    
                    # Vérifier la correspondance des ratios de distances
                    if len(ref_distances) == 3 and ref_distances[0] > 0:
                        ratios_ref = ref_distances[1:] / ref_distances[0]
                        ratios = distances[1:] / distances[0]
                        
                        if np.allclose(ratios, ratios_ref, rtol=tolerance):
                            return triangle
        
        return None
    
    def _compute_affine_transform(self, ref_points: np.ndarray, 
                                 target_points: np.ndarray) -> np.ndarray:
        """
        Calcule la transformation affine entre deux ensembles de 3 points.
        
        Parameters
        ----------
        ref_points : np.ndarray
            Points de référence (3x2)
        target_points : np.ndarray
            Points cibles (3x2)
            
        Returns
        -------
        np.ndarray
            Matrice de transformation affine (3x3)
        """
        # Utiliser une transformation affine 2D (6 paramètres)
        # x' = a*x + b*y + c
        # y' = d*x + e*y + f
        
        A = np.zeros((6, 6))
        b = np.zeros(6)
        
        for i in range(3):
            A[2*i, 0] = ref_points[i, 0]
            A[2*i, 1] = ref_points[i, 1]
            A[2*i, 2] = 1
            A[2*i+1, 3] = ref_points[i, 0]
            A[2*i+1, 4] = ref_points[i, 1]
            A[2*i+1, 5] = 1
            
            b[2*i] = target_points[i, 0]
            b[2*i+1] = target_points[i, 1]
        
        params = np.linalg.solve(A, b)
        
        transform = np.array([
            [params[0], params[1], params[2]],
            [params[3], params[4], params[5]],
            [0, 0, 1]
        ])
        
        return transform
    
    def _match_stars(self, transformed_coords: np.ndarray, stars: Table,
                    max_distance: float = 2.0) -> Table:
        """
        Matche les coordonnées transformées avec les étoiles détectées.
        
        Parameters
        ----------
        transformed_coords : np.ndarray
            Coordonnées transformées (Nx2)
        stars : Table
            Table des étoiles détectées
        max_distance : float
            Distance maximale pour un match (en pixels)
            
        Returns
        -------
        Table
            Table des étoiles matchées avec leurs coordonnées transformées
        """
        # Normaliser les colonnes si nécessaire
        if 'x' not in stars.colnames and 'xcentroid' in stars.colnames:
            stars['x'] = stars['xcentroid']
        if 'y' not in stars.colnames and 'ycentroid' in stars.colnames:
            stars['y'] = stars['ycentroid']
        star_coords = np.array([stars['x'], stars['y']]).T
        
        # Calculer les distances
        distances = cdist(transformed_coords, star_coords)
        
        matched = []
        used_stars = set()
        
        for i, trans_coord in enumerate(transformed_coords):
            # Trouver l'étoile la plus proche
            min_dist_idx = np.argmin(distances[i])
            min_dist = distances[i, min_dist_idx]
            
            if min_dist <= max_distance and min_dist_idx not in used_stars:
                matched.append({
                    'x': stars[min_dist_idx]['x'],
                    'y': stars[min_dist_idx]['y'],
                    'flux': stars[min_dist_idx]['flux'],
                    'x_ref': transformed_coords[i, 0],
                    'y_ref': transformed_coords[i, 1],
                    'distance': min_dist
                })
                used_stars.add(min_dist_idx)
        
        return Table(matched) if matched else Table()
    
    def step5_aperture_photometry(self, calibrated_images: List,
                                  matching_results: Dict[int, Dict],
                                  aperture_range: Tuple[float, float] = (1.5, 7.0),
                                  n_apertures: int = 12) -> Dict[str, np.ndarray]:
        """
        Étape 5 : Annulus aperture photometry pour toutes les étoiles.
        
        Parameters
        ----------
        calibrated_images : List[CCDData]
            Liste des images calibrées
        matching_results : Dict[int, Dict]
            Résultats du matching de l'étape 4
        aperture_range : Tuple[float, float]
            Plage de tailles d'ouverture à tester (en pixels)
        n_apertures : int
            Nombre de tailles d'ouverture à tester
            
        Returns
        -------
        Dict[str, np.ndarray]
            Dictionnaire des courbes de lumière {star_id: light_curve}
        """
        logger.info("Étape 5 : Annulus aperture photometry")
        
        # Déterminer les tailles d'ouverture optimales pour chaque étoile
        # sur l'image de référence
        ref_ccd = calibrated_images[self.reference_image_idx]
        ref_matching = matching_results[self.reference_image_idx]
        
        if len(ref_matching['matched_stars']) == 0:
            raise ValueError("Aucune étoile matchée sur l'image de référence")
        
        # Tester différentes tailles d'ouverture
        aperture_sizes = np.linspace(aperture_range[0], aperture_range[1], n_apertures)
        
        # Pour chaque étoile, trouver la taille d'ouverture optimale
        optimal_apertures = {}
        
        for i, star in enumerate(ref_matching['matched_stars']):
            x, y = star['x'], star['y']
            
            # Tester chaque taille d'ouverture
            best_snr = -np.inf
            best_aperture = aperture_sizes[0]
            
            for ap_size in aperture_sizes:
                try:
                    # Photométrie d'ouverture
                    aperture = CircularAperture((x, y), r=ap_size)
                    annulus = CircularAnnulus((x, y), r_in=ap_size*1.5, r_out=ap_size*2.0)
                    
                    phot = aperture_photometry(ref_ccd.data, aperture)
                    bkg_phot = aperture_photometry(ref_ccd.data, annulus)
                    
                    flux = phot['aperture_sum'][0]
                    bkg_mean = bkg_phot['aperture_sum'][0] / annulus.area
                    net_flux = flux - (bkg_mean * aperture.area)
                    
                    # Estimation du SNR
                    if net_flux > 0:
                        snr = net_flux / np.sqrt(net_flux + (bkg_mean * aperture.area))
                        if snr > best_snr:
                            best_snr = snr
                            best_aperture = ap_size
                except Exception:
                    continue
            
            optimal_apertures[i] = best_aperture
        
        # Extraire les courbes de lumière pour toutes les images
        light_curves = {}
        
        for star_idx in range(len(ref_matching['matched_stars'])):
            star_id = f"star_{star_idx:04d}"
            light_curve = []
            times = []
            
            for idx, ccd in enumerate(calibrated_images):
                matching = matching_results.get(idx, {})
                
                if not matching.get('success', False):
                    continue
                
                matched_stars = matching['matched_stars']
                
                if star_idx >= len(matched_stars):
                    continue
                
                star = matched_stars[star_idx]
                x, y = star['x'], star['y']
                ap_size = optimal_apertures.get(star_idx, 3.0)
                
                try:
                    # Photométrie
                    aperture = CircularAperture((x, y), r=ap_size)
                    annulus = CircularAnnulus((x, y), r_in=ap_size*1.5, r_out=ap_size*2.0)
                    
                    phot = aperture_photometry(ccd.data, aperture)
                    bkg_phot = aperture_photometry(ccd.data, annulus)
                    
                    flux = phot['aperture_sum'][0]
                    bkg_mean = bkg_phot['aperture_sum'][0] / annulus.area
                    net_flux = flux - (bkg_mean * aperture.area)
                    
                    light_curve.append(net_flux)
                    
                    # Temps d'observation
                    try:
                        header = ccd.header
                        if 'JD' in header:
                            times.append(header['JD'])
                        elif 'MJD-OBS' in header:
                            times.append(header['MJD-OBS'] + 2400000.5)
                        else:
                            times.append(idx)  # Fallback : index
                    except Exception:
                        times.append(idx)
                        
                except Exception as e:
                    logger.debug(f"Erreur photométrie étoile {star_id} image {idx}: {e}")
                    continue
            
            if len(light_curve) > 0:
                light_curves[star_id] = {
                    'flux': np.array(light_curve),
                    'time': np.array(times),
                    'aperture': optimal_apertures.get(star_idx, 3.0)
                }
        
        self.light_curves = light_curves
        logger.info(f"Photométrie terminée : {len(light_curves)} courbes de lumière extraites")
        
        return light_curves
    
    def run_full_pipeline(self, fits_files: List[Path],
                         apply_dark: bool = False,
                         apply_flat: bool = False,
                         dark_path: Optional[Path] = None,
                         flat_path: Optional[Path] = None) -> Dict[str, any]:
        """
        Exécute le pipeline complet SPDE.
        
        Parameters
        ----------
        fits_files : List[Path]
            Liste des fichiers FITS à traiter
        apply_dark : bool
            Appliquer la soustraction de dark
        apply_flat : bool
            Appliquer la correction de flat-field
        dark_path : Path, optional
            Chemin vers l'image dark
        flat_path : Path, optional
            Chemin vers l'image flat
            
        Returns
        -------
        Dict[str, any]
            Dictionnaire avec tous les résultats du pipeline
        """
        logger.info("="*60)
        logger.info("Démarrage du pipeline SPDE complet")
        logger.info("="*60)
        
        # Étape 1 : Classification
        classification = self.step1_classify_images(fits_files)
        valid_files = classification['valid']
        
        if len(valid_files) == 0:
            raise ValueError("Aucune image valide trouvée")
        
        # Étape 2 : Pre-processing
        if not CCDPROC_AVAILABLE and (apply_dark or apply_flat):
            logger.warning("ccdproc non disponible. Le pre-processing sera ignoré.")
            apply_dark = False
            apply_flat = False
        
        calibrated_images = self.step2_preprocessing(
            valid_files, apply_dark, apply_flat, dark_path, flat_path
        )
        
        # Étape 3 : Quality justification
        ref_idx = self.step3_quality_justification(calibrated_images, valid_files)
        
        # Étape 4 : Automatic matching
        matching_results = self.step4_automatic_matching(calibrated_images, valid_files)
        
        # Étape 5 : Aperture photometry
        light_curves = self.step5_aperture_photometry(calibrated_images, matching_results)
        
        results = {
            'classification': classification,
            'reference_image_idx': ref_idx,
            'reference_image': valid_files[ref_idx],
            'matching_results': matching_results,
            'light_curves': light_curves,
            'n_stars': len(light_curves),
            'n_images': len(valid_files)
        }
        
        logger.info("="*60)
        logger.info("Pipeline SPDE terminé avec succès")
        logger.info(f"  - {results['n_images']} images traitées")
        logger.info(f"  - {results['n_stars']} courbes de lumière extraites")
        logger.info("="*60)
        
        # Stocker les étoiles de référence avec leurs coordonnées
        results['reference_stars'] = self.reference_stars
        
        return results
    
    def save_annotated_reference_image(self, output_path: Path, 
                                      reference_image_path: Path,
                                      variable_star_ids: Optional[List[str]] = None,
                                      mark_all_stars: bool = True):
        """
        Sauvegarde l'image de référence avec les étoiles marquées.
        
        Parameters
        ----------
        output_path : Path
            Chemin de sortie pour l'image FITS annotée
        reference_image_path : Path
            Chemin vers l'image de référence
        variable_star_ids : List[str], optional
            Liste des IDs des étoiles variables à marquer spécialement
        mark_all_stars : bool
            Si True, marque toutes les étoiles détectées
        """
        if self.reference_image_idx is None:
            raise ValueError("Aucune image de référence sélectionnée")
        
        if self.reference_stars is None or len(self.reference_stars) == 0:
            raise ValueError("Aucune étoile de référence disponible")
        
        from astropy.io import fits
        
        with fits.open(reference_image_path) as hdul:
            # Copier les données
            data = hdul[0].data.copy()
            header = hdul[0].header.copy()
            
            # Ajouter les positions des étoiles dans le header
            header['NSTARS'] = (len(self.reference_stars), 'Nombre d\'étoiles détectées')
            
            # Marquer les étoiles variables si spécifiées
            if variable_star_ids:
                header['NVARS'] = (len(variable_star_ids), 'Nombre d\'étoiles variables')
                for i, var_id in enumerate(variable_star_ids[:50]):  # Limiter à 50
                    # Extraire l'index de l'étoile depuis l'ID (format: star_XXXX)
                    try:
                        star_idx = int(var_id.split('_')[1])
                        if star_idx < len(self.reference_stars):
                            star = self.reference_stars[star_idx]
                            header[f'VAR{i}_X'] = (float(star['x']), f'Position X étoile variable {var_id}')
                            header[f'VAR{i}_Y'] = (float(star['y']), f'Position Y étoile variable {var_id}')
                            if 'ra' in star.colnames:
                                header[f'VAR{i}_RA'] = (float(star['ra']), f'RA étoile variable {var_id} (deg)')
                            if 'dec' in star.colnames:
                                header[f'VAR{i}_DEC'] = (float(star['dec']), f'Dec étoile variable {var_id} (deg)')
                    except (ValueError, IndexError):
                        continue
            
            # Créer une table avec toutes les étoiles pour une extension FITS
            stars_table = Table()
            stars_table['X'] = self.reference_stars['x']
            stars_table['Y'] = self.reference_stars['y']
            stars_table['FLUX'] = self.reference_stars['flux']
            if 'ra' in self.reference_stars.colnames:
                stars_table['RA'] = self.reference_stars['ra']
            if 'dec' in self.reference_stars.colnames:
                stars_table['DEC'] = self.reference_stars['dec']
            
            # Ajouter un flag pour les variables
            if variable_star_ids:
                is_variable = np.zeros(len(self.reference_stars), dtype=bool)
                for var_id in variable_star_ids:
                    try:
                        star_idx = int(var_id.split('_')[1])
                        if star_idx < len(is_variable):
                            is_variable[star_idx] = True
                    except (ValueError, IndexError):
                        pass
                stars_table['IS_VAR'] = is_variable
            else:
                stars_table['IS_VAR'] = np.zeros(len(self.reference_stars), dtype=bool)
            
            # Créer un nouveau HDUList
            new_hdul = fits.HDUList()
            new_hdul.append(fits.PrimaryHDU(data=data, header=header))
            
            # Ajouter l'extension avec la table des étoiles
            stars_hdu = fits.BinTableHDU(stars_table, name='STARS')
            new_hdul.append(stars_hdu)
            
            # Sauvegarder
            new_hdul.writeto(output_path, overwrite=True)
            logger.info(f"Image de référence annotée sauvegardée : {output_path}")
