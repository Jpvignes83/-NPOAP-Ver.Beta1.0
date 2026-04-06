# core/binary_star_catalogs.py
"""
Module pour rechercher et charger les paramètres de systèmes binaires depuis des catalogues astronomiques
"""
import logging
from typing import Optional, Dict, Tuple
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import numpy as np

logger = logging.getLogger(__name__)

try:
    from astroquery.vizier import Vizier
    from astroquery.simbad import Simbad
    VIZIER_AVAILABLE = True
except ImportError:
    VIZIER_AVAILABLE = False
    logger.warning("astroquery non disponible. Les catalogues ne pourront pas être interrogés.")


class BinaryStarCatalog:
    """
    Classe pour interroger les catalogues astronomiques pour les systèmes binaires
    """
    
    def __init__(self):
        """Initialise le querier de catalogues"""
        if VIZIER_AVAILABLE:
            self.vizier = Vizier()
            self.vizier.ROW_LIMIT = 10  # Limiter les résultats
            self.simbad = Simbad()
            self.simbad.add_votable_fields('otype', 'sp', 'pmra', 'pmdec', 'parallax', 'flux(V)')
    
    def search_by_name(self, star_name: str) -> Optional[Dict]:
        """
        Recherche un système binaire par nom
        
        Parameters
        ----------
        star_name : str
            Nom de l'étoile (ex: "V1082 Tau", "Beta Lyrae", etc.)
        
        Returns
        -------
        Dict or None
            Dictionnaire avec les paramètres trouvés ou None
        """
        if not VIZIER_AVAILABLE:
            logger.error("astroquery non disponible")
            return None
        
        try:
            # Essayer Simbad d'abord pour obtenir les coordonnées
            result_table = self.simbad.query_object(star_name)
            
            if result_table is None or len(result_table) == 0:
                logger.warning(f"Système {star_name} non trouvé dans Simbad")
                return None
            
            # Extraire les coordonnées
            ra = result_table['RA'][0]
            dec = result_table['DEC'][0]
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
            
            # Chercher dans différents catalogues
            params = {
                'name': star_name,
                'ra': coord.ra.deg,
                'dec': coord.dec.deg,
                'coordinate': coord
            }
            
            # Chercher dans VSX (Variable Star Index)
            vsx_params = self._query_vsx(star_name, coord)
            if vsx_params:
                params.update(vsx_params)
            
            # Chercher dans SB9 (Spectroscopic Binaries)
            sb9_params = self._query_sb9(star_name, coord)
            if sb9_params:
                params.update(sb9_params)
            
            # Chercher dans Gaia DR3
            gaia_params = self._query_gaia_binary(coord)
            if gaia_params:
                params.update(gaia_params)
            
            # Chercher dans DEBCat (Detached Eclipsing Binaries)
            debcat_params = self._query_debcat(star_name, coord)
            if debcat_params:
                params.update(debcat_params)
            
            return params if len(params) > 3 else None  # Plus que name, ra, dec, coord
            
        except Exception as e:
            logger.error(f"Erreur recherche {star_name}: {e}")
            return None
    
    def search_by_coordinates(self, coord: SkyCoord, radius: float = 0.1) -> Optional[Dict]:
        """
        Recherche un système binaire par coordonnées
        
        Parameters
        ----------
        coord : SkyCoord
            Coordonnées du système
        radius : float
            Rayon de recherche en degrés
        
        Returns
        -------
        Dict or None
            Dictionnaire avec les paramètres trouvés
        """
        if not VIZIER_AVAILABLE:
            logger.error("astroquery non disponible")
            return None
        
        try:
            params = {
                'ra': coord.ra.deg,
                'dec': coord.dec.deg,
                'coordinate': coord
            }
            
            # Chercher dans VSX
            vsx_params = self._query_vsx_by_coord(coord, radius)
            if vsx_params:
                params.update(vsx_params)
            
            # Chercher dans SB9
            sb9_params = self._query_sb9_by_coord(coord, radius)
            if sb9_params:
                params.update(sb9_params)
            
            # Chercher dans Gaia
            gaia_params = self._query_gaia_binary(coord)
            if gaia_params:
                params.update(gaia_params)
            
            return params if len(params) > 2 else None
            
        except Exception as e:
            logger.error(f"Erreur recherche coordonnées: {e}")
            return None
    
    def _query_vsx(self, star_name: str, coord: SkyCoord) -> Optional[Dict]:
        """
        Interroge VSX (AAVSO Variable Star Index)
        Catalogue: B/vsx/vsx
        """
        try:
            v = Vizier(columns=['**'], row_limit=1)
            result = v.query_region(coord, radius=0.1*u.deg, catalog='B/vsx/vsx')
            
            if result and len(result) > 0:
                vsx = result[0]
                if len(vsx) > 0:
                    params = {}
                    
                    # Période
                    if 'Period' in vsx.colnames:
                        period = vsx['Period'][0]
                        if not np.isnan(period) and period > 0:
                            params['period'] = float(period)  # jours
                    
                    # Type de variabilité
                    if 'VarType' in vsx.colnames:
                        params['variability_type'] = str(vsx['VarType'][0])
                    
                    # Magnitude maximale/minimale
                    if 'MagMax' in vsx.colnames:
                        mag_max = vsx['MagMax'][0]
                        if not np.isnan(mag_max):
                            params['magnitude_max'] = float(mag_max)
                    
                    if 'MagMin' in vsx.colnames:
                        mag_min = vsx['MagMin'][0]
                        if not np.isnan(mag_min):
                            params['magnitude_min'] = float(mag_min)
                    
                    logger.info(f"VSX: Paramètres trouvés pour {star_name}")
                    return params
        except Exception as e:
            logger.debug(f"VSX query error: {e}")
        return None
    
    def _query_vsx_by_coord(self, coord: SkyCoord, radius: float) -> Optional[Dict]:
        """Cherche dans VSX par coordonnées"""
        try:
            v = Vizier(columns=['**'], row_limit=1)
            result = v.query_region(coord, radius=radius*u.deg, catalog='B/vsx/vsx')
            
            if result and len(result) > 0:
                vsx = result[0]
                if len(vsx) > 0:
                    params = {}
                    if 'Name' in vsx.colnames:
                        params['name'] = str(vsx['Name'][0])
                    if 'Period' in vsx.colnames:
                        period = vsx['Period'][0]
                        if not np.isnan(period) and period > 0:
                            params['period'] = float(period)
                    if 'VarType' in vsx.colnames:
                        params['variability_type'] = str(vsx['VarType'][0])
                    return params
        except:
            pass
        return None
    
    def _query_sb9(self, star_name: str, coord: SkyCoord) -> Optional[Dict]:
        """
        Interroge SB9 (9th Catalogue of Spectroscopic Binary Orbits)
        Catalogue: B/sb9/sb9
        Contient: périodes, excentricités, masses, etc.
        """
        try:
            v = Vizier(columns=['**'], row_limit=1)
            result = v.query_region(coord, radius=0.1*u.deg, catalog='B/sb9/sb9')
            
            if result and len(result) > 0:
                sb9 = result[0]
                if len(sb9) > 0:
                    params = {}
                    
                    # Période orbitale
                    if 'Per' in sb9.colnames:
                        period = sb9['Per'][0]
                        if not np.isnan(period) and period > 0:
                            params['period'] = float(period)  # jours
                    
                    # Excentricité
                    if 'e' in sb9.colnames:
                        ecc = sb9['e'][0]
                        if not np.isnan(ecc):
                            params['eccentricity'] = float(ecc)
                    
                    # Argument du périastre
                    if 'omega' in sb9.colnames:
                        omega = sb9['omega'][0]
                        if not np.isnan(omega):
                            params['omega'] = float(omega)  # degrés
                    
                    # Demi-grand axe
                    if 'a' in sb9.colnames:
                        a = sb9['a'][0]
                        if not np.isnan(a):
                            params['semi_major_axis'] = float(a)  # UA
                    
                    # Inclinaison (si disponible)
                    if 'i' in sb9.colnames:
                        inc = sb9['i'][0]
                        if not np.isnan(inc):
                            params['inclination'] = float(inc)  # degrés
                    
                    # Masses
                    if 'M1' in sb9.colnames:
                        m1 = sb9['M1'][0]
                        if not np.isnan(m1):
                            params['mass_primary'] = float(m1)  # masses solaires
                    
                    if 'M2' in sb9.colnames:
                        m2 = sb9['M2'][0]
                        if not np.isnan(m2):
                            params['mass_secondary'] = float(m2)  # masses solaires
                    
                    logger.info(f"SB9: Paramètres orbitaux trouvés")
                    return params
        except Exception as e:
            logger.debug(f"SB9 query error: {e}")
        return None
    
    def _query_sb9_by_coord(self, coord: SkyCoord, radius: float) -> Optional[Dict]:
        """Cherche dans SB9 par coordonnées"""
        return self._query_sb9("", coord)  # Même fonction, on ignore le nom
    
    def _query_gaia_binary(self, coord: SkyCoord) -> Optional[Dict]:
        """
        Interroge Gaia DR3 pour les paramètres de binaires
        Catalogue: I/355/gaiadr3
        """
        try:
            v = Vizier(columns=['**'], row_limit=1)
            result = v.query_region(coord, radius=0.05*u.deg, catalog='I/355/gaiadr3')
            
            if result and len(result) > 0:
                gaia = result[0]
                if len(gaia) > 0:
                    params = {}
                    
                    # Parallaxe (pour distance)
                    if 'Plx' in gaia.colnames:
                        plx = gaia['Plx'][0]
                        if not np.isnan(plx) and plx > 0:
                            distance_pc = 1000.0 / plx  # parsecs
                            params['distance'] = float(distance_pc)
                    
                    # Magnitude G
                    if 'Gmag' in gaia.colnames:
                        gmag = gaia['Gmag'][0]
                        if not np.isnan(gmag):
                            params['magnitude_g'] = float(gmag)
                    
                    logger.info("Gaia: Paramètres de base trouvés")
                    return params
        except Exception as e:
            logger.debug(f"Gaia query error: {e}")
        return None
    
    def _query_debcat(self, star_name: str, coord: SkyCoord) -> Optional[Dict]:
        """
        Interroge DEBCat (Detached Eclipsing Binaries Catalog)
        Catalogue: J/ApJS/232/23/debcat
        Contient des paramètres détaillés pour les binaires à éclipses détachées
        """
        try:
            v = Vizier(columns=['**'], row_limit=1)
            # DEBCat peut être dans différents catalogues Vizier
            catalogs = ['J/ApJS/232/23/debcat', 'J/A+A/608/A22']
            
            for cat in catalogs:
                try:
                    result = v.query_region(coord, radius=0.1*u.deg, catalog=cat)
                    
                    if result and len(result) > 0:
                        deb = result[0]
                        if len(deb) > 0:
                            params = {}
                            
                            # Identifier les colonnes disponibles
                            for col in deb.colnames:
                                col_lower = col.lower()
                                value = deb[col][0]
                                
                                if np.isnan(value) if isinstance(value, (int, float)) else False:
                                    continue
                                
                                # Période
                                if 'per' in col_lower or 'period' in col_lower:
                                    if value > 0:
                                        params['period'] = float(value)
                                
                                # Températures
                                if 'teff' in col_lower or 'temp' in col_lower:
                                    if '1' in col_lower or 'p' in col_lower:
                                        params['teff_primary'] = float(value)
                                    elif '2' in col_lower or 's' in col_lower:
                                        params['teff_secondary'] = float(value)
                                
                                # Rayons
                                if 'r' in col_lower and ('sol' in col_lower or 'sun' in col_lower):
                                    if '1' in col_lower or 'p' in col_lower:
                                        params['radius_primary'] = float(value)
                                    elif '2' in col_lower or 's' in col_lower:
                                        params['radius_secondary'] = float(value)
                                
                                # Masses
                                if 'm' in col_lower and ('sol' in col_lower or 'sun' in col_lower):
                                    if '1' in col_lower or 'p' in col_lower:
                                        params['mass_primary'] = float(value)
                                    elif '2' in col_lower or 's' in col_lower:
                                        params['mass_secondary'] = float(value)
                                
                                # Excentricité
                                if 'e' in col_lower and 'ec' in col_lower:
                                    params['eccentricity'] = float(value)
                                
                                # Inclinaison
                                if 'i' in col_lower and 'inc' in col_lower:
                                    params['inclination'] = float(value)
                            
                            if params:
                                logger.info("DEBCat: Paramètres détaillés trouvés")
                                return params
                except:
                    continue
        except Exception as e:
            logger.debug(f"DEBCat query error: {e}")
        return None
    
    def format_for_phoebe(self, params: Dict) -> Dict:
        """
        Formate les paramètres récupérés pour PHOEBE2
        
        Parameters
        ----------
        params : Dict
            Paramètres bruts des catalogues
        
        Returns
        -------
        Dict
            Paramètres formatés pour PHOEBE
        """
        phoebe_params = {}
        
        # Paramètres orbitaux
        if 'period' in params:
            phoebe_params['period'] = params['period']  # jours
        
        if 'eccentricity' in params:
            phoebe_params['eccentricity'] = params['eccentricity']
        elif 'e' in params:
            phoebe_params['eccentricity'] = params['e']
        else:
            phoebe_params['eccentricity'] = 0.0  # Par défaut circulaire
        
        if 'omega' in params:
            phoebe_params['argument_of_periastron'] = params['omega']  # degrés
        elif 'argument_of_periastron' in params:
            phoebe_params['argument_of_periastron'] = params['argument_of_periastron']
        
        if 'inclination' in params:
            phoebe_params['inclination'] = params['inclination']  # degrés
        elif 'i' in params:
            phoebe_params['inclination'] = params['i']
        else:
            phoebe_params['inclination'] = 90.0  # Par défaut éclipsant
        
        if 'semi_major_axis' in params:
            phoebe_params['semi_major_axis'] = params['semi_major_axis']  # UA ou R☉
        
        # Paramètres des composantes
        if 'mass_primary' in params:
            phoebe_params['mass_primary'] = params['mass_primary']  # M☉
        
        if 'mass_secondary' in params:
            phoebe_params['mass_secondary'] = params['mass_secondary']  # M☉
        
        if 'radius_primary' in params:
            phoebe_params['radius_primary'] = params['radius_primary']  # R☉
        
        if 'radius_secondary' in params:
            phoebe_params['radius_secondary'] = params['radius_secondary']  # R☉
        
        if 'teff_primary' in params:
            phoebe_params['teff_primary'] = params['teff_primary']  # K
        
        if 'teff_secondary' in params:
            phoebe_params['teff_secondary'] = params['teff_secondary']  # K
        
        # Distance
        if 'distance' in params:
            phoebe_params['distance'] = params['distance']  # parsecs
        
        return phoebe_params

