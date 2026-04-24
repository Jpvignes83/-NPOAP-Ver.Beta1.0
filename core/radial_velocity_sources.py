"""
Module pour récupérer les vitesses radiales et autres paramètres stellaires
depuis diverses bases de données astronomiques.
"""

import logging
import pandas as pd
import numpy as np
from typing import Optional, Dict, List
from pathlib import Path

logger = logging.getLogger(__name__)

# Vérifier les dépendances optionnelles
try:
    from astroquery.simbad import Simbad
    from astropy.coordinates import SkyCoord
    SIMBAD_AVAILABLE = True
except ImportError:
    SIMBAD_AVAILABLE = False
    logger.warning("astroquery.simbad non disponible. Les requêtes SIMBAD seront désactivées.")

try:
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    NASA_AVAILABLE = True
except ImportError:
    NASA_AVAILABLE = False
    logger.warning("astroquery NASA non disponible.")

try:
    import pyvo as vo
    PYVO_AVAILABLE = True
except ImportError:
    PYVO_AVAILABLE = False
    logger.warning("pyvo non disponible. Les requêtes ESO TAP seront désactivées.")


def query_simbad_coordinates(target: str) -> tuple:
    """
    Interroge SIMBAD pour obtenir les coordonnées et alias d'une cible.
    
    Parameters
    ----------
    target : str
        Nom de la cible
        
    Returns
    -------
    tuple
        (target_alias, ra_deg, dec_deg) ou (None, None, None) si échec
    """
    if not SIMBAD_AVAILABLE:
        logger.warning("SIMBAD non disponible")
        return None, None, None
    
    try:
        pos = SkyCoord.from_name(target)
        target_alias = Simbad.query_objectids(target)
        return target_alias, pos.ra.degree, pos.dec.degree
    except Exception as e:
        logger.warning(f"Erreur lors de la requête SIMBAD pour {target}: {e}")
        return None, None, None


def query_nasa_radial_velocity(target: str) -> float:
    """
    Récupère la vitesse radiale systémique (gamma) depuis NASA Exoplanet Archive.
    
    Parameters
    ----------
    target : str
        Nom de la cible (hostname)
        
    Returns
    -------
    float
        Vitesse radiale systémique en km/s, ou 0.0 si non trouvée
    """
    if not NASA_AVAILABLE:
        logger.warning("NASA Exoplanet Archive non disponible")
        return 0.0
    
    try:
        # Essayer plusieurs variantes du nom
        name_variants = [
            target,
            target.upper(),
            target.replace(' ', '-'),
            target.replace('-', ' ')
        ]
        
        for variant in name_variants:
            try:
                result_table = NasaExoplanetArchive.query_criteria(
                    table="pscomppars",
                    where=f"hostname like '{variant}%'"
                )
                
                if len(result_table) > 0:
                    st_radv = result_table[0]["st_radv"]
                    # Extraire la valeur numérique si c'est un Quantity
                    if hasattr(st_radv, 'value'):
                        gamma = float(st_radv.value)
                    else:
                        gamma = float(str(st_radv).split()[0])
                    
                    logger.info(f"Gamma récupéré depuis NASA Exoplanet Archive pour {target}: {gamma} km/s")
                    return gamma
            except Exception:
                continue
        
        logger.warning(f"Gamma non trouvé pour {target} dans NASA Exoplanet Archive")
        return 0.0
        
    except Exception as e:
        logger.warning(f"Erreur lors de la récupération du gamma pour {target}: {e}")
        return 0.0


def query_eso_tap_harps(target: str, search_radius_arcmin: float = 2.5) -> pd.DataFrame:
    """
    Interroge ESO TAP pour récupérer les données HARPS.
    
    Parameters
    ----------
    target : str
        Nom de la cible
    search_radius_arcmin : float
        Rayon de recherche en arcminutes (défaut: 2.5)
        
    Returns
    -------
    pd.DataFrame
        DataFrame avec colonnes: bjd, rv, rv_err, ins_name
    """
    if not PYVO_AVAILABLE:
        logger.warning("pyvo non disponible pour les requêtes ESO TAP")
        return pd.DataFrame()
    
    if not SIMBAD_AVAILABLE:
        logger.warning("SIMBAD non disponible pour obtenir les coordonnées")
        return pd.DataFrame()
    
    try:
        ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
        tapobs = vo.dal.TAPService(ESO_TAP_OBS)
        
        # Obtenir les coordonnées depuis SIMBAD
        target_alias, ra, dec = query_simbad_coordinates(target)
        if ra is None or dec is None:
            logger.warning(f"Impossible d'obtenir les coordonnées pour {target}")
            return pd.DataFrame()
        
        # Rayon de recherche en degrés
        sr = search_radius_arcmin / 60.0
        
        # Requête SQL pour trouver les observations
        query = f"""SELECT *
                   FROM ivoa.ObsCore
                   WHERE intersects(s_region, circle('', {ra}, {dec}, {sr}))=1"""
        
        res = tapobs.search(query=query)
        
        if len(res) == 0:
            logger.info(f"Aucune observation HARPS trouvée pour {target}")
            return pd.DataFrame()
        
        product_id_list = tuple(res['dp_id'])
        logger.info(f"Trouvé {len(product_id_list)} produits HARPS pour {target}")
        
        # Note: L'extraction complète des données HARPS nécessite de télécharger
        # les fichiers FITS depuis l'archive ESO, ce qui est complexe.
        # Pour l'instant, on retourne un DataFrame vide avec la structure attendue.
        # Une implémentation complète nécessiterait requests, tarfile, et astropy.io.fits
        
        logger.warning("Extraction complète des données HARPS non implémentée (nécessite téléchargement de fichiers)")
        return pd.DataFrame(columns=['bjd', 'rv', 'rv_err', 'ins_name'])
        
    except Exception as e:
        logger.warning(f"Erreur lors de la requête ESO TAP HARPS pour {target}: {e}")
        return pd.DataFrame()


def query_dace(target: str) -> pd.DataFrame:
    """
    Interroge DACE (Data Analysis Center for Exoplanets) pour les vitesses radiales.
    
    Parameters
    ----------
    target : str
        Nom de la cible
        
    Returns
    -------
    pd.DataFrame
        DataFrame avec colonnes: bjd, rv, rv_err, ins_name
    """
    try:
        from dace_query.spectroscopy import Spectroscopy
        
        # Obtenir les alias depuis SIMBAD
        target_alias, ra, dec = query_simbad_coordinates(target)
        if target_alias is None:
            logger.warning(f"Impossible d'obtenir les alias SIMBAD pour {target}")
            return pd.DataFrame()
        
        aliases = [(alias[0].replace(" ", ""), ra, dec) for alias in target_alias]
        
        for simbad_alias, ra, dec in aliases:
            try:
                radial_velocities_table = Spectroscopy.get_timeseries(
                    simbad_alias,
                    sorted_by_instrument=False,
                    output_format="pandas"
                )
                
                if not radial_velocities_table.empty:
                    break
            except Exception:
                continue
        else:
            logger.info(f"Aucune donnée de vitesse radiale trouvée pour {target} dans DACE")
            return pd.DataFrame()
        
        # Filtrer les données valides (erreur > 0)
        valid_data = radial_velocities_table[radial_velocities_table['rv_err'] > 0]
        
        if valid_data.empty:
            logger.info(f"Aucune donnée de vitesse radiale valide pour {target} dans DACE")
            return pd.DataFrame()
        
        # Convertir RJD en BJD
        bjd = valid_data['rjd'] + 2400000.0
        rv = valid_data['rv']
        rv_err = valid_data['rv_err']
        ins_name = valid_data['ins_name']
        
        return pd.DataFrame({
            'bjd': bjd,
            'rv': rv,
            'rv_err': rv_err,
            'ins_name': ins_name
        })
        
    except ImportError:
        logger.warning("dace_query non disponible. Installation: pip install dace-query")
        return pd.DataFrame()
    except Exception as e:
        logger.warning(f"Erreur lors de la requête DACE pour {target}: {e}")
        return pd.DataFrame()


def get_radial_velocity_data(target: str, sources: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Récupère les données de vitesses radiales depuis plusieurs sources.
    
    Parameters
    ----------
    target : str
        Nom de la cible
    sources : list, optional
        Liste des sources à interroger. Si None, utilise toutes les sources disponibles.
        Sources possibles: 'dace', 'nasa_gamma', 'eso_harps'
        
    Returns
    -------
    pd.DataFrame
        DataFrame avec colonnes: bjd, rv, rv_err, ins_name
    """
    if sources is None:
        sources = ['dace', 'nasa_gamma']
        if PYVO_AVAILABLE:
            sources.append('eso_harps')
    
    dataframes = []
    
    for source in sources:
        try:
            if source == 'dace':
                df = query_dace(target)
                if not df.empty:
                    dataframes.append(df)
            elif source == 'nasa_gamma':
                # Gamma est une valeur unique, pas une série temporelle
                # On pourrait l'ajouter comme métadonnée
                gamma = query_nasa_radial_velocity(target)
                if gamma != 0.0:
                    logger.info(f"Gamma (vitesse systémique) pour {target}: {gamma} km/s")
            elif source == 'eso_harps':
                df = query_eso_tap_harps(target)
                if not df.empty:
                    dataframes.append(df)
        except Exception as e:
            logger.warning(f"Erreur lors de la requête {source} pour {target}: {e}")
            continue
    
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.sort_values('bjd', inplace=True)
        combined_df.drop_duplicates(subset=['bjd', 'rv', 'rv_err', 'ins_name'], inplace=True)
        logger.info(f"Total de {len(combined_df)} observations de vitesse radiale récupérées pour {target}")
        return combined_df
    else:
        logger.warning(f"Aucune donnée de vitesse radiale trouvée pour {target}")
        return pd.DataFrame(columns=['bjd', 'rv', 'rv_err', 'ins_name'])


def get_stellar_parameters(target: str) -> Dict:
    """
    Récupère les paramètres stellaires depuis diverses sources.
    
    Parameters
    ----------
    target : str
        Nom de la cible
        
    Returns
    -------
    dict
        Dictionnaire contenant les paramètres stellaires disponibles
    """
    params = {}
    
    # Récupérer gamma depuis NASA
    if NASA_AVAILABLE:
        gamma = query_nasa_radial_velocity(target)
        if gamma != 0.0:
            params['gamma'] = gamma  # km/s
            params['gamma_source'] = 'NASA Exoplanet Archive'
    
    # Récupérer les coordonnées depuis SIMBAD
    if SIMBAD_AVAILABLE:
        target_alias, ra, dec = query_simbad_coordinates(target)
        if ra is not None and dec is not None:
            params['ra'] = ra  # degrés
            params['dec'] = dec  # degrés
            params['coordinates_source'] = 'SIMBAD'
            if target_alias is not None:
                params['aliases'] = [alias[0] for alias in target_alias]
    
    return params
