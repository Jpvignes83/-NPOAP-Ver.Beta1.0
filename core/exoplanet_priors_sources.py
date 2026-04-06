"""
Module pour récupérer les priors (a/R* et i) depuis plusieurs sources alternatives
pour l'ajustement des coefficients de limb-darkening avec priors.

Sources supportées :
1. NASA Exoplanet Archive (via astroquery)
2. Extrasolar Planets Encyclopaedia (via API web)
3. Open Exoplanet Catalogue (via fichier XML)
4. AAVSO Exoplanet Database (via API)
"""

import numpy as np
import logging
from typing import Optional, Dict, Tuple
from urllib.request import urlopen, Request
from urllib.parse import urlencode, quote
from urllib.error import URLError
import json
import xml.etree.ElementTree as ET
from io import StringIO
import time
from pathlib import Path

logger = logging.getLogger(__name__)

# Fichier de cache des priors (une fois pour toutes, pas de réinterrogation)
_PRIORS_CACHE_FILENAME = "priors_cache.json"


def _priors_cache_path() -> Path:
    """Chemin du fichier cache des priors (à la racine du projet NPOAP)."""
    return Path(__file__).resolve().parent.parent / _PRIORS_CACHE_FILENAME


def _normalize_planet_key(name: str) -> str:
    """Clé unique pour le cache : nom normalisé (espaces réduits, strip)."""
    return " ".join(str(name).strip().split()) if name else ""


def _load_priors_cache() -> Dict[str, dict]:
    """Charge le cache des priors depuis le fichier JSON. Retourne {} si absent ou invalide."""
    path = _priors_cache_path()
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except Exception as e:
        logger.debug("Cache priors: lecture impossible (%s): %s", path, e)
        return {}


def _save_priors_cache(cache: Dict[str, dict]) -> None:
    """Sauvegarde le cache des priors dans le fichier JSON."""
    path = _priors_cache_path()
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2, ensure_ascii=False)
        logger.debug("Cache priors enregistré: %s", path)
    except Exception as e:
        logger.warning("Impossible d'enregistrer le cache priors (%s): %s", path, e)

try:
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    NASA_AVAILABLE = True
except ImportError:
    NASA_AVAILABLE = False
    logger.warning("astroquery non disponible, NASA Exoplanet Archive désactivé")


def get_priors_from_nasa(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors (a/R* et i) depuis NASA Exoplanet Archive.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète (ex: "HD 209458 b", "WASP-12 b", "wasp-19 b")
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination' en degrés, ou None si non trouvé
    """
    if not NASA_AVAILABLE:
        return None
    
    # Normaliser le nom : essayer plusieurs variantes
    # Format standard NASA: "WASP-19 b" (majuscules, tiret, espace avant 'b')
    name_variants = []
    
    # 1. Nom original
    name_variants.append(target_name)
    
    # 2. Format standard NASA (majuscules + tiret) - PRIORITÉ
    # Convertir "wasp-19 b" ou "WASP 19 b" en "WASP-19 b"
    nasa_format = target_name.upper().replace(' ', '-')
    if nasa_format != target_name:
        name_variants.append(nasa_format)
    
    # 3. Autres variantes
    name_variants.extend([
        target_name.replace(' ', '-'),  # Espaces -> tirets
        target_name.replace('-', ' '),  # Tirets -> espaces
        target_name.upper(),  # Majuscules
        target_name.lower(),  # Minuscules
        target_name.replace(' ', ''),  # Sans espaces
    ])
    
    # Supprimer les doublons tout en préservant l'ordre
    seen = set()
    unique_variants = []
    for variant in name_variants:
        if variant not in seen:
            seen.add(variant)
            unique_variants.append(variant)
    name_variants = unique_variants
    
    try:
        # Essayer d'abord la table pscomppars (meilleures valeurs)
        table = None
        used_name = None
        for variant in name_variants:
            try:
                table = NasaExoplanetArchive.query_object(variant, table="pscomppars")
                if len(table) > 0:
                    used_name = variant
                    logger.debug(f"NASA: Données trouvées dans pscomppars avec le nom: {variant}")
                    break
            except Exception:
                continue
        
        # Si pas trouvé dans pscomppars, essayer exoplanets
        if table is None or len(table) == 0:
            for variant in name_variants:
                try:
                    table = NasaExoplanetArchive.query_object(variant, table="exoplanets")
                    if len(table) > 0:
                        used_name = variant
                        logger.debug(f"NASA: Données trouvées dans exoplanets avec le nom: {variant}")
                        break
                except Exception:
                    continue
        
        if table is None or len(table) == 0:
            logger.debug(f"NASA: Aucune donnée trouvée pour {target_name} (variantes essayées: {len(name_variants)})")
            return None
        
        row = table[0]
        logger.info(f"NASA: Données trouvées! Colonnes disponibles: {', '.join(row.colnames[:30])}")
        
        # Récupérer a/R* et inclinaison
        # Colonnes possibles selon la table
        a_rs = None
        if 'pl_ratdor' in row.colnames:
            if not np.ma.is_masked(row['pl_ratdor']):
                a_rs = float(row['pl_ratdor'])
        elif 'pl_orbsmax' in row.colnames and 'st_rad' in row.colnames:
            # Calculer a/R* si disponible
            if (not np.ma.is_masked(row['pl_orbsmax']) and 
                not np.ma.is_masked(row['st_rad']) and 
                float(row['st_rad']) > 0):
                a_rs = float(row['pl_orbsmax']) / float(row['st_rad'])
        
        inclination = None
        if 'pl_orbincl' in row.colnames:
            if not np.ma.is_masked(row['pl_orbincl']):
                inclination = float(row['pl_orbincl'])  # En degrés
                logger.debug(f"NASA: Inclinaison trouvée dans pl_orbincl: {inclination:.4f}°")
        
        # Si a/R* ou inclinaison manquent, essayer de les calculer ou chercher dans d'autres colonnes
        if a_rs is None:
            # Essayer de calculer a/R* depuis a (UA) et R* (R☉)
            if 'pl_orbsmax' in row.colnames and 'st_rad' in row.colnames:
                if (not np.ma.is_masked(row['pl_orbsmax']) and 
                    not np.ma.is_masked(row['st_rad']) and 
                    float(row['st_rad']) > 0):
                    # Convertir a de UA en R☉ : 1 UA = 215 R☉ (approximatif)
                    a_ua = float(row['pl_orbsmax'])
                    r_star = float(row['st_rad'])
                    a_rs = (a_ua * 215.0) / r_star
                    logger.debug(f"NASA: a/R* calculé depuis a={a_ua:.4f} UA et R*={r_star:.4f} R☉")
        
        if inclination is None:
            # Chercher dans d'autres colonnes possibles
            for col in ['pl_orbincl', 'pl_incl', 'inclination', 'i']:
                if col in row.colnames and not np.ma.is_masked(row[col]):
                    inclination = float(row[col])
                    logger.debug(f"NASA: Inclinaison trouvée dans la colonne {col}")
                    break
        
        if a_rs is None or inclination is None:
            logger.warning(f"NASA: Paramètres incomplets pour {used_name or target_name} (a/R*={a_rs}, i={inclination})")
            # Log les valeurs brutes pour debug
            if 'pl_ratdor' in row.colnames:
                val = row['pl_ratdor']
                logger.info(f"  → pl_ratdor (brut): {val} (masked: {np.ma.is_masked(val)})")
            if 'pl_orbsmax' in row.colnames:
                val = row['pl_orbsmax']
                logger.info(f"  → pl_orbsmax (brut): {val} (masked: {np.ma.is_masked(val)})")
            if 'st_rad' in row.colnames:
                val = row['st_rad']
                logger.info(f"  → st_rad (brut): {val} (masked: {np.ma.is_masked(val)})")
            if 'pl_orbincl' in row.colnames:
                val = row['pl_orbincl']
                logger.info(f"  → pl_orbincl (brut): {val} (masked: {np.ma.is_masked(val)})")
            return None
        
        logger.info(f"NASA: Priors récupérés pour {used_name or target_name}: a/R*={a_rs:.4f}, i={inclination:.4f}°")
        return {'a_rs': a_rs, 'inclination': inclination, 'source': 'NASA Exoplanet Archive'}
    
    except Exception as e:
        logger.debug(f"NASA: Erreur pour {target_name}: {e}")
        return None


def get_priors_from_nasa_tap(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis NASA Exoplanet Archive via l'API TAP (Table Access Protocol).
    Méthode standardisée et recommandée par la NASA.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs', 'inclination', et leurs erreurs si disponibles
    """
    try:
        import pandas as pd
        from urllib.request import urlopen, Request
        from urllib.parse import urlencode
        from io import StringIO
        
        # Normaliser le nom : essayer plusieurs variantes
        name_variants = []
        name_variants.append(target_name)
        nasa_format = target_name.upper().replace(' ', '-')
        if nasa_format != target_name:
            name_variants.append(nasa_format)
        name_variants.extend([
            target_name.replace(' ', '-'),
            target_name.replace('-', ' '),
            target_name.upper(),
            target_name.lower(),
            target_name.replace(' ', ''),
        ])
        
        # Supprimer les doublons
        seen = set()
        unique_variants = []
        for variant in name_variants:
            if variant not in seen:
                seen.add(variant)
                unique_variants.append(variant)
        name_variants = unique_variants
        
        # Utiliser astroquery au lieu de TAP direct (évite les erreurs HTTP 400)
        try:
            from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
            import numpy as np
            
            # Fonction helper pour extraire la valeur
            def extract_value(val):
                if pd.notna(val) and not (isinstance(val, float) and np.isnan(val)):
                    try:
                        if hasattr(val, 'value'):
                            return float(val.value)
                        return float(val)
                    except (ValueError, TypeError):
                        pass
                return None
            
            # Essayer pscomppars d'abord (meilleures valeurs)
            for variant in name_variants:
                try:
                    table = NasaExoplanetArchive.query_object(variant, table="pscomppars")
                    if len(table) > 0:
                        row = table[0]
                        
                        # Récupérer a/R*
                        a_rs = None
                        a_rs_err = None
                        if 'pl_ratdor' in row.colnames:
                            a_rs = extract_value(row['pl_ratdor'])
                            if a_rs is not None and 'pl_ratdor_err' in row.colnames:
                                a_rs_err = extract_value(row['pl_ratdor_err'])
                        elif 'pl_orbsmax' in row.colnames and 'st_rad' in row.colnames:
                            # Calculer a/R* depuis a (UA) et R* (R☉)
                            a_ua = extract_value(row['pl_orbsmax'])
                            r_star = extract_value(row['st_rad'])
                            if a_ua is not None and r_star is not None and r_star > 0:
                                # Convertir a de UA en R☉ : 1 UA = 215 R☉ (approximatif)
                                a_rs = (a_ua * 215.0) / r_star
                                logger.debug(f"NASA (astroquery): a/R* calculé depuis a={a_ua:.4f} UA et R*={r_star:.4f} R☉")
                        
                        # Récupérer l'inclinaison
                        inclination = None
                        inclination_err = None
                        if 'pl_orbincl' in row.colnames:
                            inclination = extract_value(row['pl_orbincl'])
                            if inclination is not None and 'pl_orbincl_err' in row.colnames:
                                inclination_err = extract_value(row['pl_orbincl_err'])
                        
                        if a_rs is not None and inclination is not None:
                            result = {
                                'a_rs': a_rs,
                                'inclination': inclination,
                                'source': 'NASA Exoplanet Archive (astroquery)'
                            }
                            if a_rs_err is not None:
                                result['a_rs_err'] = a_rs_err
                            if inclination_err is not None:
                                result['inclination_err'] = inclination_err
                            
                            logger.info(f"✅ NASA (astroquery pscomppars): Priors récupérés pour {variant}: a/R*={a_rs:.4f}, i={inclination:.4f}°")
                            return result
                except Exception as e:
                    logger.debug(f"Erreur astroquery pscomppars pour {variant}: {e}")
                    continue
        
            # Si pas trouvé dans pscomppars, essayer exoplanets
            for variant in name_variants:
                try:
                    table = NasaExoplanetArchive.query_object(variant, table="exoplanets")
                    if len(table) > 0:
                        row = table[0]
                        
                        # Même logique d'extraction
                        a_rs = None
                        if 'pl_ratdor' in row.colnames:
                            a_rs = extract_value(row['pl_ratdor'])
                        elif 'pl_orbsmax' in row.colnames and 'st_rad' in row.colnames:
                            a_ua = extract_value(row['pl_orbsmax'])
                            r_star = extract_value(row['st_rad'])
                            if a_ua is not None and r_star is not None and r_star > 0:
                                a_rs = (a_ua * 215.0) / r_star
                        
                        inclination = None
                        if 'pl_orbincl' in row.colnames:
                            inclination = extract_value(row['pl_orbincl'])
                        
                        if a_rs is not None and inclination is not None:
                            result = {
                                'a_rs': a_rs,
                                'inclination': inclination,
                                'source': 'NASA Exoplanet Archive (astroquery exoplanets)'
                            }
                            logger.info(f"✅ NASA (astroquery exoplanets): Priors récupérés pour {variant}: a/R*={a_rs:.4f}, i={inclination:.4f}°")
                            return result
                except Exception:
                    continue
            
            logger.debug("NASA TAP (astroquery): Aucune donnée trouvée avec astroquery")
            return None
                    
        except ImportError:
            logger.debug("astroquery non disponible pour nasa_tap")
            return None
        except Exception as e:
            logger.warning(f"Erreur lors de la récupération des priors via astroquery: {e}")
            return None
    
    except Exception as e:
        logger.debug(f"NASA TAP: Erreur générale pour {target_name}: {e}")
        return None


def get_priors_from_nasa_api_direct(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis NASA Exoplanet Archive via l'API REST directe.
    Alternative à astroquery pour plus de flexibilité.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination', ou None si non trouvé
    """
    try:
        import pandas as pd
        
        # Normaliser le nom
        name_variants = [
            target_name,
            target_name.upper().replace(' ', '-'),
            target_name.replace(' ', '-'),
            target_name.replace('-', ' '),
        ]
        name_variants = list(dict.fromkeys(name_variants))
        
        base_url = "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI"
        
        for variant in name_variants:
            try:
                # Essayer la table pscomppars avec wildcards
                params = {
                    'table': 'pscomppars',
                    'format': 'csv',
                    'select': 'pl_name,pl_ratdor,pl_orbincl,pl_orbsmax,st_rad',
                    'where': f"pl_name like '%{variant}%'"
                }
                url = f"{base_url}?{urlencode(params)}"
                
                req = Request(url)
                with urlopen(req, timeout=10) as response:
                    data = response.read().decode('utf-8')
                
                if not data or data.strip() == '' or 'pl_name' not in data:
                    continue
                
                df = pd.read_csv(StringIO(data))
                if len(df) == 0:
                    continue
                
                row = df.iloc[0]
                
                # Récupérer a/R* et inclinaison
                a_rs = None
                if pd.notna(row.get('pl_ratdor')):
                    a_rs = float(row['pl_ratdor'])
                elif pd.notna(row.get('pl_orbsmax')) and pd.notna(row.get('st_rad')):
                    if float(row['st_rad']) > 0:
                        a_ua = float(row['pl_orbsmax'])
                        r_star = float(row['st_rad'])
                        a_rs = (a_ua * 215.0) / r_star
                
                inclination = None
                if pd.notna(row.get('pl_orbincl')):
                    inclination = float(row['pl_orbincl'])
                
                if a_rs is not None and inclination is not None:
                    logger.info(f"✅ NASA API: Priors récupérés pour {variant}: a/R*={a_rs:.4f}, i={inclination:.4f}°")
                    return {'a_rs': a_rs, 'inclination': inclination, 'source': 'NASA Exoplanet Archive (API)'}
                else:
                    logger.warning(f"NASA API: Paramètres incomplets pour {variant} (a/R*={a_rs}, i={inclination})")
                    # Log les valeurs brutes
                    logger.info(f"  → pl_ratdor: {row.get('pl_ratdor', 'N/A')}")
                    logger.info(f"  → pl_orbsmax: {row.get('pl_orbsmax', 'N/A')}")
                    logger.info(f"  → st_rad: {row.get('st_rad', 'N/A')}")
                    logger.info(f"  → pl_orbincl: {row.get('pl_orbincl', 'N/A')}")
                
            except Exception as e:
                logger.warning(f"NASA API: Erreur pour {variant}: {e}")
                continue
        
        return None
    
    except Exception as e:
        logger.debug(f"NASA API: Erreur générale pour {target_name}: {e}")
        return None


def get_priors_from_exoplanet_eu(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis The Extrasolar Planets Encyclopaedia.
    
    Note: Cette source n'a pas d'API publique directe, mais on peut essayer
    de parser les pages web ou utiliser des données agrégées.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination', ou None si non trouvé
    """
    try:
        # L'Encyclopédie n'a pas d'API publique, mais on peut essayer de parser
        # les données depuis leur catalogue CSV ou JSON s'ils sont disponibles
        # Pour l'instant, on retourne None car il faudrait scraper les pages web
        # ce qui n'est pas idéal
        
        # TODO: Implémenter si une API ou un endpoint CSV/JSON devient disponible
        return None
    
    except Exception as e:
        logger.debug(f"Exoplanet.eu: Erreur pour {target_name}: {e}")
        return None


def get_priors_from_open_exoplanet_catalogue(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis Open Exoplanet Catalogue.
    
    Open Exoplanet Catalogue fournit un fichier XML avec toutes les exoplanètes.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination', ou None si non trouvé
    """
    try:
        # URL du fichier XML complet
        url = "https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_gzip/master/systems.xml.gz"
        
        # Pour l'instant, on ne télécharge pas le fichier complet car il est volumineux
        # On pourrait implémenter un cache local ou une recherche par API si disponible
        
        # TODO: Implémenter le parsing du XML si nécessaire
        # Le fichier est en format XML avec une structure hiérarchique
        return None
    
    except Exception as e:
        logger.debug(f"Open Exoplanet Catalogue: Erreur pour {target_name}: {e}")
        return None


def get_priors_from_aavso(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis AAVSO Exoplanet Database.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination', ou None si non trouvé
    """
    try:
        # AAVSO a une API mais elle est principalement orientée vers les observations
        # Les paramètres orbitaux peuvent être disponibles mais nécessitent un parsing spécifique
        
        # URL de l'API AAVSO (exemple, à vérifier)
        # base_url = "https://archive.aavso.org/exoplanet-section"
        
        # Pour l'instant, on retourne None car l'API exacte doit être vérifiée
        # TODO: Implémenter si l'API AAVSO fournit ces paramètres
        return None
    
    except Exception as e:
        logger.debug(f"AAVSO: Erreur pour {target_name}: {e}")
        return None


def get_priors_from_simbad(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis SIMBAD (via astroquery si disponible).
    
    SIMBAD peut contenir des informations sur les systèmes exoplanétaires,
    mais les paramètres orbitaux détaillés peuvent ne pas être disponibles.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète ou de l'étoile hôte
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination', ou None si non trouvé
    """
    try:
        # SIMBAD via astroquery
        try:
            from astroquery.simbad import Simbad
            from astroquery.exceptions import TableParseError
            
            # SIMBAD contient principalement des informations sur les étoiles
            # Les paramètres orbitaux des planètes ne sont généralement pas dans SIMBAD
            # On retourne None pour l'instant
            
            # TODO: Vérifier si SIMBAD peut fournir ces informations
            return None
            
        except ImportError:
            logger.debug("astroquery.simbad non disponible")
            return None
    
    except Exception as e:
        logger.debug(f"SIMBAD: Erreur pour {target_name}: {e}")
        return None


def get_priors_from_pylightcurve(target_name: str) -> Optional[Dict[str, float]]:
    """
    Récupère les priors depuis pylightcurve (base de données locale).
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant 'a_rs' et 'inclination', ou None si non trouvé
    """
    try:
        import pylightcurve as plc
        
        # Normaliser le nom : essayer plusieurs variantes
        name_variants = [
            target_name,  # Nom original
            target_name.lower(),  # Minuscules
            target_name.upper(),  # Majuscules
            target_name.replace(' ', '-'),  # Espaces -> tirets
            target_name.replace('-', ' '),  # Tirets -> espaces
            target_name.replace(' ', ''),  # Sans espaces
        ]
        # Supprimer les doublons
        name_variants = list(dict.fromkeys(name_variants))
        
        pl = None
        used_name = None
        for variant in name_variants:
            try:
                pl = plc.get_planet(variant)
                if pl is not None:
                    used_name = variant
                    logger.debug(f"pylightcurve: Planète trouvée avec le nom: {variant}")
                    break
            except Exception:
                continue
        
        if pl is None:
            logger.debug(f"pylightcurve: Aucune planète trouvée pour {target_name}")
            return None
        
        # Récupérer a/R* et inclinaison
        a_rs = None
        if hasattr(pl, 'a_over_rs') and pl.a_over_rs > 0:
            a_rs = float(pl.a_over_rs)
        elif hasattr(pl, 'semi_major_axis') and hasattr(pl, 'stellar_radius'):
            if pl.semi_major_axis > 0 and pl.stellar_radius > 0:
                a_rs = float(pl.semi_major_axis) / float(pl.stellar_radius)
        
        inclination = None
        if hasattr(pl, 'inclination') and pl.inclination > 0:
            inclination = float(pl.inclination)
        elif hasattr(pl, 'orbital_inclination') and pl.orbital_inclination > 0:
            inclination = float(pl.orbital_inclination)
        
        if a_rs is None or inclination is None:
            logger.debug(f"pylightcurve: Paramètres incomplets pour {used_name or target_name}")
            return None
        
        logger.info(f"pylightcurve: Priors récupérés pour {used_name or target_name}: a/R*={a_rs:.4f}, i={inclination:.4f}°")
        return {'a_rs': a_rs, 'inclination': inclination, 'source': 'pylightcurve'}
    
    except Exception as e:
        logger.debug(f"pylightcurve: Erreur pour {target_name}: {e}")
        return None


def get_priors_from_all_sources(target_name: str,
                                sources_order: Optional[list] = None,
                                use_cache: bool = True,
                                force_refresh: bool = False) -> Optional[Dict[str, float]]:
    """
    Récupère les priors (a/R* et i) depuis le cache ou plusieurs sources en cascade.
    
    Les priors trouvés sont enregistrés dans priors_cache.json pour ne pas
    réinterroger les bases aux prochains appels.
    
    Parameters
    ----------
    target_name : str
        Nom de l'exoplanète
    sources_order : list, optional
        Ordre des sources à essayer (voir défaut dans le code)
    use_cache : bool
        Utiliser le cache disque (défaut True)
    force_refresh : bool
        Si True, ignorer le cache et réinterroger les sources (défaut False)
    
    Returns
    -------
    dict ou None
        Dictionnaire contenant:
        - 'a_rs': float, demi-grand axe en unités de rayon stellaire
        - 'inclination': float, inclinaison orbitale en degrés
        - 'source': str, nom de la source utilisée
        - 'a_rs_err': float, optional, erreur sur a/R*
        - 'inclination_err': float, optional, erreur sur l'inclinaison
        ou None si aucune source n'a fourni de données
    """
    key = _normalize_planet_key(target_name)
    if not key:
        logger.warning("Nom de planète vide pour le cache priors")
        return None

    # Utiliser le cache si demandé et pas de rafraîchissement forcé
    if use_cache and not force_refresh:
        cache = _load_priors_cache()
        if key in cache:
            priors = cache[key]
            if isinstance(priors, dict) and priors.get("a_rs") is not None and priors.get("inclination") is not None:
                logger.info("Priors chargés depuis le cache pour %s (a/R*=%.4f, i=%.4f°, source=%s)",
                            target_name, priors.get("a_rs"), priors.get("inclination"), priors.get("source", "?"))
                return priors

    if sources_order is None:
        sources_order = [
            'pylightcurve',      # Rapide, local
            'nasa_tap',          # NASA via TAP (recommandé, standardisé)
            'nasa',              # Fiable, complet (astroquery)
            'nasa_api_direct',   # Alternative NASA (API REST)
            'exoplanet_eu',      # Alternative
            'aavso',             # Observations
            'simbad'             # Dernier recours
        ]
    
    source_functions = {
        'pylightcurve': get_priors_from_pylightcurve,
        'nasa_tap': get_priors_from_nasa_tap,
        'nasa': get_priors_from_nasa,
        'nasa_api_direct': get_priors_from_nasa_api_direct,
        'exoplanet_eu': get_priors_from_exoplanet_eu,
        'open_exoplanet': get_priors_from_open_exoplanet_catalogue,
        'aavso': get_priors_from_aavso,
        'simbad': get_priors_from_simbad
    }
    
    for source_name in sources_order:
        if source_name not in source_functions:
            logger.warning(f"Source inconnue: {source_name}")
            continue
        
        func = source_functions[source_name]
        logger.info(f"🔍 Essai de la source '{source_name}' pour {target_name}...")
        try:
            result = func(target_name)
            if result is not None:
                logger.info(f"  → Source {source_name} a retourné un résultat")
                if 'a_rs' in result and 'inclination' in result:
                    # Vérifier que les valeurs sont valides
                    a_rs_val = result.get('a_rs')
                    incl_val = result.get('inclination')
                    logger.info(f"  → Valeurs trouvées: a/R*={a_rs_val}, i={incl_val}°")
                    if a_rs_val is not None and incl_val is not None:
                        if a_rs_val > 0 and 0 < incl_val <= 90:
                            logger.info(f"✅ Priors récupérés depuis {result.get('source', source_name)} pour {target_name}: a/R*={a_rs_val:.4f}, i={incl_val:.4f}°")
                            if use_cache:
                                cache = _load_priors_cache()
                                safe = {}
                                for k, v in result.items():
                                    if v is None:
                                        safe[k] = None
                                    elif isinstance(v, (np.floating, np.integer)):
                                        safe[k] = float(v)
                                    elif isinstance(v, (int, float, str)):
                                        safe[k] = v
                                cache[key] = safe
                                _save_priors_cache(cache)
                            return result
                        else:
                            logger.warning(f"❌ Source {source_name}: valeurs invalides (a/R*={a_rs_val}, i={incl_val}°) - a/R* doit être > 0, i doit être entre 0 et 90°")
                    else:
                        logger.warning(f"❌ Source {source_name}: valeurs None (a/R*={a_rs_val}, i={incl_val})")
                else:
                    logger.warning(f"❌ Source {source_name}: clés manquantes dans le résultat. Clés présentes: {list(result.keys())}")
            else:
                logger.info(f"  → Source {source_name} n'a pas trouvé de données pour {target_name}")
        except Exception as e:
            logger.warning(f"❌ Erreur avec source {source_name} pour {target_name}: {e}", exc_info=True)
            continue
    
    logger.warning(f"⚠️ Aucune source n'a fourni de priors valides pour {target_name}")
    logger.info(f"Sources essayées: {', '.join(sources_order)}")
    return None


def format_priors_for_display(priors: Dict[str, float]) -> str:
    """
    Formate les priors pour l'affichage dans l'interface utilisateur.
    
    Parameters
    ----------
    priors : dict
        Dictionnaire contenant les priors
    
    Returns
    -------
    str
        Chaîne formatée pour l'affichage
    """
    if priors is None:
        return "Aucun prior disponible"
    
    a_rs = priors.get('a_rs', None)
    i = priors.get('inclination', None)
    source = priors.get('source', 'Inconnue')
    
    a_rs_err = priors.get('a_rs_err', None)
    i_err = priors.get('inclination_err', None)
    
    lines = [f"Source: {source}"]
    
    if a_rs is not None:
        if a_rs_err is not None:
            lines.append(f"a/R* = {a_rs:.4f} ± {a_rs_err:.4f}")
        else:
            lines.append(f"a/R* = {a_rs:.4f}")
    
    if i is not None:
        if i_err is not None:
            lines.append(f"i = {i:.4f}° ± {i_err:.4f}°")
        else:
            lines.append(f"i = {i:.4f}°")
    
    return "\n".join(lines)
