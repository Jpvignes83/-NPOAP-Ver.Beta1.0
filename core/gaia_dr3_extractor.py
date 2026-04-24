# core/gaia_dr3_extractor.py
"""
Module pour extraire des données du catalogue Gaia DR3 via l'API TAP de l'archive ESA.

Ce module permet de créer des catalogues filtrés par hémisphère (nord/sud) et 
répartis par bandes de longitude horaire (RA) pour optimiser les téléchargements.
"""

import logging
import os
import gzip
from pathlib import Path
from typing import Optional, Tuple
import numpy as np
from astropy.table import Table
from astropy import units as u
import time

logger = logging.getLogger(__name__)

try:
    from astroquery.gaia import Gaia
    GAIA_TAP_AVAILABLE = True
except ImportError:
    GAIA_TAP_AVAILABLE = False
    logger.warning("astroquery.gaia non disponible. Les requêtes TAP ne fonctionneront pas.")


class GaiaDR3Extractor:
    """
    Extracteur de données Gaia DR3 via l'API TAP de l'archive ESA.
    
    Permet d'extraire des catalogues filtrés par :
    - Hémisphère (nord/sud)
    - Magnitude limite (G < mag_limit)
    - Répartition par bandes de RA (longitude horaire)
    """
    
    def __init__(self, output_dir: Optional[str | Path] = None):
        """
        Initialise l'extracteur Gaia DR3.
        
        Parameters
        ----------
        output_dir : str | Path, optional
            Répertoire de sortie pour les catalogues. Si None, utilise 'gaia_catalogues/'
        """
        if not GAIA_TAP_AVAILABLE:
            raise RuntimeError(
                "astroquery.gaia n'est pas disponible.\n"
                "Installez-le avec : pip install astroquery"
            )
        
        self.output_dir = Path(output_dir) if output_dir else Path("gaia_catalogues")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Configuration de l'API Gaia
        # Utiliser l'archive ESA principale
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = -1  # Pas de limite par défaut
        
        logger.info(f"Extracteur Gaia DR3 initialisé. Répertoire de sortie : {self.output_dir}")
    
    def extract_hemisphere_catalog(
        self,
        hemisphere: str,
        mag_limit: float = 15.0,
        ra_step_deg: float = 30.0,
        output_format: str = "csv",
        combine_files: bool = False,
        filter_variables: bool = False,
        filter_galaxies: bool = False,
        skip_existing: bool = True,
        ra_min_deg: Optional[float] = None,
        ra_max_deg: Optional[float] = None,
        dec_min: Optional[float] = None,
        dec_max: Optional[float] = None
    ) -> Tuple[list[Path], int]:
        """
        Extrait un catalogue pour un hémisphère avec répartition par RA.
        
        Parameters
        ----------
        hemisphere : str
            'north' ou 'south' (hémisphère nord ou sud)
        mag_limit : float
            Magnitude limite G (par défaut 15.0)
        ra_step_deg : float
            Pas de RA en degrés pour diviser le ciel (par défaut 30° = 2h)
        output_format : str
            Format de sortie : 'csv', 'csv.gz', 'fits', ou 'votable' (par défaut 'csv')
        combine_files : bool
            Si True, combine tous les fichiers en un seul catalogue final
        
        Returns
        -------
        tuple
            (liste des fichiers créés, nombre total d'étoiles)
        """
        hemisphere = hemisphere.lower()
        if hemisphere not in ['north', 'south']:
            raise ValueError(f"hemisphere doit être 'north' ou 'south', reçu : {hemisphere}")
        
        # Conditions de déclinaison
        # Si des limites DEC sont spécifiées, les utiliser
        # Sinon, utiliser les limites par défaut pour l'hémisphère
        if dec_min is not None or dec_max is not None:
            # Vérifier et corriger si les limites sont inversées
            if dec_min is not None and dec_max is not None:
                if dec_min > dec_max:
                    # Les limites sont inversées, les échanger
                    logger.warning(f"⚠️ Limites DEC inversées détectées ({dec_min}° > {dec_max}°), inversion automatique")
                    dec_min, dec_max = dec_max, dec_min
                    logger.info(f"   Limites DEC corrigées: {dec_min}° - {dec_max}°")
            
            dec_conditions = []
            if dec_min is not None:
                dec_conditions.append(f"dec >= {dec_min}")
            if dec_max is not None:
                dec_conditions.append(f"dec < {dec_max}")
            dec_condition = " AND ".join(dec_conditions)
            # Déterminer le nom de l'hémisphère selon la plage DEC
            if dec_min is not None and dec_min >= -5:
                hem_name = "nord"
            elif dec_max is not None and dec_max < 5:
                hem_name = "sud"
            else:
                hem_name = "nord" if hemisphere == 'north' else "sud"
            logger.info(f"Plage DEC filtrée : {dec_min or -90}° - {dec_max or 90}°")
        else:
            # Utiliser les limites par défaut pour l'hémisphère
            if hemisphere == 'north':
                dec_condition = "dec >= -5"
                hem_name = "nord"
            else:
                dec_condition = "dec < 5"
                hem_name = "sud"
        
        logger.info(f"Extraction du catalogue hémisphère {hem_name} (G < {mag_limit})")
        logger.info(f"Répartition par RA : pas de {ra_step_deg}° ({ra_step_deg/15:.1f}h)")
        
        # Diviser le ciel en bandes de RA
        # Si des limites RA sont spécifiées, les utiliser
        ra_start = ra_min_deg if ra_min_deg is not None else 0.0
        ra_end = ra_max_deg if ra_max_deg is not None else 360.0
        
        # Log des limites globales appliquées
        if ra_min_deg is not None or ra_max_deg is not None:
            logger.info(f"🎯 LIMITES RA APPLIQUÉES: {ra_start:.1f}° - {ra_end:.1f}° ({ra_start/15:.1f}h - {ra_end/15:.1f}h)")
        else:
            logger.info(f"🌐 Aucune limite RA: téléchargement complet (0° - 360°)")
        
        if dec_min is not None or dec_max is not None:
            logger.info(f"🎯 LIMITES DEC APPLIQUÉES: {dec_min or -90:.1f}° - {dec_max or 90:.1f}°")
        else:
            logger.info(f"🌐 Limites DEC par défaut pour hémisphère {hem_name}: {dec_condition}")
        
        n_bands = int(np.ceil((ra_end - ra_start) / ra_step_deg))
        ra_ranges = []
        for i in range(n_bands):
            ra_min = ra_start + i * ra_step_deg
            ra_max = min(ra_start + (i + 1) * ra_step_deg, ra_end)
            ra_ranges.append((ra_min, ra_max))
        
        logger.info(f"Nombre de bandes de RA : {len(ra_ranges)}")
        
        output_files = []
        total_stars = 0
        
        for band_idx, (ra_min, ra_max) in enumerate(ra_ranges):
            # Construire le nom de fichier pour vérifier s'il existe déjà
            ra_min_h = int(ra_min / 15)
            ra_max_h = int(ra_max / 15) if ra_max < 360 else 24
            filter_suffix = ""
            if filter_variables and filter_galaxies:
                filter_suffix = "_variables_galaxies"
            elif filter_variables:
                filter_suffix = "_variables"
            elif filter_galaxies:
                filter_suffix = "_galaxies"
            filename = f"gaia_dr3_{hem_name}_ra{ra_min_h:02d}h-{ra_max_h:02d}h_mag{mag_limit:.1f}{filter_suffix}.{output_format}"
            output_path = self.output_dir / filename
            
            # Vérifier si le fichier existe déjà (AVANT la requête pour éviter de télécharger inutilement)
            if skip_existing and output_path.exists():
                logger.info(f"Traitement bande {band_idx + 1}/{len(ra_ranges)} : RA {ra_min:.1f}° - {ra_max:.1f}°")
                logger.info(f"  ⏭ Fichier existe déjà, ignoré : {filename}")
                output_files.append(output_path)
                # Compter les lignes du fichier existant (soustraire l'en-tête)
                try:
                    if str(output_path).endswith('.csv.gz'):
                        with gzip.open(output_path, 'rt', encoding='utf-8') as f:
                            n_stars = sum(1 for line in f) - 1
                    else:
                        with open(output_path, 'r', encoding='utf-8') as f:
                            n_stars = sum(1 for line in f) - 1
                    total_stars += n_stars
                except:
                    pass
                continue
            
            logger.info(f"Traitement bande {band_idx + 1}/{len(ra_ranges)} : RA {ra_min:.1f}° - {ra_max:.1f}°")
            
            # Construire la requête ADQL
            # Pour gérer les cas où ra_max > 360, on utilise modulo
            if ra_max >= 360:
                # Cas spécial : dernière bande qui peut dépasser 360°
                ra_condition = f"(ra >= {ra_min} OR ra < {ra_max % 360})"
            else:
                ra_condition = f"ra >= {ra_min} AND ra < {ra_max}"
            
            # Log des conditions appliquées pour cette bande
            logger.debug(f"  Conditions de requête pour bande {band_idx + 1}:")
            logger.debug(f"    RA: {ra_condition}")
            logger.debug(f"    DEC: {dec_condition}")
            if dec_min is not None or dec_max is not None:
                logger.info(f"  📍 Zone limitée: RA [{ra_min:.1f}°-{ra_max:.1f}°], DEC [{dec_min or -90:.1f}°-{dec_max or 90:.1f}°]")
            else:
                logger.info(f"  📍 Zone: RA [{ra_min:.1f}°-{ra_max:.1f}°], DEC [{dec_condition}]")
            
            # Construire la requête ADQL avec filtres optionnels
            # Pour les étoiles variables : JOIN avec vari_summary
            # Pour les galaxies : JOIN avec galaxy_candidates
            # Selon la documentation Gaia DR3 : https://gea.esac.esa.int/archive/documentation/GDR3/
            where_conditions = [
                ra_condition,
                dec_condition,
                f"phot_g_mean_mag < {mag_limit}",
                "phot_g_mean_mag IS NOT NULL"
            ]
            
            if filter_variables and filter_galaxies:
                # Les deux filtres : sources qui sont à la fois variables ET galaxies (rare)
                # Utiliser phot_variable_flag au lieu de JOIN avec vari_summary
                query = f"""
                SELECT gs.*
                FROM gaiadr3.gaia_source AS gs
                INNER JOIN gaiadr3.galaxy_candidates AS gc ON gc.source_id = gs.source_id
                WHERE 
                    {ra_condition}
                    AND {dec_condition}
                    AND gs.phot_g_mean_mag < {mag_limit}
                    AND gs.phot_g_mean_mag IS NOT NULL
                    AND gs.phot_variable_flag = 'VARIABLE'
                    AND gc.classlabel_dsc_joint = 'galaxy'
                """
            elif filter_variables:
                # Filtre étoiles variables uniquement
                # Utiliser phot_variable_flag au lieu de JOIN avec vari_summary pour éviter les erreurs 500
                # phot_variable_flag = 'VARIABLE' indique une étoile variable confirmée
                query = f"""
                SELECT *
                FROM gaiadr3.gaia_source
                WHERE 
                    {ra_condition}
                    AND {dec_condition}
                    AND phot_g_mean_mag < {mag_limit}
                    AND phot_g_mean_mag IS NOT NULL
                    AND phot_variable_flag = 'VARIABLE'
                """
            elif filter_galaxies:
                # Filtre galaxies uniquement
                query = f"""
                SELECT gs.*
                FROM gaiadr3.gaia_source AS gs
                INNER JOIN gaiadr3.galaxy_candidates AS gc ON gc.source_id = gs.source_id
                WHERE 
                    {ra_condition}
                    AND {dec_condition}
                    AND gs.phot_g_mean_mag < {mag_limit}
                    AND gs.phot_g_mean_mag IS NOT NULL
                    AND gc.classlabel_dsc_joint = 'galaxy'
                """
            else:
                # Pas de filtre : requête standard
                query = f"""
                SELECT *
                FROM gaiadr3.gaia_source
                WHERE 
                    {ra_condition}
                    AND {dec_condition}
                    AND phot_g_mean_mag < {mag_limit}
                    AND phot_g_mean_mag IS NOT NULL
                """
            
            # Tentatives avec retry pour gérer les erreurs 500
            max_retries = 3
            retry_delay = 5  # secondes
            success = False
            results = None
            
            for retry in range(max_retries):
                try:
                    # Exécuter la requête TAP
                    if retry > 0:
                        logger.debug(f"Tentative {retry + 1}/{max_retries} pour RA {ra_min:.1f}°-{ra_max:.1f}°")
                        time.sleep(retry_delay * retry)  # Délai progressif
                    
                    logger.debug(f"Requête ADQL : {query[:200]}...")
                    job = Gaia.launch_job_async(query, dump_to_file=False)
                    results = job.get_results()
                    
                    n_stars = len(results)
                    total_stars += n_stars
                    logger.info(f"  → {n_stars:,} étoiles récupérées")
                    success = True
                    break  # Succès, sortir de la boucle de retry
                    
                except Exception as e:
                    error_msg = str(e)
                    # Vérifier si c'est une erreur 500 (serveur)
                    if "500" in error_msg or "Error 500" in error_msg:
                        if retry < max_retries - 1:
                            logger.warning(f"  ⚠️ Erreur serveur 500 (tentative {retry + 1}/{max_retries}), nouvelle tentative dans {retry_delay * (retry + 1)}s...")
                            continue
                        else:
                            logger.error(f"  ❌ Erreur serveur 500 persistante après {max_retries} tentatives")
                    else:
                        # Autre type d'erreur, ne pas retry
                        logger.error(f"  ❌ Erreur lors de la requête RA {ra_min:.1f}°-{ra_max:.1f}° : {error_msg}")
                        break
            
            if success and results is not None and len(results) > 0:
                try:
                    # Le nom de fichier et output_path ont déjà été définis avant la requête
                    # S'assurer que le répertoire de sortie existe
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    
                    # Sauvegarder selon le format
                    if output_format.lower() == "csv":
                        results.write(str(output_path), format='csv', overwrite=True)
                    elif output_format.lower() == "csv.gz":
                        # Sauvegarder en CSV compressé
                        # Écrire d'abord en CSV temporaire, puis compresser
                        temp_csv = str(output_path).replace('.csv.gz', '.csv.tmp')
                        results.write(temp_csv, format='csv', overwrite=True)
                        # Compresser le fichier
                        with open(temp_csv, 'rb') as f_in:
                            with gzip.open(str(output_path), 'wb') as f_out:
                                f_out.writelines(f_in)
                        # Supprimer le fichier temporaire
                        os.remove(temp_csv)
                    elif output_format.lower() == "fits":
                        results.write(str(output_path), format='fits', overwrite=True)
                    elif output_format.lower() == "votable":
                        results.write(str(output_path), format='votable', overwrite=True)
                    else:
                        raise ValueError(f"Format non supporté : {output_format}")
                    
                    output_files.append(output_path)
                    logger.info(f"  → Fichier sauvegardé : {filename}")
                except Exception as e:
                    logger.error(f"  ❌ Erreur lors de la sauvegarde : {e}")
            
            # Pause entre les requêtes pour éviter de surcharger l'API
            time.sleep(2)  # Augmenté à 2 secondes pour réduire la charge
        
        logger.info(f"Catalogue {hem_name} : {total_stars:,} étoiles au total dans {len(output_files)} fichiers")
        
        # Combiner les fichiers si demandé
        if combine_files and len(output_files) > 1:
            combined_file = self._combine_tables(
                output_files, 
                hemisphere=hem_name,
                mag_limit=mag_limit,
                output_format=output_format
            )
            if combined_file:
                output_files.append(combined_file)
                logger.info(f"Fichier combiné créé : {combined_file.name}")
        
        return output_files, total_stars
    
    def _combine_tables(
        self, 
        table_files: list[Path],
        hemisphere: str,
        mag_limit: float,
        output_format: str
    ) -> Optional[Path]:
        """
        Combine plusieurs tables en une seule.
        
        Parameters
        ----------
        table_files : list[Path]
            Liste des fichiers de tables à combiner
        hemisphere : str
            Nom de l'hémisphère ('nord' ou 'sud')
        mag_limit : float
            Magnitude limite
        output_format : str
            Format de sortie
        
        Returns
        -------
        Path ou None
            Chemin du fichier combiné, ou None en cas d'erreur
        """
        try:
            logger.info(f"Combinaison de {len(table_files)} fichiers...")
            
            combined_table = None
            for file_path in table_files:
                # Déterminer le format de lecture selon l'extension
                read_format = output_format
                if str(file_path).endswith('.csv.gz'):
                    read_format = 'csv'
                elif str(file_path).endswith('.csv'):
                    read_format = 'csv'
                elif str(file_path).endswith('.fits'):
                    read_format = 'fits'
                elif str(file_path).endswith('.votable') or str(file_path).endswith('.xml'):
                    read_format = 'votable'
                
                # Lire le fichier (gérer la compression gzip)
                if str(file_path).endswith('.csv.gz'):
                    with gzip.open(str(file_path), 'rt') as f:
                        table = Table.read(f, format='csv')
                else:
                    table = Table.read(str(file_path), format=read_format)
                
                if combined_table is None:
                    combined_table = table
                else:
                    combined_table = np.concatenate([combined_table, table])
            
            if combined_table is None:
                return None
            
            # Tri par RA pour faciliter l'utilisation
            combined_table.sort('ra')
            
            # Nom du fichier combiné
            filename = f"gaia_dr3_{hemisphere}_complete_mag{mag_limit:.1f}.{output_format}"
            output_path = self.output_dir / filename
            
            # Sauvegarder
            if output_format.lower() == "csv":
                combined_table.write(str(output_path), format='csv', overwrite=True)
            elif output_format.lower() == "csv.gz":
                # Sauvegarder en CSV compressé
                temp_csv = str(output_path).replace('.csv.gz', '.csv.tmp')
                combined_table.write(temp_csv, format='csv', overwrite=True)
                # Compresser le fichier
                with open(temp_csv, 'rb') as f_in:
                    with gzip.open(str(output_path), 'wb') as f_out:
                        f_out.writelines(f_in)
                # Supprimer le fichier temporaire
                os.remove(temp_csv)
            elif output_format.lower() == "fits":
                combined_table.write(str(output_path), format='fits', overwrite=True)
            elif output_format.lower() == "votable":
                combined_table.write(str(output_path), format='votable', overwrite=True)
            
            logger.info(f"Fichier combiné : {len(combined_table):,} étoiles")
            return output_path
            
        except Exception as e:
            logger.error(f"Erreur lors de la combinaison des fichiers : {e}")
            return None
    
    def extract_both_hemispheres(
        self,
        mag_limit: float = 15.0,
        ra_step_deg: float = 30.0,
        output_format: str = "csv",
        combine_files: bool = False,
        filter_variables: bool = False,
        filter_galaxies: bool = False,
        skip_existing: bool = True,
        ra_min_deg: Optional[float] = None,
        ra_max_deg: Optional[float] = None,
        dec_min: Optional[float] = None,
        dec_max: Optional[float] = None
    ) -> Tuple[dict, int, int]:
        """
        Extrait les catalogues pour les deux hémisphères.
        
        Parameters
        ----------
        mag_limit : float
            Magnitude limite G
        ra_step_deg : float
            Pas de RA en degrés
        output_format : str
            Format de sortie : 'csv', 'csv.gz', 'fits', ou 'votable'
        combine_files : bool
            Si True, combine les fichiers par hémisphère
        
        Returns
        -------
        tuple
            (dict avec clés 'north' et 'south' contenant les listes de fichiers,
             nombre d'étoiles nord, nombre d'étoiles sud)
        """
        logger.info("=" * 80)
        logger.info("EXTRACTION COMPLÈTE DES CATALOGUES GAIA DR3")
        logger.info("=" * 80)
        
        # Extraire hémisphère nord
        logger.info("\n[1/2] Extraction hémisphère NORD")
        north_files, north_count = self.extract_hemisphere_catalog(
            hemisphere='north',
            mag_limit=mag_limit,
            ra_step_deg=ra_step_deg,
            output_format=output_format,
            combine_files=combine_files,
            filter_variables=filter_variables,
            filter_galaxies=filter_galaxies,
            skip_existing=skip_existing,
            ra_min_deg=ra_min_deg,
            ra_max_deg=ra_max_deg,
            dec_min=dec_min,
            dec_max=dec_max
        )
        
        # Extraire hémisphère sud
        logger.info("\n[2/2] Extraction hémisphère SUD")
        south_files, south_count = self.extract_hemisphere_catalog(
            hemisphere='south',
            mag_limit=mag_limit,
            ra_step_deg=ra_step_deg,
            output_format=output_format,
            combine_files=combine_files,
            filter_variables=filter_variables,
            filter_galaxies=filter_galaxies,
            skip_existing=skip_existing,
            ra_min_deg=ra_min_deg,
            ra_max_deg=ra_max_deg,
            dec_min=dec_min,
            dec_max=dec_max
        )
        
        results = {
            'north': north_files,
            'south': south_files
        }
        
        logger.info("\n" + "=" * 80)
        logger.info("RÉSUMÉ")
        logger.info("=" * 80)
        logger.info(f"Hémisphère NORD : {north_count:,} étoiles dans {len(north_files)} fichiers")
        logger.info(f"Hémisphère SUD  : {south_count:,} étoiles dans {len(south_files)} fichiers")
        logger.info(f"TOTAL           : {north_count + south_count:,} étoiles")
        logger.info(f"Répertoire de sortie : {self.output_dir}")
        logger.info("=" * 80)
        
        return results, north_count, south_count


# =============================================================================
# FONCTION D'UTILISATION DIRECTE
# =============================================================================

def extract_gaia_catalogs(
    output_dir: str | Path = "gaia_catalogues",
    mag_limit: float = 15.0,
    ra_step_deg: float = 30.0,
    hemisphere: Optional[str] = None,
    output_format: str = "csv",
    combine_files: bool = False
) -> dict:
    """
    Fonction utilitaire pour extraire les catalogues Gaia DR3.
    
    Parameters
    ----------
    output_dir : str | Path
        Répertoire de sortie
    mag_limit : float
        Magnitude limite G (par défaut 15.0)
    ra_step_deg : float
        Pas de RA en degrés pour diviser le ciel (par défaut 30° = 2h)
    hemisphere : str, optional
        'north', 'south', ou None pour les deux (par défaut None)
    output_format : str
        Format de sortie : 'csv', 'csv.gz', 'fits', ou 'votable'
    combine_files : bool
        Si True, combine les fichiers par hémisphère
        
    Returns
    -------
    dict
        Dictionnaire avec les résultats de l'extraction
        
    Example
    -------
    >>> # Extraire les deux hémisphères avec G < 15, pas de 30°
    >>> results = extract_gaia_catalogs(
    ...     output_dir="catalogues_gaia",
    ...     mag_limit=15.0,
    ...     ra_step_deg=30.0
    ... )
    """
    extractor = GaiaDR3Extractor(output_dir=output_dir)
    
    if hemisphere is None:
        results, north_count, south_count = extractor.extract_both_hemispheres(
            mag_limit=mag_limit,
            ra_step_deg=ra_step_deg,
            output_format=output_format,
            combine_files=combine_files
        )
        return {
            'files': results,
            'north_count': north_count,
            'south_count': south_count,
            'total_count': north_count + south_count
        }
    else:
        files, count = extractor.extract_hemisphere_catalog(
            hemisphere=hemisphere,
            mag_limit=mag_limit,
            ra_step_deg=ra_step_deg,
            output_format=output_format,
            combine_files=combine_files
        )
        return {
            'files': files,
            'count': count,
            'hemisphere': hemisphere
        }
