# core/gaia_catalog_overlay.py
"""
Module pour charger et afficher les catalogues Gaia DR3 (étoiles variables et galaxies)
sur les images FITS dans les fenêtres de visualisation.
"""

import logging
import gzip
from pathlib import Path
from typing import Optional, List, Tuple
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

logger = logging.getLogger(__name__)


class GaiaCatalogOverlay:
    """
    Classe pour charger et gérer les catalogues Gaia DR3 pour l'affichage sur les images.
    """
    
    def __init__(self, catalog_dir: Optional[Path] = None):
        """
        Initialise le gestionnaire de catalogues.
        
        Parameters
        ----------
        catalog_dir : Path, optional
            Répertoire contenant les catalogues Gaia DR3 (CSV.GZIP)
        """
        self.catalog_dir = Path(catalog_dir) if catalog_dir else None
        self.variables_cache = {}  # Cache pour les étoiles variables par fichier
        self.galaxies_cache = {}   # Cache pour les galaxies par fichier
    
    def find_catalog_files(
        self, 
        center_coord: SkyCoord, 
        search_radius_deg: float = 1.0,
        catalog_type: str = "both"  # "variables", "galaxies", "both"
    ) -> Tuple[List[Path], List[Path]]:
        """
        Trouve les fichiers de catalogues Gaia DR3 qui couvrent une région du ciel.
        
        Parameters
        ----------
        center_coord : SkyCoord
            Coordonnées du centre de la région
        search_radius_deg : float
            Rayon de recherche en degrés
        catalog_type : str
            Type de catalogue : "variables", "galaxies", ou "both"
        
        Returns
        -------
        tuple
            (liste fichiers variables, liste fichiers galaxies)
        """
        if self.catalog_dir is None or not self.catalog_dir.exists():
            return [], []
        
        variables_files = []
        galaxies_files = []
        
        # Calculer les limites RA/DEC de la région
        ra_center = center_coord.ra.deg
        dec_center = center_coord.dec.deg
        
        ra_min = (ra_center - search_radius_deg) % 360
        ra_max = (ra_center + search_radius_deg) % 360
        dec_min = dec_center - search_radius_deg
        dec_max = dec_center + search_radius_deg
        
        # Chercher les fichiers de catalogues
        if catalog_type in ["variables", "both"]:
            var_patterns = [
                "gaia_dr3_*_variables*.csv.gz",
                "gaia_dr3_*_variables*.csv"
            ]
            for pattern in var_patterns:
                for f in self.catalog_dir.glob(pattern):
                    if self._file_covers_region(f, ra_min, ra_max, dec_min, dec_max):
                        variables_files.append(f)
        
        if catalog_type in ["galaxies", "both"]:
            gal_patterns = [
                "gaia_dr3_*_galaxies*.csv.gz",
                "gaia_dr3_*_galaxies*.csv"
            ]
            for pattern in gal_patterns:
                for f in self.catalog_dir.glob(pattern):
                    if self._file_covers_region(f, ra_min, ra_max, dec_min, dec_max):
                        galaxies_files.append(f)
        
        return variables_files, galaxies_files
    
    def _file_covers_region(
        self, 
        file_path: Path, 
        ra_min: float, 
        ra_max: float, 
        dec_min: float, 
        dec_max: float
    ) -> bool:
        """
        Vérifie si un fichier de catalogue couvre la région spécifiée.
        Parse le nom de fichier pour extraire les plages RA/DEC.
        """
        name = file_path.stem.replace('.csv', '')
        
        # Parser le nom : gaia_dr3_{hem}_ra{min}h-{max}h_mag{mag}[_suffix]
        # ou gaia_dr3_{hem}_ra{min}h{min_m}m-{max}h{max_m}m_mag{mag}[_suffix]
        try:
            if '_ra' in name:
                ra_part = name.split('_ra')[1].split('_')[0]
                
                # Format simple : ra00h-02h
                if 'h-' in ra_part and 'm' not in ra_part:
                    ra_min_h_str, ra_max_h_str = ra_part.split('h-')
                    ra_min_h = int(ra_min_h_str)
                    ra_max_h = int(ra_max_h_str.split('h')[0])
                    file_ra_min = ra_min_h * 15.0
                    file_ra_max = ra_max_h * 15.0
                # Format avec minutes : ra00h00m-00h20m
                elif 'h' in ra_part and 'm-' in ra_part:
                    parts = ra_part.split('-')
                    if len(parts) == 2:
                        min_part = parts[0]
                        max_part = parts[1]
                        # Parser 00h00m
                        if 'h' in min_part and 'm' in min_part:
                            min_h = int(min_part.split('h')[0])
                            min_m = int(min_part.split('h')[1].split('m')[0])
                            max_h = int(max_part.split('h')[0])
                            max_m = int(max_part.split('h')[1].split('m')[0])
                            file_ra_min = (min_h + min_m / 60.0) * 15.0
                            file_ra_max = (max_h + max_m / 60.0) * 15.0
                        else:
                            return True  # Si parsing échoue, inclure le fichier
                    else:
                        return True
                else:
                    return True  # Format non reconnu, inclure par sécurité
                
                # Vérifier le chevauchement RA (gérer le cas où ra_max > 360)
                if file_ra_max < file_ra_min:  # Cas où le fichier couvre 360°
                    ra_overlap = (ra_min < file_ra_max) or (ra_max > file_ra_min)
                else:
                    ra_overlap = (ra_min < file_ra_max) and (ra_max > file_ra_min)
                
                # Pour DEC, vérifier l'hémisphère dans le nom
                if 'nord' in name or 'north' in name:
                    file_dec_min = -5.0
                    file_dec_max = 90.0
                elif 'sud' in name or 'south' in name:
                    file_dec_min = -90.0
                    file_dec_max = 5.0
                else:
                    file_dec_min = -90.0
                    file_dec_max = 90.0
                
                dec_overlap = (dec_min < file_dec_max) and (dec_max > file_dec_min)
                
                return ra_overlap and dec_overlap
        except Exception as e:
            logger.debug(f"Erreur parsing nom fichier {file_path.name}: {e}")
            return True  # En cas d'erreur, inclure le fichier par sécurité
        
        return True
    
    def load_catalog_file(self, file_path: Path) -> pd.DataFrame:
        """
        Charge un fichier de catalogue Gaia DR3 (CSV ou CSV.GZIP).
        
        Parameters
        ----------
        file_path : Path
            Chemin vers le fichier de catalogue
        
        Returns
        -------
        pd.DataFrame
            DataFrame avec les colonnes du catalogue
        """
        try:
            if file_path.suffix == '.gz' or file_path.name.endswith('.csv.gz'):
                # Fichier compressé
                with gzip.open(file_path, 'rt', encoding='utf-8') as f:
                    df = pd.read_csv(f, comment='#', low_memory=False)
            else:
                # Fichier CSV normal
                df = pd.read_csv(file_path, comment='#', low_memory=False)
            
            logger.debug(f"Fichier chargé : {file_path.name} ({len(df)} objets)")
            return df
        except Exception as e:
            logger.error(f"Erreur chargement {file_path.name}: {e}")
            return pd.DataFrame()
    
    def get_objects_in_region(
        self,
        center_coord: SkyCoord,
        search_radius_deg: float,
        catalog_type: str = "both",  # "variables", "galaxies", "both"
        wcs: Optional[WCS] = None,
        image_shape: Optional[Tuple[int, int]] = None  # (height, width) pour filtrage pixel
    ) -> Tuple[List[SkyCoord], List[float], List[str]]:
        """
        Récupère les objets (étoiles variables ou galaxies) dans une région du ciel.
        
        Parameters
        ----------
        center_coord : SkyCoord
            Coordonnées du centre
        search_radius_deg : float
            Rayon de recherche en degrés
        catalog_type : str
            Type de catalogue : "variables", "galaxies", ou "both"
        wcs : WCS, optional
            WCS de l'image pour filtrer les objets hors champ
        
        Returns
        -------
        tuple
            (liste SkyCoord, liste magnitudes G, liste types)
        """
        coords = []
        mags = []
        types = []
        
        # Trouver les fichiers de catalogues
        var_files, gal_files = self.find_catalog_files(
            center_coord, 
            search_radius_deg, 
            catalog_type
        )
        
        # Charger les étoiles variables
        if catalog_type in ["variables", "both"]:
            for var_file in var_files:
                if var_file not in self.variables_cache:
                    self.variables_cache[var_file] = self.load_catalog_file(var_file)
                
                df = self.variables_cache[var_file]
                if len(df) == 0:
                    continue
                
                # Extraire les colonnes RA/DEC et magnitude
                ra_col = None
                dec_col = None
                mag_col = None
                
                for col in df.columns:
                    col_lower = col.lower()
                    if ra_col is None and ('ra' == col_lower or col_lower.startswith('ra_') or 'ra' in col_lower):
                        ra_col = col
                    if dec_col is None and ('dec' == col_lower or col_lower.startswith('dec_') or 'dec' in col_lower):
                        dec_col = col
                    if mag_col is None and ('phot_g_mean_mag' in col_lower or 'g_mean_mag' in col_lower or 'gmag' in col_lower):
                        mag_col = col
                
                if ra_col and dec_col:
                    try:
                        # Filtrer par région
                        df_region = df.copy()
                        df_region['coord'] = SkyCoord(
                            ra=df_region[ra_col] * u.deg,
                            dec=df_region[dec_col] * u.deg
                        )
                        df_region['sep'] = center_coord.separation(df_region['coord']).deg
                        df_region = df_region[df_region['sep'] <= search_radius_deg]
                        
                        # Filtrer par champ de l'image si WCS et image_shape disponibles
                        if wcs is not None and image_shape is not None:
                            try:
                                h, w = image_shape
                                valid_rows = []
                                for _, row in df_region.iterrows():
                                    try:
                                        px, py = wcs.world_to_pixel(row['coord'])
                                        if 0 <= px < w and 0 <= py < h:
                                            valid_rows.append(row)
                                    except:
                                        continue
                                df_region = pd.DataFrame(valid_rows) if valid_rows else pd.DataFrame()
                            except Exception as e:
                                logger.debug(f"Erreur filtrage WCS: {e}")
                        
                        # Ajouter aux listes
                        for _, row in df_region.iterrows():
                            coords.append(row['coord'])
                            mag = row[mag_col] if mag_col and mag_col in row else None
                            mags.append(mag if mag is not None and not pd.isna(mag) else 99.0)
                            types.append("variable")
                    except Exception as e:
                        logger.error(f"Erreur traitement variables {var_file.name}: {e}")
        
        # Charger les galaxies
        if catalog_type in ["galaxies", "both"]:
            for gal_file in gal_files:
                if gal_file not in self.galaxies_cache:
                    self.galaxies_cache[gal_file] = self.load_catalog_file(gal_file)
                
                df = self.galaxies_cache[gal_file]
                if len(df) == 0:
                    continue
                
                # Extraire les colonnes RA/DEC et magnitude
                ra_col = None
                dec_col = None
                mag_col = None
                
                for col in df.columns:
                    col_lower = col.lower()
                    if ra_col is None and ('ra' == col_lower or col_lower.startswith('ra_') or 'ra' in col_lower):
                        ra_col = col
                    if dec_col is None and ('dec' == col_lower or col_lower.startswith('dec_') or 'dec' in col_lower):
                        dec_col = col
                    if mag_col is None and ('phot_g_mean_mag' in col_lower or 'g_mean_mag' in col_lower or 'gmag' in col_lower):
                        mag_col = col
                
                if ra_col and dec_col:
                    try:
                        # Filtrer par région
                        df_region = df.copy()
                        df_region['coord'] = SkyCoord(
                            ra=df_region[ra_col] * u.deg,
                            dec=df_region[dec_col] * u.deg
                        )
                        df_region['sep'] = center_coord.separation(df_region['coord']).deg
                        df_region = df_region[df_region['sep'] <= search_radius_deg]
                        
                        # Filtrer par champ de l'image si WCS et image_shape disponibles
                        if wcs is not None and image_shape is not None:
                            try:
                                h, w = image_shape
                                valid_rows = []
                                for _, row in df_region.iterrows():
                                    try:
                                        px, py = wcs.world_to_pixel(row['coord'])
                                        if 0 <= px < w and 0 <= py < h:
                                            valid_rows.append(row)
                                    except:
                                        continue
                                df_region = pd.DataFrame(valid_rows) if valid_rows else pd.DataFrame()
                            except Exception as e:
                                logger.debug(f"Erreur filtrage WCS: {e}")
                        
                        # Ajouter aux listes
                        for _, row in df_region.iterrows():
                            coords.append(row['coord'])
                            mag = row[mag_col] if mag_col and mag_col in row else None
                            mags.append(mag if mag is not None and not pd.isna(mag) else 99.0)
                            types.append("galaxy")
                    except Exception as e:
                        logger.error(f"Erreur traitement galaxies {gal_file.name}: {e}")
        
        logger.info(f"Objets trouvés : {len([t for t in types if t == 'variable'])} variables, {len([t for t in types if t == 'galaxy'])} galaxies")
        return coords, mags, types
