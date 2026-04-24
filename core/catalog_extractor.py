# core/catalog_extractor.py
"""
Module pour extraire des données depuis divers catalogues astronomiques via astroquery/Vizier.
"""

import gzip
import logging
import os
from pathlib import Path
from typing import Optional, List, Dict, Any
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import numpy as np

logger = logging.getLogger(__name__)

# Vérifier la disponibilité d'astroquery
try:
    from astroquery.vizier import Vizier
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False
    logger.warning("astroquery non disponible. L'extraction de catalogues ne fonctionnera pas.")


# Dictionnaire des catalogues populaires disponibles
POPULAR_CATALOGS = {
    "Gaia DR3": {
        "vizier_id": "I/355/gaiadr3",
        "description": "Catalogue Gaia Data Release 3 - Astrométrie et photométrie précises",
        "object_types": ["Étoiles"],
        "mag_column": "Gmag"
    },
    "TIC-8.2": {
        "vizier_id": "IV/39/tic82",
        "description": "TESS Input Catalog version 8.2",
        "object_types": ["Étoiles"],
        "mag_column": "Tmag"
    },
    "TESS EBS": {
        "vizier_id": None,  # TESS EBS n'est pas dans Vizier, doit être téléchargé depuis le site web
        "url": "https://tessebs.villanova.edu/",
        "description": "TESS Eclipsing Binary Stars - Catalogue d'étoiles binaires à éclipses observées par TESS",
        "object_types": ["Étoiles Binaires"],
        "mag_column": None
    },
    "Exoplanet.eu": {
        "vizier_id": None,  # Exoplanet.eu n'est pas dans Vizier, doit être téléchargé depuis le site web
        "url": "https://exoplanet.eu/catalog/",
        "description": "Catalogue d'exoplanètes du CNRS",
        "object_types": ["Exoplanètes"],
        "mag_column": None
    }
}


class CatalogExtractor:
    """
    Classe pour extraire des données depuis des catalogues astronomiques via Vizier/astroquery.
    """
    
    def __init__(self, output_dir: Optional[Path] = None):
        """
        Initialise l'extracteur de catalogues.
        
        Parameters
        ----------
        output_dir : Path, optional
            Répertoire de sortie pour les fichiers extraits
        """
        if not ASTROQUERY_AVAILABLE:
            raise ImportError("astroquery n'est pas installé. Installez-le avec: pip install astroquery")
        
        self.vizier = Vizier()
        self.vizier.ROW_LIMIT = -1  # Pas de limite par défaut
        
        if output_dir is None:
            output_dir = Path.home() / "catalogues"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def extract_by_region(
        self,
        catalog_name: str,
        center_coord: SkyCoord,
        radius: u.Quantity,
        mag_limit: Optional[float] = None,
        mag_column: Optional[str] = None
    ) -> Table:
        """
        Extrait des données depuis un catalogue dans une région circulaire.
        
        Parameters
        ----------
        catalog_name : str
            Nom du catalogue (doit être dans POPULAR_CATALOGS)
        center_coord : SkyCoord
            Coordonnées du centre de la région
        radius : u.Quantity
            Rayon de recherche (en degrés ou arcmin)
        mag_limit : float, optional
            Limite de magnitude (plus brillant)
        mag_column : str, optional
            Nom de la colonne de magnitude à filtrer
        
        Returns
        -------
        Table
            Table Astropy avec les données extraites
        """
        if catalog_name not in POPULAR_CATALOGS:
            raise ValueError(f"Catalogue '{catalog_name}' non trouvé dans POPULAR_CATALOGS")
        
        catalog_info = POPULAR_CATALOGS[catalog_name]
        vizier_id = catalog_info.get("vizier_id")
        
        if vizier_id is None:
            raise ValueError(
                f"Le catalogue '{catalog_name}' n'est pas disponible via Vizier.\n"
                f"URL: {catalog_info.get('url', 'N/A')}"
            )
        
        # Convertir le rayon en degrés si nécessaire
        if radius.unit == u.arcmin:
            radius_deg = radius.to(u.deg).value
        else:
            radius_deg = radius.to(u.deg).value
        
        logger.info(f"Extraction depuis {catalog_name} (Vizier: {vizier_id})")
        logger.info(f"Centre: {center_coord}, Rayon: {radius_deg:.4f} deg")
        
        # Requête Vizier
        try:
            tables = self.vizier.query_region(
                center_coord,
                radius=radius_deg * u.deg,
                catalog=vizier_id
            )
            
            if len(tables) == 0:
                logger.warning("Aucune table retournée par Vizier")
                return Table()
            
            # Prendre la première table (généralement la bonne)
            table = tables[0]
            
            logger.info(f"{len(table)} objets trouvés")
            
            # Filtrer par magnitude si demandé
            if mag_limit is not None:
                if mag_column is None:
                    mag_column = catalog_info.get("mag_column")
                
                if mag_column and mag_column in table.colnames:
                    # Filtrer : magnitude <= mag_limit (plus brillant)
                    mask = table[mag_column] <= mag_limit
                    table = table[mask]
                    logger.info(f"{len(table)} objets après filtrage magnitude <= {mag_limit}")
                else:
                    logger.warning(f"Colonne de magnitude '{mag_column}' non trouvée dans le catalogue")
            
            return table
            
        except Exception as e:
            logger.error(f"Erreur lors de l'extraction depuis Vizier: {e}", exc_info=True)
            raise
    
    def extract_by_criteria(
        self,
        catalog_name: str,
        ra_min: Optional[float] = None,
        ra_max: Optional[float] = None,
        dec_min: Optional[float] = None,
        dec_max: Optional[float] = None,
        mag_min: Optional[float] = None,
        mag_max: Optional[float] = None,
        mag_column: Optional[str] = None
    ) -> Table:
        """
        Extrait des données depuis un catalogue selon des critères de boîte rectangulaire.
        
        Parameters
        ----------
        catalog_name : str
            Nom du catalogue
        ra_min : float, optional
            RA minimum en degrés
        ra_max : float, optional
            RA maximum en degrés
        dec_min : float, optional
            DEC minimum en degrés
        dec_max : float, optional
            DEC maximum en degrés
        mag_min : float, optional
            Magnitude minimum (plus faible)
        mag_max : float, optional
            Magnitude maximum (plus brillant)
        mag_column : str, optional
            Nom de la colonne de magnitude
        
        Returns
        -------
        Table
            Table Astropy avec les données extraites
        """
        if catalog_name not in POPULAR_CATALOGS:
            raise ValueError(f"Catalogue '{catalog_name}' non trouvé dans POPULAR_CATALOGS")
        
        catalog_info = POPULAR_CATALOGS[catalog_name]
        vizier_id = catalog_info.get("vizier_id")
        
        if vizier_id is None:
            raise ValueError(
                f"Le catalogue '{catalog_name}' n'est pas disponible via Vizier.\n"
                f"URL: {catalog_info.get('url', 'N/A')}"
            )
        
        # Construire les contraintes
        constraints = {}
        
        if ra_min is not None or ra_max is not None:
            if ra_min is not None and ra_max is not None:
                constraints["RA"] = f"{ra_min}..{ra_max}"
            elif ra_min is not None:
                constraints["RA"] = f">{ra_min}"
            elif ra_max is not None:
                constraints["RA"] = f"<{ra_max}"
        
        if dec_min is not None or dec_max is not None:
            if dec_min is not None and dec_max is not None:
                constraints["DE"] = f"{dec_min}..{dec_max}"
            elif dec_min is not None:
                constraints["DE"] = f">{dec_min}"
            elif dec_max is not None:
                constraints["DE"] = f"<{dec_max}"
        
        if not constraints:
            raise ValueError("Au moins un critère (RA ou DEC) doit être spécifié")
        
        logger.info(f"Extraction depuis {catalog_name} avec contraintes: {constraints}")
        
        try:
            # Requête Vizier avec contraintes
            tables = self.vizier.query_constraints(catalog=vizier_id, **constraints)
            
            if len(tables) == 0:
                logger.warning("Aucune table retournée par Vizier")
                return Table()
            
            table = tables[0]
            logger.info(f"{len(table)} objets trouvés")
            
            # Filtrer par magnitude si demandé
            if mag_column is None:
                mag_column = catalog_info.get("mag_column")
            
            if mag_column and mag_column in table.colnames:
                if mag_min is not None:
                    mask_min = table[mag_column] >= mag_min
                    table = table[mask_min]
                
                if mag_max is not None:
                    mask_max = table[mag_column] <= mag_max
                    table = table[mask_max]
                
                if mag_min is not None or mag_max is not None:
                    logger.info(f"{len(table)} objets après filtrage magnitude")
            
            return table
            
        except Exception as e:
            logger.error(f"Erreur lors de l'extraction depuis Vizier: {e}", exc_info=True)
            raise
    
    def reduce_to_essential_columns(self, table: Table) -> Table:
        """
        Réduit la table aux colonnes utiles pour les catalogues d'étoiles (comme le texte explicatif).
        Conservées : ID (source_id), ra, dec, parallax, phot_g_mean_mag, pmra, pmdec.
        Toutes les colonnes d'erreur et surnuméraires sont supprimées.
        """
        if len(table) == 0:
            return table
        colnames_lower = {c.strip().lower(): c for c in table.colnames}
        
        def pick(*candidates: str):
            for name in candidates:
                low = name.lower()
                if low in colnames_lower:
                    return colnames_lower[low]
                for k, orig in colnames_lower.items():
                    if k == low or k.replace(" ", "_") == low or k.replace("-", "_") == low:
                        return orig
            return None
        
        # Colonnes standard à garder (nom de sortie -> candidats possibles dans le catalogue)
        mapping = [
            ("ra", ["ra", "RA_ICRS", "ra_icrs", "RA"]),
            ("dec", ["dec", "DE_ICRS", "dec_icrs", "DEC", "Dec"]),
            ("parallax", ["parallax", "Plx", "parallax_mas", "PM"]),
            ("source_id", ["source_id", "Source", "designation", "id", "ID"]),
            ("pmra", ["pmra", "PMRA", "pmRA"]),
            ("pmdec", ["pmdec", "PMDec", "pmDE"]),
            ("phot_g_mean_mag", ["phot_g_mean_mag", "Gmag", "G"]),
        ]
        keep_cols = []
        rename_map = {}
        for out_name, candidates in mapping:
            col = pick(*candidates)
            if col and col not in keep_cols:
                keep_cols.append(col)
                if col != out_name:
                    rename_map[col] = out_name
        if not keep_cols:
            logger.warning("Aucune colonne essentielle trouvée, table inchangée")
            return table
        out = table[keep_cols].copy()
        if rename_map:
            out.rename_columns(list(rename_map.keys()), list(rename_map.values()))
        logger.info(f"Colonnes conservées : {list(out.colnames)} (toutes les autres supprimées)")
        return out
    
    def save_table(
        self,
        table: Table,
        filename: str,
        format: str = "csv"
    ) -> Path:
        """
        Sauvegarde une table dans un fichier.
        
        Parameters
        ----------
        table : Table
            Table Astropy à sauvegarder
        filename : str
            Nom du fichier (sans extension)
        format : str
            Format de sortie : "csv", "csv.gz", "fits", "votable"
        
        Returns
        -------
        Path
            Chemin du fichier sauvegardé
        """
        if format.lower() == "csv":
            output_path = self.output_dir / f"{filename}.csv"
            table.write(str(output_path), format="ascii.csv", overwrite=True)
        elif format.lower() == "csv.gz":
            output_path = self.output_dir / f"{filename}.csv.gz"
            temp_csv = self.output_dir / f"{filename}.csv.tmp"
            table.write(str(temp_csv), format="ascii.csv", overwrite=True)
            try:
                with open(temp_csv, "rb") as f_in:
                    with gzip.open(str(output_path), "wb") as f_out:
                        f_out.writelines(f_in)
            finally:
                if temp_csv.exists():
                    os.remove(temp_csv)
        elif format.lower() == "fits":
            output_path = self.output_dir / f"{filename}.fits"
            table.write(str(output_path), format="fits", overwrite=True)
        elif format.lower() == "votable":
            output_path = self.output_dir / f"{filename}.xml"
            table.write(str(output_path), format="votable", overwrite=True)
        else:
            raise ValueError(f"Format non supporté: {format}. Utilisez 'csv', 'csv.gz', 'fits' ou 'votable'")
        
        logger.info(f"Table sauvegardée: {output_path}")
        return output_path
