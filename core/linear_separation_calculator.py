# core/linear_separation_calculator.py
"""
Module pour calculer les séparations linéaires entre étoiles dans les catalogues Gaia DR3.
Méthode basée sur Laurent (2022) pour identifier les systèmes binaires physiques.
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u

logger = logging.getLogger(__name__)


class LinearSeparationCalculator:
    """
    Calculateur de séparation linéaire pour identifier les systèmes binaires physiques.
    
    Calcule la séparation linéaire (SL) entre paires d'étoiles en utilisant :
    - La séparation angulaire (en arcsec)
    - La parallaxe Gaia (en mas)
    - La formule : SL (pc) = séparation angulaire (arcsec) / parallaxe (arcsec)
    
    Les couples avec SL < seuil (par défaut 10 pc) sont considérés comme physiques.
    """
    
    def __init__(self):
        """Initialise le calculateur de séparation linéaire."""
        self.pc_to_au = 206265.0  # 1 pc = 206265 AU
        self.mas_to_arcsec = 0.001  # 1 mas = 0.001 arcsec
    
    def calculate_angular_separation(
        self, 
        ra1: np.ndarray, 
        dec1: np.ndarray, 
        ra2: np.ndarray, 
        dec2: np.ndarray
    ) -> np.ndarray:
        """
        Calcule la séparation angulaire entre deux ensembles de coordonnées.
        
        Parameters
        ----------
        ra1, dec1 : np.ndarray
            Coordonnées du premier ensemble (en degrés)
        ra2, dec2 : np.ndarray
            Coordonnées du second ensemble (en degrés)
        
        Returns
        -------
        np.ndarray
            Séparations angulaires en arcsec
        """
        coords1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
        coords2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)
        separations = coords1.separation(coords2)
        return separations.arcsec
    
    def calculate_linear_separation(
        self,
        angular_sep_arcsec: np.ndarray,
        parallax_mas: np.ndarray
    ) -> np.ndarray:
        """
        Calcule la séparation linéaire en parsecs.
        
        Parameters
        ----------
        angular_sep_arcsec : np.ndarray
            Séparation angulaire en arcsec
        parallax_mas : np.ndarray
            Parallaxe en milliarcsecondes (mas)
        
        Returns
        -------
        np.ndarray
            Séparations linéaires en parsecs
        """
        # Convertir parallaxe en arcsec
        parallax_arcsec = parallax_mas * self.mas_to_arcsec
        
        # Éviter division par zéro
        valid_mask = parallax_arcsec > 0
        
        linear_sep = np.full_like(angular_sep_arcsec, np.nan)
        linear_sep[valid_mask] = angular_sep_arcsec[valid_mask] / parallax_arcsec[valid_mask]
        
        return linear_sep
    
    def find_pairs(
        self,
        df: pd.DataFrame,
        max_angular_separation_arcsec: float = 60.0,
        min_angular_separation_arcsec: float = 0.01,
        parallax_column: str = 'parallax'
    ) -> List[Dict]:
        """
        Trouve les paires d'étoiles à l'aide d'une recherche spatiale (search_around_sky),
        en O(n log n) au lieu de O(n²), pour permettre le traitement de grands catalogues.
        """
        if len(df) < 2:
            return []
        
        required_cols = ['ra', 'dec', parallax_column]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Colonnes manquantes dans le DataFrame : {missing_cols}")
        
        ra = np.asarray(df['ra'].values, dtype=float)
        dec = np.asarray(df['dec'].values, dtype=float)
        parallax = np.asarray(df[parallax_column].values, dtype=float)
        
        valid_parallax_mask = ~np.isnan(parallax) & (parallax > 0)
        if not np.any(valid_parallax_mask):
            logger.warning("Aucune étoile avec parallaxe valide trouvée")
            return []
        
        ra_valid = ra[valid_parallax_mask]
        dec_valid = dec[valid_parallax_mask]
        parallax_valid = parallax[valid_parallax_mask]
        indices_valid = np.where(valid_parallax_mask)[0]
        n = len(ra_valid)
        
        logger.info(f"Recherche de paires parmi {n} étoiles (recherche spatiale rapide)...")
        
        coords = SkyCoord(ra=ra_valid * u.deg, dec=dec_valid * u.deg)
        seplimit = max_angular_separation_arcsec * u.arcsec
        result = search_around_sky(coords, coords, seplimit)
        # Astropy 6+ : CoordinateSearchResult (indices_to_first_set, indices_to_second_set, angular_separation, …)
        # Anciennes versions : idx_catalog0 / idx_catalog1, sep2d ; accès tuple (0,1,2) en secours
        idx1 = getattr(result, "indices_to_first_set", None)
        idx2 = getattr(result, "indices_to_second_set", None)
        sep_angle = getattr(result, "angular_separation", None)
        if idx1 is None:
            idx1 = getattr(result, "idx_catalog0", None)
        if idx2 is None:
            idx2 = getattr(result, "idx_catalog1", None)
        if sep_angle is None:
            sep_angle = getattr(result, "sep2d", None)
        if (idx1 is None or idx2 is None) and hasattr(result, "__len__") and len(result) >= 3:
            idx1, idx2 = result[0], result[1]
        if sep_angle is None and hasattr(result, "__len__") and len(result) >= 3:
            sep_angle = result[2]
        if idx1 is None or idx2 is None:
            raise RuntimeError(
                "CoordinateSearchResult sans idx1/idx2 (ni par nom ni par indice). "
                "Mettez à jour astropy ou adaptez linear_separation_calculator."
            )
        if sep_angle is None:
            raise RuntimeError(
                "CoordinateSearchResult sans séparation angulaire (sep2d/angular_separation ou result[2]). "
                "Mettez à jour astropy ou adaptez linear_separation_calculator."
            )
        
        # Garder seulement i < j (éviter doublons et auto-paires)
        mask = idx1 < idx2
        idx1, idx2 = idx1[mask], idx2[mask]
        # .to(u.arcsec).value fonctionne pour Angle ou Quantity
        sep_arcsec = sep_angle[mask].to(u.arcsec).value
        
        # Filtrer par séparation minimale
        mask_min = sep_arcsec >= min_angular_separation_arcsec
        idx1, idx2, sep_arcsec = idx1[mask_min], idx2[mask_min], sep_arcsec[mask_min]
        
        pairs = []
        keep_cols = [c for c in df.columns if c in ('ra', 'dec', parallax_column, 'source_id', 'pmra', 'pmdec', 'phot_g_mean_mag')]
        
        for k in range(len(idx1)):
            i, j = int(idx1[k]), int(idx2[k])
            angular_sep = sep_arcsec[k]
            parallax_mean = (parallax_valid[i] + parallax_valid[j]) / 2.0
            linear_sep = self.calculate_linear_separation(
                np.array([angular_sep]), np.array([parallax_mean])
            )[0]
            ii, jj = indices_valid[i], indices_valid[j]
            pair = {
                'ra1': ra_valid[i], 'dec1': dec_valid[i],
                'ra2': ra_valid[j], 'dec2': dec_valid[j],
                'angular_separation_arcsec': angular_sep,
                'parallax1_mas': parallax_valid[i], 'parallax2_mas': parallax_valid[j],
                'parallax_mean_mas': parallax_mean,
                'linear_separation_pc': linear_sep,
                'index1': ii, 'index2': jj
            }
            for col in keep_cols:
                if col not in ('ra', 'dec', parallax_column):
                    try:
                        pair[f'{col}_1'] = df.iloc[ii][col]
                        pair[f'{col}_2'] = df.iloc[jj][col]
                    except Exception:
                        pass
            pairs.append(pair)
        
        logger.info(f"{len(pairs)} paires trouvées (sép. ang. {min_angular_separation_arcsec}-{max_angular_separation_arcsec} arcsec)")
        return pairs
    
    def filter_physical_pairs(
        self,
        pairs: List[Dict],
        threshold_pc: float = 10.0
    ) -> List[Dict]:
        """
        Filtre les paires pour ne garder que celles avec séparation linéaire < seuil.
        
        Parameters
        ----------
        pairs : List[Dict]
            Liste de dictionnaires de paires
        threshold_pc : float
            Seuil de séparation linéaire en parsecs (par défaut 10.0 pc)
        
        Returns
        -------
        List[Dict]
            Liste filtrée des paires physiques
        """
        physical_pairs = [
            pair for pair in pairs
            if not np.isnan(pair.get('linear_separation_pc', np.nan))
            and pair['linear_separation_pc'] < threshold_pc
        ]
        
        logger.info(f"{len(physical_pairs)} paires physiques (SL < {threshold_pc} pc) sur {len(pairs)} paires totales")
        
        return physical_pairs
    
    def analyze_gaia_csv_file(
        self,
        csv_path: Path,
        max_angular_separation_arcsec: float = 60.0,
        min_angular_separation_arcsec: float = 0.01,
        threshold_pc: float = 10.0,
        parallax_column: str = 'parallax'
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Analyse un fichier CSV Gaia DR3 et trouve les paires d'étoiles.
        
        Parameters
        ----------
        csv_path : Path | str
            Chemin vers le fichier Gaia DR3 (format CSV ou CSV.GZ)
        max_angular_separation_arcsec : float
            Séparation angulaire maximale en arcsec
        min_angular_separation_arcsec : float
            Séparation angulaire minimale en arcsec
        threshold_pc : float
            Seuil de séparation linéaire en parsecs pour les couples physiques
        parallax_column : str
            Nom de la colonne contenant la parallaxe
        
        Returns
        -------
        Tuple[List[Dict], List[Dict]]
            (toutes les paires, paires physiques uniquement)
        """
        csv_path = Path(csv_path)
        logger.info(f"Lecture du fichier CSV : {csv_path}")
        
        # Formats admis : .csv et .csv.gz (compression détectée automatiquement par pandas)
        read_kw = {}
        if str(csv_path).lower().endswith(".gz"):
            read_kw["compression"] = "infer"
        try:
            df = pd.read_csv(csv_path, **read_kw)
            if len(df.columns) <= 1:
                try:
                    df_semi = pd.read_csv(csv_path, sep=";", **read_kw)
                    if len(df_semi.columns) > 1:
                        df = df_semi
                except Exception:
                    pass
                if len(df.columns) <= 1:
                    try:
                        df_tab = pd.read_csv(csv_path, sep="\t", **read_kw)
                        if len(df_tab.columns) > 1:
                            df = df_tab
                    except Exception:
                        pass
        except Exception as e:
            raise ValueError(f"Erreur lors de la lecture du fichier CSV/CSV.GZ : {e}")
        
        if df.empty:
            raise ValueError("Le fichier CSV est vide.")
        
        # En-têtes : BOM UTF-8 et espaces (exports Excel / Vizier)
        df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
        
        logger.info(f"Fichier lu : {len(df)} lignes, colonnes : {list(df.columns)}")
        
        # Normaliser les noms de colonnes (casse et variantes Gaia / Vizier)
        col_lower = {c.strip().lower(): c for c in df.columns}
        col_orig = {c: c for c in df.columns}
        
        def get_standard_col(*candidates: str) -> Optional[str]:
            for name in candidates:
                if name in df.columns:
                    return name
                low = name.lower()
                for k, orig in col_lower.items():
                    if k == low or k.replace(" ", "") == low or orig.lower() == low:
                        return orig
            return None
        
        # RA : Gaia, Vizier (_RAJ2000), TOPCAT, etc.
        ra_col = get_standard_col(
            "ra", "RA_ICRS", "ra_icrs", "RA",
            "_RAJ2000", "raj2000", "RAJ2000", "RAdeg", "ra_deg",
        )
        if ra_col:
            df["ra"] = pd.to_numeric(df[ra_col], errors="coerce")
        
        # DEC
        dec_col = get_standard_col(
            "dec", "DE_ICRS", "dec_icrs", "DEC", "Dec",
            "_DEJ2000", "dej2000", "DEJ2000", "DEdeg", "de_deg",
        )
        if dec_col:
            df["dec"] = pd.to_numeric(df[dec_col], errors="coerce")
        
        missing_cols = [c for c in ["ra", "dec"] if c not in df.columns]
        if missing_cols:
            raise ValueError(
                f"Colonnes manquantes : {missing_cols}. "
                f"Colonnes disponibles : {list(df.columns)}. "
                f"Attendu : ra (ou RA_ICRS, _RAJ2000), dec (ou DE_ICRS, _DEJ2000), parallax (ou Plx)."
            )
        
        # Parallaxe : parallax, Plx, parallax_mas, PM (en mas pour Gaia)
        if parallax_column not in df.columns:
            # Ne pas utiliser « PM » : ambigu avec mouvement propre sur certains CSV
            parallax_col = get_standard_col("parallax", "Plx", "parallax_mas", "Parallax")
            if parallax_col:
                parallax_column = parallax_col
            else:
                raise ValueError(
                    f"Colonne parallaxe non trouvée. Colonnes disponibles : {list(df.columns)}. "
                    f"Attendu : parallax ou Plx (en milliarcsec pour Gaia)."
                )
        
        # Colonne parallaxe numérique (Gaia en mas)
        df["parallax"] = pd.to_numeric(df[parallax_column], errors="coerce")
        n_valid = df["parallax"].notna() & (df["parallax"] > 0)
        if n_valid.sum() == 0:
            raise ValueError(
                f"Aucune ligne avec parallaxe valide (positive). "
                f"Colonne utilisée : '{parallax_column}'. Vérifiez que les valeurs sont en milliarcsec (mas)."
            )
        logger.info(f"Étoiles avec parallaxe valide : {n_valid.sum()} / {len(df)}")
        
        # Ne conserver que : ID, coordonnées (sans erreurs), magnitudes, parallaxes, vitesses propres
        id_col = get_standard_col("source_id", "Source", "source", "id", "ID")
        pmra_col = get_standard_col("pmra", "PMRA")
        pmdec_col = get_standard_col("pmdec", "PMDec")
        mag_col = get_standard_col("phot_g_mean_mag", "Gmag", "G")
        cols_to_keep = ["ra", "dec", "parallax"]
        if id_col and id_col in df.columns:
            df["source_id"] = df[id_col].astype(str)
            cols_to_keep.append("source_id")
        if pmra_col and pmra_col in df.columns:
            df["pmra"] = pd.to_numeric(df[pmra_col], errors="coerce")
            cols_to_keep.append("pmra")
        if pmdec_col and pmdec_col in df.columns:
            df["pmdec"] = pd.to_numeric(df[pmdec_col], errors="coerce")
            cols_to_keep.append("pmdec")
        if mag_col and mag_col in df.columns:
            df["phot_g_mean_mag"] = pd.to_numeric(df[mag_col], errors="coerce")
            cols_to_keep.append("phot_g_mean_mag")
        df = df[[c for c in cols_to_keep if c in df.columns]].copy()
        logger.info(f"Colonnes conservées : {list(df.columns)} (toutes les autres supprimées)")
        
        # Supprimer les lignes avec ra/dec manquants
        df = df.dropna(subset=["ra", "dec"]).copy()
        if len(df) < 2:
            raise ValueError(
                f"Pas assez d'étoiles avec ra, dec et parallaxe valides (au moins 2 requises, trouvé {len(df)})."
            )
        
        # Trouver toutes les paires (on utilise la colonne normalisée 'parallax')
        all_pairs = self.find_pairs(
            df,
            max_angular_separation_arcsec=max_angular_separation_arcsec,
            min_angular_separation_arcsec=min_angular_separation_arcsec,
            parallax_column="parallax"
        )
        
        # Filtrer les paires physiques
        physical_pairs = self.filter_physical_pairs(all_pairs, threshold_pc=threshold_pc)
        
        return all_pairs, physical_pairs
