# gui/easy_lucky_imaging_tab.py
"""
Onglet Easy Lucky Imaging - Techniques REDUC
Basé sur REDUC (http://www.astrosurf.com/hfosaf/reduc/tutoriel.htm)
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import logging
from pathlib import Path
import threading
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import minimize
from scipy.ndimage import center_of_mass, maximum_filter, gaussian_filter

logger = logging.getLogger(__name__)

# Import du module de réduction d'images (techniques REDUC)
try:
    from core.binary_star_reduction import BinaryStarReduction
    REDUCTION_AVAILABLE = True
except ImportError as e:
    REDUCTION_AVAILABLE = False
    logger.warning(f"Module de réduction non disponible: {e}")


def get_image_orientation(wcs, center_x, center_y):
    """
    Calcule l'orientation réelle du Nord et de l'Est sur le capteur
    en utilisant la matrice locale du WCS.
    
    Cette méthode utilise proj_plane_pixel_scales pour obtenir les échelles
    exactes et projette des points célestes sur le plan image pour déterminer
    les directions réelles du Nord et de l'Est.
    
    Parameters
    ----------
    wcs : WCS
        Objet WCS valide
    center_x, center_y : float
        Coordonnées pixel du centre de l'image
        
    Returns
    -------
    dict
        Dictionnaire avec:
        - 'north_angle': angle du Nord en degrés (0° = droite, 90° = haut)
        - 'east_angle': angle de l'Est en degrés
        - 'rotation': rotation du capteur (angle entre axe Y image et Nord)
        - 'parity': "Standard" ou "Inversée (Miroir/Sud)"
    """
    # 1. Obtenir les échelles de pixels (degrés par pixel)
    pixel_scales = proj_plane_pixel_scales(wcs)
    
    # 2. Point central en coordonnées célestes
    c_center = wcs.pixel_to_world(center_x, center_y)
    
    # 3. Projeter un point vers le Nord (Dec+)
    # On utilise l'échelle pixel en Dec (pixel_scales[1])
    c_north = SkyCoord(
        ra=c_center.ra,
        dec=c_center.dec + (pixel_scales[1] * u.deg),
        frame='icrs'
    )
    x_n, y_n = wcs.world_to_pixel(c_north)
    
    # 4. Projeter un point vers l'Est (RA+)
    # On utilise l'échelle pixel en RA (pixel_scales[0])
    # Correction cos(dec) pour la conversion RA -> distance angulaire
    cos_dec = np.cos(np.radians(c_center.dec.deg))
    if cos_dec < 1e-6:
        cos_dec = 1e-6
    delta_ra_deg = pixel_scales[0] / cos_dec
    
    c_east = SkyCoord(
        ra=c_center.ra + (delta_ra_deg * u.deg),
        dec=c_center.dec,
        frame='icrs'
    )
    x_e, y_e = wcs.world_to_pixel(c_east)
    
    # 5. Calcul des angles sur le capteur (en degrés, convention mathématique)
    # Angle 0 = Axe X (droite), 90 = Axe Y (haut)
    angle_north = np.degrees(np.arctan2(y_n - center_y, x_n - center_x))
    angle_east = np.degrees(np.arctan2(y_e - center_y, x_e - center_x))
    
    # Normaliser dans [0, 360)
    if angle_north < 0:
        angle_north += 360.0
    if angle_east < 0:
        angle_east += 360.0
    
    # 6. Calcul de la déviation (Rotation du capteur)
    # C'est l'angle entre l'axe Y de l'image (haut) et le Nord Céleste
    sensor_rotation = (90.0 - angle_north + 360.0) % 360.0
    
    # 7. Vérification de la parité (Mirroring / Hémisphère)
    # On regarde si l'Est est à +90° ou -90° du Nord
    diff = (angle_east - angle_north + 360) % 360
    parity = "Standard" if diff > 180 else "Inversée (Miroir/Sud)"
    
    return {
        'north_angle': angle_north,
        'east_angle': angle_east,
        'rotation': sensor_rotation,
        'parity': parity
    }


class EasyLuckyImagingTab(ttk.Frame):
    """
    Onglet pour le traitement d'images d'étoiles binaires avec techniques REDUC
    """
    
    def __init__(self, parent_notebook, base_dir=None):
        super().__init__(parent_notebook, padding=10)
        
        if base_dir is None:
            self.base_dir = Path.home()
        else:
            self.base_dir = Path(base_dir)
        
        # Processeur de réduction
        if REDUCTION_AVAILABLE:
            self.reducer = BinaryStarReduction()
        else:
            self.reducer = None
            if REDUCTION_AVAILABLE is False:
                logger.error("[INIT] Module de réduction non disponible - Le module core/binary_star_reduction.py doit être disponible.")
        
        # Liste de travail (images retenues)
        self.work_list = []  # Liste de tuples (path, metrics)
        self.work_list_file = None  # Chemin du fichier de liste sauvegardé
        
        self.create_widgets()
        
        # Barre de progression globale
        progress_container = ttk.Frame(self)
        progress_container.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
        self.progress_label = ttk.Label(progress_container, text="", font=("Helvetica", 9))
        self.progress_label.pack(anchor="w")
        self.progress = ttk.Progressbar(progress_container, mode="indeterminate")
        self.progress.pack(side=tk.BOTTOM, fill=tk.X)
    
    def create_widgets(self):
        """Crée l'interface utilisateur"""
        
        # En-tête
        header_frame = ttk.Frame(self, padding=10)
        header_frame.pack(fill="x")
        
        title_label = ttk.Label(
            header_frame,
            text="Techniques de Réduction REDUC",
            font=("Helvetica", 14, "bold")
        )
        title_label.pack()
        
        info_label = ttk.Label(
            header_frame,
            text="Basé sur REDUC (http://www.astrosurf.com/hfosaf/reduc/tutoriel.htm)",
            font=("Helvetica", 8),
            foreground="gray"
        )
        info_label.pack()
        
        # Notebook pour les différentes fonctionnalités
        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Onglet 1: Lucky Imaging (BestOf + ELI fusionnés)
        lucky_frame = ttk.Frame(notebook, padding=10)
        notebook.add(lucky_frame, text="Lucky Imaging (BestOf + ELI)")
        self.create_lucky_imaging_tab(lucky_frame)
        
        # Onglet 2: Mesure de séparation
        measure_frame = ttk.Frame(notebook, padding=10)
        notebook.add(measure_frame, text="Mesure de séparation")
        self.create_measure_tab(measure_frame)
    
    # === Helpers spécifiques au cadre "Recherche de Binaires (Laurent)" dans Lucky Imaging ===
    def _browse_gaia_file_for_laurent(self):
        """Choix du fichier Gaia DR3 (CSV/CSV.GZ) pour la méthode de Laurent."""
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("CSV ou GZ", "*.csv *.csv.gz *.gz"), ("Tous les fichiers", "*.*")]
        )
        if path:
            self.gaia_file_var.set(path)
    
    def _browse_binaries_output_dir_for_laurent(self):
        """Choix du répertoire de sortie pour les couples binaires (méthode Laurent)."""
        directory = filedialog.askdirectory(initialdir=self.base_dir)
        if directory:
            self.binaries_output_dir_var.set(directory)
    
    def _browse_nina_csv_for_laurent(self):
        """Choix du fichier CSV à convertir en JSON NINA (méthode Laurent)."""
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("CSV", "*.csv"), ("Tous les fichiers", "*.*")]
        )
        if path:
            self.nina_csv_file_var.set(path)
    
    def _browse_nina_output_dir_for_laurent(self):
        """Choix du répertoire de sortie pour les fichiers JSON NINA."""
        directory = filedialog.askdirectory(initialdir=self.base_dir)
        if directory:
            self.nina_output_dir_var.set(directory)
    
    def _start_nina_conversion_from_lucky(self):
        """Lance la conversion CSV → JSON NINA directement depuis l'onglet Lucky Imaging."""
        import threading
        threading.Thread(target=self._run_nina_conversion_from_lucky, daemon=True).start()
    
    def _run_nina_conversion_from_lucky(self):
        """Implémentation locale de la conversion CSV → JSON NINA (inspirée de CataloguesTab)."""
        try:
            import subprocess
            import sys
            from pathlib import Path
        except Exception as e:
            logger.error(f"[LuckyImaging] Erreur import conversion NINA: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur d'initialisation de la conversion NINA : {e}")
            return
        
        try:
            csv_file = self.nina_csv_file_var.get().strip()
            if not csv_file:
                messagebox.showerror("Erreur", "Veuillez sélectionner un fichier CSV à convertir.")
                return
            
            csv_path = Path(csv_file)
            if not csv_path.exists():
                messagebox.showerror("Erreur", f"Fichier non trouvé : {csv_file}")
                return
            
            output_dir = self.nina_output_dir_var.get().strip()
            if not output_dir:
                output_dir = str(self.base_dir / "nina_json")
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            script_path = Path(__file__).parent.parent / "utils" / "convert_binaries_to_nina.py"
            if not script_path.exists():
                messagebox.showerror("Erreur", f"Script de conversion non trouvé : {script_path}")
                return
            
            python_exe = sys.executable
            cmd = [
                python_exe,
                str(script_path),
                str(csv_path),
                str(output_dir)
            ]
            
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding='utf-8',
                errors='replace'
            )
            
            output_lines = []
            for line in process.stdout:
                line = line.strip()
                if line:
                    output_lines.append(line)
                    logger.info(f"[NINA] {line}")
            
            return_code = process.wait()
            
            if return_code == 0:
                json_files = list(output_dir.glob("*.json"))
                count = len(json_files)
                messagebox.showinfo(
                    "Succès",
                    f"Conversion terminée avec succès !\n\n"
                    f"{count} fichiers JSON créés\n"
                    f"\nFichiers dans :\n{output_dir}"
                )
            else:
                error_msg = f"Erreur lors de la conversion (code {return_code})"
                if output_lines:
                    logger.error("[NINA] " + "\n".join(output_lines[-10:]))
                messagebox.showerror("Erreur", error_msg)
        
        except Exception as e:
            logger.error(f"[LuckyImaging] Erreur lors de la conversion NINA : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de la conversion NINA : {e}")
    
    def _start_laurent_analysis_from_lucky(self):
        """Lance l'analyse de séparation linéaire directement depuis l'onglet Lucky Imaging."""
        import threading
        threading.Thread(target=self._run_laurent_analysis_from_lucky, daemon=True).start()
    
    def _run_laurent_analysis_from_lucky(self):
        """Implémentation locale de la méthode de Laurent (sans passer par l'onglet Catalogues)."""
        try:
            from core.linear_separation_calculator import LinearSeparationCalculator
            import pandas as pd
            from pathlib import Path
        except Exception as e:
            logger.error(f"[LuckyImaging] Impossible d'importer LinearSeparationCalculator : {e}", exc_info=True)
            messagebox.showerror(
                "Erreur",
                "Module core.linear_separation_calculator introuvable.\n"
                "L'analyse de séparation linéaire (méthode Laurent) n'est pas disponible."
            )
            return
        
        gaia_file = self.gaia_file_var.get().strip()
        if not gaia_file:
            messagebox.showerror("Erreur", "Veuillez sélectionner un fichier CSV Gaia DR3.")
            return
        
        gaia_path = Path(gaia_file)
        if not gaia_path.exists():
            messagebox.showerror("Erreur", f"Fichier non trouvé : {gaia_file}")
            return
        
        # Paramètres
        try:
            threshold_pc = float(self.linear_sep_threshold_var.get())
        except Exception:
            threshold_pc = 10.0
        
        try:
            max_angular_sep = float(self.max_angular_sep_var.get())
        except Exception:
            max_angular_sep = 60.0
        
        # Répertoire de sortie
        output_dir = self.binaries_output_dir_var.get().strip()
        if not output_dir:
            output_dir = str(self.base_dir / "binaries_laurent")
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        
        try:
            calculator = LinearSeparationCalculator()
            all_pairs, physical_pairs = calculator.analyze_gaia_csv_file(
                gaia_path,
                max_angular_separation_arcsec=max_angular_sep,
                threshold_pc=threshold_pc,
            )
        except Exception as e:
            logger.error(f"[LuckyImaging] Erreur analyse Laurent : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de l'analyse du fichier Gaia : {e}")
            return
        
        if not all_pairs:
            messagebox.showwarning(
                "Attention",
                "Aucun couple n'a été trouvé dans le fichier.\n"
                "Essayez d'augmenter la séparation angulaire max ou vérifiez le contenu du CSV."
            )
            return
        
        # Sauvegarde des résultats (similaire à CataloguesTab)
        results_df = pd.DataFrame(all_pairs)
        physical_df = pd.DataFrame(physical_pairs) if physical_pairs else pd.DataFrame()
        
        all_path = output_dir_path / "laurent_all_pairs.csv"
        phys_path = output_dir_path / "laurent_physical_pairs.csv"
        try:
            results_df.to_csv(all_path, index=False)
            if not physical_df.empty:
                physical_df.to_csv(phys_path, index=False)
        except Exception as e:
            logger.error(f"[LuckyImaging] Erreur sauvegarde résultats Laurent : {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de la sauvegarde des résultats : {e}")
            return
        
        message = (
            f"Analyse terminée avec succès.\n\n"
            f"{len(all_pairs)} couples trouvés.\n"
            f"{len(physical_pairs)} couples physiques (SL < {threshold_pc} pc).\n\n"
            f"Résultats enregistrés dans :\n{output_dir_path}"
        )
        messagebox.showinfo("Succès", message)
    
    def create_lucky_imaging_tab(self, parent):
        """Crée l'onglet Lucky Imaging (BestOf + ELI fusionnés)"""
        
        # Colonne principale gauche pour tous les sous-cadres Lucky Imaging
        main_frame = ttk.Frame(parent)
        main_frame.pack(fill="both", expand=True)
        
        # Deux colonnes de même poids pour équilibrer la largeur des cadres
        left_col = ttk.Frame(main_frame)
        right_col = ttk.LabelFrame(main_frame, text="Recherche de Binaires (Laurent 2022)", padding=8)
        
        left_col.pack(side="left", fill="both", expand=True)
        right_col.pack(side="left", fill="both", expand=True, padx=(5, 0), pady=5)
        
        # Chargement d'images
        load_frame = ttk.LabelFrame(left_col, text="1. Chargement des images", padding=10)
        load_frame.pack(fill="x", pady=5)
        
        ttk.Label(load_frame, text="Dossier d'images:").pack(anchor="w", pady=2)
        
        self.lucky_dir_var = tk.StringVar()
        dir_frame = ttk.Frame(load_frame)
        dir_frame.pack(fill="x", pady=5)
        ttk.Entry(dir_frame, textvariable=self.lucky_dir_var, width=50).pack(side="left", fill="x", expand=True)
        ttk.Button(dir_frame, text="📁 Parcourir", command=lambda: self.browse_folder(self.lucky_dir_var)).pack(side="left", padx=(5, 0))
        
        # Section 1: Analyse et tri (BestOf)
        analysis_frame = ttk.LabelFrame(left_col, text="2. Analyse et tri des images (BestOf)", padding=10)
        analysis_frame.pack(fill="x", pady=5)
        
        ttk.Label(analysis_frame, text="Pourcentage d'images à analyser/conserver:").pack(anchor="w", pady=2)
        
        self.lucky_percent_var = tk.StringVar(value="10")
        scale_frame = ttk.Frame(analysis_frame)
        scale_frame.pack(fill="x", pady=5)
        scale = ttk.Scale(scale_frame, from_=5, to=50, orient=tk.HORIZONTAL,
                         variable=self.lucky_percent_var, length=300,
                         command=lambda v: self.lucky_percent_label.config(text=f"{float(v):.0f}%"))
        scale.pack(side="left", fill="x", expand=True)
        self.lucky_percent_label = ttk.Label(scale_frame, text="10%", width=6)
        self.lucky_percent_label.pack(side="left", padx=(10, 0))
        
        ttk.Label(analysis_frame, 
                 text="Les images seront triées par qualité (FWHM, SNR, contraste) et seules les meilleures seront utilisées pour le stacking.",
                 font=("Helvetica", 8),
                 foreground="gray",
                 wraplength=600).pack(anchor="w", pady=(5, 0))
        
        # Bouton d'analyse dans la section Analyse
        analysis_button_frame = ttk.Frame(analysis_frame)
        analysis_button_frame.pack(fill="x", pady=(10, 0))
        
        ttk.Button(
            analysis_button_frame,
            text="🔍 Analyser les images (crée la liste de travail)",
            command=self.run_analysis_only
        ).pack(pady=5)
        
        # Section 3: Résultats de l'analyse et liste de travail
        result_frame = ttk.LabelFrame(left_col, text="3. Résultats de l'analyse et liste de travail", padding=10)
        result_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Zone de texte avec scrollbar pour les résultats
        text_frame = ttk.Frame(result_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(text_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.lucky_result_text = tk.Text(text_frame, height=10, yscrollcommand=scrollbar.set, width=40)
        self.lucky_result_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.lucky_result_text.yview)
        
        # Informations sur la liste de travail
        work_list_frame = ttk.Frame(result_frame)
        work_list_frame.pack(fill="x", pady=(5, 0))
        
        self.work_list_info_label = ttk.Label(
            work_list_frame,
            text="Aucune liste de travail active",
            font=("Helvetica", 9),
            foreground="gray"
        )
        self.work_list_info_label.pack(side="left")
        
        # Boutons de gestion de la liste
        list_buttons_frame = ttk.Frame(result_frame)
        list_buttons_frame.pack(fill="x", pady=5)
        
        ttk.Button(
            list_buttons_frame,
            text="💾 Sauvegarder la liste",
            command=self.save_work_list,
            state="disabled"
        ).pack(side="left", padx=2)
        self.save_list_button = list_buttons_frame.winfo_children()[-1]
        
        ttk.Button(
            list_buttons_frame,
            text="📂 Charger une liste",
            command=self.load_work_list
        ).pack(side="left", padx=2)
        
        ttk.Button(
            list_buttons_frame,
            text="📋 Voir la liste",
            command=self.view_work_list
        ).pack(side="left", padx=2)
        
        self.gaia_file_var = tk.StringVar()
        self.max_angular_sep_var = tk.StringVar(value="60.0")
        self.linear_sep_threshold_var = tk.StringVar(value="10.0")
        self.binaries_output_dir_var = tk.StringVar()
        
        row1 = ttk.Frame(right_col)
        row1.pack(fill="x", pady=2)
        ttk.Label(row1, text="Gaia DR3 (CSV/CSV.GZ):").pack(anchor="w")
        row1b = ttk.Frame(right_col)
        row1b.pack(fill="x", pady=1)
        ttk.Entry(row1b, textvariable=self.gaia_file_var, width=28).pack(side="left", fill="x", expand=True)
        ttk.Button(row1b, text="📁", command=self._browse_gaia_file_for_laurent, width=3).pack(side="left", padx=2)
        
        row2 = ttk.Frame(right_col)
        row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Sep. angulaire max (\"):", width=20).pack(side="left")
        ttk.Entry(row2, textvariable=self.max_angular_sep_var, width=8).pack(side="left", padx=2)
        
        row3 = ttk.Frame(right_col)
        row3.pack(fill="x", pady=2)
        ttk.Label(row3, text="Seuil SL (pc):", width=20).pack(side="left")
        ttk.Entry(row3, textvariable=self.linear_sep_threshold_var, width=8).pack(side="left", padx=2)
        
        # Répertoire de sortie pour la méthode Laurent
        row4 = ttk.Frame(right_col)
        row4.pack(fill="x", pady=4)
        ttk.Label(row4, text="Répertoire de sortie:", width=20).pack(side="left")
        out_entry = ttk.Entry(row4, textvariable=self.binaries_output_dir_var, width=20)
        out_entry.pack(side="left", fill="x", expand=True)
        ttk.Button(
            row4,
            text="📁",
            width=3,
            command=self._browse_binaries_output_dir_for_laurent
        ).pack(side="left", padx=2)
        
        ttk.Button(
            right_col,
            text="🔍 Analyser couples (séparation linéaire)",
            command=self._start_laurent_analysis_from_lucky
        ).pack(pady=(6, 4), fill="x")
        
        # Export vers NINA (section compacte sous la méthode Laurent)
        nina_frame = ttk.LabelFrame(right_col, text="Export vers NINA", padding=6)
        nina_frame.pack(fill="x", pady=(6, 0))
        
        self.nina_csv_file_var = tk.StringVar()
        self.nina_output_dir_var = tk.StringVar()
        
        ttk.Label(nina_frame, text="Fichier CSV à convertir:").pack(anchor="w", pady=(0, 2))
        nina_csv_row = ttk.Frame(nina_frame)
        nina_csv_row.pack(fill="x", pady=1)
        ttk.Entry(nina_csv_row, textvariable=self.nina_csv_file_var, width=24).pack(side="left", fill="x", expand=True)
        ttk.Button(
            nina_csv_row,
            text="📁",
            width=3,
            command=self._browse_nina_csv_for_laurent
        ).pack(side="left", padx=2)
        
        ttk.Label(nina_frame, text="Répertoire de sortie JSON:").pack(anchor="w", pady=(4, 2))
        nina_out_row = ttk.Frame(nina_frame)
        nina_out_row.pack(fill="x", pady=1)
        ttk.Entry(nina_out_row, textvariable=self.nina_output_dir_var, width=24).pack(side="left", fill="x", expand=True)
        ttk.Button(
            nina_out_row,
            text="📁",
            width=3,
            command=self._browse_nina_output_dir_for_laurent
        ).pack(side="left", padx=2)
        
        ttk.Button(
            nina_frame,
            text="📤 Convertir en JSON NINA",
            command=self._start_nina_conversion_from_lucky
        ).pack(pady=(6, 0), fill="x")
        
        # Section 4: Stacking (ELI)
        stacking_frame = ttk.LabelFrame(left_col, text="4. Stacking (ELI - Easy Lucky Imaging)", padding=10)
        stacking_frame.pack(fill="x", pady=5)
        
        # Bouton pour voir l'image de référence
        ref_image_button_frame = ttk.Frame(stacking_frame)
        ref_image_button_frame.pack(fill="x", pady=(0, 10))
        
        ttk.Button(
            ref_image_button_frame,
            text="🖼️ Voir l'image de référence",
            command=self.show_reference_image
        ).pack(side="left", padx=(0, 10))
        
        ttk.Label(
            ref_image_button_frame,
            text="(Affiche l'image avec les meilleurs critères de qualité. Clic pour sélectionner la cible.)",
            font=("Helvetica", 8),
            foreground="gray"
        ).pack(side="left")
        
        ttk.Label(stacking_frame, text="Position de référence pour l'alignement sub-pixel (x, y):").pack(anchor="w", pady=2)
        ref_pos_frame = ttk.Frame(stacking_frame)
        ref_pos_frame.pack(fill="x", pady=5)
        self.ref_x_var = tk.StringVar(value="")
        self.ref_y_var = tk.StringVar(value="")
        ttk.Entry(ref_pos_frame, textvariable=self.ref_x_var, width=15).pack(side="left", padx=(0, 5))
        ttk.Label(ref_pos_frame, text=", ").pack(side="left")
        ttk.Entry(ref_pos_frame, textvariable=self.ref_y_var, width=15).pack(side="left")
        ttk.Label(ref_pos_frame, text="(en pixels - position d'une étoile brillante)", 
                 font=("Helvetica", 8), foreground="gray").pack(side="left", padx=(10, 0))
        
        ttk.Label(stacking_frame, text="Méthode de stacking:").pack(anchor="w", pady=(10, 2))
        self.lucky_method_var = tk.StringVar(value="median")
        method_frame = ttk.Frame(stacking_frame)
        method_frame.pack(fill="x", pady=5)
        for method, label in [("median", "Médiane (recommandé)"), ("mean", "Moyenne"), ("sigma_clip", "Sigma-clip")]:
            ttk.Radiobutton(method_frame, text=label, variable=self.lucky_method_var, value=method).pack(side="left", padx=(0, 20))
        
        # Option de sauvegarde automatique
        save_options_frame = ttk.Frame(stacking_frame)
        save_options_frame.pack(fill="x", pady=(10, 0))
        
        self.auto_save_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            save_options_frame,
            text="Enregistrement automatique dans le dossier source",
            variable=self.auto_save_var
        ).pack(side="left")
        
        ttk.Label(
            save_options_frame,
            text="(sinon, demande de choisir l'emplacement)",
            font=("Helvetica", 8),
            foreground="gray"
        ).pack(side="left", padx=(10, 0))
        
        # Bouton de stacking dans la section ELI
        stacking_button_frame = ttk.Frame(stacking_frame)
        stacking_button_frame.pack(fill="x", pady=(15, 5))
        
        ttk.Button(
            stacking_button_frame,
            text="✨ Créer image empilée (utilise la liste de travail)",
            command=self.run_stacking_from_work_list
        ).pack(pady=5)
        
        ttk.Label(
            stacking_button_frame,
            text="⚠️ Une liste de travail doit être créée ou chargée avant de lancer le stacking",
            font=("Helvetica", 8),
            foreground="orange"
        ).pack()
    
    def create_measure_tab(self, parent):
        """Crée l'onglet de mesure de séparation"""
        
        # Chargement d'image
        load_frame = ttk.LabelFrame(parent, text="1. Chargement de l'image", padding=10)
        load_frame.pack(fill="x", pady=5)
        
        self.measure_image_var = tk.StringVar()
        img_frame = ttk.Frame(load_frame)
        img_frame.pack(fill="x", pady=5)
        ttk.Entry(img_frame, textvariable=self.measure_image_var, width=50).pack(side="left", fill="x", expand=True)
        ttk.Button(img_frame, text="📁 Parcourir", command=self.browse_image_file).pack(side="left", padx=(5, 0))
        
        # Bouton pour visualiser l'image
        visualize_button_frame = ttk.Frame(load_frame)
        visualize_button_frame.pack(fill="x", pady=(5, 0))
        
        ttk.Button(
            visualize_button_frame,
            text="🖼️ Visualiser l'image et sélectionner les étoiles",
            command=self.show_measure_image
        ).pack(side="left", padx=(0, 10))
        
        ttk.Label(
            visualize_button_frame,
            text="(Ouvre une fenêtre interactive. Cliquez d'abord sur étoile 1, puis sur étoile 2.)",
            font=("Helvetica", 8),
            foreground="gray"
        ).pack(side="left")
        
        # Positions des étoiles
        positions_frame = ttk.LabelFrame(parent, text="2. Positions approximatives des étoiles", padding=10)
        positions_frame.pack(fill="x", pady=5)
        
        ttk.Label(positions_frame, text="Étoile 1 (x, y):").pack(anchor="w", pady=2)
        star1_frame = ttk.Frame(positions_frame)
        star1_frame.pack(fill="x", pady=5)
        self.star1_x_var = tk.StringVar(value="")
        self.star1_y_var = tk.StringVar(value="")
        ttk.Entry(star1_frame, textvariable=self.star1_x_var, width=15).pack(side="left", padx=(0, 5))
        ttk.Label(star1_frame, text=", ").pack(side="left")
        ttk.Entry(star1_frame, textvariable=self.star1_y_var, width=15).pack(side="left")
        ttk.Label(star1_frame, text="(en pixels)", font=("Helvetica", 8), foreground="gray").pack(side="left", padx=(10, 0))
        
        ttk.Label(positions_frame, text="Étoile 2 (x, y):").pack(anchor="w", pady=(10, 2))
        star2_frame = ttk.Frame(positions_frame)
        star2_frame.pack(fill="x", pady=5)
        self.star2_x_var = tk.StringVar(value="")
        self.star2_y_var = tk.StringVar(value="")
        ttk.Entry(star2_frame, textvariable=self.star2_x_var, width=15).pack(side="left", padx=(0, 5))
        ttk.Label(star2_frame, text=", ").pack(side="left")
        ttk.Entry(star2_frame, textvariable=self.star2_y_var, width=15).pack(side="left")
        ttk.Label(star2_frame, text="(en pixels)", font=("Helvetica", 8), foreground="gray").pack(side="left", padx=(10, 0))
        
        # Échelle pixel
        ttk.Label(positions_frame, text="Échelle pixel (arcsec/pixel):").pack(anchor="w", pady=(10, 2))
        self.pixel_scale_var = tk.StringVar(value="1.0")
        ttk.Entry(positions_frame, textvariable=self.pixel_scale_var, width=15).pack(anchor="w", pady=5)
        
        # Bouton Analyse
        analyze_button_frame = ttk.Frame(positions_frame)
        analyze_button_frame.pack(fill="x", pady=(10, 0))
        
        ttk.Button(
            analyze_button_frame,
            text="📊 Analyser et générer rapport",
            command=self.analyze_and_generate_report
        ).pack(pady=5)
        
        # Résultats
        result_frame = ttk.LabelFrame(parent, text="3. Résultats", padding=10)
        result_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        self.measure_result_text = tk.Text(result_frame, height=10)
        self.measure_result_text.pack(fill=tk.BOTH, expand=True)
    
    def browse_folder(self, var):
        """Ouvre un dialogue pour sélectionner un dossier"""
        folder = filedialog.askdirectory(initialdir=self.base_dir)
        if folder:
            var.set(folder)
    
    def browse_image_file(self):
        """Ouvre un dialogue pour sélectionner un fichier image"""
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("FITS", "*.fits"), ("Tous les fichiers", "*.*")]
        )
        if path:
            self.measure_image_var.set(path)
    
    def _start_analysis_progress(self, total: int):
        """Passe la barre de progression en mode déterministe pour l'analyse (appelé dans le thread GUI)."""
        self.progress.config(mode="determinate", maximum=total, value=0)
        self.progress_label.config(text=f"Analyse en cours: 0 / {total} images")

    def _update_analysis_progress(self, current: int, total: int):
        """Met à jour la barre et le label de progression (appelé dans le thread GUI)."""
        self.progress.config(value=current)
        self.progress_label.config(text=f"Analyse en cours: {current} / {total} images")

    def _stop_analysis_progress(self, total: int):
        """Remet la barre à 100 % puis repasse en mode indéterminé (appelé dans le thread GUI)."""
        self.progress.config(value=total)
        self.progress_label.config(text="Analyse terminée.")
        self.progress.config(mode="indeterminate")
        self.progress_label.after(2000, lambda: self.progress_label.config(text=""))

    def run_analysis_only(self):
        """Exécute uniquement l'analyse BestOf sans créer l'image empilée"""
        if not REDUCTION_AVAILABLE or self.reducer is None:
            logger.error("[ANALYSE] Module de réduction non disponible")
            return
        
        image_dir = self.lucky_dir_var.get()
        if not image_dir:
            messagebox.showwarning("Attention", "Sélectionnez d'abord un dossier d'images!")
            return
        
        folder = Path(image_dir)
        # Rechercher les fichiers FITS (insensible à la casse)
        fits_files = list(folder.glob("*.fits")) + list(folder.glob("*.FITS"))
        # Dédupliquer en utilisant le chemin absolu normalisé (pour Windows)
        # Utiliser resolve() pour normaliser les chemins et gérer la casse
        seen_paths = set()
        unique_files = []
        for f in fits_files:
            resolved = f.resolve()
            if resolved not in seen_paths:
                seen_paths.add(resolved)
                unique_files.append(f)
        image_files = sorted(unique_files, key=str)
        
        if not image_files:
            messagebox.showwarning("Attention", "Aucune image FITS trouvée dans le dossier!")
            return
        
        total_images = len(image_files)

        def progress_callback(current, total):
            """Appelé depuis le thread worker ; planifie la mise à jour GUI dans le thread principal."""
            self.after(0, lambda c=current, t=total: self._update_analysis_progress(c, t))

        def analysis_task():
            try:
                self.after(0, lambda: self._start_analysis_progress(total_images))
                percent = float(self.lucky_percent_var.get()) / 100.0

                results = self.reducer.bestof_sort(
                    image_files, top_percent=percent, progress_callback=progress_callback
                )

                self.after(0, lambda: self._stop_analysis_progress(total_images))
                
                # Afficher les résultats dans la zone de texte
                self.lucky_result_text.delete(1.0, tk.END)
                result_text = f"Analyse terminée!\n"
                result_text += f"{'='*60}\n\n"
                result_text += f"Total images analysées: {len(image_files)}\n"
                result_text += f"Meilleures images sélectionnées: {len(results)} ({percent*100:.0f}%)\n\n"
                result_text += f"Top {min(20, len(results))} meilleures images:\n"
                result_text += f"{'-'*60}\n"
                
                for i, (path, metrics) in enumerate(results[:20], 1):
                    score = metrics.get('score', 0.0)
                    fwhm = metrics.get('fwhm', 0.0)
                    snr = metrics.get('snr', 0.0)
                    contrast = metrics.get('contrast', 0.0)
                    n_stars = metrics.get('n_stars', 0)
                    
                    result_text += f"\n{i:3d}. {path.name}\n"
                    result_text += f"      Score: {score:.2f}\n"
                    result_text += f"      FWHM: {fwhm:.2f}\" | SNR: {snr:.2f} | Contraste: {contrast:.3f} | Étoiles: {n_stars}\n"
                
                result_text += f"\n{'='*60}\n"
                result_text += f"Pour créer l'image empilée, remplissez les paramètres de stacking ci-dessous\n"
                result_text += f"et cliquez sur 'Analyser et créer image empilée'.\n"
                
                self.lucky_result_text.insert(1.0, result_text)
                
                # Sauvegarder la liste de travail
                self.work_list = results
                self.work_list_file = None  # Réinitialiser le fichier (nouvelle analyse)
                
                # Mettre à jour l'interface
                self.update_work_list_info()
                
                logger.info(f"Analyse: {len(results)} meilleures images sur {len(image_files)}")
                
            except Exception as e:
                self.after(0, lambda: self._stop_analysis_progress(total_images))
                logger.error(f"[ANALYSE] Erreur lors de l'analyse: {e}", exc_info=True)

        threading.Thread(target=analysis_task, daemon=True).start()
    
    def update_work_list_info(self):
        """Met à jour l'affichage des informations sur la liste de travail"""
        if self.work_list:
            n_images = len(self.work_list)
            if self.work_list_file:
                file_name = Path(self.work_list_file).name
                self.work_list_info_label.config(
                    text=f"✅ Liste de travail active: {n_images} images (depuis {file_name})",
                    foreground="green"
                )
            else:
                self.work_list_info_label.config(
                    text=f"✅ Liste de travail active: {n_images} images (analysées)",
                    foreground="green"
                )
            self.save_list_button.config(state="normal")
        else:
            self.work_list_info_label.config(
                text="Aucune liste de travail active",
                foreground="gray"
            )
            self.save_list_button.config(state="disabled")
    
    def save_work_list(self):
        """Sauvegarde la liste de travail dans un fichier JSON"""
        if not self.work_list:
            messagebox.showwarning("Attention", "Aucune liste de travail à sauvegarder!")
            return
        
        # Demander le fichier de sauvegarde
        save_path = filedialog.asksaveasfilename(
            initialdir=self.base_dir,
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("Tous les fichiers", "*.*")],
            title="Sauvegarder la liste de travail"
        )
        
        if not save_path:
            return
        
        try:
            # Préparer les données à sauvegarder
            data = {
                'total_images': len(self.work_list),
                'source_directory': str(self.lucky_dir_var.get()),
                'percent_used': float(self.lucky_percent_var.get()),
                'images': []
            }
            
            for path, metrics in self.work_list:
                data['images'].append({
                    'path': str(path),
                    'filename': path.name,
                    'metrics': {
                        'score': float(metrics.get('score', 0.0)),
                        'fwhm': float(metrics.get('fwhm', 0.0)),
                        'snr': float(metrics.get('snr', 0.0)),
                        'contrast': float(metrics.get('contrast', 0.0)),
                        'n_stars': int(metrics.get('n_stars', 0))
                    }
                })
            
            # Sauvegarder en JSON
            with open(save_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            
            self.work_list_file = save_path
            self.update_work_list_info()
            
            messagebox.showinfo("Succès", f"Liste de travail sauvegardée!\n\n{save_path}\n\n{len(self.work_list)} images")
            logger.info(f"Liste de travail sauvegardée: {save_path} ({len(self.work_list)} images)")
            
        except Exception as e:
            logger.error(f"[SAUVEGARDE] Erreur lors de la sauvegarde: {e}", exc_info=True)
            messagebox.showwarning("Attention", f"Erreur lors de la sauvegarde. Vérifiez les logs pour plus de détails.")
    
    def load_work_list(self):
        """Charge une liste de travail depuis un fichier JSON"""
        load_path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("JSON", "*.json"), ("Tous les fichiers", "*.*")],
            title="Charger une liste de travail"
        )
        
        if not load_path:
            return
        
        try:
            with open(load_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # Vérifier la structure
            if 'images' not in data:
                raise ValueError("Format de fichier invalide: 'images' manquant")
            
            # Reconstruire la liste de travail (en évitant les doublons)
            self.work_list = []
            seen_paths = set()
            for item in data['images']:
                path = Path(item['path'])
                # Vérifier que le fichier existe
                if not path.exists():
                    logger.warning(f"Image non trouvée: {path}")
                    continue
                
                # Dédupliquer en utilisant le chemin résolu
                resolved = path.resolve()
                if resolved in seen_paths:
                    logger.debug(f"Fichier dupliqué ignoré: {path}")
                    continue
                seen_paths.add(resolved)
                
                metrics = item.get('metrics', {})
                self.work_list.append((path, metrics))
            
            if not self.work_list:
                messagebox.showwarning("Attention", "Aucune image valide trouvée dans la liste!")
                return
            
            self.work_list_file = load_path
            
            # Restaurer les paramètres si disponibles
            if 'source_directory' in data:
                self.lucky_dir_var.set(data['source_directory'])
            if 'percent_used' in data:
                self.lucky_percent_var.set(str(data['percent_used']))
                self.lucky_percent_label.config(text=f"{data['percent_used']:.0f}%")
            
            # Afficher les informations
            self.display_work_list_info(data)
            
            # Mettre à jour l'interface
            self.update_work_list_info()
            
            messagebox.showinfo("Succès", f"Liste de travail chargée!\n\n{len(self.work_list)} images")
            logger.info(f"Liste de travail chargée: {load_path} ({len(self.work_list)} images)")
            
        except Exception as e:
            logger.error(f"[CHARGEMENT] Erreur lors du chargement: {e}", exc_info=True)
            messagebox.showwarning("Attention", f"Erreur lors du chargement. Vérifiez les logs pour plus de détails.")
    
    def display_work_list_info(self, data=None):
        """Affiche les informations de la liste de travail dans la zone de texte"""
        if not self.work_list:
            return
        
        if data is None:
            # Reconstruire les infos depuis la liste
            result_text = f"Liste de travail chargée\n"
            result_text += f"{'='*60}\n\n"
            result_text += f"Nombre d'images: {len(self.work_list)}\n"
            if self.work_list_file:
                result_text += f"Source: {Path(self.work_list_file).name}\n"
            result_text += f"\n{'='*60}\n\n"
            result_text += f"Top {min(20, len(self.work_list))} images:\n"
            result_text += f"{'-'*60}\n"
            
            for i, (path, metrics) in enumerate(self.work_list[:20], 1):
                score = metrics.get('score', 0.0)
                fwhm = metrics.get('fwhm', 0.0)
                snr = metrics.get('snr', 0.0)
                contrast = metrics.get('contrast', 0.0)
                
                result_text += f"\n{i:3d}. {path.name}\n"
                result_text += f"      Score: {score:.2f} | FWHM: {fwhm:.2f}\" | SNR: {snr:.2f}\n"
        else:
            # Utiliser les données du fichier
            result_text = f"Liste de travail chargée\n"
            result_text += f"{'='*60}\n\n"
            result_text += f"Nombre d'images: {len(self.work_list)}\n"
            if 'source_directory' in data:
                result_text += f"Dossier source: {data['source_directory']}\n"
            if 'percent_used' in data:
                result_text += f"Pourcentage utilisé: {data['percent_used']:.0f}%\n"
            result_text += f"Fichier: {Path(self.work_list_file).name}\n"
            result_text += f"\n{'='*60}\n\n"
            result_text += f"Top {min(20, len(self.work_list))} images:\n"
            result_text += f"{'-'*60}\n"
            
            for i, item in enumerate(data['images'][:20], 1):
                filename = item.get('filename', item['path'])
                metrics = item.get('metrics', {})
                score = metrics.get('score', 0.0)
                fwhm = metrics.get('fwhm', 0.0)
                snr = metrics.get('snr', 0.0)
                
                result_text += f"\n{i:3d}. {filename}\n"
                result_text += f"      Score: {score:.2f} | FWHM: {fwhm:.2f}\" | SNR: {snr:.2f}\n"
        
        result_text += f"\n{'='*60}\n"
        result_text += f"Vous pouvez maintenant créer l'image empilée en utilisant cette liste.\n"
        
        self.lucky_result_text.delete(1.0, tk.END)
        self.lucky_result_text.insert(1.0, result_text)
    
    def view_work_list(self):
        """Affiche la liste de travail dans la zone de texte"""
        if not self.work_list:
            messagebox.showwarning("Attention", "Aucune liste de travail active!")
            return
        
        self.display_work_list_info()
    
    def run_stacking_from_work_list(self):
        """Exécute le stacking en utilisant la liste de travail"""
        if not REDUCTION_AVAILABLE or self.reducer is None:
            logger.error("[ANALYSE] Module de réduction non disponible")
            return
        
        if not self.work_list:
            messagebox.showwarning(
                "Attention",
                "Aucune liste de travail active!\n\n"
                "Veuillez d'abord analyser les images ou charger une liste sauvegardée."
            )
            return
        
        ref_x_str = self.ref_x_var.get()
        ref_y_str = self.ref_y_var.get()
        if not ref_x_str or not ref_y_str:
            messagebox.showwarning("Attention", "Entrez la position de référence (x, y) pour l'alignement!")
            return
        
        try:
            ref_x = float(ref_x_str)
            ref_y = float(ref_y_str)
        except ValueError:
            logger.error(f"[STACKING] Les positions doivent être des nombres (x={ref_x_str}, y={ref_y_str})")
            return
        
        # Déterminer le chemin de sauvegarde
        if self.auto_save_var.get():
            # Sauvegarde automatique dans le dossier source
            if self.work_list:
                source_dir = self.work_list[0][0].parent
            elif self.lucky_dir_var.get():
                source_dir = Path(self.lucky_dir_var.get())
            else:
                source_dir = self.base_dir
            
            # Générer un nom de fichier avec timestamp
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            method_name = self.lucky_method_var.get()
            n_images = len(self.work_list) if self.work_list else 0
            output_filename = f"ELI_{method_name}_{n_images}imgs_{timestamp}.fits"
            output_path = source_dir / output_filename
        else:
            # Demander à l'utilisateur de choisir l'emplacement
            default_dir = self.base_dir
            if self.work_list:
                default_dir = self.work_list[0][0].parent
            elif self.lucky_dir_var.get():
                default_dir = Path(self.lucky_dir_var.get())
            
            # Suggérer un nom de fichier
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            method_name = self.lucky_method_var.get()
            n_images = len(self.work_list) if self.work_list else 0
            default_filename = f"ELI_{method_name}_{n_images}imgs_{timestamp}.fits"
            
            output_path = filedialog.asksaveasfilename(
                initialdir=default_dir,
                initialfile=default_filename,
                defaultextension=".fits",
                filetypes=[("FITS", "*.fits"), ("Tous les fichiers", "*.*")]
            )
            
            if not output_path:
                return
            
            output_path = Path(output_path)
        
        def stacking_task():
            try:
                self.progress.start()
                
                method = self.lucky_method_var.get()
                best_image_files = [path for path, _ in self.work_list]
                
                # Afficher le statut
                self.lucky_result_text.delete(1.0, tk.END)
                result_text = f"Stacking en cours...\n"
                result_text += f"{'='*60}\n\n"
                result_text += f"Utilisation de la liste de travail: {len(best_image_files)} images\n"
                if self.work_list_file:
                    result_text += f"Source: {Path(self.work_list_file).name}\n"
                result_text += f"Méthode: {method}\n"
                result_text += f"Position référence: ({ref_x}, {ref_y})\n"
                self.lucky_result_text.insert(1.0, result_text)
                self.lucky_result_text.update()
                
                # Stacking (ELI)
                success = self.reducer.eli_lucky_imaging(
                    best_image_files,
                    output_path,
                    (ref_x, ref_y),
                    top_percent=1.0,  # Toutes les images sont déjà sélectionnées
                    method=method
                )
                
                self.progress.stop()
                
                if success:
                    result_text += f"\n{'='*60}\n"
                    result_text += f"✅ Image empilée créée avec succès!\n"
                    result_text += f"Fichier: {output_path.name}\n"
                    result_text += f"Chemin complet: {output_path}\n"
                    
                    self.lucky_result_text.delete(1.0, tk.END)
                    self.lucky_result_text.insert(1.0, result_text)
                    
                    messagebox.showinfo("Succès", f"Image empilée créée avec succès!\n\n{output_path}")
                    logger.info(f"Stacking terminé: {output_path}")
                else:
                    error_msg = "Échec de la création de l'image empilée"
                    logger.error(f"[STACKING] {error_msg}")
                    result_text += f"\n{'='*60}\n"
                    result_text += f"❌ {error_msg}\n"
                    self.lucky_result_text.delete(1.0, tk.END)
                    self.lucky_result_text.insert(1.0, result_text)
                    
            except Exception as e:
                self.progress.stop()
                logger.error(f"[STACKING] Erreur lors du stacking: {e}", exc_info=True)
                result_text = f"Erreur lors du stacking:\n{'='*60}\n{str(e)}\n"
                self.lucky_result_text.delete(1.0, tk.END)
                self.lucky_result_text.insert(1.0, result_text)
        
        threading.Thread(target=stacking_task, daemon=True).start()
    
    def refine_centroid(self, data, x, y, box_size=25, min_separation=3, exclude_positions=None):
        """
        Affine le centroïde d'une étoile autour d'une position approximative
        Peut discriminer deux étoiles très proches en détectant les pics locaux
        et en sélectionnant celui le plus proche du point de clic
        
        Parameters
        ----------
        data : np.ndarray
            Données de l'image
        x, y : float
            Position approximative (coordonnées pixel)
        box_size : int
            Taille de la région de recherche
        min_separation : float
            Séparation minimale en pixels pour considérer deux pics comme distincts
        exclude_positions : list of tuple, optional
            Liste de positions (x, y) à exclure de la détection (ex: position de l'étoile 1)
        
        Returns
        -------
        Tuple[float, float]
            Position du centroïde affiné (x, y)
        """
        try:
            h, w = data.shape
            x_int = int(round(x))
            y_int = int(round(y))
            
            # Vérifier que la position est dans l'image
            if x_int < box_size or x_int >= w - box_size or y_int < box_size or y_int >= h - box_size:
                return float('nan'), float('nan')
            
            # Extraire une région autour de la position
            x_min = max(0, x_int - box_size // 2)
            x_max = min(w, x_int + box_size // 2 + 1)
            y_min = max(0, y_int - box_size // 2)
            y_max = min(h, y_int + box_size // 2 + 1)
            
            region = data[y_min:y_max, x_min:x_max].astype(float)
            
            # Soustraire le fond (utiliser les bords de la région)
            border = np.concatenate([
                region[0, :], region[-1, :],
                region[:, 0], region[:, -1]
            ])
            background = np.median(border)
            region_sub = np.maximum(region - background, 0)
            
            # Lisser légèrement pour réduire le bruit et faciliter la détection des pics
            region_smooth = gaussian_filter(region_sub, sigma=0.5)
            
            # Détecter les pics locaux dans la région
            # Utiliser un filtre maximum pour trouver les maximums locaux
            footprint_size = max(3, int(min_separation))
            local_maxima = maximum_filter(region_smooth, size=footprint_size)
            max_value = np.max(region_smooth)
            if max_value > 0:
                peaks_mask = (region_smooth == local_maxima) & (region_smooth > max_value * 0.1)
            else:
                peaks_mask = np.zeros_like(region_smooth, dtype=bool)
            
            # Note: on n'exclut pas encore ici - on va sélectionner le pic le plus proche du clic
            # qui n'est pas dans les positions exclues (plus bas dans la logique)
            
            # Trouver les positions de tous les pics détectés
            peak_positions = np.where(peaks_mask)
            n_peaks = len(peak_positions[0])
            best_peak_idx = None
            
            if n_peaks == 0:
                # Aucun pic détecté, utiliser le maximum global
                max_pos = np.unravel_index(np.argmax(region_sub), region_sub.shape)
                best_peak_local = (max_pos[1], max_pos[0])  # (x_local, y_local)
            else:
                # Extraire les valeurs des pics détectés
                peak_values = region_smooth[peak_positions]
                
                # Trier les pics par intensité (plus brillant d'abord)
                sorted_indices = np.argsort(peak_values)[::-1]
                
                # Trouver le pic le plus proche de la position de clic dans la région locale
                # en excluant les positions déjà sélectionnées
                click_x_local = x - x_min
                click_y_local = y - y_min
                min_distance = float('inf')
                best_peak_idx = sorted_indices[0]  # Initialiser avec le plus brillant
                
                for idx in sorted_indices:
                    peak_y_local = peak_positions[0][idx]
                    peak_x_local = peak_positions[1][idx]
                    
                    # Vérifier si ce pic est trop proche d'une position exclue
                    skip_peak = False
                    if exclude_positions:
                        for excl_x, excl_y in exclude_positions:
                            excl_x_local = excl_x - x_min
                            excl_y_local = excl_y - y_min
                            distance_to_excl = np.sqrt((peak_x_local - excl_x_local)**2 + (peak_y_local - excl_y_local)**2)
                            # Si le pic est très proche d'une position exclue (< min_separation), le sauter
                            if distance_to_excl < min_separation:
                                skip_peak = True
                                break
                    
                    if skip_peak:
                        continue
                    
                    distance = np.sqrt((peak_x_local - click_x_local)**2 + (peak_y_local - click_y_local)**2)
                    
                    # Préférer les pics proches même s'ils sont moins brillants
                    # Mais si un pic est beaucoup plus brillant (>2x), le considérer aussi
                    if distance < min_distance:
                        min_distance = distance
                        best_peak_idx = idx
                    elif peak_values[idx] > peak_values[best_peak_idx] * 2.0 and distance < min_distance * 1.5:
                        # Si beaucoup plus brillant et pas trop loin, le considérer
                        min_distance = distance
                        best_peak_idx = idx
                
                # Si aucun pic valide n'a été trouvé (tous exclus), utiliser le plus proche du clic
                if min_distance == float('inf') and n_peaks > 0:
                    # Forcer la sélection du pic le plus proche du clic, même s'il est proche d'une position exclue
                    for idx in sorted_indices:
                        peak_y_local = peak_positions[0][idx]
                        peak_x_local = peak_positions[1][idx]
                        distance = np.sqrt((peak_x_local - click_x_local)**2 + (peak_y_local - click_y_local)**2)
                        if distance < min_distance:
                            min_distance = distance
                            best_peak_idx = idx
                
                best_peak_local = (peak_positions[1][best_peak_idx], peak_positions[0][best_peak_idx])
            
            # Coordonnées locales du pic sélectionné
            max_x_local, max_y_local = best_peak_local
            
            # Utiliser center_of_mass pour un centroïde plus précis
            # Utiliser une région plus petite centrée sur le pic sélectionné
            sub_size = max(7, int(min_separation * 2))
            sub_x_min = max(0, max_x_local - sub_size // 2)
            sub_x_max = min(region_sub.shape[1], max_x_local + sub_size // 2 + 1)
            sub_y_min = max(0, max_y_local - sub_size // 2)
            sub_y_max = min(region_sub.shape[0], max_y_local + sub_size // 2 + 1)
            
            sub_region = region_sub[sub_y_min:sub_y_max, sub_x_min:sub_x_max]
            
            # Masquer les autres pics proches pour éviter leur interférence
            # Créer un masque qui exclut les autres pics à plus de min_separation pixels
            mask = np.ones_like(sub_region, dtype=bool)
            if n_peaks > 0 and best_peak_idx is not None:
                for idx in range(n_peaks):
                    if idx != best_peak_idx:
                        other_peak_x = peak_positions[1][idx] - sub_x_min
                        other_peak_y = peak_positions[0][idx] - sub_y_min
                        # Si l'autre pic est dans la sous-région, le masquer
                        if 0 <= other_peak_x < sub_region.shape[1] and 0 <= other_peak_y < sub_region.shape[0]:
                            # Créer un masque circulaire autour de l'autre pic
                            y_coords, x_coords = np.ogrid[:sub_region.shape[0], :sub_region.shape[1]]
                            distance_from_other = np.sqrt((x_coords - other_peak_x)**2 + (y_coords - other_peak_y)**2)
                            mask &= (distance_from_other > min_separation / 2)
            
            sub_region_masked = sub_region.copy()
            sub_region_masked[~mask] = 0
            
            # S'assurer que les valeurs sont positives pour center_of_mass
            sub_region_masked = np.maximum(sub_region_masked, 0)
            
            if np.sum(sub_region_masked) > 0:
                cy_local, cx_local = center_of_mass(sub_region_masked)
                # Vérifier que le centroïde est valide
                if np.isfinite(cx_local) and np.isfinite(cy_local):
                    # Ajuster pour les coordonnées globales
                    ref_x = x_min + sub_x_min + cx_local
                    ref_y = y_min + sub_y_min + cy_local
                    return float(ref_x), float(ref_y)
            
            # Fallback : utiliser la position du pic si center_of_mass échoue
            ref_x = x_min + max_x_local
            ref_y = y_min + max_y_local
            return float(ref_x), float(ref_y)
            
        except Exception as e:
            logger.error(f"[CENTROID] Erreur lors du raffinage du centroïde: {e}", exc_info=True)
            return float('nan'), float('nan')
    
    def show_reference_image(self):
        """Affiche l'image de référence (meilleure image) pour sélection de cible"""
        if not self.work_list:
            messagebox.showwarning(
                "Attention",
                "Aucune liste de travail disponible!\n\n"
                "Veuillez d'abord analyser les images ou charger une liste de travail."
            )
            return
        
        # Prendre la première image (meilleure qualité)
        ref_path, ref_metrics = self.work_list[0]
        
        if not ref_path.exists():
            logger.error(f"[IMAGE_REF] Image non trouvée: {ref_path}")
            return
        
        try:
            # Charger l'image (copie explicite pour éviter données invalides après fermeture FITS)
            with fits.open(ref_path) as hdul:
                raw = hdul[0].data
                data = np.array(raw, dtype=float, copy=True)
                data = np.squeeze(data)  # Enlever dimensions (1, ny, nx) -> (ny, nx)
            
            if data.size == 0 or data.ndim != 2:
                messagebox.showerror("Erreur", f"Image invalide: dimensions {data.shape}")
                return
            
            # Créer une fenêtre de visualisation
            ref_window = tk.Toplevel(self)
            ref_window.title(f"Image de référence - {ref_path.name}")
            ref_window.geometry("1200x1050")
            
            # Frame principal
            main_frame = ttk.Frame(ref_window)
            main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=3)
            
            # Zone d'instructions
            info_frame = ttk.Frame(main_frame)
            info_frame.pack(fill="x", pady=(0, 3))
            
            ttk.Label(
                info_frame,
                text="Clic gauche sur une étoile brillante pour sélectionner la cible de référence.",
                font=("Helvetica", 10, "bold"),
                foreground="blue"
            ).pack()
            
            ttk.Label(
                info_frame,
                text=f"Image: {ref_path.name} | Score: {ref_metrics.get('score', 0):.2f} | FWHM: {ref_metrics.get('fwhm', 0):.2f}\"",
                font=("Helvetica", 9),
                foreground="gray"
            ).pack(pady=(2, 0))
            
            # Zone de visualisation matplotlib
            fig_frame = ttk.Frame(main_frame)
            fig_frame.pack(fill=tk.BOTH, expand=True, pady=(2, 0))
            
            fig = Figure(figsize=(10, 8), dpi=100)
            fig.subplots_adjust(left=0.08, bottom=0.08, right=0.92, top=0.92, wspace=0.2, hspace=0.2)
            ax = fig.add_subplot(111)
            
            # Nettoyer NaN/Inf pour ZScale
            data_finite = np.ma.masked_invalid(data)
            if data_finite.count() == 0:
                messagebox.showerror("Erreur", "L'image ne contient que des NaN/Inf.")
                return
            
            # Calculer les limites initiales avec ZScale
            try:
                interval = ZScaleInterval()
                vmin_init, vmax_init = interval.get_limits(data)
            except Exception:
                vmin_init, vmax_init = np.nanpercentile(data, [10, 90])
            
            # Stocker les valeurs min/max absolues pour les sliders
            data_min = float(np.nanmin(data))
            data_max = float(np.nanmax(data))
            if data_max <= data_min:
                data_max = data_min + 1.0
            
            # Ajuster les limites initiales pour être dans la plage valide
            vmin_init = max(data_min, min(vmin_init, vmax_init - 1e-6))
            vmax_init = min(data_max, max(vmax_init, vmin_init + 1e-6))
            if vmax_init <= vmin_init:
                vmin_init, vmax_init = data_min, data_max
            
            # Variables pour le contraste
            vmin = vmin_init
            vmax = vmax_init
            
            im = ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax, aspect='equal', interpolation='nearest')
            ax.set_xlabel("X (pixels)")
            ax.set_ylabel("Y (pixels)")
            ax.set_title("Image de référence - Cliquez sur une étoile")
            cbar = fig.colorbar(im, ax=ax, label='Intensité')
            
            canvas = FigureCanvasTkAgg(fig, master=fig_frame)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            canvas.draw()
            ref_window.update_idletasks()
            # Redessiner après affichage (évite écran blanc si layout pas encore finalisé)
            def _delayed_draw():
                try:
                    canvas.draw_idle()
                except Exception:
                    pass
            ref_window.after(200, _delayed_draw)
            
            # Barre d'outils matplotlib avec zoom et navigation - directement sous l'image
            toolbar = NavigationToolbar2Tk(canvas, fig_frame)
            toolbar.update()
            toolbar.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 0))
            
            # Frame pour les contrôles de contraste
            contrast_frame = ttk.LabelFrame(main_frame, text="Contrôle du contraste", padding=3)
            contrast_frame.pack(fill="x", pady=(3, 0))
            
            # Fonction pour mettre à jour l'affichage
            def update_image_display():
                im.set_clim(vmin=vmin, vmax=vmax)
                cbar.update_normal(im)
                canvas.draw_idle()
            
            # Slider pour vmin (minimum)
            min_frame = ttk.Frame(contrast_frame)
            min_frame.pack(fill="x", pady=(0, 2))
            ttk.Label(min_frame, text="Minimum:", width=8).pack(side=tk.LEFT)
            min_var = tk.DoubleVar(value=vmin)
            min_slider = ttk.Scale(min_frame, from_=data_min, to=data_max, 
                                  variable=min_var, orient=tk.HORIZONTAL, length=300)
            min_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(3, 3))
            min_label = ttk.Label(min_frame, text=f"{vmin:.2f}", width=8)
            min_label.pack(side=tk.LEFT)
            
            def update_min(val):
                nonlocal vmin
                vmin = float(min_var.get())
                # S'assurer que vmin < vmax
                if vmin >= vmax:
                    vmin = vmax - (data_max - data_min) * 0.01
                    min_var.set(vmin)
                min_label.config(text=f"{vmin:.2f}")
                update_image_display()
            
            min_slider.config(command=update_min)
            
            # Slider pour vmax (maximum)
            max_frame = ttk.Frame(contrast_frame)
            max_frame.pack(fill="x")
            ttk.Label(max_frame, text="Maximum:", width=8).pack(side=tk.LEFT)
            max_var = tk.DoubleVar(value=vmax)
            max_slider = ttk.Scale(max_frame, from_=data_min, to=data_max, 
                                  variable=max_var, orient=tk.HORIZONTAL, length=300)
            max_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(3, 3))
            max_label = ttk.Label(max_frame, text=f"{vmax:.2f}", width=8)
            max_label.pack(side=tk.LEFT)
            
            def update_max(val):
                nonlocal vmax
                vmax = float(max_var.get())
                # S'assurer que vmax > vmin
                if vmax <= vmin:
                    vmax = vmin + (data_max - data_min) * 0.01
                    max_var.set(vmax)
                max_label.config(text=f"{vmax:.2f}")
                update_image_display()
            
            max_slider.config(command=update_max)
            
            # Bouton reset
            reset_contrast_btn = ttk.Button(contrast_frame, text="Réinitialiser (ZScale)", 
                                           command=lambda: reset_contrast())
            reset_contrast_btn.pack(pady=(3, 0))
            
            def reset_contrast():
                nonlocal vmin, vmax
                vmin = vmin_init
                vmax = vmax_init
                min_var.set(vmin)
                max_var.set(vmax)
                min_label.config(text=f"{vmin:.2f}")
                max_label.config(text=f"{vmax:.2f}")
                update_image_display()
            
            # Variables pour stocker la position sélectionnée
            selected_pos = {'x': None, 'y': None}
            marker_point = None
            
            def on_click(event):
                """Gère le clic sur l'image"""
                # Ignorer si le clic n'est pas dans les axes ou si la toolbar est active
                if event.inaxes != ax or event.button != 1:
                    return
                
                # Ignorer si on est en mode zoom/pan de la toolbar
                # Vérifier le mode actif de la toolbar matplotlib
                toolbar_mode = getattr(toolbar, '_active', None) or getattr(toolbar, 'mode', None)
                if toolbar_mode and toolbar_mode != '' and toolbar_mode != 'NONE':
                    return
                
                # Position du clic
                click_x = event.xdata
                click_y = event.ydata
                
                if click_x is None or click_y is None:
                    return
                
                # Rafiner le centroïde
                refined_x, refined_y = self.refine_centroid(data, click_x, click_y)
                
                if np.isfinite(refined_x) and np.isfinite(refined_y):
                    selected_pos['x'] = refined_x
                    selected_pos['y'] = refined_y
                    
                    # Afficher un marqueur
                    nonlocal marker_point
                    if marker_point:
                        marker_point.remove()
                    
                    marker_point = ax.plot(refined_x, refined_y, 'r+', markersize=20, 
                                          markeredgewidth=2, label='Cible sélectionnée')[0]
                    
                    # Afficher les coordonnées
                    if hasattr(ax, 'text_target'):
                        ax.text_target.remove()
                    
                    ax.text_target = ax.text(0.02, 0.98, 
                                            f"Position: x={refined_x:.2f}, y={refined_y:.2f}",
                                            transform=ax.transAxes,
                                            fontsize=12,
                                            verticalalignment='top',
                                            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
                    
                    canvas.draw()
                    
                    # Mettre à jour les champs dans l'onglet principal
                    self.ref_x_var.set(f"{refined_x:.2f}")
                    self.ref_y_var.set(f"{refined_y:.2f}")
                    
                    logger.info(f"Cible sélectionnée: ({refined_x:.2f}, {refined_y:.2f})")
                else:
                    messagebox.showwarning("Attention", "Impossible de calculer le centroïde à cette position.")
            
            # Connecter le gestionnaire de clic
            canvas.mpl_connect('button_press_event', on_click)
            
            # Bouton de fermeture
            button_frame = ttk.Frame(main_frame)
            button_frame.pack(fill="x", pady=(3, 0))
            
            ttk.Button(
                button_frame,
                text="Fermer",
                command=ref_window.destroy
            ).pack()
            
            logger.info(f"Affichage de l'image de référence: {ref_path.name}")
            
        except Exception as e:
            logger.error(f"[IMAGE_REF] Erreur lors du chargement de l'image: {e}", exc_info=True)
    
    def show_measure_image(self):
        """Affiche l'image pour la mesure de séparation avec sélection interactive"""
        image_path = self.measure_image_var.get()
        if not image_path:
            messagebox.showwarning("Attention", "Sélectionnez d'abord une image!")
            return
        
        image_path = Path(image_path)
        if not image_path.exists():
            logger.error(f"[MESURE] Image non trouvée: {image_path}")
            return
        
        try:
            # Charger l'image
            with fits.open(image_path) as hdul:
                data = hdul[0].data.astype(float)
                header = hdul[0].header
            
            # Vérifier que WCS est présent dans le header (plate-solving requis)
            try:
                wcs = WCS(header)
                if not wcs.has_celestial:
                    raise ValueError("WCS n'est pas céleste")
            except Exception as e:
                messagebox.showerror(
                    "Erreur",
                    f"WCS (plate-solving) non trouvé dans l'image!\n\n"
                    f"L'image doit avoir été astrométriée (plate-solving)\n"
                    f"pour déterminer le Nord céleste.\n\n"
                    f"Erreur: {e}"
                )
                logger.error(f"[MESURE] WCS non trouvé: {e}")
                return
            
            # Calculer la direction du Nord céleste via WCS (plate-solving)
            # On prend un point au centre de l'image et un point légèrement au-dessus
            ny, nx = data.shape
            center_x, center_y = nx / 2.0, ny / 2.0
            offset = min(ny, nx) * 0.1  # 10% de l'image comme décalage (environ 50-100 pixels)
            
            # Calculer l'orientation réelle via WCS (méthode robuste)
            orientation = get_image_orientation(wcs, center_x, center_y)
            north_celestial_angle = orientation['north_angle']
            east_celestial_angle = orientation['east_angle']
            
            logger.info(f"[MESURE] Orientation WCS: N={north_celestial_angle:.2f}°, E={east_celestial_angle:.2f}°")
            logger.info(f"[MESURE] Rotation capteur: {orientation['rotation']:.2f}°, Parité: {orientation['parity']}")
            
            # Calculer la taille du pixel depuis le WCS
            try:
                pixel_scales = proj_plane_pixel_scales(wcs)  # En degrés/pixel
                pixel_scale_arcsec = np.mean(pixel_scales) * 3600.0  # Conversion en arcsec/pixel
                logger.info(f"[MESURE] Échelle pixel calculée depuis WCS: {pixel_scale_arcsec:.3f}\" / pixel")
            except Exception as e:
                pixel_scale_arcsec = None
                logger.warning(f"[MESURE] Impossible de calculer l'échelle pixel depuis WCS: {e}")
            
            # Créer une fenêtre de visualisation
            measure_window = tk.Toplevel(self)
            measure_window.title(f"Mesure de séparation - {image_path.name}")
            measure_window.geometry("1100x1050")
            
            # Frame principal
            main_frame = ttk.Frame(measure_window)
            main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=3)
            
            # Zone d'instructions
            info_frame = ttk.Frame(main_frame)
            info_frame.pack(fill="x", pady=(0, 3))
            
            self.measure_instruction_label = ttk.Label(
                info_frame,
                text="1️⃣ Cliquez sur l'étoile 1 (primaire)",
                font=("Helvetica", 10, "bold"),
                foreground="blue"
            )
            self.measure_instruction_label.pack()
            
            ttk.Label(
                info_frame,
                text=f"Image: {image_path.name}",
                font=("Helvetica", 9),
                foreground="gray"
            ).pack(pady=(2, 0))
            
            # Zone de visualisation matplotlib
            fig_frame = ttk.Frame(main_frame)
            fig_frame.pack(fill=tk.BOTH, expand=True, pady=(2, 0))
            
            fig = Figure(figsize=(10, 8), dpi=100)
            # Réduire fortement les marges blanches
            fig.subplots_adjust(left=0.08, bottom=0.05, right=0.90, top=0.90, wspace=0.1, hspace=0.1)
            ax = fig.add_subplot(111)
            
            # Calculer la médiane pour centrer les sliders
            data_median = float(np.nanmedian(data))
            data_min_abs = float(np.nanmin(data))
            data_max_abs = float(np.nanmax(data))
            
            # Calculer les limites initiales avec ZScale
            try:
                interval = ZScaleInterval()
                vmin_init, vmax_init = interval.get_limits(data)
            except:
                vmin_init, vmax_init = np.percentile(data, [10, 90])
            
            # Échelle centrée sur médiane avec 5000 unités
            scale_range = 2500.0  # ±2500 unités autour de la médiane = 5000 unités total
            slider_min = max(data_min_abs, data_median - scale_range)
            slider_max = min(data_max_abs, data_median + scale_range)
            
            # Limites initiales centrées sur médiane
            vmin_init = max(slider_min, data_median - 100)
            vmax_init = min(slider_max, data_median + 100)
            
            # Variables pour le contraste
            vmin = vmin_init
            vmax = vmax_init
            
            im = ax.imshow(data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax, aspect='equal')
            ax.set_xlabel("X (pixels)")
            ax.set_ylabel("Y (pixels)")
            ax.set_title("Mesure de séparation - Sélectionnez les deux étoiles")
            cbar = fig.colorbar(im, ax=ax, label='Intensité')
            
            canvas = FigureCanvasTkAgg(fig, master=fig_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            
            # Barre d'outils matplotlib avec zoom et navigation - directement sous l'image
            toolbar = NavigationToolbar2Tk(canvas, fig_frame)
            toolbar.update()
            toolbar.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 0))
            
            # Frame pour les contrôles de contraste
            contrast_frame = ttk.LabelFrame(main_frame, text="Contrôle du contraste", padding=3)
            contrast_frame.pack(fill="x", pady=(3, 0))
            
            # Fonction pour mettre à jour l'affichage
            def update_image_display():
                im.set_clim(vmin=vmin, vmax=vmax)
                cbar.update_normal(im)
                canvas.draw_idle()
            
            # Slider pour vmin (minimum) - centré sur médiane avec échelle de 500 unités
            min_frame = ttk.Frame(contrast_frame)
            min_frame.pack(fill="x", pady=(0, 2))
            ttk.Label(min_frame, text="Minimum:", width=8).pack(side=tk.LEFT)
            min_var = tk.DoubleVar(value=vmin)
            min_slider = ttk.Scale(min_frame, from_=slider_min, to=slider_max, 
                                  variable=min_var, orient=tk.HORIZONTAL, length=300)
            min_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(3, 3))
            min_label = ttk.Label(min_frame, text=f"{vmin:.2f}", width=8)
            min_label.pack(side=tk.LEFT)
            
            def update_min(val):
                nonlocal vmin
                vmin = float(min_var.get())
                # S'assurer que vmin < vmax
                if vmin >= vmax:
                    vmin = vmax - (slider_max - slider_min) * 0.01
                    min_var.set(vmin)
                min_label.config(text=f"{vmin:.2f}")
                update_image_display()
            
            min_slider.config(command=update_min)
            
            # Slider pour vmax (maximum) - centré sur médiane avec échelle de 500 unités
            max_frame = ttk.Frame(contrast_frame)
            max_frame.pack(fill="x")
            ttk.Label(max_frame, text="Maximum:", width=8).pack(side=tk.LEFT)
            max_var = tk.DoubleVar(value=vmax)
            max_slider = ttk.Scale(max_frame, from_=slider_min, to=slider_max, 
                                  variable=max_var, orient=tk.HORIZONTAL, length=300)
            max_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(3, 3))
            max_label = ttk.Label(max_frame, text=f"{vmax:.2f}", width=8)
            max_label.pack(side=tk.LEFT)
            
            def update_max(val):
                nonlocal vmax
                vmax = float(max_var.get())
                # S'assurer que vmax > vmin
                if vmax <= vmin:
                    vmax = vmin + (slider_max - slider_min) * 0.01
                    max_var.set(vmax)
                max_label.config(text=f"{vmax:.2f}")
                update_image_display()
            
            max_slider.config(command=update_max)
            
            # Bouton reset
            reset_contrast_btn = ttk.Button(contrast_frame, text="Réinitialiser (Médiane)", 
                                           command=lambda: reset_contrast())
            reset_contrast_btn.pack(pady=(3, 0))
            
            def reset_contrast():
                nonlocal vmin, vmax
                vmin = vmin_init
                vmax = vmax_init
                min_var.set(vmin)
                max_var.set(vmax)
                min_label.config(text=f"{vmin:.2f}")
                max_label.config(text=f"{vmax:.2f}")
                update_image_display()
            
            # Variables pour stocker les positions
            star1_pos = {'x': None, 'y': None, 'refined_x': None, 'refined_y': None}
            star2_pos = {'x': None, 'y': None, 'refined_x': None, 'refined_y': None}
            marker_star1 = None
            marker_star2 = None
            selection_step = 1  # 1 = étoile 1, 2 = étoile 2
            # Variables pour l'affichage de l'angle theta et des éléments graphiques
            theta_arrow_north = None
            theta_arrow_east = None
            theta_arrow_separation = None
            theta_arc = None
            theta_text = None
            compass_north = None
            compass_east = None
            separation_line = None
            rho_text = None
            
            def update_instruction():
                if selection_step == 1:
                    self.measure_instruction_label.config(
                        text="1️⃣ Cliquez sur l'étoile 1 (primaire)",
                        foreground="blue"
                    )
                elif selection_step == 2:
                    self.measure_instruction_label.config(
                        text="2️⃣ Cliquez sur l'étoile 2 (secondaire)",
                        foreground="orange"
                    )
                else:
                    self.measure_instruction_label.config(
                        text="✅ Les deux étoiles sont sélectionnées",
                        foreground="green"
                    )
            
            def on_click(event):
                """Gère le clic sur l'image"""
                # Ignorer si le clic n'est pas dans les axes ou si la toolbar est active
                if event.inaxes != ax or event.button != 1:
                    return
                
                # Ignorer si on est en mode zoom/pan de la toolbar
                # Vérifier le mode actif de la toolbar matplotlib
                toolbar_mode = getattr(toolbar, '_active', None) or getattr(toolbar, 'mode', None)
                if toolbar_mode and toolbar_mode != '' and toolbar_mode != 'NONE':
                    return
                
                # Position du clic
                click_x = event.xdata
                click_y = event.ydata
                
                if click_x is None or click_y is None:
                    return
                
                nonlocal marker_star1, marker_star2, selection_step
                
                if selection_step == 1:
                    # Sélection étoile 1
                    refined_x, refined_y = self.refine_centroid(data, click_x, click_y)
                    
                    if not np.isfinite(refined_x) or not np.isfinite(refined_y):
                        messagebox.showwarning("Attention", "Impossible de calculer le centroïde à cette position.")
                        return
                    
                    star1_pos['x'] = click_x
                    star1_pos['y'] = click_y
                    star1_pos['refined_x'] = refined_x
                    star1_pos['refined_y'] = refined_y
                    
                    # Afficher un marqueur
                    if marker_star1:
                        marker_star1.remove()
                    
                    marker_star1 = ax.plot(refined_x, refined_y, 'w+', markersize=20, 
                                          markeredgewidth=2, label='Étoile 1')[0]
                    ax.text(refined_x + 5, refined_y + 5, 'Star 1', color='white', 
                           fontsize=10, weight='bold')
                    
                    # Mettre à jour les champs
                    self.star1_x_var.set(f"{refined_x:.2f}")
                    self.star1_y_var.set(f"{refined_y:.2f}")
                    
                    selection_step = 2
                    update_instruction()
                    logger.info(f"Étoile 1 sélectionnée: ({refined_x:.2f}, {refined_y:.2f})")
                    
                elif selection_step == 2:
                    # Sélection étoile 2 - exclure la position de l'étoile 1
                    # Rafiner le centroïde en excluant la position de l'étoile 1
                    if star1_pos['refined_x'] is not None and star1_pos['refined_y'] is not None:
                        refined_x, refined_y = self.refine_centroid(
                            data, click_x, click_y, 
                            exclude_positions=[(star1_pos['refined_x'], star1_pos['refined_y'])]
                        )
                    else:
                        refined_x, refined_y = self.refine_centroid(data, click_x, click_y)
                    
                    if not np.isfinite(refined_x) or not np.isfinite(refined_y):
                        messagebox.showwarning("Attention", "Impossible de calculer le centroïde à cette position. Essayez de cliquer plus loin de l'étoile 1.")
                        return
                    
                    # Sélection étoile 2
                    star2_pos['x'] = click_x
                    star2_pos['y'] = click_y
                    star2_pos['refined_x'] = refined_x
                    star2_pos['refined_y'] = refined_y
                    
                    # Afficher un marqueur
                    if marker_star2:
                        marker_star2.remove()
                    
                    marker_star2 = ax.plot(refined_x, refined_y, 'w+', markersize=20, 
                                          markeredgewidth=2, label='Étoile 2')[0]
                    ax.text(refined_x + 5, refined_y + 5, 'Star 2', color='white', 
                           fontsize=10, weight='bold')
                    
                    # Calculer la séparation et l'angle theta
                    dx = refined_x - star1_pos['refined_x']
                    dy = refined_y - star1_pos['refined_y']
                    separation = np.sqrt(dx**2 + dy**2)
                    
                    # Position de l'étoile 1 (point de référence T1)
                    x1 = star1_pos['refined_x']
                    y1 = star1_pos['refined_y']
                    
                    # Calculer l'angle de position theta depuis le Nord céleste vers l'Est (WCS)
                    coord_t1 = wcs.pixel_to_world(x1, y1)
                    coord_t2 = wcs.pixel_to_world(refined_x, refined_y)
                    dra_deg = (coord_t2.ra.deg - coord_t1.ra.deg) * np.cos(np.radians(coord_t1.dec.deg))
                    ddec_deg = coord_t2.dec.deg - coord_t1.dec.deg
                    theta = np.degrees(np.arctan2(dra_deg, ddec_deg))
                    if theta < 0:
                        theta += 360.0
                    
                    # Calculer la séparation en secondes d'arc (rho)
                    # Utiliser le pixel scale calculé depuis le WCS, ou celui fourni par l'utilisateur
                    pixel_scale_to_use = None
                    if pixel_scale_arcsec is not None and pixel_scale_arcsec > 0:
                        pixel_scale_to_use = pixel_scale_arcsec
                        logger.debug(f"[MESURE] Utilisation de l'échelle pixel WCS: {pixel_scale_to_use:.3f}\" / pixel")
                    else:
                        # Fallback sur la valeur fournie par l'utilisateur
                        try:
                            pixel_scale_str = self.pixel_scale_var.get().strip()
                            if pixel_scale_str:
                                pixel_scale_to_use = float(pixel_scale_str)
                                if pixel_scale_to_use <= 0:
                                    pixel_scale_to_use = None
                        except (ValueError, TypeError):
                            pixel_scale_to_use = None
                    
                    if pixel_scale_to_use is not None and pixel_scale_to_use > 0:
                        separation_arcsec = separation * pixel_scale_to_use
                        rho_value = separation_arcsec
                    else:
                        rho_value = None
                    
                    # Dessiner l'affichage de l'angle theta
                    nonlocal theta_arrow_north, theta_arrow_east, theta_arrow_separation, theta_arc, theta_text
                    nonlocal compass_north, compass_east, separation_line, rho_text
                    
                    # Nettoyer les éléments précédents
                    for elem in [theta_arrow_north, theta_arrow_east, theta_arrow_separation, theta_arc, theta_text, 
                                compass_north, compass_east, separation_line, rho_text]:
                        if elem is not None:
                            try:
                                elem.remove()
                            except:
                                pass
                    
                    # Longueur des flèches pour le compas (allongées nettement)
                    compass_length = min(separation * 1.2, 300)
                    
                    # Dessiner le Nord céleste (flèche depuis T1) - en blanc
                    # north_celestial_angle est déjà dans le système matplotlib (0°=droite, anti-horaire)
                    north_angle_rad = np.radians(north_celestial_angle)
                    from matplotlib.patches import FancyArrowPatch
                    compass_north = FancyArrowPatch(
                        (x1, y1),
                        (x1 + compass_length * np.cos(north_angle_rad), 
                         y1 + compass_length * np.sin(north_angle_rad)),
                        arrowstyle='->', mutation_scale=20, color='white', lw=2, alpha=0.8
                    )
                    ax.add_patch(compass_north)
                    # Label Nord - en blanc, fond transparent
                    ax.text(x1 + compass_length * 0.7 * np.cos(north_angle_rad), 
                           y1 + compass_length * 0.7 * np.sin(north_angle_rad) + 10,
                           'N', color='white', fontsize=12, weight='bold',
                           ha='center', va='bottom')
                    
                    # Dessiner l'Est céleste (flèche depuis T1) - en blanc
                    # east_celestial_angle est déjà dans le système matplotlib (0°=droite, anti-horaire)
                    east_angle_rad = np.radians(east_celestial_angle)
                    compass_east = FancyArrowPatch(
                        (x1, y1),
                        (x1 + compass_length * np.cos(east_angle_rad), 
                         y1 + compass_length * np.sin(east_angle_rad)),
                        arrowstyle='->', mutation_scale=20, color='white', lw=2, alpha=0.8
                    )
                    ax.add_patch(compass_east)
                    # Label Est - en blanc, fond transparent
                    ax.text(x1 + compass_length * 0.7 * np.cos(east_angle_rad) + 10,
                           y1 + compass_length * 0.7 * np.sin(east_angle_rad),
                           'E', color='white', fontsize=12, weight='bold',
                           ha='left', va='center')
                    
                    # Dessiner la ligne de séparation T1-T2
                    separation_line = ax.plot([x1, refined_x], 
                                              [y1, refined_y], 
                                              'w-', linewidth=2, alpha=0.8, label='Séparation T1-T2')[0]
                    
                    # Afficher ρ (séparation) au milieu de la droite T1-T2
                    mid_x = (x1 + refined_x) / 2.0
                    mid_y = (y1 + refined_y) / 2.0
                    # Perpendiculaire à la ligne pour le placement du texte
                    # Angle de la ligne T1T2 dans le système pixel
                    t1t2_angle_pixel = np.degrees(np.arctan2(dy, dx))
                    if t1t2_angle_pixel < 0:
                        t1t2_angle_pixel += 360.0
                    perp_angle = t1t2_angle_pixel + 90.0
                    offset_rho = 15.0  # Décalage du texte
                    rho_text_x = mid_x + offset_rho * np.cos(np.radians(perp_angle))
                    rho_text_y = mid_y + offset_rho * np.sin(np.radians(perp_angle))
                    
                    if rho_value is not None:
                        rho_text = ax.text(rho_text_x, rho_text_y,
                                           f'ρ = {rho_value:.3f}"',
                                           color='white', fontsize=11, weight='bold',
                                           ha='center', va='center')
                    else:
                        rho_text = ax.text(rho_text_x, rho_text_y,
                                           f'ρ = {separation:.2f} px',
                                           color='white', fontsize=11, weight='bold',
                                           ha='center', va='center')
                    
                    # Dessiner l'arc pointillé pour theta (Nord -> Est -> séparation)
                    # Le sens est CONSTANT : toujours du Nord vers l'Est (même sens que Nord->Est)
                    # indépendamment de la valeur de theta
                    arc_radius = min(separation * 0.35, 120)
                    
                    # Angles dans le système matplotlib (en radians)
                    north_angle_rad = np.radians(north_celestial_angle)
                    east_angle_rad = np.radians(east_celestial_angle)
                    separation_angle_rad = np.radians(t1t2_angle_pixel)
                    
                    # Calculer le sens du Nord vers l'Est
                    diff_north_to_east = (east_angle_rad - north_angle_rad) % (2 * np.pi)
                    use_ccw = diff_north_to_east <= np.pi  # True si Est est dans le sens anti-horaire depuis le Nord
                    
                    # Tracer l'arc du Nord vers la séparation en utilisant le MÊME sens que Nord->Est
                    diff_north_to_sep_ccw = (separation_angle_rad - north_angle_rad) % (2 * np.pi)
                    diff_north_to_sep_cw = 2 * np.pi - diff_north_to_sep_ccw
                    
                    if use_ccw:
                        # Même sens que Nord->Est : sens anti-horaire
                        if diff_north_to_sep_ccw <= np.pi:
                            # Chemin direct anti-horaire
                            arc_angles_rad = np.linspace(north_angle_rad, separation_angle_rad, 50)
                        else:
                            # Passage par 0 en sens anti-horaire
                            arc1 = np.linspace(north_angle_rad, 2 * np.pi, 25)
                            arc2 = np.linspace(0.0, separation_angle_rad, 25)
                            arc_angles_rad = np.concatenate([arc1, arc2])
                    else:
                        # Même sens que Nord->Est : sens horaire
                        if diff_north_to_sep_cw <= np.pi:
                            # Chemin direct horaire
                            if separation_angle_rad < north_angle_rad:
                                arc_angles_rad = np.linspace(north_angle_rad, separation_angle_rad, 50)
                            else:
                                # Passage par 0 en sens horaire
                                arc1 = np.linspace(north_angle_rad, 0.0, 25)
                                arc2 = np.linspace(2 * np.pi, separation_angle_rad, 25)
                                arc_angles_rad = np.concatenate([arc1, arc2])
                        else:
                            # Chemin horaire avec passage par 0
                            arc1 = np.linspace(north_angle_rad, 0.0, 25)
                            arc2 = np.linspace(2 * np.pi, separation_angle_rad, 25)
                            arc_angles_rad = np.concatenate([arc1, arc2])
                    
                    # Tracer l'arc
                    arc_x = x1 + arc_radius * np.cos(arc_angles_rad)
                    arc_y = y1 + arc_radius * np.sin(arc_angles_rad)
                    theta_arc = ax.plot(arc_x, arc_y, color='white', lw=1, alpha=0.7, linestyle='--')[0]
                    
                    # Texte avec la valeur de theta au milieu de l'arc
                    # Calculer l'angle au milieu de l'arc tracé
                    if len(arc_angles_rad) > 0:
                        mid_idx = len(arc_angles_rad) // 2
                        mid_theta_rad = arc_angles_rad[mid_idx]
                    else:
                        mid_theta_rad = start_angle_rad + np.radians(theta / 2.0) if theta <= 180.0 else start_angle_rad - np.radians((360.0 - theta) / 2.0)
                        mid_theta_rad = mid_theta_rad % (2 * np.pi)
                    text_x = x1 + arc_radius * 1.3 * np.cos(mid_theta_rad)
                    text_y = y1 + arc_radius * 1.3 * np.sin(mid_theta_rad)
                    theta_text = ax.text(text_x, text_y,
                                        f'θ = {theta:.1f}° (N→E)',
                                        color='white', fontsize=11, weight='bold',
                                        ha='center', va='center')
                    
                    # Mettre à jour les champs
                    self.star2_x_var.set(f"{refined_x:.2f}")
                    self.star2_y_var.set(f"{refined_y:.2f}")
                    
                    selection_step = 3
                    update_instruction()
                    if rho_value is not None:
                        logger.info(f"Étoile 2 sélectionnée: ({refined_x:.2f}, {refined_y:.2f}), Angle theta: {theta:.2f}°, Rho: {separation:.2f} px = {rho_value:.3f}\"")
                    else:
                        logger.info(f"Étoile 2 sélectionnée: ({refined_x:.2f}, {refined_y:.2f}), Angle theta: {theta:.2f}°, Rho: {separation:.2f} px")
                
                canvas.draw()
            
            # Connecter le gestionnaire de clic
            canvas.mpl_connect('button_press_event', on_click)
            update_instruction()
            
            # Bouton de réinitialisation
            button_frame = ttk.Frame(main_frame)
            button_frame.pack(fill="x", pady=(3, 0))
            
            def reset_selection():
                nonlocal marker_star1, marker_star2, selection_step
                nonlocal theta_arrow_north, theta_arrow_east, theta_arrow_separation, theta_arc, theta_text
                nonlocal compass_north, compass_east, separation_line, rho_text
                if marker_star1:
                    marker_star1.remove()
                    marker_star1 = None
                if marker_star2:
                    marker_star2.remove()
                    marker_star2 = None
                # Nettoyer la ligne de séparation
                if separation_line:
                    separation_line.remove()
                    separation_line = None
                # Nettoyer les autres lignes
                for line in list(ax.lines):
                    if line.get_label() in ['Séparation T1-T2']:
                        line.remove()
                # Nettoyer les textes
                for txt in list(ax.texts):
                    txt_text = txt.get_text()
                    if txt_text in ['Star 1', 'Star 2', 'N', 'E'] or 'ρ' in txt_text or 'θ' in txt_text:
                        txt.remove()
                # Nettoyer les éléments de l'angle theta
                for elem in [theta_arrow_north, theta_arrow_east, theta_arrow_separation, theta_arc, theta_text]:
                    if elem is not None:
                        try:
                            elem.remove()
                        except:
                            pass
                        if elem == theta_arrow_north:
                            theta_arrow_north = None
                        elif elem == theta_arrow_east:
                            theta_arrow_east = None
                        elif elem == theta_arrow_separation:
                            theta_arrow_separation = None
                        elif elem == theta_arc:
                            theta_arc = None
                        elif elem == theta_text:
                            theta_text = None
                # Nettoyer les éléments du compas
                if compass_north:
                    try:
                        compass_north.remove()
                    except (ValueError, AttributeError):
                        pass
                    compass_north = None
                if compass_east:
                    try:
                        compass_east.remove()
                    except (ValueError, AttributeError):
                        pass
                    compass_east = None
                if rho_text:
                    try:
                        rho_text.remove()
                    except (ValueError, AttributeError):
                        pass
                    rho_text = None
                # Nettoyer les patches (arcs et flèches)
                from matplotlib.patches import Arc, FancyArrowPatch
                for patch in list(ax.patches):
                    if isinstance(patch, (Arc, FancyArrowPatch)):
                        patch.remove()
                star1_pos.update({'x': None, 'y': None, 'refined_x': None, 'refined_y': None})
                star2_pos.update({'x': None, 'y': None, 'refined_x': None, 'refined_y': None})
                self.star1_x_var.set("")
                self.star1_y_var.set("")
                self.star2_x_var.set("")
                self.star2_y_var.set("")
                selection_step = 1
                update_instruction()
                canvas.draw()
            
            ttk.Button(
                button_frame,
                text="🔄 Réinitialiser",
                command=reset_selection
            ).pack(side="left", padx=(0, 10))
            
            ttk.Button(
                button_frame,
                text="Fermer",
                command=measure_window.destroy
            ).pack(side="left")
            
            logger.info(f"Affichage de l'image pour mesure: {image_path.name}")
            
        except Exception as e:
            logger.error(f"[MESURE] Erreur lors du chargement de l'image: {e}", exc_info=True)
    
    def analyze_and_generate_report(self):
        """Effectue l'analyse et génère un rapport au format catalogue"""
        if not REDUCTION_AVAILABLE or self.reducer is None:
            logger.error("[ANALYSE] Module de réduction non disponible")
            return
        
        image_path = self.measure_image_var.get()
        if not image_path:
            messagebox.showwarning("Attention", "Sélectionnez d'abord une image!")
            return
        
        star1_x_str = self.star1_x_var.get()
        star1_y_str = self.star1_y_var.get()
        star2_x_str = self.star2_x_var.get()
        star2_y_str = self.star2_y_var.get()
        pixel_scale_str = self.pixel_scale_var.get()
        
        if not all([star1_x_str, star1_y_str, star2_x_str, star2_y_str]):
            messagebox.showwarning("Attention", "Sélectionnez les deux étoiles avant d'analyser!")
            return
        
        try:
            star1_pos = (float(star1_x_str), float(star1_y_str))
            star2_pos = (float(star2_x_str), float(star2_y_str))
            pixel_scale = float(pixel_scale_str)
        except ValueError:
            logger.error(f"[ANALYSE_SEP] Les valeurs doivent être des nombres (star1=({star1_x_str}, {star1_y_str}), star2=({star2_x_str}, {star2_y_str}), scale={pixel_scale_str})")
            return
        
        def analyze_task():
            try:
                self.progress.start()
                
                result = self.reducer.measure_binary_separation(
                    Path(image_path),
                    star1_pos,
                    star2_pos,
                    pixel_scale
                )
                
                self.progress.stop()
                
                # Afficher les résultats
                self.measure_result_text.delete(1.0, tk.END)
                result_text = "Résultats de la mesure:\n"
                result_text += f"{'='*60}\n\n"
                result_text += f"Étoile 1 (centroïde précis):\n"
                result_text += f"  x = {result['x1']:.3f} pixels\n"
                result_text += f"  y = {result['y1']:.3f} pixels\n\n"
                result_text += f"Étoile 2 (centroïde précis):\n"
                result_text += f"  x = {result['x2']:.3f} pixels\n"
                result_text += f"  y = {result['y2']:.3f} pixels\n\n"
                result_text += f"{'-'*60}\n\n"
                result_text += f"Séparation: {result['separation_pix']:.3f} pixels\n"
                result_text += f"Séparation: {result['separation_arcsec']:.3f} arcsec\n"
                result_text += f"Angle de position: {result['position_angle']:.2f}°\n\n"
                result_text += f"(Angle mesuré depuis le nord vers l'est)\n"
                
                self.measure_result_text.insert(1.0, result_text)
                
                # Générer le rapport
                report_path = self.generate_separation_report(
                    Path(image_path),
                    result,
                    pixel_scale,
                    star1_pos,
                    star2_pos
                )
                
                messagebox.showinfo(
                    "Analyse terminée",
                    f"Analyse terminée avec succès!\n\n"
                    f"Séparation: {result['separation_arcsec']:.3f} arcsec\n"
                    f"Angle de position: {result['position_angle']:.2f}°\n\n"
                    f"Rapport sauvegardé:\n{report_path}"
                )
                
                logger.info(f"Mesure séparation: {result['separation_arcsec']:.3f}\", PA: {result['position_angle']:.2f}°")
                
            except Exception as e:
                self.progress.stop()
                logger.error(f"[ANALYSE_SEP] Erreur lors de l'analyse: {e}", exc_info=True)
        
        threading.Thread(target=analyze_task, daemon=True).start()
    
    def generate_separation_report(self, image_path: Path, result: dict, pixel_scale: float,
                                   star1_approx: tuple, star2_approx: tuple) -> Path:
        """
        Génère un rapport de mesure de séparation au format Washington Double Star Catalog (WDS)
        """
        from datetime import datetime
        from astropy.wcs import WCS
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        
        # Charger l'en-tête FITS pour récupérer les métadonnées
        ra_center = None
        dec_center = None
        epoch_year = None
        
        try:
            with fits.open(image_path) as hdul:
                header = hdul[0].header
                # Récupérer la date d'observation si disponible
                date_obs = header.get('DATE-OBS', datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                # Extraire l'année pour l'époque
                try:
                    if isinstance(date_obs, str):
                        epoch_year = float(date_obs.split('-')[0])
                    else:
                        epoch_year = datetime.now().year
                except:
                    epoch_year = datetime.now().year
                
                # Récupérer le filtre si disponible
                filter_name = header.get('FILTER', header.get('FILT', 'Unknown'))
                # Récupérer le temps d'exposition si disponible
                exptime = header.get('EXPTIME', header.get('EXPOSURE', 'Unknown'))
                
                # Essayer d'extraire les coordonnées du WCS
                try:
                    wcs = WCS(header)
                    if wcs.has_celestial:
                        # Calculer les coordonnées du centre de l'image
                        if wcs.naxis == 2:
                            center_pix = (header.get('NAXIS1', 0) / 2, header.get('NAXIS2', 0) / 2)
                        else:
                            center_pix = (0, 0)
                        center_coord = wcs.pixel_to_world_values(*center_pix)
                        ra_center = center_coord[0]
                        dec_center = center_coord[1]
                except:
                    # Essayer d'extraire depuis l'en-tête directement
                    ra_center = header.get('CRVAL1', header.get('RA', None))
                    dec_center = header.get('CRVAL2', header.get('DEC', None))
                    
        except:
            date_obs = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            epoch_year = datetime.now().year
            filter_name = "Unknown"
            exptime = "Unknown"
        
        # Générer le nom du fichier rapport
        report_name = f"WDS_report_{image_path.stem}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        report_path = image_path.parent / report_name
        
        # Formater les coordonnées pour WDS (format strict)
        # Colonnes 1-10: HHMMSS.s+DDMMSS (format arcminute)
        # Colonnes 113-130: HHMMSS.ss+DDMMSS (format arcsecond - plus précis)
        ra_str = "Unknown"  # Pour affichage détaillé
        dec_str = "Unknown"  # Pour affichage détaillé
        ra_str_arcmin = "          "  # 10 caractères pour format WDS compact
        ra_str_arcsec = "                  "  # 18 caractères pour format WDS compact
        if ra_center is not None and dec_center is not None:
            try:
                coord = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame='icrs')
                # Format pour affichage détaillé
                ra_str = coord.ra.to_string(unit=u.hourangle, precision=1, pad=True)
                dec_str = coord.dec.to_string(unit=u.deg, precision=0, pad=True, alwayssign=True)
                
                # Format arcminute (col 1-10): HHMMSS.s+DDMMSS (10 caractères exactement)
                # Format: HHMMSS.s + DDMMSS = 7 + 1 + 6 = 14 caractères, mais limité à 10
                # En fait, le format WDS combine: HHMMSS.s + DDMMSS dans 10 caractères
                ra_hms = coord.ra.hms
                dec_dms = coord.dec.dms
                sign = '+' if dec_dms.d >= 0 else '-'
                # Format WDS arcminute: HHMMSS.s+DDMMSS (10 caractères)
                ra_part = f"{int(ra_hms.h):02d}{int(ra_hms.m):02d}{ra_hms.s:04.1f}"
                dec_part = f"{sign}{int(abs(dec_dms.d)):02d}{int(abs(dec_dms.m)):02d}{int(abs(dec_dms.s)):02d}"
                # Combiner et tronquer à 10 caractères
                combined = ra_part + dec_part
                ra_str_arcmin = combined[:10].ljust(10)
                
                # Format arcsecond (col 113-130): HHMMSS.ss+DDMMSS (18 caractères exactement)
                # Format: HHMMSS.ss + DDMMSS = 8 + 6 = 14 caractères, mais on peut mettre le signe
                ra_part_sec = f"{int(ra_hms.h):02d}{int(ra_hms.m):02d}{ra_hms.s:05.2f}"
                dec_part_sec = f"{sign}{int(abs(dec_dms.d)):02d}{int(abs(dec_dms.m)):02d}{int(abs(dec_dms.s)):02d}"
                # Combiner et pad à 18 caractères
                combined_sec = ra_part_sec + dec_part_sec
                ra_str_arcsec = combined_sec[:18].ljust(18)
            except Exception as e:
                logger.warning(f"Erreur formatage coordonnées WDS: {e}")
                pass
        
        # Générer le contenu du rapport au format WDS
        report_lines = []
        report_lines.append("="*100)
        report_lines.append("WASHINGTON DOUBLE STAR CATALOG (WDS) FORMAT REPORT")
        report_lines.append("="*100)
        report_lines.append("")
        report_lines.append(f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"Source image: {image_path.name}")
        report_lines.append("")
        report_lines.append("-"*100)
        report_lines.append("WDS STANDARD FORMAT")
        report_lines.append("-"*100)
        report_lines.append("")
        
        # Format WDS standard (colonnes fixes)
        report_lines.append("WDS Identification Format:")
        report_lines.append("")
        report_lines.append(f"  RA (J2000):     {ra_str}")
        report_lines.append(f"  Dec (J2000):    {dec_str}")
        report_lines.append(f"  Epoch:          {epoch_year:.1f}")
        report_lines.append("")
        report_lines.append("  Position Angle (theta):  {:.1f} deg".format(result['position_angle']))
        report_lines.append("  Separation (rho):         {:.3f} arcsec".format(result['separation_arcsec']))
        report_lines.append("")
        
        # Informations additionnelles
        report_lines.append("-"*100)
        report_lines.append("OBSERVATION DETAILS")
        report_lines.append("-"*100)
        report_lines.append("")
        report_lines.append(f"  Observation Date:        {date_obs}")
        report_lines.append(f"  Filter:                  {filter_name}")
        report_lines.append(f"  Exposure Time:           {exptime}")
        report_lines.append(f"  Pixel Scale:             {pixel_scale:.4f} arcsec/pixel")
        report_lines.append("")
        
        # Positions précises
        report_lines.append("-"*100)
        report_lines.append("STAR POSITIONS (pixels)")
        report_lines.append("-"*100)
        report_lines.append("")
        report_lines.append(f"  Primary Star (A):")
        report_lines.append(f"    Approximate position:  ({star1_approx[0]:.2f}, {star1_approx[1]:.2f}) px")
        report_lines.append(f"    Refined centroid:      ({result['x1']:.3f}, {result['y1']:.3f}) px")
        report_lines.append("")
        report_lines.append(f"  Secondary Star (B):")
        report_lines.append(f"    Approximate position:  ({star2_approx[0]:.2f}, {star2_approx[1]:.2f}) px")
        report_lines.append(f"    Refined centroid:      ({result['x2']:.3f}, {result['y2']:.3f}) px")
        report_lines.append("")
        
        # Mesures détaillées
        report_lines.append("-"*100)
        report_lines.append("MEASUREMENT DATA")
        report_lines.append("-"*100)
        report_lines.append("")
        report_lines.append(f"  Separation (pixels):     {result['separation_pix']:.4f} px")
        report_lines.append(f"  Separation (arcsec):     {result['separation_arcsec']:.4f} arcsec")
        report_lines.append(f"  Position Angle:          {result['position_angle']:.2f} deg")
        report_lines.append("")
        report_lines.append("  Notes:")
        report_lines.append("    - Position angle measured from North (0°) toward East (90°)")
        report_lines.append("    - Position angle follows standard astronomical convention")
        report_lines.append("    - Measurements made using sub-pixel centroid refinement")
        report_lines.append("")
        
        # Format WDS compact (format ligne strict - 130 colonnes)
        report_lines.append("-"*100)
        report_lines.append("WDS COMPACT FORMAT (Single Line - 130 columns)")
        report_lines.append("-"*100)
        report_lines.append("")
        
        # Format strict selon spécifications WDS
        # Colonnes 1-10: RA Dec (arcminute) - format HHMMSS.s+DDMMSS
        coords_arcmin = ra_str_arcmin[:10].ljust(10)
        
        # Colonnes 11-17: Discoverer & Number (max 7 caractères)
        discoverer = "Custom".ljust(7)[:7]
        
        # Colonnes 18-22: Components (max 5 caractères, ex: "AB  ")
        components = "AB   "[:5].ljust(5)
        
        # Colonnes 24-27: Date (first) - I4 (4 chiffres)
        date_first = f"{int(epoch_year):4d}"
        
        # Colonnes 29-32: Date (last) - I4 (4 chiffres)
        date_last = f"{int(epoch_year):4d}"
        
        # Colonnes 34-37: Number of Observations - I4 (max 9999)
        n_obs = "   1"[:4].rjust(4)
        
        # Colonnes 39-41: Position Angle (first) - I3 (000-359)
        pa_first = f"{int(result['position_angle']):03d}"[:3]
        
        # Colonnes 43-45: Position Angle (last) - I3 (000-359)
        pa_last = f"{int(result['position_angle']):03d}"[:3]
        
        # Colonnes 47-51: Separation (first) - F5.1 (5 caractères, 1 décimale)
        sep_first = f"{result['separation_arcsec']:5.1f}"[:5]
        
        # Colonnes 53-57: Separation (last) - F5.1 (5 caractères, 1 décimale)
        sep_last = f"{result['separation_arcsec']:5.1f}"[:5]
        
        # Colonnes 59-63: Magnitude of First Component - F5.2 (5 caractères, 2 décimales)
        mag1 = "     "[:5].ljust(5)  # Non disponible
        
        # Colonnes 65-69: Magnitude of Second Component - F5.2 (5 caractères, 2 décimales)
        mag2 = "     "[:5].ljust(5)  # Non disponible
        
        # Colonnes 71-79: Spectral Type - A9 (9 caractères)
        spec_type = "         "[:9].ljust(9)  # Non disponible
        
        # Colonnes 81-84: Primary Proper Motion (RA) - I4
        pm_ra1 = "    "[:4].rjust(4)  # Non disponible
        
        # Colonnes 85-88: Primary Proper Motion (Dec) - I4
        pm_dec1 = "    "[:4].rjust(4)  # Non disponible
        
        # Colonnes 90-93: Secondary Proper Motion (RA) - I4
        pm_ra2 = "    "[:4].rjust(4)  # Non disponible
        
        # Colonnes 94-97: Secondary Proper Motion (Dec) - I4
        pm_dec2 = "    "[:4].rjust(4)  # Non disponible
        
        # Colonnes 99-106: Durchmusterung Number - A8 (8 caractères)
        durchmusterung = "        "[:8].ljust(8)  # Non disponible
        
        # Colonnes 108-111: Notes - A4 (4 caractères)
        notes = "    "[:4].ljust(4)  # Vide
        
        # Colonnes 113-130: 2000 arcsecond coordinates - A18 (format HHMMSS.ss+DDMMSS)
        coords_arcsec = ra_str_arcsec[:18].ljust(18)
        
        # Construire la ligne WDS (130 caractères exactement)
        wds_line = (
            coords_arcmin +          # Col 1-10: 10 chars
            discoverer +              # Col 11-17: 7 chars
            components +              # Col 18-22: 5 chars
            " " +                     # Col 23: espace
            date_first +              # Col 24-27: 4 chars
            " " +                     # Col 28: espace
            date_last +               # Col 29-32: 4 chars
            " " +                     # Col 33: espace
            n_obs +                   # Col 34-37: 4 chars
            " " +                     # Col 38: espace
            pa_first +                # Col 39-41: 3 chars
            " " +                     # Col 42: espace
            pa_last +                 # Col 43-45: 3 chars
            " " +                     # Col 46: espace
            sep_first +               # Col 47-51: 5 chars
            " " +                     # Col 52: espace
            sep_last +                # Col 53-57: 5 chars
            " " +                     # Col 58: espace
            mag1 +                    # Col 59-63: 5 chars
            " " +                     # Col 64: espace
            mag2 +                    # Col 65-69: 5 chars
            " " +                     # Col 70: espace
            spec_type +               # Col 71-79: 9 chars
            " " +                     # Col 80: espace
            pm_ra1 +                  # Col 81-84: 4 chars
            pm_dec1 +                 # Col 85-88: 4 chars
            " " +                     # Col 89: espace
            pm_ra2 +                  # Col 90-93: 4 chars
            pm_dec2 +                 # Col 94-97: 4 chars
            " " +                     # Col 98: espace
            durchmusterung +          # Col 99-106: 8 chars
            " " +                     # Col 107: espace
            notes +                   # Col 108-111: 4 chars
            " " +                     # Col 112: espace
            coords_arcsec             # Col 113-130: 18 chars
        )
        
        # Vérifier la longueur (doit être exactement 130 caractères)
        if len(wds_line) != 130:
            logger.warning(f"Ligne WDS longueur incorrecte: {len(wds_line)} (attendu 130)")
            wds_line = wds_line[:130].ljust(130)
        
        report_lines.append(wds_line)
        report_lines.append("")
        report_lines.append(f"Line length: {len(wds_line)} characters (should be 130)")
        report_lines.append("")
        report_lines.append("="*100)
        
        # Écrire le rapport
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_lines))
        
        logger.info(f"Rapport WDS sauvegardé: {report_path}")
        return report_path
