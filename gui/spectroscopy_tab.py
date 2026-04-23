# gui/spectroscopy_tab.py
"""
Onglet pour l'analyse de spectres d'étoiles avec specutils
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import logging
import numpy as np
from pathlib import Path

from gui.manual_help import add_manual_help_header

logger = logging.getLogger(__name__)

# Import specutils
try:
    from specutils import Spectrum1D
    # Essayer d'importer Spectrum pour compatibilité (anciennes versions)
    try:
        from specutils import Spectrum
    except ImportError:
        Spectrum = Spectrum1D  # Alias pour compatibilité
    
    from specutils.fitting import fit_generic_continuum
    from specutils.analysis import (
        equivalent_width, line_flux, centroid, 
        fwhm, gaussian_sigma_width, gaussian_fwhm
    )
    from specutils import SpectralRegion
    from astropy import units as u
    from astropy.io import fits
    from astropy.visualization import quantity_support
    SPECUTILS_AVAILABLE = True
except ImportError as e:
    SPECUTILS_AVAILABLE = False
    logger.warning(f"specutils n'est pas installé ou erreur d'import: {e}. Utilisez: pip install specutils")

# Import Prospector (optionnel, pour analyse de galaxies)
logger.info("=" * 60)
logger.info("GUI/SPECTROSCOPY_TAB: Début import Prospector")
logger.info("=" * 60)
try:
    logger.info("[GUI] Tentative d'import depuis core.prospector_analysis...")
    from core.prospector_analysis import ProspectorAnalyzer, PROSPECTOR_AVAILABLE as PROSP_AVAIL
    PROSPECTOR_AVAILABLE = PROSP_AVAIL
    logger.info(f"[GUI] ✓ Import réussi - PROSPECTOR_AVAILABLE = {PROSPECTOR_AVAILABLE}")
    logger.info(f"[GUI] ProspectorAnalyzer disponible: {ProspectorAnalyzer is not None}")
except ImportError as e:
    PROSPECTOR_AVAILABLE = False
    logger.error(f"[GUI] ✗ Erreur ImportError lors de l'import: {e}")
    logger.error(f"[GUI] Type d'erreur: {type(e).__name__}")
    import traceback
    logger.error(f"[GUI] Traceback:\n{traceback.format_exc()}")
    ProspectorAnalyzer = None
except Exception as e:
    PROSPECTOR_AVAILABLE = False
    logger.error(f"[GUI] ✗ Erreur inattendue lors de l'import: {e}")
    logger.error(f"[GUI] Type d'erreur: {type(e).__name__}")
    import traceback
    logger.error(f"[GUI] Traceback:\n{traceback.format_exc()}")
    ProspectorAnalyzer = None
logger.info(f"[GUI] PROSPECTOR_AVAILABLE final = {PROSPECTOR_AVAILABLE}")
logger.info("=" * 60)

# Import matplotlib
try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("matplotlib n'est pas disponible")


class SpectroscopyTab(ttk.Frame):
    """
    Onglet pour l'analyse de spectres d'étoiles
    """
    
    def __init__(self, parent_notebook, base_dir=None):
        super().__init__(parent_notebook, padding=10)
        
        if base_dir is None:
            self.base_dir = Path.home()
        else:
            self.base_dir = Path(base_dir)
        
        # Données du spectre
        self.spectrum = None
        self.spectrum_normalized = None
        
        # Prospector (pour analyse de galaxies)
        logger.info("[GUI] Initialisation de l'analyseur Prospector...")
        logger.info(f"[GUI] PROSPECTOR_AVAILABLE = {PROSPECTOR_AVAILABLE}")
        logger.info(f"[GUI] ProspectorAnalyzer = {ProspectorAnalyzer}")
        self.prospector_analyzer = None
        self.sed_data = None  # Données SED pour Prospector
        if PROSPECTOR_AVAILABLE and ProspectorAnalyzer is not None:
            try:
                logger.info("[GUI] Tentative de création de ProspectorAnalyzer...")
                self.prospector_analyzer = ProspectorAnalyzer()
                logger.info("[GUI] ✓ Analyseur Prospector initialisé avec succès")
            except Exception as e:
                logger.error(f"[GUI] ✗ Impossible d'initialiser Prospector: {e}")
                logger.error(f"[GUI] Type d'erreur: {type(e).__name__}")
                import traceback
                logger.error(f"[GUI] Traceback:\n{traceback.format_exc()}")
                self.prospector_analyzer = None
        else:
            logger.warning(f"[GUI] ⚠ ProspectorAnalyzer non initialisé - PROSPECTOR_AVAILABLE={PROSPECTOR_AVAILABLE}, ProspectorAnalyzer={ProspectorAnalyzer}")
        
        # Initialiser quantity_support pour les unités sur les axes
        if SPECUTILS_AVAILABLE:
            try:
                quantity_support()
            except:
                pass
        
        self.create_widgets()
        
        if not SPECUTILS_AVAILABLE:
            self.show_specutils_install_message()
    
    def create_widgets(self):
        """Crée l'interface utilisateur"""
        add_manual_help_header(self, "11-spectroscopie")

        # En-tête
        header_frame = ttk.Frame(self)
        header_frame.pack(fill="x", pady=(0, 10))
        
        title_label = ttk.Label(
            header_frame,
            text="Analyse de Spectres d'Étoiles",
            font=("Helvetica", 14, "bold")
        )
        title_label.pack(side="left")
        
        # Frame principal avec PanedWindow pour diviser l'interface
        main_paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True)
        
        # === PARTIE GAUCHE : CONTRÔLES ===
        left_frame = ttk.Frame(main_paned, padding=10)
        main_paned.add(left_frame, weight=1)
        
        # === PARTIE DROITE : VISUALISATION ===
        right_frame = ttk.Frame(main_paned, padding=10)
        main_paned.add(right_frame, weight=2)
        
        # Construction des sections
        self.create_controls_section(left_frame)
        self.create_visualization_section(right_frame)
    
    def create_controls_section(self, parent):
        """Crée la section de contrôles à gauche"""
        
        # 1. Installation/Statut specutils
        status_frame = ttk.LabelFrame(parent, text="1. Statut", padding=10)
        status_frame.pack(fill="x", pady=5)
        
        # Statut specutils
        self.status_label = ttk.Label(
            status_frame,
            text="État: Vérification...",
            foreground="blue"
        )
        self.status_label.pack(pady=2)
        
        if SPECUTILS_AVAILABLE:
            try:
                import specutils
                version = getattr(specutils, '__version__', 'installe')
                self.status_label.config(
                    text=f"✓ specutils {version} installé",
                    foreground="green"
                )
            except:
                self.status_label.config(
                    text="✓ specutils installé",
                    foreground="green"
                )
        else:
            self.status_label.config(
                text="✗ specutils non installé",
                foreground="red"
            )
            ttk.Button(
                status_frame,
                text="Installer specutils",
                command=self.install_specutils
            ).pack(pady=2)
        
        # Statut Prospector
        self.prospector_status_label = ttk.Label(
            status_frame,
            text="",
            foreground="blue"
        )
        self.prospector_status_label.pack(pady=2)
        
        logger.info(f"[GUI] Affichage du statut Prospector - PROSPECTOR_AVAILABLE = {PROSPECTOR_AVAILABLE}")
        if PROSPECTOR_AVAILABLE:
            try:
                logger.info("[GUI] Vérification de la version de prospect...")
                import prospect
                version = getattr(prospect, '__version__', 'installé')
                logger.info(f"[GUI] Version de prospect: {version}")
                status_text = f"✓ Prospector {version} installé"
                
                # Vérifier FSPS
                try:
                    logger.info("[GUI] Vérification de FSPS_AVAILABLE...")
                    from core.prospector_analysis import FSPS_AVAILABLE
                    logger.info(f"[GUI] FSPS_AVAILABLE = {FSPS_AVAILABLE}")
                    if FSPS_AVAILABLE:
                        status_text += " (FSPS disponible)"
                        self.prospector_status_label.config(
                            text=status_text,
                            foreground="green"
                        )
                        logger.info("[GUI] ✓ Statut Prospector affiché: installé avec FSPS")
                    else:
                        status_text += " (sans FSPS - utilise fichiers stub)"
                        self.prospector_status_label.config(
                            text=status_text,
                            foreground="orange"
                        )
                        logger.info("[GUI] ⚠ Statut Prospector affiché: installé sans FSPS (fichiers stub disponibles)")
                except Exception as e:
                    logger.warning(f"[GUI] Erreur lors de la vérification de FSPS: {e}")
                    self.prospector_status_label.config(
                        text=status_text,
                        foreground="green"
                    )
                    logger.info("[GUI] ✓ Statut Prospector affiché: installé")
            except Exception as e:
                logger.error(f"[GUI] ✗ Erreur lors de la vérification de prospect: {e}")
                import traceback
                logger.error(f"[GUI] Traceback:\n{traceback.format_exc()}")
                self.prospector_status_label.config(
                    text="✓ Prospector installé",
                    foreground="green"
                )
        else:
            logger.warning("[GUI] ⚠ PROSPECTOR_AVAILABLE = False, affichage du statut 'non installé'")
            self.prospector_status_label.config(
                text="✗ Prospector non installé",
                foreground="red"
            )
            ttk.Button(
                status_frame,
                text="Installer Prospector",
                command=self.install_prospector
            ).pack(pady=2)
        
        # 2. Chargement du spectre
        load_frame = ttk.LabelFrame(parent, text="2. Charger un spectre", padding=10)
        load_frame.pack(fill="x", pady=5)
        
        self.file_var = tk.StringVar()
        
        file_frame = ttk.Frame(load_frame)
        file_frame.pack(fill="x", pady=2)
        ttk.Label(file_frame, text="Fichier:").pack(side="left", padx=(0, 5))
        ttk.Entry(file_frame, textvariable=self.file_var, width=30).pack(side="left", fill="x", expand=True)
        ttk.Button(file_frame, text="Parcourir", command=self.browse_file).pack(side="left", padx=(5, 0))
        
        ttk.Label(
            load_frame, 
            text="Formats supportés: FITS, ASCII", 
            font=("Helvetica", 8), 
            foreground="gray"
        ).pack(anchor="w", pady=2)
        
        ttk.Button(
            load_frame,
            text="Charger le spectre",
            command=self.load_spectrum
        ).pack(pady=5)
        
        # 3. Analyse du spectre
        analysis_frame = ttk.LabelFrame(parent, text="3. Analyse", padding=10)
        analysis_frame.pack(fill="x", pady=5)
        
        # Normalisation du continuum
        ttk.Button(
            analysis_frame,
            text="Normaliser le continuum",
            command=self.normalize_continuum
        ).pack(pady=2, fill="x")
        
        # Analyse de raie
        line_frame = ttk.Frame(analysis_frame)
        line_frame.pack(fill="x", pady=5)
        
        ttk.Label(line_frame, text="Région de raie (Å):").pack(anchor="w")
        
        line_input_frame = ttk.Frame(line_frame)
        line_input_frame.pack(fill="x", pady=2)
        
        ttk.Label(line_input_frame, text="Min:").pack(side="left", padx=(0, 5))
        self.line_min_var = tk.StringVar(value="6560")
        ttk.Entry(line_input_frame, textvariable=self.line_min_var, width=10).pack(side="left", padx=2)
        
        ttk.Label(line_input_frame, text="Max:").pack(side="left", padx=(10, 5))
        self.line_max_var = tk.StringVar(value="6580")
        ttk.Entry(line_input_frame, textvariable=self.line_max_var, width=10).pack(side="left", padx=2)
        
        ttk.Button(
            analysis_frame,
            text="Analyser la raie",
            command=self.analyze_line
        ).pack(pady=2, fill="x")
        
        # Séparateur
        ttk.Separator(analysis_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Section Prospector (analyse de galaxies) - Toujours affichée
        prospector_label_frame = ttk.Frame(analysis_frame)
        prospector_label_frame.pack(anchor="w", pady=(5, 2))
        
        ttk.Label(
            prospector_label_frame,
            text="Analyse de Galaxies (Prospector)",
            font=("Helvetica", 9, "bold")
        ).pack(side="left")
        
        if not PROSPECTOR_AVAILABLE:
            ttk.Label(
                prospector_label_frame,
                text=" [Non installé]",
                font=("Helvetica", 8),
                foreground="red"
            ).pack(side="left", padx=(5, 0))
        
        if PROSPECTOR_AVAILABLE:
            ttk.Button(
                analysis_frame,
                text="🌌 Inférer Propriétés Stellaires (SED)",
                command=self.analyze_galaxy_with_prospector
            ).pack(pady=2, fill="x")
            
            ttk.Button(
                analysis_frame,
                text="📊 Créer SED depuis Photométrie",
                command=self.create_sed_from_photometry
            ).pack(pady=2, fill="x")
        else:
            ttk.Label(
                analysis_frame,
                text="Installez Prospector pour activer ces fonctionnalités",
                font=("Helvetica", 8),
                foreground="gray"
            ).pack(pady=5)
        
        # 4. Informations du spectre
        info_frame = ttk.LabelFrame(parent, text="4. Informations", padding=10)
        info_frame.pack(fill="x", pady=5)
        
        self.info_text = tk.Text(info_frame, height=8, width=30, wrap=tk.WORD, font=("Courier", 9))
        info_scroll = ttk.Scrollbar(info_frame, orient="vertical", command=self.info_text.yview)
        self.info_text.config(yscrollcommand=info_scroll.set)
        
        info_frame_inner = ttk.Frame(info_frame)
        info_frame_inner.pack(fill="both", expand=True)
        self.info_text.pack(side="left", fill="both", expand=True)
        info_scroll.pack(side="right", fill="y")
        
        # 5. Export
        export_frame = ttk.LabelFrame(parent, text="5. Export", padding=10)
        export_frame.pack(fill="x", pady=5)
        
        ttk.Button(
            export_frame,
            text="Exporter le spectre",
            command=self.export_spectrum
        ).pack(pady=2, fill="x")
        
        ttk.Button(
            export_frame,
            text="Exporter le graphique",
            command=self.export_plot
        ).pack(pady=2, fill="x")
    
    def create_visualization_section(self, parent):
        """Crée la section de visualisation à droite"""
        
        if not MATPLOTLIB_AVAILABLE:
            ttk.Label(
                parent,
                text="matplotlib n'est pas disponible",
                foreground="red"
            ).pack()
            return
        
        # Figure matplotlib
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        
        # Canvas
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, parent)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def show_specutils_install_message(self):
        """Affiche un message si specutils n'est pas installé"""
        messagebox.showwarning(
            "specutils non installé",
            "specutils n'est pas installé.\n\n"
            "Pour installer:\n"
            "pip install specutils\n\n"
            "Certaines fonctionnalités seront désactivées."
        )
    
    def install_specutils(self):
        """Tente d'installer specutils"""
        import threading
        import subprocess
        
        def install():
            try:
                result = subprocess.run(
                    ["pip", "install", "specutils"],
                    capture_output=True,
                    text=True,
                    timeout=300
                )
                if result.returncode == 0:
                    messagebox.showinfo("Succès", "specutils a été installé avec succès!\nRedémarrez l'application.")
                else:
                    messagebox.showerror("Erreur", f"Erreur lors de l'installation:\n{result.stderr}")
            except Exception as e:
                messagebox.showerror("Erreur", f"Erreur lors de l'installation: {e}")
        
        thread = threading.Thread(target=install, daemon=True)
        thread.start()
        messagebox.showinfo("Installation", "Installation de specutils en cours...\nVeuillez patienter.")
    
    def install_prospector(self):
        """Tente d'installer Prospector et ses dépendances"""
        import threading
        import subprocess
        import sys
        
        def install():
            logger.info("=" * 60)
            logger.info("[INSTALL] Début de l'installation de Prospector")
            logger.info("=" * 60)
            logger.info(f"[INSTALL] Python: {sys.executable}")
            logger.info(f"[INSTALL] Version Python: {sys.version}")
            try:
                # Installer les dépendances de Prospector d'abord
                # IMPORTANT: sedpy doit être installé depuis GitHub car PyPI n'a pas le module observate
                dependencies = [
                    ("git+https://github.com/bd-j/sedpy.git", "sedpy (depuis GitHub)"),
                    "dynesty>=2.0.0",
                    "dill>=0.3.0",
                    "h5py>=3.0.0",
                    "emcee>=3.1.0",
                    "corner>=2.2.3",
                    "matplotlib>=3.8.4",
                ]
                
                messagebox.showinfo(
                    "Installation Prospector",
                    "Installation de Prospector depuis GitHub...\n\n"
                    "IMPORTANT: Prospector doit être installé depuis GitHub\n"
                    "(le package 'prospector' sur PyPI est un autre outil).\n\n"
                    "Dépendances à installer:\n"
                    "- sedpy (depuis GitHub)\n"
                    "- dynesty, dill, h5py, emcee, corner, matplotlib\n"
                    "- Prospector depuis GitHub\n\n"
                    "Sous Windows sans compilateur Fortran (gfortran), la construction "
                    "du paquet « fsps » échoue souvent : un second essai sans dépendances "
                    "pip sera lancé automatiquement (stubs FSPS, comme l’installateur .ps1).\n\n"
                    "Cela peut prendre plusieurs minutes.\n"
                    "Git doit être installé : https://git-scm.com/download/win"
                )
                
                # Installer les dépendances
                logger.info(f"[INSTALL] Installation des dépendances: {dependencies}")
                for package_item in dependencies:
                    # package_item peut être une chaîne simple ou un tuple (URL, nom_affiché)
                    if isinstance(package_item, tuple):
                        package_url, package_name = package_item
                        package_display = package_name
                    else:
                        package_url = package_item
                        package_display = package_item
                    
                    logger.info(f"[INSTALL] Installation de {package_display}...")
                    result = subprocess.run(
                        [sys.executable, "-m", "pip", "install", package_url],
                        capture_output=True,
                        text=True,
                        timeout=600
                    )
                    logger.info(f"[INSTALL] Retour code: {result.returncode}")
                    logger.info(f"[INSTALL] Sortie stdout:\n{result.stdout}")
                    if result.stderr:
                        logger.warning(f"[INSTALL] Sortie stderr:\n{result.stderr}")
                    if result.returncode != 0:
                        logger.error(f"[INSTALL] ✗ Erreur lors de l'installation de {package_display}")
                        logger.error(f"[INSTALL] stderr: {result.stderr}")
                        messagebox.showerror(
                            "Erreur",
                            f"Erreur lors de l'installation de {package_display}:\n{result.stderr}"
                        )
                        return
                    logger.info(f"[INSTALL] ✓ {package_display} installé avec succès")
                
                # Installer Prospector depuis GitHub (même logique que INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1)
                logger.info("[INSTALL] Installation de Prospector depuis GitHub...")
                prospector_url = "git+https://github.com/bd-j/prospector.git"
                logger.info(f"[INSTALL] URL: {prospector_url}")
                
                def _pip_combined_output(res):
                    parts = [res.stdout or "", res.stderr or ""]
                    return "\n".join(p for p in parts if p).strip()
                
                result = subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "pip",
                        "install",
                        prospector_url,
                        "--no-cache-dir",
                    ],
                    capture_output=True,
                    text=True,
                    timeout=600,
                )
                logger.info(f"[INSTALL] Retour code: {result.returncode}")
                logger.info(f"[INSTALL] Sortie combinée:\n{_pip_combined_output(result)}")
                
                if result.returncode != 0:
                    logger.warning(
                        "[INSTALL] Premier essai échoué (souvent fsps : CMake/Fortran absents sous Windows). "
                        "Second essai avec --no-deps…"
                    )
                    result = subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "pip",
                            "install",
                            prospector_url,
                            "--no-cache-dir",
                            "--no-deps",
                        ],
                        capture_output=True,
                        text=True,
                        timeout=600,
                    )
                    logger.info(f"[INSTALL] Retour code (--no-deps): {result.returncode}")
                    logger.info(f"[INSTALL] Sortie combinée:\n{_pip_combined_output(result)}")
                
                if result.returncode != 0:
                    logger.error("[INSTALL] ✗ Erreur lors de l'installation de Prospector")
                    combo = _pip_combined_output(result)
                    logger.error(f"[INSTALL] Sortie: {combo}")
                    messagebox.showerror(
                        "Erreur",
                        "Erreur lors de l'installation de Prospector depuis GitHub.\n\n"
                        f"{combo}\n\n"
                        "Vérifiez Git et la connexion réseau.\n"
                        "Vous pouvez aussi lancer depuis le dossier NPOAP :\n"
                        "INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat\n\n"
                        "Si le message mentionne Fortran / fsps : sous Windows, installez gfortran "
                        "(ex. MSYS2) ou utilisez Installation_fsps/prospector.bat sous WSL2.",
                    )
                    return
                
                logger.info("[INSTALL] ✓ Prospector installé depuis GitHub")
                
                # IMPORTANT: Définir SPS_HOME AVANT de vérifier l'installation
                # Prospector essaie d'utiliser SPS_HOME lors de l'import
                logger.info("[INSTALL] Configuration de SPS_HOME...")
                import os
                import pathlib
                
                # Déterminer SPS_HOME
                if 'SPS_HOME' in os.environ:
                    sps_home = pathlib.Path(os.environ['SPS_HOME'])
                    logger.info(f"[INSTALL] SPS_HOME déjà défini: {sps_home}")
                else:
                    sps_home = pathlib.Path.home() / '.local' / 'share' / 'fsps'
                    os.environ['SPS_HOME'] = str(sps_home)
                    logger.info(f"[INSTALL] SPS_HOME défini: {sps_home}")
                
                # Créer les fichiers stub nécessaires AVANT de vérifier l'installation
                logger.info("[INSTALL] Création des fichiers stub FSPS nécessaires...")
                try:
                    
                    sps_home.mkdir(parents=True, exist_ok=True)
                    (sps_home / 'dust').mkdir(parents=True, exist_ok=True)
                    (sps_home / 'sed').mkdir(parents=True, exist_ok=True)
                    
                    # Créer le fichier stub avec format Prospector
                    dust_file = sps_home / 'dust' / 'Nenkova08_y010_torusg_n10_q2.0.dat'
                    logger.info(f"[INSTALL] Vérification/création du fichier stub: {dust_file}")
                    try:
                        # Vérifier si le fichier existe et a le bon format
                        # Format Prospector: 4 lignes d'en-tête (skip_header=4) puis ~125 lignes avec 10 colonnes
                        needs_recreate = True
                        if dust_file.exists():
                            try:
                                with open(dust_file, 'r') as f:
                                    lines = f.readlines()
                                if len(lines) >= 5:  # Au moins 4 en-têtes + 1 ligne de données
                                    # Vérifier la 5ème ligne (première ligne de données après skip_header=4)
                                    data_line = lines[4].strip()
                                    if data_line and not data_line.startswith('#'):
                                        # IMPORTANT: Vérifier le séparateur (3 espaces exactement)
                                        # Prospector utilise delimiter='   ' (3 espaces)
                                        cols_3spaces = data_line.split('   ')  # Séparateur: 3 espaces
                                        if len(cols_3spaces) == 10:
                                            # Vérifier aussi la ligne suivante
                                            if len(lines) >= 6:
                                                data_line2 = lines[5].strip()
                                                if data_line2 and not data_line2.startswith('#'):
                                                    cols2_3spaces = data_line2.split('   ')  # Séparateur: 3 espaces
                                                    if len(cols2_3spaces) == 10:
                                                        # Vérifier qu'on a au moins ~125 lignes de données
                                                        if len(lines) >= 129:  # 4 en-têtes + 125 lignes
                                                            needs_recreate = False
                                                            logger.info(f"[INSTALL] ✓ Fichier stub existe et a le bon format (10 colonnes avec delimiter=3 espaces, ~{len(lines)} lignes)")
                                                        else:
                                                            logger.warning(f"[INSTALL] Fichier stub a seulement {len(lines)} lignes (besoin d'au moins 129), recréation nécessaire")
                                                    else:
                                                        logger.warning(f"[INSTALL] Fichier stub ligne 6 a {len(cols2_3spaces)} colonnes au lieu de 10 (delimiter incorrect), recréation nécessaire")
                                                else:
                                                    logger.warning(f"[INSTALL] Fichier stub ligne 6 invalide, recréation nécessaire")
                                            else:
                                                logger.warning(f"[INSTALL] Fichier stub a moins de 6 lignes, recréation nécessaire")
                                        else:
                                            # Peut-être que le fichier utilise un autre séparateur (1 espace par exemple)
                                            cols_1space = data_line.split()  # Séparateur: espace(s) quelconque(s)
                                            if len(cols_1space) == 10:
                                                logger.warning(f"[INSTALL] Fichier stub utilise un séparateur incorrect (1 espace au lieu de 3 espaces), recréation nécessaire")
                                            else:
                                                logger.warning(f"[INSTALL] Fichier stub ligne 5 a {len(cols_3spaces)} colonnes avec delimiter=3 espaces (attendu: 10), recréation nécessaire")
                                    else:
                                        logger.warning(f"[INSTALL] Fichier stub ligne 5 invalide ou commentaire, recréation nécessaire")
                                else:
                                    logger.warning(f"[INSTALL] Fichier stub a seulement {len(lines)} lignes (besoin d'au moins 5), recréation nécessaire")
                            except Exception as read_error:
                                logger.warning(f"[INSTALL] Erreur lors de la lecture du fichier stub: {read_error}, recréation nécessaire")
                                needs_recreate = True
                        
                        if needs_recreate:
                            logger.info(f"[INSTALL] Création/recréation du fichier stub avec format Prospector...")
                            # Supprimer l'ancien fichier s'il existe avec un mauvais format
                            if dust_file.exists():
                                try:
                                    dust_file.unlink()
                                    logger.info(f"[INSTALL] Ancien fichier stub supprimé (format incorrect)")
                                except Exception as del_error:
                                    logger.warning(f"[INSTALL] Impossible de supprimer l'ancien fichier: {del_error}")
                            with open(dust_file, 'w') as f:
                                # Format attendu par Prospector (fake_fsps.py ligne 191):
                                # dtype=[('wave', 'f8'), ('fnu_5', '<U20'), ..., ('fnu_150', '<U20')]
                                # delimiter='   ' (3 espaces), skip_header=4, ~125 lignes
                                
                                # Écrire 4 lignes d'en-tête (skip_header=4)
                                f.write("# Nenkova08 AGN torus dust model - Stub file\n")
                                f.write("# This is a stub file created automatically\n")
                                f.write("# Replace with real FSPS data file for full functionality\n")
                                f.write("# wave   fnu_5   fnu_10   fnu_20   fnu_30   fnu_40   fnu_60   fnu_80   fnu_100   fnu_150\n")
                                
                                # Créer ~125 lignes de données
                                for i in range(125):
                                    wave = 1.0 + i * 0.1  # Colonne 1: longueur d'onde
                                    fnu_values = [f"{(j+1)*0.001 + i*0.0001:.6f}" for j in range(9)]  # 9 valeurs fnu
                                    # Séparateur: 3 espaces (delimiter='   ')
                                    line = f"{wave:.6f}   " + "   ".join(fnu_values)
                                    f.write(line + "\n")
                            logger.info(f"[INSTALL] ✓ Fichier stub créé/recréé avec format Prospector (4 en-têtes + 125 lignes, delimiter=3 espaces)")
                    except Exception as file_error:
                        logger.warning(f"[INSTALL] ⚠ Erreur lors de la vérification/création du fichier stub: {file_error}")
                        # Essayer quand même de créer le fichier avec le format correct
                        try:
                            logger.info(f"[INSTALL] Tentative de création du fichier stub en mode récupération...")
                            with open(dust_file, 'w') as f:
                                # Format Prospector: 4 en-têtes + 125 lignes avec 10 colonnes
                                f.write("# Nenkova08 AGN torus dust model - Stub file (recovery mode)\n")
                                f.write("# This is a stub file created automatically\n")
                                f.write("# Replace with real FSPS data file for full functionality\n")
                                f.write("# wave   fnu_5   fnu_10   fnu_20   fnu_30   fnu_40   fnu_60   fnu_80   fnu_100   fnu_150\n")
                                for i in range(125):
                                    wave = 1.0 + i * 0.1
                                    fnu_values = [f"{(j+1)*0.001 + i*0.0001:.6f}" for j in range(9)]
                                    line = f"{wave:.6f}   " + "   ".join(fnu_values)
                                    f.write(line + "\n")
                            logger.info(f"[INSTALL] ✓ Fichier stub créé en mode récupération avec format Prospector")
                        except Exception as e2:
                            logger.error(f"[INSTALL] ✗ Impossible de créer le fichier stub: {e2}")
                            import traceback
                            logger.error(f"[INSTALL] Traceback:\n{traceback.format_exc()}")
                    
                except Exception as stub_error:
                    logger.warning(f"[INSTALL] ⚠ Erreur lors de la création du fichier stub: {stub_error}")
                    import traceback
                    logger.warning(f"[INSTALL] Traceback:\n{traceback.format_exc()}")
                
                # IMPORTANT: S'assurer que SPS_HOME est défini avant l'import de prospect
                # Prospector utilise SPS_HOME lors de l'import
                if 'SPS_HOME' not in os.environ:
                    os.environ['SPS_HOME'] = str(sps_home)
                    logger.info(f"[INSTALL] SPS_HOME défini pour l'import: {sps_home}")
                
                # Vérifier l'installation
                logger.info("[INSTALL] Vérification de l'installation...")
                logger.info(f"[INSTALL] SPS_HOME={os.environ.get('SPS_HOME', 'NON DÉFINI')}")
                prospect_import_ok = False
                try:
                    import prospect
                    logger.info(f"[INSTALL] ✓ Module 'prospect' importé avec succès")
                    logger.info(f"[INSTALL] Version: {getattr(prospect, '__version__', 'inconnue')}")
                    logger.info(f"[INSTALL] Emplacement: {getattr(prospect, '__file__', 'inconnu')}")
                    
                    # Essayer d'importer SpecModel
                    try:
                        from prospect.models import SpecModel
                        logger.info("[INSTALL] ✓ SpecModel importable")
                    except ImportError as e:
                        logger.warning(f"[INSTALL] ⚠ SpecModel non importable: {e}")
                    prospect_import_ok = True
                except Exception as e:
                    logger.error(f"[INSTALL] ✗ Impossible d'importer prospect après installation: {e}")
                    logger.error(f"[INSTALL] Type d'erreur: {type(e).__name__}")
                    import traceback
                    logger.error(f"[INSTALL] Traceback:\n{traceback.format_exc()}")
                
                if prospect_import_ok:
                    messagebox.showinfo(
                        "Succès",
                        "Prospector et ses dépendances ont été installés avec succès.\n\n"
                        "IMPORTANT : redémarrez l'application pour utiliser Prospector.\n\n"
                        "Sans compilateur Fortran, fsps n'est en général pas installé : "
                        "les stubs SPS_HOME suffisent pour une grande partie des usages.\n"
                        "Pour FSPS complet : gfortran sous Windows ou "
                        "Installation_fsps/prospector.bat sous WSL2.",
                    )
                    logger.info("[INSTALL] ========================================")
                    logger.info("[INSTALL] Installation terminée avec succès")
                    logger.info("[INSTALL] ========================================")
                    self.prospector_status_label.config(
                        text="✓ Prospector installé (redémarrez)",
                        foreground="green",
                    )
                else:
                    messagebox.showwarning(
                        "Installation incomplète",
                        "pip a terminé mais l'import de « prospect » a échoué dans cette session.\n\n"
                        "Consultez les logs, vérifiez SPS_HOME et les stubs, puis redémarrez l'application.\n"
                        "Sinon lancez INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat depuis le dossier NPOAP.",
                    )
                    self.prospector_status_label.config(
                        text="Prospector : vérifier logs / redémarrer",
                        foreground="orange",
                    )
                
            except Exception as e:
                logger.error(f"[INSTALL] ✗ Erreur inattendue lors de l'installation: {e}")
                logger.error(f"[INSTALL] Type d'erreur: {type(e).__name__}")
                import traceback
                logger.error(f"[INSTALL] Traceback:\n{traceback.format_exc()}")
                messagebox.showerror("Erreur", f"Erreur lors de l'installation: {e}")
        
        thread = threading.Thread(target=install, daemon=True)
        thread.start()
    
    def browse_file(self):
        """Ouvre un dialogue pour sélectionner un fichier"""
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            title="Sélectionner un fichier de spectre",
            filetypes=[
                ("FITS files", "*.fits *.fit *.fts"),
                ("ASCII files", "*.txt *.dat *.ascii"),
                ("Tous les fichiers", "*.*")
            ]
        )
        
        if path:
            self.file_var.set(path)
    
    def load_spectrum(self):
        """Charge un spectre depuis un fichier"""
        file_path = self.file_var.get()
        
        if not file_path or not Path(file_path).exists():
            messagebox.showerror("Erreur", "Veuillez sélectionner un fichier valide")
            return
        
        if not SPECUTILS_AVAILABLE:
            messagebox.showerror("Erreur", "specutils n'est pas installé")
            return
        
        try:
            file_path_obj = Path(file_path)
            
            # Détecter le type de fichier
            if file_path_obj.suffix.lower() in ['.fits', '.fit', '.fts']:
                # Charger un fichier FITS avec specutils ou astropy
                loaded = False
                
                # Méthode 1: Essayer avec Spectrum1D.read() (specutils moderne)
                try:
                    from specutils import Spectrum1D
                    self.spectrum = Spectrum1D.read(file_path)
                    loaded = True
                    logger.info(f"Spectre chargé avec Spectrum1D.read() depuis {file_path}")
                except Exception as e1:
                    logger.debug(f"Tentative Spectrum1D.read() échouée: {e1}")
                    
                    # Méthode 2: Essayer avec specutils.io.read_fits
                    try:
                        from specutils.io import read_fits
                        self.spectrum = read_fits(file_path)
                        loaded = True
                        logger.info(f"Spectre chargé avec specutils.io.read_fits depuis {file_path}")
                    except Exception as e2:
                        logger.debug(f"Tentative specutils.io.read_fits échouée: {e2}")
                
                # Méthode 3: Fallback - lire manuellement depuis FITS avec astropy
                if not loaded:
                    logger.info("Utilisation du fallback astropy pour charger le FITS")
                    with fits.open(file_path) as hdul:
                        # Essayer différentes extensions
                        for hdu in hdul:
                            if hdu.data is not None and len(hdu.data.shape) >= 1:
                                # Si on a une extension avec des données
                                header = hdu.header
                                
                                # Vérifier si c'est un tableau structuré (structured array)
                                if hdu.data.dtype.names is not None:
                                    # Tableau structuré - extraire les colonnes appropriées
                                    logger.info(f"Tableau structuré détecté avec colonnes: {hdu.data.dtype.names}")
                                    
                                    # Chercher la colonne de longueur d'onde (wave, wavelength, lambda, etc.)
                                    wave_col = None
                                    for col_name in ['wave', 'wavelength', 'lambda', 'wvl', 'wave_log']:
                                        if col_name in hdu.data.dtype.names:
                                            wave_col = col_name
                                            break
                                    
                                    # Chercher la colonne de flux (flux_norm, flux_raw, flux, flux_tac, etc.)
                                    flux_col = None
                                    for col_name in ['flux_norm', 'flux_raw', 'flux', 'flux_tac', 'flux_irc', 'flux_obs']:
                                        if col_name in hdu.data.dtype.names:
                                            flux_col = col_name
                                            break
                                    
                                    if wave_col is None or flux_col is None:
                                        logger.warning(f"Colonnes requises non trouvées. Colonnes disponibles: {hdu.data.dtype.names}")
                                        logger.warning(f"Tentative avec la première colonne pour wave et la seconde pour flux")
                                        # Utiliser les deux premières colonnes numériques
                                        if len(hdu.data.dtype.names) >= 2:
                                            wave_col = hdu.data.dtype.names[0]
                                            flux_col = hdu.data.dtype.names[1]
                                        else:
                                            continue
                                    
                                    logger.info(f"Utilisation de '{wave_col}' pour la longueur d'onde et '{flux_col}' pour le flux")
                                    
                                    # Extraire les colonnes
                                    wavelength_data = np.asarray(hdu.data[wave_col], dtype=np.float64)
                                    flux_data = np.asarray(hdu.data[flux_col], dtype=np.float64)
                                    
                                    # Créer l'axe spectral
                                    wavelength = wavelength_data * u.AA
                                    
                                else:
                                    # Tableau simple - utiliser le code existant
                                    # S'assurer que flux_data est un tableau numpy 1D avec le bon type
                                    if hdu.data.size > 1:
                                        flux_data = np.asarray(hdu.data.flatten(), dtype=np.float64)
                                    else:
                                        flux_data = np.asarray([hdu.data], dtype=np.float64)
                                    
                                    # Essayer de lire les mots-clés WCS pour créer l'axe spectral
                                    if 'CRVAL1' in header and 'CDELT1' in header:
                                        # WCS présent, utiliser CRVAL1 et CDELT1
                                        crval1 = header['CRVAL1']  # Valeur de référence (première longueur d'onde)
                                        cdelt1 = header.get('CDELT1', header.get('CD1_1', 1.0))  # Pas en longueur d'onde
                                        crpix1 = header.get('CRPIX1', 1.0)  # Pixel de référence (généralement 1)
                                        
                                        # Créer l'axe spectral: wavelength = CRVAL1 + (pixel - CRPIX1) * CDELT1
                                        n_pixels = len(flux_data)
                                        pix = np.arange(1, n_pixels + 1)  # Pixels de 1 à N
                                        wavelength = (crval1 + (pix - crpix1) * cdelt1) * u.AA
                                        
                                    elif 'WAT0_001' in header or 'WCSDIM' in header:
                                        # Essayer avec WCS astropy
                                        from astropy.wcs import WCS
                                        try:
                                            w = WCS(header)
                                            pix = np.arange(len(flux_data))
                                            wavelength = w.pixel_to_world(pix) * u.AA
                                        except:
                                            # Si WCS échoue, utiliser un axe artificiel
                                            wavelength = np.arange(len(flux_data)) * u.AA
                                            logger.warning("WCS non disponible, utilisation d'un axe artificiel")
                                    else:
                                        # Pas de WCS, essayer d'autres mots-clés communs
                                        if 'WAVELENGTH' in header:
                                            wavelength_start = header.get('WAVELENGTH', 0)
                                            wavelength_step = header.get('WAVESTEP', 1.0)
                                            wavelength = (wavelength_start + np.arange(len(flux_data)) * wavelength_step) * u.AA
                                        elif 'LAMBDA0' in header:
                                            lambda0 = header.get('LAMBDA0', 0)
                                            delta_lambda = header.get('DELTA_LAMBDA', 1.0)
                                            wavelength = (lambda0 + np.arange(len(flux_data)) * delta_lambda) * u.AA
                                        else:
                                            # Dernier recours: axe artificiel
                                            wavelength = np.arange(len(flux_data)) * u.AA
                                            logger.warning("Aucune information WCS trouvée, utilisation d'un axe artificiel")
                                
                                # Créer le spectre avec specutils
                                if len(flux_data.shape) == 1:
                                    # Spectre 1D
                                    # S'assurer que flux_data est un tableau numpy float64
                                    flux_array = np.asarray(flux_data, dtype=np.float64)
                                    
                                    # Essayer de détecter l'unité du flux
                                    flux_unit = u.Unit('erg cm-2 s-1 AA-1')  # Unité par défaut
                                    if 'BUNIT' in header:
                                        try:
                                            flux_unit = u.Unit(header['BUNIT'])
                                        except:
                                            pass
                                    
                                    # Créer le Quantity astropy avec le tableau numpy
                                    flux = flux_array * flux_unit
                                    
                                    # S'assurer que wavelength est aussi un Quantity valide
                                    if not isinstance(wavelength, u.Quantity):
                                        wavelength = wavelength * u.AA
                                    
                                    # Filtrer les valeurs non finies (NaN, Inf) avant de créer le spectre
                                    finite_mask = np.isfinite(flux_array) & np.isfinite(wavelength_data)
                                    if not np.any(finite_mask):
                                        logger.warning(f"Aucune valeur finie trouvée dans l'extension {hdu.name if hasattr(hdu, 'name') else 'PRIMARY'}")
                                        continue
                                    
                                    # Filtrer les données
                                    flux_filtered = flux_array[finite_mask]
                                    wavelength_filtered = wavelength_data[finite_mask] * u.AA
                                    flux_filtered_qty = flux_filtered * flux_unit
                                    
                                    # Créer le spectre avec les données filtrées
                                    self.spectrum = Spectrum1D(spectral_axis=wavelength_filtered, flux=flux_filtered_qty)
                                    loaded = True
                                    logger.info(f"Spectre chargé depuis extension {hdu.name if hasattr(hdu, 'name') else 'PRIMARY'} ({np.sum(finite_mask)}/{len(flux_array)} points valides)")
                                    break
                    
                    if not loaded:
                        raise ValueError(f"Format FITS non reconnu : aucune extension avec des données spectrales valides trouvée dans {file_path}")
            else:
                # Charger un fichier ASCII (format simple: lambda flux)
                # Note: specutils peut avoir besoin d'un format spécifique
                # Pour l'instant, on lit manuellement
                data = np.loadtxt(file_path)
                if data.shape[1] >= 2:
                    wavelength = data[:, 0] * u.AA
                    flux = data[:, 1] * u.Unit('erg cm-2 s-1 AA-1')
                    self.spectrum = Spectrum1D(spectral_axis=wavelength, flux=flux)
                else:
                    raise ValueError("Format ASCII non reconnu. Attendu: lambda flux")
            
            # Réinitialiser le spectre normalisé
            self.spectrum_normalized = None
            
            # Afficher les informations
            self.update_info()
            
            # Afficher le spectre
            self.plot_spectrum()
            
            messagebox.showinfo("Succès", "Spectre chargé avec succès")
            logger.info(f"Spectre chargé: {file_path}")
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors du chargement du spectre:\n{e}")
            logger.error(f"Erreur chargement spectre: {e}", exc_info=True)
    
    def plot_spectrum(self):
        """Affiche le spectre"""
        if self.spectrum is None:
            return
        
        if not MATPLOTLIB_AVAILABLE:
            return
        
        self.ax.clear()
        
        # Utiliser le spectre normalisé si disponible, sinon le spectre original
        spectrum_to_plot = self.spectrum_normalized if self.spectrum_normalized is not None else self.spectrum
        
        # Convertir en numpy arrays pour le plotting
        wavelength = spectrum_to_plot.spectral_axis.value
        flux = spectrum_to_plot.flux.value
        
        # Afficher
        self.ax.step(wavelength, flux, where='mid', linewidth=1)
        
        self.ax.set_xlabel(f'Wavelength ({spectrum_to_plot.spectral_axis.unit})')
        self.ax.set_ylabel(f'Flux ({spectrum_to_plot.flux.unit})')
        self.ax.set_title('Spectre stellaire')
        self.ax.grid(True, alpha=0.3)
        
        self.canvas.draw()
    
    def normalize_continuum(self):
        """Normalise le continuum du spectre"""
        if self.spectrum is None:
            messagebox.showerror("Erreur", "Chargez d'abord un spectre")
            return
        
        if not SPECUTILS_AVAILABLE:
            return
        
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                # Filtrer les valeurs non finies (NaN, Inf) avant de normaliser
                finite_mask = np.isfinite(self.spectrum.flux.value) & np.isfinite(self.spectrum.spectral_axis.value)
                
                if not np.any(finite_mask):
                    messagebox.showerror("Erreur", "Le spectre ne contient aucune valeur finie")
                    logger.error("Aucune valeur finie trouvée dans le spectre")
                    return
                
                # Si toutes les valeurs sont finies, utiliser directement le spectre
                if np.all(finite_mask):
                    spectrum_to_fit = self.spectrum
                    logger.info("Toutes les valeurs du spectre sont finies, normalisation directe")
                else:
                    # Filtrer le spectre pour ne garder que les valeurs finies
                    n_finite = np.sum(finite_mask)
                    n_total = len(finite_mask)
                    logger.info(f"Filtrage des valeurs non finies: {n_finite}/{n_total} points valides")
                    spectrum_to_fit = Spectrum1D(
                        spectral_axis=self.spectrum.spectral_axis[finite_mask],
                        flux=self.spectrum.flux[finite_mask]
                    )
                
                # Normaliser le continuum sur le spectre filtré
                continuum = fit_generic_continuum(spectrum_to_fit)
                continuum_flux = continuum(spectrum_to_fit.spectral_axis)
                self.spectrum_normalized = spectrum_to_fit / continuum_flux
            
            # Réafficher
            self.plot_spectrum()
            
            # Mettre à jour les informations
            self.update_info()
            
            messagebox.showinfo("Succès", "Continuum normalisé")
            logger.info("Continuum normalisé")
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de la normalisation:\n{e}")
            logger.error(f"Erreur normalisation: {e}", exc_info=True)
    
    def analyze_line(self):
        """Analyse une raie spectrale"""
        if self.spectrum is None:
            messagebox.showerror("Erreur", "Chargez d'abord un spectre")
            return
        
        if not SPECUTILS_AVAILABLE:
            return
        
        try:
            # Lire les limites de la raie
            line_min = float(self.line_min_var.get())
            line_max = float(self.line_max_var.get())
            
            if line_min >= line_max:
                messagebox.showerror("Erreur", "La valeur min doit être inférieure à la valeur max")
                return
            
            # Créer la région spectrale
            region = SpectralRegion(line_min * u.AA, line_max * u.AA)
            
            # Utiliser le spectre normalisé si disponible
            spectrum_to_analyze = self.spectrum_normalized if self.spectrum_normalized is not None else self.spectrum
            
            # Calculer les mesures
            ew = equivalent_width(spectrum_to_analyze, regions=region)
            flux = line_flux(spectrum_to_analyze, regions=region)
            cent = centroid(spectrum_to_analyze, region=region)
            
            # Afficher les résultats
            results_text = f"Analyse de la raie ({line_min:.1f} - {line_max:.1f} Å)\n"
            results_text += f"=" * 40 + "\n"
            results_text += f"Largeur équivalente: {ew:.3f}\n"
            results_text += f"Flux de la raie: {flux:.3e}\n"
            results_text += f"Centroïde: {cent:.3f}\n"
            
            messagebox.showinfo("Résultats", results_text)
            
            # Mettre à jour les informations
            self.update_info()
            self.info_text.insert(tk.END, "\n" + results_text)
            self.info_text.see(tk.END)
            
            logger.info(f"Raie analysée: EW={ew}, Flux={flux}, Centroïde={cent}")
            
        except ValueError as e:
            messagebox.showerror("Erreur", f"Valeurs invalides:\n{e}")
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de l'analyse:\n{e}")
            logger.error(f"Erreur analyse raie: {e}", exc_info=True)
    
    def update_info(self):
        """Met à jour les informations du spectre"""
        if self.spectrum is None:
            self.info_text.delete(1.0, tk.END)
            self.info_text.insert(1.0, "Aucun spectre chargé")
            return
        
        info_text = "Informations du spectre\n"
        info_text += "=" * 30 + "\n\n"
        info_text += f"Nombre de points: {len(self.spectrum.flux)}\n"
        info_text += f"Plage de longueurs d'onde:\n"
        info_text += f"  Min: {self.spectrum.spectral_axis.min():.2f}\n"
        info_text += f"  Max: {self.spectrum.spectral_axis.max():.2f}\n"
        info_text += f"\nFlux:\n"
        info_text += f"  Min: {self.spectrum.flux.min():.3e}\n"
        info_text += f"  Max: {self.spectrum.flux.max():.3e}\n"
        info_text += f"  Moyenne: {self.spectrum.flux.mean():.3e}\n"
        
        if self.spectrum_normalized is not None:
            info_text += f"\n[Continuum normalisé]\n"
        
        self.info_text.delete(1.0, tk.END)
        self.info_text.insert(1.0, info_text)
    
    def export_spectrum(self):
        """Exporte le spectre vers un fichier"""
        if self.spectrum is None:
            messagebox.showerror("Erreur", "Aucun spectre à exporter")
            return
        
        file_path = filedialog.asksaveasfilename(
            initialdir=self.base_dir,
            title="Exporter le spectre",
            defaultextension=".fits",
            filetypes=[
                ("FITS files", "*.fits"),
                ("ASCII files", "*.txt"),
                ("Tous les fichiers", "*.*")
            ]
        )
        
        if file_path:
            try:
                # Utiliser le spectre normalisé si disponible
                spectrum_to_export = self.spectrum_normalized if self.spectrum_normalized is not None else self.spectrum
                
                file_path_obj = Path(file_path)
                if file_path_obj.suffix.lower() in ['.fits', '.fit', '.fts']:
                    # Exporter en FITS (nécessite specutils)
                    from specutils.io import write_fits
                    write_fits.write_fits_spectrum1d(spectrum_to_export, file_path)
                else:
                    # Exporter en ASCII
                    data = np.column_stack([
                        spectrum_to_export.spectral_axis.value,
                        spectrum_to_export.flux.value
                    ])
                    np.savetxt(file_path, data, fmt='%.8e', header='wavelength flux')
                
                messagebox.showinfo("Succès", f"Spectre exporté: {file_path}")
                logger.info(f"Spectre exporté: {file_path}")
                
            except Exception as e:
                messagebox.showerror("Erreur", f"Erreur lors de l'export:\n{e}")
                logger.error(f"Erreur export: {e}", exc_info=True)
    
    def export_plot(self):
        """Exporte le graphique"""
        if not MATPLOTLIB_AVAILABLE or self.spectrum is None:
            messagebox.showerror("Erreur", "Aucun graphique à exporter")
            return
        
        file_path = filedialog.asksaveasfilename(
            initialdir=self.base_dir,
            title="Exporter le graphique",
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("Tous les fichiers", "*.*")
            ]
        )
        
        if file_path:
            try:
                self.fig.savefig(file_path, dpi=150, bbox_inches='tight')
                messagebox.showinfo("Succès", f"Graphique exporté: {file_path}")
                logger.info(f"Graphique exporté: {file_path}")
            except Exception as e:
                messagebox.showerror("Erreur", f"Erreur lors de l'export:\n{e}")
                logger.error(f"Erreur export graphique: {e}", exc_info=True)
    
    def create_sed_from_photometry(self):
        """Crée une SED à partir de données photométriques."""
        if not PROSPECTOR_AVAILABLE or self.prospector_analyzer is None:
            messagebox.showerror("Erreur", "Prospector n'est pas disponible")
            return
        
        # Fenêtre pour entrer les données photométriques
        dialog = tk.Toplevel(self)
        dialog.title("Créer SED depuis Photométrie")
        dialog.geometry("600x400")
        
        # Instructions
        instructions = ttk.Label(
            dialog,
            text="Entrez les données photométriques:\nFormat: longueur d'onde (Å), flux, erreur (optionnel)\nUne ligne par point.",
            wraplength=550
        )
        instructions.pack(pady=10)
        
        # Zone de texte pour les données
        text_frame = ttk.Frame(dialog)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        text_area = tk.Text(text_frame, height=15, width=70)
        scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL, command=text_area.yview)
        text_area.config(yscrollcommand=scrollbar.set)
        
        text_area.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Exemple
        example = "4000.0, 1.5e-16, 1.0e-17\n5000.0, 2.3e-16, 1.2e-17\n6000.0, 3.1e-16, 1.5e-17\n"
        text_area.insert(1.0, example)
        
        def process_photometry():
            try:
                data_text = text_area.get(1.0, tk.END).strip()
                lines = [line.strip() for line in data_text.split('\n') if line.strip()]
                
                wavelengths = []
                fluxes = []
                flux_errors = []
                
                for line in lines:
                    parts = [p.strip() for p in line.split(',')]
                    if len(parts) >= 2:
                        w = float(parts[0])
                        f = float(parts[1])
                        e = float(parts[2]) if len(parts) >= 3 else 0.0
                        
                        wavelengths.append(w)
                        fluxes.append(f)
                        flux_errors.append(e)
                
                if len(wavelengths) == 0:
                    messagebox.showerror("Erreur", "Aucune donnée valide")
                    return
                
                # Créer la SED
                self.sed_data = self.prospector_analyzer.create_sed_from_photometry(
                    np.array(wavelengths),
                    np.array(fluxes),
                    np.array(flux_errors) if any(flux_errors) else None
                )
                
                messagebox.showinfo("Succès", 
                                  f"SED créée avec {len(wavelengths)} points photométriques")
                logger.info(f"SED créée avec {len(wavelengths)} points")
                
                # Afficher dans les infos
                self.info_text.insert(tk.END, f"\n\n=== SED Photométrique ===\n")
                self.info_text.insert(tk.END, f"Points: {len(wavelengths)}\n")
                self.info_text.insert(tk.END, f"Plage: {min(wavelengths):.1f} - {max(wavelengths):.1f} Å\n")
                self.info_text.see(tk.END)
                
                dialog.destroy()
                
            except Exception as e:
                messagebox.showerror("Erreur", f"Erreur lors du traitement:\n{e}")
                logger.error(f"Erreur création SED: {e}", exc_info=True)
        
        # Boutons
        button_frame = ttk.Frame(dialog)
        button_frame.pack(pady=10)
        
        ttk.Button(button_frame, text="Créer SED", command=process_photometry).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Annuler", command=dialog.destroy).pack(side=tk.LEFT, padx=5)
    
    def analyze_galaxy_with_prospector(self):
        """Analyse un spectre de galaxie avec Prospector pour inférer les propriétés stellaires."""
        if not PROSPECTOR_AVAILABLE or self.prospector_analyzer is None:
            messagebox.showerror("Erreur", "Prospector n'est pas disponible")
            return
        
        # Vérifier si on a un spectre ou une SED
        if self.spectrum is None and self.sed_data is None:
            messagebox.showerror("Erreur", 
                               "Chargez d'abord un spectre ou créez une SED photométrique")
            return
        
        try:
            # Préparer les données
            if self.sed_data is not None:
                # Utiliser la SED photométrique
                sed_data = self.sed_data
            elif self.spectrum is not None:
                # Créer une SED depuis le spectre
                sed_data = self.prospector_analyzer.create_sed_from_spectrum(self.spectrum)
            else:
                messagebox.showerror("Erreur", "Aucune donnée disponible")
                return
            
            # Afficher un message de progression
            progress_window = tk.Toplevel(self)
            progress_window.title("Inférence Prospector")
            progress_window.geometry("400x150")
            
            progress_label = ttk.Label(
                progress_window,
                text="Inférence des propriétés stellaires en cours...\n\n"
                     "Cela peut prendre plusieurs minutes.\n"
                     "Veuillez patienter.",
                justify=tk.CENTER
            )
            progress_label.pack(pady=20)
            
            progress_window.update()
            
            # Lancer l'inférence (dans un thread pour ne pas bloquer l'UI)
            import threading
            
            def run_inference():
                try:
                    # Inférer les propriétés stellaires
                    results = self.prospector_analyzer.fit_stellar_properties(
                        sed_data,
                        n_walkers=100,  # Paramètres par défaut
                        n_steps=1000
                    )
                    
                    # Générer le résumé
                    summary = self.prospector_analyzer.get_stellar_properties_summary(results)
                    
                    progress_window.destroy()
                    
                    # Afficher les résultats
                    result_window = tk.Toplevel(self)
                    result_window.title("Résultats Inférence Prospector")
                    result_window.geometry("500x400")
                    
                    result_text = tk.Text(result_window, wrap=tk.WORD, font=("Courier", 10))
                    result_scroll = ttk.Scrollbar(result_window, orient=tk.VERTICAL, command=result_text.yview)
                    result_text.config(yscrollcommand=result_scroll.set)
                    
                    result_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
                    result_scroll.pack(side=tk.RIGHT, fill=tk.Y)
                    
                    result_text.insert(1.0, summary)
                    result_text.config(state=tk.DISABLED)
                    
                    # Mettre à jour les infos
                    self.info_text.insert(tk.END, f"\n\n{summary}\n")
                    self.info_text.see(tk.END)
                    
                    logger.info("Inférence Prospector terminée")
                    
                except Exception as e:
                    progress_window.destroy()
                    messagebox.showerror("Erreur", 
                                       f"Erreur lors de l'inférence Prospector:\n{e}")
                    logger.error(f"Erreur inférence Prospector: {e}", exc_info=True)
            
            thread = threading.Thread(target=run_inference, daemon=True)
            thread.start()
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de l'analyse:\n{e}")
            logger.error(f"Erreur analyse galaxie Prospector: {e}", exc_info=True)

