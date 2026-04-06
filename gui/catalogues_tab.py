# gui/catalogues_tab.py
"""
Onglet unifié pour extraire des données de catalogues astronomiques.
Regroupe les outils d'extraction pour : étoiles, étoiles binaires.
(MAST / exoplanètes : voir l'onglet Analyse des données → Catalogues MAST.)
"""

import logging
import re
import shutil
import threading
import subprocess
import sys
import concurrent.futures
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from pathlib import Path
from datetime import datetime
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np

logger = logging.getLogger(__name__)

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

try:
    from astroquery.mast import Observations, Catalogs
    MAST_AVAILABLE = True
except ImportError:
    MAST_AVAILABLE = False
    Catalogs = None

try:
    import lightkurve as lk
    LIGHTKURVE_AVAILABLE = True
except ImportError:
    LIGHTKURVE_AVAILABLE = False
    lk = None

try:
    from core.catalog_extractor import CatalogExtractor, POPULAR_CATALOGS
    CATALOG_EXTRACTOR_AVAILABLE = True
except ImportError as e:
    CATALOG_EXTRACTOR_AVAILABLE = False
    CatalogExtractor = None
    POPULAR_CATALOGS = {}
    logger.warning(f"CatalogExtractor non disponible: {e} (module supprimé, certaines fonctionnalités seront désactivées)")


class CataloguesTab(ttk.Frame):
    """
    Onglet unifié pour l'extraction de catalogues astronomiques.
    Organisé par type d'objet : étoiles, étoiles binaires.
    """
    
    def __init__(self, parent, exoplanets_only=False):
        super().__init__(parent, padding=10)
        
        self.extractor = None
        self.output_dir = Path.home() / "catalogues"
        self.output_dir.mkdir(exist_ok=True)
        
        # Tous les répertoires de catalogues sont choisis par l'utilisateur
        # Aucun catalogue n'est stocké dans .NPOAP
        
        if exoplanets_only:
            self.pack(fill=tk.BOTH, expand=True)
            try:
                self.create_exoplanets_tab(embed_parent=self)
            except Exception as e:
                logger.error(
                    "Erreur création panneau Exoplanètes (mode intégré) : %s",
                    e,
                    exc_info=True,
                )
                ttk.Label(
                    self,
                    text=f"Erreur interface Exoplanètes / MAST : {e}",
                    foreground="red",
                    wraplength=560,
                    justify="left",
                ).pack(padx=10, pady=10, anchor="w")
            return
        
        self.create_widgets()
    
    def create_widgets(self):
        """Crée l'interface utilisateur avec Notebook interne organisé par type d'objet."""
        
        # En-tête
        header_frame = ttk.Frame(self)
        header_frame.pack(fill="x", pady=(0, 10))
        
        title_label = ttk.Label(
            header_frame,
            text="📚 Extraction de Catalogues Astronomiques",
            font=("Helvetica", 14, "bold")
        )
        title_label.pack()
        
        subtitle_label = ttk.Label(
            header_frame,
            text="Outils d'extraction pour étoiles et étoiles binaires (MAST : Analyse des données)",
            font=("Helvetica", 9),
            foreground="gray"
        )
        subtitle_label.pack()
        
        # Notebook interne pour organiser par type d'objet
        self.inner_notebook = ttk.Notebook(self)
        self.inner_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Zone logs créée dans l'onglet Étoiles
        try:
            self.create_stars_tab()
        except Exception as e:
            logger.error(f"Erreur création onglet Étoiles : {e}", exc_info=True)
        
        self.log_message(
            "Interface d'extraction de catalogues prête.\n"
            "Téléchargement MAST / courbes exoplanètes : menu Analyse des données → A. Catalogues MAST.",
            "info",
        )
    
    def init_catalogs(self):
        """Initialise la liste des catalogues disponibles."""
        if not CATALOG_EXTRACTOR_AVAILABLE:
            # Module catalog_extractor supprimé, fonctionnalité désactivée
            # Les catalogues Gaia DR3 sont toujours disponibles via l'onglet Étoiles
            logger.debug("CatalogExtractor non disponible (module supprimé), fonctionnalité désactivée")
            return
        
        catalogs = list(POPULAR_CATALOGS.keys())
        self.catalog_combo['values'] = catalogs
        if catalogs:
            self.catalog_combo.current(0)
            self.on_catalog_selected()
    
    def on_catalog_selected(self, event=None):
        """Appelé quand un catalogue est sélectionné."""
        catalog_name = self.catalog_var.get()
        if catalog_name in POPULAR_CATALOGS:
            catalog_info = POPULAR_CATALOGS[catalog_name]
            self.catalog_desc_label.config(text=catalog_info.get("description", ""))
            
            # Mettre à jour les types d'objets
            object_types = catalog_info.get("object_types", [])
            self.object_type_combo['values'] = object_types
            if object_types:
                self.object_type_combo.current(0)
    
    def on_search_mode_changed(self, *args):
        """Change l'affichage selon le mode de recherche."""
        mode = self.search_mode_var.get()
        if mode == "region":
            self.region_frame.pack(fill="x", pady=2)
            self.box_frame.pack_forget()
        else:
            self.region_frame.pack_forget()
            self.box_frame.pack(fill="x", pady=2)
    
    def select_output_directory(self, var=None):
        """Sélectionne le répertoire de sortie."""
        initial_dir = str(self.output_dir)
        if var and var.get():
            initial_dir = var.get()
        
        dir_path = filedialog.askdirectory(
            title="Sélectionner le répertoire de sortie",
            initialdir=initial_dir
        )
        if dir_path:
            self.output_dir = Path(dir_path)
            if var:
                var.set(str(self.output_dir))
            elif hasattr(self, 'output_dir_var'):
                self.output_dir_var.set(str(self.output_dir))
    
    def log_message(self, message: str, level: str = "info"):
        """
        Ajoute un message dans la zone de logs ET dans le fichier de log principal.
        
        Parameters
        ----------
        message : str
            Message à afficher
        level : str
            Niveau de log : "info", "warning", "error", "debug"
        """
        # Toujours écrire dans le logger Python (fichier de log)
        log_level = getattr(logging, level.upper(), logging.INFO)
        logger.log(log_level, f"[Catalogues] {message}")
        
        # Afficher dans la fenêtre si disponible
        if hasattr(self, 'log_text') and self.log_text is not None:
            try:
                self.log_text.insert(tk.END, f"[{level.upper()}] {message}\n")
                self.log_text.see(tk.END)
            except Exception as e:
                # Si erreur d'affichage, juste logger
                logger.warning(f"Erreur affichage dans log_text : {e}")
    
    def start_extraction_thread(self):
        """Lance l'extraction dans un thread séparé."""
        if not CATALOG_EXTRACTOR_AVAILABLE:
            messagebox.showerror("Erreur", "Module d'extraction non disponible.")
            return
        
        thread = threading.Thread(target=self.run_extraction, daemon=True)
        thread.start()
    
    def run_extraction(self):
        """Exécute l'extraction des données."""
        try:
            # Vérifier les paramètres
            catalog_name = self.catalog_var.get()
            if not catalog_name:
                self.log_message("❌ Veuillez sélectionner un catalogue.", "error")
                return
            
            # Initialiser l'extracteur
            output_dir = Path(self.output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            self.extractor = CatalogExtractor(output_dir=output_dir)
            self.log_message(f"📚 Extraction depuis : {catalog_name}", "info")
            
            # Récupérer les paramètres
            search_mode = self.search_mode_var.get()
            
            # Récupérer les filtres de magnitude
            mag_min = None
            mag_max = None
            if self.mag_min_var.get().strip():
                try:
                    mag_min = float(self.mag_min_var.get())
                except ValueError:
                    self.log_message("⚠️ Magnitude min invalide, ignorée", "warning")
            
            if self.mag_max_var.get().strip():
                try:
                    mag_max = float(self.mag_max_var.get())
                except ValueError:
                    self.log_message("⚠️ Magnitude max invalide, ignorée", "warning")
            
            mag_column = self.mag_column_var.get().strip() or None
            
            # Exécuter l'extraction selon le mode
            if search_mode == "region":
                # Mode région
                try:
                    ra_str = self.ra_center_var.get().strip()
                    dec_str = self.dec_center_var.get().strip()
                    radius_str = self.radius_var.get().strip()
                    
                    if not ra_str or not dec_str or not radius_str:
                        self.log_message("❌ Veuillez remplir tous les champs (RA, DEC, rayon)", "error")
                        return
                    
                    # Parser les coordonnées
                    try:
                        center_coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
                    except:
                        # Essayer en degrés
                        center_coord = SkyCoord(ra=float(ra_str) * u.deg, dec=float(dec_str) * u.deg)
                    
                    radius = float(radius_str) * u.deg
                    
                    self.log_message(f"📍 Centre : {center_coord}", "info")
                    self.log_message(f"📏 Rayon : {radius}", "info")
                    
                    table = self.extractor.extract_by_region(
                        catalog_name=catalog_name,
                        center_coord=center_coord,
                        radius=radius,
                        mag_limit=mag_max,
                        mag_column=mag_column
                    )
                    
                except Exception as e:
                    self.log_message(f"❌ Erreur lors de l'extraction : {e}", "error")
                    logger.exception(e)
                    return
            
            else:
                # Mode box
                try:
                    ra_min = float(self.ra_min_var.get()) if self.ra_min_var.get().strip() else None
                    ra_max = float(self.ra_max_var.get()) if self.ra_max_var.get().strip() else None
                    dec_min = float(self.dec_min_var.get()) if self.dec_min_var.get().strip() else None
                    dec_max = float(self.dec_max_var.get()) if self.dec_max_var.get().strip() else None
                    
                    if ra_min is None and ra_max is None and dec_min is None and dec_max is None:
                        self.log_message("❌ Veuillez spécifier au moins une limite RA ou DEC", "error")
                        return
                    
                    self.log_message(f"📦 Box : RA [{ra_min or 0}, {ra_max or 360}], DEC [{dec_min or -90}, {dec_max or 90}]", "info")
                    
                    table = self.extractor.extract_by_criteria(
                        catalog_name=catalog_name,
                        ra_min=ra_min,
                        ra_max=ra_max,
                        dec_min=dec_min,
                        dec_max=dec_max,
                        mag_min=mag_min,
                        mag_max=mag_max,
                        mag_column=mag_column
                    )
                    
                except Exception as e:
                    self.log_message(f"❌ Erreur lors de l'extraction : {e}", "error")
                    logger.exception(e)
                    return
            
            # Sauvegarder les résultats
            if len(table) == 0:
                self.log_message("⚠️ Aucun résultat trouvé.", "warning")
                return
            
            self.log_message(f"✅ {len(table)} objets extraits", "info")
            
            # Générer un nom de fichier
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{catalog_name.replace(' ', '_')}_{timestamp}"
            
            output_format = self.output_format_var.get()
            output_path = self.extractor.save_table(table, filename, format=output_format)
            
            self.log_message(f"💾 Fichier sauvegardé : {output_path}", "info")
            self.log_message(f"📊 Colonnes : {', '.join(table.colnames[:10])}{'...' if len(table.colnames) > 10 else ''}", "info")
            
            messagebox.showinfo("Succès", f"Extraction terminée !\n{len(table)} objets extraits\nFichier : {output_path}")
            
        except Exception as e:
            self.log_message(f"❌ Erreur fatale : {e}", "error")
            logger.exception(e)
            messagebox.showerror("Erreur", f"Erreur lors de l'extraction : {e}")
    
    def create_logs_section(self, parent=None):
        """Crée la zone de logs. Si parent est fourni (ex: right_frame), elle est packée dedans ; sinon en bas de self."""
        target = parent if parent is not None else self
        log_frame = ttk.LabelFrame(target, text="Logs et résultats", padding=10)
        if parent is None:
            log_frame.pack(fill=tk.BOTH, expand=False, side=tk.BOTTOM, padx=5, pady=5)
        else:
            log_frame.pack(fill=tk.BOTH, expand=True, pady=(5, 0), padx=5)
        log_text_frame = ttk.Frame(log_frame)
        log_text_frame.pack(fill=tk.BOTH, expand=True)
        self.log_text = tk.Text(
            log_text_frame,
            wrap=tk.WORD,
            font=("Consolas", 9),
            bg="#1e1e1e",
            fg="#d4d4d4",
            insertbackground="#d4d4d4",
            height=8
        )
        scrollbar = ttk.Scrollbar(log_text_frame, orient="vertical", command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.pack(side="left", fill=tk.BOTH, expand=True)
        scrollbar.pack(side="right", fill="y")
    
    def create_stars_tab(self):
        """Crée l'onglet pour l'extraction de catalogues d'étoiles (Vizier)."""
        stars_frame = ttk.Frame(self.inner_notebook, padding=5)
        self.inner_notebook.add(stars_frame, text="⭐ Étoiles")
        
        # En-tête
        header = ttk.Label(
            stars_frame,
            text="Extraction depuis catalogues Vizier (CDS Strasbourg)",
            font=("Helvetica", 10, "bold")
        )
        header.pack(pady=5)
        
        # Frame principal (colonne extraction défilable + colonne logs)
        main_container = ttk.Frame(stars_frame)
        main_container.pack(fill=tk.BOTH, expand=True)
        
        # Colonne gauche (sous-onglet Étoiles) : cadre « Extraction de catalogues » défilable
        left_outer = ttk.Frame(main_container)
        left_outer.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=(0, 10), pady=3)

        left_scroll = ttk.Scrollbar(left_outer, orient=tk.VERTICAL)
        left_canvas = tk.Canvas(
            left_outer,
            highlightthickness=0,
            yscrollcommand=left_scroll.set,
        )
        left_scroll.config(command=left_canvas.yview)
        left_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        left_frame = ttk.LabelFrame(left_canvas, text="Extraction de catalogues", padding=5)
        _left_canvas_window = left_canvas.create_window((0, 0), window=left_frame, anchor=tk.NW)

        def _left_update_scrollregion(_event=None):
            left_canvas.update_idletasks()
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))

        def _left_canvas_on_configure(event):
            left_canvas.itemconfigure(_left_canvas_window, width=event.width)

        left_frame.bind("<Configure>", lambda e: _left_update_scrollregion())
        left_canvas.bind("<Configure>", _left_canvas_on_configure)

        def _left_on_mousewheel(event):
            delta = getattr(event, "delta", 0)
            if delta:
                left_canvas.yview_scroll(int(-1 * (delta / 120)), "units")
            elif getattr(event, "num", None) == 4:
                left_canvas.yview_scroll(-1, "units")
            elif getattr(event, "num", None) == 5:
                left_canvas.yview_scroll(1, "units")

        def _left_bind_wheel_recursive(widget):
            widget.bind("<MouseWheel>", _left_on_mousewheel)
            widget.bind("<Button-4>", _left_on_mousewheel)
            widget.bind("<Button-5>", _left_on_mousewheel)
            for c in widget.winfo_children():
                _left_bind_wheel_recursive(c)
        
        # Colonne droite : logs
        right_frame = ttk.Frame(main_container)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, pady=5)
        right_frame.config(width=200)  # Largeur en pixels
        
        # ----- A. Extraction Vizier (CDS) -----
        vizier_frame = ttk.LabelFrame(left_frame, text="A. Extraction Vizier (catalogues CDS)", padding=2)
        vizier_frame.pack(fill="x", pady=(0, 5))
        
        # Sélection du catalogue (cadre réduit)
        catalog_frame = ttk.LabelFrame(vizier_frame, text="1. Catalogue", padding=2)
        catalog_frame.pack(fill="x", pady=1)
        catalog_frame.config(width=120)  # Largeur réduite pour le cadre
        catalog_frame.pack_propagate(False)  # Empêche le cadre de s'étendre
        
        ttk.Label(catalog_frame, text="Catalogue :", font=("Arial", 8)).pack(anchor="w")
        self.stars_catalog_var = tk.StringVar()
        self.stars_catalog_combo = ttk.Combobox(
            catalog_frame,
            textvariable=self.stars_catalog_var,
            state="readonly",
            width=15  # Largeur réduite
        )
        self.stars_catalog_combo.pack(fill="x", pady=1)
        self.stars_catalog_combo.bind("<<ComboboxSelected>>", self.on_stars_catalog_selected)
        
        self.stars_catalog_desc_label = ttk.Label(
            catalog_frame,
            text="",
            font=("Helvetica", 7),
            foreground="gray",
            wraplength=100  # Réduit encore plus
        )
        self.stars_catalog_desc_label.pack(fill="x", pady=1)
        
        criteria_frame = ttk.LabelFrame(
            vizier_frame,
            text="2. Centre, rayon et magnitude",
            padding=6,
        )
        criteria_frame.pack(fill="x", pady=1)

        ttk.Label(
            criteria_frame,
            text="Recherche Vizier dans un disque autour du centre (h/m/s ou degrés ; rayon en degrés).",
            font=("Helvetica", 8),
            foreground="gray",
            wraplength=300,
            justify=tk.LEFT,
        ).pack(anchor="w", pady=(0, 4))

        ttk.Label(criteria_frame, text="Centre — ascension droite :").pack(anchor="w")
        self.stars_ra_center_var = tk.StringVar()
        ttk.Entry(criteria_frame, textvariable=self.stars_ra_center_var, width=18).pack(fill="x", pady=1)
        ttk.Label(criteria_frame, text="Centre — déclinaison :").pack(anchor="w")
        self.stars_dec_center_var = tk.StringVar()
        ttk.Entry(criteria_frame, textvariable=self.stars_dec_center_var, width=18).pack(fill="x", pady=1)
        ttk.Label(criteria_frame, text="Rayon (degrés) :").pack(anchor="w")
        self.stars_radius_var = tk.StringVar(value="1.0")
        ttk.Entry(criteria_frame, textvariable=self.stars_radius_var, width=18).pack(fill="x", pady=1)

        ttk.Separator(criteria_frame, orient=tk.HORIZONTAL).pack(fill="x", pady=6)
        ttk.Label(criteria_frame, text="Magnitude (optionnel)", font=("Helvetica", 9)).pack(anchor="w")
        mag_row = ttk.Frame(criteria_frame)
        mag_row.pack(fill="x", pady=2)
        self.stars_mag_min_var = tk.StringVar()
        self.stars_mag_max_var = tk.StringVar()
        ttk.Label(mag_row, text="Min (≥)", font=("Arial", 8)).pack(side=tk.LEFT)
        ttk.Entry(mag_row, textvariable=self.stars_mag_min_var, width=7).pack(side=tk.LEFT, padx=2)
        ttk.Label(mag_row, text="Max (≤)", font=("Arial", 8)).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Entry(mag_row, textvariable=self.stars_mag_max_var, width=7).pack(side=tk.LEFT, padx=2)
        ttk.Label(mag_row, text="vide = pas de filtre", font=("Helvetica", 7), foreground="gray").pack(
            side=tk.LEFT, padx=(6, 0)
        )
        
        # Répertoire de sortie pour l'extraction Vizier
        ttk.Label(vizier_frame, text="Répertoire de sortie (Vizier) :", font=("Helvetica", 9)).pack(anchor="w", pady=(2, 1))
        vizier_output_frame = ttk.Frame(vizier_frame)
        vizier_output_frame.pack(fill="x", pady=1)
        self.stars_output_dir_var = tk.StringVar(value=str(self.output_dir))
        ttk.Entry(vizier_output_frame, textvariable=self.stars_output_dir_var, width=25).pack(side="left", fill="x", expand=True)
        ttk.Button(vizier_output_frame, text="📁", command=lambda: self.select_output_directory(self.stars_output_dir_var), width=3).pack(side="left", padx=2)
        
        ttk.Button(
            vizier_frame,
            text="📥 Extraire les données (Vizier)",
            command=self.start_stars_extraction_thread,
            width=40
        ).pack(pady=3)
        
        # ----- B. Paramètres d'extraction catalogue Gaia DR3 (archive ESA) -----
        gaia_frame = ttk.LabelFrame(left_frame, text="B. Paramètres extraction catalogue Gaia DR3 (archive ESA)", padding=2)
        gaia_frame.pack(fill="x", pady=(5, 5))
        
        # Magnitude limite G
        ttk.Label(gaia_frame, text="Magnitude limite G :").pack(anchor="w", pady=(2, 1))
        self.gaia_mag_var = tk.StringVar(value="18.0")
        ttk.Entry(gaia_frame, textvariable=self.gaia_mag_var, width=15).pack(anchor="w", pady=1)
        
        # Largeur de degré horaire (RA step)
        ttk.Label(gaia_frame, text="Largeur de degré horaire (RA step) :").pack(anchor="w", pady=(2, 1))
        ra_step_frame = ttk.Frame(gaia_frame)
        ra_step_frame.pack(fill="x", pady=1)
        self.gaia_ra_step_var = tk.StringVar(value="20.0")
        ra_step_combo = ttk.Combobox(ra_step_frame, textvariable=self.gaia_ra_step_var, values=("5.0", "10.0", "15.0", "20.0", "30.0", "60.0"), width=12, state="readonly")
        ra_step_combo.pack(side="left", padx=2)
        ttk.Label(ra_step_frame, text="°", font=("Arial", 8), foreground="gray").pack(side="left", padx=5)
        
        # Hauteur de déclinaison (DEC range)
        ttk.Label(gaia_frame, text="Hauteur de déclinaison :").pack(anchor="w", pady=(2, 1))
        dec_range_frame = ttk.Frame(gaia_frame)
        dec_range_frame.pack(fill="x", pady=1)
        ttk.Label(dec_range_frame, text="DEC min :", width=8).pack(side="left")
        self.gaia_dec_min_var = tk.StringVar()
        ttk.Entry(dec_range_frame, textvariable=self.gaia_dec_min_var, width=12).pack(side="left", padx=2)
        ttk.Label(dec_range_frame, text="°", width=2).pack(side="left")
        ttk.Label(dec_range_frame, text="DEC max :", width=8).pack(side="left", padx=(10, 0))
        self.gaia_dec_max_var = tk.StringVar()
        ttk.Entry(dec_range_frame, textvariable=self.gaia_dec_max_var, width=12).pack(side="left", padx=2)
        ttk.Label(dec_range_frame, text="°", width=2).pack(side="left")
        
        # Plage RA (optionnel)
        ttk.Label(gaia_frame, text="Plage RA (optionnel) :").pack(anchor="w", pady=(2, 1))
        ra_range_frame = ttk.Frame(gaia_frame)
        ra_range_frame.pack(fill="x", pady=1)
        ttk.Label(ra_range_frame, text="RA min :", width=8).pack(side="left")
        self.gaia_ra_min_h_var = tk.StringVar()
        ttk.Entry(ra_range_frame, textvariable=self.gaia_ra_min_h_var, width=12).pack(side="left", padx=2)
        ttk.Label(ra_range_frame, text="h", width=2).pack(side="left")
        ttk.Label(ra_range_frame, text="RA max :", width=8).pack(side="left", padx=(10, 0))
        self.gaia_ra_max_h_var = tk.StringVar()
        ttk.Entry(ra_range_frame, textvariable=self.gaia_ra_max_h_var, width=12).pack(side="left", padx=2)
        ttk.Label(ra_range_frame, text="h", width=2).pack(side="left")
        
        # Hémisphère
        ttk.Label(gaia_frame, text="Hémisphère :").pack(anchor="w", pady=(2, 1))
        self.gaia_hemisphere_var = tk.StringVar(value="both")
        hem_frame = ttk.Frame(gaia_frame)
        hem_frame.pack(fill="x", pady=1)
        ttk.Radiobutton(hem_frame, text="Nord", variable=self.gaia_hemisphere_var, value="north").pack(side="left", padx=5)
        ttk.Radiobutton(hem_frame, text="Sud", variable=self.gaia_hemisphere_var, value="south").pack(side="left", padx=5)
        ttk.Radiobutton(hem_frame, text="Les deux", variable=self.gaia_hemisphere_var, value="both").pack(side="left", padx=5)
        
        # Option ignorer fichiers existants
        self.gaia_skip_existing_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            gaia_frame,
            text="Ignorer les fichiers existants",
            variable=self.gaia_skip_existing_var
        ).pack(anchor="w", pady=(2, 1))
        
        # Répertoire de sortie Gaia DR3
        ttk.Label(gaia_frame, text="Répertoire de sortie Gaia DR3 :").pack(anchor="w", pady=(2, 1))
        gaia_output_frame = ttk.Frame(gaia_frame)
        gaia_output_frame.pack(fill="x", pady=1)
        self.gaia_output_dir_var = tk.StringVar(value=str(self.output_dir / "gaia_dr3"))
        ttk.Entry(gaia_output_frame, textvariable=self.gaia_output_dir_var, width=25).pack(side="left", fill="x", expand=True)
        ttk.Button(gaia_output_frame, text="📁", command=lambda: self.select_output_directory(self.gaia_output_dir_var), width=3).pack(side="left", padx=2)
        
        # Bouton téléchargement Gaia DR3
        ttk.Button(
            gaia_frame,
            text="📥 Télécharger catalogues Gaia DR3",
            command=self.start_gaia_download_thread,
            width=40
        ).pack(pady=3)
        
        # Filtres magnitude
        #"mag_frame = ttk.LabelFrame(left_frame, text="3. Filtres magnitude", padding=5)
        #mag_frame.pack(fill="x", pady=5)
       # 
       # ttk.Label(mag_frame, text="Magnitude min :").pack(anchor="w")
       # self.stars_mag_min_var = tk.StringVar()
       # ttk.Entry(mag_frame, textvariable=self.stars_mag_min_var, width=30).pack(fill="x", pady=2)
        
       #ttk.Label(mag_frame, text="Magnitude max :").pack(anchor="w")
        #self.stars_mag_max_var = tk.StringVar()
        #ttk.Entry(mag_frame, textvariable=self.stars_mag_max_var, width=30).pack(fill="x", pady=2)
       # 
        # Options de sortie
        #output_frame = ttk.LabelFrame(left_frame, text="4. Options de sortie", padding=5)
        #output_frame.pack(fill="x", pady=5)
       # 
       # ttk.Label(output_frame, text="Répertoire de sortie :").pack(anchor="w")
       # output_path_frame = ttk.Frame(output_frame)
       # output_path_frame.pack(fill="x", pady=2)
       # self.stars_output_dir_var = tk.StringVar(value=str(self.output_dir))
       # ttk.Entry(output_path_frame, textvariable=self.stars_output_dir_var, width=25).pack(side="left", fill="x", expand=True)
       # ttk.Button(output_path_frame, text="📁", command=lambda: self.select_output_directory(self.stars_output_dir_var), width=3).pack(side="left", padx=2)
        
        # Options de sortie (format et filtres pour Gaia DR3)
        output_frame = ttk.LabelFrame(gaia_frame, text="Format et filtres (Gaia DR3)", padding=2)
        output_frame.pack(fill="x", pady=(3, 1))
        
        # Format
        ttk.Label(output_frame, text="Format :").pack(anchor="w")
        self.stars_output_format_var = tk.StringVar(value="csv.gz")
        format_frame = ttk.Frame(output_frame)
        format_frame.pack(fill="x", pady=1)
        ttk.Radiobutton(format_frame, text="CSV", variable=self.stars_output_format_var, value="csv").pack(side="left", padx=5)
        ttk.Radiobutton(format_frame, text="CSV.GZIP", variable=self.stars_output_format_var, value="csv.gz").pack(side="left", padx=5)
        
        # Filtres par type d'objet
        ttk.Label(output_frame, text="Filtres :").pack(anchor="w", pady=(3, 1))
        self.gaia_filter_variables_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(output_frame, text="Étoiles variables uniquement", variable=self.gaia_filter_variables_var).pack(anchor="w", pady=1)
        self.gaia_filter_galaxies_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(output_frame, text="Galaxies uniquement", variable=self.gaia_filter_galaxies_var).pack(anchor="w", pady=1)
        
        # Initialiser les catalogues
        self.init_stars_catalogs()

        _left_bind_wheel_recursive(left_frame)
        left_canvas.bind("<MouseWheel>", _left_on_mousewheel)
        left_canvas.bind("<Button-4>", _left_on_mousewheel)
        left_canvas.bind("<Button-5>", _left_on_mousewheel)
        left_canvas.after_idle(_left_update_scrollregion)

        self.create_logs_section(parent=right_frame)
    
    
    def _npoap_project_root(self) -> Path:
        """Racine du dépôt NPOAP (parent de gui/)."""
        return Path(__file__).resolve().parent.parent

    def _tess_lc_fits_script(self) -> Path:
        return self._npoap_project_root() / "scripts" / "tess_lc_fits_to_txt.py"

    def create_exoplanets_tab(self, embed_parent):
        """Panneau MAST + TESS (FITS édité → LcTools .txt), intégré dans ``embed_parent``.

        Utilisé depuis l'onglet Analyse des données (Catalogues MAST), pas depuis Catalogues.
        """
        exoplanets_frame = ttk.Frame(embed_parent, padding=10)
        exoplanets_frame.pack(fill=tk.BOTH, expand=True)

        ttk.Label(
            exoplanets_frame,
            text="MAST, puis courbe TESS : FITS édité → normalisation LcTools (.txt)",
            font=("Helvetica", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        exo_nb = ttk.Notebook(exoplanets_frame)
        exo_nb.pack(fill=tk.BOTH, expand=True, pady=4)

        mast_tab = ttk.Frame(exo_nb, padding=8)
        exo_nb.add(mast_tab, text="📡 Téléchargement MAST")
        # Toujours défini (liaison onglet TESS → dossier MAST même si astroquery absent)
        self.mast_output_dir_var = tk.StringVar(value=str(self.output_dir / "mast_lightcurves"))

        if not MAST_AVAILABLE:
            ttk.Label(
                mast_tab,
                text="astroquery.mast non disponible. Installez astroquery pour activer le téléchargement MAST.",
                foreground="red",
            ).pack(anchor="w", pady=8)
            mast_na_dir = ttk.LabelFrame(mast_tab, text="Dossier sortie (liaison onglet TESS → bouton MAST)", padding=6)
            mast_na_dir.pack(fill="x", pady=(0, 6))
            ttk.Label(
                mast_na_dir,
                text="Indiquez où se trouvent vos *_lc.fits téléchargés à la main (ex. archive STScI).",
                foreground="gray",
                wraplength=640,
                justify="left",
            ).pack(anchor="w", pady=(0, 4))
            row_na = ttk.Frame(mast_na_dir)
            row_na.pack(fill="x", pady=2)
            ttk.Entry(row_na, textvariable=self.mast_output_dir_var, width=50).pack(
                side="left", fill="x", expand=True, padx=(0, 4)
            )
            ttk.Button(
                row_na,
                text="📁",
                width=3,
                command=lambda: self.select_output_directory(self.mast_output_dir_var),
            ).pack(side="left")
            dialog_frame_na = ttk.LabelFrame(mast_tab, text="Dialogue MAST (inactif)", padding=6)
            dialog_frame_na.pack(fill=tk.BOTH, expand=True, pady=(6, 2))
            dtf = ttk.Frame(dialog_frame_na)
            dtf.pack(fill=tk.BOTH, expand=True)
            self.mast_dialog_text = tk.Text(
                dtf,
                wrap=tk.WORD,
                font=("Consolas", 9),
                bg="#111111",
                fg="#e0e0e0",
                height=8,
            )
            dscroll_na = ttk.Scrollbar(dtf, orient="vertical", command=self.mast_dialog_text.yview)
            self.mast_dialog_text.configure(yscrollcommand=dscroll_na.set)
            self.mast_dialog_text.pack(side="left", fill=tk.BOTH, expand=True)
            dscroll_na.pack(side="right", fill="y")
            self.mast_dialog_text.insert(tk.END, "Téléchargement MAST désactivé (installez astroquery).\n")
        else:
            ttk.Label(
                mast_tab,
                text="Courbes MAST (lightkurve en priorité, puis astroquery)",
                foreground="gray",
            ).pack(anchor="w", pady=(0, 6))

            params_frame = ttk.LabelFrame(mast_tab, text="Paramètres MAST", padding=8)
            params_frame.pack(fill="x", pady=5)

            target_row = ttk.Frame(params_frame)
            target_row.pack(fill="x", pady=2)
            ttk.Label(target_row, text="Cible :", width=18).pack(side="left")
            self.mast_target_var = tk.StringVar()
            ttk.Entry(target_row, textvariable=self.mast_target_var, width=35).pack(side="left", padx=4)
            ttk.Label(target_row, text="(ex: TIC 25155310, TOI-700)", foreground="gray").pack(side="left")

            radius_row = ttk.Frame(params_frame)
            radius_row.pack(fill="x", pady=2)
            ttk.Label(radius_row, text="Rayon (deg) :", width=18).pack(side="left")
            self.mast_radius_var = tk.StringVar(value="0.02")
            ttk.Entry(radius_row, textvariable=self.mast_radius_var, width=10).pack(side="left", padx=4)

            mission_row = ttk.Frame(params_frame)
            mission_row.pack(fill="x", pady=2)
            ttk.Label(mission_row, text="Mission (optionnel) :", width=18).pack(side="left")
            self.mast_mission_var = tk.StringVar(value="TESS")
            mission_combo = ttk.Combobox(
                mission_row,
                textvariable=self.mast_mission_var,
                state="readonly",
                values=["", "TESS", "Kepler", "K2", "HST", "JWST"],
                width=12,
            )
            mission_combo.pack(side="left", padx=4)
            ttk.Label(mission_row, text="Laisser vide = toutes missions", foreground="gray").pack(side="left")

            # Pipelines TESS (HLSP) : requêtes lightkurve séparées par auteur MAST
            tess_pipe_row = ttk.Frame(params_frame)
            tess_pipe_row.pack(fill="x", pady=2)
            ttk.Label(tess_pipe_row, text="Pipelines TESS :", width=18).pack(side="left")
            self.mast_tess_spoc_var = tk.BooleanVar(value=True)
            self.mast_tess_qlp_var = tk.BooleanVar(value=True)
            self.mast_tess_eleanor_var = tk.BooleanVar(value=True)
            ttk.Checkbutton(
                tess_pipe_row, text="SPOC", variable=self.mast_tess_spoc_var
            ).pack(side="left", padx=(0, 6))
            ttk.Checkbutton(
                tess_pipe_row, text="QLP", variable=self.mast_tess_qlp_var
            ).pack(side="left", padx=(0, 6))
            ttk.Checkbutton(
                tess_pipe_row,
                text="GSFC-ELEANOR-LITE",
                variable=self.mast_tess_eleanor_var,
            ).pack(side="left", padx=(0, 6))
            ttk.Label(
                tess_pipe_row,
                text="(mission TESS : une ou plusieurs cases = requêtes par auteur ; aucune = tous les produits)",
                foreground="gray",
                font=("Helvetica", 8),
            ).pack(side="left", padx=(4, 0))

            engine_row = ttk.Frame(params_frame)
            engine_row.pack(fill="x", pady=2)
            ttk.Label(engine_row, text="Moteur :", width=18).pack(side="left")
            self.mast_engine_var = tk.StringVar(value="auto")
            ttk.Combobox(
                engine_row,
                textvariable=self.mast_engine_var,
                state="readonly",
                values=["auto", "lightkurve", "astroquery"],
                width=14,
            ).pack(side="left", padx=4)
            ttk.Label(
                engine_row,
                text="Si timeout astroquery : pip install lightkurve",
                foreground="gray",
            ).pack(side="left")

            max_row = ttk.Frame(params_frame)
            max_row.pack(fill="x", pady=2)
            ttk.Label(max_row, text="Max produits :", width=18).pack(side="left")
            self.mast_max_products_var = tk.StringVar(value="200")
            ttk.Entry(max_row, textvariable=self.mast_max_products_var, width=10).pack(side="left", padx=4)

            out_row = ttk.Frame(params_frame)
            out_row.pack(fill="x", pady=2)
            ttk.Label(out_row, text="Répertoire sortie :", width=18).pack(side="left")
            ttk.Entry(out_row, textvariable=self.mast_output_dir_var, width=45).pack(
                side="left", fill="x", expand=True, padx=4
            )
            ttk.Button(
                out_row,
                text="📁",
                command=lambda: self.select_output_directory(self.mast_output_dir_var),
                width=3,
            ).pack(side="left")
            ttk.Label(
                params_frame,
                text="Les FITS sont regroupés à plat dans un sous-dossier du répertoire ci-dessus "
                "(ex. …\\TESS pour mission TESS ; pas d’arborescence mastDownload).",
                foreground="gray",
                font=("Helvetica", 8),
                wraplength=720,
                justify="left",
            ).pack(anchor="w", pady=(2, 0))

            ttk.Button(
                mast_tab,
                text="📥 Télécharger les courbes de lumière (MAST)",
                command=self.start_mast_lightcurves_download_thread,
                width=45,
            ).pack(pady=8)

            hint = (
                "LcTools/LcGenerator télécharge souvent en direct depuis les archives ; "
                "astroquery Observations.query_object utilise l’API Portal MAST et peut rester "
                "sans réponse derrière certains pare-feu. lightkurve contourne en général ce blocage."
            )
            ttk.Label(mast_tab, text=hint, foreground="gray", wraplength=720, justify="left").pack(
                anchor="w", pady=(0, 4)
            )

            dialog_frame = ttk.LabelFrame(mast_tab, text="Dialogue MAST", padding=6)
            dialog_frame.pack(fill=tk.BOTH, expand=True, pady=(6, 2))

            dialog_text_frame = ttk.Frame(dialog_frame)
            dialog_text_frame.pack(fill=tk.BOTH, expand=True)

            self.mast_dialog_text = tk.Text(
                dialog_text_frame,
                wrap=tk.WORD,
                font=("Consolas", 9),
                bg="#111111",
                fg="#e0e0e0",
                insertbackground="#e0e0e0",
                height=12,
            )
            dialog_scroll = ttk.Scrollbar(
                dialog_text_frame, orient="vertical", command=self.mast_dialog_text.yview
            )
            self.mast_dialog_text.configure(yscrollcommand=dialog_scroll.set)
            self.mast_dialog_text.pack(side="left", fill=tk.BOTH, expand=True)
            dialog_scroll.pack(side="right", fill="y")

            self._mast_dialog("Prêt. Configurez les paramètres puis lancez la requête MAST.")
            self._mast_dialog(f"Interpréteur Python : {sys.executable}")
            self._mast_dialog(
                f"lightkurve importable : {LIGHTKURVE_AVAILABLE} "
                "(utilisez LANCER_NPOAP_ASTROENV.bat si False alors que pip install est dans astroenv)"
            )

        fits_edit_tab = ttk.Frame(exo_nb, padding=8)
        exo_nb.add(fits_edit_tab, text="🧩 Edition Fits")
        self._create_fits_edit_tab(fits_edit_tab)

        tess_wf_tab = ttk.Frame(exo_nb, padding=8)
        exo_nb.add(tess_wf_tab, text="📊 LcTools")
        self._create_tess_workflow_tab(tess_wf_tab)

    def _fits_edit_log(self, message: str) -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {message}\n"
        logger.info(f"[Catalogues][FITS_EDIT] {message}")
        try:
            if hasattr(self, "fits_edit_log_text") and self.fits_edit_log_text is not None:
                self.fits_edit_log_text.insert(tk.END, line)
                self.fits_edit_log_text.see(tk.END)
        except tk.TclError:
            pass

    @staticmethod
    def _find_lightcurve_hdu(hdul: fits.HDUList):
        """Retourne le BinTableHDU contenant TIME, sinon None."""
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and getattr(hdu, "columns", None) is not None:
                cols = {c.upper() for c in hdu.columns.names or []}
                if "TIME" in cols:
                    return hdu
        return None

    @staticmethod
    def _find_matching_col(columns, candidates):
        cols_map = {str(c).upper(): c for c in columns}
        for cand in candidates:
            c = cols_map.get(cand.upper())
            if c is not None:
                return c
        return None

    def _extract_pdcsap_table_from_file(self, fits_path: Path) -> Table:
        """Extrait TIME, PDCSAP_FLUX, PDCSAP_FLUX_ERROR d'un FITS en table astropy."""
        with fits.open(fits_path, memmap=False) as hdul:
            lc_hdu = self._find_lightcurve_hdu(hdul)
            if lc_hdu is None:
                raise ValueError(f"Aucune table avec colonne TIME: {fits_path.name}")
            data = lc_hdu.data
            if data is None:
                raise ValueError(f"Table vide: {fits_path.name}")
            names = list(data.names or [])
            time_col = self._find_matching_col(names, ["TIME"])
            flux_col = self._find_matching_col(names, ["PDCSAP_FLUX"])
            err_col = self._find_matching_col(names, ["PDCSAP_FLUX_ERROR", "PDCSAP_FLUX_ERR"])
            if not time_col or not flux_col or not err_col:
                raise ValueError(
                    f"Colonnes manquantes dans {fits_path.name}: "
                    f"TIME={bool(time_col)}, PDCSAP_FLUX={bool(flux_col)}, PDCSAP_FLUX_ERROR/ERR={bool(err_col)}"
                )
            t = np.asarray(data[time_col], dtype=float)
            f = np.asarray(data[flux_col], dtype=float)
            e = np.asarray(data[err_col], dtype=float)
            return Table(
                {
                    "TIME": t,
                    "PDCSAP_FLUX": f,
                    "PDCSAP_FLUX_ERROR": e,
                }
            )

    def _create_fits_edit_tab(self, parent: ttk.Frame) -> None:
        ttk.Label(
            parent,
            text="Edition FITS : nettoyage des colonnes + fusion multi-fichiers",
            font=("Helvetica", 10, "bold"),
        ).pack(anchor="w", pady=(0, 6))

        # --- Bloc 1 : nettoyage FITS ---
        clean_box = ttk.LabelFrame(
            parent,
            text="1 — Nettoyer un FITS (colonnes conservées : TIME, PDCSAP_FLUX, PDCSAP_FLUX_ERROR)",
            padding=8,
        )
        clean_box.pack(fill="x", pady=4)

        in_row = ttk.Frame(clean_box)
        in_row.pack(fill="x", pady=2)
        ttk.Label(in_row, text="FITS entrée :", width=18).pack(side="left")
        self.fits_clean_input_var = tk.StringVar()
        ttk.Entry(in_row, textvariable=self.fits_clean_input_var, width=60).pack(
            side="left", fill="x", expand=True, padx=4
        )
        ttk.Button(
            in_row,
            text="📂",
            width=3,
            command=lambda: self._tess_lc_browse_fits(self.fits_clean_input_var),
        ).pack(side="left", padx=(0, 2))
        ttk.Button(
            in_row,
            text="MAST",
            width=5,
            command=self._browse_fits_clean_from_mast,
        ).pack(side="left")

        out_row = ttk.Frame(clean_box)
        out_row.pack(fill="x", pady=2)
        ttk.Label(out_row, text="FITS sortie :", width=18).pack(side="left")
        self.fits_clean_output_var = tk.StringVar()
        ttk.Entry(out_row, textvariable=self.fits_clean_output_var, width=60).pack(
            side="left", fill="x", expand=True, padx=4
        )
        ttk.Button(
            out_row,
            text="💾",
            width=3,
            command=self._browse_fits_clean_output,
        ).pack(side="left")

        ttk.Button(
            clean_box,
            text="✓ Nettoyer le FITS",
            command=self._run_clean_single_fits,
            width=28,
        ).pack(anchor="w", pady=(6, 2))

        # --- Bloc 2 : fusion FITS ---
        merge_box = ttk.LabelFrame(
            parent,
            text="2 — Fusionner plusieurs FITS (TIME trié, mêmes 3 colonnes)",
            padding=8,
        )
        merge_box.pack(fill="both", expand=True, pady=6)

        list_btns = ttk.Frame(merge_box)
        list_btns.pack(fill="x", pady=(0, 4))
        ttk.Button(list_btns, text="➕ Ajouter FITS", command=self._add_merge_fits_files).pack(side="left", padx=2)
        ttk.Button(list_btns, text="➖ Retirer sélection", command=self._remove_merge_fits_selected).pack(side="left", padx=2)
        ttk.Button(list_btns, text="🧹 Vider liste", command=self._clear_merge_fits_list).pack(side="left", padx=2)

        self.fits_merge_listbox = tk.Listbox(merge_box, height=8, selectmode=tk.EXTENDED, font=("Consolas", 9))
        self.fits_merge_listbox.pack(fill="both", expand=True, pady=2)
        self.fits_merge_paths: list[str] = []

        merge_out_row = ttk.Frame(merge_box)
        merge_out_row.pack(fill="x", pady=2)
        ttk.Label(merge_out_row, text="FITS fusionné :", width=18).pack(side="left")
        self.fits_merge_output_var = tk.StringVar()
        ttk.Entry(merge_out_row, textvariable=self.fits_merge_output_var, width=60).pack(
            side="left", fill="x", expand=True, padx=4
        )
        ttk.Button(
            merge_out_row,
            text="💾",
            width=3,
            command=self._browse_fits_merge_output,
        ).pack(side="left")

        ttk.Button(
            merge_box,
            text="🔗 Fusionner les FITS",
            command=self._run_merge_fits,
            width=28,
        ).pack(anchor="w", pady=(6, 2))

        log_box = ttk.LabelFrame(parent, text="Journal Edition Fits", padding=4)
        log_box.pack(fill="x", pady=(4, 0))
        self.fits_edit_log_text = tk.Text(
            log_box,
            wrap=tk.WORD,
            font=("Consolas", 9),
            height=8,
            bg="#1a1a1a",
            fg="#e8e8e8",
            insertbackground="#e8e8e8",
        )
        self.fits_edit_log_text.pack(fill="x", expand=False)
        self._fits_edit_log("Prêt. Nettoyage FITS et fusion multi-FITS disponibles.")

    def _browse_fits_clean_from_mast(self) -> None:
        raw = self.mast_output_dir_var.get().strip() if hasattr(self, "mast_output_dir_var") else ""
        base = Path(raw).resolve() if raw and Path(raw).is_dir() else Path.home()
        initialdir = str(base)
        for sub in ("TESS", "Kepler", "K2", "MAST"):
            cand = base / sub
            if cand.is_dir():
                initialdir = str(cand)
                break
        p = filedialog.askopenfilename(
            title="FITS entrée à nettoyer",
            filetypes=[("FITS", "*.fits *.fits.gz"), ("Tous", "*.*")],
            initialdir=initialdir,
        )
        if p:
            self.fits_clean_input_var.set(p)
            self._default_fits_clean_output_from_input()

    def _default_fits_clean_output_from_input(self) -> None:
        src = self.fits_clean_input_var.get().strip()
        if not src:
            return
        try:
            p = Path(src)
            if p.name.lower().endswith(".fits.gz"):
                base = p.name[:-8]
            elif p.name.lower().endswith(".fits"):
                base = p.name[:-5]
            else:
                base = p.stem
            self.fits_clean_output_var.set(str(p.with_name(f"{base}_pdcsap_only.fits")))
        except Exception:
            pass

    def _browse_fits_clean_output(self) -> None:
        initial = self.fits_clean_output_var.get().strip()
        path = filedialog.asksaveasfilename(
            title="Enregistrer le FITS nettoyé sous...",
            defaultextension=".fits",
            filetypes=[("FITS", "*.fits"), ("Tous", "*.*")],
            initialdir=str(Path(initial).parent) if initial else str(Path.home()),
            initialfile=Path(initial).name if initial else "",
        )
        if path:
            self.fits_clean_output_var.set(path)

    def _run_clean_single_fits(self) -> None:
        src = Path(self.fits_clean_input_var.get().strip())
        dst_raw = self.fits_clean_output_var.get().strip()
        if not src.is_file():
            messagebox.showerror("Edition Fits", "FITS entrée introuvable.")
            return
        if not dst_raw:
            self._default_fits_clean_output_from_input()
            dst_raw = self.fits_clean_output_var.get().strip()
        if not dst_raw:
            messagebox.showerror("Edition Fits", "Veuillez définir un FITS de sortie.")
            return
        dst = Path(dst_raw).expanduser().resolve()
        try:
            tbl = self._extract_pdcsap_table_from_file(src)
            primary = fits.PrimaryHDU()
            hdu = fits.BinTableHDU(Table(tbl), name="LIGHTCURVE")
            fits.HDUList([primary, hdu]).writeto(dst, overwrite=True)
            self._fits_edit_log(f"Nettoyage OK : {src.name} -> {dst}")
            messagebox.showinfo("Edition Fits", f"FITS nettoyé créé :\n{dst}")
        except Exception as e:
            self._fits_edit_log(f"Nettoyage erreur: {e}")
            messagebox.showerror("Edition Fits", f"Erreur nettoyage FITS:\n{e}")

    def _refresh_merge_fits_listbox(self) -> None:
        if not hasattr(self, "fits_merge_listbox"):
            return
        self.fits_merge_listbox.delete(0, tk.END)
        for p in self.fits_merge_paths:
            self.fits_merge_listbox.insert(tk.END, p)

    def _add_merge_fits_files(self) -> None:
        raw = self.mast_output_dir_var.get().strip() if hasattr(self, "mast_output_dir_var") else ""
        base = Path(raw).resolve() if raw and Path(raw).is_dir() else Path.home()
        paths = filedialog.askopenfilenames(
            title="Ajouter des FITS à fusionner",
            filetypes=[("FITS", "*.fits *.fits.gz"), ("Tous", "*.*")],
            initialdir=str(base),
        )
        if not paths:
            return
        for p in paths:
            ap = str(Path(p).resolve())
            if ap not in self.fits_merge_paths:
                self.fits_merge_paths.append(ap)
        self._refresh_merge_fits_listbox()

    def _remove_merge_fits_selected(self) -> None:
        if not hasattr(self, "fits_merge_listbox"):
            return
        sel = list(self.fits_merge_listbox.curselection())
        if not sel:
            return
        for idx in sorted(sel, reverse=True):
            del self.fits_merge_paths[idx]
        self._refresh_merge_fits_listbox()

    def _clear_merge_fits_list(self) -> None:
        self.fits_merge_paths = []
        self._refresh_merge_fits_listbox()

    def _browse_fits_merge_output(self) -> None:
        initial = self.fits_merge_output_var.get().strip()
        path = filedialog.asksaveasfilename(
            title="Enregistrer le FITS fusionné sous...",
            defaultextension=".fits",
            filetypes=[("FITS", "*.fits"), ("Tous", "*.*")],
            initialdir=str(Path(initial).parent) if initial else str(Path.home()),
            initialfile=Path(initial).name if initial else "merged_pdcsap.fits",
        )
        if path:
            self.fits_merge_output_var.set(path)

    def _run_merge_fits(self) -> None:
        if not self.fits_merge_paths:
            messagebox.showerror("Edition Fits", "Ajoutez au moins un FITS à fusionner.")
            return
        out_raw = self.fits_merge_output_var.get().strip()
        if not out_raw:
            messagebox.showerror("Edition Fits", "Veuillez définir le FITS de sortie fusionné.")
            return
        out_path = Path(out_raw).expanduser().resolve()

        tables = []
        ok = 0
        for p in self.fits_merge_paths:
            fp = Path(p)
            if not fp.is_file():
                self._fits_edit_log(f"Ignoré (introuvable): {fp}")
                continue
            try:
                t = self._extract_pdcsap_table_from_file(fp)
                tables.append(t)
                ok += 1
                self._fits_edit_log(f"Ajouté au merge: {fp.name} ({len(t)} lignes)")
            except Exception as e:
                self._fits_edit_log(f"Ignoré (incompatible): {fp.name} ({e})")

        if not tables:
            messagebox.showerror("Edition Fits", "Aucun FITS compatible à fusionner.")
            return

        merged = vstack(tables, join_type="exact", metadata_conflicts="silent")
        if "TIME" in merged.colnames:
            try:
                merged.sort("TIME")
            except Exception:
                pass
        primary = fits.PrimaryHDU()
        hdu = fits.BinTableHDU(Table(merged), name="LIGHTCURVE")
        fits.HDUList([primary, hdu]).writeto(out_path, overwrite=True)
        self._fits_edit_log(f"Fusion OK : {ok} fichier(s) -> {out_path} ({len(merged)} lignes)")
        messagebox.showinfo(
            "Edition Fits",
            f"Fusion terminée.\nFichiers lus: {ok}\nLignes fusionnées: {len(merged)}\nSortie:\n{out_path}",
        )

    def _tess_lc_log(self, message: str) -> None:
        """Journal TESS/LcTools ; mise à jour GUI via after(0) (threads worker)."""
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {message}\n"
        logger.info(f"[Catalogues][TESS_LC_TXT] {message}")

        def _append() -> None:
            try:
                if self.winfo_exists() and getattr(self, "tess_lc_log_text", None) is not None:
                    self.tess_lc_log_text.insert(tk.END, line)
                    self.tess_lc_log_text.see(tk.END)
            except tk.TclError:
                pass

        try:
            self.after(0, _append)
        except tk.TclError:
            pass

    def _create_tess_workflow_tab(self, parent: ttk.Frame) -> None:
        """Flux LcTools uniquement : FITS -> normalisation -> .txt."""
        script = self._tess_lc_fits_script()
        if not script.is_file():
            ttk.Label(parent, text=f"Script introuvable : {script}", foreground="red").pack(anchor="w", pady=6)
            return

        # Formulaire défilant (évite que le journal en bas masque le bouton « Valider » étape 2)
        scroll_wrap = ttk.Frame(parent)
        scroll_wrap.pack(fill=tk.BOTH, expand=True, pady=(0, 2))
        canvas = tk.Canvas(scroll_wrap, highlightthickness=0)
        vsb = ttk.Scrollbar(scroll_wrap, orient="vertical", command=canvas.yview)
        form_root = ttk.Frame(canvas)
        _canvas_win = canvas.create_window((0, 0), window=form_root, anchor="nw")

        def _tess_scroll_region(_event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _tess_canvas_width(event):
            canvas.itemconfigure(_canvas_win, width=max(1, int(event.width)))

        form_root.bind("<Configure>", _tess_scroll_region)
        canvas.bind("<Configure>", _tess_canvas_width)

        def _tess_wheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind("<Enter>", lambda _e: canvas.bind_all("<MouseWheel>", _tess_wheel))
        canvas.bind("<Leave>", lambda _e: canvas.unbind_all("<MouseWheel>"))
        canvas.configure(yscrollcommand=vsb.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        ttk.Label(
            form_root,
            text=(
                "Normalisation LcTools à partir d'un FITS de courbe de lumière. "
                "Le flux est normalisé par la moyenne des points hors transit (OOT)."
            ),
            foreground="gray",
            wraplength=780,
            justify="left",
        ).pack(anchor="w", pady=(0, 10))

        # Variables priors (utilisées par la construction de commande)
        self.tess_wf_planet_var = tk.StringVar()
        self.tess_wf_skipeu_var = tk.BooleanVar(value=False)
        self.tess_wf_skipnasa_var = tk.BooleanVar(value=False)
        self.tess_wf_t0_var = tk.StringVar()
        self.tess_wf_p_var = tk.StringVar()
        self.tess_wf_t14_var = tk.StringVar()
        self.tess_wf_btjdoff_var = tk.StringVar(value="2457000.0")
        self.tess_wf_eu_to_var = tk.StringVar(value="120")
        self.tess_wf_pick_var = tk.BooleanVar(value=False)

        pri_fr = ttk.LabelFrame(
            form_root,
            text="Priors (catalogues : exoplanet.eu + NASA ; complétés ou remplacés par la saisie manuelle)",
            padding=6,
        )
        pri_fr.pack(fill="x", pady=6)

        prow = ttk.Frame(pri_fr)
        prow.pack(fill="x", pady=2)
        ttk.Label(prow, text="Nom planète (--planet-name) :", width=26).pack(side="left")
        ttk.Entry(prow, textvariable=self.tess_wf_planet_var, width=28).pack(side="left", padx=4)
        ttk.Checkbutton(prow, text="Sans exoplanet.eu", variable=self.tess_wf_skipeu_var).pack(side="left", padx=6)
        ttk.Checkbutton(prow, text="Sans NASA", variable=self.tess_wf_skipnasa_var).pack(side="left")

        man_fr = ttk.Frame(pri_fr)
        man_fr.pack(fill="x", pady=2)
        ttk.Label(man_fr, text="Manuel (optionnel) — T₀ BTJD :", width=28).pack(side="left")
        ttk.Entry(man_fr, textvariable=self.tess_wf_t0_var, width=12).pack(side="left", padx=2)
        ttk.Label(man_fr, text="P (j) :").pack(side="left")
        ttk.Entry(man_fr, textvariable=self.tess_wf_p_var, width=12).pack(side="left", padx=2)
        ttk.Label(man_fr, text="T₁₄ (h) :").pack(side="left")
        ttk.Entry(man_fr, textvariable=self.tess_wf_t14_var, width=10).pack(side="left", padx=2)
        ttk.Label(man_fr, text="offset BTJD→BJD :", foreground="gray").pack(side="left", padx=(8, 0))
        ttk.Entry(man_fr, textvariable=self.tess_wf_btjdoff_var, width=10).pack(side="left", padx=2)

        adv_fr = ttk.Frame(pri_fr)
        adv_fr.pack(fill="x", pady=2)
        ttk.Label(adv_fr, text="Timeout TAP EU (s) :").pack(side="left", padx=4)
        ttk.Entry(adv_fr, textvariable=self.tess_wf_eu_to_var, width=8).pack(side="left", padx=4)
        ttk.Checkbutton(
            adv_fr,
            text="Ouvrir l’éditeur matplotlib des fenêtres de transit",
            variable=self.tess_wf_pick_var,
        ).pack(side="left", padx=12)

        # --- Étape unique : normalisation LcTools -> .txt ---
        sec_lc = ttk.LabelFrame(
            form_root,
            text="1 — Normalisation LcTools → .txt (BJD-TDB ; flux / moyenne hors transit uniquement)",
            padding=8,
        )
        sec_lc.pack(fill="x", pady=8)

        lc_in_row = ttk.Frame(sec_lc)
        lc_in_row.pack(fill="x", pady=2)
        ttk.Label(lc_in_row, text="Fichier FITS entrée :", width=24).pack(side="left")
        self.tess_wf_lc_fits_var = tk.StringVar()
        ttk.Entry(lc_in_row, textvariable=self.tess_wf_lc_fits_var, width=46).pack(
            side="left", fill="x", expand=True, padx=4
        )
        fit_btns = ttk.Frame(lc_in_row)
        fit_btns.pack(side="left", padx=(4, 0))
        ttk.Button(
            fit_btns,
            text="📂",
            width=3,
            command=lambda: self._tess_lc_browse_fits(self.tess_wf_lc_fits_var),
        ).pack(side="left", padx=(0, 2))
        ttk.Button(
            fit_btns,
            text="MAST",
            width=5,
            command=self._tess_wf_browse_fits_from_mast_output_dir,
        ).pack(side="left")

        lc_out_row = ttk.Frame(sec_lc)
        lc_out_row.pack(fill="x", pady=2)
        ttk.Label(lc_out_row, text="Sortie .txt (-o), optionnel :", width=24).pack(side="left")
        self.tess_wf_lc_txt_out_var = tk.StringVar()
        ttk.Entry(lc_out_row, textvariable=self.tess_wf_lc_txt_out_var, width=46).pack(
            side="left", fill="x", expand=True, padx=4
        )
        ttk.Label(
            lc_out_row,
            text="Vide = même dossier / même base + .txt",
            foreground="gray",
        ).pack(side="left")

        ttk.Label(
            sec_lc,
            text=(
                "Priors : nom planète, T₀, P, T₁₄ et catalogues exoplanet.eu/NASA. "
                "Utilisez l’onglet « Edition Fits » si vous souhaitez nettoyer/fusionner des FITS avant normalisation."
            ),
            foreground="gray",
            wraplength=700,
            justify="left",
        ).pack(anchor="w", pady=(4, 2))

        lc_opt = ttk.LabelFrame(sec_lc, text="Options de normalisation", padding=6)
        lc_opt.pack(fill="x", pady=6)

        ttk.Label(
            lc_opt,
            text=(
                "La normalisation utilise **uniquement** les points **hors transit** pour calculer la moyenne de référence "
                "(équivalent --norm-baseline oot ; pas de normalisation sur toute la série)."
            ),
            foreground="gray",
            wraplength=700,
            justify="left",
        ).pack(anchor="w", pady=(0, 6))

        lc_adv = ttk.Frame(lc_opt)
        lc_adv.pack(fill="x", pady=2)
        self.tess_wf_lc_gq_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(lc_adv, text="QUALITY==0 uniquement", variable=self.tess_wf_lc_gq_var).pack(
            side="left", padx=4
        )
        self.tess_wf_lc_pick_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            lc_adv,
            text="Ré-ouvrir l’éditeur matplotlib (affinage des bandes transit avant export)",
            variable=self.tess_wf_lc_pick_var,
        ).pack(side="left", padx=12)

        ttk.Separator(lc_opt, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(8, 6))
        ttk.Button(
            lc_opt,
            text="✓ Valider — normaliser et exporter (LcTools .txt)",
            command=self.start_tess_wf_norm_lctools_thread,
            width=42,
        ).pack(anchor="w", pady=(0, 4))
        ttk.Label(
            lc_opt,
            text="(Lance l’export .txt à partir du FITS sélectionné ci-dessus.)",
            foreground="gray",
            font=("Helvetica", 8),
            wraplength=680,
            justify="left",
        ).pack(anchor="w", pady=(0, 2))

        log_fr = ttk.LabelFrame(parent, text="Journal", padding=4)
        log_fr.pack(fill=tk.X, expand=False, pady=6)
        lf = ttk.Frame(log_fr)
        lf.pack(fill=tk.X, expand=False)
        self.tess_lc_log_text = tk.Text(
            lf,
            wrap=tk.WORD,
            font=("Consolas", 9),
            height=8,
            bg="#1a1a1a",
            fg="#e8e8e8",
            insertbackground="#e8e8e8",
        )
        sb = ttk.Scrollbar(lf, orient="vertical", command=self.tess_lc_log_text.yview)
        self.tess_lc_log_text.configure(yscrollcommand=sb.set)
        self.tess_lc_log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side="right", fill="y")
        self._tess_lc_log(f"Script : {script}")
        self._tess_lc_log(f"Python : {sys.executable}")

    def _tess_lc_browse_fits(self, var: tk.StringVar) -> None:
        path = filedialog.askopenfilename(
            title="Fichier courbe TESS *_lc.fits",
            filetypes=[("FITS", "*.fits *.fits.gz"), ("Tous", "*.*")],
            initialdir=str(Path(var.get()).parent) if var.get() else str(Path.home()),
        )
        if path:
            var.set(path)

    def _tess_wf_browse_fits_from_mast_output_dir(self) -> None:
        """Ouvre le sélecteur de fichier dans le dossier sortie MAST (onglet 📡 Téléchargement MAST)."""
        raw = self.mast_output_dir_var.get().strip() if hasattr(self, "mast_output_dir_var") else ""
        base = Path(raw).resolve() if raw and Path(raw).is_dir() else Path.home()
        initialdir = str(base)
        for sub in ("TESS", "Kepler", "K2", "MAST"):
            cand = base / sub
            if cand.is_dir():
                initialdir = str(cand)
                break
        path = filedialog.askopenfilename(
            title="Courbe *_lc.fits (dossier plat TESS / Kepler / … sous le répertoire MAST)",
            filetypes=[("FITS courbe", "*.fits *.fits.gz"), ("Tous", "*.*")],
            initialdir=initialdir,
        )
        if path:
            if hasattr(self, "tess_wf_lc_fits_var"):
                self.tess_wf_lc_fits_var.set(path)
            # Proposition automatique du .txt de sortie
            try:
                p = Path(path)
                if hasattr(self, "tess_wf_lc_txt_out_var") and not self.tess_wf_lc_txt_out_var.get().strip():
                    self.tess_wf_lc_txt_out_var.set(str(p.with_name(p.stem + "_lctools.txt")))
            except Exception:
                pass

    def _tess_wf_sync_lc_from_edited(self) -> None:
        """FITS étape 1 → champ normalisation LcTools ; propose un .txt de sortie distinct du brut."""
        p = self.tess_wf_edited_fits_out_var.get().strip()
        if not p:
            messagebox.showinfo(
                "Étape 1",
                "Indiquez d’abord le FITS de sortie (étape 1) ou enregistrez avec « Valider ».",
            )
            return
        self.tess_wf_lc_fits_var.set(p)
        try:
            fp = Path(p)
            if not self.tess_wf_lc_txt_out_var.get().strip():
                self.tess_wf_lc_txt_out_var.set(str(fp.with_name(fp.stem + "_lctools.txt")))
        except Exception:
            pass

    def _tess_wf_default_edited_out_from_input(self) -> None:
        raw_in = self.tess_wf_norm_input_var.get().strip()
        if not raw_in:
            return
        try:
            p = Path(raw_in)
            name = p.name
            if name.lower().endswith(".fits.gz"):
                stem = name[:-8]
            elif name.lower().endswith(".fits"):
                stem = name[:-5]
            else:
                stem = p.stem
            default = p.parent / f"{stem}_edited.fits"
            if not self.tess_wf_edited_fits_out_var.get().strip():
                self.tess_wf_edited_fits_out_var.set(str(default))
        except Exception:
            pass

    def _tess_wf_browse_fits_and_default_edited_out(self) -> None:
        self._tess_lc_browse_fits(self.tess_wf_norm_input_var)
        self._tess_wf_default_edited_out_from_input()

    def _tess_wf_browse_edited_fits_save(self) -> None:
        initial = self.tess_wf_edited_fits_out_var.get().strip()
        initialdir = str(Path(initial).parent) if initial else str(Path.home())
        initialfile = Path(initial).name if initial else ""
        path = filedialog.asksaveasfilename(
            title="Enregistrer le FITS édité sous…",
            defaultextension=".fits",
            filetypes=[("FITS", "*.fits"), ("Tous", "*.*")],
            initialdir=initialdir,
            initialfile=initialfile,
        )
        if path:
            self.tess_wf_edited_fits_out_var.set(path)

    def start_tess_wf_edit_fits_thread(self) -> None:
        threading.Thread(target=self.run_tess_wf_edit_fits_export, daemon=True).start()

    def start_tess_wf_norm_lctools_thread(self) -> None:
        threading.Thread(target=self.run_tess_wf_norm_lctools_export, daemon=True).start()

    def _build_tess_wf_norm_lctools_cmd(self) -> list[str] | None:
        """--normalized-lctools sur le FITS (typiquement celui de l’étape 1)."""
        script = self._tess_lc_fits_script()
        if not script.is_file():
            return None
        inp = Path(self.tess_wf_lc_fits_var.get().strip())
        if not inp.is_file():
            return None
        cmd: list[str] = [sys.executable, str(script), str(inp.resolve())]
        cmd.append("--normalized-lctools")
        # Toujours baseline hors transit uniquement (moyenne OOT pour normaliser le flux).
        out_txt = self.tess_wf_lc_txt_out_var.get().strip()
        if out_txt:
            cmd.extend(["-o", str(Path(out_txt).expanduser().resolve())])
        if self.tess_wf_lc_gq_var.get():
            cmd.append("--good-quality-only")
        planet = self.tess_wf_planet_var.get().strip()
        if planet:
            cmd.extend(["--planet-name", planet])
        if self.tess_wf_skipeu_var.get():
            cmd.append("--skip-exoplanet-eu")
        if self.tess_wf_skipnasa_var.get():
            cmd.append("--skip-nasa-ephemeris")
        t0s = self.tess_wf_t0_var.get().strip()
        if t0s:
            cmd.extend(["--t0-btjd", t0s])
        ps = self.tess_wf_p_var.get().strip()
        if ps:
            cmd.extend(["--period-days", ps])
        t14 = self.tess_wf_t14_var.get().strip()
        if t14:
            cmd.extend(["--duration-hours", t14])
        off = self.tess_wf_btjdoff_var.get().strip()
        if off and off != "2457000.0":
            try:
                float(off)
                cmd.extend(["--btjd-offset", off])
            except ValueError:
                pass
        try:
            eu_to = float(self.tess_wf_eu_to_var.get().strip() or "120")
            cmd.extend(["--exoplanet-eu-timeout", str(eu_to)])
        except ValueError:
            pass
        if self.tess_wf_lc_pick_var.get():
            cmd.append("--pick-transit-windows")
        return cmd

    def _build_tess_wf_edit_fits_cmd(self) -> list[str] | None:
        """Édition masques + --edit-save-fits (pas d’export LcTools .txt)."""
        script = self._tess_lc_fits_script()
        if not script.is_file():
            return None
        inp = Path(self.tess_wf_norm_input_var.get().strip())
        if not inp.is_file():
            return None
        out_fits = self.tess_wf_edited_fits_out_var.get().strip()
        if not out_fits:
            return None
        out_path = Path(out_fits).expanduser().resolve()
        cmd: list[str] = [
            sys.executable,
            str(script),
            str(inp.resolve()),
            "--edit-save-fits",
            str(out_path),
        ]
        planet = self.tess_wf_planet_var.get().strip()
        if planet:
            cmd.extend(["--planet-name", planet])
        if self.tess_wf_skipeu_var.get():
            cmd.append("--skip-exoplanet-eu")
        if self.tess_wf_skipnasa_var.get():
            cmd.append("--skip-nasa-ephemeris")
        t0s = self.tess_wf_t0_var.get().strip()
        if t0s:
            cmd.extend(["--t0-btjd", t0s])
        ps = self.tess_wf_p_var.get().strip()
        if ps:
            cmd.extend(["--period-days", ps])
        t14 = self.tess_wf_t14_var.get().strip()
        if t14:
            cmd.extend(["--duration-hours", t14])
        off = self.tess_wf_btjdoff_var.get().strip()
        if off and off != "2457000.0":
            try:
                float(off)
                cmd.extend(["--btjd-offset", off])
            except ValueError:
                pass
        try:
            eu_to = float(self.tess_wf_eu_to_var.get().strip() or "120")
            cmd.extend(["--exoplanet-eu-timeout", str(eu_to)])
        except ValueError:
            pass
        if self.tess_wf_pick_var.get():
            cmd.append("--pick-transit-windows")
        return cmd

    def run_tess_wf_norm_lctools_export(self) -> None:
        cmd = self._build_tess_wf_norm_lctools_cmd()
        if cmd is None:
            self._tess_lc_log("[LcTools] Erreur : FITS introuvable (étape 2).")
            self.after(
                0,
                lambda: messagebox.showerror(
                    "Normalisation LcTools",
                    "Choisissez le FITS de l’étape 1 (bouton ↪ ou parcourir).",
                ),
            )
            return
        self._tess_lc_log("[LcTools] Commande : " + " ".join(f'"{c}"' if " " in c else c for c in cmd))
        cwd = str(self._npoap_project_root())
        try:
            proc = subprocess.run(
                cmd,
                cwd=cwd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=None,
            )
        except Exception as e:
            self._tess_lc_log(f"[LcTools] Erreur : {e}")
            self.after(0, lambda: messagebox.showerror("Normalisation LcTools", str(e)))
            return
        if proc.stdout:
            for line in proc.stdout.splitlines():
                self._tess_lc_log(line)
        if proc.stderr:
            for line in proc.stderr.splitlines():
                self._tess_lc_log("[stderr] " + line)
        if proc.returncode == 0:
            self._tess_lc_log("[LcTools] Terminé avec succès.")
            self.after(0, lambda: messagebox.showinfo("Normalisation LcTools", "Export .txt normalisé terminé."))
        else:
            self.after(
                0,
                lambda: messagebox.showerror("Normalisation LcTools", f"Code retour {proc.returncode}."),
            )

    def run_tess_wf_edit_fits_export(self) -> None:
        cmd = self._build_tess_wf_edit_fits_cmd()
        if cmd is None:
            self._tess_lc_log("[FITS édité] Erreur : FITS d’entrée introuvable ou chemin de sortie vide.")
            self.after(
                0,
                lambda: messagebox.showerror(
                    "FITS édité",
                    "Choisissez un *_lc.fits d’entrée existant et un chemin de sortie pour le FITS enregistré.",
                ),
            )
            return
        self._tess_lc_log("[FITS édité] Commande : " + " ".join(f'"{c}"' if " " in c else c for c in cmd))
        cwd = str(self._npoap_project_root())
        try:
            proc = subprocess.run(
                cmd,
                cwd=cwd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=None,
            )
        except Exception as e:
            self._tess_lc_log(f"[FITS édité] Erreur : {e}")
            self.after(0, lambda: messagebox.showerror("FITS édité", str(e)))
            return
        if proc.stdout:
            for line in proc.stdout.splitlines():
                self._tess_lc_log(line)
        if proc.stderr:
            for line in proc.stderr.splitlines():
                self._tess_lc_log("[stderr] " + line)
        if proc.returncode == 0:
            self._tess_lc_log("[FITS édité] Terminé avec succès.")
            self.after(0, lambda: messagebox.showinfo("FITS édité", "FITS enregistré (transits → NaN sur les flux)."))
        else:
            self.after(
                0,
                lambda: messagebox.showerror("FITS édité", f"Code retour {proc.returncode}."),
            )

    def start_mast_lightcurves_download_thread(self):
        """Lance le téléchargement MAST dans un thread séparé."""
        thread = threading.Thread(target=self.run_mast_lightcurves_download, daemon=True)
        thread.start()

    def _mast_dialog(self, message: str):
        """Écrit un message verbeux dans la fenêtre de dialogue MAST."""
        try:
            ts = datetime.now().strftime("%H:%M:%S")
            logger.info(f"[Catalogues][MAST] {message}")
            if hasattr(self, "mast_dialog_text") and self.mast_dialog_text is not None:
                self.mast_dialog_text.insert(tk.END, f"[{ts}] {message}\n")
                self.mast_dialog_text.see(tk.END)
        except Exception:
            pass

    def _call_with_timeout(self, func, timeout_s: int, step_name: str, *args, **kwargs):
        """
        Exécute un appel potentiellement bloquant avec timeout.
        Retourne le résultat ou lève TimeoutError/Exception.
        """
        self._mast_dialog(f"Étape '{step_name}': délai max {timeout_s}s.")
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        future = executor.submit(func, *args, **kwargs)
        try:
            return future.result(timeout=timeout_s)
        except concurrent.futures.TimeoutError:
            # Important: ne pas attendre la fin de la tâche bloquée
            future.cancel()
            executor.shutdown(wait=False, cancel_futures=True)
            # from None : évite une traceback en chaîne confuse dans les logs
            raise TimeoutError(
                f"Timeout sur l'étape '{step_name}' après {timeout_s}s."
            ) from None
        finally:
            # Si non-timeout, libère proprement le worker
            if not future.cancelled():
                executor.shutdown(wait=False, cancel_futures=True)

    @staticmethod
    def _normalize_mast_target(raw: str) -> str:
        """Normalise TOI/TIC pour de meilleures résolutions MAST/lightkurve."""
        s = " ".join(str(raw).strip().split())
        if not s:
            return s
        low = s.lower()
        m = re.match(r"^toi[\s-]+(\d+(?:\.\d+)?)$", low)
        if m:
            return f"TOI {m.group(1)}"
        m = re.match(r"^tic[\s-]+(\d+)$", low)
        if m:
            return f"TIC {m.group(1)}"
        if re.match(r"^\d{8,}$", low):
            return f"TIC {low}"
        return s

    def _mast_connectivity_probe(self) -> None:
        """Test rapide HTTPS vers MAST (diagnostic réseau / pare-feu)."""
        if not REQUESTS_AVAILABLE:
            self._mast_dialog("Test connectivité: requests non disponible, ignoré.")
            return
        try:
            r = requests.head("https://mast.stsci.edu", timeout=10, allow_redirects=True)
            self._mast_dialog(f"Test connectivité mast.stsci.edu: HTTP {r.status_code}")
        except Exception as ex:
            self._mast_dialog(f"Test connectivité mast.stsci.edu: ÉCHEC ({ex})")

    @staticmethod
    def _lightkurve_mission_kw(mission: str):
        if mission == "TESS":
            return "TESS"
        if mission == "Kepler":
            return "Kepler"
        if mission == "K2":
            return "K2"
        return None

    @staticmethod
    def _lightkurve_supported_mission(mission: str) -> bool:
        return mission in ("", "TESS", "Kepler", "K2")

    @staticmethod
    def _mast_is_fits_product(path: Path) -> bool:
        """True si fichier FITS courbe / produit MAST typique (.fits ou .fits.gz)."""
        if not path.is_file():
            return False
        name = path.name.lower()
        return name.endswith(".fits") or name.endswith(".fits.gz")

    @staticmethod
    def _mast_mission_flat_subdir(mission: str) -> str:
        """Sous-dossier plat pour les FITS après téléchargement (pas d’arborescence mastDownload)."""
        m = (mission or "").strip().upper()
        if m == "TESS":
            return "TESS"
        if m == "KEPLER":
            return "Kepler"
        if m == "K2":
            return "K2"
        return "MAST"

    def _mast_unique_flat_destination(self, dest_dir: Path, src: Path) -> Path:
        """Nom de fichier unique dans dest_dir (évite collisions entre secteurs / pipelines)."""
        name = src.name
        candidate = dest_dir / name
        if not candidate.exists():
            return candidate
        stem = src.stem
        suffixes = "".join(src.suffixes)
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", src.parent.name)[:48]
        for n in range(2, 100000):
            cand = dest_dir / f"{slug}__{stem}_{n}{suffixes}"
            if not cand.exists():
                return cand
        return dest_dir / f"{slug}__{stem}_overflow{suffixes}"

    def _mast_flatten_fits_to_folder(self, source_root: Path, dest_dir: Path) -> int:
        """
        Déplace tous les FITS trouvés sous source_root vers dest_dir (structure plate).
        Retourne le nombre de fichiers déplacés.
        """
        if not source_root.is_dir():
            return 0
        dest_dir.mkdir(parents=True, exist_ok=True)
        moved = 0
        # Parcours : ne pas prendre les fichiers déjà dans dest_dir
        try:
            dest_resolved = dest_dir.resolve()
        except OSError:
            dest_resolved = dest_dir
        move_errors = 0
        for p in sorted(source_root.rglob("*")):
            if not self._mast_is_fits_product(p):
                continue
            try:
                if dest_resolved in p.resolve().parents or p.resolve().parent == dest_resolved:
                    continue
            except OSError:
                pass
            target = self._mast_unique_flat_destination(dest_dir, p)
            try:
                shutil.move(str(p), str(target))
                moved += 1
            except OSError as ex:
                move_errors += 1
                self._mast_dialog(f"Aplatissement: impossible de déplacer {p.name} → {target.name}: {ex}")
                logger.warning("MAST flatten move failed: %s", ex)
        # Supprimer le staging seulement si tout a été déplacé (évite perte de FITS en cas d’erreur)
        if move_errors == 0:
            try:
                shutil.rmtree(source_root, ignore_errors=True)
            except Exception:
                pass
        else:
            self._mast_dialog(
                f"Aplatissement partiel : {move_errors} erreur(s), dossier staging conservé : {source_root}"
            )
        self._mast_dialog(f"Aplatissement: {moved} fichier(s) FITS → {dest_dir}")
        return moved

    def _mast_reset_staging_dir(self, staging_dir: Path) -> None:
        """Répertoire temporaire pour téléchargements MAST (arborescence type mastDownload), puis vidé après aplatissement."""
        if staging_dir.exists():
            shutil.rmtree(staging_dir, ignore_errors=True)
        staging_dir.mkdir(parents=True, exist_ok=True)

    def _mast_tess_lightkurve_authors(self) -> list[str] | None:
        """
        Auteurs MAST passés à lightkurve (author=...) pour la mission TESS.
        None = ne pas filtrer par auteur (recherche large, tous produits retournés par MAST).
        Liste non vide = une requête par auteur (SPOC, QLP, GSFC-ELEANOR-LITE, …).
        """
        if not hasattr(self, "mast_tess_spoc_var"):
            return None
        authors: list[str] = []
        if self.mast_tess_spoc_var.get():
            authors.append("SPOC")
        if self.mast_tess_qlp_var.get():
            authors.append("QLP")
        if self.mast_tess_eleanor_var.get():
            authors.append("GSFC-ELEANOR-LITE")
        if not authors:
            return None
        return authors

    def _should_try_tic_catalog_path(self, target_norm: str) -> bool:
        t = target_norm.lower()
        if "tic" in t:
            return True
        return bool(re.match(r"^\d{8,}$", target_norm.strip()))

    def _tic_resolve_to_skycoord(self, target_norm: str, radius_deg: float, query_timeout_s: int):
        """
        Résout un TIC via Catalogs uniquement (rapide). Ne passe pas par Observations.query_region.
        Retourne SkyCoord ou None.
        """
        if Catalogs is None:
            return None
        tic_query = target_norm
        if re.match(r"^\d{8,}$", target_norm.strip()):
            tic_query = f"TIC {target_norm.strip()}"
        self._mast_dialog(f"Catalogs.query_object('{tic_query}', catalog='TIC')...")
        tic_table = self._call_with_timeout(
            lambda: Catalogs.query_object(tic_query, catalog="TIC", radius=f"{radius_deg} deg"),
            query_timeout_s,
            "Catalogs.query_object(TIC)",
        )
        if tic_table is None or len(tic_table) == 0:
            self._mast_dialog("TIC: aucune ligne catalogue.")
            return None
        row = tic_table[0]
        ra_col = "ra" if "ra" in tic_table.colnames else ("RA" if "RA" in tic_table.colnames else None)
        dec_col = "dec" if "dec" in tic_table.colnames else ("DEC" if "DEC" in tic_table.colnames else None)
        if not ra_col or not dec_col:
            self._mast_dialog(f"TIC: colonnes ra/dec introuvables ({tic_table.colnames[:12]})")
            return None
        ra_v = row[ra_col]
        dec_v = row[dec_col]
        if hasattr(ra_v, "value"):
            ra_v = ra_v.value
        if hasattr(dec_v, "value"):
            dec_v = dec_v.value
        coord = SkyCoord(ra=float(ra_v) * u.deg, dec=float(dec_v) * u.deg)
        self._mast_dialog(f"TIC résolu -> RA={float(ra_v):.6f}°, Dec={float(dec_v):.6f}°")
        return coord

    def _observations_via_tic_coords(self, target_norm: str, radius_deg: float, query_timeout_s: int):
        """Secours: TIC via Catalogs puis Observations.query_region (souvent bloquant / timeout)."""
        coord = self._tic_resolve_to_skycoord(target_norm, radius_deg, query_timeout_s)
        if coord is None:
            return None
        self._mast_dialog("Observations.query_region(coord, ...) — peut rester sans réponse (API Portal MAST).")
        return self._call_with_timeout(
            lambda: Observations.query_region(coord, radius=f"{radius_deg} deg"),
            query_timeout_s,
            "Observations.query_region",
        )

    def _lightkurve_download_chunk(
        self,
        search,
        output_dir: Path,
        n_rows: int,
        log_prefix: str,
    ) -> list:
        """
        Télécharge les n_rows premières lignes d'un SearchResult lightkurve (itératif, tolérant aux erreurs).
        """
        results: list = []
        if search is None or n_rows <= 0:
            return results
        failed = 0
        self._mast_dialog(f"{log_prefix} téléchargement itératif ({n_rows} entrée(s))")
        for i in range(n_rows):
            one = search[i : i + 1]
            try:
                if hasattr(one, "download_all"):
                    results.append(one.download_all(download_dir=str(output_dir)))
                else:
                    results.append(one.download(download_dir=str(output_dir)))
            except Exception as ex:
                failed += 1
                results.append(None)
                self._mast_dialog(f"{log_prefix} entrée {i + 1}/{n_rows} ignorée ({ex})")
        self._mast_dialog(f"{log_prefix} fin chunk, échecs ignorés={failed}")
        return results

    def _download_with_lightkurve(
        self,
        target_norm: str,
        mission: str,
        output_dir: Path,
        max_products: int,
        search_timeout_s: int,
        download_timeout_s: int,
        tess_authors: list[str] | None = None,
    ) -> bool:
        """Télécharge via lightkurve (souvent plus fiable que Observations.query_object seul)."""
        if not LIGHTKURVE_AVAILABLE or lk is None:
            return False
        m_kw = self._lightkurve_mission_kw(mission)
        use_per_author = mission == "TESS" and bool(tess_authors)

        def do_download():
            if use_per_author and tess_authors:
                all_results: list = []
                rem = max_products if max_products > 0 else 10**9
                for auth in tess_authors:
                    if rem <= 0:
                        break

                    def do_search(author=auth):
                        return lk.search_lightcurve(target_norm, mission=m_kw, author=author)

                    self._mast_dialog(
                        f"lightkurve.search_lightcurve('{target_norm}', mission={m_kw!r}, author={auth!r})..."
                    )
                    try:
                        search = self._call_with_timeout(
                            do_search, search_timeout_s, f"lightkurve_search_{auth}"
                        )
                    except Exception as ex:
                        self._mast_dialog(f"lightkurve author={auth}: échec recherche ({ex})")
                        logger.warning("lightkurve author=%s: %s", auth, ex)
                        continue
                    if search is None or len(search) == 0:
                        self._mast_dialog(f"lightkurve author={auth}: 0 résultat.")
                        continue
                    self._mast_dialog(
                        f"lightkurve author={auth}: {len(search)} jeu(x) de courbes disponible(s)."
                    )
                    n = min(len(search), rem)
                    chunk = self._lightkurve_download_chunk(
                        search, output_dir, n, f"lightkurve[{auth}]"
                    )
                    all_results.extend(chunk)
                    rem -= n
                return all_results

            def do_search():
                return lk.search_lightcurve(target_norm, mission=m_kw)

            self._mast_dialog(f"lightkurve.search_lightcurve('{target_norm}', mission={m_kw!r})...")
            search = self._call_with_timeout(do_search, search_timeout_s, "lightkurve_search_lightcurve")
            if search is None or len(search) == 0:
                self._mast_dialog("lightkurve: 0 résultat.")
                return []
            self._mast_dialog(f"lightkurve: {len(search)} jeu(x) de courbes disponible(s).")
            n = min(len(search), max_products) if max_products > 0 else len(search)
            return self._lightkurve_download_chunk(search, output_dir, n, "lightkurve")

        self._mast_dialog(f"lightkurve: téléchargement vers {output_dir}...")
        downloaded_results = self._call_with_timeout(do_download, download_timeout_s, "lightkurve_download")
        success_count = sum(1 for r in (downloaded_results or []) if r is not None)
        if success_count == 0:
            self._mast_dialog("lightkurve: aucune entrée valide téléchargée.")
            return False
        self._mast_dialog(f"lightkurve: entrées valides téléchargées={success_count}")
        self._mast_dialog("lightkurve: téléchargement terminé (voir dossier de sortie).")
        return True

    def _download_with_lightkurve_coords(
        self,
        coord: SkyCoord,
        radius_deg: float,
        mission: str,
        output_dir: Path,
        max_products: int,
        search_timeout_s: int,
        download_timeout_s: int,
        tess_authors: list[str] | None = None,
    ) -> bool:
        """
        Recherche lightkurve par position (évite Observations.query_region qui bloque souvent).
        """
        if not LIGHTKURVE_AVAILABLE or lk is None:
            return False
        m_kw = self._lightkurve_mission_kw(mission)
        rad = float(radius_deg) * u.deg
        use_per_author = mission == "TESS" and bool(tess_authors)

        def do_download():
            if use_per_author and tess_authors:
                all_results: list = []
                rem = max_products if max_products > 0 else 10**9
                for auth in tess_authors:
                    if rem <= 0:
                        break

                    def do_search(author=auth):
                        return lk.search_lightcurve(coord, radius=rad, mission=m_kw, author=author)

                    self._mast_dialog(
                        f"lightkurve.search_lightcurve(SkyCoord, r={radius_deg}°, mission={m_kw!r}, author={auth!r})..."
                    )
                    try:
                        search = self._call_with_timeout(
                            do_search, search_timeout_s, f"lightkurve_search_coords_{auth}"
                        )
                    except Exception as ex:
                        self._mast_dialog(f"lightkurve (coords) author={auth}: échec recherche ({ex})")
                        logger.warning("lightkurve coords author=%s: %s", auth, ex)
                        continue
                    if search is None or len(search) == 0:
                        self._mast_dialog(f"lightkurve (coords) author={auth}: 0 résultat.")
                        continue
                    self._mast_dialog(
                        f"lightkurve (coords) author={auth}: {len(search)} jeu(x) disponible(s)."
                    )
                    n = min(len(search), rem)
                    chunk = self._lightkurve_download_chunk(
                        search, output_dir, n, f"lightkurve(coords)[{auth}]"
                    )
                    all_results.extend(chunk)
                    rem -= n
                return all_results

            def do_search():
                return lk.search_lightcurve(coord, radius=rad, mission=m_kw)

            self._mast_dialog(
                f"lightkurve.search_lightcurve(SkyCoord, radius={radius_deg} deg, mission={m_kw!r}) "
                f"(sans query_region Portal)..."
            )
            search = self._call_with_timeout(do_search, search_timeout_s, "lightkurve_search_coords")
            if search is None or len(search) == 0:
                self._mast_dialog("lightkurve (coords): 0 résultat.")
                return []
            self._mast_dialog(f"lightkurve (coords): {len(search)} jeu(x) de courbes disponible(s).")
            n = min(len(search), max_products) if max_products > 0 else len(search)
            return self._lightkurve_download_chunk(search, output_dir, n, "lightkurve(coords)")

        self._mast_dialog(f"lightkurve (coords): téléchargement vers {output_dir}...")
        downloaded_results = self._call_with_timeout(do_download, download_timeout_s, "lightkurve_download_coords")
        success_count = sum(1 for r in (downloaded_results or []) if r is not None)
        if success_count == 0:
            self._mast_dialog("lightkurve (coords): aucune entrée valide téléchargée.")
            return False
        self._mast_dialog(f"lightkurve (coords): entrées valides téléchargées={success_count}")
        self._mast_dialog("lightkurve (coords): téléchargement terminé.")
        return True

    def run_mast_lightcurves_download(self):
        """Télécharge des courbes de lumière depuis MAST (lightkurve puis astroquery)."""
        try:
            if not MAST_AVAILABLE:
                self.log_message("❌ astroquery.mast non disponible.", "error")
                messagebox.showerror("Erreur", "astroquery.mast non disponible. Installez astroquery.")
                return

            target_raw = self.mast_target_var.get().strip()
            if not target_raw:
                self.log_message("❌ Veuillez saisir une cible.", "error")
                return

            target = self._normalize_mast_target(target_raw)
            engine = getattr(self, "mast_engine_var", None)
            engine = (engine.get().strip().lower() if engine else "auto")

            try:
                radius_deg = float(self.mast_radius_var.get().strip() or "0.02")
            except ValueError:
                self.log_message("⚠️ Rayon invalide, valeur 0.02 deg utilisée.", "warning")
                radius_deg = 0.02

            try:
                max_products = int(self.mast_max_products_var.get().strip() or "200")
            except ValueError:
                self.log_message("⚠️ Max produits invalide, valeur 200 utilisée.", "warning")
                max_products = 200

            mission = self.mast_mission_var.get().strip()
            output_dir = Path(self.mast_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            flat_sub = self._mast_mission_flat_subdir(mission)
            dest_flat = output_dir / flat_sub
            staging_dir = output_dir / ".mast_staging"
            self._mast_reset_staging_dir(staging_dir)
            self._mast_dialog(
                f"Sortie : FITS regroupés à plat dans « {dest_flat.name} » (staging temporaire : .mast_staging)."
            )
            # Snapshot avant téléchargement pour mesurer les fichiers réellement écrits.
            before_fits = {p.resolve() for p in output_dir.rglob("*") if p.is_file() and self._mast_is_fits_product(p)}
            query_timeout_s = 120
            products_timeout_s = 120
            download_timeout_s = 900
            # Recherches MAST via lightkurve peuvent être lentes (réseau / serveurs)
            lk_search_timeout_s = 300

            self._mast_dialog("Nouvelle session de requête démarrée.")
            self._mast_dialog(
                f"astroquery.mast={MAST_AVAILABLE}, lightkurve={LIGHTKURVE_AVAILABLE}, moteur={engine}"
            )
            self._mast_connectivity_probe()
            if target_raw != target:
                self._mast_dialog(f"Cible normalisée: '{target_raw}' -> '{target}'")
            self._mast_dialog(f"Paramètres -> radius='{radius_deg} deg', mission='{mission or 'ALL'}'")

            tess_authors = self._mast_tess_lightkurve_authors() if mission == "TESS" else None
            if mission == "TESS":
                if tess_authors:
                    self._mast_dialog(f"Pipelines TESS (lightkurve author) : {', '.join(tess_authors)}")
                else:
                    self._mast_dialog("Pipelines TESS : aucune case cochée → recherche sans filtre auteur (tous produits)")

            # --- lightkurve (prioritaire en auto pour TESS/Kepler/K2) ---
            lightkurve_attempted = False
            if engine == "lightkurve" and not LIGHTKURVE_AVAILABLE:
                messagebox.showerror(
                    "lightkurve manquant",
                    "Installez lightkurve pour utiliser ce moteur :\n  pip install lightkurve",
                )
                return
            if engine in ("auto", "lightkurve") and LIGHTKURVE_AVAILABLE and self._lightkurve_supported_mission(mission):
                lightkurve_attempted = True
                ok = False
                # TIC : d'abord par coordonnées ; un timeout ne doit pas empêcher l'essai par nom.
                if self._should_try_tic_catalog_path(target):
                    try:
                        coord = self._tic_resolve_to_skycoord(target, radius_deg, query_timeout_s)
                        if coord is not None:
                            try:
                                ok = self._download_with_lightkurve_coords(
                                    coord,
                                    radius_deg,
                                    mission,
                                    staging_dir,
                                    max_products,
                                    search_timeout_s=lk_search_timeout_s,
                                    download_timeout_s=600,
                                    tess_authors=tess_authors,
                                )
                            except TimeoutError as te:
                                self._mast_dialog(f"lightkurve (coords): timeout — {te} ; essai par nom…")
                                logger.warning("lightkurve coords: %s", te)
                    except Exception as tic_lk_err:
                        self._mast_dialog(f"lightkurve (coords): erreur — {tic_lk_err}")
                        logger.exception(tic_lk_err)
                        if engine == "lightkurve":
                            raise
                if not ok:
                    try:
                        ok = self._download_with_lightkurve(
                            target,
                            mission,
                            staging_dir,
                            max_products,
                            search_timeout_s=lk_search_timeout_s,
                            download_timeout_s=600,
                            tess_authors=tess_authors,
                        )
                    except TimeoutError as te:
                        self._mast_dialog(f"lightkurve (nom): timeout — {te}")
                        logger.warning("lightkurve string search: %s", te)
                    except Exception as lk_err:
                        self._mast_dialog(f"lightkurve (nom): erreur — {lk_err}")
                        logger.exception(lk_err)
                        if engine == "lightkurve":
                            raise
                if ok:
                    self._mast_flatten_fits_to_folder(staging_dir, dest_flat)
                    after_fits = {p.resolve() for p in output_dir.rglob("*") if p.is_file() and self._mast_is_fits_product(p)}
                    new_fits = after_fits - before_fits
                    self.log_message(
                        f"🧾 Fichiers FITS réellement écrits sur disque: {len(new_fits)} (total détecté: {len(after_fits)})",
                        "info",
                    )
                    self.log_message(f"✅ Téléchargement lightkurve terminé -> {dest_flat}", "info")
                    messagebox.showinfo(
                        "Succès (lightkurve)",
                        f"Téléchargement via lightkurve terminé.\n"
                        f"Cible: {target}\n"
                        f"Fichiers FITS réellement écrits: {len(new_fits)}\n"
                        f"FITS détectés sous le répertoire de sortie: {len(after_fits)}\n"
                        f"Dossier plat (FITS) :\n{dest_flat}",
                    )
                    return

            if engine == "lightkurve":
                shutil.rmtree(staging_dir, ignore_errors=True)
                messagebox.showwarning(
                    "lightkurve",
                    "Aucune courbe téléchargée avec lightkurve.\n"
                    "Vérifiez la cible ou essayez le moteur « astroquery ».",
                )
                return

            # Mode auto : si lightkurve a été essayé sans succès, ne pas enchaîner sur le Portal
            # (chez vous query_region + query_object = 2×120 s de timeout inutiles).
            if (
                engine == "auto"
                and lightkurve_attempted
                and self._lightkurve_supported_mission(mission)
            ):
                self._mast_dialog(
                    "Arrêt (moteur auto) : lightkurve n'a rien téléchargé (0 résultat ou timeout). "
                    "Les services MAST derrière lightkurve / astroquery sont lents ou filtrés — pas de nouvelle tentative Portal automatique."
                )
                messagebox.showwarning(
                    "MAST — lightkurve sans résultat",
                    "lightkurve n'a pas abouti (aucun fichier ou recherche expirée).\n\n"
                    "Les requêtes MAST (lightkurve comme astroquery) peuvent rester sans réponse selon le réseau.\n"
                    "En mode « auto », on ne relance pas le Portal astroquery (~4 min de timeouts).\n\n"
                    "Options :\n"
                    "• Réessayer plus tard, autre réseau ou VPN.\n"
                    "• Moteur « astroquery » pour forcer query_region / query_object.\n"
                    "• Téléchargement manuel sur https://archive.stsci.edu (TESS).",
                )
                shutil.rmtree(staging_dir, ignore_errors=True)
                return

            # --- astroquery Observations ---
            self._mast_reset_staging_dir(staging_dir)
            if (
                not LIGHTKURVE_AVAILABLE
                and self._should_try_tic_catalog_path(target)
                and engine in ("auto", "astroquery")
            ):
                self._mast_dialog(
                    "⚠ lightkurve=False : Observations.query_region / query_object utilisent l’API Portal MAST "
                    "et peuvent rester bloqués >120 s (mast.stsci.edu peut répondre en HTTP sans que le Portal réponde). "
                    "Installez lightkurve avec le MÊME interpréteur que NPOAP, ex. : py -3 -m pip install lightkurve"
                )
            self.log_message(f"🔎 Recherche MAST (astroquery) pour cible: {target}", "info")
            self.log_message(f"📏 Rayon: {radius_deg} deg", "info")
            if mission:
                self.log_message(f"🛰️ Mission filtrée: {mission}", "info")

            observations = None
            if engine in ("auto", "astroquery") and self._should_try_tic_catalog_path(target):
                try:
                    observations = self._observations_via_tic_coords(target, radius_deg, query_timeout_s)
                    if observations is not None and len(observations) > 0:
                        self._mast_dialog(f"query_region: {len(observations)} observation(s).")
                except Exception as tic_ex:
                    self._mast_dialog(f"Voie TIC+query_region échouée: {tic_ex}")

            if observations is None or len(observations) == 0:
                self._mast_dialog("Appel MAST: Observations.query_object(...)")
                observations = self._call_with_timeout(
                    Observations.query_object,
                    query_timeout_s,
                    "query_object",
                    target,
                    radius=f"{radius_deg} deg",
                )
            if observations is None or len(observations) == 0:
                self._mast_dialog("Réponse MAST: 0 observation retournée.")
                self.log_message("⚠️ Aucune observation trouvée.", "warning")
                messagebox.showwarning("Aucun résultat", "Aucune observation trouvée sur MAST.")
                return

            self._mast_dialog(f"Réponse MAST: {len(observations)} observations reçues.")
            if "obs_collection" in observations.colnames:
                collections = sorted({str(x) for x in observations["obs_collection"]})
                self._mast_dialog(f"Collections observées: {', '.join(collections[:8])}{' ...' if len(collections) > 8 else ''}")
            if "target_name" in observations.colnames:
                targets = sorted({str(x) for x in observations["target_name"]})
                self._mast_dialog(f"Cibles résolues: {', '.join(targets[:6])}{' ...' if len(targets) > 6 else ''}")

            self.log_message(f"✅ Observations trouvées: {len(observations)}", "info")

            if mission and "obs_collection" in observations.colnames:
                before = len(observations)
                observations = observations[observations["obs_collection"] == mission]
                self._mast_dialog(f"Filtre mission '{mission}': {before} -> {len(observations)} observations")
                self.log_message(f"✅ Après filtre mission: {len(observations)} observations", "info")

            if len(observations) == 0:
                self._mast_dialog("Aucune observation après filtre mission.")
                self.log_message("⚠️ Aucune observation après filtre mission.", "warning")
                messagebox.showwarning("Aucun résultat", "Aucune observation après filtre mission.")
                return

            self._mast_dialog("Appel MAST: Observations.get_product_list(...)")
            products = self._call_with_timeout(
                Observations.get_product_list,
                products_timeout_s,
                "get_product_list",
                observations,
            )
            if products is None or len(products) == 0:
                self._mast_dialog("Réponse produits: 0 produit.")
                self.log_message("⚠️ Aucun produit disponible.", "warning")
                messagebox.showwarning("Aucun produit", "Aucun produit de données trouvé.")
                return

            self._mast_dialog(f"Produits bruts reçus: {len(products)}")
            if "productType" in products.colnames:
                ptypes = sorted({str(x) for x in products["productType"]})
                self._mast_dialog(f"Types de produits: {', '.join(ptypes[:8])}{' ...' if len(ptypes) > 8 else ''}")
            if "dataproduct_type" in products.colnames:
                dtypes = sorted({str(x) for x in products["dataproduct_type"]})
                self._mast_dialog(f"Data product types: {', '.join(dtypes[:8])}{' ...' if len(dtypes) > 8 else ''}")

            self._mast_dialog("Filtrage local: productType=SCIENCE, extension=fits|fits.gz")
            filtered = Observations.filter_products(
                products,
                productType=["SCIENCE"],
                extension=["fits", "fits.gz"],
            )
            if filtered is None or len(filtered) == 0:
                self._mast_dialog("Aucun produit après filtrage SCIENCE/FITS.")
                self.log_message("⚠️ Aucun produit SCIENCE FITS trouvé.", "warning")
                messagebox.showwarning("Aucun produit", "Aucun produit SCIENCE FITS trouvé.")
                return

            # Prioriser les séries temporelles quand possible
            if "dataproduct_type" in filtered.colnames:
                timeseries = filtered[filtered["dataproduct_type"] == "timeseries"]
                if len(timeseries) > 0:
                    self._mast_dialog(f"Priorisation timeseries: {len(filtered)} -> {len(timeseries)}")
                    filtered = timeseries
                    self.log_message(f"✅ Produits 'timeseries': {len(filtered)}", "info")
                else:
                    self._mast_dialog("Aucun dataproduct_type='timeseries' trouvé, conservation du filtrage courant.")

            if max_products > 0 and len(filtered) > max_products:
                self._mast_dialog(f"Limitation max_products={max_products}: {len(filtered)} -> {max_products}")
                filtered = filtered[:max_products]
                self.log_message(f"ℹ️ Limitation à {len(filtered)} produits", "info")

            self._mast_dialog("Début téléchargement: Observations.download_products(...)")
            self._mast_dialog(f"Dossier staging: {staging_dir} → puis aplatissement vers {dest_flat}")
            self.log_message(f"📥 Téléchargement de {len(filtered)} fichiers...", "info")
            manifest = self._call_with_timeout(
                Observations.download_products,
                download_timeout_s,
                "download_products",
                filtered,
                download_dir=str(staging_dir),
                cache=True,
            )

            downloaded = 0
            if manifest is not None and "Status" in manifest.colnames:
                downloaded = sum(1 for s in manifest["Status"] if str(s).upper() == "COMPLETE")
                self._mast_dialog(f"Manifest reçu: {len(manifest)} lignes, COMPLETE={downloaded}")
                if "Local Path" in manifest.colnames:
                    local_paths = [str(p) for p in manifest["Local Path"][:5]]
                    self._mast_dialog("Exemples de fichiers:")
                    for p in local_paths:
                        self._mast_dialog(f"  - {p}")
            else:
                self._mast_dialog("Manifest absent ou sans colonne Status.")

            self._mast_flatten_fits_to_folder(staging_dir, dest_flat)
            after_fits = {p.resolve() for p in output_dir.rglob("*") if p.is_file() and self._mast_is_fits_product(p)}
            new_fits = after_fits - before_fits
            self._mast_dialog(
                f"Comptage disque: nouveaux FITS/FITS.GZ={len(new_fits)}, total sous la sortie={len(after_fits)}"
            )
            self.log_message(
                f"🧾 Fichiers FITS réellement écrits sur disque: {len(new_fits)} (total détecté: {len(after_fits)})",
                "info",
            )
            self._mast_dialog("Session terminée avec succès.")
            self.log_message(f"✅ Téléchargement terminé: {downloaded}/{len(filtered)} complets", "info")
            messagebox.showinfo(
                "Succès",
                f"Téléchargement MAST terminé.\n"
                f"Cible: {target}\n"
                f"Produits demandés: {len(filtered)}\n"
                f"Téléchargements complets: {downloaded}\n"
                f"Fichiers FITS réellement écrits: {len(new_fits)}\n"
                f"FITS détectés sous le répertoire de sortie: {len(after_fits)}\n"
                f"Dossier plat (FITS) :\n{dest_flat}"
            )

        except TimeoutError as e:
            self._mast_dialog(f"ERREUR TIMEOUT: {e}")
            self.log_message(f"❌ Timeout MAST: {e}", "error")
            logger.error("Timeout MAST (astroquery Portal): %s", e)
            messagebox.showerror(
                "Timeout MAST (astroquery Portal)",
                f"{e}\n\n"
                "L’API Portal MAST (Observations.query_object / query_region) n’a pas répondu à temps.\n"
                "LcTools/LcGenerator utilise souvent des téléchargements directs depuis les archives.\n\n"
                "Essayez :\n"
                "  pip install lightkurve\n"
                "Puis dans NPOAP : moteur « auto » ou « lightkurve ».",
            )
        except Exception as e:
            self._mast_dialog(f"ERREUR: {e}")
            self.log_message(f"❌ Erreur MAST: {e}", "error")
            logger.exception(e)
            messagebox.showerror("Erreur", f"Erreur lors du téléchargement MAST:\n{e}")
    
    def select_gaia_file(self):
        """Sélectionne un fichier CSV ou CSV.GZ Gaia DR3 à analyser (séparation linéaire)."""
        file_path = filedialog.askopenfilename(
            title="Sélectionner un fichier Gaia DR3 (CSV ou CSV.GZ)",
            filetypes=[("CSV / CSV.GZ", "*.csv;*.csv.gz"), ("CSV", "*.csv"), ("CSV compressé", "*.csv.gz"), ("Tous", "*.*")]
        )
        if file_path:
            self.gaia_file_var.set(file_path)
    
    def select_nina_csv_file(self):
        """Sélectionne un fichier CSV à convertir en JSON NINA."""
        file_path = filedialog.askopenfilename(
            title="Sélectionner un fichier CSV d'étoiles binaires",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if file_path:
            self.nina_csv_file_var.set(file_path)
    
    def start_nina_conversion_thread(self):
        """Lance la conversion CSV vers JSON NINA dans un thread séparé."""
        thread = threading.Thread(target=self.run_nina_conversion, daemon=True)
        thread.start()
    
    def run_nina_conversion(self):
        """Exécute la conversion CSV vers JSON NINA."""
        try:
            csv_file = self.nina_csv_file_var.get().strip()
            if not csv_file:
                self.log_message("❌ Veuillez sélectionner un fichier CSV", "error")
                messagebox.showerror("Erreur", "Veuillez sélectionner un fichier CSV à convertir")
                return
            
            csv_path = Path(csv_file)
            if not csv_path.exists():
                self.log_message(f"❌ Fichier non trouvé : {csv_file}", "error")
                messagebox.showerror("Erreur", f"Fichier non trouvé : {csv_file}")
                return
            
            output_dir = Path(self.nina_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            self.log_message(f"📤 Conversion CSV vers JSON NINA...", "info")
            self.log_message(f"📁 Fichier CSV : {csv_path.name}", "info")
            self.log_message(f"📁 Répertoire de sortie : {output_dir}", "info")
            
            # Chemin du script Python
            script_path = Path(__file__).parent.parent / "utils" / "convert_binaries_to_nina.py"
            if not script_path.exists():
                error_msg = f"Script de conversion non trouvé : {script_path}"
                self.log_message(f"❌ {error_msg}", "error")
                messagebox.showerror("Erreur", error_msg)
                return
            
            # Exécuter le script Python
            python_exe = sys.executable
            cmd = [
                python_exe,
                str(script_path),
                str(csv_path),
                str(output_dir)
            ]
            
            self.log_message(f"🔧 Commande : {' '.join(cmd)}", "info")
            
            # Exécuter dans un processus
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding='utf-8',
                errors='replace'
            )
            
            # Lire la sortie en temps réel
            output_lines = []
            for line in process.stdout:
                line = line.strip()
                if line:
                    output_lines.append(line)
                    self.log_message(line, "info")
            
            # Attendre la fin du processus
            return_code = process.wait()
            
            if return_code == 0:
                # Compter les fichiers JSON créés
                json_files = list(output_dir.glob("*.json"))
                count = len(json_files)
                
                self.log_message(f"✅ Conversion terminée avec succès !", "info")
                self.log_message(f"📊 {count} fichiers JSON créés", "info")
                
                messagebox.showinfo(
                    "Succès",
                    f"Conversion terminée avec succès !\n\n"
                    f"{count} fichiers JSON créés\n"
                    f"\nFichiers dans :\n{output_dir}"
                )
            else:
                error_msg = f"Erreur lors de la conversion (code {return_code})"
                self.log_message(f"❌ {error_msg}", "error")
                if output_lines:
                    self.log_message("Détails de l'erreur :", "error")
                    for line in output_lines[-10:]:  # Afficher les 10 dernières lignes
                        self.log_message(f"  {line}", "error")
                messagebox.showerror("Erreur", error_msg)
                
        except Exception as e:
            error_msg = f"Erreur lors de la conversion : {e}"
            logger.error(error_msg, exc_info=True)
            self.log_message(f"❌ {error_msg}", "error")
            messagebox.showerror("Erreur", error_msg)
    
    def start_linear_separation_analysis_thread(self):
        """Lance l'analyse de séparation linéaire dans un thread séparé."""
        thread = threading.Thread(target=self.run_linear_separation_analysis, daemon=True)
        thread.start()
    
    def run_linear_separation_analysis(self):
        """Exécute l'analyse de séparation linéaire selon la méthode de Laurent (2022)."""
        try:
            from core.linear_separation_calculator import LinearSeparationCalculator
            import pandas as pd
            
            gaia_file = self.gaia_file_var.get().strip()
            if not gaia_file:
                self.log_message("❌ Veuillez sélectionner un fichier CSV Gaia DR3", "error")
                messagebox.showerror("Erreur", "Veuillez sélectionner un fichier CSV Gaia DR3")
                return
            
            gaia_path = Path(gaia_file)
            if not gaia_path.exists():
                self.log_message(f"❌ Fichier non trouvé : {gaia_file}", "error")
                messagebox.showerror("Erreur", f"Fichier non trouvé : {gaia_file}")
                return
            
            # Paramètres
            threshold_pc = 10.0
            try:
                threshold_pc = float(self.linear_sep_threshold_var.get())
            except ValueError:
                self.log_message("⚠️ Seuil invalide, utilisation de 10.0 pc par défaut", "warning")
            
            max_angular_sep = 60.0
            try:
                max_angular_sep = float(self.max_angular_sep_var.get())
            except ValueError:
                self.log_message("⚠️ Séparation angulaire max invalide, utilisation de 60.0 arcsec par défaut", "warning")
            
            self.log_message(f"📚 Analyse de séparation linéaire (seuil: {threshold_pc} pc, sép. ang. max: {max_angular_sep} arcsec)...", "info")
            self.log_message(f"📁 Fichier Gaia : {gaia_path.name}", "info")
            
            # Initialiser le calculateur
            calculator = LinearSeparationCalculator()
            
            # Analyser le fichier Gaia
            try:
                all_pairs, physical_pairs = calculator.analyze_gaia_csv_file(
                    gaia_path,
                    max_angular_separation_arcsec=max_angular_sep,
                    threshold_pc=threshold_pc,
                )
            except Exception as e:
                error_msg = f"Erreur lors de l'analyse du fichier Gaia : {e}"
                logger.error(error_msg, exc_info=True)
                self.log_message(f"❌ {error_msg}", "error")
                messagebox.showerror("Erreur", error_msg)
                return
            
            if not all_pairs:
                self.log_message("⚠️ Aucun couple n'a été trouvé dans le fichier", "warning")
                messagebox.showwarning(
                    "Attention",
                    "Aucun couple n'a été trouvé dans le fichier.\n\n"
                    "Conseils :\n"
                    "• Augmentez la « Séparation angulaire max » (ex. 120 arcsec)\n"
                    "• Vérifiez que le CSV contient les colonnes ra, dec, parallax (ou RA_ICRS, DE_ICRS, Plx)\n"
                    "• Le fichier doit contenir au moins 2 étoiles avec parallaxe > 0"
                )
                return
            
            # Créer des DataFrames
            results_df = pd.DataFrame(all_pairs)
            physical_df = pd.DataFrame(physical_pairs) if physical_pairs else pd.DataFrame()
            
            # Réorganiser les colonnes : ra1, dec1, ra2, dec2 en premier
            if not results_df.empty:
                # Colonnes prioritaires
                priority_cols = ['ra1', 'dec1', 'ra2', 'dec2']
                # Autres colonnes dans l'ordre d'apparition
                other_cols = [col for col in results_df.columns if col not in priority_cols]
                # Nouvel ordre des colonnes
                new_column_order = priority_cols + other_cols
                # Réorganiser (garder seulement les colonnes qui existent)
                available_cols = [col for col in new_column_order if col in results_df.columns]
                results_df = results_df[available_cols]
            
            if not physical_df.empty:
                # Même réorganisation pour le DataFrame des couples physiques
                priority_cols = ['ra1', 'dec1', 'ra2', 'dec2']
                other_cols = [col for col in physical_df.columns if col not in priority_cols]
                new_column_order = priority_cols + other_cols
                available_cols = [col for col in new_column_order if col in physical_df.columns]
                physical_df = physical_df[available_cols]
            
            self.log_message(f"✅ {len(all_pairs)} couples analysés avec succès", "info")
            self.log_message(f"✅ {len(physical_pairs)} couples physiques (SL < {threshold_pc} pc)", "info")
            
            # Sauvegarder les résultats
            output_dir = Path(self.binaries_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Fichier complet (tous les couples)
            output_file_all = output_dir / f"linear_separation_all_{timestamp}.csv"
            results_df.to_csv(output_file_all, index=False)
            self.log_message(f"💾 Résultats complets : {output_file_all.name}", "info")
            
            # Fichier results.csv : couples physiques uniquement (SL < seuil pc)
            results_csv = output_dir / "results.csv"
            if len(physical_pairs) > 0:
                physical_df.to_csv(results_csv, index=False)
                self.log_message(f"💾 Couples physiques (SL < {threshold_pc} pc) : {results_csv.name}", "info")
            else:
                # Fichier vide ou avec en-têtes seulement si aucun couple physique
                physical_df.to_csv(results_csv, index=False)
                self.log_message(f"💾 Aucun couple physique ; {results_csv.name} créé (vide)", "info")
            
            messagebox.showinfo(
                "Succès",
                f"Analyse terminée !\n\n"
                f"{len(all_pairs)} couples analysés avec succès\n"
                f"{len(physical_pairs)} couples physiques (SL < {threshold_pc} pc)\n"
                f"\nCouples physiques écrits dans :\n{results_csv}\n"
                f"Répertoire : {output_dir}"
            )
            
        except ImportError as e:
            self.log_message(f"❌ Module requis manquant : {e}", "error")
            messagebox.showerror("Erreur", f"Module requis manquant : {e}\n\nInstallez : pip install pandas")
        except Exception as e:
            error_msg = f"Erreur lors de l'analyse : {e}"
            logger.error(error_msg, exc_info=True)
            self.log_message(f"❌ {error_msg}", "error")
            messagebox.showerror("Erreur", error_msg)
            logger.exception(e)
            messagebox.showerror("Erreur", f"Erreur lors de l'analyse : {e}")
    
    # ============================================================
    # MÉTHODES POUR ONGLET ÉTOILES
    # ============================================================
    
    def init_stars_catalogs(self):
        """Initialise la liste des catalogues disponibles pour les étoiles."""
        if not CATALOG_EXTRACTOR_AVAILABLE:
            self.log_message("⚠️ Module d'extraction non disponible. Installez astroquery.", "error")
            return
        
        # Filtrer les catalogues pertinents pour les étoiles
        star_catalogs = []
        for name, info in POPULAR_CATALOGS.items():
            if "Étoiles" in info.get("object_types", []) or name in ["Gaia DR3", "Gaia EDR3", "Hipparcos", "Tycho-2", "2MASS", "USNO-B1.0", "SDSS DR16", "WISE", "TIC-8.2"]:
                star_catalogs.append(name)
        
        self.stars_catalog_combo['values'] = star_catalogs
        if star_catalogs:
            self.stars_catalog_combo.current(0)
            self.on_stars_catalog_selected()
    
    def on_stars_catalog_selected(self, event=None):
        """Appelé quand un catalogue est sélectionné dans l'onglet Étoiles."""
        catalog_name = self.stars_catalog_var.get()
        if catalog_name in POPULAR_CATALOGS:
            catalog_info = POPULAR_CATALOGS[catalog_name]
            self.stars_catalog_desc_label.config(text=catalog_info.get("description", ""))
    
    def start_stars_extraction_thread(self):
        """Lance l'extraction dans un thread séparé pour l'onglet Étoiles."""
        if not CATALOG_EXTRACTOR_AVAILABLE:
            messagebox.showerror("Erreur", "Module d'extraction non disponible.")
            return
        
        thread = threading.Thread(target=self.run_stars_extraction, daemon=True)
        thread.start()
    
    def run_stars_extraction(self):
        """Exécute l'extraction des données pour l'onglet Étoiles."""
        try:
            catalog_name = self.stars_catalog_var.get()
            if not catalog_name:
                self.log_message("❌ Veuillez sélectionner un catalogue.", "error")
                return
            
            output_dir = Path(self.stars_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            self.extractor = CatalogExtractor(output_dir=output_dir)
            self.log_message(f"📚 Extraction depuis : {catalog_name}", "info")

            mag_min = None
            mag_max = None
            if self.stars_mag_min_var.get().strip():
                try:
                    mag_min = float(self.stars_mag_min_var.get())
                except ValueError:
                    self.log_message("⚠️ Magnitude min invalide, ignorée", "warning")

            if self.stars_mag_max_var.get().strip():
                try:
                    mag_max = float(self.stars_mag_max_var.get())
                except ValueError:
                    self.log_message("⚠️ Magnitude max invalide, ignorée", "warning")

            ra_str = self.stars_ra_center_var.get().strip()
            dec_str = self.stars_dec_center_var.get().strip()
            radius_str = self.stars_radius_var.get().strip()

            if not ra_str or not dec_str or not radius_str:
                self.log_message("❌ Veuillez remplir tous les champs (RA, DEC, rayon)", "error")
                return

            try:
                center_coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
            except Exception:
                center_coord = SkyCoord(ra=float(ra_str) * u.deg, dec=float(dec_str) * u.deg)

            radius = float(radius_str) * u.deg

            self.log_message(f"📍 Centre : {center_coord}", "info")
            self.log_message(f"📏 Rayon : {radius}", "info")

            table = self.extractor.extract_by_region(
                catalog_name=catalog_name,
                center_coord=center_coord,
                radius=radius,
                mag_limit=mag_max,
                mag_column=None,
            )
            if mag_min is not None and len(table) > 0 and catalog_name in POPULAR_CATALOGS:
                mag_col = POPULAR_CATALOGS[catalog_name].get("mag_column")
                if mag_col and mag_col in table.colnames:
                    table = table[table[mag_col] >= mag_min]
                    self.log_message(f"🔽 Après filtre magnitude ≥ {mag_min} : {len(table)} objet(s)", "info")

            # Sauvegarder les résultats
            if len(table) == 0:
                self.log_message("⚠️ Aucun résultat trouvé.", "warning")
                return
            
            self.log_message(f"✅ {len(table)} objets extraits", "info")
            
            # Réduire aux colonnes utiles (ra, dec, parallax, source_id, pmra, pmdec, phot_g_mean_mag)
            table = self.extractor.reduce_to_essential_columns(table)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{catalog_name.replace(' ', '_')}_{timestamp}"
            
            output_format = self.stars_output_format_var.get()
            output_path = self.extractor.save_table(table, filename, format=output_format)
            
            self.log_message(f"💾 Fichier sauvegardé : {output_path}", "info")
            self.log_message(f"📊 Colonnes : {', '.join(table.colnames)}", "info")
            
            messagebox.showinfo("Succès", f"Extraction terminée !\n{len(table)} objets extraits\nFichier : {output_path}")
            
        except Exception as e:
            self.log_message(f"❌ Erreur fatale : {e}", "error")
            logger.exception(e)
            messagebox.showerror("Erreur", f"Erreur lors de l'extraction : {e}")
    
    # ============================================================
    # MÉTHODES POUR EXTRACTION GAIA DR3
    # ============================================================
    
    def start_gaia_download_thread(self):
        """Lance le téléchargement Gaia DR3 dans un thread séparé."""
        thread = threading.Thread(target=self.run_gaia_download, daemon=True)
        thread.start()
    
    def run_gaia_download(self):
        """Exécute le téléchargement des catalogues Gaia DR3."""
        try:
            from core.gaia_dr3_extractor import GaiaDR3Extractor
            
            output_dir = Path(self.gaia_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            mag_limit = float(self.gaia_mag_var.get()) if self.gaia_mag_var.get() else 15.0
            ra_step_deg = float(self.gaia_ra_step_var.get()) if self.gaia_ra_step_var.get() else 30.0
            
            ra_min_h = None
            ra_max_h = None
            if self.gaia_ra_min_h_var.get().strip():
                ra_min_h = float(self.gaia_ra_min_h_var.get())
            if self.gaia_ra_max_h_var.get().strip():
                ra_max_h = float(self.gaia_ra_max_h_var.get())
            
            # Récupérer les limites de déclinaison
            dec_min = None
            dec_max = None
            if self.gaia_dec_min_var.get().strip():
                dec_min = float(self.gaia_dec_min_var.get())
            if self.gaia_dec_max_var.get().strip():
                dec_max = float(self.gaia_dec_max_var.get())
            
            hemisphere = self.gaia_hemisphere_var.get()
            skip_existing = self.gaia_skip_existing_var.get()
            filter_variables = self.gaia_filter_variables_var.get()
            filter_galaxies = self.gaia_filter_galaxies_var.get()
            
            self.log_message(f"📥 Début téléchargement catalogues Gaia DR3...", "info")
            if filter_variables:
                self.log_message(f"   Filtre : Étoiles variables uniquement", "info")
            if filter_galaxies:
                self.log_message(f"   Filtre : Galaxies uniquement", "info")
            if filter_variables and filter_galaxies:
                self.log_message(f"   ⚠️ Attention : Les deux filtres sont activés (peut donner peu de résultats)", "warning")
            self.log_message(f"   Répertoire: {output_dir}", "info")
            self.log_message(f"   Magnitude max G: {mag_limit}", "info")
            self.log_message(f"   Pas RA: {ra_step_deg}°", "info")
            if ra_min_h is not None or ra_max_h is not None:
                self.log_message(f"   RA: {ra_min_h or 0}h - {ra_max_h or 24}h", "info")
            if dec_min is not None or dec_max is not None:
                self.log_message(f"   DEC: {dec_min or -90}° - {dec_max or 90}°", "info")
            self.log_message(f"   Hémisphère: {hemisphere}", "info")
            self.log_message(f"   Ignorer fichiers existants: {skip_existing}", "info")
            
            # Compter les fichiers existants
            existing_files = list(output_dir.glob("gaia_dr3_*.csv.gz")) + list(output_dir.glob("gaia_dr3_*.csv"))
            if existing_files and skip_existing:
                self.log_message(f"   {len(existing_files)} fichiers existants détectés (seront ignorés)", "info")
            
            extractor = GaiaDR3Extractor(output_dir=output_dir)
            
            if hemisphere == "both":
                # Convertir ra_min_h et ra_max_h en degrés
                ra_min_deg = None
                ra_max_deg = None
                if ra_min_h is not None:
                    ra_min_deg = ra_min_h * 15.0  # Convertir heures en degrés
                if ra_max_h is not None:
                    ra_max_deg = ra_max_h * 15.0  # Convertir heures en degrés
                
                results, north_count, south_count = extractor.extract_both_hemispheres(
                    mag_limit=mag_limit,
                    ra_step_deg=ra_step_deg,
                    output_format="csv.gz",
                    skip_existing=skip_existing,
                    ra_min_deg=ra_min_deg,
                    ra_max_deg=ra_max_deg,
                    dec_min=dec_min,
                    dec_max=dec_max,
                    filter_variables=filter_variables,
                    filter_galaxies=filter_galaxies
                )
                self.log_message(f"✅ Téléchargement terminé: {north_count + south_count} étoiles au total", "info")
                self.log_message(f"   Hémisphère nord: {north_count} étoiles", "info")
                self.log_message(f"   Hémisphère sud: {south_count} étoiles", "info")
            else:
                # Convertir ra_min_h et ra_max_h en degrés
                ra_min_deg = None
                ra_max_deg = None
                if ra_min_h is not None:
                    ra_min_deg = ra_min_h * 15.0  # Convertir heures en degrés
                if ra_max_h is not None:
                    ra_max_deg = ra_max_h * 15.0  # Convertir heures en degrés
                
                files, count = extractor.extract_hemisphere_catalog(
                    hemisphere=hemisphere,
                    mag_limit=mag_limit,
                    ra_step_deg=ra_step_deg,
                    output_format="csv.gz",
                    skip_existing=skip_existing,
                    ra_min_deg=ra_min_deg,
                    ra_max_deg=ra_max_deg,
                    dec_min=dec_min,
                    dec_max=dec_max,
                    filter_variables=filter_variables,
                    filter_galaxies=filter_galaxies
                )
                self.log_message(f"✅ Téléchargement terminé: {count} étoiles", "info")
            
            self.log_message(f"✅ Téléchargement catalogues Gaia DR3 terminé", "info")
            messagebox.showinfo("Succès", "Téléchargement des catalogues Gaia DR3 terminé !")
            
        except Exception as e:
            error_msg = f"Erreur lors du téléchargement Gaia DR3: {e}"
            self.log_message(f"❌ {error_msg}", "error")
            logger.exception(e)
            messagebox.showerror("Erreur", error_msg)