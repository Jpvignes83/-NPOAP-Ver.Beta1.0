# ============================================================
# CODE À AJOUTER À gui/catalogues_tab.py
# ============================================================
# 
# Ce fichier contient uniquement le code à ajouter/modifier.
# Instructions :
# 1. Remplacer la méthode create_binary_stars_tab (lignes ~478-491)
# 2. Remplacer la méthode create_exoplanets_tab (lignes ~493-506)  
# 3. Ajouter toutes les méthodes ci-dessous à la fin du fichier (après run_stars_extraction)
# ============================================================

# ============================================================
# MÉTHODE create_binary_stars_tab (REMPLACER)
# ============================================================

def create_binary_stars_tab_REPLACE(self):
    """Crée l'onglet pour l'extraction de catalogues d'étoiles binaires."""
    binaries_frame = ttk.Frame(self.inner_notebook, padding=10)
    self.inner_notebook.add(binaries_frame, text="⭐ Étoiles Binaires")
    
    header = ttk.Label(
        binaries_frame,
        text="Extraction depuis catalogues d'étoiles binaires",
        font=("Helvetica", 10, "bold")
    )
    header.pack(pady=5)
    
    # Frame principal
    main_container = ttk.Frame(binaries_frame)
    main_container.pack(fill=tk.BOTH, expand=True)
    
    # Colonne gauche : Catalogues
    left_frame = ttk.LabelFrame(main_container, text="Catalogues disponibles", padding=10)
    left_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=(0, 10), pady=5)
    
    # AAVSO
    aavso_frame = ttk.LabelFrame(left_frame, text="AAVSO Target Tool", padding=5)
    aavso_frame.pack(fill="x", pady=5)
    ttk.Label(aavso_frame, text="Étoiles binaires à éclipses depuis AAVSO", 
             font=("Helvetica", 8), foreground="gray").pack(anchor="w")
    ttk.Button(aavso_frame, text="📥 Télécharger AAVSO (Binaires)", 
              command=self.download_aavso_binaries).pack(pady=5)
    self.aavso_binaries_status = ttk.Label(aavso_frame, text="Non téléchargé", foreground="gray")
    self.aavso_binaries_status.pack()
    
    # VSX (Vizier)
    vsx_frame = ttk.LabelFrame(left_frame, text="VSX (Variable Star Index)", padding=5)
    vsx_frame.pack(fill="x", pady=5)
    ttk.Label(vsx_frame, text="Via Vizier - Extraction par région", 
             font=("Helvetica", 8), foreground="gray").pack(anchor="w")
    ttk.Button(vsx_frame, text="📥 Extraire VSX", 
              command=self.start_vsx_extraction_thread).pack(pady=5)
    
    # SB9 (Vizier)
    sb9_frame = ttk.LabelFrame(left_frame, text="SB9 (Spectroscopic Binaries)", padding=5)
    sb9_frame.pack(fill="x", pady=5)
    ttk.Label(sb9_frame, text="Via Vizier - Extraction par région", 
             font=("Helvetica", 8), foreground="gray").pack(anchor="w")
    ttk.Button(sb9_frame, text="📥 Extraire SB9", 
              command=self.start_sb9_extraction_thread).pack(pady=5)
    
    # TESS EBS
    tess_frame = ttk.LabelFrame(left_frame, text="TESS EBS", padding=5)
    tess_frame.pack(fill="x", pady=5)
    ttk.Label(tess_frame, text="TESS Eclipsing Binary Stars Catalog", 
             font=("Helvetica", 8), foreground="gray").pack(anchor="w")
    ttk.Button(tess_frame, text="📥 Télécharger TESS EBS", 
              command=self.start_tess_ebs_extraction_thread).pack(pady=5)
    
    # Options de sortie
    output_frame = ttk.LabelFrame(left_frame, text="Options de sortie", padding=5)
    output_frame.pack(fill="x", pady=5)
    
    ttk.Label(output_frame, text="Répertoire :").pack(anchor="w")
    output_path_frame = ttk.Frame(output_frame)
    output_path_frame.pack(fill="x", pady=2)
    self.binaries_output_dir_var = tk.StringVar(value=str(self.output_dir))
    ttk.Entry(output_path_frame, textvariable=self.binaries_output_dir_var, width=25).pack(side="left", fill="x", expand=True)
    ttk.Button(output_path_frame, text="📁", command=lambda: self.select_output_directory(self.binaries_output_dir_var), width=3).pack(side="left", padx=2)
    
    # Mettre à jour le statut AAVSO
    self.update_aavso_binaries_status()


# ============================================================
# MÉTHODE create_exoplanets_tab (REMPLACER)
# ============================================================

def create_exoplanets_tab_REPLACE(self):
    """Crée l'onglet pour l'extraction de catalogues d'exoplanètes."""
    exoplanets_frame = ttk.Frame(self.inner_notebook, padding=10)
    self.inner_notebook.add(exoplanets_frame, text="🔭 Exoplanètes")
    
    header = ttk.Label(
        exoplanets_frame,
        text="Extraction depuis catalogues d'exoplanètes",
        font=("Helvetica", 10, "bold")
    )
    header.pack(pady=5)
    
    # Frame principal
    main_container = ttk.Frame(exoplanets_frame)
    main_container.pack(fill=tk.BOTH, expand=True)
    
    # Colonne gauche : Catalogues
    left_frame = ttk.LabelFrame(main_container, text="Catalogues disponibles", padding=10)
    left_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=(0, 10), pady=5)
    
    # AAVSO
    aavso_frame = ttk.LabelFrame(left_frame, text="AAVSO Target Tool", padding=5)
    aavso_frame.pack(fill="x", pady=5)
    ttk.Label(aavso_frame, text="Exoplanètes depuis AAVSO", 
             font=("Helvetica", 8), foreground="gray").pack(anchor="w")
    ttk.Button(aavso_frame, text="📥 Télécharger AAVSO (Exoplanètes)", 
              command=self.download_aavso_exoplanets).pack(pady=5)
    self.aavso_exoplanets_status = ttk.Label(aavso_frame, text="Non téléchargé", foreground="gray")
    self.aavso_exoplanets_status.pack()
    
    # exoplanet.eu
    exoplanet_eu_frame = ttk.LabelFrame(left_frame, text="exoplanet.eu", padding=5)
    exoplanet_eu_frame.pack(fill="x", pady=5)
    ttk.Label(exoplanet_eu_frame, text="Catalogue complet des exoplanètes confirmées", 
             font=("Helvetica", 8), foreground="gray").pack(anchor="w")
    
    format_frame = ttk.Frame(exoplanet_eu_frame)
    format_frame.pack(fill="x", pady=2)
    self.exoplanet_eu_format_var = tk.StringVar(value="csv")
    ttk.Radiobutton(format_frame, text="CSV", variable=self.exoplanet_eu_format_var, value="csv").pack(side="left", padx=5)
    ttk.Radiobutton(format_frame, text="VOTable", variable=self.exoplanet_eu_format_var, value="votable").pack(side="left", padx=5)
    
    ttk.Button(exoplanet_eu_frame, text="📥 Télécharger exoplanet.eu", 
              command=self.start_exoplanet_eu_extraction_thread).pack(pady=5)
    
    # Options de sortie
    output_frame = ttk.LabelFrame(left_frame, text="Options de sortie", padding=5)
    output_frame.pack(fill="x", pady=5)
    
    ttk.Label(output_frame, text="Répertoire :").pack(anchor="w")
    output_path_frame = ttk.Frame(output_frame)
    output_path_frame.pack(fill="x", pady=2)
    self.exoplanets_output_dir_var = tk.StringVar(value=str(self.output_dir))
    ttk.Entry(output_path_frame, textvariable=self.exoplanets_output_dir_var, width=25).pack(side="left", fill="x", expand=True)
    ttk.Button(output_path_frame, text="📁", command=lambda: self.select_output_directory(self.exoplanets_output_dir_var), width=3).pack(side="left", padx=2)
    
    # Mettre à jour le statut AAVSO
    self.update_aavso_exoplanets_status()


# ============================================================
# MÉTHODES À AJOUTER À LA FIN DU FICHIER
# ============================================================
# Ajouter tout le code ci-dessous à la fin de catalogues_tab.py
# (après la méthode run_stars_extraction, vers la ligne 733)
# ============================================================

    # ============================================================
    # MÉTHODES POUR ONGLET ÉTOILES BINAIRES
    # ============================================================
    
    def download_aavso_binaries(self):
        """Télécharge le catalogue AAVSO pour les étoiles binaires à éclipses."""
        try:
            if not REQUESTS_AVAILABLE:
                raise RuntimeError("requests n'est pas disponible. Installez-le avec: pip install requests")
            
            self.aavso_binaries_status.config(text="Téléchargement en cours...", foreground="blue")
            self.update()
            
            from utils.fetch_aavso_targets import fetch_aavso_targets
            
            output_dir = Path(self.binaries_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = output_dir / f"aavso_binaries_{timestamp}.csv"
            
            self.log_message("Téléchargement du catalogue AAVSO (étoiles binaires)...", "info")
            result = fetch_aavso_targets(
                output_file=str(output_file),
                exoplanets=False,
                eclipsing_binaries=True
            )
            
            if result:
                self.log_message(f"✅ Catalogue AAVSO téléchargé : {result}", "info")
                self.update_aavso_binaries_status()
                messagebox.showinfo("Succès", f"Catalogue AAVSO (étoiles binaires) téléchargé\nFichier : {result}")
            else:
                raise RuntimeError("Échec du téléchargement")
                
        except Exception as e:
            self.log_message(f"❌ Erreur téléchargement AAVSO binaires: {e}", "error")
            logger.exception(e)
            self.aavso_binaries_status.config(text="Erreur", foreground="red")
            messagebox.showerror("Erreur", f"Erreur lors du téléchargement : {e}")
    
    def update_aavso_binaries_status(self):
        """Met à jour le statut du catalogue AAVSO binaires."""
        try:
            output_dir = Path(self.binaries_output_dir_var.get())
            aavso_files = list(output_dir.glob("aavso_binaries_*.csv"))
            if aavso_files:
                latest = max(aavso_files, key=lambda p: p.stat().st_mtime)
                date = datetime.fromtimestamp(latest.stat().st_mtime)
                self.aavso_binaries_status.config(
                    text=f"Téléchargé ({date.strftime('%Y-%m-%d %H:%M')})",
                    foreground="green"
                )
            else:
                self.aavso_binaries_status.config(text="Non téléchargé", foreground="gray")
        except Exception as e:
            logger.debug(f"Erreur mise à jour statut AAVSO binaires: {e}")
    
    def start_vsx_extraction_thread(self):
        """Lance l'extraction VSX dans un thread."""
        if not CATALOG_EXTRACTOR_AVAILABLE:
            messagebox.showerror("Erreur", "Module d'extraction non disponible.")
            return
        messagebox.showinfo("Info", "Pour extraire VSX, utilisez l'onglet '⭐ Étoiles' et sélectionnez le catalogue 'VSX'.")
    
    def start_sb9_extraction_thread(self):
        """Lance l'extraction SB9 dans un thread."""
        if not CATALOG_EXTRACTOR_AVAILABLE:
            messagebox.showerror("Erreur", "Module d'extraction non disponible.")
            return
        messagebox.showinfo("Info", "Pour extraire SB9, utilisez l'onglet '⭐ Étoiles' et sélectionnez le catalogue 'SB9'.")
    
    def start_tess_ebs_extraction_thread(self):
        """Lance l'extraction TESS EBS dans un thread."""
        try:
            if not CATALOG_EXTRACTOR_AVAILABLE:
                messagebox.showerror("Erreur", "Module d'extraction non disponible.")
                return
            
            thread = threading.Thread(target=self.run_tess_ebs_extraction, daemon=True)
            thread.start()
        except Exception as e:
            self.log_message(f"❌ Erreur : {e}", "error")
            logger.exception(e)
    
    def run_tess_ebs_extraction(self):
        """Exécute l'extraction TESS EBS."""
        try:
            output_dir = Path(self.binaries_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            self.extractor = CatalogExtractor(output_dir=output_dir)
            self.log_message("📚 Extraction TESS EBS...", "info")
            
            table = self.extractor.extract_tess_ebs()
            
            if len(table) == 0:
                self.log_message("⚠️ Aucun résultat trouvé.", "warning")
                return
            
            self.log_message(f"✅ {len(table)} objets extraits", "info")
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"TESS_EBS_{timestamp}"
            output_path = self.extractor.save_table(table, filename, format="csv")
            
            self.log_message(f"💾 Fichier sauvegardé : {output_path}", "info")
            messagebox.showinfo("Succès", f"Extraction terminée !\n{len(table)} objets extraits\nFichier : {output_path}")
            
        except Exception as e:
            self.log_message(f"❌ Erreur fatale : {e}", "error")
            logger.exception(e)
            messagebox.showerror("Erreur", f"Erreur lors de l'extraction : {e}")
    
    # ============================================================
    # MÉTHODES POUR ONGLET EXOPLANÈTES
    # ============================================================
    
    def download_aavso_exoplanets(self):
        """Télécharge le catalogue AAVSO pour les exoplanètes."""
        try:
            if not REQUESTS_AVAILABLE:
                raise RuntimeError("requests n'est pas disponible. Installez-le avec: pip install requests")
            
            self.aavso_exoplanets_status.config(text="Téléchargement en cours...", foreground="blue")
            self.update()
            
            from utils.fetch_aavso_targets import fetch_aavso_targets
            
            output_dir = Path(self.exoplanets_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = output_dir / f"aavso_exoplanets_{timestamp}.csv"
            
            self.log_message("Téléchargement du catalogue AAVSO (exoplanètes)...", "info")
            result = fetch_aavso_targets(
                output_file=str(output_file),
                exoplanets=True,
                eclipsing_binaries=False
            )
            
            if result:
                self.log_message(f"✅ Catalogue AAVSO téléchargé : {result}", "info")
                self.update_aavso_exoplanets_status()
                messagebox.showinfo("Succès", f"Catalogue AAVSO (exoplanètes) téléchargé\nFichier : {result}")
            else:
                raise RuntimeError("Échec du téléchargement")
                
        except Exception as e:
            self.log_message(f"❌ Erreur téléchargement AAVSO exoplanètes: {e}", "error")
            logger.exception(e)
            self.aavso_exoplanets_status.config(text="Erreur", foreground="red")
            messagebox.showerror("Erreur", f"Erreur lors du téléchargement : {e}")
    
    def update_aavso_exoplanets_status(self):
        """Met à jour le statut du catalogue AAVSO exoplanètes."""
        try:
            output_dir = Path(self.exoplanets_output_dir_var.get())
            aavso_files = list(output_dir.glob("aavso_exoplanets_*.csv"))
            if aavso_files:
                latest = max(aavso_files, key=lambda p: p.stat().st_mtime)
                date = datetime.fromtimestamp(latest.stat().st_mtime)
                self.aavso_exoplanets_status.config(
                    text=f"Téléchargé ({date.strftime('%Y-%m-%d %H:%M')})",
                    foreground="green"
                )
            else:
                self.aavso_exoplanets_status.config(text="Non téléchargé", foreground="gray")
        except Exception as e:
            logger.debug(f"Erreur mise à jour statut AAVSO exoplanètes: {e}")
    
    def start_exoplanet_eu_extraction_thread(self):
        """Lance l'extraction exoplanet.eu dans un thread."""
        try:
            if not CATALOG_EXTRACTOR_AVAILABLE:
                messagebox.showerror("Erreur", "Module d'extraction non disponible.")
                return
            
            thread = threading.Thread(target=self.run_exoplanet_eu_extraction, daemon=True)
            thread.start()
        except Exception as e:
            self.log_message(f"❌ Erreur : {e}", "error")
            logger.exception(e)
    
    def run_exoplanet_eu_extraction(self):
        """Exécute l'extraction exoplanet.eu."""
        try:
            output_dir = Path(self.exoplanets_output_dir_var.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            
            self.extractor = CatalogExtractor(output_dir=output_dir)
            self.log_message("📚 Extraction exoplanet.eu...", "info")
            
            format_type = self.exoplanet_eu_format_var.get()
            table = self.extractor.extract_exoplanet_eu(format=format_type)
            
            if len(table) == 0:
                self.log_message("⚠️ Aucun résultat trouvé.", "warning")
                return
            
            self.log_message(f"✅ {len(table)} exoplanètes extraites", "info")
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"exoplanet_eu_{timestamp}"
            output_path = self.extractor.save_table(table, filename, format=format_type)
            
            self.log_message(f"💾 Fichier sauvegardé : {output_path}", "info")
            messagebox.showinfo("Succès", f"Extraction terminée !\n{len(table)} exoplanètes extraites\nFichier : {output_path}")
            
        except Exception as e:
            self.log_message(f"❌ Erreur fatale : {e}", "error")
            logger.exception(e)
            messagebox.showerror("Erreur", f"Erreur lors de l'extraction : {e}")
