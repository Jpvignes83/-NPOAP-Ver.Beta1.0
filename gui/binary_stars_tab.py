# gui/binary_stars_tab.py
"""
Onglet pour la modélisation d'étoiles binaires à éclipses avec PHOEBE2
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import logging
import numpy as np
from pathlib import Path
import threading

from gui.manual_help import add_manual_help_header

logger = logging.getLogger(__name__)

# Import de la fenêtre de visualisation
try:
    from gui.binary_stars_viewer import BinaryStarsViewer
    VIEWER_AVAILABLE = True
except ImportError:
    VIEWER_AVAILABLE = False

# Périodogramme (analyse de données)
PERIODOGRAM_AVAILABLE = False
run_lomb_scargle = run_bls = run_plavchan = None
try:
    from core.periodogram_tools import run_lomb_scargle, run_bls, run_plavchan
    PERIODOGRAM_AVAILABLE = True
    logger.info("[Périodogramme] Module core.periodogram_tools chargé.")
except ImportError as e:
    logger.warning("[Périodogramme] Import core.periodogram_tools échoué: %s", e, exc_info=True)
    import sys
    from pathlib import Path as _Path
    _parent = _Path(__file__).resolve().parent.parent
    if str(_parent) not in sys.path:
        sys.path.insert(0, str(_parent))
    try:
        from core.periodogram_tools import run_lomb_scargle, run_bls, run_plavchan
        PERIODOGRAM_AVAILABLE = True
        logger.info("[Périodogramme] Module chargé après ajout de %s à sys.path.", _parent)
    except ImportError as e2:
        logger.warning("[Périodogramme] Import échoué après fallback: %s", e2, exc_info=True)

# Import PHOEBE2 avec gestion d'erreur si non installé
try:
    # Patch pour Windows : créer un module readline factice si nécessaire
    import sys
    if sys.platform == 'win32' and 'readline' not in sys.modules:
        class MockReadline:
            def add_history(self, *args, **kwargs): pass
            def read_history_file(self, *args, **kwargs): pass
            def write_history_file(self, *args, **kwargs): pass
            def set_history_length(self, *args, **kwargs): pass
            def get_history_length(self, *args, **kwargs): return 0
            def set_completer(self, *args, **kwargs): pass
            def set_completer_delims(self, *args, **kwargs): pass
            def set_completion_display_matches_hook(self, *args, **kwargs): pass
            def parse_and_bind(self, *args, **kwargs): pass
        sys.modules['readline'] = MockReadline()
    
    import phoebe
    PHOEBE_AVAILABLE = True
    logger.debug(f"[PHOEBE] Import réussi (readline mock actif sur Windows si besoin)")
except (ImportError, Exception) as e:
    PHOEBE_AVAILABLE = False
    logger.warning(f"PHOEBE2 n'est pas installé ou erreur d'import: {e}. Utilisez: pip install phoebe", exc_info=True)


class BinaryStarsTab(ttk.Frame):
    """
    Onglet pour la modélisation d'étoiles binaires avec PHOEBE2
    """
    
    def __init__(self, parent_notebook, base_dir=None):
        super().__init__(parent_notebook, padding=10)
        logger.debug("[PHOEBE] BinaryStarsTab __init__, base_dir=%s", base_dir)
        if base_dir is None:
            self.base_dir = Path.home()
        else:
            self.base_dir = Path(base_dir)
        self.b = None
        self.lc_time = None
        self.lc_flux = None
        self.lc_error = None
        self.create_widgets()
        if not PHOEBE_AVAILABLE:
            logger.warning("[PHOEBE] PHOEBE2 non disponible à l'initialisation")
            self.show_phoebe_install_message()
        logger.info("[PHOEBE] Onglet modélisation binaires initialisé (périodogramme=%s)", PERIODOGRAM_AVAILABLE)
    
    def create_widgets(self):
        """Crée l'interface utilisateur"""
        add_manual_help_header(self, "8-étoiles-binaires")

        # En-tête
        header_frame = ttk.Frame(self)
        header_frame.pack(fill="x", pady=(0, 10))
        
        title_label = ttk.Label(
            header_frame,
            text="Modélisation d'Étoiles Binaires à Éclipses",
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
        
        # 1. Statut PHOEBE2 (ligne compacte)
        status_frame = ttk.Frame(parent)
        status_frame.pack(fill="x", pady=(0, 2))
        
        self.phoebe_status_label = ttk.Label(
            status_frame,
            text="État: Vérification...",
            foreground="blue"
        )
        self.phoebe_status_label.pack(anchor="w", padx=2, pady=1)
        
        if PHOEBE_AVAILABLE:
            try:
                version = phoebe.__version__
                self.phoebe_status_label.config(
                    text=f"✓ PHOEBE2 {version} installé",
                    foreground="green"
                )
            except:
                self.phoebe_status_label.config(
                    text="✓ PHOEBE2 installé",
                    foreground="green"
                )
        else:
            self.phoebe_status_label.config(
                text="✗ PHOEBE2 non installé",
                foreground="red"
            )
        
        # 2. Charger des données observées
        data_frame = ttk.LabelFrame(parent, text="2. Charger des données", padding=6)
        data_frame.pack(fill="x", pady=3)
        
        self.data_file_var = tk.StringVar()
        
        file_frame = ttk.Frame(data_frame)
        file_frame.pack(fill="x", pady=1)
        ttk.Label(file_frame, text="Fichier (CSV):").pack(side="left", padx=(0, 5))
        ttk.Entry(file_frame, textvariable=self.data_file_var, width=30).pack(side="left", fill="x", expand=True)
        ttk.Button(file_frame, text="Parcourir", command=self.browse_data_file).pack(side="left", padx=(5, 0))
        
        ttk.Label(data_frame, text="Format: time(s), flux(es), error/ferrs (optionnel). Temps en BJD-TDB ou JD.", font=("Helvetica", 8), foreground="gray").pack(anchor="w", pady=1)
        
        ttk.Button(
            data_frame,
            text="Charger les données",
            command=self.load_observed_data
        ).pack(pady=2)
        
        # 4. Périodogramme (estimer la période à partir des données chargées)
        perio_frame = ttk.LabelFrame(parent, text="4. Périodogramme", padding=6)
        perio_frame.pack(fill="x", pady=3)
        ttk.Label(perio_frame, text="Utilise les données chargées ci-dessus. Min/Max période (jours):", font=("Helvetica", 8), foreground="gray").pack(anchor="w", pady=(0, 2))
        perio_params = ttk.Frame(perio_frame)
        perio_params.pack(fill="x", pady=1)
        ttk.Label(perio_params, text="Min P (j):").pack(side="left", padx=(0, 2))
        self.perio_min_p_var = tk.StringVar(value="0.5")
        ttk.Entry(perio_params, textvariable=self.perio_min_p_var, width=6).pack(side="left", padx=(0, 8))
        ttk.Label(perio_params, text="Max P (j):").pack(side="left", padx=(0, 2))
        self.perio_max_p_var = tk.StringVar(value="10.0")
        ttk.Entry(perio_params, textvariable=self.perio_max_p_var, width=6).pack(side="left", padx=(0, 8))
        perio_btns = ttk.Frame(perio_frame)
        perio_btns.pack(fill="x", pady=2)
        ttk.Button(perio_btns, text="Lomb-Scargle", command=lambda: self._run_periodogram("LS")).pack(side="left", padx=2)
        ttk.Button(perio_btns, text="BLS (Transit)", command=lambda: self._run_periodogram("BLS")).pack(side="left", padx=2)
        ttk.Button(perio_btns, text="Plavchan", command=lambda: self._run_periodogram("PLAV")).pack(side="left", padx=2)
        
        # 5. Créer un nouveau système (bundle PHOEBE pour modéliser les données chargées)
        bundle_frame = ttk.LabelFrame(parent, text="5. Créer un nouveau système", padding=6)
        bundle_frame.pack(fill="x", pady=3)
        
        ttk.Label(bundle_frame, text="Type de système:").pack(anchor="w", pady=2)
        
        self.system_type_var = tk.StringVar(value="binary")
        system_types = [
            ("Binaire à éclipses", "binary"),
            ("Système de contact", "contact"),
        ]
        
        for text, value in system_types:
            ttk.Radiobutton(
                bundle_frame,
                text=text,
                variable=self.system_type_var,
                value=value
            ).pack(anchor="w", padx=20)
        
        ttk.Button(
            bundle_frame,
            text="Créer un nouveau Bundle",
            command=self.create_bundle
        ).pack(pady=5)
        
        # 5. Paramètres du modèle
        params_frame = ttk.LabelFrame(parent, text="5. Paramètres du modèle", padding=10)
        params_frame.pack(fill="x", pady=5)
        
        # Paramètres de base
        self.period_var = tk.StringVar(value="1.0")
        self.t0_var = tk.StringVar(value="0.0")
        self.incl_var = tk.StringVar(value="90.0")
        
        param_grid = [
            ("Période (jours):", self.period_var),
            ("t0 (JD):", self.t0_var),
            ("Inclinaison (°):", self.incl_var),
        ]
        
        for label, var in param_grid:
            frame = ttk.Frame(params_frame)
            frame.pack(fill="x", pady=1)
            ttk.Label(frame, text=label, width=15).pack(side="left")
            ttk.Entry(frame, textvariable=var, width=12).pack(side="left", padx=3)
        
        # 6. Calcul et visualisation
        compute_frame = ttk.LabelFrame(parent, text="6. Calcul", padding=6)
        compute_frame.pack(fill="x", pady=3)
        
        ttk.Button(
            compute_frame,
            text="Calculer le modèle",
            command=self.compute_model
        ).pack(pady=2)
        
        ttk.Button(
            compute_frame,
            text="Ajuster les paramètres",
            command=self.fit_parameters
        ).pack(pady=2)
        
        ttk.Button(
            compute_frame,
            text="🎬 Visualisation 3D",
            command=self.open_viewer
        ).pack(pady=2)
        
        # Barre de progression
        self.progress = ttk.Progressbar(compute_frame, mode="indeterminate")
        self.progress.pack(fill="x", pady=2)
        
        # 7. Sauvegarder/Charger
        save_frame = ttk.LabelFrame(parent, text="7. Sauvegarder/Charger", padding=6)
        save_frame.pack(fill="x", pady=3)
        
        ttk.Button(
            save_frame,
            text="Sauvegarder le Bundle",
            command=self.save_bundle
        ).pack(pady=2)
        
        ttk.Button(
            save_frame,
            text="Charger un Bundle",
            command=self.load_bundle
        ).pack(pady=2)
    
    def create_visualization_section(self, parent):
        """Crée la section de visualisation à droite"""
        
        # Zone pour les graphiques matplotlib : ordre vertical = courbe de lumière, périodogramme, orbite (en bas)
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.gridspec import GridSpec
        
        self.fig = plt.figure(figsize=(8, 11))
        gs = GridSpec(3, 1, figure=self.fig, height_ratios=[2, 1.2, 2], hspace=0.35)
        self.ax_lc = self.fig.add_subplot(gs[0, 0])   # en haut : courbe de lumière
        self.ax_perio = self.fig.add_subplot(gs[1, 0])  # au milieu : périodogramme
        self.ax_orb = self.fig.add_subplot(gs[2, 0])    # en bas : orbite
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        toolbar = NavigationToolbar2Tk(self.canvas, parent)
        toolbar.update()
        
        # Initialisation des axes
        self.ax_lc.set_xlabel("Temps (JD)")
        self.ax_lc.set_ylabel("Flux normalisé")
        self.ax_lc.set_title("Courbe de lumière")
        self.ax_lc.grid(True, alpha=0.3)
        
        self.ax_perio.set_xlabel("Période (jours)")
        self.ax_perio.set_ylabel("Puissance")
        self.ax_perio.set_title("Périodogramme")
        self.ax_perio.grid(True, linestyle=':', alpha=0.6)
        
        self.ax_orb.set_xlabel("X (R☉)")
        self.ax_orb.set_ylabel("Y (R☉)")
        self.ax_orb.set_title("Orbite")
        self.ax_orb.set_aspect('equal')
        self.ax_orb.grid(True, alpha=0.3)
    
    def show_phoebe_install_message(self):
        """Affiche un message si PHOEBE2 n'est pas installé"""
        messagebox.showwarning(
            "PHOEBE2 non installé",
            "PHOEBE2 n'est pas installé.\n\n"
            "Pour l'installer, exécutez dans un terminal:\n"
            "pip install phoebe\n\n"
            "L'onglet fonctionnera de manière limitée jusqu'à l'installation."
        )
    
    def install_phoebe(self):
        """Tente d'installer PHOEBE2"""
        logger.info("[PHOEBE] install_phoebe demandé")
        import subprocess
        import sys
        
        def install():
            try:
                self.progress.start()
                logger.debug("[PHOEBE] Lancement pip install phoebe")
                result = subprocess.run(
                    [sys.executable, "-m", "pip", "install", "phoebe"],
                    capture_output=True,
                    text=True
                )
                self.progress.stop()
                if result.returncode == 0:
                    logger.info("[PHOEBE] pip install phoebe réussi")
                    messagebox.showinfo("Succès", "PHOEBE2 a été installé avec succès!\nVeuillez redémarrer l'application.")
                else:
                    logger.warning("[PHOEBE] pip install phoebe échec: %s", result.stderr)
                    messagebox.showerror("Erreur", f"Échec de l'installation:\n{result.stderr}")
            except Exception as e:
                self.progress.stop()
                logger.error("[PHOEBE] Erreur install_phoebe: %s", e, exc_info=True)
                messagebox.showerror("Erreur", f"Erreur lors de l'installation: {e}")
        threading.Thread(target=install, daemon=True).start()
    
    def create_bundle(self):
        """Crée un nouveau Bundle PHOEBE2"""
        logger.debug("[PHOEBE] create_bundle appelé")
        if not PHOEBE_AVAILABLE:
            logger.warning("[PHOEBE] Bundle non créé: PHOEBE2 non disponible")
            messagebox.showerror("Erreur", "PHOEBE2 n'est pas installé!")
            return
        
        try:
            system_type = self.system_type_var.get()
            logger.info(f"[PHOEBE] Création du bundle, type demandé: {system_type}")
            
            if system_type == "binary":
                logger.debug("[PHOEBE] Appel phoebe.default_binary()")
                self.b = phoebe.default_binary()
            elif system_type == "contact":
                logger.debug("[PHOEBE] Appel phoebe.default_binary(contact_binary=True)")
                self.b = phoebe.default_binary(contact_binary=True)
            else:
                logger.debug("[PHOEBE] Appel phoebe.default_binary() (fallback)")
                self.b = phoebe.default_binary()
            
            logger.info(f"[PHOEBE] Bundle créé avec succès: {system_type}")
            messagebox.showinfo("Succès", f"Bundle créé ({system_type})")
            
        except Exception as e:
            logger.error(f"[PHOEBE] Erreur lors de la création du Bundle: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de la création du Bundle: {e}")
    
    def browse_data_file(self):
        """Ouvre un dialogue pour sélectionner un fichier de données"""
        logger.debug("[PHOEBE] browse_data_file, base_dir=%s", self.base_dir)
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("CSV", "*.csv"), ("Tous les fichiers", "*.*")]
        )
        if path:
            self.data_file_var.set(path)
            logger.info("[PHOEBE] Fichier sélectionné: %s", path)
    
    def load_observed_data(self):
        """Charge des données observées depuis un fichier CSV"""
        file_path = self.data_file_var.get()
        logger.debug(f"[PHOEBE] load_observed_data appelé, fichier: {file_path!r}")
        if not file_path:
            messagebox.showwarning("Attention", "Sélectionnez d'abord un fichier!")
            return
        
        try:
            import pandas as pd
            
            logger.info(f"[PHOEBE] Chargement CSV: {file_path}")
            df = pd.read_csv(file_path)
            logger.debug(f"[PHOEBE] Colonnes détectées: {list(df.columns)}")
            
            # Noms de colonnes acceptés (time/times, flux/fluxes, error/ferrs)
            time_col = 'time' if 'time' in df.columns else ('times' if 'times' in df.columns else None)
            flux_col = 'flux' if 'flux' in df.columns else ('fluxes' if 'fluxes' in df.columns else None)
            if time_col is None or flux_col is None:
                messagebox.showerror(
                    "Erreur",
                    "Le fichier doit contenir des colonnes de temps et de flux.\n"
                    "Accepté: time ou times, flux ou fluxes.\n"
                    f"Colonnes trouvées: {', '.join(df.columns)}"
                )
                return
            
            self.lc_time = df[time_col].values.astype(float)
            self.lc_flux = df[flux_col].values.astype(float)
            
            err_col = 'error' if 'error' in df.columns else ('ferrs' if 'ferrs' in df.columns else None)
            if err_col is not None:
                self.lc_error = df[err_col].values.astype(float)
            else:
                self.lc_error = None
            
            # Afficher les données
            self.ax_lc.clear()
            if self.lc_error is not None:
                self.ax_lc.errorbar(
                    self.lc_time, self.lc_flux, yerr=self.lc_error,
                    fmt='o', markersize=3, alpha=0.6, label="Observations"
                )
            else:
                self.ax_lc.plot(
                    self.lc_time, self.lc_flux, 'o', markersize=3, alpha=0.6, label="Observations"
                )
            self.ax_lc.set_xlabel("Temps (JD)")
            self.ax_lc.set_ylabel("Flux normalisé")
            self.ax_lc.set_title("Courbe de lumière observée")
            self.ax_lc.legend()
            self.ax_lc.grid(True, alpha=0.3)
            self.canvas.draw()
            
            messagebox.showinfo("Succès", f"Données chargées: {len(self.lc_time)} points")
            logger.info(f"Données chargées depuis {file_path}: {len(self.lc_time)} points")
            
        except Exception as e:
            logger.error(f"[PHOEBE] Erreur chargement données: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors du chargement: {e}")
    
    def _run_periodogram(self, algo):
        """Lance l'algorithme de périodogramme (Lomb-Scargle, BLS ou Plavchan) sur les données chargées."""
        logger.info("[Périodogramme] _run_periodogram algo=%s", algo)
        if not PERIODOGRAM_AVAILABLE or run_lomb_scargle is None:
            logger.error("[Périodogramme] Module non disponible (PERIODOGRAM_AVAILABLE=%s, run_lomb_scargle=%s)",
                         PERIODOGRAM_AVAILABLE, run_lomb_scargle is not None)
            messagebox.showerror(
                "Erreur",
                "Module périodogramme non disponible (core.periodogram_tools).\nConsultez les logs pour la cause."
            )
            return
        if self.lc_time is None or self.lc_flux is None:
            logger.warning("[Périodogramme] Données non chargées (lc_time=%s, lc_flux=%s)",
                           self.lc_time is not None, self.lc_flux is not None)
            messagebox.showwarning("Attention", "Chargez d'abord des données (étape 2).")
            return
        try:
            min_p = float(self.perio_min_p_var.get())
            max_p = float(self.perio_max_p_var.get())
            if min_p >= max_p:
                raise ValueError("Min P doit être < Max P")
        except Exception as e:
            logger.warning("[Périodogramme] Paramètres invalides: %s", e)
            messagebox.showerror("Erreur", "Vérifiez Min P et Max P (nombres, Min < Max).")
            return
        logger.debug("[Périodogramme] Lancement thread min_p=%s max_p=%s", min_p, max_p)
        threading.Thread(target=self._worker_periodogram, args=(algo, min_p, max_p), daemon=True).start()
    
    def _worker_periodogram(self, algo, min_p, max_p):
        """Exécute le calcul en arrière-plan et met à jour l'interface via after()."""
        res = None
        name = "Périodogramme"
        logger.info("[Périodogramme] _worker_periodogram algo=%s min_p=%s max_p=%s n_points=%s",
                    algo, min_p, max_p, len(self.lc_time) if self.lc_time is not None else 0)
        try:
            if algo == "LS":
                res = run_lomb_scargle(self.lc_time, self.lc_flux, min_period=min_p, max_period=max_p)
                name = "Lomb-Scargle"
            elif algo == "BLS":
                res = run_bls(self.lc_time, self.lc_flux, min_period=min_p, max_period=max_p)
                name = "BLS (Transit)"
            elif algo == "PLAV":
                res = run_plavchan(self.lc_time, self.lc_flux, min_period=min_p, max_period=max_p)
                name = "Plavchan"
            if res is not None:
                logger.info("[Périodogramme] Calcul terminé %s, meilleure période=%.6f j", name, res[2])
                self.after(0, lambda r=res, n=name: self._on_periodogram_done(r[0], r[1], r[2], n))
        except Exception as e:
            logger.error("[Périodogramme] Erreur calcul: %s", e, exc_info=True)
            self.after(0, lambda msg=str(e): messagebox.showerror("Erreur calcul", msg))
    
    def _on_periodogram_done(self, periods, powers, best_period, name):
        """Affiche le périodogramme dans le 3e graphique et propose de copier la période."""
        logger.info("[Périodogramme] _on_periodogram_done %s best=%.6f j", name, best_period)
        self.ax_perio.clear()
        self.ax_perio.plot(periods, powers, color="black", linewidth=1)
        self.ax_perio.set_xlabel("Période (jours)")
        self.ax_perio.set_ylabel("Puissance")
        self.ax_perio.set_title(f"{name} — Meilleur pic: {best_period:.5f} j")
        self.ax_perio.axvline(best_period, color="red", linestyle="--", linewidth=1.5)
        self.ax_perio.grid(True, linestyle=':', alpha=0.6)
        self.canvas.draw()
        self.period_var.set(f"{best_period:.6f}")
        logger.debug("[Périodogramme] Période mise dans champ Période (jours)")
        messagebox.showinfo("Périodogramme", f"Période détectée : {best_period:.6f} j\nElle a été mise dans « Période (jours) ».")
    
    def compute_model(self):
        """Calcule le modèle PHOEBE2"""
        logger.debug("[PHOEBE] compute_model appelé")
        if not PHOEBE_AVAILABLE:
            messagebox.showerror("Erreur", "PHOEBE2 n'est pas installé!")
            return
        
        if self.b is None:
            logger.warning("[PHOEBE] compute_model: pas de bundle (créer d'abord un Bundle)")
            messagebox.showwarning("Attention", "Créez d'abord un Bundle!")
            return
        
        try:
            # Mise à jour des paramètres depuis l'interface
            period = float(self.period_var.get())
            t0 = float(self.t0_var.get())
            incl = float(self.incl_var.get())
            logger.info(f"[PHOEBE] Paramètres: period={period}, t0={t0}, incl={incl}")
            
            self.b['period@binary'] = period
            self.b['t0@system'] = t0
            self.b['incl@binary'] = incl
            
            # Vérifier si un dataset existe, sinon en créer un
            if len(self.b.datasets) == 0:
                # Créer un dataset de courbe de lumière
                if self.lc_time is not None:
                    # Utiliser les temps observés
                    logger.debug(f"[PHOEBE] Ajout dataset lc01 avec {len(self.lc_time)} points observés")
                    self.b.add_dataset('lc', times=self.lc_time, dataset='lc01')
                else:
                    # Créer une grille de temps par défaut
                    times = np.linspace(0, period, 100)
                    logger.debug("[PHOEBE] Ajout dataset lc01 avec grille par défaut")
                    self.b.add_dataset('lc', times=times, dataset='lc01')
            
            # Calcul du modèle
            self.progress.start()
            import time
            t0 = time.perf_counter()
            logger.info("[PHOEBE] Lancement run_compute(phoebe01)...")
            self.b.run_compute(compute='phoebe01')
            elapsed = time.perf_counter() - t0
            self.progress.stop()
            logger.info("[PHOEBE] run_compute terminé en %.2f s", elapsed)
            
            # Récupération des résultats
            dataset_labels = self.b.datasets
            if len(dataset_labels) > 0:
                dataset_label = dataset_labels[0]
                try:
                    times = self.b['value@times@{}@model'.format(dataset_label)]
                    fluxes = self.b['value@fluxes@{}@model'.format(dataset_label)]
                except:
                    # Fallback sur lc01 si le format est différent
                    times = self.b['value@times@lc01@model']
                    fluxes = self.b['value@fluxes@lc01@model']
            else:
                raise ValueError("Aucun dataset trouvé après le calcul")
            
            # Normaliser le flux du modèle (PHOEBE renvoie un flux physique) pour comparer aux observations (flux normalisé ~1)
            flux_max = np.nanmax(fluxes)
            if flux_max > 0:
                fluxes = fluxes / flux_max
            
            # Visualisation
            self.ax_lc.clear()
            
            # Afficher les observations si disponibles
            if self.lc_time is not None:
                if self.lc_error is not None:
                    self.ax_lc.errorbar(
                        self.lc_time, self.lc_flux, yerr=self.lc_error,
                        fmt='o', markersize=3, alpha=0.6, label="Observations"
                    )
                else:
                    self.ax_lc.plot(
                        self.lc_time, self.lc_flux, 'o', markersize=3, alpha=0.6, label="Observations"
                    )
            
            # Afficher le modèle (déjà normalisé)
            self.ax_lc.plot(times, fluxes, '-', linewidth=2, label="Modèle")
            self.ax_lc.set_xlabel("Temps (JD)")
            self.ax_lc.set_ylabel("Flux normalisé")
            self.ax_lc.set_title("Courbe de lumière: Observations vs Modèle")
            self.ax_lc.legend()
            self.ax_lc.grid(True, alpha=0.3)
            
            self.canvas.draw()
            messagebox.showinfo("Succès", "Modèle calculé avec succès!\n\nTemps de calcul : %.2f s" % elapsed)
            logger.info("Modèle PHOEBE2 calculé (%.2f s)", elapsed)
            
        except Exception as e:
            self.progress.stop()
            messagebox.showerror("Erreur", f"Erreur lors du calcul: {e}")
            logger.error(f"Erreur calcul modèle: {e}", exc_info=True)
    
    def fit_parameters(self):
        """Ajuste période, t0 et inclinaison sur la courbe de lumière observée (minimisation chi²)."""
        logger.info("[PHOEBE] fit_parameters demandé")
        if not PHOEBE_AVAILABLE:
            messagebox.showerror("Erreur", "PHOEBE2 n'est pas installé!")
            return
        if self.b is None:
            messagebox.showwarning("Attention", "Créez d'abord un Bundle (étape 4).")
            return
        if self.lc_time is None or self.lc_flux is None:
            messagebox.showwarning("Attention", "Chargez des données observées (étape 2).")
            return
        period0 = float(self.period_var.get())
        t0_0 = float(self.t0_var.get())
        incl0 = float(self.incl_var.get())
        n_pts = len(self.lc_time)
        # S'assurer que le dataset existe aux temps observés
        if len(self.b.datasets) == 0:
            self.b.add_dataset('lc', times=self.lc_time, dataset='lc01')
        self.progress.config(mode="determinate", maximum=100, value=0)
        self.progress_label = getattr(self, 'progress_label', None)
        if self.progress_label is not None:
            self.progress_label.config(text="Ajustement en cours...")
        self._fit_done = False
        self._fit_error = None
        self._fit_result = None

        def objective(x):
            period, t0, incl = x[0], x[1], x[2]
            try:
                self.b['period@binary'] = period
                self.b['t0@system'] = t0
                self.b['incl@binary'] = incl
                self.b.run_compute(compute='phoebe01')
                fmod = np.asarray(self.b['value@fluxes@lc01@model'], dtype=float)
            except Exception:
                return 1e30
            fmax = np.nanmax(fmod)
            if fmax <= 0:
                return 1e30
            fmod_norm = fmod / fmax
            if self.lc_error is not None and np.all(self.lc_error > 0):
                chi2 = np.sum(((self.lc_flux - fmod_norm) / self.lc_error) ** 2)
            else:
                chi2 = np.sum((self.lc_flux - fmod_norm) ** 2)
            return chi2

        def run_fit():
            from scipy.optimize import minimize
            try:
                t_min, t_max = np.min(self.lc_time), np.max(self.lc_time)
                bounds = [
                    (max(0.1, period0 * 0.5), min(100, period0 * 2)),
                    (t_min - 2 * period0, t_max + 2 * period0),
                    (50.0, 90.0),
                ]
                res = minimize(
                    objective,
                    x0=[period0, t0_0, incl0],
                    method='L-BFGS-B',
                    bounds=bounds,
                    options=dict(maxfun=80, maxiter=40),
                )
                if res.success:
                    self._fit_result = (res.x[0], res.x[1], res.x[2], res.fun)
                else:
                    self._fit_result = (res.x[0], res.x[1], res.x[2], res.fun)
            except Exception as e:
                logger.error("[PHOEBE] Erreur fit_parameters: %s", e, exc_info=True)
                self._fit_error = e
            self._fit_done = True

        def on_fit_done():
            self.progress.config(value=100)
            if getattr(self, 'progress_label', None) is not None:
                self.progress_label.config(text="Ajustement terminé.")
                self.progress_label.after(2000, lambda: self.progress_label.config(text=""))
            self.progress.config(mode="indeterminate")
            if self._fit_error is not None:
                messagebox.showerror("Erreur", "Erreur lors de l'ajustement: %s" % self._fit_error)
                return
            p, t0, incl, chi2 = self._fit_result
            self.period_var.set("%.6f" % p)
            self.t0_var.set("%.6f" % t0)
            self.incl_var.set("%.2f" % incl)
            logger.info("[PHOEBE] Ajustement: period=%.6f t0=%.6f incl=%.2f chi2=%.2f", p, t0, incl, chi2)
            messagebox.showinfo(
                "Ajustement terminé",
                "Paramètres ajustés:\n\nPériode = %.6f j\nt0 = %.6f\nInclinaison = %.2f°\n\nChi² = %.2f\n\nCliquez sur « Calculer le modèle » pour afficher la courbe."
                % (p, t0, incl, chi2)
            )
            self.compute_model()

        def poll():
            if self._fit_done:
                self.after(0, on_fit_done)
                return
            self.after(200, poll)

        threading.Thread(target=run_fit, daemon=True).start()
        self.after(200, poll)
    
    def save_bundle(self):
        """Sauvegarde le Bundle PHOEBE2"""
        logger.debug("[PHOEBE] save_bundle")
        if not PHOEBE_AVAILABLE:
            messagebox.showerror("Erreur", "PHOEBE2 n'est pas installé!")
            return
        if self.b is None:
            logger.warning("[PHOEBE] save_bundle: pas de bundle")
            messagebox.showwarning("Attention", "Créez d'abord un Bundle!")
            return
        path = filedialog.asksaveasfilename(
            initialdir=self.base_dir,
            defaultextension=".phoebe",
            filetypes=[("PHOEBE", "*.phoebe"), ("Tous les fichiers", "*.*")]
        )
        if path:
            try:
                self.b.save(path)
                logger.info("[PHOEBE] Bundle sauvegardé: %s", path)
                messagebox.showinfo("Succès", f"Bundle sauvegardé: {path}")
            except Exception as e:
                logger.error("[PHOEBE] Erreur sauvegarde: %s", e, exc_info=True)
                messagebox.showerror("Erreur", f"Erreur lors de la sauvegarde: {e}")
    
    def load_bundle(self):
        """Charge un Bundle PHOEBE2"""
        logger.debug("[PHOEBE] load_bundle")
        if not PHOEBE_AVAILABLE:
            messagebox.showerror("Erreur", "PHOEBE2 n'est pas installé!")
            return
        path = filedialog.askopenfilename(
            initialdir=self.base_dir,
            filetypes=[("PHOEBE", "*.phoebe"), ("Tous les fichiers", "*.*")]
        )
        if path:
            try:
                self.b = phoebe.open(path)
                logger.info("[PHOEBE] Bundle chargé: %s", path)
                messagebox.showinfo("Succès", f"Bundle chargé: {path}")
            except Exception as e:
                logger.error("[PHOEBE] Erreur chargement bundle: %s", e, exc_info=True)
                messagebox.showerror("Erreur", f"Erreur lors du chargement: {e}")
    
    def open_viewer(self):
        """Ouvre la fenêtre de visualisation 3D"""
        logger.debug("[PHOEBE] open_viewer")
        if not VIEWER_AVAILABLE:
            logger.warning("[PHOEBE] Viewer non disponible")
            messagebox.showerror(
                "Erreur",
                "Le module de visualisation n'est pas disponible.\n"
                "Assurez-vous que gui/binary_stars_viewer.py existe."
            )
            return
        try:
            viewer = BinaryStarsViewer(self, bundle=self.b)
            viewer.protocol("WM_DELETE_WINDOW", viewer.on_close)
            logger.info("[PHOEBE] Fenêtre de visualisation 3D ouverte")
        except Exception as e:
            logger.error("[PHOEBE] Erreur ouverture viewer: %s", e, exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de l'ouverture de la visualisation: {e}")

