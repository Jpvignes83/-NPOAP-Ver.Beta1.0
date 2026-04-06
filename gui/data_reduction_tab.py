import logging
import threading
import queue
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, Toplevel
from pathlib import Path
import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.io import fits
from astropy.visualization import ZScaleInterval

from core.image_processor import ImageProcessor
from core.astrometry import AstrometrySolverNova

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class CCDProcGUI:
    def __init__(self, parent):
        self.frame = ttk.Frame(parent)
        self.frame.pack(fill=tk.BOTH, expand=True)

        self.processor = None
        self.files = []
        self.bias_files = []
        self.dark_files = []
        self.flat_files = []
        self.output_dir = None
        self.sep_catalog_folder = ""
        self.progress_prefix = "📊 Progression :"
        self.astrometry_thread = None
        self.astrometry_cancelled = False
        self.astrometry_process = None
        self.viewer_directory = None  # Répertoire pour la visualisation (section 5)
        self._ui_queue = queue.Queue()
        self.create_widgets()
        self._start_ui_queue_poller()

        localappdata = os.environ.get("LOCALAPPDATA")
        if not localappdata:
            localappdata = os.path.join(os.environ.get("USERPROFILE"), "Local Settings", "Application Data")
        self.bash_path = os.path.join(localappdata, "cygwin_ansvr", "bin", "bash.exe")

    # ------------------------------------------------------------------
    # GUI
    # ------------------------------------------------------------------
    def create_widgets(self):
        # --- Panneau Gauche (Contrôles) ---
        left_frame = ttk.Frame(self.frame)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        # --- Panneau Droit (Logs) ---
        right_frame = ttk.Frame(self.frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # ============================================================
        # SECTION 1 : FICHIERS & CALIBRATION
        # ============================================================
        calib_frame = ttk.LabelFrame(left_frame, text="1. Fichiers & Calibration")
        calib_frame.pack(fill="x", padx=5, pady=5)

        calib_buttons = [
            ("📁 Définir Répertoire", self.set_output_directory),
            ("📂 Charger Lights", self.load_files),
            ("📂 Charger Bias", self.load_bias),
            ("📂 Charger Darks", self.load_darks),
            ("📂 Charger Flats", self.load_flats),
        ]

        for text, command in calib_buttons:
            btn = ttk.Button(calib_frame, text=text, command=command, width=30)
            btn.pack(pady=2, padx=5, fill="x")

        # Option scaling darks (extrapolation temps d'exposition, type AstroImageJ)
        self.scale_darks_var = tk.BooleanVar(value=False)
        scale_darks_cb = ttk.Checkbutton(
            calib_frame,
            text="Scaler les darks au temps d'exposition des lights (si différent)",
            variable=self.scale_darks_var,
        )
        scale_darks_cb.pack(pady=4, padx=5, anchor="w")

        ttk.Button(
            calib_frame, text="🚀 Lancer Calibration", command=self.run_calibration, width=30
        ).pack(pady=2, padx=5, fill="x")

        # ============================================================
        # SECTION 2 : ASTROMÉTRIE
        # ============================================================
        astro_frame = ttk.LabelFrame(left_frame, text="2. Astrométrie (Plate Solving)")
        astro_frame.pack(fill="x", padx=5, pady=10)

        ttk.Button(
            astro_frame, 
            text="🌐 Via Astrometry.net (NOVA)", 
            command=self.run_astrometry_nova
        ).pack(pady=2, padx=5, fill="x")

        ttk.Button(
            astro_frame, 
            text="🖥️ Astrométrie Locale (WSL)", 
            command=self.run_astrometry_local
        ).pack(pady=2, padx=5, fill="x")
        
        # Répertoire catalogues pour LOCAL_G
        # catalog_frame = ttk.Frame(astro_frame)
        # catalog_frame.pack(fill="x", padx=5, pady=5)
        # ttk.Label(catalog_frame, text="Répertoire catalogues").pack(anchor="w")
        # catalog_path_frame = ttk.Frame(catalog_frame)
        # catalog_path_frame.pack(fill="x", pady=2)
        # self.catalog_dir_var = tk.StringVar()
        # ttk.Entry(catalog_path_frame, textvariable=self.catalog_dir_var, width=25).pack(side="left", fill="x", expand=True)
        # ttk.Button(
            # catalog_path_frame,
            # text="📁",
            # command=self.select_catalog_directory,
            # width=3
        # ).pack(side="left", padx=2)
        
        # ttk.Button(
            # astro_frame, 
            # text="⭐ Astrométrie Locale G (Catalogues binaires)", 
            # command=self.run_astrometry_local_g
        # ).pack(pady=2, padx=5, fill="x")

        # ============================================================
        # SECTION 3 : ALIGNEMENT & EMPILEMENT
        # ============================================================
        post_frame = ttk.LabelFrame(left_frame, text="4. Post-Traitement")
        post_frame.pack(fill="x", padx=5, pady=10)

        # Bouton Alignement
        ttk.Button(
            post_frame,
            text="📐 Aligner Images (WCS)",
            command=self.start_alignment_thread
        ).pack(pady=2, padx=5, fill="x")

        # Bouton Empilement
        ttk.Button(
            post_frame,
            text="📚 Empiler Images (Stack)",
            command=self.start_stacking_thread
        ).pack(pady=2, padx=5, fill="x")

        # ============================================================
        # SECTION 5 : VISUALISATION
        # ============================================================
        viz_frame = ttk.LabelFrame(left_frame, text="5. Visualisation")
        viz_frame.pack(fill="x", padx=5, pady=10)

        ttk.Button(
            viz_frame,
            text="📁 Définir répertoire images",
            command=self._set_viewer_directory,
            width=30
        ).pack(pady=2, padx=5, fill="x")

        ttk.Button(
            viz_frame,
            text="🖼️ Ouvrir visualisation",
            command=self._open_image_viewer,
            width=30
        ).pack(pady=2, padx=5, fill="x")

        # ============================================================
        # LOGS & PROGRESSION (Panneau Droit)
        # ============================================================
        log_label = ttk.Label(
            right_frame,
            text="📜 Journal des événements",
            font=("Arial", 10, "bold"),
        )
        log_label.pack(fill=tk.X)

        self.log_text = tk.Text(right_frame, height=15, width=80, state="disabled")
        self.log_text.pack(fill=tk.BOTH, expand=True)

        self.progress_label = ttk.Label(right_frame, text="📊 Progression : 0%")
        self.progress_label.pack(fill=tk.X, pady=5)

        self.progress_bar = ttk.Progressbar(
            right_frame,
            length=200,
            mode="determinate",
            maximum=100,
        )
        self.progress_bar.pack(fill=tk.X, padx=5)

    # ------------------------------------------------------------------
    # Exécution UI thread-safe
    # ------------------------------------------------------------------
    def _start_ui_queue_poller(self):
        def _poll():
            while True:
                try:
                    func, args, kwargs = self._ui_queue.get_nowait()
                except queue.Empty:
                    break
                try:
                    func(*args, **kwargs)
                except Exception as e:
                    logging.error(f"Erreur UI queue: {e}")
            self.frame.after(50, _poll)

        self.frame.after(50, _poll)

    def _call_on_ui_thread(self, func, *args, **kwargs):
        if threading.current_thread() is threading.main_thread():
            func(*args, **kwargs)
        else:
            self._ui_queue.put((func, args, kwargs))

    def _showerror(self, title, message):
        self._call_on_ui_thread(messagebox.showerror, title, message)

    def _showwarning(self, title, message):
        self._call_on_ui_thread(messagebox.showwarning, title, message)

    def _showinfo(self, title, message):
        self._call_on_ui_thread(messagebox.showinfo, title, message)

    # ------------------------------------------------------------------
    # Log
    # ------------------------------------------------------------------
    def log_message(self, message, level="info"):
        def _append_log():
            self.log_text.configure(state="normal")
            self.log_text.insert(tk.END, message + "\n")
            self.log_text.configure(state="disabled")
            self.log_text.yview(tk.END)

        self._call_on_ui_thread(_append_log)

        if level == "info":
            logging.info(message)
        elif level == "warning":
            logging.warning(message)
        elif level == "error":
            logging.error(message)

    # ------------------------------------------------------------------
    # Sélection des fichiers
    # ------------------------------------------------------------------
    def load_files(self):
        files = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits")])
        if files:
            self.files = list(files)
            self.log_message(f"📂 {len(self.files)} fichiers d'images chargés.")

    def load_bias(self):
        files = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits")])
        if files:
            self.bias_files = list(files)
            self.log_message(f"📂 {len(self.bias_files)} fichiers Bias chargés.")

    def load_darks(self):
        files = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits")])
        if files:
            self.dark_files = list(files)
            self.log_message(f"📂 {len(self.dark_files)} fichiers Dark chargés.")

    def load_flats(self):
        files = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits")])
        if files:
            self.flat_files = list(files)
            self.log_message(f"📂 {len(self.flat_files)} fichiers Flat chargés.")

    # ------------------------------------------------------------------
    # Répertoire de travail
    # ------------------------------------------------------------------
    def set_output_directory(self):
        directory = filedialog.askdirectory(title="Sélectionner le répertoire de travail")
        if directory:
            self.output_dir = Path(directory)
            self.processor = ImageProcessor(base_dir=self.output_dir)
            self.log_message(f"📁 Répertoire de travail défini : {self.output_dir}")
    
    def select_catalog_directory(self):
        """Sélectionne le répertoire contenant les catalogues binaires."""
        initial_dir = self.catalog_dir_var.get() if self.catalog_dir_var.get() else str(Path.home())
        directory = filedialog.askdirectory(
            title="Sélectionner le répertoire contenant les catalogues binaires (.1476)",
            initialdir=initial_dir
        )
        if directory:
            self.catalog_dir_var.set(directory)
            self.log_message(f"📁 Répertoire catalogues défini : {directory}")

            # Mise à jour de photometry_exoplanets_tab.base_dir
            try:
                notebook = self.frame.master       # ttk.Notebook
                main_window = notebook.master      # MainWindow (normalement)
                science_dir = self.output_dir / "science"
                
                # Vérifie que c'est bien une instance de MainWindow
                if hasattr(main_window, "photometry_exoplanets_tab"):
                    main_window.photometry_exoplanets_tab.base_dir = str(science_dir)
                    logging.info(f"📂 Onglet exoplanètes mis à jour : {science_dir}")
                else:
                    logging.warning("⚠️ MainWindow ne contient pas photometry_exoplanets_tab.")
            except Exception as e:
                logging.warning(f"Impossible de mettre à jour PhotometryExoplanetsTab.base_dir : {e}")


    # ------------------------------------------------------------------
    # Calibration
    # ------------------------------------------------------------------
    def run_calibration(self):
        if not self.files or not self.output_dir:
            self._showerror(
                "❌ Erreur",
                "Aucune image ou répertoire de travail non spécifié.",
            )
            return

        self.progress_prefix = "📊 Calibration :"
        self.update_progress(0)

        thread = threading.Thread(target=self.calibration_task, daemon=True)
        thread.start()

    def calibration_task(self):
        self.log_message("🚀 Début de la calibration...")
        try:
            self.processor.process_calibration(
                self.files,
                self.bias_files,
                self.dark_files,
                self.flat_files,
                progress_callback=self.update_progress,
                scale_darks=self.scale_darks_var.get(),
            )
            self.update_progress(100)
            self.log_message("✅ Calibration terminée. Lancez l'astrométrie maintenant.")
        except Exception as e:
            self.log_message(f"❌ Erreur durant la calibration : {e}", "error")
            self.update_progress(0)

    # ------------------------------------------------------------------
    # Astrométrie locale / NOVA
    # ------------------------------------------------------------------
    def run_astrometry_local(self):
        if self.processor is None:
            self._showerror(
                "❌ Erreur",
                "Aucun répertoire de travail défini.\n"
                "Veuillez d'abord choisir le répertoire et lancer la calibration.",
            )
            return

        calibrated_dir = self.processor.calibrated_dir
        if not list(calibrated_dir.glob("*.fits")):
            self._showerror(
                "❌ Erreur",
                f"Aucune image calibrée trouvée dans le dossier :\n{calibrated_dir}",
            )
            return

        self.progress_prefix = "📊 Astrométrie locale :"
        self.update_progress(0)

        thread = threading.Thread(
            target=self.astrometry_task,
            args=("LOCAL",),
            daemon=True,
        )
        thread.start()

    def run_astrometry_local_g(self):
        """Lance l'astrométrie avec le solveur local utilisant les catalogues binaires (magnitude G)."""
        if self.processor is None:
            self._showerror(
                "❌ Erreur",
                "Aucun répertoire de travail défini.\n"
                "Veuillez d'abord choisir le répertoire et lancer la calibration.",
            )
            return

        calibrated_dir = self.processor.calibrated_dir
        if not list(calibrated_dir.glob("*.fits")):
            self._showerror(
                "❌ Erreur",
                f"Aucune image calibrée trouvée dans le dossier :\n{calibrated_dir}",
            )
            return

        self.progress_prefix = "📊 Astrométrie locale G :"
        self.update_progress(0)

        thread = threading.Thread(
            target=self.astrometry_task,
            args=("LOCAL_G",),
            daemon=True,
        )
        thread.start()

    def run_astrometry_nova(self):
        if self.processor is None:
            self._showerror(
                "❌ Erreur",
                "Aucun répertoire de travail défini.\n"
                "Veuillez d'abord choisir le répertoire et lancer la calibration.",
            )
            return

        calibrated_dir = self.processor.calibrated_dir
        if not list(calibrated_dir.glob("*.fits")):
            self._showerror(
                "❌ Erreur",
                f"Aucune image calibrée trouvée dans le dossier :\n{calibrated_dir}",
            )
            return

        self.progress_prefix = "📊 Astrométrie NOVA :"
        self.update_progress(0)

        thread = threading.Thread(
            target=self.astrometry_task,
            args=("NOVA",),
            daemon=True,
        )
        thread.start()
    
    def astrometry_task(self, method="LOCAL"):
        self.log_message(f"🚀 Démarrage de l'astrométrie méthode {method}...")

        if self.processor is None:
            self.log_message("❌ ImageProcessor non initialisé.", "error")
            self.update_progress(0)
            return

        calibrated_dir = self.processor.calibrated_dir
        output_fits = sorted(calibrated_dir.glob("*.fits"))
        if not output_fits:
            self.log_message(
                f"⚠️ Aucun fichier calibré trouvé dans le dossier {calibrated_dir}.",
                "warning",
            )
            self.update_progress(0)
            return

        method = method.upper()

        if method == "LOCAL":
            self.processor.bash_path = self.bash_path
            try:
                self.processor.process_astrometry(
                    method="LOCAL",
                    progress_callback=self.update_progress,
                )
                self.log_message("✅ Astrométrie locale terminée.")
            except Exception as e:
                self.log_message(f"❌ Erreur durant l'astrométrie locale : {e}", "error")
                self.update_progress(0)

        elif method == "LOCAL_G":
            try:
                # Récupérer le répertoire de catalogues
                catalog_dir_str = self.catalog_dir_var.get().strip()
                if not catalog_dir_str:
                    self._showerror(
                        "Erreur",
                        "Veuillez sélectionner le répertoire contenant les catalogues binaires (.1476)"
                    )
                    self.update_progress(0)
                    return
                
                catalog_dir = Path(catalog_dir_str)
                if not catalog_dir.exists():
                    self._showerror(
                        "Erreur",
                        f"Répertoire introuvable : {catalog_dir}\nVeuillez vérifier le chemin."
                    )
                    self.update_progress(0)
                    return
                
                self.log_message("🔄 Démarrage astrométrie locale G (catalogues binaires)...")
                self.log_message(f"   Répertoire catalogues : {catalog_dir}")
                self.log_message("   Utilise directement la magnitude G de Gaia")
                self.processor.process_astrometry(
                    method="LOCAL_G",
                    progress_callback=self.update_progress,
                    catalog_dir=catalog_dir
                )
                self.log_message("✅ Astrométrie locale G terminée.")
            except Exception as e:
                self.log_message(f"❌ Erreur durant l'astrométrie locale G : {e}", "error")
                self._showerror("Erreur", f"Erreur durant l'astrométrie locale G : {e}")
                self.update_progress(0)
        
        elif method == "NOVA":
            try:
                solver = AstrometrySolverNova(
                    output_dir=self.processor.astrometry_dir,
                )
                solver.solve_directory(
                    calibrated_dir,
                    progress_callback=self.update_progress,
                )
                self.log_message("✅ Astrométrie NOVA terminée.")
            except Exception as e:
                self.log_message(f"❌ Erreur durant l'astrométrie NOVA : {e}", "error")
                self.update_progress(0)

        else:
            self.log_message(f"❌ Méthode inconnue : {method}", "error")
            self.update_progress(0)
            return

        self.log_message("🎯 Astrométrie terminée.")
            
    def start_alignment_thread(self):
        threading.Thread(target=self.run_alignment, daemon=True).start()

    def run_alignment(self):
        if not self.processor: return
            
        # Source : Science | Destination : Science/Aligned (défini dans processor.aligned_dir)
        input_dir = self.processor.science_dir
        output_dir = self.processor.aligned_dir 

        if not list(input_dir.glob("*.fits")):
            self._showwarning("Attention", f"Aucune image dans {input_dir.name}")
            return

        self.progress_prefix = "📐 Alignement :"
        try:
            self.processor.process_alignment_wcs(
                input_dir=input_dir,
                output_dir=output_dir, # Force l'écriture dans science/aligned
                progress_callback=self.update_progress
            )
            self.log_message(f"✅ Images alignées dans : {output_dir}")
            self._showinfo("Succès", "Alignement terminé.")
        except Exception as e:
            self.log_message(f"❌ Erreur alignement : {e}", "error")
        finally:
            self.update_progress(0)

    def start_stacking_thread(self):
        threading.Thread(target=self.run_stacking, daemon=True).start()

    def run_stacking(self):
        if not self.processor: return
            
        # On ouvre par défaut le dossier 'science/aligned'
        initial_dir = self.processor.aligned_dir
        if not initial_dir.exists():
            initial_dir = self.processor.science_dir

        files = filedialog.askopenfilenames(
            title="Sélectionnez les images à empiler",
            initialdir=initial_dir,
            filetypes=[("FITS", "*.fits")]
        )
        if not files: return

        save_path = filedialog.asksaveasfilename(
            title="Enregistrer le Master sous...",
            initialdir=self.processor.science_dir,
            initialfile="Master_Stacked.fits",
            defaultextension=".fits"
        )
        if not save_path: return

        self.progress_prefix = "📚 Empilement :"
        try:
            self.processor.process_stacking(
                input_files=[Path(f) for f in files],
                output_path=Path(save_path),
                progress_callback=self.update_progress
            )
            self.log_message(f"✅ Master créé : {Path(save_path).name}")
            self._showinfo("Succès", "Empilement terminé.")
        except Exception as e:
            self.log_message(f"❌ Erreur empilement : {e}", "error")
        finally:
            self.update_progress(0)

    # ------------------------------------------------------------------
    # Section 5 : Visualisation
    # ------------------------------------------------------------------
    def _set_viewer_directory(self):
        """Définit le répertoire contenant les images à visualiser."""
        directory = filedialog.askdirectory(title="Répertoire des images à visualiser")
        if directory:
            self.viewer_directory = Path(directory)
            self.log_message(f"📁 Répertoire visualisation : {self.viewer_directory}")

    def _open_image_viewer(self):
        """Ouvre la fenêtre de visualisation des images (navigation image par image)."""
        ImageViewerWindow(self)

    # ------------------------------------------------------------------
    # Progression (thread-safe)
    # ------------------------------------------------------------------
    def update_progress(self, percent):
        try:
            p = max(0, min(100, float(percent)))
        except Exception:
            p = 0.0

        def _do_update():
            self.progress_bar["value"] = p
            self.progress_label.config(
                text=f"{self.progress_prefix} {p:.0f}%"
            )
        self._call_on_ui_thread(_do_update)


# =============================================================================
# Fenêtre de visualisation d'images (section 5 Réduction) — inspirée de l'onglet Astéroïdes
# =============================================================================
class ImageViewerWindow(Toplevel):
    """Fenêtre de visualisation des images FITS d'un répertoire avec navigation."""

    def __init__(self, parent_tab):
        super().__init__()
        self.parent_tab = parent_tab
        self.title("Visualisation des images")
        # Taille augmentée de 30 % (900x700 → 1170x910)
        self.geometry("1170x910")

        # Répertoire : celui défini dans l'onglet ou demande à l'utilisateur
        directory = parent_tab.viewer_directory
        if directory is None or not Path(directory).exists():
            directory = filedialog.askdirectory(title="Choisir le répertoire des images")
            if not directory:
                self.destroy()
                return
            parent_tab.viewer_directory = Path(directory)

        self.directory = Path(directory)
        self.image_files = sorted([str(f) for f in self.directory.glob("*.fits")])
        if not self.image_files:
            messagebox.showwarning(
                "Aucune image",
                f"Aucun fichier FITS trouvé dans :\n{self.directory}",
                parent=self
            )
            self.destroy()
            return

        self.current_index = 0
        self.current_data = None
        self.auto_play_active = False
        self.auto_play_job = None
        self.auto_play_delay = 25  # ms

        # Contenu
        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        # Label image courante
        self.info_label = ttk.Label(main, text="", font=("", 10))
        self.info_label.pack(fill=tk.X, pady=(0, 5))

        # Figure matplotlib
        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=main)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, pady=5)

        # Boutons de navigation (barre centrée sous les images, boutons élargis de 50 %)
        nav_wrapper = ttk.Frame(main)
        nav_wrapper.pack(fill=tk.X, pady=5)
        nav = ttk.Frame(nav_wrapper)
        nav.pack(anchor=tk.CENTER)

        ttk.Button(nav, text="⏮ Première", command=self._first, width=18).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav, text="◀ Précédent", command=self._previous, width=18).pack(side=tk.LEFT, padx=2)
        self.play_btn = ttk.Button(nav, text="▶ Défilement auto", command=self._toggle_auto_play, width=21)
        self.play_btn.pack(side=tk.LEFT, padx=2)
        ttk.Button(nav, text="Suivant ▶", command=self._next, width=18).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav, text="Dernière ⏭", command=self._last, width=18).pack(side=tk.LEFT, padx=2)

        # Vitesse du défilement (ms) — centrée
        speed_wrapper = ttk.Frame(main)
        speed_wrapper.pack(fill=tk.X)
        speed_frame = ttk.Frame(speed_wrapper)
        speed_frame.pack(anchor=tk.CENTER)
        ttk.Label(speed_frame, text="Vitesse défilement (ms):").pack(side=tk.LEFT, padx=(0, 5))
        self.delay_var = tk.IntVar(value=self.auto_play_delay)
        ttk.Entry(speed_frame, textvariable=self.delay_var, width=6).pack(side=tk.LEFT, padx=2)

        # Charger et afficher la première image
        self._load_image_at(self.current_index)
        self._update_info()
        self._refresh_display()

    def _load_image_at(self, index):
        """Charge l'image à l'index donné."""
        if not 0 <= index < len(self.image_files):
            return
        self.current_index = index
        path = self.image_files[index]
        try:
            with fits.open(path) as hdul:
                self.current_data = hdul[0].data.astype(float)
        except Exception as e:
            logging.warning(f"Impossible de charger {path}: {e}")
            self.current_data = None
        self._update_info()
        self._refresh_display()

    def _update_info(self):
        """Met à jour le label d'information."""
        n = len(self.image_files)
        if n == 0:
            self.info_label.config(text="Aucune image")
            return
        name = Path(self.image_files[self.current_index]).name
        self.info_label.config(
            text=f"Image {self.current_index + 1} / {n} — {name}"
        )

    def _refresh_display(self):
        """Affiche l'image courante avec ZScale."""
        if self.current_data is None:
            self.ax.clear()
            self.canvas.draw()
            return
        self.ax.clear()
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(self.current_data)
        self.ax.imshow(self.current_data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        ny, nx = self.current_data.shape
        self.ax.set_xlim(0, nx)
        self.ax.set_ylim(0, ny)
        self.canvas.draw()

    def _first(self):
        """Aller à la première image."""
        self._stop_auto_play()
        self._load_image_at(0)

    def _previous(self):
        """Image précédente (ou dernière si on est à la première)."""
        self._stop_auto_play()
        if self.current_index > 0:
            self._load_image_at(self.current_index - 1)
        else:
            self._load_image_at(len(self.image_files) - 1)

    def _next(self):
        """Image suivante (ou première si on est à la dernière)."""
        self._stop_auto_play()
        if self.current_index < len(self.image_files) - 1:
            self._load_image_at(self.current_index + 1)
        else:
            self._load_image_at(0)

    def _last(self):
        """Aller à la dernière image."""
        self._stop_auto_play()
        self._load_image_at(len(self.image_files) - 1)

    def _toggle_auto_play(self):
        """Démarre ou arrête le défilement automatique."""
        if not self.image_files:
            return
        if self.auto_play_active:
            self._stop_auto_play()
        else:
            self._start_auto_play()

    def _start_auto_play(self):
        """Démarre le défilement automatique."""
        if not self.image_files:
            return
        self.auto_play_active = True
        self.play_btn.config(text="⏸ Arrêter défilement")
        self.auto_play_delay = max(25, self.delay_var.get())
        self.delay_var.set(self.auto_play_delay)
        self._auto_play_step()

    def _stop_auto_play(self):
        """Arrête le défilement automatique."""
        self.auto_play_active = False
        if hasattr(self, "play_btn"):
            self.play_btn.config(text="▶ Défilement auto")
        if self.auto_play_job:
            self.after_cancel(self.auto_play_job)
            self.auto_play_job = None

    def _auto_play_step(self):
        """Une étape du défilement automatique."""
        if not self.auto_play_active or not self.image_files:
            self._stop_auto_play()
            return
        next_index = (self.current_index + 1) % len(self.image_files)
        self._load_image_at(next_index)
        self.auto_play_job = self.after(self.auto_play_delay, self._auto_play_step)

    def destroy(self):
        """À la fermeture, arrêter le défilement auto."""
        self._stop_auto_play()
        super().destroy()