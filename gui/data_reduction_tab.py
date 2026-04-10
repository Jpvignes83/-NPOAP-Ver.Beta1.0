import logging
import threading
import queue
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, Toplevel
from pathlib import Path
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from scipy.ndimage import maximum_filter

from core.image_processor import ImageProcessor, PipelineControl
from gui.asteroid_photometry_tab import estimate_fwhm_marginal

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
        self.progress_prefix = "📊 Progression :"
        self.astrometry_thread = None
        self.astrometry_cancelled = False
        self.astrometry_process = None
        self.viewer_directory = None  # Répertoire pour la visualisation (section 4)
        self.resize_directory = None  # Répertoire source pour redimensionnement (section 5)
        self.quality_rows = []        # Données de qualité images (section Qualité)
        self._ui_queue = queue.Queue()
        self.pipeline_control = PipelineControl()
        self._progress_floor = 0.0
        self._progress_monotone = False
        self.resize_width_var = tk.StringVar(value="1024")
        self.resize_height_var = tk.StringVar(value="1024")
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

        astro_btns = ttk.Frame(astro_frame)
        astro_btns.pack(fill="x", padx=5, pady=2)
        self._pack_action_with_pipeline_row(
            astro_btns,
            text="🌐 Via Astrometry.net (NOVA)",
            command=self.run_astrometry_nova,
        )
        self._pack_action_with_pipeline_row(
            astro_btns,
            text="🖥️ Astrométrie locale (WSL)",
            command=self.run_astrometry_local,
        )

        # ============================================================
        # SECTION 3 : ALIGNEMENT, EMPILEMENT & QUALITÉ
        # ============================================================
        post_frame = ttk.LabelFrame(left_frame, text="3. Post-Traitement")
        post_frame.pack(fill="both", padx=5, pady=10, expand=False)

        post_btns = ttk.Frame(post_frame)
        post_btns.pack(fill="x", padx=5, pady=2)
        self._pack_action_with_pipeline_row(
            post_btns,
            text="📐 Aligner images (WCS)",
            command=self.start_alignment_thread,
        )
        self._pack_action_with_pipeline_row(
            post_btns,
            text="📚 Empiler images (Stack)",
            command=self.start_stacking_thread,
        )

        # Sous-section Qualité des images (liste triable)
        quality_frame = ttk.LabelFrame(post_frame, text="Qualité des images (science ou calibrées)")
        quality_frame.pack(fill="both", padx=5, pady=(8, 2), expand=True)

        ttk.Button(
            quality_frame,
            text="📊 Analyser la qualité des images",
            command=self.analyze_image_quality
        ).pack(pady=2, padx=5, fill="x")

        # Tableau de résultats
        cols = ("name", "fwhm", "ellipticity", "background", "snr", "stars")
        self.quality_tree = ttk.Treeview(quality_frame, columns=cols, show="headings", height=6)
        self.quality_tree.pack(fill="both", expand=True, pady=(4, 0))

        headers = {
            "name": "Nom",
            "fwhm": "FWHM [px]",
            "ellipticity": "Ellipticité",
            "background": "Fond médian",
            "snr": "SNR pic",
            "stars": "Stars level",
        }
        widths = {
            "name": 170,
            "fwhm": 80,
            "ellipticity": 90,
            "background": 90,
            "snr": 80,
            "stars": 80,
        }
        for cid in cols:
            self.quality_tree.heading(cid, text=headers[cid], command=lambda c=cid: self.sort_quality_by(c, False))
            self.quality_tree.column(cid, width=widths[cid], anchor="center")

        # ============================================================
        # SECTION 4 : VISUALISATION
        # ============================================================
        viz_frame = ttk.LabelFrame(left_frame, text="4. Visualisation")
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
        # SECTION 5 : ÉDITION DES IMAGES
        # ============================================================
        edit_frame = ttk.LabelFrame(left_frame, text="5. Édition des images")
        edit_frame.pack(fill="x", padx=5, pady=10)

        ttk.Button(
            edit_frame,
            text="📁 Définir dossier à redimensionner",
            command=self._set_resize_directory,
            width=30,
        ).pack(pady=2, padx=5, fill="x")

        size_frame = ttk.Frame(edit_frame)
        size_frame.pack(fill="x", padx=5, pady=4)
        ttk.Label(size_frame, text="Largeur (px) :").grid(row=0, column=0, sticky="w", padx=(0, 4))
        ttk.Entry(size_frame, textvariable=self.resize_width_var, width=8).grid(row=0, column=1, sticky="w", padx=(0, 10))
        ttk.Label(size_frame, text="Hauteur (px) :").grid(row=0, column=2, sticky="w", padx=(0, 4))
        ttk.Entry(size_frame, textvariable=self.resize_height_var, width=8).grid(row=0, column=3, sticky="w")

        ttk.Button(
            edit_frame,
            text="✂️ Redimensionner en batch (centré)",
            command=self._start_resize_batch_thread,
            width=30,
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

    def _pack_action_with_pipeline_row(
        self,
        parent: ttk.Frame,
        *,
        text: str,
        command,
    ) -> None:
        """Bouton d'action et Pause | Reprise | Stop sur la même ligne."""
        row = ttk.Frame(parent)
        row.pack(fill="x", pady=2)
        ctrl = ttk.Frame(row)
        ctrl.pack(side=tk.RIGHT, padx=(6, 0))
        ttk.Button(ctrl, text="Pause", command=self.pipeline_pause, width=8).pack(
            side=tk.LEFT, padx=1
        )
        ttk.Button(ctrl, text="Reprise", command=self.pipeline_resume, width=8).pack(
            side=tk.LEFT, padx=1
        )
        ttk.Button(ctrl, text="Stop", command=self.pipeline_stop, width=8).pack(
            side=tk.LEFT, padx=1
        )
        ttk.Button(row, text=text, command=command).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 2)
        )

    def _reset_pipeline_for_new_task(self) -> None:
        self.pipeline_control.reset()
        self._progress_floor = 0.0
        self._progress_monotone = False

    def pipeline_pause(self) -> None:
        self.pipeline_control.pause()
        self.log_message("⏸ Pipeline en pause (reprise possible).")

    def pipeline_resume(self) -> None:
        self.pipeline_control.resume()
        self.log_message("▶ Pipeline repris.")

    def pipeline_stop(self) -> None:
        self.pipeline_control.stop()
        self.log_message("⏹ Arrêt demandé : fin de l'étape en cours, puis arrêt entre fichiers.")

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
            filter_warn = self.processor.process_calibration(
                self.files,
                self.bias_files,
                self.dark_files,
                self.flat_files,
                progress_callback=self.update_progress,
                scale_darks=self.scale_darks_var.get(),
            )
            self.update_progress(100)
            if filter_warn:
                self.log_message(filter_warn, "warning")
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
        self._reset_pipeline_for_new_task()
        self._progress_monotone = True
        self.update_progress(0)

        thread = threading.Thread(
            target=self.astrometry_task,
            args=("LOCAL",),
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
        self._reset_pipeline_for_new_task()
        self._progress_monotone = True
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
            self._progress_monotone = False
            self.update_progress(0)
            return

        calibrated_dir = self.processor.calibrated_dir
        output_fits = sorted(calibrated_dir.glob("*.fits"))
        if not output_fits:
            self.log_message(
                f"⚠️ Aucun fichier calibré trouvé dans le dossier {calibrated_dir}.",
                "warning",
            )
            self._progress_monotone = False
            self.update_progress(0)
            return

        method = method.upper()
        success = False

        try:
            if method == "LOCAL":
                self.processor.bash_path = self.bash_path
                try:
                    self.processor.process_astrometry(
                        method="LOCAL",
                        progress_callback=self.update_progress,
                        pipeline_control=self.pipeline_control,
                    )
                    if self.pipeline_control.should_stop():
                        self.log_message("⚠️ Astrométrie locale interrompue (Stop).", "warning")
                    else:
                        self.log_message("✅ Astrométrie locale terminée.")
                        success = True
                except Exception as e:
                    self.log_message(f"❌ Erreur durant l'astrométrie locale : {e}", "error")
                    self.update_progress(0)

            elif method == "NOVA":
                try:
                    self.processor.process_astrometry(
                        method="NOVA",
                        progress_callback=self.update_progress,
                        pipeline_control=self.pipeline_control,
                    )
                    if self.pipeline_control.should_stop():
                        self.log_message("⚠️ Astrométrie NOVA interrompue (Stop).", "warning")
                    else:
                        self.log_message("✅ Astrométrie NOVA terminée.")
                        success = True
                except Exception as e:
                    self.log_message(f"❌ Erreur durant l'astrométrie NOVA : {e}", "error")
                    self.update_progress(0)

            else:
                self.log_message(f"❌ Méthode inconnue : {method}", "error")
                self.update_progress(0)
                return
        finally:
            self._progress_monotone = False
            self._progress_floor = 0.0

        if success:
            self.log_message("🎯 Astrométrie terminée.")
            self.update_progress(100)
            
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
        self._reset_pipeline_for_new_task()
        try:
            self.processor.process_alignment_wcs(
                input_dir=input_dir,
                output_dir=output_dir,
                progress_callback=self.update_progress,
                pipeline_control=self.pipeline_control,
            )
            if self.pipeline_control.should_stop():
                self.log_message("⚠️ Alignement interrompu (Stop).", "warning")
            else:
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
        self._reset_pipeline_for_new_task()
        try:
            self.processor.process_stacking(
                input_files=[Path(f) for f in files],
                output_path=Path(save_path),
                progress_callback=self.update_progress,
                pipeline_control=self.pipeline_control,
            )
            self.log_message(f"✅ Master créé : {Path(save_path).name}")
            self._showinfo("Succès", "Empilement terminé.")
        except RuntimeError as e:
            msg = str(e)
            if "interrompu" in msg.lower():
                self.log_message(f"⚠️ {msg}", "warning")
            else:
                self.log_message(f"❌ Erreur empilement : {e}", "error")
        except Exception as e:
            self.log_message(f"❌ Erreur empilement : {e}", "error")
        finally:
            self.update_progress(0)

    # ------------------------------------------------------------------
    # Section 4 : Visualisation
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
    # Section 5 : Édition des images (redimensionnement batch centré)
    # ------------------------------------------------------------------
    def _set_resize_directory(self):
        directory = filedialog.askdirectory(title="Répertoire des images à redimensionner")
        if directory:
            self.resize_directory = Path(directory)
            self.log_message(f"📁 Répertoire redimensionnement : {self.resize_directory}")

    def _start_resize_batch_thread(self):
        thread = threading.Thread(target=self._run_resize_batch, daemon=True)
        thread.start()

    def _run_resize_batch(self):
        if self.resize_directory is None or not self.resize_directory.exists():
            self._showwarning("Redimensionnement", "Veuillez d'abord définir un dossier source.")
            return

        try:
            target_width = int(self.resize_width_var.get())
            target_height = int(self.resize_height_var.get())
            if target_width <= 0 or target_height <= 0:
                raise ValueError
        except Exception:
            self._showerror("Redimensionnement", "Largeur et hauteur doivent être des entiers strictement positifs.")
            return

        fits_files = sorted(
            [
                *self.resize_directory.glob("*.fits"),
                *self.resize_directory.glob("*.fit"),
                *self.resize_directory.glob("*.fts"),
            ]
        )
        if not fits_files:
            self._showwarning("Redimensionnement", f"Aucun fichier FITS trouvé dans : {self.resize_directory}")
            return

        output_dir = self.resize_directory / "resized"
        output_dir.mkdir(parents=True, exist_ok=True)

        self.progress_prefix = "✂️ Redimensionnement :"
        self.update_progress(0)
        self.log_message(
            f"🚀 Redimensionnement batch de {len(fits_files)} images vers {target_width}x{target_height} px (sortie: {output_dir})"
        )

        done = 0
        skipped = 0
        total = len(fits_files)

        for idx, src_path in enumerate(fits_files, start=1):
            try:
                with fits.open(src_path) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header.copy()

                if data is None:
                    raise ValueError("données absentes")
                data = np.asarray(data)
                if data.ndim != 2:
                    raise ValueError(f"dimension non supportée ({data.ndim}D)")

                ny, nx = data.shape
                if target_width > nx or target_height > ny:
                    raise ValueError(
                        f"taille demandée {target_width}x{target_height} supérieure à {nx}x{ny}"
                    )

                # Centre défini par le pixel central (moitié des distances en pixels).
                cx = nx // 2
                cy = ny // 2
                x0 = cx - (target_width // 2)
                y0 = cy - (target_height // 2)
                x1 = x0 + target_width
                y1 = y0 + target_height

                cropped = data[y0:y1, x0:x1]
                if cropped.shape != (target_height, target_width):
                    raise ValueError("échec du découpage centré")

                header["NAXIS1"] = target_width
                header["NAXIS2"] = target_height

                dst_path = output_dir / src_path.name
                fits.PrimaryHDU(data=cropped, header=header).writeto(dst_path, overwrite=True)
                done += 1
            except Exception as e:
                skipped += 1
                self.log_message(f"⚠️ Ignorée {src_path.name} : {e}", "warning")

            self.update_progress(100.0 * idx / total)

        self.log_message(
            f"✅ Redimensionnement terminé : {done} image(s) écrite(s), {skipped} ignorée(s), dossier : {output_dir}"
        )
        self.update_progress(0)

    # ------------------------------------------------------------------
    # Qualité des images (sous-section du post-traitement, section 3)
    # ------------------------------------------------------------------
    def _image_quality_source_dir(self) -> Path | None:
        """
        Retourne le dossier à analyser pour la qualité :
        priorité au dossier science/aligned s'il existe, sinon science,
        sinon le répertoire de visualisation si défini.
        """
        if self.processor is None:
            return self.viewer_directory
        # priorité aux images alignées puis science
        aligned = getattr(self.processor, "aligned_dir", None)
        science = getattr(self.processor, "science_dir", None)
        if aligned is not None and aligned.exists() and list(aligned.glob("*.fits")):
            return aligned
        if science is not None and science.exists() and list(science.glob("*.fits")):
            return science
        return self.viewer_directory

    def analyze_image_quality(self):
        """Lance l'analyse de qualité des images et remplit le tableau Qualité (section 3)."""
        directory = self._image_quality_source_dir()
        if directory is None or not Path(directory).exists():
            self._showwarning("Qualité images", "Aucun répertoire d'images trouvé (science/aligned ou visualisation).")
            return

        files = sorted(Path(directory).glob("*.fits"))
        if not files:
            self._showwarning("Qualité images", f"Aucune image FITS trouvée dans : {directory}")
            return

        self.log_message(f"🔎 Analyse qualité sur {len(files)} images dans {directory}")
        self.quality_rows = []
        # Nettoyer le tableau
        for item in self.quality_tree.get_children():
            self.quality_tree.delete(item)

        def task():
            total = len(files)
            for idx, fpath in enumerate(files, start=1):
                row = self._compute_image_quality_metrics(fpath)
                if row is not None:
                    self.quality_rows.append(row)
                    self._call_on_ui_thread(self._insert_quality_row, row)
                # mise à jour progression globale
                percent = 100.0 * idx / total
                self.update_progress(percent)

            self.log_message("✅ Analyse de qualité terminée.")
            self.update_progress(0)

        thread = threading.Thread(target=task, daemon=True)
        thread.start()

    def _estimate_ellipticity(self, data: np.ndarray, x: int, y: int, box_size: int = 25) -> float:
        """
        Estime une ellipticité simple autour d'un maximum local.
        e = 1 - b/a avec a ≥ b les demi-axes principaux.
        """
        half = box_size // 2
        ny, nx = data.shape
        if y - half < 0 or y + half + 1 > ny or x - half < 0 or x + half + 1 > nx:
            return float("nan")

        sub = data[y - half:y + half + 1, x - half:x + half + 1].astype(float)
        if sub.size == 0 or not np.isfinite(sub).any():
            return float("nan")

        sub = sub - np.nanmin(sub)
        mask = np.isfinite(sub) & (sub > 0)
        if not mask.any():
            return float("nan")

        yy, xx = np.indices(sub.shape)
        w = sub[mask]
        xw = xx[mask]
        yw = yy[mask]

        w_sum = float(np.sum(w))
        if w_sum <= 0:
            return float("nan")

        x_mean = float(np.sum(xw * w) / w_sum)
        y_mean = float(np.sum(yw * w) / w_sum)

        dx = xw - x_mean
        dy = yw - y_mean
        cov_xx = float(np.sum(w * dx * dx) / w_sum)
        cov_yy = float(np.sum(w * dy * dy) / w_sum)
        cov_xy = float(np.sum(w * dx * dy) / w_sum)

        cov = np.array([[cov_xx, cov_xy], [cov_xy, cov_yy]])
        try:
            vals, _ = np.linalg.eigh(cov)
        except np.linalg.LinAlgError:
            return float("nan")

        vals = np.sort(vals)[::-1]
        if vals[0] <= 0 or vals[1] < 0:
            return float("nan")

        a = np.sqrt(vals[0])
        b = np.sqrt(max(vals[1], 0.0))
        if a <= 0:
            return float("nan")
        if b > a:
            a, b = b, a
        e = 1.0 - (b / a)
        return float(e)

    def _compute_image_quality_metrics(self, path: Path):
        """Calcule quelques métriques simples de qualité pour une image FITS."""
        try:
            with fits.open(path) as hdul:
                data = hdul[0].data
        except Exception as e:
            logging.warning(f"Lecture FITS impossible pour {path}: {e}")
            return None

        if data is None:
            return None

        data = np.asarray(data, dtype=float)
        if data.size == 0 or not np.isfinite(data).any():
            return None

        finite = np.isfinite(data)
        vals = data[finite]
        background = float(np.median(vals))
        noise = float(np.std(vals))
        peak = float(np.max(vals))
        snr_peak = float((peak - background) / noise) if noise > 0 else float("nan")

        # FWHM / ellipticité : candidat = pixel le plus brillant
        try:
            ny, nx = data.shape
            flat_index = int(np.nanargmax(data))
            y = flat_index // nx
            x = flat_index % nx
            fwhm_val, _, _ = estimate_fwhm_marginal(data, x, y, box_size=25)
            fwhm = float(fwhm_val) if fwhm_val is not None else float("nan")
            ellipticity = self._estimate_ellipticity(data, x, y, box_size=25)
        except Exception:
            fwhm = float("nan")
            ellipticity = float("nan")

        # Estimation très simple du "niveau d'étoiles" : nombre de maxima locaux significatifs
        try:
            thresh = background + 3.0 * noise
            sig = np.where(np.isfinite(data), data, background)
            data_max = maximum_filter(sig, size=5, mode="nearest")
            peaks = (sig == data_max) & (sig > thresh)
            stars = int(np.clip(peaks.sum(), 0, 9999))
        except Exception:
            stars = 0

        return {
            "name": path.name,
            "fwhm": fwhm,
            "ellipticity": ellipticity,
            "background": background,
            "snr": snr_peak,
            "stars": stars,
        }

    def _insert_quality_row(self, row: dict):
        """Insère une ligne dans le tableau de qualité."""
        def fmt(v):
            if isinstance(v, (int, np.integer)):
                return str(v)
            return "" if np.isnan(v) else f"{v:.2f}"

        self.quality_tree.insert(
            "",
            "end",
            values=(
                row["name"],
                fmt(row["fwhm"]),
                fmt(row["ellipticity"]),
                fmt(row["background"]),
                fmt(row["snr"]),
                fmt(row["stars"]),
            ),
        )

    def sort_quality_by(self, col: str, descending: bool):
        """
        Trie les lignes de la section qualité par colonne.
        Pour le nom, tri alphabétique ; pour les autres, tri numérique.
        """
        if not self.quality_rows:
            return

        if col == "name":
            self.quality_rows.sort(key=lambda r: r["name"], reverse=descending)
        else:
            self.quality_rows.sort(key=lambda r: (np.isnan(r[col]), r[col]), reverse=descending)

        # Re-remplir le tableau
        for item in self.quality_tree.get_children():
            self.quality_tree.delete(item)
        for row in self.quality_rows:
            self._insert_quality_row(row)

        # inverser le sens pour le prochain clic
        self.quality_tree.heading(col, command=lambda c=col: self.sort_quality_by(c, not descending))

    # ------------------------------------------------------------------
    # Progression (thread-safe)
    # ------------------------------------------------------------------
    def update_progress(self, percent):
        try:
            p = max(0, min(100, float(percent)))
        except Exception:
            p = 0.0
        if getattr(self, "_progress_monotone", False):
            p = max(self._progress_floor, p)
            self._progress_floor = p

        def _do_update():
            self.progress_bar["value"] = p
            self.progress_label.config(
                text=f"{self.progress_prefix} {p:.0f}%"
            )
        self._call_on_ui_thread(_do_update)


# =============================================================================
# Fenêtre de visualisation d'images (section 4 Réduction) — inspirée de l'onglet Astéroïdes
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
        # Contraste : 1.0 = plage ZScale par défaut ; >1 resserre la plage (plus de contraste)
        self.contrast_var = tk.DoubleVar(value=1.0)
        # Luminosité : 0 = centre ZScale ; décale la fenêtre d'affichage (fraction de la demi‑plage ZScale)
        self.brightness_var = tk.DoubleVar(value=0.0)

        # Contenu
        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        # Label image courante
        self.info_label = ttk.Label(main, text="", font=("", 10))
        self.info_label.pack(fill=tk.X, pady=(0, 5))

        # Réglage du contraste (relatif au stretch ZScale de l'image courante)
        contrast_row = ttk.Frame(main)
        contrast_row.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(contrast_row, text="Contraste").pack(side=tk.LEFT, padx=(0, 6))
        self.contrast_scale = ttk.Scale(
            contrast_row,
            from_=0.25,
            to=4.0,
            variable=self.contrast_var,
            command=lambda _s: self._on_stretch_changed(),
        )
        self.contrast_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 6))
        self.contrast_value_label = ttk.Label(contrast_row, text="1.00", width=5)
        self.contrast_value_label.pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(contrast_row, text="Défaut", command=self._reset_contrast, width=8).pack(
            side=tk.LEFT
        )

        brightness_row = ttk.Frame(main)
        brightness_row.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(brightness_row, text="Luminosité").pack(side=tk.LEFT, padx=(0, 6))
        self.brightness_scale = ttk.Scale(
            brightness_row,
            from_=-1.0,
            to=1.0,
            variable=self.brightness_var,
            command=lambda _s: self._on_stretch_changed(),
        )
        self.brightness_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 6))
        self.brightness_value_label = ttk.Label(brightness_row, text="+0.00", width=5)
        self.brightness_value_label.pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(brightness_row, text="Zéro", command=self._reset_brightness, width=8).pack(
            side=tk.LEFT
        )

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

    def _on_stretch_changed(self):
        self._update_contrast_label()
        self._update_brightness_label()
        self._refresh_display()

    def _update_contrast_label(self):
        try:
            c = float(self.contrast_var.get())
            self.contrast_value_label.config(text=f"{c:.2f}")
        except (tk.TclError, ValueError, TypeError):
            self.contrast_value_label.config(text="—")

    def _update_brightness_label(self):
        try:
            b = float(self.brightness_var.get())
            self.brightness_value_label.config(text=f"{b:+.2f}")
        except (tk.TclError, ValueError, TypeError):
            self.brightness_value_label.config(text="—")

    def _reset_contrast(self):
        self.contrast_var.set(1.0)
        self._on_stretch_changed()

    def _reset_brightness(self):
        self.brightness_var.set(0.0)
        self._on_stretch_changed()

    def _display_limits_from_zscale(self, data: np.ndarray) -> tuple[float, float]:
        """vmin, vmax : ZScale, puis luminosité (décalage du centre) et contraste (largeur de plage)."""
        interval = ZScaleInterval()
        vmin0, vmax0 = interval.get_limits(data)
        if not np.isfinite(vmin0) or not np.isfinite(vmax0) or vmin0 >= vmax0:
            return vmin0, vmax0
        try:
            c = float(self.contrast_var.get())
        except (tk.TclError, ValueError, TypeError):
            c = 1.0
        c = max(0.05, min(50.0, c))
        try:
            b = float(self.brightness_var.get())
        except (tk.TclError, ValueError, TypeError):
            b = 0.0
        b = max(-2.0, min(2.0, b))
        half = 0.5 * (float(vmax0) - float(vmin0))
        mid = 0.5 * (float(vmin0) + float(vmax0))
        mid_adj = mid + b * half
        half_adj = half / c
        return mid_adj - half_adj, mid_adj + half_adj

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
        """Affiche l'image courante avec ZScale et facteur de contraste utilisateur."""
        if self.current_data is None:
            self.ax.clear()
            self.canvas.draw()
            return
        self.ax.clear()
        vmin, vmax = self._display_limits_from_zscale(self.current_data)
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