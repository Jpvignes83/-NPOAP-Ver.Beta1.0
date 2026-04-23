# gui/transient_photometry_tab.py
"""
Onglet de photométrie des transitoires utilisant STDPipe.
Fonctionnalités:
- Astrométrie automatique
- Téléchargement d'images de référence (Pan-STARRS, SDSS, etc.)
- Soustraction d'images
- Détection de transitoires
- Photométrie calibrée
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from threading import Thread
import json
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table
import astropy.units as u
from astropy.visualization import ZScaleInterval
import numpy as np
import os
import logging
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import Circle

# Import du wrapper STDPipe
try:
    from core.stdpipe_wrapper import STDPipeWrapper, STDPIPE_AVAILABLE, STDPIPE_ERROR
except ImportError:
    STDPIPE_AVAILABLE = False
    STDPIPE_ERROR = "Impossible d'importer le module stdpipe_wrapper"

# Import du client Astro-COLIBRI (remplace TNS pour la recherche de transitoires)
try:
    from core.astro_colibri_client import AstroColibriClient
    ASTRO_COLIBRI_AVAILABLE = True
except ImportError:
    ASTRO_COLIBRI_AVAILABLE = False
    AstroColibriClient = None

logger = logging.getLogger(__name__)


class TransientPhotometryTab(ttk.Frame):
    """
    Onglet pour l'analyse de transitoires utilisant STDPipe.
    """
    
    def __init__(self, parent, base_dir=None):
        super().__init__(parent, padding=10)

        if base_dir is None:
            self.base_dir = Path.home()
        else:
            self.base_dir = Path(base_dir)

        # Variables
        self.target_coord = None
        self.science_image_path = None
        self.reference_image_path = None
        self.subtracted_image_path = None
        self.transients_table = None
        self.stdpipe_wrapper = None
        
        # Recherche transitoires (Astro-COLIBRI)
        self.transient_search_results = []
        
        # Initialiser STDPipe si disponible
        if STDPIPE_AVAILABLE:
            try:
                self.stdpipe_wrapper = STDPipeWrapper()
                logger.info("STDPipe wrapper initialisé")
            except Exception as e:
                logger.error(f"Erreur initialisation STDPipe: {e}")
                messagebox.showwarning(
                    "STDPipe",
                    f"STDPipe n'est pas correctement installé:\n{e}\n\n"
                    "Certaines fonctionnalités seront limitées."
                )
        
        self.create_widgets()

    def create_widgets(self):
        # Layout principal (uniquement la partie STDPipe à gauche)
        container = ttk.Frame(self, padding=5)
        container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        sne_fr = ttk.LabelFrame(container, text="Courbes SN Ia (optionnel, sncosmo)", padding=8)
        sne_fr.pack(fill=tk.X, pady=(0, 6))
        ttk.Label(
            sne_fr,
            text="Ajustement SALT2 / SALT3 sur CSV (MJD, filtre Gaia G / G_Bp / G_Rp, mag, erreur). "
            "Profil d’installation : requirements-cosmology-sne.txt",
            wraplength=760,
            justify=tk.LEFT,
        ).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8))
        ttk.Button(sne_fr, text="Ouvrir la boîte de dialogue…", command=self._open_sncosmo_dialog).pack(
            side=tk.RIGHT
        )

        # Frame principal (workflow STDPipe existant)
        left_frame = ttk.Frame(container, padding=5)
        left_frame.pack(fill=tk.BOTH, expand=True)
        
        # ========== PARTIE GAUCHE ==========
        # Titre
        title = ttk.Label(left_frame, text="Photométrie Transitoires (STDPipe)", 
                         font=("Helvetica", 12, "bold"))
        title.grid(row=0, column=0, columnspan=3, pady=10)
        
        # Vérification disponibilité STDPipe
        if not STDPIPE_AVAILABLE:
            error_msg = "⚠️ STDPipe n'est pas disponible.\n\n"
            if 'STDPIPE_ERROR' in globals() and STDPIPE_ERROR:
                error_msg += f"Erreur détectée: {STDPIPE_ERROR}\n\n"
            error_msg += "Pour installer:\n"
            error_msg += "  conda activate astroenv\n"
            error_msg += "  pip install stdpipe\n\n"
            error_msg += "Assurez-vous d'utiliser le même environnement Python que NPOAP."
            
            warning = ttk.Label(left_frame, 
                               text=error_msg,
                               foreground="red", font=("Helvetica", 9),
                               justify=tk.LEFT)
            warning.grid(row=1, column=0, columnspan=3, pady=5, padx=10, sticky="w")
            
            # Ajouter un bouton pour tester l'import
            def test_stdpipe():
                try:
                    import stdpipe
                    version = getattr(stdpipe, '__version__', 'version inconnue')
                    messagebox.showinfo("Test STDPipe", 
                                       f"STDPipe est importable !\nVersion: {version}\n\n"
                                       f"Mais l'application ne le détecte pas.\n"
                                       f"Vérifiez que vous utilisez le même environnement Python.\n\n"
                                       f"Erreur détectée: {STDPIPE_ERROR if 'STDPIPE_ERROR' in globals() else 'Aucune'}")
                except ImportError as e:
                    messagebox.showerror("Test STDPipe", 
                                        f"STDPipe n'est pas importable:\n{e}\n\n"
                                        f"Installez-le avec:\n"
                                        f"  conda activate astroenv\n"
                                        f"  pip install stdpipe")
                except Exception as e:
                    messagebox.showerror("Test STDPipe", f"Erreur: {e}")
            
            test_btn = ttk.Button(left_frame, text="🔍 Tester l'import STDPipe", command=test_stdpipe)
            test_btn.grid(row=2, column=0, columnspan=3, pady=5)
            return

        # Section 1: Chargement image science
        section1 = ttk.LabelFrame(left_frame, text="1. Image Science", padding=10)
        section1.grid(row=2, column=0, columnspan=3, sticky="ew", pady=5, padx=5)
        
        ttk.Label(section1, text="Image FITS calibrée:").grid(row=0, column=0, sticky="e", padx=5)
        self.science_fits_var = tk.StringVar()
        ttk.Entry(section1, textvariable=self.science_fits_var, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(section1, text="Parcourir", command=self.browse_science_fits).grid(row=0, column=2, padx=5)
        
        # Statut astrométrie
        self.astrometry_status_var = tk.StringVar(value="Non vérifié")
        self.astrometry_status_label = ttk.Label(section1, textvariable=self.astrometry_status_var, 
                                                 foreground="gray")
        self.astrometry_status_label.grid(row=1, column=0, columnspan=2, sticky="w", padx=5)
        
        ttk.Button(section1, text="🔭 Résoudre/Ré-astrométrie", 
                  command=self.solve_astrometry).grid(row=1, column=2, pady=5)

        # Section 2: Image de référence
        section2 = ttk.LabelFrame(left_frame, text="2. Image de Référence", padding=10)
        section2.grid(row=3, column=0, columnspan=3, sticky="ew", pady=5, padx=5)
        
        ttk.Label(section2, text="Catalogue:").grid(row=0, column=0, sticky="e", padx=5)
        self.catalog_var = tk.StringVar(value="panstarrs")
        catalog_combo = ttk.Combobox(section2, textvariable=self.catalog_var,
                                     values=["panstarrs", "sdss", "des"], state="readonly", width=15)
        catalog_combo.grid(row=0, column=1, padx=5)
        
        
        ttk.Label(section2, text="Filtre:").grid(row=0, column=2, sticky="e", padx=5)
        self.filter_var = tk.StringVar(value="r")
        filter_combo = ttk.Combobox(section2, textvariable=self.filter_var,
                                   values=["g", "r", "i", "z"], state="readonly", width=10)
        filter_combo.grid(row=0, column=3, padx=5)
        
        ttk.Button(section2, text="⬇️ Télécharger Image Référence", 
                  command=self.download_reference).grid(row=1, column=0, columnspan=4, pady=5)
        
        ttk.Label(section2, text="Ou charger image locale:").grid(row=2, column=0, sticky="e", padx=5)
        self.reference_fits_var = tk.StringVar()
        ttk.Entry(section2, textvariable=self.reference_fits_var, width=40).grid(row=2, column=1, columnspan=2, padx=5)
        ttk.Button(section2, text="Parcourir", command=self.browse_reference_fits).grid(row=2, column=3, padx=5)

        # Section 3: Soustraction
        section3 = ttk.LabelFrame(left_frame, text="3. Soustraction d'Images", padding=10)
        section3.grid(row=4, column=0, columnspan=3, sticky="ew", pady=5, padx=5)
        
        ttk.Label(section3, text="Méthode:").grid(row=0, column=0, sticky="e", padx=5)
        self.subtraction_method_var = tk.StringVar(value="simple")
        method_combo = ttk.Combobox(section3, textvariable=self.subtraction_method_var,
                                   values=["simple", "hotpants", "zogy", "alardlupton"], state="readonly", width=15)
        method_combo.grid(row=0, column=1, padx=5)

        # Options Alard-Lupton (utilisées uniquement si méthode = alardlupton)
        self.alard_poisson_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(section3, text="Pondération Poisson (Alard-Lupton)", variable=self.alard_poisson_var
                       ).grid(row=1, column=0, columnspan=2, sticky="w", padx=5)
        ttk.Label(section3, text="Gain (ADU/photon):").grid(row=2, column=0, sticky="e", padx=5)
        self.alard_gain_var = tk.DoubleVar(value=1.0)
        ttk.Entry(section3, textvariable=self.alard_gain_var, width=8).grid(row=2, column=1, sticky="w", padx=5)

        ttk.Button(section3, text="➖ Lancer Soustraction",
                  command=self.run_subtraction).grid(row=3, column=0, columnspan=2, pady=5)

        ttk.Button(section3, text="👁️ Visualiser Images",
                  command=self.show_subtraction_images).grid(row=4, column=0, columnspan=2, pady=5)

        # Section 4: Détection de transitoires
        section4 = ttk.LabelFrame(left_frame, text="4. Détection de Transitoires", padding=10)
        section4.grid(row=5, column=0, columnspan=3, sticky="ew", pady=5, padx=5)
        
        # Méthode de détection
        ttk.Label(section4, text="Méthode:").grid(row=0, column=0, sticky="e", padx=5)
        self.detection_method_var = tk.StringVar(value="photutils_segmentation")
        detection_method_combo = ttk.Combobox(section4, textvariable=self.detection_method_var,
                                             values=["photutils_segmentation", "photutils_dao", "photutils_iraf"],
                                             state="readonly", width=25)
        detection_method_combo.grid(row=0, column=1, padx=5)
        
        # FWHM (pour méthodes dao/iraf)
        ttk.Label(section4, text="FWHM (pixels):").grid(row=0, column=2, sticky="e", padx=5)
        self.fwhm_var = tk.DoubleVar(value=3.0)
        ttk.Entry(section4, textvariable=self.fwhm_var, width=10).grid(row=0, column=3, padx=5)
        
        # Seuil
        ttk.Label(section4, text="Seuil (sigma):").grid(row=1, column=0, sticky="e", padx=5)
        self.threshold_var = tk.DoubleVar(value=5.0)
        ttk.Entry(section4, textvariable=self.threshold_var, width=10).grid(row=1, column=1, padx=5)
        
        # Options avancées (déblending pour segmentation)
        self.deblend_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(section4, text="Déblending (segmentation)", 
                       variable=self.deblend_var).grid(row=1, column=2, columnspan=2, sticky="w", padx=5)
        
        ttk.Button(section4, text="🔍 Détecter Transitoires", 
                  command=self.detect_transients).grid(row=2, column=0, columnspan=4, pady=5)

        # Section 5: Photométrie
        section5 = ttk.LabelFrame(left_frame, text="5. Photométrie", padding=10)
        section5.grid(row=6, column=0, columnspan=3, sticky="ew", pady=5, padx=5)
        
        ttk.Label(section5, text="Catalogue calib.:").grid(row=0, column=0, sticky="e", padx=5)
        self.phot_catalog_var = tk.StringVar(value="gaia")
        phot_catalog_combo = ttk.Combobox(section5, textvariable=self.phot_catalog_var,
                                         values=["gaia", "panstarrs"], state="readonly", width=15)
        phot_catalog_combo.grid(row=0, column=1, padx=5)
        
        ttk.Label(section5, text="Filtre:").grid(row=0, column=2, sticky="e", padx=5)
        self.phot_filter_var = tk.StringVar(value="G")
        phot_filter_combo = ttk.Combobox(section5, textvariable=self.phot_filter_var,
                                        values=["G", "G_Bp", "G_Rp", "U", "B", "V", "R", "I"],
                                        state="readonly", width=12)
        phot_filter_combo.grid(row=0, column=3, padx=5)
        
        ttk.Label(section5, text="Méthode:").grid(row=1, column=0, sticky="e", padx=5)
        self.phot_method_var = tk.StringVar(value="aperture")
        phot_method_combo = ttk.Combobox(section5, textvariable=self.phot_method_var,
                                        values=["aperture", "psf"], state="readonly", width=15)
        phot_method_combo.grid(row=1, column=1, padx=5)
        
        ttk.Label(section5, text="FWHM (PSF, px):").grid(row=1, column=2, sticky="e", padx=5)
        self.phot_fwhm_var = tk.DoubleVar(value=3.0)
        ttk.Entry(section5, textvariable=self.phot_fwhm_var, width=10).grid(row=1, column=3, padx=5)

        self.use_npoap_trans_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            section5,
            text="Utiliser le coefficient de transformation (priorité : « Coefficient médian » du CSV NPOAP)",
            variable=self.use_npoap_trans_var,
        ).grid(row=2, column=0, columnspan=4, sticky="w", pady=2, padx=2)
        ttk.Label(section5, text="Filtre observateur (si vide : clé FILTER du FITS science) :").grid(
            row=3, column=0, columnspan=2, sticky="e", padx=5, pady=2
        )
        self.npoap_filter_override_var = tk.StringVar(value="")
        ttk.Entry(section5, textvariable=self.npoap_filter_override_var, width=20).grid(
            row=3, column=2, columnspan=2, sticky="w", padx=5, pady=2
        )

        ttk.Button(
            section5, text="📊 Effectuer Photométrie", command=self.perform_photometry
        ).grid(row=4, column=0, columnspan=4, pady=5)

        ttk.Button(
            section5, text="💾 Exporter Résultats", command=self.export_results
        ).grid(row=5, column=0, columnspan=4, pady=5)

        # Barre de progression
        self.progress = ttk.Progressbar(left_frame, length=400, mode="determinate")
        self.progress.grid(row=7, column=0, columnspan=3, pady=10, sticky="ew")
        
        # (Ancienne partie droite Astro-COLIBRI supprimée de l'interface)

    def _open_sncosmo_dialog(self):
        """Ajustement SN Ia (SALT2/SALT3) via sncosmo — dépendance optionnelle."""
        try:
            from gui.transient_sncosmo_dialog import open_transient_sncosmo_dialog

            top = self.winfo_toplevel()
            open_transient_sncosmo_dialog(top)
        except Exception as e:
            logger.exception("Ouverture dialogue sncosmo")
            messagebox.showerror("sncosmo", f"Impossible d'ouvrir l'outil SN Ia :\n{e}", parent=self.winfo_toplevel())

    def browse_science_fits(self):
        initial_dir = self.base_dir
        current_value = self.science_fits_var.get()
        if current_value:
            try:
                p = Path(current_value).parent
                if p.exists():
                    initial_dir = p
            except Exception:
                pass

        path = filedialog.askopenfilename(
            initialdir=initial_dir,
            filetypes=[("FITS", "*.fits"), ("FITS", "*.fit")]
        )
        if path:
            self.science_fits_var.set(path)
            self.science_image_path = path
            self.load_header(path)

    def browse_reference_fits(self):
        initial_dir = self.base_dir
        path = filedialog.askopenfilename(
            initialdir=initial_dir,
            filetypes=[("FITS", "*.fits"), ("FITS", "*.fit")]
        )
        if path:
            self.reference_fits_var.set(path)
            self.reference_image_path = path

    def load_header(self, path: str):
        try:
            with fits.open(path) as hdul:
                header = hdul[0].header
                wcs = WCS(header)
                if wcs.is_celestial:
                    ny, nx = hdul[0].data.shape
                    ra_dec = wcs.pixel_to_world(nx // 2, ny // 2)
                    self.target_coord = SkyCoord(ra=ra_dec.ra, dec=ra_dec.dec, frame="icrs")
                    
                    # Vérifier la qualité du WCS (coordonnées finies et déclinaison dans [-90, 90])
                    try:
                        test_coords = wcs.pixel_to_world([0, nx//2, nx], [0, ny//2, ny])
                        ok = all(
                            np.isfinite(c.ra.deg) and np.isfinite(c.dec.deg) and -90 <= c.dec.deg <= 90
                            for c in test_coords
                        )
                        if ok:
                            self.astrometry_status_var.set("✅ Astrométrie valide")
                            self.astrometry_status_label.config(foreground="green")
                        else:
                            self.astrometry_status_var.set("⚠️ WCS présent mais douteux")
                            self.astrometry_status_label.config(foreground="orange")
                    except Exception:
                        self.astrometry_status_var.set("⚠️ WCS présent mais invalide")
                        self.astrometry_status_label.config(foreground="orange")
                    
                    logger.info(f"Coordonnées du champ: RA={self.target_coord.ra.deg:.6f}°, "
                               f"Dec={self.target_coord.dec.deg:.6f}°")
                else:
                    self.astrometry_status_var.set("❌ Pas de WCS valide")
                    self.astrometry_status_label.config(foreground="red")
                    self.target_coord = None
        except Exception as e:
            logger.error(f"Erreur lecture header: {e}")
            self.astrometry_status_var.set("❌ Erreur lecture")
            self.astrometry_status_label.config(foreground="red")
            messagebox.showerror("Erreur", f"Erreur lors de la lecture de l'image:\n{e}")

    def solve_astrometry(self):
        if not self.stdpipe_wrapper:
            messagebox.showerror("Erreur", "STDPipe n'est pas disponible")
            return
        
        science_path = self.science_fits_var.get()
        if not science_path or not Path(science_path).exists():
            messagebox.showerror("Erreur", "Sélectionnez d'abord une image science")
            return
        
        # Vérifier si l'image a déjà un WCS valide
        try:
            with fits.open(science_path) as hdul:
                header = hdul[0].header
                wcs = WCS(header)
                if wcs.is_celestial:
                    # Vérifier que le WCS donne des coordonnées valides (finies, Dec dans [-90,90])
                    try:
                        ny, nx = hdul[0].data.shape
                        test_coords = wcs.pixel_to_world([nx//2], [ny//2])
                        c = test_coords[0]
                        if np.isfinite(c.ra.deg) and np.isfinite(c.dec.deg) and -90 <= c.dec.deg <= 90:
                            # Demander confirmation pour ré-astrométrie
                            result = messagebox.askyesno(
                                "Astrométrie existante",
                                "Cette image a déjà un WCS astrométrique valide.\n\n"
                                "Voulez-vous quand même ré-astrométriser l'image ?\n\n"
                                "Cliquez sur 'Non' pour utiliser le WCS existant."
                            )
                            if not result:
                                # Mettre à jour le statut et les coordonnées
                                self.load_header(science_path)
                                messagebox.showinfo("Info", "Utilisation du WCS existant.")
                                return
                    except Exception:
                        pass  # WCS invalide, continuer avec l'astrométrie
        except Exception as e:
            logger.debug(f"Erreur vérification WCS: {e}")
            # Continuer avec l'astrométrie
        
        self.progress["value"] = 0
        Thread(target=self._solve_astrometry_task, args=(science_path,), daemon=True).start()

    def _solve_astrometry_task(self, image_path: str):
        try:
            self.progress["value"] = 10
            logger.info(f"Résolution astrométrique de {Path(image_path).name}...")
            
            wcs = self.stdpipe_wrapper.solve_astrometry(image_path)
            
            if wcs:
                self.progress["value"] = 100
                messagebox.showinfo("Succès", "Astrométrie résolue avec succès!")
                # Mettre à jour les coordonnées et le statut
                self.load_header(image_path)
            else:
                self.progress["value"] = 0
                self.astrometry_status_var.set("❌ Échec astrométrie")
                self.astrometry_status_label.config(foreground="red")
                messagebox.showerror("Erreur", "Échec de la résolution astrométrique")
                
        except Exception as e:
            self.progress["value"] = 0
            logger.error(f"Erreur astrométrie: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de l'astrométrie:\n{e}")

    def download_reference(self):
        if not self.stdpipe_wrapper:
            messagebox.showerror("Erreur", "STDPipe n'est pas disponible")
            return
        
        # Vérifier si on a des coordonnées (depuis WCS existant ou astrométrie)
        if not self.target_coord:
            # Essayer de relire le header pour vérifier
            science_path = self.science_fits_var.get()
            if science_path and Path(science_path).exists():
                self.load_header(science_path)
            
            if not self.target_coord:
                messagebox.showwarning(
                    "Erreur", 
                    "Aucune coordonnée disponible.\n\n"
                    "Chargez une image avec WCS valide ou résolvez l'astrométrie d'abord."
                )
                return

        catalog = self.catalog_var.get()
        filter_name = self.filter_var.get()
        
        self.progress["value"] = 0
        Thread(target=self._download_task, args=(self.target_coord, catalog, filter_name), daemon=True).start()

    def _download_task(self, coord: SkyCoord, catalog: str, filter_name: str):
        try:
            self.progress["value"] = 20
            logger.info(f"Téléchargement image référence {catalog} filtre {filter_name}...")
            
            ref_image = self.stdpipe_wrapper.download_reference_image(
                coord, catalog=catalog, filter_name=filter_name
            )
            
            if ref_image:
                self.reference_image_path = ref_image
                self.reference_fits_var.set(ref_image)
                self.progress["value"] = 100
                messagebox.showinfo("Succès", f"Image de référence téléchargée:\n{ref_image}")
            else:
                self.progress["value"] = 0
                messagebox.showerror("Erreur", "Échec du téléchargement de l'image de référence")
                
        except Exception as e:
            self.progress["value"] = 0
            logger.error(f"Erreur téléchargement: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors du téléchargement:\n{e}")

    def run_subtraction(self):
        if not self.stdpipe_wrapper:
            messagebox.showerror("Erreur", "STDPipe n'est pas disponible")
            return

        science_path = self.science_fits_var.get()
        ref_path = self.reference_fits_var.get()
        
        if not science_path or not Path(science_path).exists():
            messagebox.showerror("Erreur", "Sélectionnez une image science valide")
            return

        if not ref_path or not Path(ref_path).exists():
            messagebox.showerror("Erreur", "Sélectionnez une image de référence valide")
            return

        # Chemin de sortie pour l'image soustraite
        science_dir = Path(science_path).parent
        output_path = science_dir / f"subtracted_{Path(science_path).name}"
        
        self.progress["value"] = 0
        Thread(target=self._subtraction_task, 
               args=(science_path, ref_path, str(output_path)), daemon=True).start()

    def _subtraction_task(self, science_path: str, ref_path: str, output_path: str):
        try:
            self.progress["value"] = 20
            logger.info(f"Soustraction: {Path(science_path).name} - {Path(ref_path).name}")
            
            method = self.subtraction_method_var.get()
            
            # Vérifier la compatibilité avant soustraction
            self.progress["value"] = 10
            logger.info("Vérification de la compatibilité des images...")
            compatibility_info = self.stdpipe_wrapper.check_image_compatibility(science_path, ref_path)
            
            if not compatibility_info['compatible']:
                error_msg = "Images incompatibles:\n" + "\n".join(compatibility_info['errors'])
                self.progress["value"] = 0
                messagebox.showerror("Erreur", error_msg)
                return
            
            # Afficher les avertissements si présents
            if compatibility_info['warnings']:
                warnings_msg = "Avertissements de compatibilité:\n\n" + "\n".join(compatibility_info['warnings'])
                warnings_msg += "\n\nSTDPipe effectuera automatiquement la mise à l'échelle et l'alignement.\n"
                warnings_msg += "Voulez-vous continuer ?"
                if not messagebox.askyesno("Avertissements", warnings_msg):
                    self.progress["value"] = 0
                    return
            
            # Informations détaillées sur les échelles si disponibles
            if compatibility_info.get('pixel_scale1') and compatibility_info.get('pixel_scale2'):
                logger.info(f"Échelle pixel science: {compatibility_info['pixel_scale1']:.3f}\"/pixel")
                logger.info(f"Échelle pixel référence: {compatibility_info['pixel_scale2']:.3f}\"/pixel")
                if compatibility_info.get('scale_ratio'):
                    logger.info(f"Ratio d'échelle: {compatibility_info['scale_ratio']:.3f}")
            
            self.progress["value"] = 20
            kwargs = {"check_compatibility": False}
            if method == "alardlupton":
                try:
                    kwargs["alard_lupton_use_poisson_weights"] = self.alard_poisson_var.get()
                    kwargs["alard_lupton_gain"] = float(self.alard_gain_var.get())
                except (tk.TclError, ValueError):
                    kwargs["alard_lupton_use_poisson_weights"] = True
                    kwargs["alard_lupton_gain"] = 1.0
            success, _ = self.stdpipe_wrapper.subtract_images(
                science_path, ref_path, output_path, method=method, **kwargs
            )
            
            if success:
                self.subtracted_image_path = output_path
                self.progress["value"] = 100
                # Ouvrir automatiquement la fenêtre de visualisation
                self.after(0, self.show_subtraction_images)
                self.after(0, lambda: messagebox.showinfo("Succès", f"Image soustraite sauvegardée:\n{output_path}"))
            else:
                self.progress["value"] = 0
                messagebox.showerror("Erreur", "Échec de la soustraction d'images")
                
        except Exception as e:
            self.progress["value"] = 0
            logger.error(f"Erreur soustraction: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de la soustraction:\n{e}")

    def detect_transients(self):
        if not self.stdpipe_wrapper:
            messagebox.showerror("Erreur", "STDPipe n'est pas disponible")
            return
        
        if not self.subtracted_image_path or not Path(self.subtracted_image_path).exists():
            messagebox.showerror("Erreur", "Lancez d'abord la soustraction d'images")
            return
        
        threshold = self.threshold_var.get()
        method = self.detection_method_var.get()
        fwhm = self.fwhm_var.get()
        deblend = self.deblend_var.get()
        
        self.progress["value"] = 0
        Thread(target=self._detect_transients_task, 
               args=(self.subtracted_image_path, threshold, method, fwhm, deblend), daemon=True).start()

    def _detect_transients_task(self, subtracted_image: str, threshold: float, 
                                method: str, fwhm: float, deblend: bool):
        try:
            self.progress["value"] = 20
            logger.info(f"Détection de transitoires: méthode={method}, seuil={threshold} sigma, FWHM={fwhm}")
            
            # Charger l'image
            with fits.open(subtracted_image) as hdul:
                data = hdul[0].data
                header = hdul[0].header
            
            self.progress["value"] = 40
            
            # Utiliser la méthode de détection choisie
            transients = self.stdpipe_wrapper.detect_sources(
                data,
                header=header,
                threshold_sigma=threshold,
                fwhm=fwhm if method in ['photutils_dao', 'photutils_iraf'] else None,
                method=method,
                deblend=deblend if method == 'photutils_segmentation' else False,
                aper=3.0,
                edge=10,
                minarea=5
            )
            
            self.progress["value"] = 80
            
            if transients and len(transients) > 0:
                self.transients_table = transients
                self.progress["value"] = 100
                messagebox.showinfo("Succès", f"{len(transients)} transitoire(s) détecté(s)!")
            else:
                self.progress["value"] = 0
                messagebox.showinfo("Info", "Aucun transitoire détecté avec ces paramètres")
                
        except Exception as e:
            self.progress["value"] = 0
            logger.error(f"Erreur détection: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de la détection:\n{e}")

    def perform_photometry(self):
        if not self.stdpipe_wrapper:
            messagebox.showerror("Erreur", "STDPipe n'est pas disponible")
            return
        
        if not self.transients_table:
            messagebox.showerror("Erreur", "Détectez d'abord des transitoires")
            return
        
        if not self.science_image_path or not Path(self.science_image_path).exists():
            messagebox.showerror("Erreur", "Image science manquante")
            return
        
        catalog = self.phot_catalog_var.get()
        filter_name = self.phot_filter_var.get()
        method = self.phot_method_var.get()
        fwhm = self.phot_fwhm_var.get()
        
        use_npoap = self.use_npoap_trans_var.get()
        npo_filt = (self.npoap_filter_override_var.get() or "").strip()
        self.progress["value"] = 0
        Thread(
            target=self._photometry_task,
            args=(
                self.science_image_path,
                self.transients_table,
                catalog,
                filter_name,
                method,
                fwhm,
                use_npoap,
                npo_filt,
            ),
            daemon=True,
        ).start()

    def _photometry_task(
        self,
        image_path: str,
        sources: Table,
        catalog: str,
        filter_name: str,
        method: str,
        fwhm: float,
        use_npoap: bool = True,
        npoap_obs_filter: str = "",
    ):
        try:
            self.progress["value"] = 20
            logger.info(
                f"Photométrie avec catalogue {catalog}, filtre {filter_name}, méthode {method}..."
            )

            photometry = self.stdpipe_wrapper.perform_photometry(
                image_path,
                sources,
                catalog=catalog,
                filter_name=filter_name,
                photometry_method=method,
                fwhm=fwhm if method == "psf" else None,
                use_npoap_transformation=use_npoap,
                npoap_obs_filter=npoap_obs_filter or None,
            )
            
            if photometry:
                self.transients_table = photometry  # Remplacer avec photométrie calibrée
                self.progress["value"] = 100
                messagebox.showinfo("Succès", "Photométrie effectuée avec succès!")
            else:
                self.progress["value"] = 0
                messagebox.showerror("Erreur", "Échec de la photométrie")
                
        except Exception as e:
            self.progress["value"] = 0
            logger.error(f"Erreur photométrie: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de la photométrie:\n{e}")

    def export_results(self):
        if not self.transients_table:
            messagebox.showerror("Erreur", "Aucun résultat à exporter")
            return

        output_dir = Path(self.science_image_path).parent if self.science_image_path else self.base_dir
        output_path = filedialog.asksaveasfilename(
            initialdir=output_dir,
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("FITS Table", "*.fits")]
        )
        
        if not output_path:
            return
        
        try:
            if output_path.endswith('.csv'):
                self.transients_table.write(output_path, format='csv', overwrite=True)
            else:
                self.transients_table.write(output_path, format='fits', overwrite=True)
            
            messagebox.showinfo("Succès", f"Résultats exportés:\n{output_path}")
            
        except Exception as e:
            logger.error(f"Erreur export: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors de l'export:\n{e}")
    
    def show_subtraction_images(self):
        """Affiche une fenêtre avec les images science, référence et soustraite côte à côte."""
        science_path = self.science_fits_var.get()
        ref_path = self.reference_fits_var.get()
        subtracted_path = self.subtracted_image_path
        
        # Vérifier qu'on a au moins l'image soustraite
        if not subtracted_path or not Path(subtracted_path).exists():
            messagebox.showwarning("Attention", 
                                 "Aucune image soustraite disponible.\n\n"
                                 "Lancez d'abord la soustraction d'images.")
            return
        
        # Créer une fenêtre Toplevel
        viz_window = tk.Toplevel(self)
        viz_window.title("Visualisation - Images Avant/Après Soustraction")
        viz_window.geometry("1200x500")
        
        # Créer la figure matplotlib avec 3 sous-graphiques
        fig = plt.Figure(figsize=(12, 4))
        
        # Charger et afficher les images
        try:
            # Image science
            if science_path and Path(science_path).exists():
                with fits.open(science_path) as hdul:
                    data_sci = hdul[0].data.astype(float)
                ax1 = fig.add_subplot(131)
                interval = ZScaleInterval()
                vmin, vmax = interval.get_limits(data_sci)
                ax1.imshow(data_sci, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
                ax1.set_title("Image Science", fontsize=10, fontweight='bold')
                ax1.set_xlabel("X (pixels)")
                ax1.set_ylabel("Y (pixels)")
                ax1.grid(True, alpha=0.3)
            else:
                ax1 = fig.add_subplot(131)
                ax1.text(0.5, 0.5, "Image Science\nnon disponible", 
                        ha='center', va='center', transform=ax1.transAxes, fontsize=12)
                ax1.set_title("Image Science", fontsize=10)
            
            # Image de référence
            if ref_path and Path(ref_path).exists():
                with fits.open(ref_path) as hdul:
                    data_ref = hdul[0].data.astype(float)
                ax2 = fig.add_subplot(132)
                interval = ZScaleInterval()
                vmin, vmax = interval.get_limits(data_ref)
                ax2.imshow(data_ref, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
                ax2.set_title("Image de Référence", fontsize=10, fontweight='bold')
                ax2.set_xlabel("X (pixels)")
                ax2.set_ylabel("Y (pixels)")
                ax2.grid(True, alpha=0.3)
            else:
                ax2 = fig.add_subplot(132)
                ax2.text(0.5, 0.5, "Image de Référence\nnon disponible", 
                        ha='center', va='center', transform=ax2.transAxes, fontsize=12)
                ax2.set_title("Image de Référence", fontsize=10)
            
            # Image soustraite
            with fits.open(subtracted_path) as hdul:
                data_sub = hdul[0].data.astype(float)
            ax3 = fig.add_subplot(133)
            interval = ZScaleInterval()
            vmin, vmax = interval.get_limits(data_sub)
            ax3.imshow(data_sub, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
            ax3.set_title("Image Soustraite", fontsize=10, fontweight='bold', color='green')
            ax3.set_xlabel("X (pixels)")
            ax3.set_ylabel("Y (pixels)")
            ax3.grid(True, alpha=0.3)
            
            # Ajuster l'espacement entre les sous-graphiques
            fig.tight_layout(pad=2.0)
            
            # Intégrer dans Tkinter
            canvas = FigureCanvasTkAgg(fig, viz_window)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            
            # Ajouter la barre d'outils matplotlib
            toolbar = NavigationToolbar2Tk(canvas, viz_window)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            
        except Exception as e:
            logger.error(f"Erreur lors de la visualisation des images: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Impossible d'afficher les images:\n{e}")
            viz_window.destroy()
    
    # ------------------------------------------------------------------
    # Recherche transitoires (Astro-COLIBRI, remplace TNS)
    # ------------------------------------------------------------------
    def create_transient_search_frame(self, parent):
        """Cadre de recherche de transitoires via l'API Astro-COLIBRI."""
        ttk.Label(parent, text="🔍 Recherche Astro-COLIBRI",
                  font=("Helvetica", 12, "bold")).pack(pady=(0, 10), anchor="w")

        config_frame = ttk.LabelFrame(parent, text="Configuration", padding=10)
        config_frame.pack(fill="x", pady=5)
        ttk.Label(config_frame, text="UID (requis, 100 req/jour):").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        self.colibri_uid_var = tk.StringVar()
        ttk.Entry(config_frame, textvariable=self.colibri_uid_var, width=30).grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(config_frame, text="💾 Sauvegarder", command=self.save_colibri_config).grid(row=1, column=0, columnspan=2, pady=5)
        self.load_colibri_config()

        search_frame = ttk.LabelFrame(parent, text="Recherche par dates", padding=10)
        search_frame.pack(fill="x", pady=5)
        ttk.Label(search_frame, text="Date début (YYYY-MM-DD):").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        self.colibri_date_start_var = tk.StringVar()
        ttk.Entry(search_frame, textvariable=self.colibri_date_start_var, width=25).grid(row=0, column=1, padx=5, pady=2)
        ttk.Label(search_frame, text="Date fin (YYYY-MM-DD):").grid(row=1, column=0, sticky="e", padx=5, pady=2)
        self.colibri_date_end_var = tk.StringVar()
        ttk.Entry(search_frame, textvariable=self.colibri_date_end_var, width=25).grid(row=1, column=1, padx=5, pady=2)
        ttk.Button(search_frame, text="🔍 Rechercher", command=self.search_colibri).grid(row=2, column=0, columnspan=2, pady=10)

        results_frame = ttk.LabelFrame(parent, text="Résultats", padding=10)
        results_frame.pack(fill="both", expand=True, pady=5)
        results_scroll = ttk.Scrollbar(results_frame)
        results_scroll.pack(side="right", fill="y")
        # Augmenter légèrement la hauteur de la zone de résultats (~+10 %)
        self.transient_results_listbox = tk.Listbox(results_frame, yscrollcommand=results_scroll.set, height=17)
        self.transient_results_listbox.pack(side="left", fill="both", expand=True)
        results_scroll.config(command=self.transient_results_listbox.yview)
        self.transient_results_listbox.bind("<Double-Button-1>", self.on_transient_result_selected)
        btn_frame = ttk.Frame(results_frame)
        btn_frame.pack(fill="x", pady=5)
        # Boutons empilés verticalement
        ttk.Button(btn_frame, text="📋 Détails", command=self.get_transient_details).pack(side="top", padx=2, pady=1, fill="x")
        ttk.Button(btn_frame, text="🔄 Actualiser", command=self.refresh_transient_results).pack(side="top", padx=2, pady=1, fill="x")
        ttk.Button(btn_frame, text="📁 NINA.json", command=self.export_transient_to_nina).pack(side="top", padx=2, pady=1, fill="x")

    def _ra_dec_deg_to_nina_input(self, ra_deg, dec_deg):
        """Convertit RA/Dec en degrés vers le format NINA InputCoordinates (sexagésimal)."""
        ra_h = float(ra_deg) % 360.0 / 15.0
        ra_hours = int(ra_h)
        ra_m = (ra_h - ra_hours) * 60.0
        ra_minutes = int(ra_m)
        ra_seconds = (ra_m - ra_minutes) * 60.0
        dec = float(dec_deg)
        negative_dec = dec < 0
        ad = abs(dec)
        dec_degrees = int(ad)
        dec_m = (ad - dec_degrees) * 60.0
        dec_minutes = int(dec_m)
        dec_seconds = (dec_m - dec_minutes) * 60.0
        return {
            "RAHours": ra_hours,
            "RAMinutes": ra_minutes,
            "RASeconds": round(ra_seconds, 4),
            "NegativeDec": negative_dec,
            "DecDegrees": -dec_degrees if negative_dec else dec_degrees,
            "DecMinutes": dec_minutes,
            "DecSeconds": round(dec_seconds, 4),
        }

    def export_transient_to_nina(self):
        """Crée un fichier NINA.json (format DeepSkyObjectContainer) pour la cible sélectionnée dans Documents\\N.I.N.A\\Targets."""
        sel = self.transient_results_listbox.curselection()
        if not sel or sel[0] >= len(self.transient_search_results):
            messagebox.showwarning("Attention", "Sélectionnez une cible dans la liste")
            return
        evt = self.transient_search_results[sel[0]]
        # Même nom que dans la liste résultats (source_name puis trigger_id), pas celui de la fenêtre détail
        name = evt.get("source_name") or evt.get("trigger_id") or "Transient"
        ra = evt.get("ra")
        dec = evt.get("dec")
        if ra is None or dec is None:
            messagebox.showwarning("Attention", "Coordonnées (RA/Dec) manquantes pour cette cible")
            return
        try:
            ra = float(ra)
            dec = float(dec)
        except (TypeError, ValueError):
            messagebox.showerror("Erreur", "Coordonnées invalides")
            return
        template_path = Path(__file__).parent.parent / "templates" / "nina_target_template.json"
        if not template_path.exists():
            template_path = Path.home() / "Downloads" / "AT 2026fbz.json"
        if not template_path.exists():
            messagebox.showerror("Erreur", "Template NINA introuvable.\nPlacez le fichier dans:\n%s" % (Path(__file__).parent.parent / "templates" / "nina_target_template.json"))
            return
        try:
            with open(template_path, "r", encoding="utf-8") as f:
                nina_data = json.load(f)
        except (OSError, ValueError) as e:
            messagebox.showerror("Erreur", "Impossible de lire le template NINA:\n%s" % e)
            return
        coords = self._ra_dec_deg_to_nina_input(ra, dec)
        if "Target" in nina_data and isinstance(nina_data["Target"], dict):
            nina_data["Target"]["TargetName"] = name
            if "InputCoordinates" in nina_data["Target"]:
                for k, v in coords.items():
                    nina_data["Target"]["InputCoordinates"][k] = v
        nina_data["Name"] = name
        items = nina_data.get("Items") or {}
        values = items.get("$values") or []
        if values and isinstance(values[0], dict) and "Coordinates" in values[0]:
            for k, v in coords.items():
                values[0]["Coordinates"][k] = v
        nina_dir = (Path.home() / "Documents" / "N.I.N.A" / "Targets").resolve()
        try:
            nina_dir.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            messagebox.showerror("Erreur", "Impossible de créer le répertoire NINA Targets:\n%s" % e)
            return
        safe_name = "".join(c if c.isalnum() or c in " ._-" else "_" for c in str(name)).strip() or "target"
        safe_name = safe_name.replace(" ", "_")
        filepath = (nina_dir / ("%s.json" % safe_name)).resolve()
        try:
            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(nina_data, f, indent=2, ensure_ascii=False)
            try:
                if os.name == "nt":
                    os.startfile(nina_dir)
                elif os.name == "posix":
                    import subprocess
                    subprocess.run(["xdg-open", str(nina_dir)], check=False)
            except Exception:
                pass
            messagebox.showinfo("Succès", "Fichier créé.\n\nDossier ouvert dans l'Explorateur.\n\nChemin complet:\n%s" % filepath)
        except OSError as e:
            messagebox.showerror("Erreur", "Impossible d'écrire le fichier:\n%s" % e)
            logger.error("Export NINA: %s", e)

    def load_colibri_config(self):
        """Charge l'UID Astro-COLIBRI depuis config."""
        try:
            import config
            uid = getattr(config, "ASTRO_COLIBRI_UID", "") or ""
            self.colibri_uid_var.set(uid)
        except Exception as e:
            logger.error("Erreur chargement config Astro-COLIBRI: %s", e)

    def save_colibri_config(self):
        """Sauvegarde l'UID Astro-COLIBRI dans config.py (persistant au redémarrage)."""
        try:
            import config
            import re
            uid = self.colibri_uid_var.get().strip()
            config_path = Path(__file__).parent.parent / "config.py"
            with open(config_path, "r", encoding="utf-8") as f:
                content = f.read()
            # Remplacer la ligne ASTRO_COLIBRI_UID = "..."
            pattern = r"^(ASTRO_COLIBRI_UID\s*=\s*).*$"
            replacement = r"\g<1>" + repr(uid)
            new_content = re.sub(pattern, replacement, content, count=1, flags=re.MULTILINE)
            if new_content == content:
                logger.warning("Ligne ASTRO_COLIBRI_UID non trouvée dans config.py")
            with open(config_path, "w", encoding="utf-8") as f:
                f.write(new_content)
            config.ASTRO_COLIBRI_UID = uid
            messagebox.showinfo("Succès", "UID Astro-COLIBRI sauvegardé dans config.py (rappelé à l'ouverture)")
        except Exception as e:
            logger.error("Erreur sauvegarde config Astro-COLIBRI: %s", e)
            messagebox.showerror("Erreur", "Impossible de sauvegarder:\n%s" % e)
    
    def search_colibri(self):
        """Lance une recherche Astro-COLIBRI par plage de dates (latest_transients)."""
        if not ASTRO_COLIBRI_AVAILABLE or not AstroColibriClient:
            messagebox.showerror("Erreur", "Module Astro-COLIBRI indisponible")
            return
        from datetime import datetime
        date_start = self.colibri_date_start_var.get().strip() or None
        date_end = self.colibri_date_end_var.get().strip() or None
        if not date_start or not date_end:
            messagebox.showwarning("Attention", "Indiquez la date de début et la date de fin (YYYY-MM-DD)")
            return
        try:
            datetime.strptime(date_start, "%Y-%m-%d")
            datetime.strptime(date_end, "%Y-%m-%d")
        except ValueError:
            messagebox.showerror("Erreur", "Format de date invalide. Utilisez YYYY-MM-DD")
            return
        uid = self.colibri_uid_var.get().strip() or None
        if not uid:
            messagebox.showwarning(
                "Attention",
                "Un UID Astro-COLIBRI est requis pour la recherche par dates (inscription gratuite sur astro-colibri.com, 100 req/jour)."
            )
            return
        Thread(
            target=self._search_colibri_task,
            args=(uid, date_start, date_end),
            daemon=True,
        ).start()

    def _search_colibri_task(self, uid, date_start, date_end):
        """Recherche Astro-COLIBRI par dates (latest_transients) en arrière-plan."""
        try:
            client = AstroColibriClient(uid=uid)
            time_min = "%sT00:00:00" % date_start
            time_max = "%sT23:59:59" % date_end
            results = client.latest_transients(time_min=time_min, time_max=time_max) or []
            self.transient_search_results = results
            self.after(0, self._display_colibri_results)
            if not results:
                self.after(0, lambda: messagebox.showinfo("Info", "Aucun résultat trouvé"))
        except Exception as e:
            logger.error("Erreur recherche Astro-COLIBRI: %s", e, exc_info=True)
            self.after(0, lambda: messagebox.showerror("Erreur", "Erreur recherche Astro-COLIBRI:\n%s" % e))

    def _display_colibri_results(self):
        """Affiche les résultats Astro-COLIBRI."""
        self.transient_results_listbox.delete(0, tk.END)
        for evt in self.transient_search_results:
            name = evt.get("source_name") 
            ra_val = evt.get("ra")
            dec_val = evt.get("dec")
            if ra_val is not None and dec_val is not None:
                coord_str = "RA:%.4f Dec:%.4f" % (float(ra_val), float(dec_val))
            else:
                coord_str = "RA:N/A Dec:N/A"
            t = evt.get("time") or evt.get("timestamp") or "N/A"
            if isinstance(t, (int, float)) and t > 0:
                try:
                    from datetime import datetime
                    t = datetime.utcfromtimestamp(t / 1000.0).strftime("%Y-%m-%d %H:%M") if t > 1e12 else datetime.utcfromtimestamp(t).strftime("%Y-%m-%d %H:%M")
                except Exception:
                    t = str(t)
            self.transient_results_listbox.insert(tk.END, "%s | %s | %s" % (name, coord_str, t))

    def refresh_transient_results(self):
        if self.transient_search_results:
            self._display_colibri_results()

    def on_transient_result_selected(self, event):
        self.get_transient_details()

    def get_transient_details(self):
        sel = self.transient_results_listbox.curselection()
        if not sel or sel[0] >= len(self.transient_search_results):
            messagebox.showwarning("Attention", "Sélectionnez un objet dans la liste")
            return
        evt = self.transient_search_results[sel[0]]
        tid = evt.get("trigger_id")
        name = evt.get("source_name")
        Thread(target=self._get_colibri_details_task, args=(tid, name), daemon=True).start()

    def _get_colibri_details_task(self, trigger_id, source_name):
        try:
            client = AstroColibriClient(uid=self.colibri_uid_var.get().strip() or None)
            details = client.get_event(trigger_id=trigger_id, source_name=source_name)
            if details:
                self.after(0, lambda: self._show_transient_details(details))
            else:
                self.after(0, lambda: messagebox.showinfo("Info", "Aucun détail trouvé"))
        except Exception as e:
            logger.error("Erreur détails Astro-COLIBRI: %s", e, exc_info=True)
            self.after(0, lambda: messagebox.showerror("Erreur", "Erreur récupération détails:\n%s" % e))

    def _deg_to_ra_hms(self, ra_deg):
        """Convertit RA (degrés) en chaîne hh:mm:ss."""
        if ra_deg is None:
            return "—"
        try:
            ra = float(ra_deg) % 360.0
            h = int(ra / 15.0)
            m = int((ra / 15.0 - h) * 60.0)
            s = ((ra / 15.0 - h) * 60.0 - m) * 60.0
            return "%02d:%02d:%06.3f" % (h, m, s)
        except (TypeError, ValueError):
            return str(ra_deg)

    def _deg_to_dec_dms(self, dec_deg):
        """Convertit Dec (degrés) en chaîne ±dd:mm:ss."""
        if dec_deg is None:
            return "—"
        try:
            d = float(dec_deg)
            sign = "+" if d >= 0 else "-"
            d = abs(d)
            dd = int(d)
            m = int((d - dd) * 60.0)
            s = ((d - dd) * 60.0 - m) * 60.0
            return "%s%02d:%02d:%05.2f" % (sign, dd, m, s)
        except (TypeError, ValueError):
            return str(dec_deg)

    def _time_from_obs(self, obs):
        """Extrait une chaîne temps d'un dict observation (time, date, mjd, timestamp)."""
        if not isinstance(obs, dict):
            return None
        for key in ("time", "date", "mjd", "timestamp", "obs_time", "obs_date"):
            val = obs.get(key)
            if val is not None and val != "" and val != -1:
                return str(val)
        return None

    def _mag_from_obs(self, obs):
        """Extrait la magnitude d'un dict observation (mag, magnitude, mag_value, etc.)."""
        if not isinstance(obs, dict):
            return None
        for key in ("mag", "magnitude", "mag_value", "mag_abs"):
            val = obs.get(key)
            if val is not None and val != "" and val != -1:
                try:
                    return "%s" % float(val)
                except (TypeError, ValueError):
                    return str(val)
        return None

    def _parse_time_for_comparison(self, t):
        """Retourne une valeur comparable (MJD ou timestamp) pour trier les temps, ou None."""
        if t is None:
            return None
        if isinstance(t, (int, float)) and t > 0:
            return float(t)  # mjd ou timestamp
        s = str(t).strip()
        if not s:
            return None
        # ISO / date string -> on garde tel quel pour affichage, tri approximatif par string
        if "T" in s or "-" in s:
            return s
        try:
            return float(s)
        except ValueError:
            return s

    def _photometry_last_magnitude_time(self, photometry):
        """
        Extrait le temps de la dernière magnitude (Astro-COLIBRI: first/peak/last par filtre).
        Structure typique: photometry[survey] = {"first": {...}, "peak": {...}, "last": {"time": ..., "mag": ...}}
        ou photometry = {"last": {"time": ...}}. On prend le temps le plus récent parmi tous les "last".
        """
        time_val = self._photometry_last_time_and_mag(photometry)
        return time_val[0] if time_val else "—"

    def _photometry_last_magnitude_value(self, photometry):
        """Extrait la valeur de la dernière magnitude (celle associée au dernier temps)."""
        time_val = self._photometry_last_time_and_mag(photometry)
        return time_val[1] if time_val else "—"

    def _photometry_last_time_and_mag(self, photometry):
        """
        Retourne (time_str, mag_str) pour l'observation "last" la plus récente, ou None.
        """
        if not photometry or not isinstance(photometry, dict):
            return None
        collected = []
        # Clé top-level "last" (une seule série)
        top_last = photometry.get("last")
        t = self._time_from_obs(top_last) if isinstance(top_last, dict) else None
        m = self._mag_from_obs(top_last) if isinstance(top_last, dict) else None
        if t:
            collected.append((self._parse_time_for_comparison(t), t, m or "—"))
        # Par survey/filtre: photometry[ZTF], photometry[ATLAS], etc. -> chacun peut avoir "last"
        for key, val in photometry.items():
            if key in ("first", "peak", "last") and val is top_last:
                continue
            if not isinstance(val, dict):
                continue
            last_obs = val.get("last")
            t = self._time_from_obs(last_obs) if isinstance(last_obs, dict) else None
            m = self._mag_from_obs(last_obs) if isinstance(last_obs, dict) else None
            if t:
                collected.append((self._parse_time_for_comparison(t), t, m or "—"))
        # Liste de points: dernier point
        if not collected and isinstance(photometry.get("points"), list) and photometry["points"]:
            last_pt = photometry["points"][-1]
            t = self._time_from_obs(last_pt)
            m = self._mag_from_obs(last_pt)
            if t:
                collected.append((self._parse_time_for_comparison(t), t, m or "—"))
        # dates/times + magnitudes list
        if not collected and isinstance(photometry.get("dates"), list) and photometry["dates"]:
            times = photometry["dates"]
            mags = photometry.get("magnitudes") or photometry.get("mags") or photometry.get("mag") or []
            if isinstance(mags, list) and mags and len(mags) >= len(times):
                last_val = times[-1]
                last_mag = mags[-1] if len(mags) > 0 else "—"
                collected.append((self._parse_time_for_comparison(last_val), str(last_val), str(last_mag)))
            elif times:
                last_val = times[-1]
                collected.append((self._parse_time_for_comparison(last_val), str(last_val), "—"))
        if not collected:
            return None
        def sort_key(item):
            comp, disp, _ = item
            if isinstance(comp, (int, float)):
                return (1, comp)
            return (0, str(comp))
        collected.sort(key=sort_key)
        best = collected[-1]
        return (best[1], best[2])

    def _show_transient_details(self, details):
        """Affiche une fenêtre détails : identification, classification, RA (hms), Dec (dms), photometry last magnitude time."""
        win = tk.Toplevel(self)
        name = details.get("source_name") or details.get("trigger_id") or "Transient"
        win.title("Détails Astro-COLIBRI: %s" % name)
        win.geometry("480x260")
        text_frame = ttk.Frame(win)
        text_frame.pack(fill="both", expand=True, padx=10, pady=10)
        identification = details.get("astro_colibri_id") or details.get("identifier") or details.get("source_name") or details.get("trigger_id") or "—"
        classification = details.get("classification") or "—"
        ra_str = self._deg_to_ra_hms(details.get("ra"))
        dec_str = self._deg_to_dec_dms(details.get("dec"))
        last_mag_time = self._photometry_last_magnitude_time(details.get("photometry"))
        last_mag_value = self._photometry_last_magnitude_value(details.get("photometry"))
        lines = [
            "Identification: %s" % identification,
            "Classification: %s" % classification,
            "RA (hh:mm:ss): %s" % ra_str,
            "Dec (dd:mm:ss): %s" % dec_str,
            "Dernière magnitude: %s" % last_mag_value,
            "Photometry last magnitude time: %s" % last_mag_time,
        ]
        text_widget = tk.Text(text_frame, height=9, wrap=tk.WORD, font=("Consolas", 10))
        text_widget.pack(fill="both", expand=True)
        text_widget.insert("1.0", "\n".join(lines))
        text_widget.config(state="disabled")