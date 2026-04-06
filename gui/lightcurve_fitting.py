# gui/lightcurve_fitting.py
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import logging
import pylightcurve as plc
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from pathlib import Path
from scipy.stats import shapiro
from statsmodels.tsa.stattools import acf
from scipy.optimize import minimize
from scipy.interpolate import interp1d

# Importer le module pour la loi power-2
try:
    from core.limb_darkening_power2 import (
        transit_lightcurve_power2,
        transit_lightcurve_quadratic,
        transit_lightcurve_square_root,
        fit_power2_coefficients,
        power2_intensity,
    )
    POWER2_AVAILABLE = True
except ImportError as e:
    POWER2_AVAILABLE = False
    logging.warning(f"Module limb_darkening_power2 non disponible: {e}")

try:
    import astroquery  # noqa: F401
    from core.claret_ld_vizier import nearest_power2_gaia, format_lookup_summary

    CLARET_2022_VIZIER_AVAILABLE = True
except ImportError:
    CLARET_2022_VIZIER_AVAILABLE = False

# Importer le module pour la récupération des priors
try:
    from core.exoplanet_priors_sources import (
        get_priors_from_all_sources,
        format_priors_for_display
    )
    PRIORS_AVAILABLE = True
except ImportError as e:
    PRIORS_AVAILABLE = False
    logging.warning(f"Module exoplanet_priors_sources non disponible: {e}")

# Importer le module pour les diagnostics de qualité
try:
    from core.quality_diagnostics import QualityDiagnostics
    QUALITY_DIAGNOSTICS_AVAILABLE = True
except ImportError as e:
    QUALITY_DIAGNOSTICS_AVAILABLE = False
    logging.warning(f"Module quality_diagnostics non disponible: {e}")

# Importer le module de binning temporel
try:
    from core.temporal_binning import bin_lightcurve, optimal_bin_time
    TEMPORAL_BINNING_AVAILABLE = True
except ImportError as e:
    TEMPORAL_BINNING_AVAILABLE = False
    logging.warning(f"Module temporal_binning non disponible: {e}")

from core.gaia_pylightcurve_support import (
    is_gaia_pylightcurve_filter,
    prepare_gaia_pylightcurve_transit,
)

try:
    from core.exotethys_ldc import (
        EXOTETHYS_AVAILABLE,
        EXOTETHYS_STELLAR_MODELS,
        inject_exotethys_claret4_into_planet,
        prepare_planet_exotethys_claret4,
        run_exotethys_ldc_claret4,
        run_exotethys_ldc_power2,
    )
except ImportError:
    EXOTETHYS_AVAILABLE = False
    EXOTETHYS_STELLAR_MODELS = ()
    run_exotethys_ldc_power2 = None
    run_exotethys_ldc_claret4 = None
    prepare_planet_exotethys_claret4 = None
    inject_exotethys_claret4_into_planet = None

logger = logging.getLogger(__name__)

# Classe ToolTip pour les tooltips
class ToolTip:
    """Classe pour créer des tooltips."""
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0
        self._id1 = self.widget.bind("<Enter>", self.enter)
        self._id2 = self.widget.bind("<Leave>", self.leave)
        self._id3 = self.widget.bind("<ButtonPress>", self.leave)

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(500, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert") if hasattr(self.widget, 'bbox') else (0, 0, 0, 0)
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                        background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                        font=("TkDefaultFont", 9), wraplength=300)
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

# --- Fonction de calcul mathématique locale ---
def perform_detrending(flux, vectors_dict):
    if not vectors_dict or len(flux) == 0:
        med = np.median(flux) if len(flux) > 0 else 1.0
        if med == 0: med = 1.0
        return flux / med, np.ones_like(flux) * med

    try:
        X_list = [np.ones_like(flux)]
        for vec in vectors_dict.values():
            if len(vec) == len(flux):
                X_list.append(vec)
        
        X = np.column_stack(X_list)
        beta, _, _, _ = np.linalg.lstsq(X, flux, rcond=None)
        model = X @ beta
        model[model == 0] = 1.0
        return flux / model, model
    except Exception:
        med = np.median(flux)
        return flux / med, np.ones_like(flux) * med

class LightcurveFitting(tk.Toplevel):
    def __init__(self, parent, cached_priors=None):
        super().__init__(parent)
        self.title("Analyse Courbe de Lumière & Modélisation Exoplanètes")
        self.geometry("1500x1000")

        self.data = None
        self.processed_data = None
        self.model_flux = None  # Modèle théorique (depuis BDD pylightcurve)
        self.fitted_model_flux = None  # Modèle ajusté aux données observées (après fitting)
        self.planet_obj = None  # Objet Planet pylightcurve pour extraire Rp/Rs
        self.excluded_indices = set()  # Indices des points exclus par l'utilisateur
        self._index_mapping = None  # Mapping: index processed_data -> index data original
        self.t0_reference = None  # T0 de référence (origine de l'époque) pour calcul O-C
        self.epoch_reference = None  # Époque de référence fixée pour éviter les sauts lors du détrend
        self.t0_fitted = None  # T0 ajusté après fitting (pour calcul O-C comme dans HOPS)
        self.rprs_fitted = None  # Rp/Rs ajusté après fitting (pour recalculer le modèle ajusté)
        self.mid_obs_reference = None  # mid_obs initial (avant ajustement) pour calcul O-C
        self.rprs_theoretical = None  # Rp/Rs théorique depuis la BDD (valeur fixe)
        self.manual_model_active = False  # True si le dernier modèle vient d'une saisie manuelle (hors catalogue)

        # Catalogue pylightcurve / NASA (décocher = planète non répertoriée, orbite saisie à la main)
        self.use_catalog_for_model = tk.BooleanVar(value=True)
        self.manual_rp_rs = tk.DoubleVar(value=0.1)
        self.manual_a_rs = tk.DoubleVar(value=8.0)
        self.manual_inclination_deg = tk.DoubleVar(value=87.0)

        # Variables de contrôle
        self.detrend_vars = {
            'AIRMASS': tk.BooleanVar(value=False),
            'FWHM_Mean': tk.BooleanVar(value=False),
            'Sky/Pixel_T1': tk.BooleanVar(value=False),
            'X_T1': tk.BooleanVar(value=False),
            'Y_T1': tk.BooleanVar(value=False)
        }
        
        self.sigma_clip_active = tk.BooleanVar(value=False)
        self.sigma_clip_value = tk.DoubleVar(value=0.0) 
        
        # Sliders
        self.obs_start = tk.DoubleVar(value=0.0)
        self.obs_end = tk.DoubleVar(value=0.0)
        self.transit_start = tk.DoubleVar(value=0.0)
        self.transit_end = tk.DoubleVar(value=0.0)
        self.transit_mid_calc = tk.StringVar(value="")

        # Modèle & Paramètres
        self.planet_name = tk.StringVar(value="")
        self.filter_name = tk.StringVar(value="")
        self.t0_var = tk.DoubleVar(value=0.0)
        self.period_var = tk.DoubleVar(value=0.0)
        
        self.teff_var = tk.DoubleVar(value=0.0) 
        self.logg_var = tk.DoubleVar(value=0.0)
        self.met_var = tk.DoubleVar(value=0.0)
        
        # Loi de limb-darkening
        self.limb_darkening_law = tk.StringVar(value="pylightcurve")  # Par défaut, utiliser pylightcurve
        self.power2_c = tk.DoubleVar(value=0.5)  # Coefficient c de la loi power-2
        self.power2_alpha = tk.DoubleVar(value=0.5)  # Coefficient α de la loi power-2
        self.claret_gaia_passband = tk.StringVar(value="G")  # G, G_BP, G_RP → Claret+ 2022 table1
        self.fit_limb_darkening = tk.BooleanVar(value=False)  # Ajuster les coefficients empiriquement
        # Source des LDC pour le transit pylightcurve : grille intégrée ou ExoTETHyS
        self.ldc_source_var = tk.StringVar(value="grid")
        self.exotethys_model_var = tk.StringVar(value="Phoenix_2018")
        self._exotethys_ldc_cache_key = None
        self._exotethys_ldc_cache_coeffs = None
        
        # Priors pour l'ajustement claret-4
        self.use_priors = tk.BooleanVar(value=False)  # Si True : a/R* et i des champs priors priment sur planet_obj (power-2, etc.)
        self.prior_a_rs = tk.DoubleVar(value=0.0)  # Prior sur a/R*
        self.prior_a_rs_err = tk.DoubleVar(value=0.0)  # Erreur sur a/R*
        self.prior_inclination = tk.DoubleVar(value=0.0)  # Prior sur l'inclinaison
        self.prior_inclination_err = tk.DoubleVar(value=0.0)  # Erreur sur l'inclinaison
        self.priors_source = tk.StringVar(value="")  # Source des priors
        self.priors_data = None  # Dictionnaire contenant les priors récupérés

        # Variables Qualité
        self.quality_vars = {
            "depth": tk.StringVar(value=""),
            "sigma_oot": tk.StringVar(value=""),
            "snr": tk.StringVar(value=""),
            "shapiro": tk.StringVar(value=""),
            "oc": tk.StringVar(value=""),
            "rprs": tk.StringVar(value=""),
            "beta": tk.StringVar(value=""),
            "chi2": tk.StringVar(value=""),
            "duration": tk.StringVar(value="")
        }

        self.stats_cache = {}
        
        # Variables pour le binning temporel
        self.apply_binning = tk.BooleanVar(value=False)
        self.binning_auto = tk.BooleanVar(value=True)
        self.exposure_time = tk.DoubleVar(value=60.0)  # secondes
        self.bin_time_manual = tk.DoubleVar(value=60.0)  # secondes
        self.binning_method = tk.StringVar(value="mean")
        self.preserve_transit_shape = tk.BooleanVar(value=True)
        
        # Diagnostics de qualité automatiques
        if QUALITY_DIAGNOSTICS_AVAILABLE:
            self.quality_diagnostics = QualityDiagnostics()
        else:
            self.quality_diagnostics = None
        
        self.create_widgets()

    def create_widgets(self):
        # Barre d'outils
        top = ttk.Frame(self, padding=5)
        top.pack(fill="x", side="top")
        ttk.Button(top, text="Charger CSV", command=self.load_csv).pack(side="left", padx=5)
        ttk.Button(top, text="Valider & Exporter (Report + PNG)", command=self.valider_modele).pack(side="left", padx=5, fill="x")

        # Panneau de droite avec scrollbar
        side_container = ttk.Frame(self, width=380)
        side_container.pack(side="right", fill="y", expand=False)
        side_container.pack_propagate(False)
        
        # Canvas et Scrollbar pour le panneau scrollable
        canvas_side = tk.Canvas(side_container, highlightthickness=0)
        scrollbar_side = ttk.Scrollbar(side_container, orient="vertical", command=canvas_side.yview)
        side = ttk.Frame(canvas_side, padding=12)
        # Définir une largeur minimale pour le frame interne
        side.update_idletasks()
        
        side_window = canvas_side.create_window((0, 0), window=side, anchor="nw")
        
        def configure_scroll_region(event):
            # Mettre à jour la région de scroll
            canvas_side.configure(scrollregion=canvas_side.bbox("all"))
            # Ajuster la largeur du frame interne à celle du canvas
            canvas_width = canvas_side.winfo_width()
            if canvas_width > 1:
                canvas_side.itemconfig(side_window, width=canvas_width)
        
        def configure_canvas_width(event):
            # Ajuster la largeur du frame interne à celle du canvas
            canvas_width = event.width
            canvas_side.itemconfig(side_window, width=canvas_width)
            # Mettre à jour la région de scroll
            canvas_side.configure(scrollregion=canvas_side.bbox("all"))
        
        side.bind("<Configure>", configure_scroll_region)
        canvas_side.bind("<Configure>", configure_canvas_width)
        canvas_side.configure(yscrollcommand=scrollbar_side.set)
        
        scrollbar_side.pack(side="right", fill="y")
        canvas_side.pack(side="left", fill="both", expand=True)
        
        # Activer le scroll avec la molette de la souris (Windows)
        def on_mousewheel(event):
            canvas_side.yview_scroll(int(-1 * (event.delta / 120)), "units")
        def bind_mousewheel(event):
            canvas_side.bind_all("<MouseWheel>", on_mousewheel)
        def unbind_mousewheel(event):
            canvas_side.unbind_all("<MouseWheel>")
        canvas_side.bind('<Enter>', bind_mousewheel)
        canvas_side.bind('<Leave>', unbind_mousewheel)

        # 0. Binning Temporel (NOUVEAU)
        if TEMPORAL_BINNING_AVAILABLE:
            lf_binning = ttk.LabelFrame(side, text="0. Binning Temporel (Optionnel)", padding=10)
            lf_binning.pack(fill="x", pady=(0, 10))
            
            ttk.Checkbutton(lf_binning, text="Appliquer binning temporel", 
                          variable=self.apply_binning).pack(anchor="w", pady=2)
            
            # Mode automatique ou manuel
            binning_mode_frame = ttk.Frame(lf_binning)
            binning_mode_frame.pack(fill="x", pady=2)
            ttk.Radiobutton(binning_mode_frame, text="Automatique", 
                          variable=self.binning_auto, value=True).pack(side="left", padx=5)
            ttk.Radiobutton(binning_mode_frame, text="Manuel", 
                          variable=self.binning_auto, value=False).pack(side="left", padx=5)
            
            # Temps d'exposition (pour automatique)
            exp_frame = ttk.Frame(lf_binning)
            exp_frame.pack(fill="x", pady=2)
            ttk.Label(exp_frame, text="Temps d'exposition (s):").pack(side="left")
            exp_entry = ttk.Entry(exp_frame, textvariable=self.exposure_time, width=10)
            exp_entry.pack(side="right")
            ToolTip(exp_entry, 
                   "Temps d'exposition en secondes.\n"
                   "Utilisé pour calculer automatiquement le binning optimal.\n"
                   "Recommandé: temps réel d'exposition de vos images.")
            
            # Temps de binning manuel
            bin_time_frame = ttk.Frame(lf_binning)
            bin_time_frame.pack(fill="x", pady=2)
            ttk.Label(bin_time_frame, text="Temps de binning (s):").pack(side="left")
            bin_entry = ttk.Entry(bin_time_frame, textvariable=self.bin_time_manual, width=10)
            bin_entry.pack(side="right")
            ToolTip(bin_entry,
                   "Temps de binning en secondes pour le mode manuel.\n"
                   "Typique: 60 s (comme Kepler). Plus long = moins de points mais meilleur SNR.")
            
            # Méthode de binning
            method_frame = ttk.Frame(lf_binning)
            method_frame.pack(fill="x", pady=2)
            ttk.Label(method_frame, text="Méthode:").pack(side="left")
            method_combo = ttk.Combobox(method_frame, textvariable=self.binning_method,
                        values=["mean", "median", "weighted"], width=12, state="readonly")
            method_combo.pack(side="right")
            ToolTip(method_combo,
                   "Méthode de binning:\n"
                   "- mean: moyenne (recommandé, rapide)\n"
                   "- median: médiane (plus robuste aux outliers)\n"
                   "- weighted: pondéré par erreurs (nécessite flux_err)")
            
            # Préserver la forme du transit
            ttk.Checkbutton(lf_binning, text="Binning adaptatif (fin pendant transit)", 
                          variable=self.preserve_transit_shape).pack(anchor="w", pady=2)
            
            # Info tooltip
            info_label = ttk.Label(lf_binning, 
                                  text="Binning optimal ~1 min (Kepler): réduit bruit\n"
                                       "sans affecter la précision des paramètres",
                                  font=("TkDefaultFont", 7), foreground="gray")
            info_label.pack(anchor="w", pady=2)
        
        # 1. Modélisation
        lf_model = ttk.LabelFrame(side, text="1. Modélisation", padding=10)
        lf_model.pack(fill="x", pady=(0, 10))
        
        ttk.Checkbutton(
            lf_model,
            text="Récupérer les paramètres depuis pylightcurve / NASA",
            variable=self.use_catalog_for_model,
            command=self._on_use_catalog_model_toggle,
        ).pack(anchor="w", pady=(0, 4))
        ttk.Label(
            lf_model,
            text="Décochez pour une planète hors catalogue : Rp/R*, a/R* et i sont saisis ci‑dessous ; "
                 "T₀ et la période viennent de la section « Infos Transit ».",
            font=("Arial", 8),
            foreground="gray",
            wraplength=320,
            justify="left",
        ).pack(anchor="w", pady=(0, 6))
        
        f_p = ttk.Frame(lf_model)
        f_p.pack(fill="x", pady=2)
        self.lbl_planet_name = ttk.Label(f_p, text="Planète :")
        self.lbl_planet_name.pack(side="left", padx=(0, 5))
        self.planet_entry = ttk.Entry(f_p, textvariable=self.planet_name, width=15)
        self.planet_entry.pack(side="right", fill="x", expand=True)
        
        f_f = ttk.Frame(lf_model)
        f_f.pack(fill="x", pady=2)
        ttk.Label(f_f, text="Filtre :").pack(side="left")
        ttk.Combobox(f_f, textvariable=self.filter_name, 
                     values=["COUSINS_V", "COUSINS_R", "COUSINS_I", "TESS", "SLOAN_g", "SLOAN_r", "SLOAN_i",
                             "G", "G_BP", "G_RP"],
                     width=12, state="readonly").pack(side="right")
        
        self.manual_orbit_frame = ttk.LabelFrame(
            lf_model, text="Orbite (saisie manuelle, hors catalogue)", padding=6
        )
        for lbl, var in [
            ("Rp/R* :", self.manual_rp_rs),
            ("a/R* :", self.manual_a_rs),
            ("Inclinaison (°) :", self.manual_inclination_deg),
        ]:
            row = ttk.Frame(self.manual_orbit_frame)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=lbl).pack(side="left")
            ttk.Entry(row, textvariable=var, width=12).pack(side="right")
        
        # Paramètres Stellaires
        lf_stellar = ttk.LabelFrame(lf_model, text="Paramètres Stellaires (Requis)", padding=5)
        self.lf_stellar_frame = lf_stellar
        lf_stellar.pack(fill="x", pady=5)
        
        for lbl, var in [("Teff (K):", self.teff_var), ("Log g:", self.logg_var), ("[Fe/H]:", self.met_var)]:
            f_s = ttk.Frame(lf_stellar)
            f_s.pack(fill="x")
            ttk.Label(f_s, text=lbl).pack(side="left")
            ttk.Entry(f_s, textvariable=var, width=8).pack(side="right")
        
        # Loi de Limb-Darkening
        lf_ld = ttk.LabelFrame(lf_model, text="Loi de Limb-Darkening", padding=5)
        lf_ld.pack(fill="x", pady=5)
        
        f_ld = ttk.Frame(lf_ld)
        f_ld.pack(fill="x", pady=2)
        ttk.Label(f_ld, text="Loi:").pack(side="left")
        ld_combo = ttk.Combobox(f_ld, textvariable=self.limb_darkening_law,
                               values=["pylightcurve", "power-2", "quadratique", "square-root"],
                               width=15, state="readonly")
        ld_combo.pack(side="right")
        ld_combo.set("pylightcurve")
        
        f_ld_src = ttk.Frame(lf_ld)
        f_ld_src.pack(fill="x", pady=2)
        ttk.Label(f_ld_src, text="LDC transit pylightcurve :").pack(side="left")
        self.ldc_source_combo = ttk.Combobox(
            f_ld_src,
            textvariable=self.ldc_source_var,
            values=("grid", "exotethys"),
            width=14,
            state="readonly",
        )
        self.ldc_source_combo.pack(side="right")
        self.ldc_source_combo.set("grid")
        ToolTip(
            self.ldc_source_combo,
            "grille = coefficients tabulés pylightcurve (rapide).\n"
            "exotethys = ExoTETHyS (Phoenix/ATLAS) pour le filtre choisi (plus lent, 1er lancement télécharge des données).",
        )
        f_exo_m = ttk.Frame(lf_ld)
        f_exo_m.pack(fill="x", pady=2)
        ttk.Label(f_exo_m, text="Modèle ExoTETHyS :").pack(side="left")
        self.exotethys_model_combo = ttk.Combobox(
            f_exo_m,
            textvariable=self.exotethys_model_var,
            values=list(EXOTETHYS_STELLAR_MODELS) if EXOTETHYS_STELLAR_MODELS else ["Phoenix_2018"],
            width=22,
            state="readonly" if EXOTETHYS_AVAILABLE else "disabled",
        )
        self.exotethys_model_combo.pack(side="right")
        if EXOTETHYS_STELLAR_MODELS:
            self.exotethys_model_combo.set("Phoenix_2018")
        
        # Paramètres Power-2 (affichés seulement si power-2 est sélectionné)
        self.power2_frame = ttk.Frame(lf_ld)
        self.power2_frame.pack(fill="x", pady=2)
        
        f_c = ttk.Frame(self.power2_frame)
        f_c.pack(fill="x")
        ttk.Label(f_c, text="c:").pack(side="left")
        ttk.Entry(f_c, textvariable=self.power2_c, width=8).pack(side="right")
        
        f_alpha = ttk.Frame(self.power2_frame)
        f_alpha.pack(fill="x", pady=2)
        ttk.Label(f_alpha, text="α:").pack(side="left")
        ttk.Entry(f_alpha, textvariable=self.power2_alpha, width=8).pack(side="right")

        if CLARET_2022_VIZIER_AVAILABLE:
            f_claret = ttk.Frame(self.power2_frame)
            f_claret.pack(fill="x", pady=4)
            ttk.Label(f_claret, text="Gaia (Claret+2022) :").pack(side="left")
            ttk.Combobox(
                f_claret,
                textvariable=self.claret_gaia_passband,
                values=["G", "G_BP", "G_RP"],
                width=8,
                state="readonly",
            ).pack(side="left", padx=4)
            ttk.Button(
                f_claret,
                text="Remplir c, α depuis VizieR",
                command=self.apply_claret_2022_power2_from_vizier,
            ).pack(side="left", padx=2)
            if EXOTETHYS_AVAILABLE and run_exotethys_ldc_power2 is not None:
                ttk.Button(
                    f_claret,
                    text="Remplir c, α depuis ExoTETHyS",
                    command=self.apply_exotethys_power2_from_filter,
                ).pack(side="left", padx=2)
            ttk.Label(
                f_claret,
                text="Teff, log g, [Fe/H] ci-dessus ; vt=2 km/s",
                font=("Arial", 7),
                foreground="gray",
            ).pack(side="left", padx=6)
        
        # Ajustement empirique des coeff. power-2 : visible uniquement si loi = power-2
        self.power2_empirical_frame = ttk.Frame(lf_ld)
        ttk.Checkbutton(
            self.power2_empirical_frame,
            text="Ajuster c et α empiriquement sur les données",
            variable=self.fit_limb_darkening,
        ).pack(anchor="w")
        ttk.Label(
            self.power2_empirical_frame,
            text=(
                "Uniquement pour la loi power-2 : le bouton « Ajuster Rp/Rs & T0 » optimise alors "
                "c et α en plus de Rp/Rs et T₀. Aucun effet avec pylightcurve, quadratique ou square-root."
            ),
            font=("Arial", 8),
            foreground="gray",
            wraplength=320,
            justify="left",
        ).pack(anchor="w", pady=(2, 0))
        
        # Bouton pour comparer les lois de limb-darkening (réf. pour réinsérer power-2 au bon endroit)
        self.ld_compare_laws_btn = ttk.Button(
            lf_ld, text="📊 Comparer les Lois",
            command=self.compare_limb_darkening_laws, width=20,
        )
        self.ld_compare_laws_btn.pack(fill="x", pady=5)
        
        # Affichage des coefficients ajustés
        self.ld_coefficients_label = ttk.Label(lf_ld, text="Coefficients: -", 
                                              foreground="gray", font=("Arial", 8))
        self.ld_coefficients_label.pack(anchor="w", pady=2)
        
        # Masquer le frame power-2 et l'option d'ajustement empirique par défaut
        self.power2_frame.pack_forget()
        
        # Callback pour afficher/masquer les paramètres power-2
        def on_ld_law_change(event=None):
            if self.limb_darkening_law.get() == "power-2":
                self.power2_frame.pack(fill="x", pady=2, before=self.ld_compare_laws_btn)
                self.power2_empirical_frame.pack(
                    fill="x", pady=(2, 4), before=self.ld_compare_laws_btn
                )
            else:
                self.power2_frame.pack_forget()
                self.power2_empirical_frame.pack_forget()
                self.fit_limb_darkening.set(False)
            # Mettre à jour l'affichage des coefficients
            self.update_ld_coefficients_display()
        
        ld_combo.bind("<<ComboboxSelected>>", on_ld_law_change)

        def on_ldc_source_change(event=None):
            self._exotethys_ldc_cache_key = None
            self.update_ld_coefficients_display()

        self.ldc_source_combo.bind("<<ComboboxSelected>>", on_ldc_source_change)
        self.exotethys_model_combo.bind("<<ComboboxSelected>>", on_ldc_source_change)
        
        # Géométrie orbitale (a/R*, i) — uniquement en mode catalogue : surcharge optionnelle pour power-2
        self.lf_geometry_transit_frame = None
        if PRIORS_AVAILABLE:
            lf_priors = ttk.LabelFrame(
                lf_model,
                text="Géométrie du transit : a/R* et inclinaison (optionnel)",
                padding=5,
            )
            self.lf_geometry_transit_frame = lf_priors
            lf_priors.pack(fill="x", pady=5)
            ttk.Label(
                lf_priors,
                text=(
                    "Réservé au mode catalogue : surcharge a/R* et i pour la loi power-2 et « Comparer les lois », "
                    "si la case ci-dessous est cochée (priorité sur le modèle planète). "
                    "En mode hors catalogue, utilisez uniquement le cadre « Orbite (saisie manuelle) » — "
                    "ce bloc est alors masqué. Indépendant du menu « Loi de Limb-Darkening » (limbe ≠ trajectoire)."
                ),
                font=("Arial", 8),
                foreground="gray",
                wraplength=320,
                justify="left",
            ).pack(anchor="w", pady=(0, 6))

            ttk.Checkbutton(
                lf_priors,
                text="Appliquer a/R* et i ci-dessous (priorité sur le modèle planète)",
                variable=self.use_priors,
            ).pack(anchor="w", pady=2)

            f_priors_btn = ttk.Frame(lf_priors)
            f_priors_btn.pack(fill="x", pady=2)
            ttk.Button(
                f_priors_btn,
                text="📥 Remplir depuis les archives",
                command=self.fetch_priors_from_sources,
                width=22,
            ).pack(side="left", padx=2)
            ttk.Button(
                f_priors_btn,
                text="✏️ Saisie manuelle a/R*, i",
                command=self.show_manual_priors_dialog,
                width=20,
            ).pack(side="left", padx=2)

            self.priors_label = ttk.Label(
                lf_priors,
                text="Aucune valeur — cochez « Appliquer… » après remplissage pour power-2",
                foreground="gray",
                font=("Arial", 8),
            )
            self.priors_label.pack(anchor="w", pady=2)

        self.btn_generate_model = ttk.Button(
            lf_model, text="Générer Modèle", command=self.calculate_model
        )
        self.btn_generate_model.pack(fill="x", pady=5)
        self.btn_fit_rprs_t0 = ttk.Button(
            lf_model, text="🔧 Ajuster Rp/Rs & T0", command=self.fit_model_parameters
        )
        self.btn_fit_rprs_t0.pack(fill="x", pady=5)

        # 2. Paramètres manuels
        lf_info = ttk.LabelFrame(side, text="2. Infos Transit", padding=10)
        lf_info.pack(fill="x", pady=10)

        f_t0 = ttk.Frame(lf_info)
        f_t0.pack(fill="x")
        ttk.Label(f_t0, text="T₀ (JD):").pack(side="left")
        ttk.Entry(f_t0, textvariable=self.t0_var, width=15).pack(side="right")

        f_per = ttk.Frame(lf_info)
        f_per.pack(fill="x", pady=2)
        ttk.Label(f_per, text="Période:").pack(side="left")
        ttk.Entry(f_per, textvariable=self.period_var, width=15).pack(side="right")

        ttk.Label(lf_info, textvariable=self.transit_mid_calc, foreground="blue", font=("Arial", 9, "bold")).pack(pady=5)

        # 3. Qualité
        lf_quality = ttk.LabelFrame(side, text="3. Qualité (Temps réel)", padding=10)
        lf_quality.pack(fill="x", pady=10)

        grid_f = ttk.Frame(lf_quality)
        grid_f.pack(fill="x")
        
        # Configurer les colonnes pour permettre l'expansion
        grid_f.columnconfigure(0, weight=0, minsize=90)
        grid_f.columnconfigure(1, weight=0, minsize=70)
        grid_f.columnconfigure(2, weight=1, minsize=80)
        
        # Normes attendues pour chaque critère
        norms = {
            "depth": "> 0.0001",
            "sigma_oot": "< 0.01",
            "snr": "> 5",
            "rprs": "> 0.01",
            "duration": "-",
            "beta": "< 1.2",
            "chi2": "~ 1.0",
            "shapiro": "> 0.05",
            "oc": "< 5 min"
        }
        
        labels = [
            ("Profondeur:", "depth"),
            ("σ OOT:", "sigma_oot"),
            ("SNR:", "snr"),
            ("Rp/Rs:", "rprs"),
            ("Durée (min):", "duration"),
            ("β (Red Noise):", "beta"),
            ("χ² réduit:", "chi2"),
            ("Shapiro p:", "shapiro"),
            ("O–C (min):", "oc")
        ]
        
        for i, (txt, var) in enumerate(labels):
            ttk.Label(grid_f, text=txt).grid(row=i, column=0, sticky="w", pady=1)
            ttk.Label(grid_f, textvariable=self.quality_vars[var], foreground="#333333").grid(row=i, column=1, sticky="e", padx=5)
            # Ajouter la norme attendue
            norm_text = norms.get(var, "-")
            ttk.Label(grid_f, text=f"({norm_text})", foreground="gray", font=("Arial", 8)).grid(row=i, column=2, sticky="w", padx=(5, 0))

        # 4. Diagnostics de Qualité Automatiques
        if QUALITY_DIAGNOSTICS_AVAILABLE:
            lf_diagnostics = ttk.LabelFrame(side, text="4. Diagnostics Automatiques", padding=10)
            lf_diagnostics.pack(fill="x", pady=10)
            
            # Zone de texte scrollable pour les diagnostics
            from tkinter.scrolledtext import ScrolledText
            self.diagnostics_text = ScrolledText(lf_diagnostics, height=8, width=40, 
                                                 wrap=tk.WORD, font=("Courier", 9),
                                                 bg="#f8f8f8", relief=tk.SUNKEN, borderwidth=1)
            self.diagnostics_text.pack(fill="x", pady=2)
            self.diagnostics_text.config(state=tk.DISABLED)  # Lecture seule
            
            # Bouton pour rafraîchir les diagnostics
            ttk.Button(lf_diagnostics, text="🔄 Rafraîchir Diagnostics", 
                      command=self.refresh_diagnostics_display).pack(fill="x", pady=2)
        else:
            # Si les diagnostics ne sont pas disponibles, afficher un message
            lf_diagnostics = ttk.LabelFrame(side, text="4. Diagnostics Automatiques", padding=10)
            lf_diagnostics.pack(fill="x", pady=10)
            ttk.Label(lf_diagnostics, text="Module quality_diagnostics non disponible", 
                     foreground="gray", font=("Arial", 8)).pack()

        # Détrendage
        lf_det = ttk.LabelFrame(side, text="Détrendage & Nettoyage", padding=8)
        lf_det.pack(fill="x", pady=10)
        
        for k, v in self.detrend_vars.items():
            ttk.Checkbutton(lf_det, text=k, variable=v, command=self.update_processing).pack(anchor="w")

        f_sig = ttk.Frame(lf_det)
        f_sig.pack(fill="x", pady=5)
        ttk.Checkbutton(f_sig, text="Sigma Clip", variable=self.sigma_clip_active, command=self.update_processing).pack(side="left")
        e_sig = ttk.Entry(f_sig, textvariable=self.sigma_clip_value, width=5)
        e_sig.pack(side="left", padx=5)
        e_sig.bind("<Return>", lambda e: self.update_processing())
        ttk.Label(f_sig, text="σ").pack(side="left")
        ttk.Button(f_sig, text="Go", width=4, command=self.update_processing).pack(side="left", padx=5)
        ttk.Button(f_sig, text="↩️ Reset points", command=self.reset_excluded_points).pack(side="left", padx=5)
        
        # Note sur la suppression de points
        info_label = ttk.Label(lf_det, text="💡 Astuce: Cliquez sur un point du graphique\npour le supprimer de l'analyse", 
                              font=("Arial", 8), foreground="gray", justify="left")
        info_label.pack(anchor="w", pady=(5, 0))

        # Zone Graphique
        center = ttk.Frame(self)
        center.pack(side="left", fill="both", expand=True)

        slider_frame = ttk.LabelFrame(center, text="Ajustement Zone Transit (In/Out)", padding=5)
        slider_frame.pack(fill="x", padx=10, pady=5)
        self.scale_start = tk.Scale(slider_frame, orient="horizontal", label="Début Transit", command=lambda v: self.on_slider_move())
        self.scale_start.pack(fill="x")
        self.scale_end = tk.Scale(slider_frame, orient="horizontal", label="Fin Transit", command=lambda v: self.on_slider_move())
        self.scale_end.pack(fill="x")
        
        obs_frame = ttk.LabelFrame(center, text="Fenêtre d'Observation (Zoom)", padding=5)
        obs_frame.pack(fill="x", padx=10, pady=5)
        self.obs_start_scale = tk.Scale(obs_frame, orient="horizontal", command=lambda v: self.on_obs_slider_move())
        self.obs_start_scale.pack(fill="x")
        self.obs_end_scale = tk.Scale(obs_frame, orient="horizontal", command=lambda v: self.on_obs_slider_move())
        self.obs_end_scale.pack(fill="x")
        
        # --- CONFIGURATION GRAPHIQUE AVEC SOUS-FENÊTRE (RESIDUALS) ---
        self.fig = plt.Figure(figsize=(9, 7), dpi=100)
        # Création de 2 axes partageant l'axe X : Main (haut) et Residuals (bas)
        gs = self.fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3, 1], hspace=0)
        self.ax_main = self.fig.add_subplot(gs[0])
        self.ax_res = self.fig.add_subplot(gs[1], sharex=self.ax_main)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=center)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)
        NavigationToolbar2Tk(self.canvas, center).pack(side="bottom", fill="x")

    def load_from_path(self, csv_path):
        try:
            self.data_path = Path(csv_path)
            self.data = pd.read_csv(csv_path)

            # Colonne temps : accepter JD-UTC, JD_UTC ou BJD-TDB et forcer le type float (évite erreur d'échelle)
            time_col = None
            for c in ('JD-UTC', 'JD_UTC', 'BJD-TDB'):
                if c in self.data.columns:
                    time_col = c
                    break
            if time_col is None:
                raise ValueError("Colonne de temps manquante (JD-UTC, JD_UTC ou BJD-TDB)")
            self.data['JD-UTC'] = pd.to_numeric(self.data[time_col], errors='coerce')
            if time_col != 'JD-UTC':
                self.data.drop(columns=[time_col], inplace=True, errors='ignore')

            if 'rel_flux_T1_fn' in self.data.columns:
                self.data['FLUX_RAW'] = self.data['rel_flux_T1_fn']
            elif 'rel_flux_T1' in self.data.columns:
                med = self.data['rel_flux_T1'].median() or 1.0
                self.data['FLUX_RAW'] = self.data['rel_flux_T1'] / med
            elif 'rel_flux' in self.data.columns:
                 self.data['FLUX_RAW'] = self.data['rel_flux']
            else:
                raise ValueError("Colonne de flux manquante")

            if 'rel_flux_err_T1' not in self.data.columns:
                self.data['rel_flux_err_T1'] = 0.0

            self.data.dropna(subset=['FLUX_RAW', 'JD-UTC'], inplace=True)
            if self.data.empty: raise ValueError("Fichier vide après nettoyage")
            
            # Réinitialiser l'index pour faciliter la gestion des exclusions
            self.data.reset_index(drop=True, inplace=True)
            # Réinitialiser les indices exclus lors du chargement d'un nouveau fichier
            self.excluded_indices = set()
            self.manual_model_active = False

            jd_min, jd_max = float(self.data['JD-UTC'].min()), float(self.data['JD-UTC'].max())
            
            self.transit_start.set(jd_min)
            self.transit_end.set(jd_max)
            self.obs_start.set(jd_min)
            self.obs_end.set(jd_max)
            
            for s in [self.scale_start, self.scale_end, self.obs_start_scale, self.obs_end_scale]:
                s.config(from_=jd_min, to=jd_max, resolution=0.0001)
            
            self.scale_start.set(jd_min)
            self.scale_end.set(jd_max)
            self.obs_start_scale.set(jd_min)
            self.obs_end_scale.set(jd_max)

            self.update_processing()
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de charger :\n{e}")

    def load_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV", "*.csv")])
        if path: self.load_from_path(path)

    def on_slider_move(self):
        self.transit_start.set(self.scale_start.get())
        self.transit_end.set(self.scale_end.get())
        self.update_processing()

    def on_obs_slider_move(self):
        self.obs_start.set(self.obs_start_scale.get())
        self.obs_end.set(self.obs_end_scale.get())
        self.update_processing()

    def calculate_quality_indicators(self):
        if self.processed_data is None or self.processed_data.empty: return

        df = self.processed_data
        t = df['JD-UTC'].values
        f = df['FLUX_FINAL'].values
        
        t_start, t_end = self.transit_start.get(), self.transit_end.get()
        if t_start < t_end:
            oot_mask = (t < t_start) | (t > t_end)
            int_mask = ~oot_mask
            # Calcul du vrai mid-time observé (méthode Kipping 2010: minimum du flux)
            if np.sum(int_mask) > 0:
                # Trouver l'index du point avec le flux minimum dans le transit
                int_indices = np.where(int_mask)[0]
                t_in_transit = t[int_indices]
                f_in_transit = f[int_indices]
                min_idx = np.argmin(f_in_transit)
                
                # Améliorer la précision par interpolation quadratique autour du minimum
                # Utiliser 3 points autour du minimum pour une interpolation plus précise
                if len(t_in_transit) >= 3 and 0 < min_idx < len(t_in_transit) - 1:
                    # Points autour du minimum (min_idx-1, min_idx, min_idx+1)
                    idx_start = max(0, min_idx - 1)
                    idx_end = min(len(t_in_transit), min_idx + 2)
                    t_local = t_in_transit[idx_start:idx_end]
                    f_local = f_in_transit[idx_start:idx_end]
                    
                    # Interpolation quadratique pour trouver le minimum exact
                    # f(t) = a*t^2 + b*t + c
                    # Le minimum est à t_min = -b/(2*a)
                    try:
                        # Fit polynomial d'ordre 2
                        coeffs = np.polyfit(t_local, f_local, 2)
                        # Le minimum d'un polynôme ax^2 + bx + c est à -b/(2a)
                        if abs(coeffs[0]) > 1e-10:  # Vérifier que a != 0
                            mid_obs = -coeffs[1] / (2 * coeffs[0])
                            # Vérifier que le résultat est dans la plage raisonnable
                            if mid_obs < t_local.min() or mid_obs > t_local.max():
                                # Si hors limites, utiliser le point minimum
                                mid_obs = t_in_transit[min_idx]
                        else:
                            mid_obs = t_in_transit[min_idx]
                    except:
                        # Si l'interpolation échoue, utiliser le point minimum
                        mid_obs = t_in_transit[min_idx]
                else:
                    # Pas assez de points pour interpolation, utiliser le point minimum
                    mid_obs = t_in_transit[min_idx]
            else:
                mid_obs = (t_start + t_end)/2
        else:
            oot_mask = np.ones_like(t, dtype=bool)
            int_mask = np.zeros_like(t, dtype=bool)
            mid_obs = np.mean(t)

        if np.sum(oot_mask) > 1:
            sigma_oot = np.std(f[oot_mask])
        else:
            sigma_oot = 0.0

        baseline = np.median(f[oot_mask]) if np.sum(oot_mask) > 0 else 1.0
        depth = baseline - np.min(f[int_mask]) if np.sum(int_mask) > 0 else 0.0
        depth = max(0.0, depth)
        snr = depth / sigma_oot if sigma_oot > 0 else 0.0
        
        # Rp/Rs théorique (depuis la BDD) - valeur fixe
        if self.planet_obj is not None and hasattr(self.planet_obj, 'rp_over_rs'):
            if self.rprs_theoretical is None:
                # Stocker la valeur théorique la première fois
                self.rprs_theoretical = float(self.planet_obj.rp_over_rs)
            rprs_theoretical = self.rprs_theoretical
        else:
            rprs_theoretical = None
        
        # Rp/Rs observé : calculé à partir de la profondeur du transit observé
        # Relation : depth ≈ (Rp/Rs)² pour un transit circulaire
        # Donc : Rp/Rs ≈ sqrt(depth)
        rprs_observed = np.sqrt(depth) if depth > 0 else 0.0
        
        # Utiliser la valeur observée pour l'affichage (change avec détrending et sélection)
        rprs = rprs_observed

        chi2 = np.nan
        beta = np.nan
        shapiro_p = np.nan
        oc = np.nan
        duration_min = np.nan
        acf_1 = np.nan
        res_mean = np.nan
        res_std = np.nan

        if self.model_flux is not None and len(self.model_flux) == len(f):
            resid = f - self.model_flux
            res_mean = np.mean(resid)
            
            # Residuals STD: Calculer seulement sur les points OOT (comme AIJ)
            if np.sum(oot_mask) > 1:
                res_std = np.std(resid[oot_mask])
            else:
                res_std = np.std(resid)
            
            # Extraire Rp/Rs du modèle si possible (via ajustement ou paramètre du modèle)
            # Pour l'instant, on utilise l'approximation, mais on pourrait améliorer avec un fit
            
            if sigma_oot > 0:
                # Chi2: utiliser seulement les points OOT pour le calcul de l'erreur
                chi2 = np.sum((resid / sigma_oot)**2) / max(1, len(resid) - 4)
            
            # Shapiro: Calculer seulement sur les résidus OOT (comme AIJ)
            try:
                if np.sum(oot_mask) > 3:
                    _, shapiro_p = shapiro(resid[oot_mask])
                elif len(resid) > 3:
                    _, shapiro_p = shapiro(resid)
            except: pass
            
            try:
                # ACF: Calculer sur tous les résidus
                acf_vals = acf(resid, nlags=1)
                if len(acf_vals) > 1: acf_1 = acf_vals[1]
            except: pass

            try:
                # Utiliser le T0 de référence (origine de l'époque) pour calculer O-C correctement
                # Méthode HOPS: comparer le T0 ajusté (ou mid_obs observé) avec la prédiction théorique
                if self.t0_reference is not None:
                    T0_ref = self.t0_reference
                else:
                    # Fallback sur t0_var si t0_reference n'est pas défini
                    T0_ref = self.t0_var.get()
                
                P = self.period_var.get()
                if P > 0 and T0_ref is not None:
                    # O-C = T0(JD) affiché dans l'interface - Mid affiché dans l'interface
                    # T0(JD) affiché = self.t0_var.get()
                    # Mid affiché = valeur extraite de transit_mid_calc
                    t0_displayed = self.t0_var.get()
                    
                    # Extraire la valeur Mid de transit_mid_calc (format: "Mid: X.XXXXX")
                    mid_displayed_str = self.transit_mid_calc.get()
                    try:
                        # Extraire le nombre après "Mid: "
                        mid_displayed = float(mid_displayed_str.replace("Mid: ", "").strip())
                    except (ValueError, AttributeError):
                        # Si extraction échoue, utiliser mid_obs comme fallback
                        mid_displayed = mid_obs
                    
                    # O-C = T0(JD) affiché - Mid affiché
                    oc = (t0_displayed - mid_displayed) * 24 * 60
                    
                    logger.info(f"O-C (calculé depuis valeurs affichées): T0(JD)={t0_displayed:.10f}, Mid={mid_displayed:.10f}, O-C={oc:.2f} min")
                    
                    # Durée: utiliser le modèle pour déterminer les limites du transit
                    in_transit_model = self.model_flux < 0.9999
                    if np.sum(in_transit_model) > 1:
                        t_mod_in = t[in_transit_model]
                        duration_min = (t_mod_in.max() - t_mod_in.min()) * 24 * 60
                    
                    # Essayer d'extraire Rp/Rs du modèle pylightcurve si disponible
                    # Pour l'instant, on garde l'approximation, mais on pourrait améliorer
                    # en ajustant le modèle pour trouver le meilleur Rp/Rs
            except: pass
            
            try:
                dt = np.median(np.diff(np.sort(t))) * 86400
                bin_size = int(600 / dt) 
                if bin_size > 1 and len(resid) > 2*bin_size:
                    n_bins = len(resid) // bin_size
                    binned = resid[:n_bins*bin_size].reshape(-1, bin_size).mean(axis=1)
                    sigma_n = np.std(binned)
                    sigma_1 = np.std(resid)
                    beta = sigma_n / (sigma_1 / np.sqrt(bin_size))
            except: pass

        self.quality_vars["depth"].set(f"{depth:.5f}")
        self.quality_vars["sigma_oot"].set(f"{sigma_oot:.6f}")
        self.quality_vars["snr"].set(f"{snr:.1f}")
        # Afficher Rp/Rs observé (change avec détrending et sélection)
        # Si valeur théorique disponible, l'afficher aussi pour comparaison
        if rprs_theoretical is not None:
            self.quality_vars["rprs"].set(f"{rprs:.4f} (théorique: {rprs_theoretical:.4f})")
        else:
            self.quality_vars["rprs"].set(f"{rprs:.4f}")
        self.quality_vars["duration"].set(f"{duration_min:.1f}" if not np.isnan(duration_min) else "-")
        self.quality_vars["chi2"].set(f"{chi2:.2f}" if not np.isnan(chi2) else "-")
        self.quality_vars["shapiro"].set(f"{shapiro_p:.3f}" if not np.isnan(shapiro_p) else "-")
        self.quality_vars["beta"].set(f"{beta:.2f}" if not np.isnan(beta) else "-")
        self.quality_vars["oc"].set(f"{oc:.1f}" if not np.isnan(oc) else "-")

        self.stats_cache = {
            "depth": depth, "snr": snr, "sigma_oot": sigma_oot, 
            "chi2": chi2, "beta": beta, "oc": oc, "rprs": rprs,
            "duration": duration_min, "shapiro_p": shapiro_p, "acf_1": acf_1,
            "res_mean": res_mean, "res_std": res_std
        }
        
        # Diagnostics de qualité automatiques
        if self.quality_diagnostics is not None:
            self.quality_diagnostics.clear()
            
            # Valider la profondeur
            self.quality_diagnostics.validate_transit_depth(depth)
            
            # Valider les paramètres orbitaux si disponibles
            if self.planet_obj is not None:
                period = self.period_var.get()
                if period > 0:
                    # Récupérer a/R* et inclinaison depuis l'objet Planet
                    a_rs = None
                    if hasattr(self.planet_obj, 'a_over_rs'):
                        a_rs = float(self.planet_obj.a_over_rs)
                    elif hasattr(self.planet_obj, 'semi_major_axis') and hasattr(self.planet_obj, 'stellar_radius'):
                        if self.planet_obj.semi_major_axis > 0 and self.planet_obj.stellar_radius > 0:
                            a_rs = float(self.planet_obj.semi_major_axis) / float(self.planet_obj.stellar_radius)
                    
                    inclination = None
                    if hasattr(self.planet_obj, 'inclination'):
                        inclination = float(self.planet_obj.inclination)
                    elif hasattr(self.planet_obj, 'orbital_inclination'):
                        inclination = float(self.planet_obj.orbital_inclination)
                    
                    if a_rs is not None and inclination is not None:
                        self.quality_diagnostics.validate_orbital_parameters(period, a_rs, inclination)
            
            # Vérifier la qualité des résidus
            if self.model_flux is not None and len(self.model_flux) == len(f):
                resid = f - self.model_flux
                self.quality_diagnostics.check_residuals_quality(resid, sigma_oot)
            
            # Vérifier le chi2
            if not np.isnan(chi2):
                n_data = len(f)
                n_params = 4  # Rp/Rs, T0, et éventuellement coefficients limb-darkening
                self.quality_diagnostics.check_chi2_quality(chi2, n_data, n_params)
            
            # Mettre à jour l'affichage graphique des diagnostics
            self.update_diagnostics_display()
            
            # Afficher les diagnostics dans les logs si des problèmes sont détectés
            if self.quality_diagnostics.warnings or self.quality_diagnostics.errors:
                report = self.quality_diagnostics.generate_report()
                logger.warning(f"Diagnostics de qualité:\n{report}")

    def update_processing(self):
        if self.data is None or self.data.empty: return

        df = self.data.copy()
        
        # Construire le mapping des indices: après chaque étape de filtrage
        # index_mapping[new_idx] = original_idx_in_data
        index_mapping = {i: i for i in range(len(df))}
        
        # Supprimer les points exclus par l'utilisateur
        if self.excluded_indices:
            valid_mask = ~df.index.isin(list(self.excluded_indices))
            # Créer nouveau mapping après exclusion
            old_mapping = index_mapping.copy()
            index_mapping = {}
            new_idx = 0
            for old_idx, orig_idx in old_mapping.items():
                if valid_mask[old_idx]:
                    index_mapping[new_idx] = orig_idx
                    new_idx += 1
            df = df[valid_mask].copy()
            df.reset_index(drop=True, inplace=True)
        
        jd_full = df['JD-UTC'].values
        obs_s, obs_e = self.obs_start.get(), self.obs_end.get()
        if obs_s < obs_e:
            mask_obs = (jd_full >= obs_s) & (jd_full <= obs_e)
            # Mettre à jour le mapping après filtrage par fenêtre d'observation
            old_mapping = index_mapping.copy()
            index_mapping = {}
            new_idx = 0
            for old_idx, orig_idx in old_mapping.items():
                if mask_obs[old_idx]:
                    index_mapping[new_idx] = orig_idx
                    new_idx += 1
            df = df[mask_obs].copy()
            df.reset_index(drop=True, inplace=True)
        
        if df.empty: return

        jd = df['JD-UTC'].values
        flux = df['FLUX_RAW'].values
        
        t_s, t_e = self.transit_start.get(), self.transit_end.get()
        oot_mask = (jd < t_s) | (jd > t_e) if t_s < t_e else np.ones_like(jd, dtype=bool)

        if self.sigma_clip_active.get() and np.sum(oot_mask) > 5:
            sigma_val = self.sigma_clip_value.get()
            if sigma_val > 0:
                med = np.median(flux[oot_mask])
                std = np.std(flux[oot_mask])
                if std > 0:
                    mask = np.abs(flux - med) < (sigma_val * std)
                    # Mettre à jour le mapping après sigma clip
                    old_mapping = index_mapping.copy()
                    index_mapping = {}
                    new_idx = 0
                    for old_idx, orig_idx in old_mapping.items():
                        if mask[old_idx]:
                            index_mapping[new_idx] = orig_idx
                            new_idx += 1
                    df = df[mask].copy()
                    df.reset_index(drop=True, inplace=True)
                    jd = df['JD-UTC'].values
                    flux = df['FLUX_RAW'].values
                    oot_mask = (jd < t_s) | (jd > t_e) if t_s < t_e else np.ones_like(jd, dtype=bool)

        # BINNING TEMPOREL (si activé) - appliqué après filtrage mais avant détrending
        if TEMPORAL_BINNING_AVAILABLE and self.apply_binning.get():
            try:
                # Déterminer le temps de binning
                if self.binning_auto.get():
                    exposure = self.exposure_time.get()
                    bin_time = optimal_bin_time(exposure, cadence=None, target_bin_time=60.0)
                    logger.info(f"Binning automatique: temps d'exposition={exposure}s, binning optimal={bin_time:.1f}s")
                else:
                    bin_time = self.bin_time_manual.get()
                    logger.info(f"Binning manuel: {bin_time}s")
                
                # Récupérer les erreurs si disponibles
                flux_err = None
                if 'rel_flux_err_T1' in df.columns:
                    flux_err = df['rel_flux_err_T1'].values
                
                # Durée du transit pour binning adaptatif (si disponible)
                transit_duration = None
                if self.preserve_transit_shape.get() and t_s < t_e:
                    transit_duration = (t_e - t_s)  # Déjà en jours (JD)
                    logger.info(f"Binning adaptatif: durée transit={(transit_duration * 86400.0):.1f}s")
                
                # Appliquer le binning
                jd_binned, flux_binned, flux_err_binned = bin_lightcurve(
                    jd, flux, flux_err,
                    bin_time=bin_time,
                    method=self.binning_method.get(),
                    preserve_transit=self.preserve_transit_shape.get() and transit_duration is not None,
                    transit_duration=transit_duration if transit_duration else None
                )
                
                # Mettre à jour les données
                jd = jd_binned
                flux = flux_binned
                
                # Recalculer oot_mask avec les nouvelles données
                oot_mask = (jd < t_s) | (jd > t_e) if t_s < t_e else np.ones_like(jd, dtype=bool)
                
                # Mettre à jour le DataFrame
                df_binned = pd.DataFrame({
                    'JD-UTC': jd,
                    'FLUX_RAW': flux
                })
                
                # Mettre à jour les erreurs si disponibles
                if flux_err_binned is not None:
                    df_binned['rel_flux_err_T1'] = flux_err_binned
                elif 'rel_flux_err_T1' in df.columns:
                    # Réestimer les erreurs depuis les données binnées
                    if len(flux) > 1:
                        df_binned['rel_flux_err_T1'] = np.std(flux) / np.sqrt(len(flux))
                    else:
                        df_binned['rel_flux_err_T1'] = 0.0
                
                # Copier/interpoler les autres colonnes nécessaires
                for col in df.columns:
                    if col not in df_binned.columns:
                        # Colonnes continues: interpoler
                        if col in ['AIRMASS', 'X_IMAGE', 'Y_IMAGE', 'FWHM']:
                            from scipy.interpolate import interp1d
                            if len(df[col].dropna()) > 1:
                                try:
                                    f_interp = interp1d(df['JD-UTC'].values, df[col].values, 
                                                       kind='linear', fill_value='extrapolate', 
                                                       bounds_error=False)
                                    df_binned[col] = f_interp(jd)
                                except:
                                    # Si interpolation échoue, prendre la valeur moyenne
                                    df_binned[col] = df[col].mean()
                            else:
                                df_binned[col] = df[col].iloc[0] if len(df) > 0 else np.nan
                        else:
                            # Autres colonnes: prendre la première valeur du bin
                            df_binned[col] = df[col].iloc[0] if len(df) > 0 else np.nan
                
                df = df_binned.copy()
                df.reset_index(drop=True, inplace=True)
                
                # Mettre à jour le mapping après binning (simplifié)
                n_original = len(index_mapping)
                n_binned = len(df)
                if n_binned > 0:
                    # Mapping approximatif: distribuer les indices originaux sur les binnés
                    index_mapping = {i: int(i * n_original / n_binned) for i in range(n_binned)}
                else:
                    index_mapping = {}
                
                logger.info(f"Binning appliqué: {n_binned} points (réduction de {n_original} à {n_binned})")
                
            except Exception as e:
                logger.warning(f"Erreur lors du binning temporel: {e}", exc_info=True)
                # Continuer avec les données non binnées en cas d'erreur

        vectors = {k: df[k].values for k, v in self.detrend_vars.items() if v.get() and k in df.columns}
        flux_clean = flux
        detrend_divisor = np.ones_like(flux, dtype=float)
        if vectors and np.sum(oot_mask) > 5:
            try:
                X_oot = np.column_stack([np.ones_like(flux[oot_mask])] + [v[oot_mask] for v in vectors.values()])
                beta = np.linalg.lstsq(X_oot, flux[oot_mask], rcond=None)[0]
                X_all = np.column_stack([np.ones_like(flux)] + list(vectors.values()))
                M = X_all @ beta
                M[np.abs(M) < 1e-15] = 1.0
                detrend_divisor = M
                flux_clean = flux / M
            except Exception:
                pass

        df['FLUX_FINAL'] = flux_clean
        err_raw = (
            df["rel_flux_err_T1"].values
            if "rel_flux_err_T1" in df.columns
            else np.zeros(len(flux), dtype=float)
        )
        # Détrending multiplicatif f' = f / M : σ_f' ≈ σ_f / |M| (M supposé sans incertitude)
        df["rel_flux_fn_detrend_err"] = err_raw / np.maximum(np.abs(detrend_divisor), 1e-15)
        self.processed_data = df
        # Stocker le mapping final: index dans processed_data -> index dans data original
        self._index_mapping = index_mapping
        
        # Recalculer le modèle théorique si nécessaire (pour correspondre aux nouvelles données)
        if self.planet_obj is not None and self.t0_reference is not None:
            try:
                flt = self.filter_name.get()
                if flt:
                    t_data = df['JD-UTC'].values
                    pl = self.planet_obj
                    period = getattr(pl, 'period', getattr(pl, 'orbital_period', None))
                    
                    # Recalculer le modèle théorique avec les paramètres originaux de la BDD
                    depth_theo = None
                    rprs_theo = None
                    if period is not None and self.t0_reference is not None:
                        # Utiliser Rp/Rs théorique si disponible, sinon utiliser la valeur actuelle de planet_obj
                        rprs_theo = self.rprs_theoretical if self.rprs_theoretical is not None else getattr(pl, 'rp_over_rs', 0.1)
                        
                        # Calculer le modèle théorique avec les paramètres de la BDD
                        pl.rp_over_rs = rprs_theo
                        logger.info(f"[RECALCUL] Modèle théorique: Rp/Rs={rprs_theo:.5f} (théorique depuis BDD)")
                        
                        epoch = round((t_data.min() - self.t0_reference) / period)
                        local_t0 = self.t0_reference + epoch * period
                        for attr in ['mid_time', 'transit_mid_time', 't0', 'tmid', 'transit_time']:
                            if hasattr(pl, attr):
                                setattr(pl, attr, local_t0)
                        
                        if not self._ldc_preparation_available(flt):
                            logger.warning(
                                "Recalcul du modèle ignoré : LDC indisponible pour ce filtre "
                                "(Gaia sans VizieR ni ExoTETHyS)."
                            )
                        else:
                            self._prepare_ldc_for_transit(pl, flt)
                            self.model_flux = pl.transit_integrated(
                                time=t_data,
                                time_format='JD_UTC',
                                exp_time=120,
                                time_stamp='mid',
                                max_sub_exp_time=1,
                                filter_name=flt
                            )
                            depth_theo = 1.0 - np.min(self.model_flux)
                            logger.info(
                                f"[RECALCUL] Modèle théorique calculé: Rp/Rs={pl.rp_over_rs:.5f}, "
                                f"T0={local_t0:.8f}, profondeur={depth_theo:.6f}"
                            )
                    
                    # Recalculer le modèle ajusté si un fitting a été fait
                    if self.rprs_fitted is not None and self.t0_fitted is not None:
                        pl.rp_over_rs = self.rprs_fitted
                        logger.info(f"[RECALCUL] Modèle ajusté: Rp/Rs={self.rprs_fitted:.5f} (ajusté), T0={self.t0_fitted:.8f}")
                        for attr in ['mid_time', 'transit_mid_time', 't0', 'tmid', 'transit_time']:
                            if hasattr(pl, attr):
                                setattr(pl, attr, self.t0_fitted)
                        if not self._ldc_preparation_available(flt):
                            logger.warning(
                                "Recalcul du modèle ajusté ignoré : LDC indisponible pour ce filtre."
                            )
                        else:
                            self._prepare_ldc_for_transit(pl, flt)
                            self.fitted_model_flux = pl.transit_integrated(
                                time=t_data,
                                time_format='JD_UTC',
                                exp_time=120,
                                time_stamp='mid',
                                max_sub_exp_time=1,
                                filter_name=flt
                            )
                            depth_fitted = 1.0 - np.min(self.fitted_model_flux)
                            logger.info(f"[RECALCUL] Modèle ajusté calculé: Rp/Rs={pl.rp_over_rs:.5f}, profondeur={depth_fitted:.6f}")
                            if depth_theo is not None and rprs_theo is not None:
                                logger.info(
                                    f"[COMPARAISON] Théorique: Rp/Rs={rprs_theo:.5f}, profondeur={depth_theo:.6f} "
                                    f"| Ajusté: Rp/Rs={self.rprs_fitted:.5f}, profondeur={depth_fitted:.6f}"
                                )
                    else:
                        self.fitted_model_flux = None
                    
            except Exception as e:
                logger.warning(f"Erreur lors du recalcul des modèles: {e}")
                import traceback
                logger.warning(traceback.format_exc())
                self.fitted_model_flux = None
        else:
            logger.debug(f"Recalcul modèles ignoré: planet_obj={self.planet_obj is not None}, t0_reference={self.t0_reference is not None}")
        
        self.calculate_quality_indicators()
        self.draw_plot()

    def draw_plot(self):
        if self.processed_data is None: return
        
        # Nettoyage
        self.ax_main.clear()
        self.ax_res.clear()
        
        t = self.processed_data['JD-UTC'].values
        f = self.processed_data['FLUX_FINAL'].values
        e = self.processed_data['rel_flux_err_T1'].values
        
        t_s, t_e = self.transit_start.get(), self.transit_end.get()
        in_transit = (t >= t_s) & (t <= t_e) if t_s < t_e else np.zeros_like(t, dtype=bool)

        # Stocker les indices pour permettre la suppression par clic
        # On stocke les scatter plots pour pouvoir les interroger
        self.scatter_main_oot = None
        self.scatter_main_in = None
        self.scatter_res = None
        self.data_indices = {}  # Dictionnaire pour mapper les indices des points affichés
        
        # --- PANEL HAUT : COURBE DE LUMIERE ---
        if np.any(~in_transit):
            oot_indices = np.where(~in_transit)[0]
            self.scatter_main_oot = self.ax_main.scatter(t[~in_transit], f[~in_transit], c='blue', s=15, alpha=0.6, picker=True, pickradius=5, label='Hors Transit')
            # Mapper les indices des points scatter aux indices du DataFrame original (self.data)
            scatter_dict = {}
            for i, local_idx in enumerate(oot_indices):
                # local_idx est l'index dans processed_data (après filtres)
                # On le mappe vers l'index dans data original via _index_mapping
                if hasattr(self, '_index_mapping') and local_idx in self._index_mapping:
                    scatter_dict[i] = self._index_mapping[local_idx]
                else:
                    scatter_dict[i] = local_idx  # Fallback
            self.data_indices[id(self.scatter_main_oot)] = scatter_dict
        
        if np.any(in_transit):
            in_indices = np.where(in_transit)[0]
            self.scatter_main_in = self.ax_main.scatter(t[in_transit], f[in_transit], c='red', s=15, alpha=0.8, picker=True, pickradius=5, label='Transit')
            # Mapper les indices
            scatter_dict = {}
            for i, local_idx in enumerate(in_indices):
                if hasattr(self, '_index_mapping') and local_idx in self._index_mapping:
                    scatter_dict[i] = self._index_mapping[local_idx]
                else:
                    scatter_dict[i] = local_idx
            self.data_indices[id(self.scatter_main_in)] = scatter_dict

        # Tracer la courbe lissée des observations
        if len(t) > 3:
            # Trier les données par temps
            idx_sorted = np.argsort(t)
            t_sorted = t[idx_sorted]
            f_sorted = f[idx_sorted]
            
            # Ne plus tracer de ligne interpolée bleue pointillée
            # Le modèle ajusté (rouge) sera affiché après le fitting

        # Modèle théorique (catalogue ou saisie manuelle) en bleu
        if self.model_flux is not None and len(self.model_flux) == len(t):
            idx = np.argsort(t)
            theo_lbl = (
                "Modèle théorique (manuel)"
                if getattr(self, "manual_model_active", False)
                else "Modèle théorique (catalogue)"
            )
            self.ax_main.plot(t[idx], self.model_flux[idx], "b-", lw=2, label=theo_lbl)
        
        # Modèle ajusté (après fitting) en rouge
        if self.fitted_model_flux is not None and len(self.fitted_model_flux) == len(t):
            idx = np.argsort(t)
            self.ax_main.plot(t[idx], self.fitted_model_flux[idx], 'r-', lw=2, label='Modèle ajusté')
        
        # --- PANEL BAS : RESIDUS ---
        # Utiliser le modèle ajusté s'il existe, sinon le modèle théorique
        if self.model_flux is not None and len(self.model_flux) == len(t):
            model_for_residuals = self.fitted_model_flux if (self.fitted_model_flux is not None and len(self.fitted_model_flux) == len(t)) else self.model_flux
            residuals = f - model_for_residuals
            self.scatter_res = self.ax_res.scatter(t, residuals, c='black', s=15, alpha=0.6, picker=True, pickradius=5)
            # Mapper tous les indices pour les résidus
            if hasattr(self, '_index_mapping'):
                scatter_dict = {i: self._index_mapping.get(i, i) for i in range(len(t))}
            else:
                scatter_dict = {i: i for i in range(len(t))}
            self.data_indices[id(self.scatter_res)] = scatter_dict
            self.ax_res.axhline(0, color='red', linestyle='--', alpha=0.5)

        # Zone Visuelle Transit
        if t_s < t_e:
            mid = (t_s + t_e)/2
            self.transit_mid_calc.set(f"Mid: {mid:.5f}")
            self.ax_main.axvspan(t_s, t_e, color='red', alpha=0.1)
            self.ax_main.axvline(mid, color='red', ls=':', alpha=0.5)

        # Styles
        self.ax_main.axhline(1.0, color='black', ls='--', alpha=0.5)
        
        pn = self.planet_name.get().strip()
        if pn:
            plt_title = f"Lightcurve: {pn}"
        elif getattr(self, "manual_model_active", False):
            plt_title = "Lightcurve (modèle manuel)"
        else:
            plt_title = "Lightcurve"
        self.ax_main.set_title(plt_title, fontsize=12, fontweight="bold")
        self.ax_main.set_ylabel("Norm. Flux")
        self.ax_main.legend(loc='best')
        
        # Axe X partagé
        plt.setp(self.ax_main.get_xticklabels(), visible=False)
        
        self.ax_res.set_xlabel("JD")
        self.ax_res.set_ylabel("Residuals")
        self.ax_res.grid(True, alpha=0.3)
        
        # Connecter l'événement de clic
        self.fig.canvas.mpl_connect('pick_event', self.on_point_pick)

        self.fig.tight_layout()
        self.canvas.draw()
    
    def on_point_pick(self, event):
        """Gère le clic sur un point pour le supprimer"""
        if event.artist not in [self.scatter_main_oot, self.scatter_main_in, self.scatter_res]:
            return
        
        if event.mouseevent.button != 1:  # Seulement clic gauche
            return
        
        # Trouver l'index du point cliqué
        scatter_id = id(event.artist)
        if scatter_id not in self.data_indices:
            return
        
        ind = event.ind[0]  # Premier point sélectionné
        if ind not in self.data_indices[scatter_id]:
            return
        
        # Obtenir l'index original dans le DataFrame
        original_idx = self.data_indices[scatter_id][ind]
        
        # Ajouter à la liste des exclus
        if original_idx not in self.excluded_indices:
            self.excluded_indices.add(original_idx)
            # Recalculer avec les nouveaux points exclus
            self.update_processing()
            logger.info(f"Point {original_idx} exclu. Total exclus: {len(self.excluded_indices)}")

    def reset_excluded_points(self):
        """Restaure tous les points exclus par l'utilisateur"""
        if self.excluded_indices:
            count = len(self.excluded_indices)
            self.excluded_indices.clear()
            self.update_processing()
            logger.info(f"{count} point(s) restauré(s)")
            messagebox.showinfo("Points restaurés", f"{count} point(s) restauré(s)")
        else:
            messagebox.showinfo("Info", "Aucun point exclu à restaurer")
    
    def update_diagnostics_display(self):
        """Met à jour l'affichage graphique des diagnostics de qualité."""
        if not QUALITY_DIAGNOSTICS_AVAILABLE or self.quality_diagnostics is None:
            return
        
        if not hasattr(self, 'diagnostics_text'):
            return
        
        # Générer le rapport
        report = self.quality_diagnostics.generate_report()
        
        # Activer l'édition, effacer et insérer le nouveau texte
        self.diagnostics_text.config(state=tk.NORMAL)
        self.diagnostics_text.delete(1.0, tk.END)
        
        # Insérer le rapport avec formatage
        self.diagnostics_text.insert(tk.END, report)
        
        # Appliquer des couleurs selon le type de message
        # Chercher les lignes avec des avertissements ou erreurs
        content = self.diagnostics_text.get(1.0, tk.END)
        lines = content.split('\n')
        
        # Configurer les tags pour les couleurs
        self.diagnostics_text.tag_config("info", foreground="green")
        self.diagnostics_text.tag_config("warning", foreground="orange")
        self.diagnostics_text.tag_config("error", foreground="red")
        self.diagnostics_text.tag_config("header", foreground="blue", font=("Courier", 9, "bold"))
        
        # Réinsérer avec formatage
        self.diagnostics_text.delete(1.0, tk.END)
        
        for line in lines:
            if "RAPPORT DE DIAGNOSTIC" in line or "="*60 in line:
                self.diagnostics_text.insert(tk.END, line + "\n", "header")
            elif "✓ INFORMATIONS:" in line or (line.startswith("  •") and "✓" in line):
                self.diagnostics_text.insert(tk.END, line + "\n", "info")
            elif "⚠ AVERTISSEMENTS:" in line or (line.startswith("  •") and "⚠" in line):
                self.diagnostics_text.insert(tk.END, line + "\n", "warning")
            elif "❌ ERREURS:" in line or (line.startswith("  •") and "❌" in line):
                self.diagnostics_text.insert(tk.END, line + "\n", "error")
            else:
                self.diagnostics_text.insert(tk.END, line + "\n")
        
        # Désactiver l'édition
        self.diagnostics_text.config(state=tk.DISABLED)
    
    def refresh_diagnostics_display(self):
        """Rafraîchit l'affichage des diagnostics (bouton manuel)."""
        if self.quality_diagnostics is not None:
            # Les diagnostics sont déjà calculés dans calculate_quality_indicators
            # On peut juste rafraîchir l'affichage
            self.update_diagnostics_display()
        else:
            if hasattr(self, 'diagnostics_text'):
                self.diagnostics_text.config(state=tk.NORMAL)
                self.diagnostics_text.delete(1.0, tk.END)
                self.diagnostics_text.insert(tk.END, "Aucun diagnostic disponible.\nChargez des données et générez un modèle.")
                self.diagnostics_text.config(state=tk.DISABLED)

    def get_planet_from_exoclock(self, target_name):
        """
        Récupère les paramètres de la planète depuis ExoClock.
        ExoClock fournit des éphémérides mises à jour régulièrement.
        
        Note: ExoClock n'a pas d'API publique directe, mais les données
        sont accessibles via NASA Exoplanet Archive qui intègre parfois
        des données ExoClock. Cette fonction est une structure pour
        future intégration si une API devient disponible.
        """
        try:
            # Tentative 1: Via astroquery si ExoClock est supporté
            # (pour l'instant, ExoClock n'est pas directement dans astroquery)
            
            # Tentative 2: Via NASA Archive qui peut inclure des données ExoClock
            # Pour l'instant, on retourne None et on laisse NASA Archive gérer
            # Cette fonction est une place-holder pour future intégration
            
            # TODO: Si ExoClock fournit une API ou un endpoint, l'utiliser ici
            # Exemple de structure future:
            # import requests
            # response = requests.get(f"https://exoclock.space/api/planet/{target_name}")
            # data = response.json()
            # ...
            
            logger.debug(f"ExoClock: pas d'API directe disponible pour {target_name}")
            return None
        except Exception as e:
            logger.debug(f"Erreur lors de la récupération ExoClock pour {target_name}: {e}")
            return None

    def get_planet_from_nasa(self, target_name):
        """
        Récupère les paramètres de la planète depuis NASA Exoplanet Archive
        en utilisant la table pscomppars qui contient les valeurs les plus précises.
        
        La table pscomppars (Planetary Systems Composite Parameters) est la table
        composite qui contient déjà les meilleures valeurs disponibles, équivalent
        au flag "most precise" du service Transit et Éphémérides.
        """
        try:
            # Utiliser la table pscomppars qui contient les valeurs les plus précises
            # Cette table est équivalente au flag "most precise" du service web
            table = NasaExoplanetArchive.query_object(target_name, table="pscomppars")
            if len(table) == 0: 
                logger.warning(f"Aucune donnée trouvée dans pscomppars pour {target_name}")
                return None
            
            row = table[0]
            
            # Vérifier que les paramètres essentiels sont disponibles
            if np.ma.is_masked(row['pl_orbper']) or np.ma.is_masked(row['pl_tranmid']): 
                logger.warning(f"Période ou T0 manquants pour {target_name}")
                return None
            
            per = float(row['pl_orbper'])
            t0 = float(row['pl_tranmid'])
            
            # Récupérer les paramètres avec des valeurs par défaut si manquants
            rp_rs = float(row['pl_ratror']) if not np.ma.is_masked(row['pl_ratror']) else 0.0
            a_rs = float(row['pl_ratdor']) if not np.ma.is_masked(row['pl_ratdor']) else 0.0
            inc = float(row['pl_orbincl']) if not np.ma.is_masked(row['pl_orbincl']) else 0.0
            teff = float(row['st_teff']) if not np.ma.is_masked(row['st_teff']) else 0.0
            logg = float(row['st_logg']) if not np.ma.is_masked(row['st_logg']) else 0.0
            met = float(row['st_met']) if not np.ma.is_masked(row['st_met']) else 0.0

            logger.info(f"Paramètres récupérés depuis NASA Archive (pscomppars - most precise) pour {target_name}: "
                       f"P={per:.6f}, T0={t0:.6f}, Rp/Rs={rp_rs:.4f}")

            pl = plc.Planet(
                name=target_name, ra=plc.Degrees(0), dec=plc.Degrees(0), 
                period=per, mid_time=t0,
                rp_over_rs=rp_rs, sma_over_rs=a_rs, inclination=inc,
                eccentricity=0, periastron=0, 
                stellar_temperature=teff, 
                stellar_logg=logg, 
                stellar_metallicity=met,
                mid_time_format='JD_UTC'
            )
            pl.period = per
            pl.mid_time = t0
            return pl
        except Exception as e:
            logger.error(f"Erreur lors de la récupération depuis NASA Archive pour {target_name}: {e}")
            return None

    def _on_use_catalog_model_toggle(self):
        """Affiche ou masque la saisie d'orbite manuelle selon le mode catalogue."""
        if self.use_catalog_for_model.get():
            self.manual_orbit_frame.pack_forget()
            self.lbl_planet_name.config(text="Planète :")
            if getattr(self, "lf_geometry_transit_frame", None) is not None:
                self.lf_geometry_transit_frame.pack(
                    fill="x", pady=5, before=self.btn_generate_model
                )
        else:
            self.manual_orbit_frame.pack(fill="x", pady=6, before=self.lf_stellar_frame)
            self.lbl_planet_name.config(text="Nom d'affichage (optionnel) :")
            if getattr(self, "lf_geometry_transit_frame", None) is not None:
                self.lf_geometry_transit_frame.pack_forget()

    def _fetch_planet_from_catalogs(self):
        """Construit un objet pylightcurve.Planet depuis pylightcurve, ExoClock (placeholder) ou NASA."""
        tgt = self.planet_name.get().strip()
        if not tgt:
            messagebox.showwarning(
                "Attention",
                "Entrez le nom de la planète pour la recherche catalogue,\n"
                "ou décochez « Récupérer les paramètres depuis pylightcurve / NASA » pour un modèle manuel.",
            )
            return None
        pl = None
        try:
            pl = plc.get_planet(tgt)
            logger.info(f"Paramètres récupérés depuis pylightcurve pour {tgt}")
        except Exception:
            try:
                pl = self.get_planet_from_exoclock(tgt)
                if pl:
                    logger.info(f"Paramètres récupérés depuis ExoClock pour {tgt}")
            except Exception:
                pass
            if not pl:
                try:
                    pl = self.get_planet_from_nasa(tgt)
                    if pl:
                        logger.info(f"Paramètres récupérés depuis NASA Exoplanet Archive pour {tgt}")
                except Exception:
                    pass
        if not pl:
            messagebox.showerror(
                "Planète introuvable",
                f"Aucune entrée pour « {tgt} » dans pylightcurve ni la NASA Archive.\n\n"
                "Décochez « Récupérer les paramètres depuis pylightcurve / NASA » pour saisir "
                "Rp/R*, a/R*, i, T₀ et la période à la main.",
            )
        return pl

    def _build_planet_manual_model(self):
        """Construit un Planet pylightcurve entièrement depuis l'interface (cible non cataloguée)."""
        period = float(self.period_var.get())
        t0_ref = float(self.t0_var.get())
        if period <= 0:
            messagebox.showerror("Erreur", "Indiquez une période > 0 dans « Infos Transit ».")
            return None
        if t0_ref < 2_000_000:
            messagebox.showerror(
                "Erreur",
                "Indiquez un T₀ (JD) réaliste dans « Infos Transit » (ex. ~2460000…), "
                "aligné sur l'éphéméride du transit pour cette fenêtre d'observation.",
            )
            return None
        teff = float(self.teff_var.get())
        logg = float(self.logg_var.get())
        met = float(self.met_var.get())
        if teff <= 0 or logg <= 0:
            messagebox.showerror(
                "Erreur",
                "Teff (K) et log g doivent être valides pour le calcul du limbe (pylightcurve).",
            )
            return None
        try:
            rp_rs = float(self.manual_rp_rs.get())
            a_rs = float(self.manual_a_rs.get())
            inc = float(self.manual_inclination_deg.get())
        except (ValueError, TypeError):
            messagebox.showerror("Erreur", "Rp/R*, a/R* et inclinaison doivent être numériques.")
            return None
        if not (0 < rp_rs < 1) or a_rs <= 0 or not (0 < inc < 90.01):
            messagebox.showerror(
                "Erreur",
                "Contraintes : 0 < Rp/R* < 1, a/R* > 0, inclinaison dans ]0 ; 90]° (degrés).",
            )
            return None
        name = self.planet_name.get().strip() or "Cible_manuelle"
        try:
            pl = plc.Planet(
                name=name,
                ra=plc.Degrees(0),
                dec=plc.Degrees(0),
                period=period,
                mid_time=t0_ref,
                rp_over_rs=rp_rs,
                sma_over_rs=a_rs,
                inclination=inc,
                eccentricity=0,
                periastron=0,
                stellar_temperature=teff,
                stellar_logg=logg,
                stellar_metallicity=met,
                mid_time_format="JD_UTC",
            )
            pl.period = period
            pl.mid_time = t0_ref
            logger.info(
                f"Modèle manuel : {name} P={period} T0={t0_ref} Rp/R*={rp_rs} a/R*={a_rs} i={inc}°"
            )
            return pl
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de construire le modèle pylightcurve :\n{e}")
            logger.exception("_build_planet_manual_model")
            return None

    def _ldc_preparation_available(self, flt: str) -> bool:
        """Indique si un recalcul transit pylightcurve peut obtenir des LDC pour ce filtre."""
        if EXOTETHYS_AVAILABLE and self.ldc_source_var.get() == "exotethys":
            return True
        if is_gaia_pylightcurve_filter(flt):
            return bool(CLARET_2022_VIZIER_AVAILABLE)
        return True

    def _prepare_ldc_for_transit(self, pl, flt: str) -> None:
        """Prépare passbands + LDC avant ``transit_integrated`` (grille, Gaia/VizieR ou ExoTETHyS)."""
        if pl is None or not flt:
            return
        if not hasattr(pl, "_npoap_ldc_original"):
            pl._npoap_ldc_original = pl.ldc_method

        if self.ldc_source_var.get() == "exotethys":
            if (
                not EXOTETHYS_AVAILABLE
                or run_exotethys_ldc_claret4 is None
                or inject_exotethys_claret4_into_planet is None
            ):
                raise RuntimeError(
                    "ExoTETHyS n'est pas installé. Exécutez : pip install exotethys h5py click"
                )
            teff = float(self.teff_var.get())
            logg = float(self.logg_var.get())
            met = float(self.met_var.get())
            model = self.exotethys_model_var.get()
            key = (flt, model, teff, logg, met)
            if self._exotethys_ldc_cache_key != key:
                self._exotethys_ldc_cache_coeffs = run_exotethys_ldc_claret4(
                    flt, teff, logg, met, stellar_model=model
                )
                self._exotethys_ldc_cache_key = key
            inject_exotethys_claret4_into_planet(pl, flt, self._exotethys_ldc_cache_coeffs)
            return

        if is_gaia_pylightcurve_filter(flt) and CLARET_2022_VIZIER_AVAILABLE:
            prepare_gaia_pylightcurve_transit(pl, flt)
            return

        pl.filters.pop(flt, None)
        pl.ldc_method = getattr(pl, "_npoap_ldc_original", pl.ldc_method)

    def calculate_model(self):
        if self.processed_data is None:
            messagebox.showwarning("Info", "Chargez un CSV d'abord")
            return

        flt = self.filter_name.get()
        if not flt:
            messagebox.showwarning("Attention", "Spécifiez un filtre.")
            return
        if self.ldc_source_var.get() == "exotethys" and not EXOTETHYS_AVAILABLE:
            messagebox.showerror(
                "ExoTETHyS",
                "Le mode « exotethys » est sélectionné mais le paquet n'est pas installé.\n"
                "Installez : pip install exotethys h5py click",
            )
            return
        if is_gaia_pylightcurve_filter(flt) and not self._ldc_preparation_available(flt):
            messagebox.showerror(
                "Filtre Gaia",
                "Filtres G / G_BP / G_RP : choisissez « exotethys » comme source LDC, ou installez "
                "astroquery pour VizieR (Claret & Southworth 2022).",
            )
            return

        try:
            self.config(cursor="watch")
            if self.use_catalog_for_model.get():
                self.manual_model_active = False
                pl = self._fetch_planet_from_catalogs()
                if pl is None:
                    return
                if hasattr(pl, "stellar_temperature") and pl.stellar_temperature > 0:
                    self.teff_var.set(pl.stellar_temperature)
                if hasattr(pl, "stellar_logg") and pl.stellar_logg > 0:
                    self.logg_var.set(pl.stellar_logg)
                if hasattr(pl, "stellar_metallicity"):
                    self.met_var.set(pl.stellar_metallicity)
            else:
                self.manual_model_active = True
                pl = self._build_planet_manual_model()
                if pl is None:
                    return

            pl.stellar_temperature = self.teff_var.get()
            pl.stellar_logg = self.logg_var.get()
            pl.stellar_metallicity = self.met_var.get()

            t0 = None
            for attr in ["mid_time", "transit_mid_time", "t0", "tmid", "transit_time"]:
                if hasattr(pl, attr):
                    t0 = getattr(pl, attr)
                    break
            if t0 is None:
                raise ValueError("Attribut T0 introuvable sur l'objet planète")

            period = getattr(pl, "period", getattr(pl, "orbital_period", None))
            if period is None:
                raise ValueError("Période introuvable")

            self.period_var.set(period)

            self.t0_reference = float(t0)
            self.epoch_reference = None
            self.t0_fitted = None
            self.rprs_fitted = None
            self.mid_obs_reference = None
            self.rprs_theoretical = None
            self.fitted_model_flux = None

            t_data = self.processed_data["JD-UTC"].values
            epoch = round((t_data.min() - t0) / period)
            local_t0 = t0 + epoch * period
            self.t0_var.set(local_t0)

            self._prepare_ldc_for_transit(pl, flt)
            self.model_flux = pl.transit_integrated(
                time=t_data,
                time_format="JD_UTC",
                exp_time=120,
                time_stamp="mid",
                max_sub_exp_time=1,
                filter_name=flt,
            )

            self.planet_obj = pl
            if hasattr(pl, "rp_over_rs"):
                self.rprs_theoretical = float(pl.rp_over_rs)

            self.update_processing()
            tgt = self.planet_name.get().strip()
            if self.manual_model_active:
                msg = f"Modèle manuel généré ({tgt or 'sans nom'})"
            else:
                msg = f"Modèle {tgt} généré (catalogue)"
            messagebox.showinfo("OK", msg)

        except Exception as e:
            messagebox.showerror("Erreur", str(e))
        finally:
            self.config(cursor="")

    def _get_power2_orbital_params(self):
        """Retourne a/R* et inclination pour power-2 (et comparaison des lois)."""
        a_rs = None
        inclination = None

        # Mode hors catalogue : une seule source — le cadre « Orbite (saisie manuelle) »
        if not self.use_catalog_for_model.get():
            try:
                a_m = float(self.manual_a_rs.get())
                i_m = float(self.manual_inclination_deg.get())
                if a_m > 0 and np.isfinite(a_m) and np.isfinite(i_m) and 0 < i_m <= 90.01:
                    logger.info(
                        f"[power-2] a/R* et i depuis saisie manuelle (hors catalogue): "
                        f"a/R*={a_m:.4f}, i={i_m:.4f}°"
                    )
                    return a_m, i_m
            except (tk.TclError, TypeError, ValueError):
                pass

        if self.use_priors.get():
            try:
                prior_a = float(self.prior_a_rs.get())
                if prior_a > 0 and np.isfinite(prior_a):
                    a_rs = prior_a
            except Exception:
                pass
            try:
                prior_i = float(self.prior_inclination.get())
                if prior_i > 0 and np.isfinite(prior_i):
                    inclination = prior_i
            except Exception:
                pass

        if a_rs is None or not np.isfinite(a_rs) or a_rs <= 0:
            if hasattr(self.planet_obj, 'a_over_rs') and self.planet_obj.a_over_rs:
                a_rs = float(self.planet_obj.a_over_rs)
            elif hasattr(self.planet_obj, 'semi_major_axis') and hasattr(self.planet_obj, 'stellar_radius'):
                try:
                    a_rs = float(self.planet_obj.semi_major_axis) / float(self.planet_obj.stellar_radius)
                except Exception:
                    a_rs = None

        if inclination is None or not np.isfinite(inclination) or inclination <= 0:
            if hasattr(self.planet_obj, 'inclination') and self.planet_obj.inclination:
                inclination = float(self.planet_obj.inclination)
            elif hasattr(self.planet_obj, 'orbital_inclination') and self.planet_obj.orbital_inclination:
                inclination = float(self.planet_obj.orbital_inclination)

        # Heuristique: si l'inclinaison est en radians, convertir en degrés
        if inclination is not None and np.isfinite(inclination) and 0 < inclination <= (np.pi + 0.01):
            inclination = np.degrees(inclination)

        if a_rs is None or not np.isfinite(a_rs) or a_rs <= 0:
            a_rs = 10.0
        if inclination is None or not np.isfinite(inclination) or inclination <= 0:
            inclination = 90.0

        logger.info(f"[power-2] Paramètres orbitaux: a/R*={a_rs:.4f}, i={inclination:.4f}°")
        return a_rs, inclination

    def fit_model_parameters(self):
        """Ajuste automatiquement Rp/Rs et T0 aux données observées"""
        if self.processed_data is None or self.model_flux is None or self.planet_obj is None:
            messagebox.showwarning("Erreur", "Générez d'abord un modèle avec 'Générer Modèle'.")
            return
        
        if self.filter_name.get() == "":
            messagebox.showwarning("Erreur", "Spécifiez un filtre d'abord.")
            return
        
        flt = self.filter_name.get()
        if is_gaia_pylightcurve_filter(flt) and not self._ldc_preparation_available(flt):
            messagebox.showerror(
                "Filtre Gaia",
                "Ajustement : activez ExoTETHyS comme source LDC ou installez astroquery (VizieR).",
            )
            return

        try:
            self.config(cursor="watch")
            
            # Données observées
            t_obs = self.processed_data['JD-UTC'].values
            f_obs = self.processed_data['FLUX_FINAL'].values
            # Erreurs si disponibles
            if 'rel_flux_err_T1' in self.processed_data.columns:
                f_err = self.processed_data['rel_flux_err_T1'].values
                f_err = np.where(f_err > 0, f_err, np.std(f_obs) * 0.001)  # Éviter les erreurs nulles
            else:
                # Estimation des erreurs à partir de la dispersion OOT
                oot_mask = (t_obs < self.transit_start.get()) | (t_obs > self.transit_end.get())
                if np.sum(oot_mask) > 1:
                    f_err = np.ones_like(f_obs) * np.std(f_obs[oot_mask])
                else:
                    f_err = np.ones_like(f_obs) * np.std(f_obs) * 0.001
            period = self.period_var.get()
            
            if period <= 0:
                messagebox.showerror("Erreur", "Période invalide.")
                return
            
            # Valeurs initiales
            rprs_initial = float(self.planet_obj.rp_over_rs) if hasattr(self.planet_obj, 'rp_over_rs') else 0.1
            t0_initial = self.t0_var.get()
            
            # Créer une copie de l'objet Planet pour l'ajustement
            pl = self.planet_obj
            
            # Vérifier si on utilise la loi power-2
            ld_law = self.limb_darkening_law.get()
            use_power2 = (ld_law == "power-2" and POWER2_AVAILABLE)
            
            # Récupérer les paramètres orbitaux pour power-2 si nécessaire
            if use_power2:
                a_rs, inclination = self._get_power2_orbital_params()
                c = self.power2_c.get()
                alpha = self.power2_alpha.get()
            
            def chi2_function(params):
                """Fonction objectif : chi2 entre modèle et données"""
                if use_power2 and self.fit_limb_darkening.get():
                    # Ajuster rp_rs, t0, c, et alpha
                    rprs, t0, c_fit, alpha_fit = params
                    
                    # Limites physiques
                    if rprs <= 0 or rprs >= 1 or c_fit < 0 or c_fit > 1 or alpha_fit < 0.1 or alpha_fit > 10:
                        return 1e10
                    
                    try:
                        # Générer le modèle avec power-2
                        model = transit_lightcurve_power2(
                            time=t_obs,
                            period=period,
                            t0=t0,
                            rp_rs=rprs,
                            a_rs=a_rs,
                            inclination=inclination,
                            c=c_fit,
                            alpha=alpha_fit,
                            n_annuli=5000  # Réduire pour la vitesse lors de l'ajustement
                        )
                        
                        # Calculer le chi2
                        residuals = f_obs - model
                        chi2 = np.sum((residuals / f_err)**2)
                        return chi2
                    except:
                        return 1e10
                else:
                    # Ajuster seulement rp_rs et t0
                    rprs, t0 = params
                    
                    # Limites physiques
                    if rprs <= 0 or rprs >= 1:
                        return 1e10
                    
                    try:
                        if use_power2:
                            # Utiliser power-2 avec coefficients fixes
                            model = transit_lightcurve_power2(
                                time=t_obs,
                                period=period,
                                t0=t0,
                                rp_rs=rprs,
                                a_rs=a_rs,
                                inclination=inclination,
                                c=c,
                                alpha=alpha,
                                n_annuli=5000
                            )
                        else:
                            # Utiliser pylightcurve
                            # Mettre à jour les paramètres
                            pl.rp_over_rs = rprs
                            for attr in ['mid_time', 'transit_mid_time', 't0', 'tmid', 'transit_time']:
                                if hasattr(pl, attr):
                                    setattr(pl, attr, t0)
                            
                            # Générer le modèle
                            self._prepare_ldc_for_transit(pl, flt)
                            model = pl.transit_integrated(
                                time=t_obs,
                                time_format='JD_UTC',
                                exp_time=120,
                                time_stamp='mid',
                                max_sub_exp_time=1,
                                filter_name=flt
                            )
                        
                        # Calculer le chi2 pondéré par les erreurs
                        residuals = f_obs - model
                        chi2 = np.sum((residuals / f_err)**2)
                        
                        return chi2
                    except:
                        return 1e10
            
            # Ajustement avec scipy.optimize.minimize
            if use_power2 and self.fit_limb_darkening.get():
                # Ajuster rp_rs, t0, c, et alpha
                c_initial = self.power2_c.get()
                alpha_initial = self.power2_alpha.get()
                initial_params = [rprs_initial, t0_initial, c_initial, alpha_initial]
                
                bounds = [
                    (max(0.01, rprs_initial * 0.5), min(0.5, rprs_initial * 2.0)),
                    (t0_initial - 0.1, t0_initial + 0.1),
                    (0.0, 1.0),  # c
                    (0.1, 10.0)  # alpha
                ]
            else:
                # Ajuster seulement rp_rs et t0
                initial_params = [rprs_initial, t0_initial]
                
                # Bornes : Rp/Rs entre 0.01 et 0.5, T0 autour de ±0.1 jour
                bounds = [(max(0.01, rprs_initial * 0.5), min(0.5, rprs_initial * 2.0)),
                         (t0_initial - 0.1, t0_initial + 0.1)]
            
            result = minimize(chi2_function, initial_params, method='L-BFGS-B', bounds=bounds, 
                            options={'maxiter': 100, 'ftol': 1e-6})
            
            if result.success:
                if use_power2 and self.fit_limb_darkening.get():
                    rprs_fit, t0_fit, c_fit, alpha_fit = result.x
                    # Mettre à jour les coefficients power-2
                    self.power2_c.set(c_fit)
                    self.power2_alpha.set(alpha_fit)
                    logger.info(f"Coefficients power-2 ajustés: c={c_fit:.4f}, α={alpha_fit:.4f}")
                    # Mettre à jour l'affichage des coefficients
                    self.update_ld_coefficients_display()
                else:
                    rprs_fit, t0_fit = result.x
                
                # Stocker les paramètres ajustés pour recalculer le modèle après changements de données
                self.t0_fitted = float(t0_fit)
                self.rprs_fitted = float(rprs_fit)
                logger.info(f"[FITTING] Paramètres ajustés stockés: Rp/Rs={self.rprs_fitted:.5f}, T0={self.t0_fitted:.8f}")
                logger.info(f"[FITTING] Valeurs initiales: Rp/Rs={rprs_initial:.5f}, T0={t0_initial:.8f}")
                
                # Mettre à jour l'objet Planet
                pl.rp_over_rs = rprs_fit
                for attr in ['mid_time', 'transit_mid_time', 't0', 'tmid', 'transit_time']:
                    if hasattr(pl, attr):
                        setattr(pl, attr, t0_fit)
                
                # Mettre à jour les variables
                self.t0_var.set(t0_fit)
                
                # Ne pas recalculer le modèle ajusté ici - update_processing() s'en chargera
                # avec les bons temps et les bons paramètres stockés dans self.rprs_fitted et self.t0_fitted
                logger.info(f"[FITTING] Avant update_processing: self.rprs_fitted={self.rprs_fitted:.5f}, self.t0_fitted={self.t0_fitted:.8f}")
                
                # Recalculer les métriques (ce qui recalculera aussi les modèles)
                self.update_processing()
                
                # Message de succès
                msg = (f"Rp/Rs ajusté: {rprs_fit:.5f} (initial: {rprs_initial:.5f})\n"
                       f"T0 ajusté: {t0_fit:.8f} (initial: {t0_initial:.8f})\n"
                       f"Chi2 final: {result.fun:.2f}")
                
                if use_power2 and self.fit_limb_darkening.get():
                    msg += f"\n\nCoefficients power-2 ajustés:\n"
                    msg += f"c = {c_fit:.4f}\n"
                    msg += f"α = {alpha_fit:.4f}"
                
                messagebox.showinfo("Ajustement réussi", msg)
            else:
                messagebox.showerror("Erreur", f"L'ajustement a échoué: {result.message}")
        
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de l'ajustement:\n{e}")
            logger.exception("Erreur ajustement")
        finally:
            self.config(cursor="")

    def valider_modele(self):
        """Génère final-report.txt (Stats + Data) et report.png"""
        if self.processed_data is None or self.model_flux is None:
            messagebox.showwarning("Erreur", "Données ou modèle manquant.")
            return
        
        # --- 1. PREPARATION DES DONNEES ---
        s = self.stats_cache
        df = self.processed_data.copy()
        
        # Préparation du DataFrame pour l'export (PARTIE 2)
        export_df = pd.DataFrame()
        # Formater JD-UTC avec 6 décimales en format standard (pas scientifique)
        export_df['JD-UTC'] = df['JD-UTC'].apply(lambda x: f"{x:.6f}")
        export_df['rel_flux_fn'] = df['FLUX_RAW']
        
        err_col = df['rel_flux_err_T1'] if 'rel_flux_err_T1' in df.columns else np.zeros(len(df))
        export_df['rel_flux_fn_err'] = err_col
        # Résiduel et Detrend
        export_df['rel_flux_fn_residuals'] = df['FLUX_FINAL'] - self.model_flux
        export_df['rel_flux_fn_detrend'] = df['FLUX_FINAL']
        if 'rel_flux_fn_detrend_err' in df.columns:
            export_df['rel_flux_fn_detrend_err'] = df['rel_flux_fn_detrend_err']
        else:
            export_df['rel_flux_fn_detrend_err'] = err_col

        # Interprétation stats
        shapiro_p = s.get('shapiro_p', 0)
        chi2 = s.get('chi2', 1)
        acf_1 = s.get('acf_1', 0)

        qual_shapiro = "Normal" if shapiro_p > 0.05 else "Rejet H0 (Non-Gaussian)"
        qual_chi2 = "Bon" if 0.8 <= chi2 <= 1.2 else ("Sous-estimé" if chi2 < 0.8 else "Bruit/Mauvais fit")
        qual_acf = "Bon (Bruit Blanc)" if abs(acf_1) < 0.2 else "Corrélé (Bruit Rouge)"

        # --- 2. GENERATION DU TEXTE (PARTIE 1 + PARTIE 2) ---
        _rep_name = self.planet_name.get().strip() or (
            "cible_manuelle" if getattr(self, "manual_model_active", False) else "(sans_nom)"
        )
        report_txt = (
            f"=== RAPPORT ANALYSE : {_rep_name} ===\n\n"
            f"--- 1. CRITÈRES STATISTIQUES & RÉSULTATS ---\n"
            f"Tmid (JD)         : {self.transit_mid_calc.get().replace('Mid: ', '')}\n"
            f"Rp/Rs             : {s.get('rprs', 0):.5f}\n"
            f"O-C (min)         : {s.get('oc', 0):.2f}\n"
            f"Transit Duration  : {s.get('duration', 0):.2f} min\n\n"
            f"Residuals Mean    : {s.get('res_mean', 0):.2e}\n"
            f"Residuals Std     : {s.get('res_std', 0):.2e}\n\n"
            f"Shapiro-Wilk (p)  : {shapiro_p:.4f}  [{qual_shapiro}]\n"
            f"   (Borne: > 0.05 pour distribution Normale)\n"
            f"ACF (Lag 1)       : {acf_1:.4f}      [{qual_acf}]\n"
            f"   (Borne: Proche de 0 pour Bruit Blanc)\n"
            f"Chi2 Réduit       : {chi2:.2f}       [{qual_chi2}]\n"
            f"   (Borne: ~1.0 idéalement)\n\n"
            f"--- 2. TABLEAU DE DONNÉES ---\n"
        )
        report_txt += export_df.to_string(index=False, col_space=15)

        # --- 3. EXPORT DES 2 FICHIERS ---
        try:
            out_dir = Path(self.data_path).parent
            
            # Fichier 1 : final-report.txt
            file_txt = out_dir / "final-report.txt"
            with open(file_txt, "w", encoding="utf-8") as f:
                f.write(report_txt)
            
            # Fichier 2 : report.png (Ajout Text Box sur plot existant)
            file_png = out_dir / "report.png"
            
            # Infos O-C / Tmid sur le graphique
            info_str = (f"Tmid: {self.transit_mid_calc.get().replace('Mid: ', '')}\n"
                        f"O-C: {s.get('oc', 0):.1f} min")
            
            txt_box = self.ax_main.text(0.02, 0.95, info_str, transform=self.ax_main.transAxes, 
                                   fontsize=10, verticalalignment='top', 
                                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            self.fig.savefig(file_png)
            txt_box.remove() # Nettoyage
            self.canvas.draw()

            messagebox.showinfo("Succès", f"Export terminé !\n1. {file_txt}\n2. {file_png}")
            
        except Exception as e:
            messagebox.showerror("Erreur Export", str(e))
    
    def fetch_priors_from_sources(self):
        """Remplit a/R* et i depuis les archives (géométrie ; utilisé si « Appliquer… » est coché)."""
        if not PRIORS_AVAILABLE:
            messagebox.showwarning("Non disponible", 
                "Le module de récupération des priors n'est pas disponible.")
            return
        
        planet_name = self.planet_name.get().strip()
        if not planet_name:
            messagebox.showwarning("Attention", "Entrez d'abord le nom de l'exoplanète.")
            return
        
        try:
            self.config(cursor="watch")
            logger.info(f"Récupération des priors pour {planet_name}...")
            
            # Récupérer les priors depuis toutes les sources disponibles
            priors = get_priors_from_all_sources(planet_name)
            
            if priors is None:
                messagebox.showwarning("Non trouvé", 
                    f"Aucun prior trouvé pour '{planet_name}'.\n"
                    f"Essayez la saisie manuelle ou vérifiez le nom de la planète.")
                self.priors_data = None
                self.priors_label.config(
                    text="Aucune valeur dans les archives — saisie manuelle possible",
                    foreground="gray",
                )
                return
            
            # Stocker les priors
            self.priors_data = priors
            self.prior_a_rs.set(priors['a_rs'])
            self.prior_inclination.set(priors['inclination'])
            
            # Erreurs si disponibles
            if 'a_rs_err' in priors:
                self.prior_a_rs_err.set(priors['a_rs_err'])
            else:
                # Estimation d'erreur par défaut (1% de la valeur)
                self.prior_a_rs_err.set(priors['a_rs'] * 0.01)
            
            if 'inclination_err' in priors:
                self.prior_inclination_err.set(priors['inclination_err'])
            else:
                # Estimation d'erreur par défaut (0.1°)
                self.prior_inclination_err.set(0.1)
            
            # Afficher la source
            source = priors.get('source', 'Inconnue')
            self.priors_source.set(source)

            # Stocker les priors récupérés séparément des priors calculés
            self.database_priors = priors
            
            # Mettre à jour l'affichage
            priors_text = format_priors_for_display(priors)
            self.priors_label.config(text=priors_text, foreground="green")
            
            messagebox.showinfo(
                "Valeurs récupérées",
                f"Source : {source}\n\n"
                f"a/R* = {priors['a_rs']:.4f} ± {self.prior_a_rs_err.get():.4f}\n"
                f"i = {priors['inclination']:.4f}° ± {self.prior_inclination_err.get():.4f}°\n\n"
                "Cochez « Appliquer a/R* et i ci-dessous… » pour les utiliser avec la loi power-2 "
                "et l'outil « Comparer les lois » (sinon la géométrie vient du modèle planète).",
            )
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de la récupération des priors:\n{e}")
            logger.exception("Erreur récupération priors")
            self.priors_data = None
            self.priors_label.config(text="Erreur de récupération", foreground="red")
        finally:
            self.config(cursor="")
    
    def show_manual_priors_dialog(self):
        """Saisie manuelle de a/R* et i (géométrie du transit, pas la loi de limbe)."""
        dialog = tk.Toplevel(self)
        dialog.title("Géométrie : a/R* et inclinaison")
        dialog.geometry("400x250")
        dialog.transient(self)
        dialog.grab_set()
        
        # Variables temporaires
        a_rs_var = tk.DoubleVar(value=self.prior_a_rs.get())
        a_rs_err_var = tk.DoubleVar(value=self.prior_a_rs_err.get())
        i_var = tk.DoubleVar(value=self.prior_inclination.get())
        i_err_var = tk.DoubleVar(value=self.prior_inclination_err.get())
        
        # Frame principal
        main_frame = ttk.Frame(dialog, padding=10)
        main_frame.pack(fill="both", expand=True)
        
        # a/R*
        f_a_rs = ttk.Frame(main_frame)
        f_a_rs.pack(fill="x", pady=5)
        ttk.Label(f_a_rs, text="a/R*:").pack(side="left", padx=5)
        ttk.Entry(f_a_rs, textvariable=a_rs_var, width=15).pack(side="left", padx=5)
        ttk.Label(f_a_rs, text="±").pack(side="left", padx=5)
        ttk.Entry(f_a_rs, textvariable=a_rs_err_var, width=10).pack(side="left", padx=5)
        
        # Inclinaison
        f_i = ttk.Frame(main_frame)
        f_i.pack(fill="x", pady=5)
        ttk.Label(f_i, text="Inclinaison (°):").pack(side="left", padx=5)
        ttk.Entry(f_i, textvariable=i_var, width=15).pack(side="left", padx=5)
        ttk.Label(f_i, text="±").pack(side="left", padx=5)
        ttk.Entry(f_i, textvariable=i_err_var, width=10).pack(side="left", padx=5)
        
        # Boutons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill="x", pady=10)
        
        def save_priors():
            # Valider les valeurs
            a_rs = a_rs_var.get()
            a_rs_err = a_rs_err_var.get()
            i = i_var.get()
            i_err = i_err_var.get()
            
            if a_rs <= 0 or i <= 0 or i > 90:
                messagebox.showerror("Erreur", "Valeurs invalides:\n"
                    "- a/R* doit être > 0\n"
                    "- Inclinaison doit être entre 0° et 90°")
                return
            
            # Sauvegarder
            self.prior_a_rs.set(a_rs)
            self.prior_a_rs_err.set(a_rs_err)
            self.prior_inclination.set(i)
            self.prior_inclination_err.set(i_err)
            
            # Créer un dictionnaire de priors
            self.priors_data = {
                'a_rs': a_rs,
                'a_rs_err': a_rs_err,
                'inclination': i,
                'inclination_err': i_err,
                'source': 'Saisie manuelle'
            }
            self.priors_source.set('Saisie manuelle')
            
            # Mettre à jour l'affichage
            priors_text = format_priors_for_display(self.priors_data)
            self.priors_label.config(text=priors_text, foreground="blue")
            
            dialog.destroy()
            messagebox.showinfo(
                "Valeurs enregistrées",
                f"a/R* = {a_rs:.4f} ± {a_rs_err:.4f}\n"
                f"i = {i:.4f}° ± {i_err:.4f}°\n\n"
                "Cochez « Appliquer a/R* et i ci-dessous… » pour les utiliser avec power-2 "
                "et « Comparer les lois ».",
            )
        
        ttk.Button(btn_frame, text="Enregistrer", command=save_priors).pack(side="right", padx=5)
        ttk.Button(btn_frame, text="Annuler", command=dialog.destroy).pack(side="right", padx=5)
    
    def apply_claret_2022_power2_from_vizier(self):
        """
        Remplit c et α (loi power-2) depuis Claret & Southworth (2022), VizieR J/A+A/664/A128 table1,
        pour les bandes Gaia G, G_BP ou G_RP (colonnes gG/hG, gBP/hBP, gRP/hRP).
        """
        if not CLARET_2022_VIZIER_AVAILABLE:
            messagebox.showwarning(
                "Indisponible",
                "astroquery n'est pas installé ou le module claret_ld_vizier est introuvable.",
            )
            return
        try:
            teff = float(self.teff_var.get())
            logg = float(self.logg_var.get())
            met = float(self.met_var.get())
        except (tk.TclError, ValueError, TypeError):
            messagebox.showerror("Erreur", "Teff, log g et [Fe/H] doivent être numériques.")
            return
        if teff <= 0 or logg <= 0:
            messagebox.showwarning(
                "Attention",
                "Renseignez des valeurs stellaires valides (Teff, log g) dans la section modélisation.",
            )
            return
        raw = self.claret_gaia_passband.get().strip()
        band_map = {"G": "G", "G_BP": "BP", "G_RP": "RP"}
        band_key = band_map.get(raw, "G")
        try:
            self.config(cursor="watch")
            self.update_idletasks()
            info = nearest_power2_gaia(teff, logg, met, 2.0, band=band_key)
            self.power2_c.set(info["c"])
            self.power2_alpha.set(info["alpha"])
            self.update_ld_coefficients_display()
            messagebox.showinfo(
                "Claret & Southworth (2022), VizieR",
                format_lookup_summary(info)
                + "\n\nLes coefficients du catalogue (g, h) correspondent à c et α dans NPOAP.\n"
                "Microturbulence fixée à 2 km/s (valeur courante ; grille aussi 0,1,4,8).",
            )
        except Exception as e:
            logger.exception("Claret 2022 VizieR")
            messagebox.showerror(
                "Erreur",
                f"Impossible de récupérer les coefficients (réseau / VizieR) :\n{e}",
            )
        finally:
            self.config(cursor="")

    def apply_exotethys_power2_from_filter(self):
        """Remplit c et α (power-2) via ExoTETHyS pour le filtre Modélisation (passband + Phoenix/ATLAS)."""
        if not EXOTETHYS_AVAILABLE or run_exotethys_ldc_power2 is None:
            messagebox.showwarning("ExoTETHyS", "Installez exotethys, h5py et click (pip).")
            return
        flt = self.filter_name.get().strip()
        if not flt:
            messagebox.showwarning("Filtre", "Choisissez un filtre dans la section Modélisation.")
            return
        try:
            teff = float(self.teff_var.get())
            logg = float(self.logg_var.get())
            met = float(self.met_var.get())
        except (tk.TclError, ValueError, TypeError):
            messagebox.showerror("Erreur", "Teff, log g et [Fe/H] doivent être numériques.")
            return
        if teff <= 0 or logg <= 0:
            messagebox.showwarning("Attention", "Renseignez Teff et log g valides.")
            return
        try:
            self.config(cursor="watch")
            self.update_idletasks()
            model = self.exotethys_model_var.get()
            ca = run_exotethys_ldc_power2(flt, teff, logg, met, stellar_model=model)
            ca = np.asarray(ca, dtype=float).ravel()
            if len(ca) < 2:
                raise ValueError("ExoTETHyS n'a pas renvoyé deux coefficients power-2.")
            self.power2_c.set(float(ca[0]))
            self.power2_alpha.set(float(ca[1]))
            self.update_ld_coefficients_display()
            messagebox.showinfo(
                "ExoTETHyS (power-2)",
                f"Filtre : {flt}\nModèle : {model}\n"
                f"c = {ca[0]:.6f}\nα = {ca[1]:.6f}\n\n"
                "Premier lancement : téléchargement possible des grilles stellaires (long).",
            )
        except Exception as e:
            logger.exception("ExoTETHyS power-2 UI")
            messagebox.showerror("Erreur", str(e))
        finally:
            self.config(cursor="")

    def update_ld_coefficients_display(self):
        """Met à jour l'affichage des coefficients de limb-darkening."""
        ld_law = self.limb_darkening_law.get()
        if ld_law == "power-2":
            c = self.power2_c.get()
            alpha = self.power2_alpha.get()
            self.ld_coefficients_label.config(
                text=f"Coefficients power-2: c={c:.4f}, α={alpha:.4f}",
                foreground="blue"
            )
        elif ld_law == "pylightcurve":
            if self.ldc_source_var.get() == "exotethys":
                self.ld_coefficients_label.config(
                    text="LDC transit : ExoTETHyS (loi claret4, pylightcurve)",
                    foreground="blue",
                )
            else:
                self.ld_coefficients_label.config(
                    text="LDC transit : grille intégrée pylightcurve",
                    foreground="gray",
                )
        elif ld_law in ("quadratique", "square-root"):
            self.ld_coefficients_label.config(
                text="Courbe analytique : utiliser « Comparer les lois » (pas le bouton Générer modèle)",
                foreground="gray",
            )
        else:
            self.ld_coefficients_label.config(
                text=f"Coefficients {ld_law}: non branché sur Générer modèle",
                foreground="gray",
            )
    
    def _heuristic_ld_u1_u2_for_compare(self):
        """
        u₁, u₂ approximatifs pour les lois quadratique et square-root dans « Comparer les lois ».
        Dépend légèrement de Teff (pas des tables Claret complètes).
        """
        try:
            teff = float(self.teff_var.get())
        except (tk.TclError, TypeError, ValueError):
            teff = 5778.0
        teff = float(np.clip(teff, 3500.0, 8000.0))
        u1 = 0.70 - (teff - 3500.0) * (0.28 / 4500.0)
        u2 = 0.45 - (teff - 3500.0) * (0.20 / 4500.0)
        return float(np.clip(u1, 0.25, 0.65)), float(np.clip(u2, 0.12, 0.42))

    def compare_limb_darkening_laws(self):
        """Affiche une fenêtre de comparaison des différentes lois de limb-darkening."""
        if self.processed_data is None or self.model_flux is None or self.planet_obj is None:
            messagebox.showwarning("Erreur", "Générez d'abord un modèle avec 'Générer Modèle'.")
            return
        
        try:
            # Créer une nouvelle fenêtre
            compare_window = tk.Toplevel(self)
            compare_window.title("Comparaison des Lois de Limb-Darkening")
            compare_window.geometry("1200x800")
            
            # Frame pour les graphiques
            fig_compare = plt.Figure(figsize=(12, 8), dpi=100)
            canvas_compare = FigureCanvasTkAgg(fig_compare, master=compare_window)
            canvas_compare.get_tk_widget().pack(fill="both", expand=True)
            
            # Créer une grille 2x2 pour les graphiques
            gs = fig_compare.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
            ax1 = fig_compare.add_subplot(gs[0, 0])  # Courbes de lumière
            ax2 = fig_compare.add_subplot(gs[0, 1])  # Résidus
            ax3 = fig_compare.add_subplot(gs[1, 0])  # Histogramme résidus
            ax4 = fig_compare.add_subplot(gs[1, 1])  # Autocorrélation
            
            # Données
            t_data = self.processed_data['JD-UTC'].values
            f_obs = self.processed_data['FLUX_FINAL'].values
            flt = self.filter_name.get()
            period = self.period_var.get()
            t0 = self.t0_var.get()
            rprs = float(self.planet_obj.rp_over_rs) if hasattr(self.planet_obj, 'rp_over_rs') else 0.1
            
            colors = {
                "pylightcurve": "blue",
                "power-2": "red",
                "quadratique": "green",
                "square-root": "orange",
            }
            models = {}
            n_ann = 5000
            u1_h, u2_h = self._heuristic_ld_u1_u2_for_compare()
            a_rs, inclination = None, None
            if POWER2_AVAILABLE:
                a_rs, inclination = self._get_power2_orbital_params()

            law_specs = [
                ("pylightcurve", "pylightcurve"),
                ("power-2", "power-2"),
                ("quadratic", "quadratique"),
                ("square-root", "square-root"),
            ]

            for law, legend_label in law_specs:
                try:
                    if law == "pylightcurve":
                        pl = self.planet_obj
                        pl.rp_over_rs = rprs
                        for attr in ['mid_time', 'transit_mid_time', 't0', 'tmid', 'transit_time']:
                            if hasattr(pl, attr):
                                setattr(pl, attr, t0)
                        if not self._ldc_preparation_available(flt):
                            logger.warning("Comparaison pylightcurve : LDC indisponible pour ce filtre — ignorée.")
                            continue
                        self._prepare_ldc_for_transit(pl, flt)
                        model = pl.transit_integrated(
                            time=t_data, time_format='JD_UTC', exp_time=120,
                            time_stamp='mid', max_sub_exp_time=1, filter_name=flt
                        )
                        models[legend_label] = model
                    elif law == "power-2":
                        if not POWER2_AVAILABLE:
                            continue
                        c = self.power2_c.get()
                        alpha = self.power2_alpha.get()
                        model = transit_lightcurve_power2(
                            time=t_data, period=period, t0=t0,
                            rp_rs=rprs, a_rs=a_rs, inclination=inclination,
                            c=c, alpha=alpha, n_annuli=n_ann
                        )
                        models[legend_label] = model
                    elif law == "quadratic":
                        if not POWER2_AVAILABLE:
                            continue
                        model = transit_lightcurve_quadratic(
                            time=t_data, period=period, t0=t0,
                            rp_rs=rprs, a_rs=a_rs, inclination=inclination,
                            u1=u1_h, u2=u2_h, n_annuli=n_ann
                        )
                        models[legend_label] = model
                    elif law == "square-root":
                        if not POWER2_AVAILABLE:
                            continue
                        model = transit_lightcurve_square_root(
                            time=t_data, period=period, t0=t0,
                            rp_rs=rprs, a_rs=a_rs, inclination=inclination,
                            u1=u1_h, u2=u2_h, n_annuli=n_ann
                        )
                        models[legend_label] = model
                except Exception as e:
                    logger.warning(f"Erreur calcul modèle {law} ({legend_label}): {e}")

            if POWER2_AVAILABLE:
                fig_compare.text(
                    0.5,
                    0.99,
                    f"Quadratique & square-root : u₁={u1_h:.3f}, u₂={u2_h:.3f} "
                    f"(heuristique depuis Teff ; pas des tables Claret)",
                    ha="center",
                    va="top",
                    fontsize=8,
                    transform=fig_compare.transFigure,
                )
            
            # Graphique 1: Courbes de lumière comparées
            ax1.plot(t_data, f_obs, 'k.', alpha=0.3, markersize=2, label='Données')
            for law, model in models.items():
                idx = np.argsort(t_data)
                ax1.plot(t_data[idx], model[idx], '-', lw=2, 
                        color=colors.get(law, "black"), label=law, alpha=0.7)
            ax1.set_xlabel("JD-UTC")
            ax1.set_ylabel("Flux normalisé")
            ax1.set_title("Comparaison des Courbes de Lumière")
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Graphique 2: Résidus
            for law, model in models.items():
                residuals = f_obs - model
                ax2.plot(t_data, residuals, '.', alpha=0.5, markersize=2,
                        color=colors.get(law, "black"), label=f"Résidus {law}")
            ax2.axhline(0, color='red', linestyle='--', alpha=0.5)
            ax2.set_xlabel("JD-UTC")
            ax2.set_ylabel("Résidus")
            ax2.set_title("Résidus par Loi")
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            # Graphique 3: Histogramme des résidus
            for law, model in models.items():
                residuals = f_obs - model
                ax3.hist(residuals, bins=50, alpha=0.5, label=law,
                        color=colors.get(law, "black"), density=True)
            ax3.set_xlabel("Résidus")
            ax3.set_ylabel("Densité")
            ax3.set_title("Distribution des Résidus")
            ax3.legend()
            ax3.grid(True, alpha=0.3, axis='y')
            
            # Graphique 4: Autocorrélation
            for law, model in models.items():
                residuals = f_obs - model
                try:
                    acf_vals = acf(residuals, nlags=20)
                    lags = np.arange(len(acf_vals))
                    ax4.plot(lags, acf_vals, 'o-', lw=2, markersize=4,
                            color=colors.get(law, "black"), label=law, alpha=0.7)
                except:
                    pass
            ax4.axhline(0, color='black', linestyle='-', alpha=0.3)
            ax4.axhline(0.2, color='red', linestyle='--', alpha=0.5, label='Seuil ±0.2')
            ax4.axhline(-0.2, color='red', linestyle='--', alpha=0.5)
            ax4.set_xlabel("Lag")
            ax4.set_ylabel("ACF")
            ax4.set_title("Autocorrélation des Résidus")
            ax4.legend()
            ax4.grid(True, alpha=0.3)
            
            if POWER2_AVAILABLE:
                fig_compare.tight_layout(rect=[0, 0, 1, 0.96])
            else:
                fig_compare.tight_layout()
            canvas_compare.draw()
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de la comparaison:\n{e}")
            logger.exception("Erreur comparaison lois")
    
    def show_diagnostic_plots(self):
        """Affiche une fenêtre avec des graphiques de diagnostic détaillés."""
        if self.processed_data is None or self.model_flux is None:
            messagebox.showwarning("Erreur", "Générez d'abord un modèle avec 'Générer Modèle'.")
            return
        
        try:
            # Créer une nouvelle fenêtre
            diag_window = tk.Toplevel(self)
            diag_window.title("Graphiques de Diagnostic")
            diag_window.geometry("1200x900")
            
            # Frame pour les graphiques
            fig_diag = plt.Figure(figsize=(12, 9), dpi=100)
            canvas_diag = FigureCanvasTkAgg(fig_diag, master=diag_window)
            canvas_diag.get_tk_widget().pack(fill="both", expand=True)
            
            # Créer une grille 2x3 pour les graphiques
            gs = fig_diag.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
            ax1 = fig_diag.add_subplot(gs[0, 0])  # Résidus vs temps
            ax2 = fig_diag.add_subplot(gs[0, 1])  # Histogramme résidus
            ax3 = fig_diag.add_subplot(gs[0, 2])  # Q-Q plot
            ax4 = fig_diag.add_subplot(gs[1, 0])  # Autocorrélation
            ax5 = fig_diag.add_subplot(gs[1, 1])  # Résidus vs modèle
            ax6 = fig_diag.add_subplot(gs[1, 2])  # Résidus vs airmass (si disponible)
            
            # Données
            t_data = self.processed_data['JD-UTC'].values
            f_obs = self.processed_data['FLUX_FINAL'].values
            model_for_residuals = self.fitted_model_flux if (self.fitted_model_flux is not None and len(self.fitted_model_flux) == len(f_obs)) else self.model_flux
            
            if model_for_residuals is None or len(model_for_residuals) != len(f_obs):
                messagebox.showwarning("Erreur", "Modèle non disponible ou incompatible.")
                return
            
            residuals = f_obs - model_for_residuals
            
            # Graphique 1: Résidus vs temps
            ax1.plot(t_data, residuals, 'k.', alpha=0.5, markersize=3)
            ax1.axhline(0, color='red', linestyle='--', alpha=0.5)
            ax1.set_xlabel("JD-UTC")
            ax1.set_ylabel("Résidus")
            ax1.set_title("Résidus vs Temps")
            ax1.grid(True, alpha=0.3)
            
            # Graphique 2: Histogramme des résidus
            ax2.hist(residuals, bins=50, alpha=0.7, edgecolor='black', density=True)
            # Ajuster une gaussienne
            from scipy.stats import norm
            mu, sigma = np.mean(residuals), np.std(residuals)
            x_gauss = np.linspace(residuals.min(), residuals.max(), 100)
            ax2.plot(x_gauss, norm.pdf(x_gauss, mu, sigma), 'r-', lw=2, 
                    label=f'Gaussienne (μ={mu:.2e}, σ={sigma:.2e})')
            ax2.set_xlabel("Résidus")
            ax2.set_ylabel("Densité")
            ax2.set_title("Distribution des Résidus")
            ax2.legend()
            ax2.grid(True, alpha=0.3, axis='y')
            
            # Graphique 3: Q-Q plot
            from scipy.stats import probplot
            probplot(residuals, dist="norm", plot=ax3)
            ax3.set_title("Q-Q Plot (Normalité)")
            ax3.grid(True, alpha=0.3)
            
            # Graphique 4: Autocorrélation
            try:
                acf_vals = acf(residuals, nlags=30)
                lags = np.arange(len(acf_vals))
                ax4.plot(lags, acf_vals, 'o-', lw=2, markersize=4, color='blue')
                ax4.axhline(0, color='black', linestyle='-', alpha=0.3)
                ax4.axhline(0.2, color='red', linestyle='--', alpha=0.5, label='Seuil ±0.2')
                ax4.axhline(-0.2, color='red', linestyle='--', alpha=0.5)
                ax4.fill_between(lags, -0.2, 0.2, alpha=0.1, color='green', label='Zone acceptable')
                ax4.set_xlabel("Lag")
                ax4.set_ylabel("ACF")
                ax4.set_title("Autocorrélation des Résidus")
                ax4.legend()
                ax4.grid(True, alpha=0.3)
            except:
                ax4.text(0.5, 0.5, "Erreur calcul ACF", ha='center', va='center', transform=ax4.transAxes)
            
            # Graphique 5: Résidus vs modèle
            ax5.scatter(model_for_residuals, residuals, alpha=0.5, s=10)
            ax5.axhline(0, color='red', linestyle='--', alpha=0.5)
            ax5.set_xlabel("Flux Modèle")
            ax5.set_ylabel("Résidus")
            ax5.set_title("Résidus vs Modèle")
            ax5.grid(True, alpha=0.3)
            
            # Graphique 6: Résidus vs airmass (si disponible)
            if 'AIRMASS' in self.processed_data.columns:
                airmass = self.processed_data['AIRMASS'].values
                ax6.scatter(airmass, residuals, alpha=0.5, s=10)
                ax6.axhline(0, color='red', linestyle='--', alpha=0.5)
                ax6.set_xlabel("AIRMASS")
                ax6.set_ylabel("Résidus")
                ax6.set_title("Résidus vs AIRMASS")
                ax6.grid(True, alpha=0.3)
            else:
                ax6.text(0.5, 0.5, "AIRMASS non disponible", ha='center', va='center', transform=ax6.transAxes)
                ax6.set_title("Résidus vs AIRMASS")
            
            fig_diag.tight_layout()
            canvas_diag.draw()
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de l'affichage des diagnostics:\n{e}")
            logger.exception("Erreur graphiques diagnostic")