# gui/cluster_analysis_tab.py
"""
Onglet pour l'analyse d'amas d'etoiles : age et distance a partir
d'images ou de photometrie dans les filtres Gaia G, G_BP, G_RP.
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

import numpy as np

try:
    from astropy.io import fits
    from astropy.table import Table
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    import astropy.units as u
    from astropy.stats import sigma_clipped_stats
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False
    WCS = None

try:
    from astroquery.gaia import Gaia
    from astroquery.vizier import Vizier
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False

try:
    from photutils.detection import DAOStarFinder
    from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
    PHOTUTILS_AVAILABLE = True
except ImportError:
    PHOTUTILS_AVAILABLE = False

try:
    from core.parsec_isochrones import PARSEC_AVAILABLE, fetch_parsec_turnoff_grid
except ImportError:
    PARSEC_AVAILABLE = False
    fetch_parsec_turnoff_grid = None

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Amas de reference (nom, rayon approximatif en deg, distance connue en pc)
# Nord : Hyades, Pleiades, Praesepe, M67. Sud : IC 2602, IC 2391, NGC 2516, NGC 3532, Blanco 1.
REFERENCE_CLUSTERS = [
    ("Hyades", 5.0, 46.3),
    ("Pleiades (M45)", 2.0, 136.2),
    ("Praesepe (M44)", 1.5, 187.0),
    ("M67", 0.5, 908.0),
    ("IC 2602 (hem. sud)", 1.5, 147.0),
    ("IC 2391 (hem. sud)", 1.2, 153.0),
    ("NGC 2516 (hem. sud)", 1.0, 409.0),
    ("NGC 3532 (hem. sud)", 1.0, 485.0),
    ("Blanco 1 (hem. sud)", 1.0, 265.0),
]


class ClusterAnalysisTab(ttk.Frame):
    """
    Onglet Analyse d'amas : diagramme couleur-magnitude (G, G_BP, G_RP),
    estimation de la distance par comparaison a un amas de reference,
    et estimation de l'age.
    """

    def __init__(self, parent_notebook):
        super().__init__(parent_notebook, padding=10)
        self.parent_notebook = parent_notebook
        self.data_table = None  # Table avec G, G_BP, G_RP (magnitudes)
        self.ref_table = None  # Amas de reference
        self._data_from_images = False
        self._wcs = None
        self._distance_pc = None
        self._distance_modulus = None
        self._age_gyr = None
        self._age_model = None  # "PARSEC" ou "grille empirique"
        from gui.manual_help import add_manual_help_header

        add_manual_help_header(self, "analyse-amas-open-clusters")
        self._build_ui()

    def _build_ui(self):
        main = ttk.Frame(self)
        main.pack(fill=tk.BOTH, expand=True)

        # ---- Gauche : controles ----
        left = ttk.LabelFrame(main, text="Donnees et parametres", padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))

        ttk.Label(left, text="Source des donnees").pack(anchor=tk.W)
        self.mode_var = tk.StringVar(value="catalog")
        ttk.Radiobutton(left, text="Catalogue Gaia (nom ou coordonnees)", variable=self.mode_var, value="catalog", command=self._on_mode_change).pack(anchor=tk.W)
        ttk.Radiobutton(left, text="3 images FITS (G, G_BP, G_RP)", variable=self.mode_var, value="images", command=self._on_mode_change).pack(anchor=tk.W)

        # Mode catalogue
        self.frm_catalog = ttk.Frame(left)
        self.frm_catalog.pack(fill=tk.X, pady=4)
        ttk.Label(self.frm_catalog, text="Nom amas ou RA,Dec (deg):").pack(anchor=tk.W)
        self.entry_query = ttk.Entry(self.frm_catalog, width=35)
        self.entry_query.pack(fill=tk.X, pady=2)
        self.entry_query.insert(0, "Hyades")
        ttk.Label(self.frm_catalog, text="Rayon (deg):").pack(anchor=tk.W)
        self.entry_radius = ttk.Entry(self.frm_catalog, width=10)
        self.entry_radius.pack(fill=tk.X, pady=2)
        self.entry_radius.insert(0, "3.0")
        ttk.Button(self.frm_catalog, text="Charger depuis Gaia", command=self._load_from_gaia).pack(pady=4)

        # Mode images
        self.frm_images = ttk.Frame(left)
        self.frm_images.pack(fill=tk.X, pady=4)
        ttk.Label(self.frm_images, text="FITS G:").pack(anchor=tk.W)
        self.entry_g = ttk.Entry(self.frm_images, width=32)
        self.entry_g.pack(fill=tk.X, pady=1)
        ttk.Button(self.frm_images, text="Parcourir...", command=lambda: self._browse_fits("G")).pack(anchor=tk.W)
        ttk.Label(self.frm_images, text="FITS G_BP:").pack(anchor=tk.W)
        self.entry_gbp = ttk.Entry(self.frm_images, width=32)
        self.entry_gbp.pack(fill=tk.X, pady=1)
        ttk.Button(self.frm_images, text="Parcourir...", command=lambda: self._browse_fits("G_BP")).pack(anchor=tk.W)
        ttk.Label(self.frm_images, text="FITS G_RP:").pack(anchor=tk.W)
        self.entry_grp = ttk.Entry(self.frm_images, width=32)
        self.entry_grp.pack(fill=tk.X, pady=1)
        ttk.Button(self.frm_images, text="Parcourir...", command=lambda: self._browse_fits("G_RP")).pack(anchor=tk.W)
        ttk.Button(self.frm_images, text="Photometrie et associer", command=self._load_from_images).pack(pady=4)
        ttk.Separator(self.frm_images, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)
        ttk.Label(self.frm_images, text="Calibration (images uniquement):").pack(anchor=tk.W)
        ttk.Button(self.frm_images, text="Calibrer avec etoiles de ref. Gaia", command=self._calibrate_with_gaia_reference).pack(pady=2)
        self.label_calib_status = ttk.Label(self.frm_images, text="", foreground="gray")
        self.label_calib_status.pack(anchor=tk.W)

        self._on_mode_change()

        # Amas de reference
        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)
        ttk.Label(left, text="Amas de reference (distance connue):").pack(anchor=tk.W)
        self.combo_ref = ttk.Combobox(left, values=[r[0] for r in REFERENCE_CLUSTERS], state="readonly", width=22)
        self.combo_ref.pack(fill=tk.X, pady=2)
        self.combo_ref.current(0)
        ttk.Button(left, text="Charger amas de reference", command=self._load_reference_cluster).pack(pady=4)

        ttk.Separator(left, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)
        ttk.Button(left, text="Tracer diagramme CMD", command=self._plot_cmd).pack(pady=2)
        ttk.Button(left, text="Estimer distance (ajustement MS)", command=self._estimate_distance).pack(pady=2)
        ttk.Button(left, text="Estimer age (tour de courbe)", command=self._estimate_age).pack(pady=2)
        ttk.Button(left, text="Exporter tableau CSV", command=self._export_csv).pack(pady=2)

        # ---- Droite : graphique ----
        right = ttk.Frame(main)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.fig = Figure(figsize=(8, 6), dpi=100) if MATPLOTLIB_AVAILABLE else None
        self.ax = self.fig.add_subplot(111) if self.fig else None
        self.ax.set_xlabel("(G_BP - G_RP) [mag]") if self.ax else None
        self.ax.set_ylabel("G [mag]") if self.ax else None
        self.ax.invert_yaxis() if self.ax else None
        self.canvas = FigureCanvasTkAgg(self.fig, master=right) if self.fig else None
        if self.canvas:
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self.toolbar = NavigationToolbar2Tk(self.canvas, right)
            self.toolbar.update()

        # Zone resultats
        self.frm_results = ttk.LabelFrame(main, text="Resultats (distance / age)", padding=6)
        self.frm_results.pack(side=tk.BOTTOM, fill=tk.X, pady=8)
        self.label_results = ttk.Label(self.frm_results, text="Chargez des donnees puis estimez la distance.")
        self.label_results.pack(anchor=tk.W)

    def _on_mode_change(self):
        if self.mode_var.get() == "catalog":
            self.frm_catalog.pack(fill=tk.X, pady=4)
            self.frm_images.pack_forget()
        else:
            self.frm_catalog.pack_forget()
            self.frm_images.pack(fill=tk.X, pady=4)

    def _browse_fits(self, band):
        path = filedialog.askopenfilename(
            title=f"Image FITS {band}",
            filetypes=[("FITS", "*.fits *.fit"), ("Tous", "*.*")]
        )
        if not path:
            return
        if band == "G":
            self.entry_g.delete(0, tk.END)
            self.entry_g.insert(0, path)
        elif band == "G_BP":
            self.entry_gbp.delete(0, tk.END)
            self.entry_gbp.insert(0, path)
        else:
            self.entry_grp.delete(0, tk.END)
            self.entry_grp.insert(0, path)

    def _load_from_gaia(self):
        if not ASTROQUERY_AVAILABLE:
            messagebox.showerror("Erreur", "astroquery n'est pas installe. pip install astroquery")
            return
        query = self.entry_query.get().strip()
        try:
            radius_deg = float(self.entry_radius.get().strip().replace(",", "."))
        except ValueError:
            radius_deg = 3.0
        if not query:
            messagebox.showwarning("Attention", "Saisissez un nom d'amas ou RA,Dec.")
            return
        try:
            job = Gaia.launch_job_async(
                f"""
                SELECT source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
                FROM gaiadr3.gaia_source
                WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {query})) = 1
                AND phot_g_mean_mag IS NOT NULL AND phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL
                AND phot_g_mean_mag < 21
                """
            )
        except Exception as e:
            # Requete par nom via ADQL ou par coordonnees
            try:
                from astroquery.simbad import Simbad
                result = Simbad.query_object(query)
                if result and len(result) > 0:
                    ra = result["RA"][0].replace(" ", ":")
                    dec = result["DEC"][0].replace(" ", ":")
                    from astropy.coordinates import SkyCoord
                    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
                    ra_deg, dec_deg = c.ra.deg, c.dec.deg
                else:
                    ra_deg, dec_deg = float(query.split(",")[0].strip()), float(query.split(",")[1].strip())
            except Exception:
                messagebox.showerror("Erreur", "Impossible de resoudre la position. Utilisez 'RA,Dec' en degres (ex: 66.5, 15.9).")
                return
            try:
                job = Gaia.launch_job_async(
                    f"""
                    SELECT source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
                    FROM gaiadr3.gaia_source
                    WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg})) = 1
                    AND phot_g_mean_mag IS NOT NULL AND phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL
                    AND phot_g_mean_mag < 21
                    """
                )
            except Exception as e2:
                messagebox.showerror("Erreur", f"echec requete Gaia: {e2}")
                return
        try:
            table = job.get_results()
        except Exception as e:
            messagebox.showerror("Erreur", f"echec recuperation resultats Gaia: {e}")
            return
        self.data_table = table
        self._data_from_images = False
        self._wcs = None
        if hasattr(self, "label_calib_status"):
            self.label_calib_status.config(text="")
        self.label_results.config(text=f"Amas charge: {len(table)} etoiles (Gaia).")
        messagebox.showinfo("OK", f"{len(table)} etoiles chargees depuis Gaia.")

    def _load_from_images(self):
        if not PHOTUTILS_AVAILABLE or not ASTROPY_AVAILABLE:
            messagebox.showerror("Erreur", "photutils et astropy requis pour la photometrie sur images.")
            return
        path_g = self.entry_g.get().strip()
        path_gbp = self.entry_gbp.get().strip()
        path_grp = self.entry_grp.get().strip()
        if not path_g or not path_gbp or not path_grp:
            messagebox.showwarning("Attention", "Indiquez les 3 fichiers FITS (G, G_BP, G_RP).")
            return
        try:
            data_g, _, wcs_g = self._load_fits_and_photometry(path_g)
            data_gbp, _, _ = self._load_fits_and_photometry(path_gbp)
            data_grp, _, _ = self._load_fits_and_photometry(path_grp)
        except Exception as e:
            messagebox.showerror("Erreur", f"Chargement/photometrie: {e}")
            return
        if len(data_g) == 0:
            messagebox.showwarning("Attention", "Aucune source detectee dans l'image G.")
            return
        from scipy.spatial import cKDTree
        def match_xy(pos_ref, pos_other, mag_other, max_sep=5.0):
            tree = cKDTree(pos_other)
            dists, idx = tree.query(pos_ref, k=1, distance_upper_bound=max_sep)
            out = np.full(len(pos_ref), np.nan)
            valid = np.isfinite(dists)
            out[valid] = mag_other[idx[valid]]
            return out
        pos_g = np.column_stack((data_g["x"], data_g["y"]))
        pos_gbp = np.column_stack((data_gbp["x"], data_gbp["y"])) if len(data_gbp) > 0 else np.zeros((0, 2))
        pos_grp = np.column_stack((data_grp["x"], data_grp["y"])) if len(data_grp) > 0 else np.zeros((0, 2))
        mag_g = data_g["mag"] if "mag" in data_g.dtype.names else -2.5 * np.log10(data_g["flux"] + 1e-30)
        mag_gbp = data_gbp["mag"] if len(data_gbp) and "mag" in data_gbp.dtype.names else np.full(len(data_gbp), np.nan)
        if len(data_gbp) and "mag" not in data_gbp.dtype.names:
            mag_gbp = -2.5 * np.log10(data_gbp["flux"] + 1e-30)
        mag_grp = data_grp["mag"] if len(data_grp) and "mag" in data_grp.dtype.names else np.full(len(data_grp), np.nan)
        if len(data_grp) and "mag" not in data_grp.dtype.names:
            mag_grp = -2.5 * np.log10(data_grp["flux"] + 1e-30)
        gbp_matched = match_xy(pos_g, pos_gbp, mag_gbp)
        grp_matched = match_xy(pos_g, pos_grp, mag_grp)
        n = len(data_g)
        self.data_table = Table({
            "x_pix": pos_g[:, 0],
            "y_pix": pos_g[:, 1],
            "phot_g_mean_mag": mag_g,
            "phot_bp_mean_mag": gbp_matched,
            "phot_rp_mean_mag": grp_matched,
        })
        self._wcs = wcs_g
        self._data_from_images = True
        n_ok = int(np.sum(np.isfinite(gbp_matched) & np.isfinite(grp_matched)))
        calib_hint = " WCS present: calibration Gaia possible." if self._wcs is not None else " Pas de WCS: ajoutez-en pour calibrer."
        if hasattr(self, "label_calib_status"):
            self.label_calib_status.config(text=calib_hint)
        self.label_results.config(text=f"Images: {n} sources ({n_ok} avec 3 bandes).{calib_hint}")
        messagebox.showinfo("OK", f"{n} sources; {n_ok} avec G, G_BP, G_RP.{calib_hint}")

    def _load_fits_and_photometry(self, path):
        with fits.open(path) as hdul:
            data = hdul[0].data
            header = hdul[0].header
        wcs = None
        if ASTROPY_AVAILABLE and WCS is not None:
            try:
                wcs = WCS(header)
                if not wcs.has_celestial:
                    wcs = None
            except Exception:
                wcs = None
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        finder = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
        sources = finder(data - median)
        if sources is None or len(sources) == 0:
            return np.array([]), data, wcs
        positions = np.transpose((sources["xcentroid"], sources["ycentroid"]))
        apertures = CircularAperture(positions, r=4.0)
        annuli = CircularAnnulus(positions, r_in=8.0, r_out=12.0)
        ap_flux = aperture_photometry(data - median, apertures)
        bkg_flux = aperture_photometry(data - median, annuli)
        n_pix_ap = np.pi * 4.0 ** 2
        n_pix_ann = np.pi * (12.0**2 - 8.0**2)
        bkg_mean = bkg_flux["aperture_sum"] / n_pix_ann
        flux = ap_flux["aperture_sum"] - bkg_mean * n_pix_ap
        mag = -2.5 * np.log10(flux + 1e-30)
        out = np.array(list(zip(positions[:, 0], positions[:, 1], flux, mag)), dtype=[("x", float), ("y", float), ("flux", float), ("mag", float)])
        return out, data, wcs

    def _calibrate_with_gaia_reference(self):
        """Calibre les magnitudes instrumentales en systeme Gaia via association aux etoiles Gaia du champ."""
        if not self._data_from_images or self.data_table is None:
            messagebox.showwarning("Attention", "Calibration possible uniquement pour des donnees chargees depuis des images.")
            return
        if self._wcs is None:
            messagebox.showwarning("Attention", "Les en-tetes FITS doivent contenir un WCS pour associer aux etoiles Gaia.")
            return
        if "x_pix" not in self.data_table.colnames or "y_pix" not in self.data_table.colnames:
            messagebox.showwarning("Attention", "Donnees images sans positions pixel.")
            return
        if not ASTROQUERY_AVAILABLE:
            messagebox.showerror("Erreur", "astroquery requis pour interroger Gaia.")
            return
        try:
            x = self.data_table["x_pix"]
            y = self.data_table["y_pix"]
            coords = self._wcs.pixel_to_world(x, y)
            if hasattr(coords, 'ra'):
                ra_deg = coords.ra.deg
                dec_deg = coords.dec.deg
            else:
                ra_deg = np.array([c.ra.deg for c in coords])
                dec_deg = np.array([c.dec.deg for c in coords])
            center_ra = np.mean(ra_deg)
            center_dec = np.mean(dec_deg)
            radius_deg = max(0.05, 1.5 * np.max([np.ptp(ra_deg), np.ptp(dec_deg)]) * np.cos(np.radians(center_dec)))
            job = Gaia.launch_job_async(
                f"""
                SELECT source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
                FROM gaiadr3.gaia_source
                WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {center_ra}, {center_dec}, {radius_deg})) = 1
                AND phot_g_mean_mag IS NOT NULL AND phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL
                AND phot_g_mean_mag < 21 AND phot_g_mean_mag > 10
                """
            )
            gaia_t = job.get_results()
            gaia_coords = SkyCoord(ra=gaia_t["ra"] * u.deg, dec=gaia_t["dec"] * u.deg)
            our_coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
            idx_gaia, sep2d, _ = our_coords.match_to_catalog_sky(gaia_coords)
            max_sep_arcsec = 2.0
            matched = sep2d.arcsec <= max_sep_arcsec
            n_matched = int(np.sum(matched))
            if n_matched < 5:
                messagebox.showwarning("Attention", f"Seulement {n_matched} associations (< 5). Verifiez le WCS ou le champ.")
                return
            g_inst = self.data_table["phot_g_mean_mag"]
            bp_inst = self.data_table["phot_bp_mean_mag"]
            rp_inst = self.data_table["phot_rp_mean_mag"]
            g_gaia = gaia_t["phot_g_mean_mag"][idx_gaia]
            bp_gaia = gaia_t["phot_bp_mean_mag"][idx_gaia]
            rp_gaia = gaia_t["phot_rp_mean_mag"][idx_gaia]
            valid = matched & np.isfinite(g_inst) & np.isfinite(bp_inst) & np.isfinite(rp_inst)
            valid = valid & np.isfinite(g_gaia) & np.isfinite(bp_gaia) & np.isfinite(rp_gaia)
            zp_g = np.median(np.asarray(g_inst)[valid] - np.asarray(g_gaia)[valid])
            zp_bp = np.median(np.asarray(bp_inst)[valid] - np.asarray(bp_gaia)[valid])
            zp_rp = np.median(np.asarray(rp_inst)[valid] - np.asarray(rp_gaia)[valid])
            self.data_table["phot_g_mean_mag"] = g_inst - zp_g
            self.data_table["phot_bp_mean_mag"] = np.where(np.isfinite(bp_inst), bp_inst - zp_bp, np.nan)
            self.data_table["phot_rp_mean_mag"] = np.where(np.isfinite(rp_inst), rp_inst - zp_rp, np.nan)
            self.label_calib_status.config(text=f"Calibre (ZP_G={zp_g:.2f}, {n_matched} ref.)")
            self.label_results.config(text=f"Calibration Gaia appliquee ({n_matched} etoiles de ref., ZP_G={zp_g:.2f}).")
            messagebox.showinfo("OK", f"Calibration appliquee avec {n_matched} etoiles Gaia.\nZP_G={zp_g:.2f}, ZP_BP={zp_bp:.2f}, ZP_RP={zp_rp:.2f}")
        except Exception as e:
            logger.exception(e)
            messagebox.showerror("Erreur", f"Calibration: {e}")

    def _load_reference_cluster(self):
        if not ASTROQUERY_AVAILABLE:
            messagebox.showerror("Erreur", "astroquery requis.")
            return
        idx = self.combo_ref.current()
        name, radius, dist_pc = REFERENCE_CLUSTERS[idx]
        try:
            Vizier.ROW_LIMIT = 5000
            v = Vizier(columns=["*"], catalog="I/355/gaiadr3")
            if "Hyades" in name:
                result = v.query_object("Hyades")
            elif "Pleiades" in name or "M45" in name:
                result = v.query_object("M45")
            elif "Praesepe" in name or "M44" in name:
                result = v.query_object("M44")
            elif "M67" in name:
                result = v.query_object("M67")
            elif "IC 2602" in name:
                result = v.query_object("IC 2602")
            elif "IC 2391" in name:
                result = v.query_object("IC 2391")
            elif "NGC 2516" in name:
                result = v.query_object("NGC 2516")
            elif "NGC 3532" in name:
                result = v.query_object("NGC 3532")
            elif "Blanco 1" in name:
                result = v.query_object("Blanco 1")
            else:
                result = v.query_object(name)
            if not result or len(result) == 0:
                messagebox.showwarning("Attention", f"Aucune donnee pour {name}.")
                return
            t = result[0]
            if "Gmag" in t.colnames and "BPmag" in t.colnames and "RPmag" in t.colnames:
                self.ref_table = Table({
                    "phot_g_mean_mag": t["Gmag"],
                    "phot_bp_mean_mag": t["BPmag"],
                    "phot_rp_mean_mag": t["RPmag"],
                })
            else:
                self.ref_table = None
                messagebox.showwarning("Attention", f"Colonnes G/BP/RP non trouvees pour {name}.")
                return
            self.label_results.config(text=f"Reference: {name} ({len(self.ref_table)} etoiles), distance {dist_pc} pc.")
        except Exception as e:
            messagebox.showerror("Erreur", f"Chargement amas de reference: {e}")
            logger.exception(e)

    def _plot_cmd(self):
        if self.data_table is None or len(self.data_table) == 0:
            messagebox.showwarning("Attention", "Chargez d'abord les donnees de l'amas.")
            return
        if not MATPLOTLIB_AVAILABLE or self.ax is None:
            return
        self.ax.clear()
        g = self.data_table["phot_g_mean_mag"]
        bp = self.data_table["phot_bp_mean_mag"]
        rp = self.data_table["phot_rp_mean_mag"]
        color = bp - rp
        valid = np.isfinite(g) & np.isfinite(color)
        self.ax.scatter(color[valid], g[valid], s=5, alpha=0.7, label="Amas", color="blue")
        if self.ref_table is not None and len(self.ref_table) > 0:
            gr = self.ref_table["phot_g_mean_mag"]
            br = self.ref_table["phot_bp_mean_mag"]
            rr = self.ref_table["phot_rp_mean_mag"]
            cr = br - rr
            valid_r = np.isfinite(gr) & np.isfinite(cr)
            if self._distance_modulus is not None and self.combo_ref.current() >= 0:
                _, _, d_ref_pc = REFERENCE_CLUSTERS[self.combo_ref.current()]
                d_ref_mod = 5.0 * np.log10(d_ref_pc) - 5.0
                shift = self._distance_modulus - d_ref_mod
                self.ax.scatter(cr[valid_r], np.asarray(gr[valid_r]) + shift, s=3, alpha=0.5, label="Ref. (ajustee)", color="green")
            else:
                self.ax.scatter(cr[valid_r], gr[valid_r], s=3, alpha=0.4, label="Reference", color="gray")
        title = "CMD"
        if self._distance_pc is not None:
            title += f"  |  d = {self._distance_pc:.0f} pc"
        if self._age_gyr is not None:
            title += f"  |  age = {self._age_gyr:.2f} Gyr"
        if self._age_model is not None:
            title += f"  ({self._age_model})"
        self.ax.set_title(title)
        self.ax.set_xlabel("(G_BP - G_RP) [mag]")
        self.ax.set_ylabel("G [mag]")
        self.ax.invert_yaxis()
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw()

    def _estimate_distance(self):
        if self.data_table is None or self.ref_table is None:
            messagebox.showwarning("Attention", "Chargez l'amas et un amas de reference.")
            return
        idx = self.combo_ref.current()
        _, _, d_ref_pc = REFERENCE_CLUSTERS[idx]
        # Ajustement grossier : decalage en magnitude (module de distance)
        g = self.data_table["phot_g_mean_mag"]
        color = self.data_table["phot_bp_mean_mag"] - self.data_table["phot_rp_mean_mag"]
        gr = self.ref_table["phot_g_mean_mag"]
        cr = self.ref_table["phot_bp_mean_mag"] - self.ref_table["phot_rp_mean_mag"]
        valid = np.isfinite(g) & np.isfinite(color) & (color > 0.2) & (color < 2.5) & (g > 10) & (g < 18)
        valid_r = np.isfinite(gr) & np.isfinite(cr) & (cr > 0.2) & (cr < 2.5) & (gr > 10) & (gr < 18)
        if np.sum(valid) < 5 or np.sum(valid_r) < 5:
            messagebox.showwarning("Attention", "Pas assez de points en sequence principale pour l'ajustement.")
            return
        # Bins en couleur, magnitude mediane pour chaque bin
        from scipy.stats import binned_statistic
        bins_c = np.linspace(0.3, 2.2, 15)
        ms_obs, _, _ = binned_statistic(color[valid], g[valid], statistic="median", bins=bins_c)
        ms_ref, be, _ = binned_statistic(cr[valid_r], gr[valid_r], statistic="median", bins=bins_c)
        bc = (be[:-1] + be[1:]) / 2
        ok = np.isfinite(ms_obs) & np.isfinite(ms_ref)
        if np.sum(ok) < 3:
            messagebox.showwarning("Attention", "Ajustement impossible (bins vides).")
            return
        dm = np.median(ms_obs[ok] - ms_ref[ok])
        # Module de distance: m - M = 5*log10(d) - 5  => d = 10^((m-M+5)/5)
        d_ref_mod = 5.0 * np.log10(d_ref_pc) - 5.0
        d_amas_mod = d_ref_mod + dm
        d_amas_pc = 10.0 ** ((d_amas_mod + 5.0) / 5.0)
        self._distance_pc = float(d_amas_pc)
        self._distance_modulus = float(d_amas_mod)
        txt = f"Distance: {d_amas_pc:.1f} pc (module {d_amas_mod:.2f})."
        if self._age_gyr is not None:
            txt += f" Age: {self._age_gyr:.2f} Gyr."
        self.label_results.config(text=txt)
        self._plot_cmd()
        messagebox.showinfo("Distance", f"Distance estimee: {d_amas_pc:.1f} pc\n(module de distance m-M = {d_amas_mod:.2f} mag)\n\nVous pouvez maintenant estimer l'age.")

    def _estimate_age(self):
        """Estime l'age de l'amas par la position du tour de courbe (turn-off) par rapport a des isochrones."""
        if self.data_table is None or len(self.data_table) == 0:
            messagebox.showwarning("Attention", "Chargez d'abord les donnees de l'amas.")
            return
        if self._distance_modulus is None:
            messagebox.showwarning("Attention", "Estimez d'abord la distance (necessaire pour passer en magnitude absolue).")
            return
        g = np.asarray(self.data_table["phot_g_mean_mag"])
        bp = np.asarray(self.data_table["phot_bp_mean_mag"])
        rp = np.asarray(self.data_table["phot_rp_mean_mag"])
        color = bp - rp
        valid = np.isfinite(g) & np.isfinite(color) & (color >= 0.2) & (color <= 2.2) & (g <= 18)
        if np.sum(valid) < 10:
            messagebox.showwarning("Attention", "Pas assez de points valides pour estimer l'age.")
            return
        M_G = g[valid] - self._distance_modulus
        color_v = color[valid]
        from scipy.stats import binned_statistic
        bins_c = np.linspace(0.3, 2.0, 18)
        ms_bright, be, _ = binned_statistic(color_v, M_G, statistic=lambda x: np.percentile(x, 10), bins=bins_c)
        bc = (be[:-1] + be[1:]) / 2
        ok = np.isfinite(ms_bright)
        if np.sum(ok) < 3:
            messagebox.showwarning("Attention", "Impossible de determiner le tour de courbe.")
            return
        idx_min = np.nanargmin(ms_bright)
        color_to = float(bc[idx_min])
        M_G_to = float(ms_bright[idx_min])
        from scipy.interpolate import interp1d
        model_name = "PARSEC"
        if PARSEC_AVAILABLE and fetch_parsec_turnoff_grid is not None:
            grid = fetch_parsec_turnoff_grid(logage_min=8.0, logage_max=10.2, step=0.2, MH=0.0)
            if grid and len(grid) >= 3:
                log_age_yr = np.array([g[0] for g in grid])
                color_iso = np.array([g[1] for g in grid])
                M_G_iso = np.array([g[2] for g in grid])
                ok = np.isfinite(color_iso) & np.isfinite(M_G_iso)
                if np.sum(ok) >= 3:
                    log_age_yr = log_age_yr[ok]
                    color_iso = color_iso[ok]
                    M_G_iso = M_G_iso[ok]
                    d2 = (color_to - color_iso) ** 2 + (M_G_to - M_G_iso) ** 2
                    idx_best = np.argmin(d2)
                    if idx_best == 0 or idx_best == len(log_age_yr) - 1:
                        log_age_est = log_age_yr[idx_best]
                    else:
                        f_log = interp1d(color_iso, log_age_yr, kind="linear", fill_value="extrapolate")
                        log_age_est = float(f_log(color_to))
                    age_gyr = 10.0 ** (log_age_est - 9.0)
                    self._age_gyr = age_gyr
                    self._age_model = "PARSEC"
                    txt = f"Distance: {self._distance_pc:.1f} pc. Age: {age_gyr:.2f} Gyr (PARSEC)."
                    self.label_results.config(text=txt)
                    self._plot_cmd()
                    messagebox.showinfo("Age", f"Age estime: {age_gyr:.2f} Gyr\n(log(age/an) ~ {log_age_est:.1f})\n\nModele: PARSEC (ezpadova), metallicite solaire.")
                    return
        log_age_yr = np.array([8.0, 8.5, 9.0, 9.3, 9.5, 9.7, 10.0, 10.2])
        color_iso = np.array([0.25, 0.45, 0.62, 0.78, 0.88, 0.98, 1.08, 1.15])
        M_G_iso = np.array([1.8, 2.6, 3.2, 3.7, 4.0, 4.25, 4.5, 4.7])
        d2 = (color_to - color_iso) ** 2 + (M_G_to - M_G_iso) ** 2
        idx_best = np.argmin(d2)
        if idx_best == 0 or idx_best == len(log_age_yr) - 1:
            log_age_est = log_age_yr[idx_best]
        else:
            f_log = interp1d(color_iso, log_age_yr, kind="linear", fill_value="extrapolate")
            log_age_est = float(f_log(color_to))
        age_gyr = 10.0 ** (log_age_est - 9.0)
        self._age_gyr = age_gyr
        self._age_model = "grille empirique"
        model_name = "grille empirique (PARSEC: pip install ezpadova)"
        txt = f"Distance: {self._distance_pc:.1f} pc. Age: {age_gyr:.2f} Gyr ({model_name})."
        self.label_results.config(text=txt)
        self._plot_cmd()
        messagebox.showinfo("Age", f"Age estime: {age_gyr:.2f} Gyr\n(log(age/an) ~ {log_age_est:.1f})\n\nModele: {model_name}")
        return

    def _export_csv(self):
        if self.data_table is None or len(self.data_table) == 0:
            messagebox.showwarning("Attention", "Aucune donnee a exporter.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv"), ("Tous", "*.*")])
        if not path:
            return
        try:
            self.data_table.write(path, format="csv", overwrite=True)
            messagebox.showinfo("OK", f"Tableau exporte vers {path}")
        except Exception as e:
            messagebox.showerror("Erreur", str(e))
