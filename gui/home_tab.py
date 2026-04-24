import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import json
from pathlib import Path
import importlib

import config
from core.api_key_dialog import APIKeyDialog, API_KEY_PATH
from gui.manual_help import add_manual_help_header_grid
from config import OBSERVATORY, EQUIPMENT_OBSERVATION


class HomeTab(ttk.Frame):
    def __init__(self, parent, base_dir=None):
        super().__init__(parent, padding=10)

        if base_dir is None:
            # racine du projet: NPOAP/
            home_tab_dir = Path(__file__).resolve().parent  # gui/
            self.base_dir = home_tab_dir.parent             # NPOAP/
        else:
            self.base_dir = Path(base_dir)

        self.config_file = self.base_dir / "config.json"
        self.config_data = {}

        # Variables observatoire
        self.obs_name_var = tk.StringVar()
        self.obs_lat_var = tk.DoubleVar()
        self.obs_lon_var = tk.DoubleVar()
        self.obs_elev_var = tk.DoubleVar()

        # Clé API (affichage masqué)
        self.api_var = tk.StringVar()

        # Configuration matériel (équipement)
        self.equip_entries = {}
        self.equip_computed_vars = {}  # FOV et échelle pixel (affichage calculé)

        # Affiliations
        self.affiliations = []
        self.aff_widgets = []

        self.load_config()
        self.create_widgets()
        self.load_values()

    # ------------------------------------------------------------------
    # Chargement / sauvegarde de la config
    # ------------------------------------------------------------------
    def load_config(self):
        """
        Charge la configuration en combinant :
        - les valeurs de config.OBSERVATORY (source de vérité)
        - les éventuelles valeurs surchargées depuis config.json
        """
        # Base : config.py
        eq = EQUIPMENT_OBSERVATION
        self.config_data = {
            "observatory": {
                "name": OBSERVATORY.get("name", ""),
                "latitude": float(OBSERVATORY.get("latitude", OBSERVATORY.get("lat", 0.0))),
                "longitude": float(OBSERVATORY.get("longitude", OBSERVATORY.get("lon", 0.0))),
                "elevation": float(OBSERVATORY.get("elevation", OBSERVATORY.get("elev", 0.0))),
                "timezone": OBSERVATORY.get("timezone", "UTC"),
            },
            "equipment": {
                "telescope_diameter_mm": eq.get("telescope_diameter_mm", ""),
                "focal_length_mm": eq.get("focal_length_mm", ""),
                "sensor_width_mm": eq.get("sensor_width_mm", ""),
                "sensor_height_mm": eq.get("sensor_height_mm", ""),
                "binning": eq.get("binning", 1),  # 1, 2, 3 ou 4 (1x1, 2x2, 3x3, 4x4)
                "pixel_size_um": eq.get("pixel_size_um", ""),
            },
            "astrometry_api_key": "",
            "affiliations": [],  # Liste de dicts: [{"text": "...", "selected": bool}]
        }

        # Surcharge par config.json si présent
        if self.config_file.exists():
            try:
                with open(self.config_file, "r", encoding="utf-8") as f:
                    saved_cfg = json.load(f)

                if "observatory" in saved_cfg:
                    self.config_data["observatory"].update(saved_cfg["observatory"])

                if "equipment" in saved_cfg:
                    self.config_data["equipment"].update(saved_cfg["equipment"])

                if "astrometry_api_key" in saved_cfg:
                    self.config_data["astrometry_api_key"] = saved_cfg["astrometry_api_key"]
                
                if "affiliations" in saved_cfg:
                    self.config_data["affiliations"] = saved_cfg["affiliations"]

            except Exception as e:
                messagebox.showwarning("Config", f"Erreur lecture config.json : {e}")

    def save_config(self):
        """
        Sauvegarde la config dans config.json
        et met à jour config.OBSERVATORY en mémoire.
        """
        # Met à jour le bloc observatory depuis les champs GUI
        obs = self.config_data.get("observatory", {})
        obs["name"] = self.entries["name"].get().strip()

        try:
            lat_str = self.entries["lat"].get().strip().replace(",", ".") or "0"
            lon_str = self.entries["lon"].get().strip().replace(",", ".") or "0"
            elev_str = self.entries["elev"].get().strip().replace(",", ".") or "0"
            obs["latitude"] = float(lat_str)
            obs["longitude"] = float(lon_str)
            obs["elevation"] = float(elev_str)
        except ValueError:
            messagebox.showerror("Erreur", "Latitude / Longitude / Élévation doivent être numériques.")
            return

        self.config_data["observatory"] = obs

        # Bloc équipement (chaînes pour affichage, numériques pour config)
        eq_cfg = self.config_data.get("equipment", {})
        for key in ("telescope_diameter_mm", "focal_length_mm", "sensor_width_mm", "sensor_height_mm", "pixel_size_um"):
            if key in self.equip_entries:
                val = self.equip_entries[key].get().strip().replace(",", ".")
                eq_cfg[key] = val if val else ""
        # Binning
        binning_str = (getattr(self, "binning_combo", None) and self.binning_combo.get()) or "1x1"
        try:
            binning = int(binning_str.split("x")[0]) if "x" in binning_str else 1
            binning = max(1, min(4, binning))
        except (ValueError, TypeError):
            binning = 1
        eq_cfg["binning"] = binning
        self.config_data["equipment"] = eq_cfg

        # Sauvegarder les affiliations
        self.config_data["affiliations"] = self.affiliations

        try:
            with open(self.config_file, "w", encoding="utf-8") as f:
                json.dump(self.config_data, f, indent=4)

            # Met à jour config.OBSERVATORY (en mémoire)
            config.OBSERVATORY = self.config_data["observatory"]
            # Met à jour config.EQUIPMENT_OBSERVATION (valeurs numériques + échelle calculée)
            eq = self.config_data.get("equipment", {})
            try:
                fl = float(eq.get("focal_length_mm") or 0)
                pu_native = float(eq.get("pixel_size_um") or 0)
                binning = int(eq.get("binning", 1))
                binning = max(1, min(4, binning))
                pu = pu_native * binning  # pixel size effectif pour l'échelle
                config.EQUIPMENT_OBSERVATION["focal_length_mm"] = fl
                config.EQUIPMENT_OBSERVATION["pixel_size_um"] = pu
                if fl > 0 and pu > 0:
                    config.EQUIPMENT_OBSERVATION["pixel_scale_arcsec"] = 206.265 * pu / fl
                for k in ("telescope_diameter_mm", "sensor_width_mm", "sensor_height_mm"):
                    v = eq.get(k, "")
                    try:
                        config.EQUIPMENT_OBSERVATION[k] = float(v) if (v and str(v).strip()) else EQUIPMENT_OBSERVATION.get(k, 0)
                    except (ValueError, TypeError):
                        config.EQUIPMENT_OBSERVATION[k] = EQUIPMENT_OBSERVATION.get(k, 0)
            except (ValueError, TypeError):
                pass
            # Met à jour config.AFFILIATIONS (en mémoire)
            config.AFFILIATIONS = self.config_data.get("affiliations", [])
            importlib.reload(config)

            messagebox.showinfo("Succès", "Configuration sauvegardée.")
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible de sauvegarder : {e}")

    # ------------------------------------------------------------------
    # GUI
    # ------------------------------------------------------------------
    def create_widgets(self):
        add_manual_help_header_grid(self, "2-accueil", columnspan=2, row=0)
        # Frame principal avec deux colonnes
        left_frame = ttk.Frame(self)
        left_frame.grid(row=1, column=0, sticky="nsew", padx=(0, 10))
        
        right_frame = ttk.Frame(self)
        right_frame.grid(row=1, column=1, sticky="nsew")
        
        # Configuration des colonnes pour l'expansion
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)
        
        # ========== COLONNE GAUCHE ==========
        # Observatoire
        ttk.Label(
            left_frame,
            text="Observatoire",
            font=("Helvetica", 12, "bold")
        ).grid(row=0, column=0, columnspan=2, pady=(0, 10), sticky="w")

        fields = [
            ("Nom", "name"),
            ("Latitude (°)", "lat"),
            ("Longitude (°)", "lon"),
            ("Élévation (m)", "elev"),
        ]
        self.entries = {}
        for i, (label, key) in enumerate(fields, 1):
            ttk.Label(left_frame, text=label + ":").grid(
                row=i, column=0, sticky="e", padx=(0, 5), pady=2
            )
            entry = ttk.Entry(left_frame, width=20)
            entry.grid(row=i, column=1, sticky="w", pady=2)
            self.entries[key] = entry

        # ---------- Configuration du matériel ----------
        ttk.Separator(left_frame, orient="horizontal").grid(
            row=5, column=0, columnspan=2, sticky="ew", pady=15
        )
        equip_lf = ttk.LabelFrame(left_frame, text="Configuration du matériel", padding=5)
        equip_lf.grid(row=6, column=0, columnspan=2, sticky="ew", pady=(0, 5))

        def _bind_equip_update(*args):
            self._update_equipment_computed()

        equip_fields = [
            ("Diamètre télescope (mm)", "telescope_diameter_mm"),
            ("Longueur focale (mm)", "focal_length_mm"),
            ("Largeur capteur (mm)", "sensor_width_mm"),
            ("Hauteur capteur (mm)", "sensor_height_mm"),
        ]
        for i, (label, key) in enumerate(equip_fields):
            ttk.Label(equip_lf, text=label + ":").grid(row=i, column=0, sticky="e", padx=(0, 5), pady=2)
            e = ttk.Entry(equip_lf, width=12)
            e.grid(row=i, column=1, sticky="w", pady=2)
            e.bind("<KeyRelease>", _bind_equip_update)
            self.equip_entries[key] = e

        # Binning (avant taille pixel) : 1x1, 2x2, 3x3, 4x4
        row_binning = len(equip_fields)
        ttk.Label(equip_lf, text="Binning:").grid(row=row_binning, column=0, sticky="e", padx=(0, 5), pady=2)
        self.binning_combo = ttk.Combobox(
            equip_lf, width=10, values=("1x1", "2x2", "3x3", "4x4"), state="readonly"
        )
        self.binning_combo.grid(row=row_binning, column=1, sticky="w", pady=2)
        self.binning_combo.set("1x1")
        self.binning_combo.bind("<<ComboboxSelected>>", lambda e: _bind_equip_update())

        row_pixel = row_binning + 1
        ttk.Label(equip_lf, text="Taille pixel (µm):").grid(row=row_pixel, column=0, sticky="e", padx=(0, 5), pady=2)
        e_pixel = ttk.Entry(equip_lf, width=12)
        e_pixel.grid(row=row_pixel, column=1, sticky="w", pady=2)
        e_pixel.bind("<KeyRelease>", _bind_equip_update)
        self.equip_entries["pixel_size_um"] = e_pixel

        row_fov = row_pixel + 1
        ttk.Label(equip_lf, text="Champ de vue (calculé):").grid(row=row_fov, column=0, sticky="e", padx=(0, 5), pady=2)
        self.fov_label = ttk.Label(equip_lf, text="— ° × — °", foreground="gray")
        self.fov_label.grid(row=row_fov, column=1, sticky="w", pady=2)

        ttk.Label(equip_lf, text="Échelle pixel (″/px):").grid(row=row_fov + 1, column=0, sticky="e", padx=(0, 5), pady=2)
        self.pixel_scale_label = ttk.Label(equip_lf, text="—", foreground="gray")
        self.pixel_scale_label.grid(row=row_fov + 1, column=1, sticky="w", pady=2)

        # Clé API
        ttk.Separator(left_frame, orient="horizontal").grid(
            row=7, column=0, columnspan=2, sticky="ew", pady=15
        )
        ttk.Label(
            left_frame,
            text="Clé API Astrometry.net",
            font=("Helvetica", 12, "bold")
        ).grid(row=8, column=0, columnspan=2, sticky="w")

        api_entry = ttk.Entry(left_frame, textvariable=self.api_var, width=40, show="*")
        api_entry.grid(row=9, column=0, columnspan=2, pady=5)

        api_btns = ttk.Frame(left_frame)
        api_btns.grid(row=10, column=0, columnspan=2, pady=5, sticky="ew")
        ttk.Button(api_btns, text="Générer/Modifier", command=self.set_api_key).pack(
            side=tk.LEFT, padx=(0, 8)
        )
        ttk.Button(api_btns, text="Voir / vérifier", command=self.verify_api_key_view).pack(
            side=tk.LEFT
        )

        # Affiliations
        ttk.Separator(left_frame, orient="horizontal").grid(
            row=11, column=0, columnspan=2, sticky="ew", pady=15
        )
        ttk.Label(
            left_frame,
            text="Affiliations",
            font=("Helvetica", 12, "bold")
        ).grid(row=12, column=0, columnspan=2, sticky="w")

        # Frame pour la liste des affiliations avec scrollbar
        aff_frame = ttk.Frame(left_frame)
        aff_frame.grid(row=13, column=0, columnspan=2, sticky="ew", pady=5)
        
        # Canvas avec scrollbar pour la liste des affiliations
        canvas_aff = tk.Canvas(aff_frame, height=150)
        scrollbar_aff = ttk.Scrollbar(aff_frame, orient="vertical", command=canvas_aff.yview)
        self.aff_scrollable_frame = ttk.Frame(canvas_aff)
        
        self.aff_scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas_aff.configure(scrollregion=canvas_aff.bbox("all"))
        )
        
        canvas_aff.create_window((0, 0), window=self.aff_scrollable_frame, anchor="nw")
        canvas_aff.configure(yscrollcommand=scrollbar_aff.set)
        
        canvas_aff.pack(side="left", fill="both", expand=True)
        scrollbar_aff.pack(side="right", fill="y")
        
        # Bind mousewheel to canvas
        def _on_mousewheel(event):
            canvas_aff.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas_aff.bind_all("<MouseWheel>", _on_mousewheel)
        
        # Liste des affiliations (stockée comme liste de dicts: {"text": "...", "selected": bool})
        self.affiliations = []
        self.aff_widgets = []  # Pour stocker les widgets (checkbox + entry + button)
        
        # Frame pour ajouter une nouvelle affiliation
        add_aff_frame = ttk.Frame(left_frame)
        add_aff_frame.grid(row=14, column=0, columnspan=2, sticky="ew", pady=5)
        
        ttk.Label(add_aff_frame, text="Nouvelle affiliation:").pack(side="left", padx=5)
        self.new_aff_entry = ttk.Entry(add_aff_frame, width=50)
        self.new_aff_entry.pack(side="left", padx=5, fill="x", expand=True)
        ttk.Button(
            add_aff_frame,
            text="Ajouter",
            command=self.add_affiliation
        ).pack(side="left", padx=5)

        # Boutons bas
        btn_frame = ttk.Frame(left_frame)
        btn_frame.grid(row=15, column=0, columnspan=2, pady=20)

        ttk.Button(
            btn_frame,
            text="Sauvegarder",
            command=self.save_config
        ).pack(side="left", padx=5)

        ttk.Button(
            btn_frame,
            text="Quitter",
            command=self.quit_app
        ).pack(side="left", padx=5)
        
        # ========== COLONNE DROITE ==========
        # Bibliographie
        ttk.Label(
            right_frame,
            text="Bibliographie",
            font=("Helvetica", 12, "bold")
        ).grid(row=0, column=0, sticky="w", pady=(0, 10))
        
        # Boutons actions bibliographie
        bib_actions = ttk.Frame(right_frame)
        bib_actions.grid(row=1, column=0, sticky="ew", pady=5)
        bib_actions.columnconfigure(0, weight=1)
        bib_actions.columnconfigure(1, weight=0)
        ttk.Button(
            bib_actions,
            text="📁 Ouvrir le répertoire Bibliographie",
            command=self.open_bibliography_folder
        ).grid(row=0, column=0, sticky="ew", padx=(0, 5))
        ttk.Button(
            bib_actions,
            text="🔄 Rafraîchir",
            command=self.scan_bibliography
        ).grid(row=0, column=1, sticky="e")
        
        # Séparateur
        ttk.Separator(right_frame, orient="horizontal").grid(
            row=2, column=0, sticky="ew", pady=10
        )
        
        # Label pour la sélection d'articles
        ttk.Label(
            right_frame,
            text="Articles disponibles:",
            font=("Helvetica", 10, "bold")
        ).grid(row=3, column=0, sticky="w", pady=(0, 5))
        
        # Frame pour la liste des articles avec scrollbar
        bib_frame = ttk.Frame(right_frame)
        bib_frame.grid(row=4, column=0, sticky="nsew", pady=5)
        right_frame.rowconfigure(4, weight=1)
        right_frame.columnconfigure(0, weight=1)
        
        # Canvas avec scrollbar pour la liste des articles
        canvas_bib = tk.Canvas(bib_frame, height=400)
        scrollbar_bib = ttk.Scrollbar(bib_frame, orient="vertical", command=canvas_bib.yview)
        self.bib_scrollable_frame = ttk.Frame(canvas_bib)
        
        self.bib_scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas_bib.configure(scrollregion=canvas_bib.bbox("all"))
        )
        
        canvas_bib.create_window((0, 0), window=self.bib_scrollable_frame, anchor="nw")
        canvas_bib.configure(yscrollcommand=scrollbar_bib.set)
        
        canvas_bib.pack(side="left", fill="both", expand=False)
        scrollbar_bib.pack(side="right", fill="y")
        
        # Bind mousewheel to canvas
        def _on_mousewheel_bib(event):
            canvas_bib.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas_bib.bind_all("<MouseWheel>", _on_mousewheel_bib)
        
        # Liste des articles (sera remplie par scan_bibliography)
        self.bibliography_articles = []
        self.bib_widgets = []
        self.category_frames = {}  # Pour stocker les frames des catégories et leur état (réduit/développé)
        
        # Scanner le répertoire bibliographie au démarrage
        self.scan_bibliography()

    def _update_equipment_computed(self):
        """Met à jour les champs calculés (champ de vue, échelle pixel) à partir des entrées matériel."""
        if not hasattr(self, "fov_label") or not self.fov_label.winfo_exists():
            return
        e = self.equip_entries
        if not e:
            return
        def _val(key):
            try:
                w = e.get(key)
                s = (w.get() if w else "").strip().replace(",", ".") or "0"
                return float(s)
            except (ValueError, TypeError, AttributeError):
                return 0.0
        try:
            fl = _val("focal_length_mm")
            sw = _val("sensor_width_mm")
            sh = _val("sensor_height_mm")
            pu = _val("pixel_size_um")
            # Binning : 1x1 -> 1, 2x2 -> 2, etc.
            binning_str = (getattr(self, "binning_combo", None) and self.binning_combo.get()) or "1x1"
            binning = int(binning_str.split("x")[0]) if "x" in binning_str else 1
            binning = max(1, min(4, binning))
            pu_effective = pu * binning  # taille de pixel effective pour l'échelle
        except Exception:
            self.fov_label.config(text="— ° × — °")
            self.pixel_scale_label.config(text="—")
            return
        deg_per_rad = 57.29577951308232
        if fl > 0:
            fov_w = (sw / fl) * deg_per_rad if sw > 0 else 0.0
            fov_h = (sh / fl) * deg_per_rad if sh > 0 else 0.0
            self.fov_label.config(text=f"{fov_w:.4f}° × {fov_h:.4f}°" if (sw and sh) else "— ° × — °")
            if pu_effective > 0:
                scale = 206.265 * pu_effective / fl  # µm effectif, fl en mm → ″/px
                self.pixel_scale_label.config(text=f"{scale:.4f} ″/px")
            else:
                self.pixel_scale_label.config(text="—")
        else:
            self.fov_label.config(text="— ° × — °")
            self.pixel_scale_label.config(text="—")

    def load_values(self):
        """
        Remplit les champs GUI à partir de self.config_data.
        """
        obs = self.config_data.get("observatory", {})

        self.obs_name_var.set(obs.get("name", ""))
        self.obs_lat_var.set(obs.get("latitude", 0.0))
        self.obs_lon_var.set(obs.get("longitude", 0.0))
        self.obs_elev_var.set(obs.get("elevation", 0.0))

        # Entries texte
        self.entries["name"].delete(0, tk.END)
        self.entries["name"].insert(0, self.obs_name_var.get())

        self.entries["lat"].delete(0, tk.END)
        self.entries["lat"].insert(0, str(self.obs_lat_var.get()))

        self.entries["lon"].delete(0, tk.END)
        self.entries["lon"].insert(0, str(self.obs_lon_var.get()))

        self.entries["elev"].delete(0, tk.END)
        self.entries["elev"].insert(0, str(self.obs_elev_var.get()))

        # Clé API (affichage masqué)
        key = self.config_data.get("astrometry_api_key", "")
        if key:
            self.api_var.set(key[:10] + "...")
        else:
            self.api_var.set("")

        # Configuration matériel
        eq = self.config_data.get("equipment", {})
        for key, entry in (self.equip_entries or {}).items():
            val = eq.get(key, "")
            entry.delete(0, tk.END)
            entry.insert(0, str(val) if val not in (None, "") else "")
        # Binning
        if getattr(self, "binning_combo", None):
            b = eq.get("binning", 1)
            try:
                b = int(b)
                b = max(1, min(4, b))
            except (ValueError, TypeError):
                b = 1
            self.binning_combo.set(f"{b}x{b}")
        self._update_equipment_computed()
        
        # Charger les affiliations
        self.affiliations = self.config_data.get("affiliations", [])
        if hasattr(self, 'aff_scrollable_frame'):
            self.refresh_affiliations_display()

    # ------------------------------------------------------------------
    # Clé API Astrometry.net
    # ------------------------------------------------------------------
    def set_api_key(self):
        dialog = APIKeyDialog(self)
        api_key = getattr(dialog, "api_key", None)
        if api_key:
            self.config_data["astrometry_api_key"] = api_key
            self.api_var.set(api_key[:10] + "...")
            # Sauvegarde dans le fichier dédié
            try:
                with open(API_KEY_PATH, "w", encoding="utf-8") as f:
                    f.write(api_key.strip())
            except Exception as e:
                messagebox.showwarning(
                    "Clé API",
                    f"Impossible d'écrire la clé API dans {API_KEY_PATH} : {e}"
                )
            # Sauvegarde globale config.json
            self.save_config()

    def verify_api_key_view(self):
        """
        Affiche le texte réellement utilisé : ~/.astrometry_api_key si non vide,
        sinon repli sur astrometry_api_key de config.json.
        """
        raw = ""
        source = "(aucune)"
        err_read = None
        if API_KEY_PATH.exists():
            try:
                raw = API_KEY_PATH.read_text(encoding="utf-8")
                if raw.strip():
                    source = str(API_KEY_PATH)
                else:
                    raw = ""
            except OSError as e:
                err_read = str(e)
                raw = ""
        if not raw.strip():
            raw = self.config_data.get("astrometry_api_key", "") or ""
            source = str(self.config_file) if raw.strip() else source

        win = tk.Toplevel(self)
        win.title("Clé API Astrometry.net — contenu enregistré")
        win.geometry("520x360")

        ttk.Label(win, text=f"Source : {source}", font=("", 9)).pack(anchor="w", padx=10, pady=(10, 4))
        if err_read:
            ttk.Label(win, text=f"Lecture fichier clé : {err_read}", foreground="orange").pack(
                anchor="w", padx=10, pady=2
            )

        ttk.Label(win, text="Texte enregistré (tel quel) :").pack(anchor="w", padx=10, pady=(8, 2))
        txt = tk.Text(win, height=6, width=64, wrap="none", font=("Consolas", 10))
        txt.pack(fill=tk.BOTH, expand=True, padx=10, pady=4)
        txt.insert("1.0", raw if raw else "(vide)")
        txt.configure(state="disabled")

        n = len(raw)
        nstrip = len(raw.strip())
        ascii_ok = raw.isascii() if raw else True
        extra_ws = (raw != raw.strip()) if raw else False
        checks = ttk.Frame(win)
        checks.pack(fill=tk.X, padx=10, pady=6)
        ttk.Label(checks, text=f"Longueur (caractères) : {n}").pack(anchor="w")
        ttk.Label(checks, text=f"Longueur après strip() : {nstrip}").pack(anchor="w")
        ttk.Label(
            checks,
            text=f"ASCII uniquement : {'oui' if ascii_ok else 'non'}",
            foreground=("gray" if ascii_ok else "darkred"),
        ).pack(anchor="w")
        ttk.Label(
            checks,
            text=f"Espaces / retours en tête ou fin : {'oui (à nettoyer)' if extra_ws else 'non'}",
            foreground=("darkred" if extra_ws else "gray"),
        ).pack(anchor="w")

        ttk.Button(win, text="Fermer", command=win.destroy).pack(pady=10)

    # ------------------------------------------------------------------
    # Affiliations
    # ------------------------------------------------------------------
    def add_affiliation(self):
        """Ajoute une nouvelle affiliation à la liste."""
        text = self.new_aff_entry.get().strip()
        if not text:
            messagebox.showwarning("Attention", "Veuillez entrer une affiliation.")
            return
        
        # Vérifier si elle n'existe pas déjà
        if any(aff.get("text", "") == text for aff in self.affiliations):
            messagebox.showwarning("Attention", "Cette affiliation existe déjà.")
            return
        
        # Ajouter avec selected=False par défaut
        self.affiliations.append({"text": text, "selected": False})
        self.new_aff_entry.delete(0, tk.END)
        self.refresh_affiliations_display()
    
    def remove_affiliation(self, index):
        """Supprime une affiliation de la liste."""
        if 0 <= index < len(self.affiliations):
            self.affiliations.pop(index)
            self.refresh_affiliations_display()
    
    def toggle_affiliation_selection(self, index):
        """Change l'état de sélection d'une affiliation."""
        if 0 <= index < len(self.affiliations):
            # Désélectionner toutes les autres
            for i, aff in enumerate(self.affiliations):
                if i == index:
                    self.affiliations[i]["selected"] = True
                else:
                    self.affiliations[i]["selected"] = False
            self.refresh_affiliations_display()
    
    def refresh_affiliations_display(self):
        """Rafraîchit l'affichage de la liste des affiliations."""
        # Supprimer tous les widgets existants
        for widget in self.aff_widgets:
            widget.destroy()
        self.aff_widgets.clear()
        
        # Recréer les widgets pour chaque affiliation
        for i, aff in enumerate(self.affiliations):
            row_frame = ttk.Frame(self.aff_scrollable_frame)
            row_frame.pack(fill="x", padx=5, pady=2)
            
            # Checkbox pour sélectionner
            var = tk.BooleanVar(value=aff.get("selected", False))
            checkbox = ttk.Checkbutton(
                row_frame,
                variable=var,
                command=lambda idx=i: self.toggle_affiliation_selection(idx)
            )
            checkbox.pack(side="left", padx=5)
            
            # Label avec le texte de l'affiliation
            label = ttk.Label(row_frame, text=aff.get("text", ""), width=60, anchor="w")
            label.pack(side="left", padx=5, fill="x", expand=True)
            
            # Bouton supprimer
            btn_remove = ttk.Button(
                row_frame,
                text="✕",
                width=3,
                command=lambda idx=i: self.remove_affiliation(idx)
            )
            btn_remove.pack(side="right", padx=5)
            
            self.aff_widgets.extend([row_frame, checkbox, label, btn_remove])
    
    def get_selected_affiliation(self):
        """Retourne l'affiliation sélectionnée (ou None)."""
        for aff in self.affiliations:
            if aff.get("selected", False):
                return aff.get("text", "")
        return None

    # ------------------------------------------------------------------
    # Bibliographie
    # ------------------------------------------------------------------
    def open_bibliography_folder(self):
        """Ouvre le répertoire Bibliographie dans l'explorateur de fichiers."""
        bib_dir = self.base_dir / "Bibliographie"
        if bib_dir.exists():
            try:
                if os.name == 'nt':  # Windows
                    os.startfile(str(bib_dir))
                elif os.name == 'posix':  # macOS et Linux
                    import subprocess
                    subprocess.Popen(['xdg-open' if os.uname().sysname == 'Linux' else 'open', str(bib_dir)])
            except Exception as e:
                messagebox.showerror("Erreur", f"Impossible d'ouvrir le répertoire : {e}")
        else:
            messagebox.showwarning("Attention", f"Le répertoire Bibliographie n'existe pas : {bib_dir}")
    
    def scan_bibliography(self):
        """Scanne le répertoire Bibliographie et liste tous les fichiers PDF."""
        bib_dir = self.base_dir / "Bibliographie"
        self.bibliography_articles = []
        
        if not bib_dir.exists():
            return
        
        # Parcourir récursivement tous les sous-dossiers
        for pdf_file in bib_dir.rglob("*.pdf"):
            # Obtenir le chemin relatif depuis Bibliographie
            rel_path = pdf_file.relative_to(bib_dir)
            # Obtenir le dossier parent (catégorie)
            category = rel_path.parent.name if rel_path.parent != Path(".") else "Racine"
            # Nom du fichier sans extension
            name = pdf_file.stem
            
            self.bibliography_articles.append({
                "name": name,
                "category": category,
                "path": str(pdf_file),
                "selected": False
            })
        
        # Trier par catégorie puis par nom
        self.bibliography_articles.sort(key=lambda x: (x["category"], x["name"]))
        
        # Rafraîchir l'affichage
        self.refresh_bibliography_display()
    
    def refresh_bibliography_display(self):
        """Rafraîchit l'affichage de la liste des articles."""
        # Supprimer tous les widgets existants
        for widget in self.bib_widgets:
            widget.destroy()
        self.bib_widgets.clear()
        self.category_frames.clear()
        
        if not self.bibliography_articles:
            no_articles_label = ttk.Label(
                self.bib_scrollable_frame,
                text="Aucun article PDF trouvé dans le répertoire Bibliographie.",
                foreground="gray"
            )
            no_articles_label.pack(pady=20)
            self.bib_widgets.append(no_articles_label)
            return
        
        # Grouper par catégorie
        articles_by_category = {}
        for article_idx, article in enumerate(self.bibliography_articles):
            category = article["category"]
            if category not in articles_by_category:
                articles_by_category[category] = []
            articles_by_category[category].append((article_idx, article))
        
        # Créer l'affichage pour chaque catégorie avec possibilité de réduire/développer
        for category in sorted(articles_by_category.keys()):
            # Frame pour la catégorie (header + contenu)
            category_container = ttk.Frame(self.bib_scrollable_frame)
            category_container.pack(fill="x", padx=5, pady=2)
            
            # Header de la catégorie (cliquable pour réduire/développer)
            category_header = ttk.Frame(category_container)
            category_header.pack(fill="x")
            
            # Bouton pour réduire/développer (état initial: réduit)
            is_expanded = self.category_frames.get(category, {}).get("expanded", False)
            expand_btn = ttk.Button(
                category_header,
                text="▼" if is_expanded else "▶",
                width=3,
                command=lambda cat=category: self.toggle_category(cat)
            )
            expand_btn.pack(side="left", padx=2)
            
            # Label avec le nom de la catégorie
            category_label = ttk.Label(
                category_header,
                text=f"📁 {category} ({len(articles_by_category[category])} article{'s' if len(articles_by_category[category]) > 1 else ''})",
                font=("Helvetica", 10, "bold"),
                foreground="blue",
                cursor="hand2"
            )
            category_label.pack(side="left", padx=5)
            # Clic sur le label pour aussi réduire/développer
            category_label.bind("<Button-1>", lambda e, cat=category: self.toggle_category(cat))
            
            # Frame pour le contenu de la catégorie (articles)
            category_content = ttk.Frame(category_container)
            if is_expanded:
                category_content.pack(fill="x", padx=(20, 0))
            else:
                category_content.pack_forget()  # Caché si réduit
            
            # Stocker les références pour pouvoir les afficher/cacher
            self.category_frames[category] = {
                "container": category_container,
                "header": category_header,
                "content": category_content,
                "expand_btn": expand_btn,
                "expanded": is_expanded
            }
            
            # Ajouter les articles de cette catégorie
            for article_idx, article in articles_by_category[category]:
                # Frame pour chaque article
                article_frame = ttk.Frame(category_content)
                article_frame.pack(fill="x", padx=5, pady=2)
                
                # Checkbox pour sélectionner
                var = tk.BooleanVar(value=article.get("selected", False))
                checkbox = ttk.Checkbutton(
                    article_frame,
                    variable=var,
                    command=lambda idx=article_idx: self.toggle_article_selection(idx)
                )
                checkbox.pack(side="left", padx=5)
                
                # Label avec le nom de l'article (tronqué si trop long)
                article_name = article["name"]
                if len(article_name) > 55:
                    article_name = article_name[:52] + "..."
                label = ttk.Label(article_frame, text=article_name, anchor="w", cursor="hand2")
                label.pack(side="left", padx=5, fill="x", expand=True)
                # Double-clic pour ouvrir le PDF
                label.bind("<Double-Button-1>", lambda e, path=article["path"]: self.open_pdf(path))
                
                # Bouton pour ouvrir le PDF
                btn_open = ttk.Button(
                    article_frame,
                    text="📄",
                    width=3,
                    command=lambda path=article["path"]: self.open_pdf(path)
                )
                btn_open.pack(side="right", padx=5)
                
                self.bib_widgets.extend([article_frame, checkbox, label, btn_open])
            
            self.bib_widgets.extend([category_container, category_header, category_content])
    
    def toggle_category(self, category):
        """Réduit ou développe une catégorie."""
        if category not in self.category_frames:
            return
        
        category_info = self.category_frames[category]
        is_expanded = category_info["expanded"]
        
        # Inverser l'état
        category_info["expanded"] = not is_expanded
        
        # Mettre à jour le bouton
        category_info["expand_btn"].config(text="▼" if not is_expanded else "▶")
        
        # Afficher/cacher le contenu
        if not is_expanded:
            category_info["content"].pack(fill="x", padx=(20, 0))
        else:
            category_info["content"].pack_forget()
    
    def toggle_article_selection(self, article_index):
        """Change l'état de sélection d'un article."""
        if 0 <= article_index < len(self.bibliography_articles):
            self.bibliography_articles[article_index]["selected"] = not self.bibliography_articles[article_index].get("selected", False)
            self.refresh_bibliography_display()
    
    def open_pdf(self, pdf_path):
        """Ouvre un fichier PDF avec l'application par défaut."""
        try:
            if os.name == 'nt':  # Windows
                os.startfile(pdf_path)
            elif os.name == 'posix':  # macOS et Linux
                import subprocess
                subprocess.Popen(['xdg-open' if os.uname().sysname == 'Linux' else 'open', pdf_path])
        except Exception as e:
            messagebox.showerror("Erreur", f"Impossible d'ouvrir le fichier PDF : {e}")
    
    def get_selected_articles(self):
        """Retourne la liste des articles sélectionnés."""
        return [art for art in self.bibliography_articles if art.get("selected", False)]

    # ------------------------------------------------------------------
    # Quitter
    # ------------------------------------------------------------------
    def quit_app(self):
        if messagebox.askokcancel("Quitter", "Voulez-vous quitter NPOAP ?"):
            self.winfo_toplevel().destroy()
            os._exit(0)