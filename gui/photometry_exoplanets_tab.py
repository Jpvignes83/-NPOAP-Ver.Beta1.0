# gui/photometry_exoplanets_tab.py
import json
import os
import sys
import shutil
import zipfile
import subprocess
import importlib.util
import tempfile
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from threading import Thread
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u
from io import StringIO


def _read_text_file_robust(file_path):
    """
    Lit un fichier texte en essayant plusieurs encodages (utf-8, utf-16, utf-8-sig, cp1252).
    Évite l'erreur 'utf-8' codec can't decode byte 0xff in position 0 lorsque le fichier
    a un BOM ou un autre encodage (ex. final-report.txt généré ailleurs).
    """
    encodings = ('utf-8', 'utf-16', 'utf-16-le', 'utf-16-be', 'utf-8-sig', 'cp1252', 'latin-1')
    for enc in encodings:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                return f.read()
        except (UnicodeDecodeError, UnicodeError):
            continue
    # Dernier recours : lire en binaire et décoder en ignorant les erreurs
    with open(file_path, 'rb') as f:
        return f.read().decode('utf-8', errors='replace')


def _to_utf8_no_bom(s):
    """Retourne les octets UTF-8 sans BOM pour upload ExoClock (évite erreur 0xff / UTF-16)."""
    if isinstance(s, bytes):
        s = s.decode('utf-8', errors='replace')
    if s.startswith('\ufeff'):
        s = s[1:]
    return s.encode('utf-8', errors='replace')


def _read_csv_robust(file_path):
    """
    Lit un CSV avec plusieurs encodages (évite 'utf-8' codec can't decode byte 0xff).
    Force les colonnes de temps (JD-UTC, JD_UTC, BJD-TDB) en float pour une échelle correcte.
    Retourne un DataFrame pandas.
    """
    content = _read_text_file_robust(file_path)
    df = pd.read_csv(StringIO(content))
    for col in ('JD-UTC', 'JD_UTC', 'BJD-TDB'):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


# --- Imports du pipeline ---
from core.photometry_pipeline import PhotometryPipeline
from gui.target_selector import launch_target_selection
import config

# --- NOUVEAU : On utilise UNIQUEMENT le visualiseur complet ---
from gui.lightcurve_fitting import LightcurveFitting

class PhotometryExoplanetsTab(ttk.Frame):
    """
    Onglet de photométrie exoplanètes.
    Workflow en 6 étapes (référence, cibles, ouvertures, batch, analyse, rapports).
    """

    def __init__(self, parent, base_dir=None):
        super().__init__(parent, padding=0)
        self.base_dir = str(base_dir) if base_dir is not None else None
        
        # Instanciation du pipeline
        self.pipeline = PhotometryPipeline()

        self.ref_image = None
        self.target_coord = None
        self.comp_coords = []
        self.current_selections = None

        self.create_widgets()

    def create_widgets(self):
        right_frame = ttk.Frame(self, padding=10)
        right_frame.pack(fill=tk.BOTH, expand=True)

        hops_frame = ttk.LabelFrame(right_frame, text="HOPS", padding=8)
        hops_frame.pack(fill="both", expand=True, pady=(0, 2))

        hops_btn_row = ttk.Frame(hops_frame)
        hops_btn_row.pack(fill="x")
        ttk.Button(
            hops_btn_row,
            text="Installer / Réinstaller HOPS (ZIP)",
            command=self.install_hops_from_zip,
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            hops_btn_row,
            text="Lancer HOPS",
            command=self.launch_hops,
        ).pack(side="left")
        ttk.Button(
            hops_btn_row,
            text="Réinitialiser HOPS",
            command=self.reset_hops,
        ).pack(side="left", padx=(6, 0))

        self.hops_embed_host = ttk.Frame(hops_frame)
        self.hops_embed_host.pack(fill="both", expand=True, pady=(8, 0))
        # Zone à défilement : HOPS (grille + figures) peut dépasser la hauteur du panneau.
        self.hops_canvas = None
        self.hops_scrollbar = None
        self.hops_canvas_window_id = None
        self.hops_isolated_container = None
        self.hops_app = None
        self._build_hops_scrollable_embed()

        # Conservé pour compatibilité interne, mais non affiché.
        self.hops_status_var = tk.StringVar(value=self._hops_status_text())
        ttk.Label(
            hops_frame,
            text="MIT License\nCopyright (c) 2017 Angelos Tsiaras and Konstantinos Karpouzas\nHOPS-modified integration: J.P Vignes <jeanpascal.vignes@gmail.com>",
            justify="left",
        ).pack(anchor="w", pady=(6, 0))
        
        # Variables conservées pour compatibilité avec les fonctions batch appelées ailleurs.
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(hops_frame, variable=self.progress_var, maximum=100)
        self.lbl_progress = ttk.Label(hops_frame, text="En attente...", font=("Arial", 8))

    def _build_hops_scrollable_embed(self) -> None:
        scroll_wrap = ttk.Frame(self.hops_embed_host)
        scroll_wrap.pack(fill="both", expand=True)
        try:
            bg = scroll_wrap.cget("background")
        except tk.TclError:
            bg = None
        if not bg:
            bg = "#f5f5f5" if sys.platform == "win32" else "#ececec"
        self.hops_canvas = tk.Canvas(scroll_wrap, highlightthickness=0, borderwidth=0, bg=bg)
        self.hops_scrollbar = ttk.Scrollbar(scroll_wrap, orient="vertical", command=self.hops_canvas.yview)
        self.hops_canvas.configure(yscrollcommand=self.hops_scrollbar.set)
        self.hops_scrollbar.pack(side="right", fill="y")
        self.hops_canvas.pack(side="left", fill="both", expand=True)

        self.hops_isolated_container = ttk.Frame(self.hops_canvas)
        self.hops_canvas_window_id = self.hops_canvas.create_window(
            (0, 0), window=self.hops_isolated_container, anchor="nw"
        )
        self.hops_isolated_container.bind("<Configure>", self._on_hops_embed_inner_configure)
        self.hops_canvas.bind("<Configure>", self._on_hops_embed_canvas_configure)
        self._bind_hops_embed_mousewheel(self.hops_canvas)
        self._bind_hops_embed_mousewheel(self.hops_isolated_container)

    def _on_hops_embed_inner_configure(self, event) -> None:
        if self.hops_canvas is None:
            return
        try:
            self.hops_canvas.configure(scrollregion=self.hops_canvas.bbox("all"))
        except tk.TclError:
            pass

    def _on_hops_embed_canvas_configure(self, event) -> None:
        if self.hops_canvas is None or self.hops_canvas_window_id is None:
            return
        try:
            self.hops_canvas.itemconfigure(self.hops_canvas_window_id, width=event.width)
        except tk.TclError:
            pass

    def _on_hops_embed_mousewheel(self, event) -> None:
        if self.hops_canvas is None:
            return
        if getattr(event, "delta", 0):
            self.hops_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        elif event.num == 4:
            self.hops_canvas.yview_scroll(-1, "units")
        elif event.num == 5:
            self.hops_canvas.yview_scroll(1, "units")

    def _bind_hops_embed_mousewheel(self, widget) -> None:
        widget.bind("<MouseWheel>", self._on_hops_embed_mousewheel)
        widget.bind("<Button-4>", self._on_hops_embed_mousewheel)
        widget.bind("<Button-5>", self._on_hops_embed_mousewheel)

    def _recreate_hops_isolated_container(self) -> None:
        """Recrée la surface Tk accueillie dans le canvas (après fermeture HOPS ou réinitialisation)."""
        if self.hops_canvas is None:
            self._build_hops_scrollable_embed()
            return
        if self.hops_app is not None:
            try:
                self.hops_app.close()
            except Exception:
                pass
            self.hops_app = None
        if self.hops_canvas_window_id is not None:
            try:
                self.hops_canvas.delete(self.hops_canvas_window_id)
            except tk.TclError:
                pass
            self.hops_canvas_window_id = None
        self.hops_isolated_container = ttk.Frame(self.hops_canvas)
        self.hops_canvas_window_id = self.hops_canvas.create_window(
            (0, 0), window=self.hops_isolated_container, anchor="nw"
        )
        self.hops_isolated_container.bind("<Configure>", self._on_hops_embed_inner_configure)
        self._bind_hops_embed_mousewheel(self.hops_isolated_container)

    def browse_ref(self):
        initial = self.base_dir or os.getcwd()
        path = filedialog.askopenfilename(
            title="Sélectionner l'image de référence (Solved)",
            initialdir=initial,
            filetypes=[("FITS files", "*.fits *.fit *.fts")]
        )
        if path:
            self.ref_image = path
            self.lbl_ref.config(text=Path(path).name, foreground="green")
            self.btn_select.config(state="normal")
            logging.info(f"Image de référence : {path}")

    def launch_selector(self):
        if not self.ref_image: return
        
        def read_header_target_coord(fits_path):
            try:
                with fits.open(fits_path) as hdul:
                    header = hdul[0].header
                ra_h = header.get("OBJCTRA") or header.get("RA")
                dec_h = header.get("OBJCTDEC") or header.get("DEC")
                if ra_h and dec_h:
                    try:
                        ra_val = float(ra_h)
                        dec_val = float(dec_h)
                        if ra_val > 24.0:
                            return SkyCoord(ra_val, dec_val, unit=(u.deg, u.deg))
                    except Exception:
                        pass
                    try:
                        return SkyCoord(ra_h, dec_h, unit=(u.hourangle, u.deg))
                    except Exception:
                        return SkyCoord(float(ra_h), float(dec_h), unit=(u.deg, u.deg))
            except Exception as e:
                logging.warning(f"Coordonnées du header indisponibles: {e}")
            return None

        def on_done(fits_path, target_data, comps_list):
            # target_data peut être un dict avec 'coord' et optionnellement 'fwhm', ou directement une SkyCoord
            if isinstance(target_data, dict):
                t1 = target_data['coord']
            else:
                t1 = target_data
            
            self.target_coord = t1
            self.comp_coords = comps_list 
            self.current_selections = []
            self.current_selections.append({'label': 'T1', 'coord': t1, 'r_ap': 8.0, 'r_in': 12.0, 'r_out': 18.0})
            for i, c in enumerate(comps_list):
                self.current_selections.append({'label': f'C{i+1}', 'coord': c, 'r_ap': 8.0, 'r_in': 12.0, 'r_out': 18.0})
            self.lbl_targets.config(text=f"T1: OK\nComps: {len(comps_list)}")
            self.btn_aperture.config(state="normal")
            self.lbl_aperture_status.config(text="Prêt à régler", foreground="orange")
            logging.info("Sélection terminée → Étape 3 disponible")
            self.open_aperture_dialog()  # Fluidité

        header_target_coord = read_header_target_coord(self.ref_image)
        launch_target_selection(self.ref_image, on_selection_done=on_done, header_target_coord=header_target_coord)

    def open_aperture_dialog(self):
        if not self.ref_image or not self.target_coord: return
        try:
            self.pipeline.launch_photometry_aperture(
                fits_path=self.ref_image,
                target_coord=self.target_coord,
                comp_coords=self.comp_coords,
                on_finish=self.on_aperture_confirmed
            )
        except Exception as e:
            logging.error(f"Erreur ouverture interface aperture: {e}")
            messagebox.showerror("Erreur", f"Interface aperture: {e}")

    def on_aperture_confirmed(self, final_selections):
        if not final_selections:
            return
        self.current_selections = final_selections
        logging.info(f"Ouvertures validées pour {len(final_selections)} étoiles")
        
        t1 = next((s for s in final_selections if s['label'] == 'T1'), None)
        if not t1:
            messagebox.showwarning("Attention", "T1 est désactivée ! Réactivez-la.")
            self.lbl_aperture_status.config(text="T1 manquante", foreground="red")
            self.btn_batch.config(state="disabled")
            return
            
        self.target_coord = t1['coord']
        self.comp_coords = [s['coord'] for s in final_selections if s['label'] != 'T1']
        self.lbl_aperture_status.config(text="Validé", foreground="green")
        self.btn_batch.config(state="normal")

    def launch_batch(self):
        if not self.current_selections:
            messagebox.showwarning("Erreur", "Validez d'abord les ouvertures (Étape 3)")
            return

        initial = self.base_dir or getattr(config, "CALIBRATED_DIR", os.getcwd())
        science_dir = filedialog.askdirectory(title="Dossier images science", initialdir=initial)
        if not science_dir: return

        logging.info(f"Batch démarré sur {science_dir}")
        self.progress_var.set(0)
        self.lbl_progress.config(text="Traitement en cours...")

        Thread(target=self.run_batch, args=(science_dir,), daemon=True).start()

    def run_batch(self, science_dir):
        try:
            self.pipeline.process_photometry_series(
                folder=science_dir,
                target_coord=self.target_coord,
                comps=self.comp_coords,
                ref_image=self.ref_image,
                selections=self.current_selections,
                progress_callback=self.update_progress_ui
            )
            res_csv = Path(science_dir) / "photometrie" / "results.csv"
            self.after(0, lambda: self.on_batch_success(res_csv))
        except Exception as e:
            logging.error(f"Batch échoué: {e}")
            self.after(0, lambda: messagebox.showerror("Erreur Batch", str(e)))

    def update_progress_ui(self, percent):
        self.after(0, lambda: self.progress_var.set(percent))
        self.after(0, lambda: self.lbl_progress.config(text=f"{percent:.1f}%"))

    def on_batch_success(self, csv_path):
        self.lbl_progress.config(text="Terminé !")
        if messagebox.askyesno("Succès", "Photométrie terminée !\nOuvrir l'analyseur avancé ?"):
            self.open_advanced_analysis(csv_path)

    def open_advanced_analysis(self, csv_path=None):
        """
        Ouvre le visualiseur complet (detrending + modélisation)
        Si csv_path=None → demande à l'utilisateur
        """
        if csv_path is None:
            initial = self.base_dir or os.getcwd()
            csv_path = filedialog.askopenfilename(
                title="Ouvrir results.csv",
                initialdir=initial,
                filetypes=[("CSV", "*.csv")]
            )
        if not csv_path or not Path(csv_path).exists():
            return

        try:
            window = LightcurveFitting(self)
            window.load_from_path(csv_path)
            logging.info(f"Analyse avancée ouverte : {csv_path}")
        except Exception as e:
            logging.error(f"Erreur ouverture analyse : {e}")
            messagebox.showerror("Erreur", f"Impossible d'ouvrir le fichier :\n{e}")

    def save_equipment_config(self):
        """Sauvegarde les spécifications matériel dans config.py"""
        try:
            import config
            
            # Préparer le dictionnaire
            equipment_data = {
                "obs_code": self.equip_obs_code.get().strip(),
                "camera": self.equip_camera.get().strip(),
                "binning": self.equip_binning.get().strip(),
                "delim": self.equip_delim.get().strip()
            }
            
            # Lire le fichier config.py
            config_path = Path(__file__).parent.parent / "config.py"
            with open(config_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Trouver et remplacer la section EQUIPMENT_OBSERVATION
            import re
            pattern = r'(EQUIPMENT_OBSERVATION\s*=\s*\{)(.*?)(\})'
            
            new_section = (
                "EQUIPMENT_OBSERVATION = {\n"
                f'    "obs_code": {json.dumps(equipment_data["obs_code"])},           # Code observateur AAVSO (5 caractères max)\n'
                f'    "camera": {json.dumps(equipment_data["camera"])},             # Nom de la caméra\n'
                f'    "binning": {json.dumps(equipment_data["binning"])},         # Binning (1x1, 2x2, 3x3, 4x4)\n'
                f'    "delim": {json.dumps(equipment_data["delim"])},             # Délimiteur pour rapports (, ; | : ! / ? ou tab)\n'
                "}"
            )
            
            if re.search(pattern, content, re.DOTALL):
                content = re.sub(pattern, new_section, content, flags=re.DOTALL)
            else:
                # Ajouter avant la dernière ligne si pas trouvé
                content += "\n\n" + new_section + "\n"
            
            # Écrire le fichier
            with open(config_path, 'w', encoding='utf-8') as f:
                f.write(content)
            
            # Mettre à jour le module config
            config.EQUIPMENT_OBSERVATION = equipment_data
            
            messagebox.showinfo("Succès", "Spécifications matériel sauvegardées dans config.py")
            logging.info("Spécifications matériel sauvegardées")
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors de la sauvegarde:\n{e}")
            logging.error(f"Erreur sauvegarde équipement: {e}")

    def _npoap_root(self) -> Path:
        return Path(__file__).resolve().parent.parent

    def _hops_base_dir(self) -> Path:
        return self._npoap_root() / "external_apps" / "hops"

    def _hops_source_dir(self) -> Path:
        return self._hops_base_dir() / "hops-master"

    def _hops_entry_file(self) -> Path:
        return self._hops_source_dir() / "hops" / "__main__.py"

    def _hops_status_text(self) -> str:
        if self._hops_entry_file().is_file():
            return f"HOPS prêt : {self._hops_source_dir()}"
        return "HOPS non installé (utilisez le bouton d'installation depuis ZIP)."

    def _default_hops_zip(self) -> Path:
        # Priorité à l'archive fournie dans la distribution NPOAP.
        bundled = self._hops_base_dir() / "HOPS-modified.zip"
        if bundled.is_file():
            return bundled
        return Path.home() / "Downloads" / "hops-master.zip"

    def _ensure_hops_dependencies(self) -> None:
        """Installe les dépendances HOPS manquantes dans l'environnement Python courant."""
        # Le premier import qui casse est généralement `exoclock`.
        if importlib.util.find_spec("exoclock") is not None:
            return

        req = self._hops_source_dir() / "requirements.txt"
        if not req.is_file():
            raise FileNotFoundError(f"requirements.txt introuvable: {req}")

        self.hops_status_var.set("Installation dépendances HOPS (exoclock, etc.)...")
        self.update_idletasks()

        progress_dialog = tk.Toplevel(self)
        progress_dialog.title("HOPS - Installation dépendances")
        progress_dialog.geometry("460x150")
        progress_dialog.transient(self.winfo_toplevel())
        progress_dialog.grab_set()
        progress_dialog.resizable(False, False)

        progress_label = ttk.Label(
            progress_dialog,
            text="Installation des dépendances HOPS en cours...\nMerci de patienter.",
            justify="center",
        )
        progress_label.pack(pady=(18, 10), padx=12)

        progress_bar = ttk.Progressbar(progress_dialog, mode="indeterminate", length=360)
        progress_bar.pack(pady=(0, 10))
        progress_bar.start(10)

        ttk.Label(
            progress_dialog,
            text="(L'installation peut durer plusieurs minutes.)",
            foreground="gray",
        ).pack()

        progress_dialog.update_idletasks()
        logging.info(f"HOPS: installation dépendances depuis {req}")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", str(req)])
        finally:
            progress_bar.stop()
            try:
                progress_dialog.grab_release()
            except Exception:
                pass
            progress_dialog.destroy()

    def _register_gaia_filters_from_resources(self) -> None:
        """
        Enregistre automatiquement les passbands Gaia de resources/filters
        dans hops.pylightcurve41 au lancement de HOPS.
        """
        filters_dir = self._npoap_root() / "resources" / "filters"
        mapping = {
            "GAIA_G": filters_dir / "GAIA_G.txt",
            "GAIA_BP": filters_dir / "GAIA_BP.txt",
            "GAIA_RP": filters_dir / "GAIA_RP.txt",
        }
        try:
            import hops.pylightcurve41 as plc
        except Exception as e:
            logging.warning(f"Impossible d'importer hops.pylightcurve41 pour enregistrer les filtres Gaia: {e}")
            return

        try:
            existing = set(plc.all_filters())
        except Exception as e:
            logging.warning(f"Impossible de lister les filtres pylightcurve (all_filters): {e}")
            existing = set()

        for filter_name, filter_path in mapping.items():
            if not filter_path.is_file():
                logging.warning(f"Passband Gaia manquant: {filter_path}")
                continue
            if filter_name in existing:
                logging.info(f"Passband Gaia deja disponible: {filter_name}")
                continue
            try:
                plc.add_filter(filter_name, str(filter_path))
                logging.info(f"Passband Gaia enregistre: {filter_name} <- {filter_path}")
            except Exception as e:
                logging.warning(f"Echec enregistrement passband Gaia {filter_name}: {e}")

    def _enable_hops_embedded_mode_source(self) -> None:
        """Patche la copie locale de HOPS pour supporter un mode embarqué."""
        app_windows = self._hops_source_dir() / "hops" / "application_windows.py"
        app_file = self._hops_source_dir() / "hops" / "application.py"
        if not app_windows.is_file() or not app_file.is_file():
            raise FileNotFoundError("Fichiers HOPS introuvables pour le patch embarqué.")

        txt = app_windows.read_text(encoding="utf-8", errors="replace")
        if "embedded=False" not in txt:
            txt = txt.replace(
                "def __init__(self, log=HOPSLog, name='HOPS - window', sizex=None, sizey=None, position=5):",
                "def __init__(self, log=HOPSLog, name='HOPS - window', sizex=None, sizey=None, position=5, embedded=False, host_root=None, host_frame=None):",
            )
            txt = txt.replace(
                "        self.root = Tk()\n        self.hide()\n        self.main_frame = self.root\n        self.root.protocol('WM_DELETE_WINDOW', self.close)",
                "        self.embedded = bool(embedded)\n\n        if self.embedded:\n            if host_root is None or host_frame is None:\n                raise ValueError(\"Embedded mode requires host_root and host_frame\")\n            self.root = host_root\n            self.main_frame = host_frame\n        else:\n            self.root = Tk()\n            self.main_frame = self.root\n            self.root.protocol('WM_DELETE_WINDOW', self.close)\n            self.hide()",
            )
            txt = txt.replace(
                "        self.root.wm_title(name)\n        if sizex and sizey:\n            self.root.geometry('{0}x{1}'.format(int(self.root.winfo_screenwidth() / sizex),\n                                                int(self.root.winfo_screenheight() / sizey)))",
                "        if not self.embedded:\n            self.root.wm_title(name)\n            if sizex and sizey:\n                self.root.geometry('{0}x{1}'.format(int(self.root.winfo_screenwidth() / sizex),\n                                                    int(self.root.winfo_screenheight() / sizey)))",
            )
            txt = txt.replace(
                "class MainWindow(HOPSWindow):\n\n    def __init__(self, log, name, sizex=None, sizey=None, position=5):\n\n        HOPSWindow.__init__(self, log, name, sizex, sizey, position)",
                "class MainWindow(HOPSWindow):\n\n    def __init__(self, log, name, sizex=None, sizey=None, position=5, embedded=False, host_root=None, host_frame=None):\n\n        HOPSWindow.__init__(self, log, name, sizex, sizey, position, embedded=embedded, host_root=host_root, host_frame=host_frame)",
            )
            app_windows.write_text(txt, encoding="utf-8")
        # Correctif géométrie Tk: tous les widgets doivent cibler main_frame en mode embarqué.
        txt = app_windows.read_text(encoding="utf-8", errors="replace")
        txt = txt.replace("widget = Label(window.root, textvar=self.variable)", "widget = Label(window.main_frame, textvar=self.variable)")
        txt = txt.replace("widget = Frame(window.root)", "widget = Frame(window.main_frame)")
        app_windows.write_text(txt, encoding="utf-8")

        txt2 = app_file.read_text(encoding="utf-8", errors="replace")
        if "def __init__(self, embedded=False, host_root=None, host_frame=None):" not in txt2:
            txt2 = txt2.replace(
                "    def __init__(self):",
                "    def __init__(self, embedded=False, host_root=None, host_frame=None):",
            )
            txt2 = txt2.replace(
                "        MainWindow.__init__(self, HOPSLog(), name='HOPS-modified (A. Tsiaras, K. Karpouzas) - J.P Vignes <jeanpascal.vignes@gmail.com>', position=1)",
                "        MainWindow.__init__(self, HOPSLog(), name='HOPS-modified (A. Tsiaras, K. Karpouzas) - J.P Vignes <jeanpascal.vignes@gmail.com>', position=1, embedded=embedded, host_root=host_root, host_frame=host_frame)",
            )
            app_file.write_text(txt2, encoding="utf-8")

        txt2 = app_file.read_text(encoding="utf-8", errors="replace")
        if "NPOAP: embedded close callback" not in txt2:
            old_sig = (
                "    def __init__(self, embedded=False, host_root=None, host_frame=None):\n\n"
                "        MainWindow.__init__(self, HOPSLog(), name='HOPS-modified (A. Tsiaras, K. Karpouzas) - J.P Vignes <jeanpascal.vignes@gmail.com>', position=1, "
                "embedded=embedded, host_root=host_root, host_frame=host_frame)"
            )
            new_sig = (
                "    def __init__(self, embedded=False, host_root=None, host_frame=None, on_embedded_close=None):  "
                "# NPOAP: embedded close callback\n\n"
                "        self._on_embedded_close = on_embedded_close if embedded else None\n"
                "        MainWindow.__init__(self, HOPSLog(), name='HOPS-modified (A. Tsiaras, K. Karpouzas) - J.P Vignes <jeanpascal.vignes@gmail.com>', position=1, "
                "embedded=embedded, host_root=host_root, host_frame=host_frame)"
            )
            if old_sig in txt2:
                txt2 = txt2.replace(old_sig, new_sig, 1)
            old_block = "        self.update_window()\n\n    def open_updates(self):"
            new_block = (
                "        self.update_window()\n\n"
                "    def close(self):\n"
                "        MainWindow.close(self)\n"
                "        if self.embedded and getattr(self, '_on_embedded_close', None):\n"
                "            cb = self._on_embedded_close\n"
                "            self._on_embedded_close = None\n"
                "            try:\n"
                "                cb()\n"
                "            except Exception:\n"
                "                traceback.print_exc()\n\n"
                "    def open_updates(self):"
            )
            if old_block in txt2:
                txt2 = txt2.replace(old_block, new_block, 1)
            app_file.write_text(txt2, encoding="utf-8")

        self._patch_hops_embedded_frame_safe(app_windows)

    def _patch_hops_embedded_frame_safe(self, app_windows: Path) -> None:
        """Complète le patch embarqué : Frame n'a pas geometry/protocol/mainloop."""
        txt = app_windows.read_text(encoding="utf-8", errors="replace")
        if "NPOAP: embedded frame-safe" in txt:
            return
        if "host_frame" not in txt or "self.embedded" not in txt:
            return

        old = """    def reposition(self):

        self.root.update_idletasks()

        if self.position == 1:"""
        new = """    def reposition(self):

        self.root.update_idletasks()

        # NPOAP: embedded frame-safe (Frame host has no geometry/protocol)
        if self.embedded:
            try:
                self.main_frame.update_idletasks()
            except TclError:
                pass
            return

        if self.position == 1:"""
        if old not in txt:
            return
        txt = txt.replace(old, new, 1)

        old = """    def show(self):

        self.reposition()

        self.root.wm_attributes("-topmost", 1)"""
        new = """    def show(self):

        self.reposition()

        if self.embedded:
            self.update_idletasks()
            return

        # Garder les fenêtres HOPS auxiliaires visibles au-dessus de NPOAP.
        self.root.wm_attributes("-topmost", 1)
        self.root.deiconify()
        self.root.lift()
        try:
            self.root.focus_force()
        except TclError:
            pass
        self.update_idletasks()
        return"""
        txt = txt.replace(old, new, 1)

        old = """    def hide(self):

        self.root.withdraw()"""
        new = """    def hide(self):

        if self.embedded:
            return

        self.root.withdraw()"""
        txt = txt.replace(old, new, 1)

        old = """        self.after(f_after)

        self.root.mainloop()"""
        new = """        self.after(f_after)

        if self.embedded:
            return

        self.root.mainloop()"""
        txt = txt.replace(old, new, 1)

        old = """    def def_close(self):

        if self.mainloop_on:
            self.root.quit()

        for job in self.jobs:
            self.root.after_cancel(job)

        self.root.destroy()"""
        new = """    def def_close(self):

        if not self.embedded:
            if self.mainloop_on:
                self.root.quit()

        for job in self.jobs:
            try:
                self.root.after_cancel(job)
            except TclError:
                pass
        self.jobs.clear()

        if self.embedded:
            self.mainloop_on = False
            try:
                for w in list(self.main_frame.winfo_children()):
                    w.destroy()
            except TclError:
                pass
            return

        self.root.destroy()"""
        txt = txt.replace(old, new, 1)

        old = """    def disable(self):

        self.root.protocol('WM_DELETE_WINDOW', self.no_action)"""
        new = """    def disable(self):

        if not self.embedded:
            self.root.protocol('WM_DELETE_WINDOW', self.no_action)"""
        txt = txt.replace(old, new, 1)

        old = """    def activate(self):

        self.root.protocol('WM_DELETE_WINDOW', self.close)"""
        new = """    def activate(self):

        if not self.embedded:
            self.root.protocol('WM_DELETE_WINDOW', self.close)"""
        txt = txt.replace(old, new, 1)

        old = """    def set_close_button_function(self, function):
        self.root.protocol('WM_DELETE_WINDOW', function)"""
        new = """    def set_close_button_function(self, function):
        if self.embedded:
            return
        self.root.protocol('WM_DELETE_WINDOW', function)"""
        txt = txt.replace(old, new, 1)

        # Evite de contaminer le theme ttk global de NPOAP via theme_use().
        txt = txt.replace(
            "            combostyle = Style()\n"
            "            combostyle.theme_create('combostyle', parent='alt',\n"
            "                                    settings={'TCombobox': {'configure':\n"
            "                                                                {'selectbackground': 'white',\n"
            "                                                                 'fieldbackground': 'white',\n"
            "                                                                 'background': 'white'}}})\n"
            "            combostyle.theme_use('combostyle')",
            "            style_name = 'HOPS.TCombobox'\n"
            "            combostyle = Style()\n"
            "            combostyle.configure(style_name,\n"
            "                                selectbackground='white',\n"
            "                                fieldbackground='white',\n"
            "                                background='white')"
        )
        txt = txt.replace(
            "        if width:\n"
            "            widget = Combobox(window.main_frame, textvariable=self.variable, state='readonly', width=width)\n"
            "        else:\n"
            "            widget = Combobox(window.main_frame, textvariable=self.variable, state='readonly', width=int(window.log.entries_width * 0.8))",
            "        if width:\n"
            "            widget = Combobox(window.main_frame, textvariable=self.variable, state='readonly', width=width, style=style_name)\n"
            "        else:\n"
            "            widget = Combobox(window.main_frame, textvariable=self.variable, state='readonly', width=int(window.log.entries_width * 0.8), style=style_name)"
        )

        app_windows.write_text(txt, encoding="utf-8")

    def install_hops_from_zip(self) -> None:
        initial_zip = self._default_hops_zip()
        # Si l'archive de distribution est présente, on l'utilise directement.
        if initial_zip.is_file():
            path = str(initial_zip)
        else:
            path = filedialog.askopenfilename(
                title="Sélectionner HOPS-modified.zip / hops-master.zip",
                initialdir=str(initial_zip.parent if initial_zip.parent.exists() else Path.home()),
                filetypes=[("Archive ZIP", "*.zip"), ("Tous", "*.*")],
            )
            if not path:
                return
        try:
            zip_path = Path(path)
            if not zip_path.is_file():
                raise FileNotFoundError("Archive ZIP introuvable.")

            base = self._hops_base_dir()
            source = self._hops_source_dir()
            base.mkdir(parents=True, exist_ok=True)
            if source.exists():
                shutil.rmtree(source)
            # Extraction robuste : accepte soit un ZIP racine "hops-master/", soit un ZIP
            # contenant directement "hops/" et les scripts d'installation.
            with tempfile.TemporaryDirectory() as td:
                tmp_dir = Path(td)
                with zipfile.ZipFile(zip_path, "r") as zf:
                    zf.extractall(tmp_dir)

                # Cherche le dossier racine qui contient hops/__main__.py
                candidate_root = None
                if (tmp_dir / "hops" / "__main__.py").is_file():
                    candidate_root = tmp_dir
                else:
                    for p in tmp_dir.rglob("__main__.py"):
                        if p.parent.name == "hops":
                            candidate_root = p.parent.parent
                            break

                if candidate_root is None:
                    raise RuntimeError("Structure ZIP invalide : impossible de trouver hops/__main__.py")

                shutil.copytree(candidate_root, source)
            self._enable_hops_embedded_mode_source()

            if not self._hops_entry_file().is_file():
                raise RuntimeError("Structure HOPS inattendue : __main__.py introuvable.")

            self.hops_app = None
            self.hops_status_var.set(self._hops_status_text())
            messagebox.showinfo("HOPS", f"Installation terminée dans :\n{source}")
            logging.info(f"HOPS installé depuis {zip_path} vers {source}")
        except Exception as e:
            logging.error(f"Installation HOPS impossible: {e}", exc_info=True)
            messagebox.showerror("HOPS", f"Échec installation HOPS:\n{e}")

    def launch_hops(self) -> None:
        source = self._hops_source_dir()
        if not self._hops_entry_file().is_file():
            messagebox.showwarning(
                "HOPS",
                "HOPS n'est pas installé.\nUtilisez d'abord « Installer / Réinstaller HOPS (ZIP) ».",
            )
            return
        try:
            # Intégration en mode embarqué dans le cadre HOPS (pas de fenêtre externe).
            added_path = False
            if str(source) not in sys.path:
                sys.path.insert(0, str(source))
                added_path = True
            self._ensure_hops_dependencies()
            self._register_gaia_filters_from_resources()

            from hops.application import HOPS
            if added_path:
                try:
                    sys.path.remove(str(source))
                except ValueError:
                    pass

            self._recreate_hops_isolated_container()

            self.hops_app = HOPS(
                embedded=True,
                # Root dédié au conteneur pour éviter toute prise de contrôle du root principal.
                host_root=self.hops_isolated_container,
                host_frame=self.hops_isolated_container,
                on_embedded_close=self._on_hops_embedded_closed,
            )
            self.hops_app.run(f_after=self.hops_app.log.check_for_updates)
            self.hops_status_var.set(f"HOPS intégré dans le cadre depuis : {source}")
            logging.info(f"HOPS intégré dans l'onglet Exoplanètes: source={source}")
        except Exception as e:
            logging.error(f"Lancement HOPS impossible: {e}", exc_info=True)
            messagebox.showerror("HOPS", f"Impossible de lancer HOPS:\n{e}")

    def _on_hops_embedded_closed(self) -> None:
        """Appelé par HOPS après EXIT en mode embarqué (cadre vidé, référence à nettoyer)."""
        self.hops_app = None
        self.hops_status_var.set(
            "HOPS fermé (EXIT). Cliquez sur « Lancer HOPS » pour rouvrir."
        )

    def reset_hops(self) -> None:
        """Réinitialise HOPS embarqué et recharge la fenêtre principale dans le cadre."""
        try:
            self.hops_status_var.set("HOPS réinitialisé. Rechargement en cours...")
            self.launch_hops()
        except Exception as e:
            logging.error(f"Réinitialisation HOPS impossible: {e}", exc_info=True)
            messagebox.showerror("HOPS", f"Impossible de réinitialiser HOPS:\n{e}")
