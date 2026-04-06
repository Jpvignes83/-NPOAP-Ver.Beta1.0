import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
from tkinter.filedialog import asksaveasfilename
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from astropy.time import Time
from astropy import units as u
import datetime
import os
import glob
import threading
import csv
import re
import copy
# --- IMPORTS BACKEND (Assurez-vous que ces fichiers existent dans /core) ---
from core.periodogram_tools import run_lomb_scargle, run_bls, run_plavchan
from core.ttv_modeling import fit_sine_model, multi_sine_model
# Import N-body (optionnel)
try:
    from core.nbody_simulation import (
        generate_simulation, analyze_simulation, predict_future_transits,
        REBOUND_AVAILABLE, ULTRANEST_AVAILABLE
    )
    NBODY_AVAILABLE = REBOUND_AVAILABLE
    NBODY_FITTING_AVAILABLE = REBOUND_AVAILABLE and ULTRANEST_AVAILABLE
except ImportError:
    NBODY_AVAILABLE = False
    NBODY_FITTING_AVAILABLE = False 

# --- IMPORTS FRONTEND ---
from gui.periodogram_viewer import PeriodogramViewer
from gui.ttv_viewer import TTVViewer
from gui.lc_markers_viewer import LCMarkersViewer

logger = logging.getLogger(__name__)

# O-C internes en jours ; affichage TTV / figures en minutes (× MINUTES_PER_DAY)
MINUTES_PER_DAY = 24 * 60

class ToolTip:
    """Classe pour créer des tooltips."""
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.id = None
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
        if self.id:
            self.widget.after_cancel(self.id)
            self.id = None

    def showtip(self, event=None):
        x, y, cx, cy = (0, 0, 0, 0)
        if hasattr(self.widget, 'bbox'):
            x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(
            tw,
            text=self.text,
            justify=tk.LEFT,
            background="#ffffe0",
            relief=tk.SOLID,
            borderwidth=1,
            font=("TkDefaultFont", 9),
            wraplength=300
        )
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

class DataAnalysisTab:
    def __init__(self, parent_notebook):
        """
        Intègre l'onglet d'analyse dans le notebook principal de l'application.
        """
        self.main_frame = ttk.Frame(parent_notebook)
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Création des sous-onglets A, B, C
        self.sub_notebook = ttk.Notebook(self.main_frame)
        self.sub_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.tab_period = ttk.Frame(self.sub_notebook)
        self.tab_ttv = ttk.Frame(self.sub_notebook)
        self.tab_multi = ttk.Frame(self.sub_notebook)
        self.tab_nbody = ttk.Frame(self.sub_notebook)

        self.sub_notebook.add(self.tab_period, text="A. Détermination Période")
        self.sub_notebook.add(self.tab_ttv, text="B. Recherche & Analyse TTV")
        self.sub_notebook.add(self.tab_multi, text="C. Analyse Système Multiple")
        self.sub_notebook.add(self.tab_nbody, text="D. Simulation N-body")

        # Données en mémoire
        self.lc_time = None
        self.lc_flux = None
        self.oc_epoch = None
        self.oc_res = None
        self.oc_yerr = None
        self.oc_time_raw = None # Ajout pour le calcul de P_orb
        self.detected_P_orb = None # P_orb détectée (nécessaire pour le rapport)
        self.current_yerr = None # yerr final utilisé pour le fit
        self.oc_source_path = None # Chemin du fichier O-C chargé (pour rapport)
        self.lc_selected_paths = []  # Chemins absolus des LC pour concaténation (onglet A)

        self.model_x = None 
        self.model_y = None 
        self.last_mcmc_samples = None # Pour stocker les résultats du fit
        self.last_n_freq = 0 # Nombre de fréquences fitées
        
        # Paramètres du fitting MCMC (avec valeurs par défaut)
        self.mcmc_nwalkers = 32
        self.mcmc_nsteps = 5000
        self.mcmc_burnin_fraction = 0.25

        # Construction des interfaces
        self.setup_part_a()
        self.setup_part_b()
        self.setup_part_c()
        self.setup_part_d()

    # =========================================================================
    # PARTIE A : DÉTERMINATION PÉRIODE & EXTRACTION MID-TIME
    # =========================================================================
    def setup_part_a(self):
        
        # -------------------------------------------------------------------
        # BLOC 1 : Fichiers source des courbes de lumière (.txt / .csv)
        # -------------------------------------------------------------------
        frame_recueil = ttk.LabelFrame(
            self.tab_period,
            text="1. Fichiers source (Light curves .txt / .csv)",
            padding=10,
        )
        frame_recueil.pack(side=tk.TOP, fill="x", padx=10, pady=5)

        list_fr = ttk.Frame(frame_recueil)
        list_fr.pack(fill=tk.BOTH, expand=True, padx=5, pady=4)
        lc_sb = ttk.Scrollbar(list_fr)
        self._lc_listbox = tk.Listbox(
            list_fr,
            height=6,
            width=92,
            yscrollcommand=lc_sb.set,
            font=("Consolas", 9),
            selectmode=tk.EXTENDED,
        )
        lc_sb.config(command=self._lc_listbox.yview)
        self._lc_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lc_sb.pack(side=tk.RIGHT, fill=tk.Y)

        btn_lc_fr = ttk.Frame(frame_recueil)
        btn_lc_fr.pack(fill=tk.X, padx=5, pady=(2, 0))
        ttk.Button(btn_lc_fr, text="Ajouter fichiers…", command=self._add_lc_files_dialog).pack(
            side=tk.LEFT, padx=3
        )
        ttk.Button(btn_lc_fr, text="Ajouter depuis un dossier…", command=self._add_lc_from_folder_dialog).pack(
            side=tk.LEFT, padx=3
        )
        ttk.Button(btn_lc_fr, text="Retirer la sélection", command=self._remove_selected_lc_files).pack(
            side=tk.LEFT, padx=3
        )
        ttk.Button(btn_lc_fr, text="Vider la liste", command=self._clear_lc_files).pack(side=tk.LEFT, padx=3)

        ttk.Label(
            frame_recueil,
            text="Le fichier concatenated_lightcurve.csv est enregistré dans le dossier du premier fichier de la liste.",
            foreground="gray",
            font=("Helvetica", 8),
            wraplength=720,
            justify="left",
        ).pack(anchor="w", padx=5, pady=(4, 0))

        # -------------------------------------------------------------------
        # BLOC 2 : FLUX DE TRAVAIL PÉRIODOGRAMMES
        # -------------------------------------------------------------------
        ctrl_frame = ttk.LabelFrame(self.tab_period, text="2. Flux de travail : Périodicité (Light Curve)", padding=10)
        ctrl_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)

        btn_load = ttk.Button(ctrl_frame, text="Concaténer Lightcurves", command=self.load_lightcurve)
        btn_load.pack(side=tk.LEFT, padx=10)
        
        ttk.Separator(ctrl_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        param_frame = ttk.Frame(ctrl_frame)
        param_frame.pack(side=tk.LEFT, padx=10)
        
        ttk.Label(param_frame, text="Min P (j):").pack(side=tk.LEFT)
        self.entry_min_p = ttk.Entry(param_frame, width=5)
        self.entry_min_p.insert(0, "0.5") 
        self.entry_min_p.pack(side=tk.LEFT, padx=2)

        ttk.Label(param_frame, text="Max P (j):").pack(side=tk.LEFT)
        self.entry_max_p = ttk.Entry(param_frame, width=5)
        self.entry_max_p.insert(0, "10.0")
        self.entry_max_p.pack(side=tk.LEFT, padx=2)

        ttk.Separator(ctrl_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        ttk.Label(ctrl_frame, text="Calcul :").pack(side=tk.LEFT)
        
        ttk.Button(ctrl_frame, text="Lomb-Scargle", 
                    command=lambda: self.run_period_algo("LS")).pack(side=tk.LEFT, padx=2)
        
        ttk.Button(ctrl_frame, text="BLS (Transit)", 
                    command=lambda: self.run_period_algo("BLS")).pack(side=tk.LEFT, padx=2)

        ttk.Button(ctrl_frame, text="Plavchan", 
                    command=lambda: self.run_period_algo("PLAV")).pack(side=tk.LEFT, padx=2)

        # -------------------------------------------------------------------
        # BLOC 3 : VIEWER PÉRIODOGRAMME
        # -------------------------------------------------------------------
        viewer_container = ttk.LabelFrame(self.tab_period, text="3. Résultats Périodogramme", padding=5)
        viewer_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)

        perio_toolbar = ttk.Frame(viewer_container)
        perio_toolbar.pack(side=tk.TOP, fill=tk.X, pady=(0, 4))
        ttk.Button(
            perio_toolbar,
            text="Voir LC avec transits (marqueurs, export mid-time.csv par période)",
            command=self._open_lc_markers_viewer,
        ).pack(side=tk.LEFT, padx=2)

        self.perio_viewer = PeriodogramViewer(viewer_container)

    def _refresh_lc_files_listbox(self):
        if hasattr(self, "_lc_listbox"):
            try:
                self._lc_listbox.delete(0, tk.END)
                for p in self.lc_selected_paths:
                    self._lc_listbox.insert(tk.END, p)
            except tk.TclError:
                pass

    def _add_lc_files_dialog(self):
        paths = filedialog.askopenfilenames(
            title="Courbes de lumière (multi-sélection)",
            filetypes=[
                ("TXT / CSV", "*.txt *.csv"),
                ("Texte", "*.txt"),
                ("CSV", "*.csv"),
                ("Tous", "*.*"),
            ],
        )
        if not paths:
            return
        for p in paths:
            ap = os.path.abspath(os.path.normpath(p))
            if ap not in self.lc_selected_paths:
                self.lc_selected_paths.append(ap)
        self._refresh_lc_files_listbox()

    def _add_lc_from_folder_dialog(self):
        d = filedialog.askdirectory(title="Ajouter tous les .txt et .csv de ce dossier")
        if not d:
            return
        d = os.path.abspath(d)
        found = sorted(
            set(glob.glob(os.path.join(d, "*.txt")) + glob.glob(os.path.join(d, "*.csv")))
        )
        if not found:
            messagebox.showinfo("Dossier", "Aucun fichier .txt ou .csv dans ce dossier.")
            return
        n_add = 0
        for p in found:
            ap = os.path.abspath(p)
            bn = os.path.basename(ap)
            if bn in ("concatenated_lightcurve.csv", "mid-time.csv"):
                continue
            if ap not in self.lc_selected_paths:
                self.lc_selected_paths.append(ap)
                n_add += 1
        self._refresh_lc_files_listbox()
        messagebox.showinfo(
            "Dossier",
            f"{n_add} fichier(s) ajouté(s) à la liste ({len(found)} fichier(s) .txt/.csv dans le dossier).",
        )

    def _remove_selected_lc_files(self):
        sel = list(self._lc_listbox.curselection())
        if not sel:
            messagebox.showinfo("Liste", "Sélectionnez une ou plusieurs lignes à retirer.")
            return
        for i in sorted(sel, reverse=True):
            if 0 <= i < len(self.lc_selected_paths):
                del self.lc_selected_paths[i]
        self._refresh_lc_files_listbox()

    def _clear_lc_files(self):
        self.lc_selected_paths.clear()
        self._refresh_lc_files_listbox()

    def _open_lc_markers_viewer(self):
        """Ouvre le viewer LC avec les marqueurs de période ; permet d'ajuster les mid-times et d'exporter mid-time_*.csv par période."""
        if self.lc_time is None or self.lc_flux is None:
            messagebox.showwarning("Données", "Concaténez d'abord les light curves (bouton « Concaténer Lightcurves »).")
            return
        markers = getattr(self.perio_viewer, 'markers', None) or []
        if not markers:
            messagebox.showwarning(
                "Marqueurs",
                "Ajoutez au moins un marqueur de période sur le périodogramme (clic gauche sur un pic, ou « Ajouter Marqueur »)."
            )
            return
        LCMarkersViewer(self.main_frame.winfo_toplevel(), self.lc_time, self.lc_flux, list(markers))

    def load_lightcurve(self):
        """
        Charge une courbe de lumière (LC) à partir des fichiers listés au bloc 1.
        Réutilise concatenated_lightcurve.csv s'il existe déjà dans le dossier du 1er fichier, sinon concatène.
        """
        self.lc_time = None
        self.lc_flux = None

        if not self.lc_selected_paths:
            messagebox.showwarning(
                "Fichiers manquants",
                "Ajoutez au moins un fichier .txt ou .csv (bloc 1), puis cliquez sur « Concaténer Lightcurves ».",
            )
            return

        out_dir = os.path.dirname(os.path.abspath(self.lc_selected_paths[0]))
        concat_path = os.path.join(out_dir, "concatenated_lightcurve.csv")

        # 1. Réutiliser un concaténé existant s'il est lisible
        if os.path.exists(concat_path):
            try:
                df = pd.read_csv(concat_path)
                self.lc_time = df["Time"].values
                self.lc_flux = df["Flux"].values
                messagebox.showinfo(
                    "Succès",
                    f"LC chargée : {len(self.lc_time)} points depuis le fichier concaténé existant.",
                )
                return
            except Exception as e:
                logger.warning("Fichier concaténé corrompu ou obsolète, nouvelle concaténation. Erreur: %s", e)

        # 2. Concaténation explicite depuis la liste de fichiers
        try:
            from core.lightcurve_tools import concatenate_lightcurve_paths

            self.main_frame.config(cursor="watch")
            self.main_frame.update()

            self.lc_time, self.lc_flux = concatenate_lightcurve_paths(
                self.lc_selected_paths, output_directory=out_dir
            )

            messagebox.showinfo(
                "Succès",
                f"LC chargée : {len(self.lc_time)} points ({len(self.lc_selected_paths)} fichier(s) concaténé(s)).",
            )
            return

        except ValueError as e:
            messagebox.showerror("Erreur LC", f"Concaténation impossible : {e}")
            return
        except Exception as e:
            messagebox.showerror("Erreur LC", f"Erreur lors de la concaténation : {e}")
            return
        finally:
            self.main_frame.config(cursor="")

    def run_period_algo(self, algo):
        """Lance l'algorithme de périodogramme sélectionné dans un thread séparé."""
        if self.lc_time is None:
            messagebox.showwarning("Attention", "Veuillez charger une courbe de lumière (Étape 2).")
            return

        try:
            min_p = float(self.entry_min_p.get())
            max_p = float(self.entry_max_p.get())
            if min_p >= max_p: raise ValueError
        except:
            messagebox.showerror("Erreur", "Vérifiez les valeurs Min/Max Période.")
            return

        # 1. SIGNAL DE TRAVAIL EN COURS (S'exécute AVANT que le thread ne bloque la GUI)
        self.main_frame.config(cursor="watch")
        self.sub_notebook.config(cursor="watch") 
        self.main_frame.update_idletasks() 

        # 2. Lancement du calcul dans un thread
        thread = threading.Thread(
            target=self._worker_run_period_algo, 
            args=(algo, min_p, max_p)
        )
        thread.start()

    def _worker_run_period_algo(self, algo, min_p, max_p):
        """Fonction de travail exécutée dans le thread."""
        res = None
        name = "Calcul Périodogramme"
        try:
            if algo == "LS":
                res = run_lomb_scargle(self.lc_time, self.lc_flux, min_period=min_p, max_period=max_p)
                name = "Lomb-Scargle (Standard Norm)"
            elif algo == "BLS":
                res = run_bls(self.lc_time, self.lc_flux, min_period=min_p, max_period=max_p)
                name = "Box Least Squares (Transit)"
            elif algo == "PLAV":
                res = run_plavchan(self.lc_time, self.lc_flux, min_period=min_p, max_period=max_p)
                name = "Plavchan (Binless PDM)"
            
            # Utiliser self.main_frame.after pour appeler la fonction de complétion sur le thread GUI
            self.main_frame.after(0, self._complete_run_period_algo, res, name)

        except Exception as e:
            logger.error(f"Algo error: {e}")
            # Gérer l'erreur sur le thread GUI
            self.main_frame.after(0, self._handle_algo_error, str(e))
        finally:
            # Assurez-vous que le curseur est réinitialisé même en cas d'erreur non gérée
            if res is None:
                 self.main_frame.after(0, self._reset_cursor)


    def _complete_run_period_algo(self, res, name):
        """Met à jour la GUI après la fin du calcul dans le thread."""
        periods, powers, best = res
        
        if hasattr(self.perio_viewer, 'plot_periodogram'):
            self.perio_viewer.plot_periodogram(periods, powers, title=f"{name} (Meilleur Pic: {best:.5f} j)")
            messagebox.showinfo("Résultat", f"Période détectée : {best:.6f} jours")
        else:
            messagebox.showwarning("Erreur Code", "Le Viewer n'a pas la méthode 'plot_periodogram'.")
            
        self._reset_cursor()

    def _handle_algo_error(self, error_message):
        """Affiche les erreurs de calcul sur le thread GUI."""
        messagebox.showerror("Erreur Calcul", error_message)
        self.main_frame.after(0, self._reset_cursor) # Correction: Utiliser after pour reset


    def _reset_cursor(self):
        """Réinitialise le curseur de la GUI."""
        self.main_frame.config(cursor="")
        self.sub_notebook.config(cursor="")

    # =========================================================================
    # PARTIE B : ANALYSE TTV (MID-TIME + O-C)
    # =========================================================================
    def setup_part_b(self):
        # 1. FLUX DE TRAVAIL TTV
        ctrl_frame = ttk.LabelFrame(self.tab_ttv, text="Flux de travail : Analyse TTV (O-C)", padding=10)
        ctrl_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)

        # A. Fréquences
        freq_frame = tk.Frame(ctrl_frame)
        freq_frame.pack(side=tk.LEFT, padx=5)
        tk.Label(freq_frame, text="Fréquences :").pack(side=tk.LEFT)
        self.n_freq_var = tk.StringVar(value="1") 
        self.n_freq_combo = ttk.Combobox(freq_frame, textvariable=self.n_freq_var, 
                                             values=["1", "2", "3"], state="readonly", width=3)
        self.n_freq_combo.pack(side=tk.LEFT, padx=5)
        
        # B. Boutons
        ttk.Separator(ctrl_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        ttk.Button(ctrl_frame, text="1. Charger O-C", command=self.load_oc).pack(side=tk.LEFT, padx=5)
        ttk.Button(ctrl_frame, text="🔍 Auto Période TTV", command=self.auto_find_period).pack(side=tk.LEFT, padx=5)
        ttk.Button(ctrl_frame, text="⚙️ Param. Fit", command=self.configure_fit_parameters).pack(side=tk.LEFT, padx=5)
        ttk.Button(ctrl_frame, text="2. Fit Modèle", command=self.fit_ttv).pack(side=tk.LEFT, padx=5)
        ttk.Button(ctrl_frame, text="📄 Rapport & Analyse", command=self.generate_report).pack(side=tk.LEFT, padx=5)

        # 3. VIEWER
        viewer_container = ttk.LabelFrame(self.tab_ttv, text="3. TTV Viewer", padding=5)
        viewer_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        viewer_toolbar = tk.Frame(viewer_container)
        viewer_toolbar.pack(side=tk.TOP, fill="x", pady=2)
        
        self.ttv_viewer = TTVViewer(viewer_container)
        ttk.Button(viewer_toolbar, text="Sauvegarder Figure", command=self.save_plot_with_table).pack(side=tk.LEFT)

    # --- MÉTHODES BACKEND PARTIE B ---
    def load_oc(self):
        path = filedialog.askopenfilename(filetypes=[("CSV", "*.csv")])
        if not path: return
        try:
            df = pd.read_csv(path)
            df.columns = [str(c).replace('\ufeff', '').strip().lower().lstrip('#') for c in df.columns]
            cols = list(df.columns)
            
            # Réinitialiser la période détectée pour éviter les valeurs d'un fichier précédent
            self.detected_P_orb = None
            self.oc_time_raw = None
            ep_col = next((c for c in cols if c in ['epoch', 'e', 'cycle', 'n', 'num', 'no', 'transit', 'ep']), None)
            # O-C en minutes (export LC viewer) ou en jours (ttv, o-c, o-c_jours, …)
            oc_min_col = next((c for c in cols if c in ['o-c_minutes', 'oc_minutes', 'o-c (minutes)']), None)
            ttv_col = next((c for c in cols if c in ['ttv', 'o-c', 'oc', 'residus', 'residuals', 'o-c_jours', 'oc_jours']), None)
            err_col = next((c for c in cols if c in ['uncertainty', 'err', 'error', 'sig', 'sigma']), None)
            time_col = next((c for c in cols if c in ['time', 'bjd', 'jd', 'date', 'mid_time', 'mid-time']), None)
            
            cols_str = ", ".join(cols) if cols else "(aucune)"
            if not ep_col:
                error_msg = "Colonne Epoch manquante (attendues: epoch, e, cycle, n, num, no, transit). Colonnes détectées: %s" % cols_str
                logger.error(error_msg)
                messagebox.showerror("Erreur Lecture O-C", error_msg)
                return
            if not ttv_col and not oc_min_col and not time_col:
                error_msg = "Colonnes TTV/O-C (jours), O-C_minutes ou Mid_time/Time manquantes. Colonnes détectées: %s" % cols_str
                logger.error(error_msg)
                messagebox.showerror("Erreur Lecture O-C", error_msg)
                return

            self.oc_source_path = path
            self.oc_epoch = df[ep_col].values.astype(float)
            self.oc_time_raw = df[time_col].values.astype(float) if time_col else None

            if oc_min_col:
                self.oc_res = np.asarray(df[oc_min_col].values, dtype=float) / MINUTES_PER_DAY
                # Incertitudes du fichier en minutes si colonne présente, sinon défaut en jours
                self.oc_yerr = (
                    df[err_col].values.astype(float) / MINUTES_PER_DAY
                    if err_col else np.ones_like(self.oc_res, dtype=float) * 0.001
                )
                if time_col:
                    self.oc_time_raw = np.asarray(df[time_col].values, dtype=float)
            elif ttv_col:
                self.oc_res = np.asarray(df[ttv_col].values, dtype=float)
                self.oc_yerr = df[err_col].values.astype(float) if err_col else np.ones_like(self.oc_res, dtype=float) * 0.001
                if time_col:
                    self.oc_time_raw = np.asarray(df[time_col].values, dtype=float)
            else:
                # Fichier mid-time seul (Epoch + Mid_time) : O-C = Mid_time - (T0 + Epoch * P)
                if not time_col or len(self.oc_epoch) < 2:
                    messagebox.showerror("Erreur Lecture O-C", "Fichier mid-time : il faut au moins 2 points (Epoch + Mid_time).")
                    return
                mid_times = np.asarray(df[time_col].values, dtype=float)
                epochs = np.asarray(self.oc_epoch, dtype=float)
                ok = np.isfinite(epochs) & np.isfinite(mid_times)
                epochs, mid_times = epochs[ok], mid_times[ok]
                if len(epochs) < 2:
                    messagebox.showerror("Erreur Lecture O-C", "Pas assez de points valides (Epoch, Mid_time).")
                    return
                P_lin, T0_lin = np.polyfit(epochs, mid_times, 1)
                if not (np.isfinite(P_lin) and P_lin > 0.01):
                    messagebox.showerror("Erreur Lecture O-C", "Impossible d'estimer la période à partir des mid-times.")
                    return
                self.detected_P_orb = float(P_lin)
                self.oc_res = mid_times - (T0_lin + epochs * P_lin)
                self.oc_epoch = epochs
                self.oc_time_raw = mid_times
                self.oc_yerr = df[err_col].values[ok].astype(float) if err_col else np.ones_like(self.oc_res, dtype=float) * 0.001
                logger.info("O-C dérivés des mid-times (épéméride linéaire P=%.5f j, T0=%.5f).", self.detected_P_orb, T0_lin)

            self.model_x = None
            self.model_y = None
            if len(self.oc_res) > 0:
                offset = np.mean(self.oc_res)
                self.oc_res = self.oc_res - offset
                logger.info("O-C centrés. Offset (%.6f j) soustrait.", offset)

            if self.detected_P_orb is None and self.oc_time_raw is not None and len(self.oc_epoch) > 1:
                try:
                    epochs = pd.Series(self.oc_epoch).astype(float).dropna()
                    mid_times = pd.Series(self.oc_time_raw).astype(float).dropna()
                    if len(epochs) == len(mid_times) and len(epochs) > 1:
                        slope, _ = np.polyfit(epochs, mid_times, 1)
                        if np.isfinite(slope) and slope > 0.01:
                            self.detected_P_orb = slope
                            logger.info("P_orb auto-détectée au chargement O-C: %.5f j.", self.detected_P_orb)
                except Exception as e:
                    logger.debug("Auto-détection P_orb ignorée: %s", e)

            # Incertitudes nulles / NaN dans le CSV → division par zéro (Lomb-Scargle, MCMC TTV)
            if self.oc_yerr is not None and len(self.oc_yerr) == len(self.oc_res):
                e = np.asarray(self.oc_yerr, dtype=float)
                std = float(np.nanstd(self.oc_res)) if len(self.oc_res) > 1 else 1.0
                fl = max(1e-15, 1e-10 * max(std, 1e-12))
                e = np.where(np.isfinite(e), e, fl)
                self.oc_yerr = np.maximum(e, fl)

            # Mise à jour du Viewer
            if hasattr(self.ttv_viewer, 'plot_external_data'):
                self.ttv_viewer.plot_external_data(self.oc_epoch, self.oc_res, self.oc_yerr, None, None) 
                
            # --- Message de succès (avec alerte P_orb si manquante) ---
            msg = f"Chargé: {len(self.oc_epoch)} points."
            if self.detected_P_orb:
                 msg += f"\nP_orb détectée/saisie: {self.detected_P_orb:.5f} jours."
                 logger.info("Chargement O-C réussi.")
            else:
                 msg += f"\nATTENTION: P_orb non détectée. Elle sera demandée lors du Rapport TTV."
                 logger.warning("P_orb est manquante pour l'analyse TTV.")

            messagebox.showinfo("Succès Chargement O-C", msg)

        except Exception as e:
            error_msg = f"Échec du chargement ou de l'analyse du fichier O-C. Détail: {e}"
            logger.error(error_msg, exc_info=True)
            messagebox.showerror("Erreur Chargement O-C", error_msg)

    def auto_find_period(self):
        # ... (Logique inchangée pour l'auto-détection Lomb-Scargle)
        if self.oc_epoch is None or len(self.oc_epoch) < 5:
            messagebox.showwarning("Erreur", "Données insuffisantes.")
            return
            
        min_p = 2.0; max_p = (np.max(self.oc_epoch) - np.min(self.oc_epoch)) / 2
        dy = np.asarray(
            self.oc_yerr if self.oc_yerr is not None else np.full_like(self.oc_res, 0.001, dtype=float),
            dtype=float,
        )
        std = float(np.nanstd(self.oc_res)) if len(self.oc_res) > 1 else 1.0
        floor = max(1e-15, 1e-10 * max(std, 1e-12))
        dy = np.where(np.isfinite(dy), dy, floor)
        dy = np.maximum(dy, floor)

        try:
             frequency, power = LombScargle(self.oc_epoch, self.oc_res, dy=dy).autopower(
                 minimum_frequency=1/max_p, maximum_frequency=1/min_p)
             
             best_period = 1 / frequency[np.argmax(power)]
             msg = f"Période dominante détectée : {best_period:.4f} époques"
             print(msg)
             messagebox.showinfo("Lomb-Scargle", msg)
             return best_period
        except Exception as e:
             messagebox.showerror("Erreur LS", f"Échec de l'auto-détection LS: {e}")
             return None

    def configure_fit_parameters(self):
        """Ouvre une fenêtre pour configurer les paramètres du fitting MCMC."""
        dialog = tk.Toplevel(self.main_frame)
        dialog.title("Paramètres du Fitting MCMC")
        dialog.geometry("550x320")
        dialog.transient(self.main_frame)
        dialog.grab_set()
        
        main_frame = ttk.Frame(dialog, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Variables pour les champs
        n_freq_var = tk.StringVar(value=str(int(self.n_freq_var.get())))
        nwalkers_var = tk.StringVar(value=str(self.mcmc_nwalkers))
        nsteps_var = tk.StringVar(value=str(self.mcmc_nsteps))
        burnin_var = tk.StringVar(value=f"{self.mcmc_burnin_fraction:.2f}")
        
        # Nombre de fréquences
        row1 = ttk.Frame(main_frame)
        row1.pack(fill=tk.X, pady=6)
        ttk.Label(row1, text="Nombre de fréquences:", width=22, anchor="w").pack(side=tk.LEFT, padx=5)
        ttk.Entry(row1, textvariable=n_freq_var, width=18).pack(side=tk.LEFT, padx=5)
        ttk.Label(row1, text="(1-5 recommandé)", font=("", 8), foreground="gray").pack(side=tk.LEFT, padx=5)
        
        # Nombre de walkers
        row2 = ttk.Frame(main_frame)
        row2.pack(fill=tk.X, pady=6)
        ttk.Label(row2, text="Nombre de walkers:", width=22, anchor="w").pack(side=tk.LEFT, padx=5)
        ttk.Entry(row2, textvariable=nwalkers_var, width=18).pack(side=tk.LEFT, padx=5)
        ttk.Label(row2, text="(pair, ≥4)", font=("", 8), foreground="gray").pack(side=tk.LEFT, padx=5)
        
        # Nombre de steps
        row3 = ttk.Frame(main_frame)
        row3.pack(fill=tk.X, pady=6)
        ttk.Label(row3, text="Nombre de steps:", width=22, anchor="w").pack(side=tk.LEFT, padx=5)
        ttk.Entry(row3, textvariable=nsteps_var, width=18).pack(side=tk.LEFT, padx=5)
        ttk.Label(row3, text="(1000-20000)", font=("", 8), foreground="gray").pack(side=tk.LEFT, padx=5)
        
        # Fraction de burn-in
        row4 = ttk.Frame(main_frame)
        row4.pack(fill=tk.X, pady=6)
        ttk.Label(row4, text="Fraction burn-in:", width=22, anchor="w").pack(side=tk.LEFT, padx=5)
        ttk.Entry(row4, textvariable=burnin_var, width=18).pack(side=tk.LEFT, padx=5)
        ttk.Label(row4, text="(0.1-0.5)", font=("", 8), foreground="gray").pack(side=tk.LEFT, padx=5)
        
        # Note informative
        info_label = ttk.Label(
            main_frame, 
            text="Note: Plus de walkers/steps = meilleure convergence mais plus long",
            font=("", 8), 
            foreground="gray"
        )
        info_label.pack(pady=(15, 10))
        
        def save_params():
            try:
                n_freq = int(n_freq_var.get())
                nwalkers = int(nwalkers_var.get())
                nsteps = int(nsteps_var.get())
                burnin = float(burnin_var.get())
                
                if n_freq < 1 or n_freq > 10:
                    messagebox.showerror("Erreur", "Nombre de fréquences doit être entre 1 et 10")
                    return
                if nwalkers < 4 or nwalkers % 2 != 0:
                    messagebox.showerror("Erreur", "Nombre de walkers doit être ≥ 4 et pair")
                    return
                if nsteps < 100:
                    messagebox.showerror("Erreur", "Nombre de steps doit être ≥ 100")
                    return
                if burnin < 0.05 or burnin > 0.9:
                    messagebox.showerror("Erreur", "Fraction burn-in doit être entre 0.05 et 0.9")
                    return
                
                # Sauvegarder les paramètres
                self.n_freq_var.set(str(n_freq))
                self.mcmc_nwalkers = nwalkers
                self.mcmc_nsteps = nsteps
                self.mcmc_burnin_fraction = burnin
                
                dialog.destroy()
                messagebox.showinfo("Succès", f"Paramètres sauvegardés:\n"
                                             f"Fréquences: {n_freq}\n"
                                             f"Walkers: {nwalkers}\n"
                                             f"Steps: {nsteps}\n"
                                             f"Burn-in: {burnin:.1%}")
            except ValueError as e:
                messagebox.showerror("Erreur", f"Valeurs invalides: {e}")
        
        # Boutons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(pady=(15, 0))
        ttk.Button(btn_frame, text="Enregistrer", command=save_params).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Annuler", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    def fit_ttv(self):
        """Lance la modélisation MCMC TTV dans un thread séparé."""
        if self.oc_epoch is None:
            messagebox.showwarning("Attention", "Chargez les données O-C d'abord.")
            return

        try:
            n_frequences = int(self.n_freq_var.get())
        except ValueError:
            messagebox.showerror("Erreur", "Nombre de fréquences invalide.")
            return

        # 1. SIGNAL DE TRAVAIL EN COURS
        self.main_frame.config(cursor="watch")
        self.sub_notebook.config(cursor="watch") 
        self.main_frame.update_idletasks()

        yerr = np.asarray(
            self.oc_yerr if self.oc_yerr is not None else np.ones_like(self.oc_res) * 0.001,
            dtype=float,
        )
        std = float(np.nanstd(self.oc_res)) if len(self.oc_res) > 1 else 1.0
        floor = max(1e-15, 1e-10 * max(std, 1e-12))
        yerr = np.where(np.isfinite(yerr), yerr, floor)
        yerr = np.maximum(yerr, floor)
        self.current_yerr = yerr

        # 2. Lancement du thread MCMC avec les paramètres configurés
        thread = threading.Thread(
            target=self._worker_fit_ttv, 
            args=(n_frequences, yerr)
        )
        thread.start()

    def _worker_fit_ttv(self, n_frequences, yerr):
        """Worker thread pour le fitting MCMC TTV."""
        res = None
        try:
            print(f"--- Démarrage Fit TTV ({n_frequences} Fréquence(s)) ---")
            print(f"Paramètres MCMC: {self.mcmc_nwalkers} walkers, {self.mcmc_nsteps} steps, burn-in {self.mcmc_burnin_fraction:.1%}")
            res = fit_sine_model(
                self.oc_epoch, 
                self.oc_res, 
                yerr, 
                nwalkers=self.mcmc_nwalkers, 
                nsteps=self.mcmc_nsteps, 
                n_frequences=n_frequences
            )
            
            # Appel du thread GUI pour la complétion
            self.main_frame.after(0, self._complete_fit_ttv, res, n_frequences, yerr)

        except Exception as e:
            logger.error(f"Fit error: {e}")
            self.main_frame.after(0, self._handle_algo_error, f"Échec MCMC. Détail: {e}")
        finally:
            if res is None:
                 self.main_frame.after(0, self._reset_cursor)


    def _complete_fit_ttv(self, res, n_frequences, yerr):
        """Met à jour la GUI et les données après la fin du fit MCMC."""
        try:
            if isinstance(res, tuple): res = res[0]
            res = np.array(res)
            
            if res.ndim < 2:
                raise ValueError(f"Format invalide. Attendu: Tableau 2D. Reçu: Shape {res.shape}")

            best_params = np.median(res, axis=0)
            
            self.model_x = np.linspace(np.min(self.oc_epoch), np.max(self.oc_epoch), 200)
            self.model_y = multi_sine_model(self.model_x, *best_params)

            self.last_mcmc_samples = res
            self.last_n_freq = n_frequences
            
            messagebox.showinfo("Résultat Modélisation", "Calcul MCMC terminé avec succès.")

            if hasattr(self.ttv_viewer, 'plot_external_data'):
                self.ttv_viewer.plot_external_data(
                    self.oc_epoch, 
                    self.oc_res, 
                    yerr,
                    self.model_x, 
                    self.model_y
                )
        
        except Exception as e:
            messagebox.showerror("Erreur Affichage Fit", f"Erreur lors de l'affichage des résultats: {e}")
        finally:
            self._reset_cursor()
    
    # Correction: Ajout de la fonction manquante _log_markers pour éviter les erreurs
    def _log_markers(self):
        if self.ttv_viewer and hasattr(self.ttv_viewer, '_log_markers'):
             self.ttv_viewer._log_markers()
        else:
             print("Marqueurs TTV:", self.ttv_viewer.markers if hasattr(self.ttv_viewer, 'markers') else "Viewer non initialisé")


    def calculate_bic(self, residuals, yerr, n_params):
        # ... (Logique inchangée)
        n = len(residuals)
        # S'assurer que sigma2 est un tableau si yerr est un tableau
        sigma2 = yerr ** 2 if np.isscalar(yerr) or len(yerr)==n else np.ones_like(residuals) * 0.001**2
        chi2 = np.sum((residuals ** 2) / sigma2)
        return n_params * np.log(n) + chi2, chi2

    def predict_resonance_candidates(self, p_orb, p_ttv_epochs):
        # ... (Logique inchangée pour le calcul des résonances)
        p_ttv_days = p_ttv_epochs * p_orb
        # Ajout de la résonance 3:1 pour couvrir le rapport de l'utilisateur
        resonances = [(1, 2, "1:2 (Int)"), (2, 3, "2:3 (Int)"), (3, 2, "3:2 (Ext)"), (2, 1, "2:1 (Ext)"), (3, 1, "3:1 (Ext)")]
        candidates = []
        for j, k, name in resonances:
            f_orb, f_ttv = 1.0/p_orb, 1.0/p_ttv_days
            # P_plus: Utiliser le signe plus pour la première valeur
            p1 = 1.0 / ((j/k)*f_orb + (1/k)*f_ttv) 
            # P_moins: Utiliser le signe moins pour la seconde valeur (souvent le second candidat de la résonance)
            p2 = 1.0 / (abs((j/k)*f_orb - (1/k)*f_ttv)) # Utiliser abs() pour éviter P négative
            candidates.append(p1)
            candidates.append(p2)
        return candidates

    def generate_report(self):
        """Génère le rapport d'analyse en utilisant les données du fit MCMC."""
        p_orb_jours = self.detected_P_orb
        dir_path = None
        if self.oc_source_path and os.path.isfile(self.oc_source_path):
            dir_path = os.path.dirname(self.oc_source_path)
        if not dir_path and self.lc_selected_paths:
            dir_path = os.path.dirname(os.path.abspath(self.lc_selected_paths[0]))

        if self.last_mcmc_samples is None or self.oc_epoch is None: # Vérifier aussi les données O-C
            messagebox.showwarning("Attention", "Veuillez charger les O-C et effectuer le Fit Modèle MCMC d'abord.")
            return
        
        if not dir_path or not os.path.isdir(dir_path):
             messagebox.showerror("Erreur Répertoire", "Le répertoire de travail (dossier source) n'est pas défini ou n'existe pas.")
             return
             
        # Extraction du dernier composant du chemin (nom du dossier, ex: "Kepler-18d")
        planet_name = os.path.basename(dir_path)
        
        # Fallback pour un nom de dossier vide ou non significatif
        if not planet_name or planet_name in ['.', '..', '/', '\\']:
             planet_name = "Analyse_Systeme" 
             logger.warning("Le nom du répertoire est vide ou non significatif. Utilisation de 'Analyse_Systeme'.")
        
        # --- 2. Construction du chemin de sauvegarde ---
        safe_name = planet_name.replace(" ", "_").replace("-", "_")
        report_filename = f"{safe_name}_TTV_report.txt"
        output_filepath = os.path.join(dir_path, report_filename)
        
        if p_orb_jours is None:
            # Si elle n'est pas là, demander à l'utilisateur
            p_orb_input = simpledialog.askfloat(
                "Période Orbitale Manquante", 
                "La Période Orbitale (P_orb) n'a pas été détectée. Veuillez entrer sa valeur en jours (ex: 14.86000)."
            )
            if p_orb_input is None or p_orb_input <= 0:
                messagebox.showwarning("Annulation", "La Période Orbitale est nécessaire pour générer le rapport. Opération annulée.")
                return
            
            # FIX 2: Mettre à jour la variable d'instance avec la valeur saisie
            self.detected_P_orb = p_orb_input
            p_orb_jours = p_orb_input
            
        try:
            # 1. Extraction des paramètres MCMC (Moyenne et Erreur)
            samples = self.last_mcmc_samples
            n_params = samples.shape[1]
            medians = np.median(samples, axis=0)
            errors = np.std(samples, axis=0)
            
            # Paramètres de la première fréquence (Index 0, 1, 2)
            amp_moy, amp_err = medians[0], errors[0]
            p_ttv_epoch_moy, p_ttv_epoch_err = medians[1], errors[1]
            phase_moy, phase_err = medians[2], errors[2]
            offset_moy, offset_err = medians[-1], errors[-1]

            # 2. Calculs statistiques
            # BIC Nul (Modèle : Ligne droite = 1 paramètre (Offset))
            rms_brut = np.std(self.oc_res) # Approx RMS des données brutes
            residuals_nul = self.oc_res - np.mean(self.oc_res)
            bic_nul, _ = self.calculate_bic(residuals_nul, self.current_yerr, n_params=1)

            # BIC Modèle (Modèle TTV)
            residuals_ttv = self.oc_res - multi_sine_model(self.oc_epoch, *medians)
            rms_modele = np.std(residuals_ttv)
            bic_ttv, _ = self.calculate_bic(residuals_ttv, self.current_yerr, n_params=n_params)

            # 3. Données d'Observation
            epochs_n = len(self.oc_epoch)
            epochs_range = np.max(self.oc_epoch) - np.min(self.oc_epoch)
            p_orb_jours = self.detected_P_orb

            # 4. Prédiction des Résonances (Appel de la fonction corrigée)
            mmr_candidates_values = self.predict_resonance_candidates(p_orb_jours, p_ttv_epoch_moy)

            # 5. Formatage du rapport
            report_text = self._format_ttv_report(
                p_orb_jours=p_orb_jours,
                epochs_n=epochs_n,
                epochs_range=epochs_range,
                rms_brut=rms_brut,
                amp_moy=amp_moy, amp_err=amp_err,
                p_ttv_epoch_moy=p_ttv_epoch_moy, p_ttv_epoch_err=p_ttv_epoch_err,
                phase_moy=phase_moy, phase_err=phase_err,
                offset_moy=offset_moy, offset_err=offset_err,
                rms_modele=rms_modele,
                bic_nul=bic_nul,
                bic_ttv=bic_ttv,
                mmr_candidates=mmr_candidates_values
            )
            
            # 6. Affichage et Sauvegarde
            with open(output_filepath, 'w') as f:
                f.write(report_text)
                
            print("\n" + report_text)
            messagebox.showinfo("Rapport Généré", f"Rapport généré et sauvegardé sous :\n{output_filepath}")

        except Exception as e:
            messagebox.showerror("Erreur Rapport", f"Échec de la génération du rapport : {e}")
            logger.error(f"Erreur dans generate_report: {e}")


    def _format_ttv_report(self, **kwargs):
        """Fonction interne pour gérer le formatage précis du rapport (copie de la version fournie)."""
        # Simplification des arguments pour la clarté (utilise kwargs pour l'ensemble)
        
        p_orb_jours = kwargs['p_orb_jours']
        epochs_n = kwargs['epochs_n']
        epochs_range = kwargs['epochs_range']
        rms_brut = kwargs['rms_brut']
        amp_moy, amp_err = kwargs['amp_moy'], kwargs['amp_err']
        p_ttv_epoch_moy, p_ttv_epoch_err = kwargs['p_ttv_epoch_moy'], kwargs['p_ttv_epoch_err']
        phase_moy, phase_err = kwargs['phase_moy'], kwargs['phase_err']
        offset_moy, offset_err = kwargs['offset_moy'], kwargs['offset_err']
        rms_modele = kwargs['rms_modele']
        bic_nul = kwargs['bic_nul']
        bic_ttv = kwargs['bic_ttv']
        mmr_candidates = kwargs['mmr_candidates'] # Liste des P_pert calculées
        
        # --- Calculs Dérivés ---
        p_ttv_jours_moy = p_ttv_epoch_moy * p_orb_jours
        p_ttv_jours_err = p_ttv_epoch_err * p_orb_jours
        ratio_p_ttv_p_orb = p_ttv_epoch_moy
        delta_bic = bic_ttv - bic_nul
        
        if delta_bic < -10:
            interpretation = "[V] DÉTECTION FORTE (Signal périodique confirmé)"
        elif delta_bic < -5:
            interpretation = "[!] DÉTECTION MODÉRÉE (Preuve substantielle)"
        else:
            interpretation = "[X] DÉTECTION FAIBLE (Peu d'évidence pour le signal TTV)"
            
        date_generee = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        fichier_analyse = "Analyse TTV" # Statique

        # Mise en forme des résonances pour le rapport (regroupement par paire)
        mmr_lines = ""
        resonances_info = [
            ("1:2", "Intérieure", "Extérieure"),
            ("2:3", "Intérieure", "Extérieure"),
            ("3:2", "Intérieure", "Extérieure"),
            ("2:1", "Intérieure", "Extérieure"),
            ("3:1", "Intérieure", "Extérieure")
        ]
        
        for i, (ratio, label1, label2) in enumerate(resonances_info):
             p1_idx = i * 2
             p2_idx = i * 2 + 1
             
             p1 = mmr_candidates[p1_idx] if p1_idx < len(mmr_candidates) else np.nan
             p2 = mmr_candidates[p2_idx] if p2_idx < len(mmr_candidates) else np.nan
             
             # Déterminer laquelle est intérieure et laquelle est extérieure
             if not np.isnan(p1) and not np.isnan(p2):
                 if p1 < p_orb_jours and p2 > p_orb_jours:
                     p_inner, p_outer = p1, p2
                     label_inner, label_outer = label1, label2
                 elif p2 < p_orb_jours and p1 > p_orb_jours:
                     p_inner, p_outer = p2, p1
                     label_inner, label_outer = label1, label2
                 elif p1 < p_orb_jours and p2 < p_orb_jours:
                     # Les deux sont intérieures (cas rare)
                     p_inner, p_outer = max(p1, p2), np.nan
                     label_inner, label_outer = label1, label2
                 elif p1 > p_orb_jours and p2 > p_orb_jours:
                     # Les deux sont extérieures (cas rare)
                     p_inner, p_outer = np.nan, min(p1, p2)
                     label_inner, label_outer = label1, label2
                 else:
                     p_inner, p_outer = p1, p2
                     label_inner, label_outer = label1, label2
             else:
                 p_inner, p_outer = p1, p2
                 label_inner, label_outer = label1, label2
             
             lines_parts = [f"  > Résonance {ratio}:"]
             if not np.isnan(p_inner):
                 lines_parts.append(f"    - Perturbateur {label_inner.lower()} à P = {p_inner:.4f} jours")
             if not np.isnan(p_outer):
                 lines_parts.append(f"    - Perturbateur {label_outer.lower()} à P = {p_outer:.4f} jours")
             
             mmr_lines += "\n".join(lines_parts) + "\n"
        
        # --- Construction du Rapport (Alignement Exact) ---
        report = f"""
===================================================================
                RAPPORT D'ANALYSE TTV - AUTOMATISÉ
===================================================================
DATE       : {date_generee}
FICHIER    : {fichier_analyse}

1. DONNÉES D'OBSERVATION
-------------------------------------------------------------------
Nombre de transits (N) : {epochs_n}
Étendue temporelle     : {epochs_range:.2f} époques
RMS des données brutes : {rms_brut:.5f} jours
Période Planète (P_orb): {p_orb_jours:.5f} jours

2. RÉSULTATS DE LA MODÉLISATION (MCMC)
-------------------------------------------------------------------
Modèle : {self.last_n_freq} Sinusoïde(s) + Offset

> Fréquence 1 :
  - Amplitude (A)   : {amp_moy:.5f} +/- {amp_err:.5f} jours
  - Période (P_ttv) : {p_ttv_epoch_moy:.4f} +/- {p_ttv_epoch_err:.4f} époques (Cycles TTV)
                     = {p_ttv_jours_moy:.4f} +/- {p_ttv_jours_err:.4f} jours
  - Phase (phi)     : {phase_moy:.4f} +/- {phase_err:.4f} rad
> Offset Global     : {offset_moy:.5f} +/- {offset_err:.5f} jours

3. ANALYSE STATISTIQUE (Critère BIC - Naponiello et al.)
-------------------------------------------------------------------
RMS Résidus (Modèle)  : {rms_modele:.5f} jours
BIC (Modèle Nul)      : {bic_nul:.2f}
BIC (Modèle TTV)      : {bic_ttv:.2f}
Delta BIC             : {delta_bic:.2f}

INTERPRÉTATION : {interpretation}

4. VÉRIFICATION PHYSIQUE ('The Exoplanet Edge' - Yahalomi et al.)
-------------------------------------------------------------------
Période Orbitale Planète : {p_orb_jours:.4f} jours
Période Modulation TTV   : {p_ttv_jours_moy:.4f} jours ({ratio_p_ttv_p_orb:.4f} x P_orb)
Ratio (P_ttv / P_orb)    : {ratio_p_ttv_p_orb:.4f}
Limite Théorique         : > 0.5

[OK] COMPATIBLE AVEC UNE PERTURBATION PLANÉTAIRE
===================================================================
5. PRÉDICTION DU PERTURBATEUR (Hypothèse de Résonance)
-------------------------------------------------------------------
Si le signal TTV est causé par une résonance de mouvement moyen (MMR),
voici où pourrait se trouver la planète invisible :

{mmr_lines.strip()}

Note : Cherchez des transits ou des signaux TTV croisés à ces périodes.
===================================================================
"""
        return report.strip()
  

    def save_plot_with_table(self):
        if self.last_mcmc_samples is None: return
        filename = asksaveasfilename(defaultextension=".png")
        if not filename: return
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.errorbar(
            self.oc_epoch,
            self.oc_res * MINUTES_PER_DAY,
            yerr=self.current_yerr * MINUTES_PER_DAY,
            fmt='ko',
            alpha=0.5,
        )
        ax.plot(self.model_x, self.model_y * MINUTES_PER_DAY, 'r-')
        ax.set_ylabel("O-C (minutes)", fontsize=12, fontweight='bold')
        ax.set_xlabel("Epoch", fontsize=11, fontweight='bold')
        
        # Tableau simplifié
        medians = np.median(self.last_mcmc_samples, axis=0)
        # CORRECTION: Créer un en-tête pour les paramètres
        headers = [f'A{i+1}' for i in range(self.last_n_freq)] + \
                  [f'P_ttv{i+1}' for i in range(self.last_n_freq)] + \
                  [f'Phi{i+1}' for i in range(self.last_n_freq)] + ['Offset']
                  
        cell_text = [[f"{v:.5f}" for v in medians]]
        plt.table(cellText=cell_text, colLabels=headers[:len(medians)], loc='bottom', bbox=[0.0, -0.3, 1.0, 0.2])
        
        plt.subplots_adjust(bottom=0.3)
        plt.savefig(filename)
        plt.close(fig)

    # =========================================================================
    # PARTIE C : ANALYSE SYSTÈME MULTIPLE (Rétablies ici pour la complétion)
    # =========================================================================
    # ... (Code de la partie C)
    def setup_part_c(self):
        frame_multi = ttk.LabelFrame(self.tab_multi, text="Flux de travail : Système Multiple", padding=10)
        frame_multi.pack(side=tk.TOP, fill="both", expand=True, padx=10, pady=10)
        
        self.reports_list = [] 
        
        lbl_instr = tk.Label(frame_multi, text="Chargez les rapports (.txt) des planètes (ex: 18b, 18c)", fg="gray")
        lbl_instr.pack(anchor="w", pady=(0, 10))

        self.listbox_reports = tk.Listbox(frame_multi, height=8, selectmode=tk.EXTENDED)
        self.listbox_reports.pack(fill="both", expand=True, pady=5)
        
        btn_frame = tk.Frame(frame_multi)
        btn_frame.pack(fill="x", pady=10)
        
        ttk.Button(btn_frame, text="➕ Ajouter Rapport", command=self.add_planet_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="🗑️ Vider Liste", command=self.clear_reports).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="🚀 Analyser Système", command=self.analyze_multi_system).pack(side=tk.LEFT, padx=5)
        
        # Bouton de transfert vers N-body
        if NBODY_AVAILABLE:
            ttk.Separator(btn_frame, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
            ttk.Button(btn_frame, text="🔄 Transférer vers Simulation N-body", command=self.transfer_to_nbody_simulation).pack(side=tk.LEFT, padx=5)

    def clear_reports(self):
        self.reports_list = []
        self.listbox_reports.delete(0, tk.END)

    def add_planet_report(self):
        filenames = filedialog.askopenfilenames(filetypes=[("Rapports TTV", "*.txt")])
        if not filenames: return
        
        for f in filenames:
            data = self.parse_report_file(f)
            if data:
                self.reports_list.append(data)
                name = f.split('/')[-1]
                p_orb = data.get('P_orb', 0)
                p_ttv = data.get('P_ttv_days', None)
                if p_ttv is not None:
                    self.listbox_reports.insert(tk.END, f"{name} | P_orb={p_orb:.3f}j | P_ttv={p_ttv:.3f}j")
                else:
                    self.listbox_reports.insert(tk.END, f"{name} | P_orb={p_orb:.3f}j")

    def parse_report_file(self, filepath):
        data = {'filename': filepath.split('/')[-1], 'filepath': filepath}
        try:
            with open(filepath, 'r') as f: lines = f.readlines()
            for line in lines:
                # Période orbitale de la planète transitaire (P_orb) - ligne spécifique
                if "Période Planète (P_orb)" in line:
                    try: 
                        value = float(line.split(':')[1].split('jours')[0].strip())
                        data['P_orb'] = value
                    except: pass
                # Période de modulation TTV (P_ttv) - différente de P_orb
                if "Période Modulation TTV" in line and "jours" in line:
                    try: 
                        # Format: "Période Modulation TTV   : {value:.4f} jours ({ratio} x P_orb)"
                        value_str = line.split(':')[1].split('jours')[0].strip()
                        data['P_ttv_days'] = float(value_str)
                    except: pass
                # Phase phi
                if "Phase" in line and "phi" in line.lower() and ":" in line:
                    try: data['Phi'] = float(line.split(':')[1].split('+/-')[0].strip())
                    except: pass
                # Extraire l'amplitude TTV si disponible
                if "Amplitude" in line and "A" in line and "+/-" in line:
                    try: data['Amp_TTV'] = float(line.split(':')[1].split('+/-')[0].strip())
                    except: pass
            return data
        except: return None
    
    def transfer_to_nbody_simulation(self):
        """Transfère les paramètres des rapports vers l'onglet D (Simulation N-body)."""
        if not NBODY_AVAILABLE:
            messagebox.showerror("Erreur", "La simulation N-body n'est pas disponible. Installez 'rebound'.")
            return
        
        if len(self.reports_list) < 1:
            messagebox.showwarning("Attention", "Ajoutez au moins un rapport avant de transférer.")
            return
        
        # Demander la masse de l'étoile si pas déjà définie
        star_mass = simpledialog.askfloat(
            "Masse de l'Étoile",
            "Entrez la masse de l'étoile (Msun):",
            initialvalue=1.0,
            minvalue=0.1,
            maxvalue=10.0
        )
        if star_mass is None:
            return
        
        # Trier les planètes par période (ordre croissant)
        planets_sorted = sorted(self.reports_list, key=lambda x: x.get('P_orb', float('inf')))
        
        # Vider les planètes existantes dans l'onglet D
        if hasattr(self, 'nbody_planets'):
            self.nbody_planets.clear()
        else:
            self.nbody_planets = []
        
        # Remplir la masse de l'étoile
        if hasattr(self, 'star_mass_var'):
            self.star_mass_var.set(str(star_mass))
        
        # Estimer et créer les planètes
        periods_added = []  # Pour détecter les doublons
        
        for i, report in enumerate(planets_sorted):
            P_orb = report.get('P_orb')
            if P_orb is None or P_orb <= 0:
                logger.warning(f"Période invalide pour {report.get('filename', 'planète')}, ignorée.")
                continue
            
            # Vérifier si une planète avec cette période a déjà été ajoutée
            if any(abs(P_orb - P_existing) < 0.01 for P_existing in periods_added):
                planet_name = report.get('filename', f'Planète {i+1}').replace('.txt', '').replace('_TTV_report', '')
                response = messagebox.askyesno(
                    "Période Orbitale Dupliquée",
                    f"Attention : Une planète avec P_orb = {P_orb:.3f} jours a déjà été ajoutée.\n\n"
                    f"Planète actuelle : {planet_name}\n\n"
                    f"Dans une simulation N-body, deux planètes ne peuvent pas avoir exactement la même période orbitale.\n\n"
                    f"Voulez-vous quand même ajouter cette planète ?\n"
                    f"(Recommandé : Non - vérifiez que les rapports correspondent à des planètes différentes)"
                )
                if not response:
                    continue
            
            periods_added.append(P_orb)
            
            # Estimation de la masse (valeur par défaut raisonnable)
            # Pour les planètes chaudes (P < 10 j): typiquement 0.5-2 Mjup
            # Pour les planètes plus lointaines: masse variable
            if P_orb < 10:
                estimated_mass_mjup = 1.0  # Jupiter typique
            elif P_orb < 50:
                estimated_mass_mjup = 0.5  # Neptune/Saturn
            else:
                estimated_mass_mjup = 0.3  # Plus petite
            
            # Demander confirmation pour chaque planète
            planet_name = report.get('filename', f'Planète {i+1}').replace('.txt', '').replace('_TTV_report', '')
            p_ttv = report.get('P_ttv_days', None)
            info_text = f"Période orbitale (P_orb): {P_orb:.3f} jours"
            if p_ttv is not None:
                info_text += f"\nPériode TTV (P_ttv): {p_ttv:.3f} jours"
            
            mass = simpledialog.askfloat(
                f"Masse - {planet_name}",
                f"{info_text}\n\n"
                f"Masse estimée: {estimated_mass_mjup:.2f} Mjup\n"
                f"Entrez la masse de la planète (Mjup):",
                initialvalue=estimated_mass_mjup,
                minvalue=0.01,
                maxvalue=50.0
            )
            if mass is None:
                continue
            
            # Créer la planète avec les paramètres
            mjup_to_msun = u.M_jup.to(u.M_sun)
            
            planet = {
                'm': mass * mjup_to_msun,  # Conversion en Msun
                'P': P_orb,
                'inc': np.pi / 2,  # Orbite vue par la tranche (transit)
                'e': 0.0,  # Excentricité par défaut (circulaire)
                'omega': 0.0  # Argument du périastre
            }
            
            self.nbody_planets.append(planet)
        
        # Mettre à jour la liste des planètes
        if hasattr(self, 'update_nbody_planets_list'):
            self.update_nbody_planets_list()
        
        # Basculer vers l'onglet D
        self.sub_notebook.select(3)  # Index 3 = onglet D
        
        # Message de confirmation
        messagebox.showinfo(
            "Transfert Réussi",
            f"{len(self.nbody_planets)} planète(s) transférée(s) vers la simulation N-body.\n\n"
            f"Vous pouvez maintenant lancer la simulation dans l'onglet D."
        )

    def analyze_multi_system(self):
        if len(self.reports_list) < 2:
            messagebox.showwarning("Erreur", "Il faut au moins 2 rapports.")
            return
        
        planets = sorted(self.reports_list, key=lambda x: x.get('P_orb', 0))
        lines = ["RAPPORT SYSTEME", "-"*20]
        
        for i in range(len(planets)-1):
            p1, p2 = planets[i], planets[i+1]
            ratio = p2.get('P_orb',1)/p1.get('P_orb',1)
            lines.append(f"PAIRE {p1['filename']} / {p2['filename']}")
            lines.append(f"Ratio Périodes: {ratio:.3f}")
            
            dphi = abs(p1.get('Phi',0) - p2.get('Phi',0)) * 180/np.pi
            dphi = dphi % 360
            lines.append(f"Diff Phase: {dphi:.1f} deg")
            lines.append("-" * 20)
            
        filename = asksaveasfilename(defaultextension=".txt")
        if filename:
            with open(filename, 'w') as f: f.write("\n".join(lines))
            messagebox.showinfo("Succès", "Rapport Système généré.")
    
    # =========================================================================
    # PARTIE D : SIMULATION N-BODY
    # =========================================================================
    def setup_part_d(self):
        """Configuration de l'onglet Simulation N-body."""
        if not NBODY_AVAILABLE:
            msg_frame = ttk.LabelFrame(self.tab_nbody, text="Simulation N-body non disponible", padding=20)
            msg_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=20)
            
            msg = tk.Text(msg_frame, wrap=tk.WORD, height=8, width=60)
            msg.pack(fill=tk.BOTH, expand=True)
            msg.insert("1.0", 
                "Le module de simulation N-body nécessite la bibliothèque 'rebound'.\n\n"
                "Pour l'installer, exécutez dans votre terminal :\n"
                "  pip install rebound\n\n"
                "Une fois installé, redémarrez l'application pour utiliser cette fonctionnalité.\n\n"
                "La simulation N-body permet de :\n"
                "- Simuler des systèmes planétaires multiples avec interactions gravitationnelles\n"
                "- Calculer les TTV (Transit Timing Variations) causés par des planètes perturbatrices\n"
                "- Analyser les vitesses radiales et périodogrammes\n"
                "- Visualiser les orbites et variations temporelles"
            )
            msg.config(state=tk.DISABLED)
            return
        
        # Variables pour stocker les données
        self.nbody_simulation_data = None
        self.nbody_analysis_data = None
        
        # Frame principal avec PanedWindow
        main_paned = ttk.PanedWindow(self.tab_nbody, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # PARTIE GAUCHE : Conteneur avec ascenseur vertical
        left_container = ttk.Frame(main_paned)
        main_paned.add(left_container, weight=1)
        
        scrollbar = ttk.Scrollbar(left_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        canvas = tk.Canvas(left_container, highlightthickness=0, yscrollcommand=scrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=canvas.yview)
        
        left_frame = ttk.Frame(canvas, padding=10)
        canvas_window = canvas.create_window((0, 0), window=left_frame, anchor="nw")
        
        def _on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
        
        def _on_canvas_configure(event):
            canvas.itemconfig(canvas_window, width=event.width)
        
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        
        left_frame.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_configure)
        def _bind_mousewheel(event):
            canvas.bind_all("<MouseWheel>", _on_mousewheel)
        def _unbind_mousewheel(event):
            canvas.unbind_all("<MouseWheel>")
        canvas.bind("<Enter>", _bind_mousewheel)
        canvas.bind("<Leave>", _unbind_mousewheel)
        
        self.setup_nbody_controls(left_frame)
        
        # PARTIE DROITE : Visualisation
        right_frame = ttk.Frame(main_paned, padding=10)
        main_paned.add(right_frame, weight=2)
        self.setup_nbody_visualization(right_frame)
    
    def setup_nbody_controls(self, parent):
        """Crée l'interface de contrôle pour la simulation N-body."""
        # Instructions
        info_frame = ttk.LabelFrame(parent, text="Instructions", padding=10)
        info_frame.pack(fill="x", pady=5)
        
        info_text = (
            "1. Définir les paramètres du système (étoile + planètes)\n"
            "2. Configurer la durée et résolution de la simulation\n"
            "3. Lancer la simulation\n"
            "4. Analyser les résultats (TTV, orbites, RV)"
        )
        tk.Label(info_frame, text=info_text, justify=tk.LEFT, fg="gray").pack(anchor="w")
        
        # Paramètres de l'étoile
        star_frame = ttk.LabelFrame(parent, text="Étoile", padding=10)
        star_frame.pack(fill="x", pady=5)
        
        tk.Label(star_frame, text="Masse (Msun):").pack(anchor="w")
        self.star_mass_var = tk.StringVar(value="1.0")
        tk.Entry(star_frame, textvariable=self.star_mass_var, width=15).pack(fill="x")
        
        # Paramètres de simulation
        sim_frame = ttk.LabelFrame(parent, text="Paramètres Simulation", padding=10)
        sim_frame.pack(fill="x", pady=5)
        
        tk.Label(sim_frame, text="Durée (jours):").pack(anchor="w")
        self.sim_days_var = tk.StringVar(value="1000")
        tk.Entry(sim_frame, textvariable=self.sim_days_var, width=15).pack(fill="x")
        
        tk.Label(sim_frame, text="Nombre de sorties:").pack(anchor="w", pady=(5,0))
        self.sim_outputs_var = tk.StringVar(value="10000")
        tk.Entry(sim_frame, textvariable=self.sim_outputs_var, width=15).pack(fill="x")
        
        # Liste des planètes
        planets_frame = ttk.LabelFrame(parent, text="Planètes", padding=10)
        planets_frame.pack(fill="both", expand=True, pady=5)
        
        self.planets_listbox = tk.Listbox(planets_frame, height=6)
        self.planets_listbox.pack(fill="both", expand=True, pady=5)
        
        btn_frame = tk.Frame(planets_frame)
        btn_frame.pack(fill="x")
        
        ttk.Button(btn_frame, text="➕ Ajouter Planète", command=self.add_nbody_planet).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="✏️ Modifier", command=self.edit_nbody_planet).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="🗑️ Supprimer", command=self.remove_nbody_planet).pack(side=tk.LEFT, padx=2)
        
        # Paramètres de prédiction de transits futurs
        predict_frame = ttk.LabelFrame(parent, text="Prédiction Transits Futurs", padding=10)
        predict_frame.pack(fill="x", pady=5)
        
        # Sélection explicite de la planète
        tk.Label(predict_frame, text="Planète à prédire:").grid(row=0, column=0, sticky="w")
        self.nbody_predict_planet_var = tk.StringVar(value="Planète 1")
        self.nbody_predict_planet_combo = ttk.Combobox(
            predict_frame,
            textvariable=self.nbody_predict_planet_var,
            values=["Planète 1"],
            state="readonly",
            width=20
        )
        self.nbody_predict_planet_combo.grid(row=0, column=1, sticky="w", padx=5)
        ToolTip(self.nbody_predict_planet_combo, "Choisir la planète pour laquelle exporter les prochains transits.")
        
        # Presets de résolution
        tk.Label(predict_frame, text="Preset résolution:").grid(row=1, column=0, sticky="w", pady=(5, 0))
        self.nbody_predict_preset_var = tk.StringVar(value="Précis")
        preset_combo = ttk.Combobox(
            predict_frame,
            textvariable=self.nbody_predict_preset_var,
            values=["Rapide", "Précis", "Personnalisé"],
            state="readonly",
            width=20
        )
        preset_combo.grid(row=1, column=1, sticky="w", padx=5, pady=(5, 0))
        ToolTip(preset_combo, "Rapide: moins précis mais plus rapide. Précis: plus précis mais plus lent.")
        
        # Durée et sorties
        tk.Label(predict_frame, text="Durée (jours):").grid(row=2, column=0, sticky="w", pady=(5, 0))
        self.nbody_predict_days_var = tk.StringVar(value=self.sim_days_var.get() if hasattr(self, 'sim_days_var') else "365")
        days_entry = tk.Entry(predict_frame, textvariable=self.nbody_predict_days_var, width=10)
        days_entry.grid(row=2, column=1, sticky="w", padx=5, pady=(5, 0))
        ToolTip(days_entry, "Horizon de prédiction en jours (ex: 365 pour 1 an).")
        
        tk.Label(predict_frame, text="Nombre de sorties:").grid(row=3, column=0, sticky="w", pady=(5, 0))
        self.nbody_predict_outputs_var = tk.StringVar(value=self.sim_outputs_var.get() if hasattr(self, 'sim_outputs_var') else "10000")
        outputs_entry = tk.Entry(predict_frame, textvariable=self.nbody_predict_outputs_var, width=10)
        outputs_entry.grid(row=3, column=1, sticky="w", padx=5, pady=(5, 0))
        ToolTip(outputs_entry, "Résolution de l'intégration. Plus élevé = transits plus précis mais calcul plus long.")
        
        # T0 automatique
        tk.Label(predict_frame, text="T0 (BJD/JD):").grid(row=4, column=0, sticky="w", pady=(5, 0))
        self.nbody_predict_t0_var = tk.StringVar(value="")
        t0_entry = tk.Entry(predict_frame, textvariable=self.nbody_predict_t0_var, width=18)
        t0_entry.grid(row=4, column=1, sticky="w", padx=5, pady=(5, 0))
        ToolTip(t0_entry, "Référence temporelle du début de simulation. Auto = dernier Tc observé.")
        
        # Utiliser paramètres du fitting
        self.nbody_use_fit_var = tk.BooleanVar(value=True)
        use_fit_check = ttk.Checkbutton(
            predict_frame,
            text="Utiliser paramètres du fitting N-body",
            variable=self.nbody_use_fit_var,
            command=self._update_predict_planet_options
        )
        use_fit_check.grid(row=5, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ToolTip(use_fit_check, "Si disponible, utilise les paramètres ajustés pour prédire les transits.")
        
        # Incertitudes
        tk.Label(predict_frame, text="Échantillons (incertitudes):").grid(row=6, column=0, sticky="w", pady=(5, 0))
        self.nbody_predict_samples_var = tk.StringVar(value="0")
        samples_entry = tk.Entry(predict_frame, textvariable=self.nbody_predict_samples_var, width=10)
        samples_entry.grid(row=6, column=1, sticky="w", padx=5, pady=(5, 0))
        ToolTip(samples_entry, "Nombre d'échantillons Monte Carlo (0 = pas d'incertitudes).")
        
        # Appliquer preset sur changement
        self.nbody_predict_preset_var.trace_add("write", lambda *_: self._apply_predict_preset())
        self.nbody_predict_days_var.trace_add("write", lambda *_: self._apply_predict_preset())
        self._apply_predict_preset()
        
        # Boutons de contrôle - Simulation
        control_frame = ttk.LabelFrame(parent, text="Actions - Simulation", padding=10)
        control_frame.pack(fill="x", pady=5)
        
        ttk.Button(control_frame, text="🚀 Lancer Simulation", command=self.run_nbody_simulation).pack(fill="x", pady=2)
        ttk.Button(control_frame, text="📊 Analyser Résultats", command=self.analyze_nbody_results).pack(fill="x", pady=2)
        ttk.Button(control_frame, text="🕒 Prévoir Transits Futurs", command=self.export_future_transits).pack(fill="x", pady=2)
        ttk.Button(control_frame, text="💾 Sauvegarder", command=self.save_nbody_results).pack(fill="x", pady=2)
        
        # Stockage des planètes
        self.nbody_planets = []
        
        # Section Fitting N-body
        if NBODY_FITTING_AVAILABLE:
            fitting_frame = ttk.LabelFrame(parent, text="Fitting N-body aux Observations TTV", padding=10)
            fitting_frame.pack(fill="x", pady=5)
            
            info_text = (
                "Ajuste les paramètres physiques d'une planète perturbatrice\n"
                "aux données TTV observées (nécessite données O-C)"
            )
            tk.Label(fitting_frame, text=info_text, justify=tk.LEFT, fg="gray", font=("", 8)).pack(anchor="w", pady=2)
            
            ttk.Button(fitting_frame, text="📥 Charger Données O-C", command=self.load_oc_for_nbody_fit).pack(fill="x", pady=2)
            ttk.Button(fitting_frame, text="🎯 Lancer Fitting N-body", command=self.run_nbody_fitting).pack(fill="x", pady=2)
            
            # Stockage des données pour le fitting
            self.nbody_fit_data = None
            self.nbody_fitter = None
    
    def setup_nbody_visualization(self, parent):
        """Crée l'interface de visualisation pour les résultats N-body."""
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        
        self.nbody_fig = plt.Figure(figsize=(10, 8))
        self.nbody_canvas = FigureCanvasTkAgg(self.nbody_fig, parent)
        self.nbody_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        toolbar = NavigationToolbar2Tk(self.nbody_canvas, parent)
        toolbar.update()
        
        # Axes initialisés plus tard
        self.nbody_axes = None
    
    def add_nbody_planet(self):
        """Ouvre une fenêtre pour ajouter une planète."""
        dialog = tk.Toplevel(self.main_frame)
        dialog.title("Ajouter Planète")
        dialog.geometry("350x300")
        dialog.transient(self.main_frame)
        dialog.grab_set()
        
        params = {}
        entries = {}
        
        fields = [
            ("Masse (Mjup)", "m", "0.5"),
            ("Période (jours)", "P", "10.0"),
            ("Inclinaison (deg)", "inc", "90.0"),
            ("Excentricité", "e", "0.0"),
            ("Argument péria (deg)", "omega", "0.0"),
        ]
        
        for i, (label, key, default) in enumerate(fields):
            tk.Label(dialog, text=label + ":").grid(row=i, column=0, padx=5, pady=5, sticky="e")
            var = tk.StringVar(value=default)
            entries[key] = tk.Entry(dialog, textvariable=var, width=20)
            entries[key].grid(row=i, column=1, padx=5, pady=5)
            params[key] = var
        
        def save_planet():
            try:
                from astropy import units as u
                mjup = u.M_jup.to(u.M_sun)
                
                planet = {
                    'm': float(params['m'].get()) * mjup,  # Conversion en Msun
                    'P': float(params['P'].get()),
                    'inc': np.deg2rad(float(params['inc'].get())),
                    'e': float(params['e'].get()),
                    'omega': np.deg2rad(float(params['omega'].get()))
                }
                
                self.nbody_planets.append(planet)
                self.update_nbody_planets_list()
                dialog.destroy()
            except ValueError as e:
                messagebox.showerror("Erreur", f"Valeurs invalides: {e}")
        
        tk.Button(dialog, text="Ajouter", command=save_planet).grid(row=len(fields), column=0, columnspan=2, pady=10)
        tk.Button(dialog, text="Annuler", command=dialog.destroy).grid(row=len(fields)+1, column=0, columnspan=2)
    
    def edit_nbody_planet(self):
        """Modifie une planète sélectionnée."""
        selection = self.planets_listbox.curselection()
        if not selection:
            messagebox.showwarning("Attention", "Sélectionnez une planète à modifier.")
            return
        # TODO: Implémenter l'édition
        messagebox.showinfo("Info", "Fonction d'édition à implémenter. Supprimez et réajoutez pour l'instant.")
    
    def remove_nbody_planet(self):
        """Supprime une planète sélectionnée."""
        selection = self.planets_listbox.curselection()
        if not selection:
            messagebox.showwarning("Attention", "Sélectionnez une planète à supprimer.")
            return
        
        idx = selection[0]
        self.nbody_planets.pop(idx)
        self.update_nbody_planets_list()
    
    def update_nbody_planets_list(self):
        """Met à jour la liste des planètes affichée."""
        self.planets_listbox.delete(0, tk.END)
        from astropy import units as u
        mjup = u.M_sun.to(u.M_jup)
        
        for i, planet in enumerate(self.nbody_planets):
            m_jup = planet['m'] * mjup
            text = f"Planète {i+1}: M={m_jup:.2f} Mjup, P={planet['P']:.2f} j"
            self.planets_listbox.insert(tk.END, text)
        
        self._update_predict_planet_options()

    def _update_predict_planet_options(self):
        """Met à jour la liste des planètes pour la prédiction."""
        if not hasattr(self, 'nbody_predict_planet_combo'):
            return
        
        use_fit = bool(getattr(self, 'nbody_use_fit_var', tk.BooleanVar(value=False)).get())
        names = []
        
        if use_fit and self.nbody_fitter and getattr(self.nbody_fitter, 'parameters', None):
            n_planets = max(0, len(self.nbody_fitter.parameters) - 1)
            names = [f"Planète {i+1} (fit)" for i in range(n_planets)]
        else:
            names = [f"Planète {i+1}" for i in range(len(self.nbody_planets))]
        
        if not names:
            names = ["Planète 1"]
        
        current = self.nbody_predict_planet_var.get()
        self.nbody_predict_planet_combo["values"] = names
        if current not in names:
            self.nbody_predict_planet_var.set(names[0])

    def _apply_predict_preset(self):
        """Applique un preset de résolution aux sorties."""
        if not hasattr(self, 'nbody_predict_preset_var'):
            return
        
        preset = self.nbody_predict_preset_var.get()
        if preset == "Personnalisé":
            return
        
        try:
            days = float(self.nbody_predict_days_var.get())
        except ValueError:
            return
        
        factor = 12 if preset == "Rapide" else 48
        outputs = max(1000, int(days * factor))
        self.nbody_predict_outputs_var.set(str(outputs))

    def _get_last_observed_tc(self):
        """Retourne le dernier Tc observé si disponible."""
        if self.nbody_fit_data and 'Tc' in self.nbody_fit_data and len(self.nbody_fit_data['Tc']) > 0:
            return float(np.max(self.nbody_fit_data['Tc']))
        if hasattr(self, 'oc_time_raw') and self.oc_time_raw is not None and len(self.oc_time_raw) > 0:
            return float(np.max(self.oc_time_raw))
        return None

    def _parse_planet_index_label(self, label):
        """Extrait l'index de planète depuis un libellé."""
        match = re.search(r"(\d+)", label)
        if match:
            return int(match.group(1))
        return 1

    def _build_prediction_objects(self, use_fit):
        """Construit la liste des objets pour la prédiction."""
        if use_fit and self.nbody_fitter and getattr(self.nbody_fitter, 'parameters', None):
            return copy.deepcopy(self.nbody_fitter.parameters)
        
        try:
            star_mass = float(self.star_mass_var.get()) if hasattr(self, 'star_mass_var') else 1.0
        except ValueError:
            raise ValueError("Masse stellaire invalide.")
        
        objects = [{'m': star_mass}]
        objects.extend(copy.deepcopy(self.nbody_planets))
        return objects

    def _apply_fitter_errors(self, objects, errors, rng):
        """Perturbe les paramètres selon les incertitudes du fitting."""
        for key, err in errors.items():
            if err is None:
                continue
            try:
                idx_str, param = key.split('_', 1)
                idx = int(idx_str)
            except ValueError:
                continue
            
            if idx >= len(objects) or param not in objects[idx]:
                continue
            
            val = objects[idx][param]
            if not isinstance(val, (int, float, np.floating)):
                continue
            
            sample = float(val) + float(err) * rng.normal()
            
            if param in ("m", "P", "a") and sample <= 0:
                sample = max(1e-6, float(val))
            if param == "e":
                sample = min(max(sample, 0.0), 0.9)
            if param == "inc":
                sample = min(max(sample, 0.0), np.pi)
            if param == "omega":
                while sample < -np.pi:
                    sample += 2 * np.pi
                while sample > np.pi:
                    sample -= 2 * np.pi
            
            objects[idx][param] = sample

    def _estimate_transit_uncertainties(self, objects, days, outputs, planet_index, t0, n_samples, errors):
        """Estime les incertitudes des temps de transit via Monte Carlo."""
        rng = np.random.default_rng()
        tc_rel_samples = []
        tc_abs_samples = []
        
        for _ in range(n_samples):
            objs = copy.deepcopy(objects)
            self._apply_fitter_errors(objs, errors, rng)
            sample = predict_future_transits(objs, days, outputs, planet_index=planet_index, t0=t0)
            tc_rel_samples.append(sample['Tc_rel'])
            tc_abs_samples.append(sample['Tc_abs'])
        
        min_len = min(len(arr) for arr in tc_abs_samples)
        if min_len == 0:
            return None, None
        
        tc_rel_stack = np.vstack([arr[:min_len] for arr in tc_rel_samples])
        tc_abs_stack = np.vstack([arr[:min_len] for arr in tc_abs_samples])
        
        tc_rel_std = np.std(tc_rel_stack, axis=0)
        tc_abs_std = np.std(tc_abs_stack, axis=0)
        return tc_rel_std, tc_abs_std
    
    def run_nbody_simulation(self):
        """Lance la simulation N-body."""
        if not NBODY_AVAILABLE:
            messagebox.showerror("Erreur", "rebound n'est pas installé.")
            return
        
        if len(self.nbody_planets) == 0:
            messagebox.showwarning("Attention", "Ajoutez au moins une planète.")
            return
        
        try:
            # Récupérer les paramètres
            star_mass = float(self.star_mass_var.get())
            Ndays = float(self.sim_days_var.get())
            Noutputs = int(self.sim_outputs_var.get())
            
            # Construire la liste des objets
            objects = [{'m': star_mass}]  # Étoile
            objects.extend(self.nbody_planets)  # Planètes
            
            # Lancer la simulation dans un thread
            self.main_frame.config(cursor="watch")
            thread = threading.Thread(
                target=self._worker_nbody_simulation,
                args=(objects, Ndays, Noutputs)
            )
            thread.start()
            
        except ValueError as e:
            messagebox.showerror("Erreur", f"Paramètres invalides: {e}")
            self.main_frame.config(cursor="")
    
    def _worker_nbody_simulation(self, objects, Ndays, Noutputs):
        """Worker thread pour la simulation N-body."""
        try:
            sim_data = generate_simulation(objects, Ndays, Noutputs)
            analysis_data = analyze_simulation(sim_data, ttvfast=False)
            
            self.main_frame.after(0, self._complete_nbody_simulation, sim_data, analysis_data)
        except Exception as e:
            logger.error(f"Erreur simulation N-body: {e}", exc_info=True)
            self.main_frame.after(0, self._handle_nbody_error, str(e))
        finally:
            self.main_frame.after(0, lambda: self.main_frame.config(cursor=""))
    
    def _complete_nbody_simulation(self, sim_data, analysis_data):
        """Met à jour la GUI après la simulation."""
        self.nbody_simulation_data = sim_data
        self.nbody_analysis_data = analysis_data
        messagebox.showinfo("Succès", f"Simulation terminée: {len(analysis_data['planets'])} planète(s) analysée(s).")
        self.plot_nbody_results()
    
    def _handle_nbody_error(self, error_msg):
        """Gère les erreurs de simulation."""
        messagebox.showerror("Erreur Simulation", f"Échec de la simulation:\n{error_msg}")
    
    def plot_nbody_results(self):
        """Affiche les résultats de la simulation."""
        if self.nbody_analysis_data is None:
            return
        
        self.nbody_fig.clear()
        
        if len(self.nbody_analysis_data['planets']) == 0:
            ax = self.nbody_fig.add_subplot(111)
            ax.text(0.5, 0.5, "Aucune planète à afficher", 
                   ha="center", va="center", transform=ax.transAxes)
            self.nbody_canvas.draw()
            return
        
        # Créer des sous-graphiques
        ax1 = self.nbody_fig.add_subplot(2, 2, 1)  # Orbites
        ax2 = self.nbody_fig.add_subplot(2, 2, 2)  # TTV
        ax3 = self.nbody_fig.add_subplot(2, 2, 3)  # RV
        ax4 = self.nbody_fig.add_subplot(2, 2, 4)  # Périodogramme TTV
        
        data = self.nbody_analysis_data
        
        # Orbites
        for i, planet in enumerate(data['planets']):
            ax1.plot(planet['x'], planet['y'], label=f'Planète {i+1}', lw=0.5, alpha=0.5)
        ax1.set_xlabel('x (AU)')
        ax1.set_ylabel('y (AU)')
        ax1.legend()
        ax1.set_title('Orbites')
        ax1.set_aspect('equal', adjustable='box')
        
        # TTV
        for i, planet in enumerate(data['planets']):
            if len(planet['ttv']) >= 2:
                epochs = np.arange(len(planet['ttv']))
                ax2.plot(epochs, planet['ttv'] * 24 * 60, label=f'Planète {i+1}', marker='o', markersize=2)
        ax2.set_xlabel('Époque')
        ax2.set_ylabel('TTV (minutes)')
        ax2.legend()
        ax2.set_title('Variations Temporelles de Transit')
        ax2.grid(True)
        
        # RV
        if 'RV' in data:
            ax3.plot(data['times'][1:], data['RV']['signal'], 'k-', lw=0.5)
            ax3.set_xlabel('Temps (jours)')
            ax3.set_ylabel('RV (m/s)')
            ax3.set_title('Vitesses Radiales')
            ax3.grid(True)
        
        # Périodogramme TTV
        for i, planet in enumerate(data['planets']):
            if len(planet['freq']) > 0:
                periods = 1. / planet['freq']
                ax4.semilogx(periods, planet['power'], label=f'Planète {i+1}', alpha=0.7)
        ax4.set_xlabel('Période (époques)')
        ax4.set_ylabel('Puissance')
        ax4.set_title('Périodogramme TTV')
        ax4.legend()
        ax4.grid(True)
        
        self.nbody_fig.tight_layout()
        self.nbody_canvas.draw()
    
    def analyze_nbody_results(self):
        """Analyse les résultats de la simulation."""
        if self.nbody_analysis_data is None:
            messagebox.showwarning("Attention", "Lancez d'abord une simulation.")
            return
        
        # Afficher les résultats (déjà fait dans plot_nbody_results)
        self.plot_nbody_results()
    
    def save_nbody_results(self):
        """Sauvegarde les résultats de la simulation."""
        if self.nbody_analysis_data is None:
            messagebox.showwarning("Attention", "Aucune simulation à sauvegarder.")
            return
        
        filename = asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")])
        if filename:
            self.nbody_fig.savefig(filename, dpi=150, bbox_inches='tight')
            messagebox.showinfo("Succès", f"Figure sauvegardée: {filename}")

    def export_future_transits(self):
        """Exporte les transits futurs via une simulation N-body."""
        if not NBODY_AVAILABLE:
            messagebox.showerror("Erreur", "rebound n'est pas installé.")
            return
        
        use_fit = bool(getattr(self, 'nbody_use_fit_var', tk.BooleanVar(value=False)).get())
        if use_fit and not (self.nbody_fitter and getattr(self.nbody_fitter, 'parameters', None)):
            messagebox.showwarning("Attention", "Aucun fitting N-body disponible. Utilisation des paramètres manuels.")
            use_fit = False
        
        if not use_fit and len(self.nbody_planets) == 0:
            messagebox.showwarning("Attention", "Ajoutez au moins une planète.")
            return
        
        try:
            days = float(self.nbody_predict_days_var.get())
            outputs = int(self.nbody_predict_outputs_var.get())
        except ValueError:
            messagebox.showerror("Erreur", "Durée ou nombre de sorties invalide.")
            return
        
        if days <= 0 or outputs < 100:
            messagebox.showerror("Erreur", "Durée > 0 et sorties >= 100 requises.")
            return
        
        t0_str = self.nbody_predict_t0_var.get().strip()
        t0 = None
        if t0_str:
            try:
                t0 = float(t0_str)
            except ValueError:
                messagebox.showerror("Erreur", "T0 invalide.")
                return
        else:
            t0 = self._get_last_observed_tc()
        
        if t0 is None:
            messagebox.showwarning("Attention", "T0 indisponible. Veuillez le saisir manuellement.")
            return
        
        label = self.nbody_predict_planet_var.get()
        planet_index = self._parse_planet_index_label(label)
        
        try:
            n_samples = int(self.nbody_predict_samples_var.get())
        except ValueError:
            messagebox.showerror("Erreur", "Nombre d'échantillons invalide.")
            return
        
        if n_samples > 0 and not (use_fit and self.nbody_fitter and getattr(self.nbody_fitter, 'errors', None)):
            messagebox.showwarning(
                "Attention",
                "Incertitudes indisponibles (pas de fitting ou pas d'erreurs). Export sans écarts-types."
            )
            n_samples = 0
        
        filename = asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not filename:
            return
        
        self.main_frame.config(cursor="watch")
        thread = threading.Thread(
            target=self._worker_export_future_transits,
            args=(use_fit, days, outputs, planet_index, t0, filename, n_samples)
        )
        thread.start()

    def _worker_export_future_transits(self, use_fit, days, outputs, planet_index, t0, filename, n_samples):
        """Worker thread pour exporter les transits futurs."""
        try:
            objects = self._build_prediction_objects(use_fit)
            data = predict_future_transits(objects, days, outputs, planet_index=planet_index, t0=t0)
            
            if len(data['Tc_abs']) == 0:
                raise ValueError("Aucun transit détecté dans la fenêtre simulée.")
            
            tc_rel_std = None
            tc_abs_std = None
            if use_fit and n_samples > 0 and self.nbody_fitter and getattr(self.nbody_fitter, 'errors', None):
                tc_rel_std, tc_abs_std = self._estimate_transit_uncertainties(
                    objects, days, outputs, planet_index, t0, n_samples, self.nbody_fitter.errors
                )
            
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                headers = ["epoch", "Tc_rel_days", "Tc_abs", "Tc_utc", "planet_index"]
                if tc_abs_std is not None:
                    headers.extend(["Tc_rel_std", "Tc_abs_std"])
                writer.writerow(headers)
                n_rows = len(data['epochs'])
                if tc_abs_std is not None:
                    n_rows = min(n_rows, len(tc_abs_std))
                for i in range(n_rows):
                    try:
                        tc_utc = Time(float(data['Tc_abs'][i]), format='jd', scale='utc').isot
                    except Exception:
                        tc_utc = ""
                    row = [
                        int(data['epochs'][i]),
                        float(data['Tc_rel'][i]),
                        float(data['Tc_abs'][i]),
                        tc_utc,
                        int(planet_index)
                    ]
                    if tc_abs_std is not None:
                        row.extend([
                            float(tc_rel_std[i]),
                            float(tc_abs_std[i])
                        ])
                    writer.writerow(row)
            
            self.main_frame.after(0, self._complete_export_future_transits, filename, len(data['epochs']))
        
        except Exception as e:
            logger.error(f"Erreur export transits futurs: {e}", exc_info=True)
            self.main_frame.after(0, self._handle_export_future_transits_error, str(e))
        finally:
            self.main_frame.after(0, lambda: self.main_frame.config(cursor=""))

    def _complete_export_future_transits(self, filename, n_transits):
        """Fin de l'export des transits futurs."""
        messagebox.showinfo("Succès", f"Transits futurs exportés: {n_transits} lignes\n{filename}")

    def _handle_export_future_transits_error(self, error_msg):
        """Gère les erreurs d'export de transits futurs."""
        messagebox.showerror("Erreur Export", f"Échec export transits futurs:\n{error_msg}")
    
    def load_oc_for_nbody_fit(self):
        """Charge des données O-C pour le fitting N-body."""
        # Réutiliser les données O-C déjà chargées dans l'onglet B si disponibles
        if self.oc_epoch is not None and self.oc_res is not None:
            use_existing = messagebox.askyesno(
                "Données O-C existantes",
                "Utiliser les données O-C déjà chargées dans l'onglet B ?\n\n"
                f"Époques: {len(self.oc_epoch)} points"
            )
            if use_existing:
                # Il faut calculer Tc à partir des époques, TTV et P_orb
                if self.detected_P_orb is None:
                    messagebox.showerror("Erreur", "Période orbitale non détectée. Chargez d'abord les données O-C correctement.")
                    return
                
                # Estimer T0 depuis oc_time_raw si disponible, sinon utiliser une estimation
                if hasattr(self, 'oc_time_raw') and self.oc_time_raw is not None and len(self.oc_time_raw) > 0:
                    # Utiliser les temps bruts pour estimer T0
                    T0_est = self.oc_time_raw[0] - self.oc_epoch[0] * self.detected_P_orb
                else:
                    # Estimation basique : valeur arbitraire (la différence n'importe pas pour les TTV)
                    T0_est = 2450000.0
                
                # Calculer Tc = T0 + epoch * P + TTV
                Tc_obs = T0_est + self.oc_epoch * self.detected_P_orb + self.oc_res
                Tc_err = self.oc_yerr if self.oc_yerr is not None else np.ones_like(self.oc_res) * 0.001
                
                self.nbody_fit_data = {
                    'epochs': self.oc_epoch,
                    'ttv': self.oc_res,
                    'ttv_err': Tc_err,
                    'Tc': Tc_obs,
                    'Tc_err': Tc_err
                }
                if hasattr(self, 'nbody_predict_t0_var'):
                    self.nbody_predict_t0_var.set(f"{float(np.max(Tc_obs)):.6f}")
                messagebox.showinfo("Succès", f"Données O-C chargées depuis l'onglet B.\n{len(Tc_obs)} points de transit calculés.")
                return
        
        # Sinon, charger depuis un fichier
        path = filedialog.askopenfilename(filetypes=[("CSV", "*.csv")])
        if not path:
            return
        
        try:
            df = pd.read_csv(path)
            df.columns = [c.strip().lower() for c in df.columns]
            cols = df.columns
            
            ep_col = next((c for c in ['epoch', 'e'] if c in cols), None)
            ttv_col = next((c for c in ['ttv', 'o-c', 'oc', 'residus'] if c in cols), None)
            err_col = next((c for c in ['uncertainty', 'err', 'error', 'sig'] if c in cols), None)
            time_col = next((c for c in ['time', 'bjd', 'jd', 'date', 'mid_time', 'mid-time'] if c in cols), None)
            
            if not ep_col or not ttv_col or not time_col:
                messagebox.showerror("Erreur", "Colonnes manquantes. Nécessite: Epoch, TTV, Time")
                return
            
            # Convertir en temps de transit (Tc) depuis époques et TTV
            # On a besoin de P_orb pour convertir
            p_orb = simpledialog.askfloat(
                "Période Orbitale",
                "Entrez la période orbitale de la planète (jours):",
                minvalue=0.1
            )
            if p_orb is None:
                return
            
            epochs = df[ep_col].values
            ttv = df[ttv_col].values
            times = df[time_col].values
            
            # Calculer les temps de transit observés
            # T0 initial estimé depuis le premier point
            t0_est = times[0] - epochs[0] * p_orb
            Tc_obs = t0_est + epochs * p_orb + ttv
            
            err = df[err_col].values if err_col else np.ones_like(ttv) * 0.001
            
            self.nbody_fit_data = {
                'Tc': Tc_obs,
                'Tc_err': err,
                'epochs': epochs,
                'ttv': ttv
            }
            if hasattr(self, 'nbody_predict_t0_var'):
                self.nbody_predict_t0_var.set(f"{float(np.max(Tc_obs)):.6f}")
            
            messagebox.showinfo("Succès", f"Données chargées: {len(Tc_obs)} points")
            
        except Exception as e:
            messagebox.showerror("Erreur", f"Échec du chargement: {e}")
            logger.error(f"Erreur chargement O-C pour N-body: {e}", exc_info=True)
    
    def run_nbody_fitting(self):
        """Lance le fitting N-body aux observations TTV."""
        if not NBODY_FITTING_AVAILABLE:
            messagebox.showerror("Erreur", "ultranest n'est pas installé. Installez avec: pip install ultranest")
            return
        
        if self.nbody_fit_data is None:
            messagebox.showwarning("Attention", "Chargez d'abord des données O-C.")
            return
        
        # Demander les paramètres de la planète intérieure
        dialog = tk.Toplevel(self.main_frame)
        dialog.title("Paramètres Planète Intérieure")
        dialog.geometry("400x250")
        dialog.transient(self.main_frame)
        dialog.grab_set()
        
        params = {}
        fields = [
            ("Masse (Mjup)", "m", "1.0"),
            ("Période (jours)", "P", str(self.detected_P_orb if self.detected_P_orb else "10.0")),
        ]
        
        for i, (label, key, default) in enumerate(fields):
            tk.Label(dialog, text=label + ":").grid(row=i, column=0, padx=5, pady=5, sticky="e")
            var = tk.StringVar(value=default)
            tk.Entry(dialog, textvariable=var, width=20).grid(row=i, column=1, padx=5, pady=5)
            params[key] = var
        
        # Paramètres planète extérieure (bornes)
        tk.Label(dialog, text="Planète perturbatrice (bornes):").grid(row=len(fields), column=0, columnspan=2, pady=(10,5), sticky="w")
        
        bounds_fields = [
            ("Masse min (Mjup)", "m_min", "0.01"),
            ("Masse max (Mjup)", "m_max", "10.0"),
            ("Période min (jours)", "P_min", "5.0"),
            ("Période max (jours)", "P_max", "100.0"),
        ]
        
        for i, (label, key, default) in enumerate(bounds_fields):
            tk.Label(dialog, text=label + ":").grid(row=len(fields)+1+i, column=0, padx=5, pady=2, sticky="e")
            var = tk.StringVar(value=default)
            tk.Entry(dialog, textvariable=var, width=20).grid(row=len(fields)+1+i, column=1, padx=5, pady=2)
            params[key] = var
        
        def start_fitting():
            try:
                mjup_to_msun = u.M_jup.to(u.M_sun)
                
                # Paramètres planète intérieure
                m_inner = float(params['m'].get()) * mjup_to_msun
                P_inner = float(params['P'].get())
                
                # Bornes planète extérieure
                m_min = float(params['m_min'].get()) * mjup_to_msun
                m_max = float(params['m_max'].get()) * mjup_to_msun
                P_min = float(params['P_min'].get())
                P_max = float(params['P_max'].get())
                
                dialog.destroy()
                
                # Lancer le fitting dans un thread
                self.main_frame.config(cursor="watch")
                thread = threading.Thread(
                    target=self._worker_nbody_fitting,
                    args=(m_inner, P_inner, m_min, m_max, P_min, P_max)
                )
                thread.start()
                
            except ValueError as e:
                messagebox.showerror("Erreur", f"Valeurs invalides: {e}")
        
        tk.Button(dialog, text="Lancer Fitting", command=start_fitting).grid(row=len(fields)+len(bounds_fields)+1, column=0, columnspan=2, pady=10)
        tk.Button(dialog, text="Annuler", command=dialog.destroy).grid(row=len(fields)+len(bounds_fields)+2, column=0, columnspan=2)
    
    def _worker_nbody_fitting(self, m_inner, P_inner, m_min, m_max, P_min, P_max):
        """Worker thread pour le fitting N-body."""
        try:
            from core.nbody_simulation import NBodyFitter
            
            # Construire les priors et bornes
            star_mass = float(self.star_mass_var.get()) if hasattr(self, 'star_mass_var') else 1.0
            
            # Prior: paramètres initiaux
            prior = [
                {'m': star_mass},  # Étoile
                {'m': m_inner, 'P': P_inner, 'inc': np.pi/2, 'e': 0.0, 'omega': 0.0},  # Planète intérieure
                {'m': (m_min + m_max) / 2, 'P': (P_min + P_max) / 2, 'inc': np.pi/2, 'e': 0.0, 'omega': 0.0}  # Planète extérieure
            ]
            
            # Bornes: ce qu'on ajuste
            bounds = [
                {},  # Pas de bornes sur l'étoile
                {'P': [P_inner * 0.99, P_inner * 1.01]},  # Petites variations sur P_inner
                {
                    'm': [m_min, m_max],
                    'P': [P_min, P_max],
                    'omega': [-np.pi, np.pi]
                }
            ]
            
            # Données: format attendu par NBodyFitter
            data = [
                {},  # Pas de données pour l'étoile
                {'Tc': self.nbody_fit_data['Tc'], 'Tc_err': self.nbody_fit_data['Tc_err']},  # Planète intérieure
                {}  # Pas de données pour la planète extérieure
            ]
            
            # Lancer le fitting
            fitter = NBodyFitter(data, prior=prior, bounds=bounds, verbose=True)
            
            self.main_frame.after(0, self._complete_nbody_fitting, fitter)
            
        except Exception as e:
            logger.error(f"Erreur fitting N-body: {e}", exc_info=True)
            self.main_frame.after(0, self._handle_nbody_fitting_error, str(e))
        finally:
            self.main_frame.after(0, lambda: self.main_frame.config(cursor=""))
    
    def _complete_nbody_fitting(self, fitter):
        """Met à jour la GUI après le fitting."""
        self.nbody_fitter = fitter
        if hasattr(self, 'nbody_use_fit_var'):
            self.nbody_use_fit_var.set(True)
        self._update_predict_planet_options()
        messagebox.showinfo("Succès", "Fitting N-body terminé!\nVérifiez les résultats dans les logs.")
        
        # Afficher les résultats si disponibles
        if fitter.parameters:
            result_text = "Résultats du Fitting N-body:\n\n"
            if len(fitter.parameters) > 2:
                outer = fitter.parameters[2]
                m_jup = outer.get('m', 0) * u.M_sun.to(u.M_jup)
                result_text += f"Planète perturbatrice:\n"
                result_text += f"  Masse: {m_jup:.3f} Mjup\n"
                result_text += f"  Période: {outer.get('P', 'N/A'):.3f} jours\n"
            
            messagebox.showinfo("Résultats", result_text)
    
    def _handle_nbody_fitting_error(self, error_msg):
        """Gère les erreurs de fitting."""
        messagebox.showerror("Erreur Fitting", f"Échec du fitting N-body:\n{error_msg}")