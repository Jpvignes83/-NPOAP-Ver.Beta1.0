# gui/lc_markers_viewer.py
"""
Viewer de la courbe de lumière concaténée avec les marqueurs de période du périodogramme
reportés sur la courbe : une couleur par période, transits ajustables à la main.
Après validation : création d'un fichier mid-time.csv par période/couleur.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd
import logging
from pathlib import Path

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.ticker import MaxNLocator, ScalarFormatter

from gui.transit_fit_dialog import TransitFitDialog

logger = logging.getLogger(__name__)

PERIOD_COLORS = [
    '#e74c3c', '#27ae60', '#3498db', '#9b59b6', '#f39c12', '#1abc9c',
]

MAX_DISPLAY_POINTS = 80_000


def _downsample(time, flux, max_pts=MAX_DISPLAY_POINTS):
    n = len(time)
    if n <= max_pts:
        return time, flux
    idx = np.linspace(0, n - 1, max_pts, dtype=int)
    return time[idx], flux[idx]


def _transit_windows(time, period, t0):
    t_min, t_max = time.min(), time.max()
    if period <= 0:
        return []
    e0 = int(np.floor((t_min - t0) / period))
    e1 = int(np.ceil((t_max - t0) / period)) + 1
    return [(ep, t0 + ep * period) for ep in range(e0, e1 + 1) if t_min <= t0 + ep * period <= t_max]


def _pick_transit_mid_time(t, f, tc_calc, half_window):
    """
    Temps du minimum de flux autour de tc_calc, en évitant les creux parasites :
    on ne garde que les points du quartile inférieur du flux dans la fenêtre,
    puis le minimum parmi eux ; sinon minimum global dans la fenêtre.

    Retourne None s'il n'existe aucune mesure (temps + flux finis) dans la fenêtre
    (trou dans la courbe, pas de données) — pas de marqueur / mid-time dans ce cas.
    """
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    mask = np.abs(t - tc_calc) <= half_window
    if not np.any(mask):
        return None
    lt = t[mask]
    lf = f[mask]
    ok = np.isfinite(lt) & np.isfinite(lf)
    lt, lf = lt[ok], lf[ok]
    if len(lt) == 0:
        return None
    if len(lt) == 1:
        return float(lt[0])
    q25 = float(np.percentile(lf, 25))
    cand = lf <= q25
    if np.count_nonzero(cand) >= 2:
        idx = np.where(cand)[0]
        j = idx[int(np.argmin(lf[cand]))]
    else:
        j = int(np.argmin(lf))
    return float(lt[j])


class LCMarkersViewer(tk.Toplevel):
    """
    Fenêtre affichant la LC concaténée avec les périodes (marqueurs du périodogramme)
    reportées en couleurs. Ajustement manuel des mid-times par clic.
    Export : un mid-time.csv par période (couleur).
    """

    def __init__(self, parent, time, flux, periods):
        super().__init__(parent)
        self.title("Courbe de lumière – transits par période (marqueurs)")
        self.geometry("1100x750")

        time = np.asarray(time, dtype=float)
        flux = np.asarray(flux, dtype=float)
        self.time_full = time
        self.flux_full = flux
        self.time, self.flux = _downsample(time, flux)
        self.t0_display = float(np.min(time))

        # Liste de {period, t0, depth, duration, color, mid_times}
        self.periods_data = []
        for i, P in enumerate(periods):
            P = float(P)
            # Pas d'estimation de profondeur: depth fixé à 0.0
            t0 = float(np.nanmin(time)) if len(time) else 0.0
            depth = 0.0
            duration = max(0.01, P * 0.05)
            self.periods_data.append({
                'period': P,
                't0': t0,
                'depth': depth,
                'duration': duration,
                'color': PERIOD_COLORS[i % len(PERIOD_COLORS)],
                'mid_times': {},
            })

        self.selected_index = tk.IntVar(value=0)
        self.display_relative = tk.BooleanVar(value=True)
        self._build_ui()

    @staticmethod
    def _unpack_mid_time_value(value):
        """
        Normalise un enregistrement mid-time.
        Compatibilité:
          - ancien format: (tc, err)
          - nouveau format: (tc, err, source) avec source in {"auto","manual"}
        """
        if isinstance(value, (tuple, list)):
            if len(value) >= 3:
                return float(value[0]), float(value[1]), str(value[2])
            if len(value) == 2:
                return float(value[0]), float(value[1]), "auto"
        return float(value), 0.0, "auto"

    def _build_ui(self):
        main = ttk.Frame(self, padding=5)
        main.pack(fill=tk.BOTH, expand=True)

        toolbar = ttk.Frame(main)
        toolbar.pack(fill=tk.X, pady=(0, 5))
        ttk.Button(toolbar, text="✓ Valider tous les transits (créer tous les mid-times à l’épéméride)", command=self._fill_all_mid_times).pack(side=tk.LEFT, padx=2)
        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Button(toolbar, text="💾 Exporter sélection (manuel uniquement)", command=self._export_selected).pack(side=tk.LEFT, padx=2)
        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Button(toolbar, text="💾 Exporter mid-time.csv par période", command=self._export).pack(side=tk.LEFT, padx=2)
        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Button(toolbar, text="📐 Affiner Tc (fit trapèze)", command=self._open_transit_fit_dialog).pack(side=tk.LEFT, padx=2)
        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Checkbutton(
            toolbar,
            text="Abscisses : jours depuis le début",
            variable=self.display_relative,
            command=self._redraw,
        ).pack(side=tk.LEFT, padx=2)

        left = ttk.LabelFrame(main, text="Visualisation périodes", padding=6)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        ttk.Label(left, text="Périodes (sélection = couleur pour clic)", font=("", 9, "bold")).pack(anchor=tk.W)
        self.listbox = tk.Listbox(left, height=12, width=32, selectmode=tk.SINGLE, font=("Consolas", 9))
        self.listbox.pack(fill=tk.Y, pady=2)
        self.listbox.bind("<<ListboxSelect>>", self._on_select)
        ttk.Label(left, text="Largeur liste:", font=("", 8)).pack(anchor=tk.W, pady=(6, 0))
        self.slider_list_width = tk.Scale(
            left,
            from_=20,
            to=50,
            orient=tk.HORIZONTAL,
            length=180,
            command=self._on_slider_periods,
        )
        self.slider_list_width.set(32)
        self.slider_list_width.pack(fill=tk.X, pady=2)
        ttk.Label(left, text="Clic sur la courbe : fixe ou modifie le mid-time (période sélectionnée).", font=("", 8), foreground="gray", wraplength=240).pack(anchor=tk.W, pady=4)
        ttk.Label(
            left,
            text="Fit trapèze : validez d’abord des mid-times, puis affinez-les un par un (période sélectionnée).",
            font=("", 8),
            foreground="gray",
            wraplength=240,
        ).pack(anchor=tk.W, pady=(2, 0))

        for i, d in enumerate(self.periods_data):
            n = len(d.get('mid_times') or {})
            depth = d.get('depth') or 0
            self.listbox.insert(tk.END, f"  P={d['period']:.5f} j  prof.={depth:.4f}  ({n} TTV)")
        self.listbox.selection_set(0)

        fig_frame = ttk.Frame(main)
        fig_frame.pack(fill=tk.BOTH, expand=True)
        self.fig = plt.Figure(figsize=(9, 5), facecolor='white')
        self.ax = self.fig.add_subplot(111, facecolor='white')
        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.nav_toolbar = NavigationToolbar2Tk(self.canvas, fig_frame)
        self.nav_toolbar.update()
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self._redraw()

    def _on_select(self, event):
        sel = self.listbox.curselection()
        if sel:
            self.selected_index.set(sel[0])

    def _on_slider_periods(self, value):
        """Ajuste la largeur de la liste des périodes (nombre de caractères)."""
        try:
            w = int(float(value))
            self.listbox.configure(width=max(20, min(50, w)))
        except (ValueError, TypeError):
            pass

    def _refresh_listbox(self):
        """Met à jour la liste des périodes avec le nombre actuel de mid-times."""
        self.listbox.delete(0, tk.END)
        for i, d in enumerate(self.periods_data):
            n = len(d.get('mid_times') or {})
            depth = d.get('depth') or 0
            self.listbox.insert(tk.END, f"  P={d['period']:.5f} j  prof.={depth:.4f}  ({n} TTV)")
        self.listbox.selection_set(min(self.selected_index.get(), len(self.periods_data) - 1))

    def _fill_all_mid_times(self):
        """Recalcule tous les mid-times automatiquement via minimum local autour de chaque transit théorique."""
        if self.time_full is None or len(self.time_full) == 0:
            messagebox.showinfo("Valider", "Aucune courbe chargée.")
            return
        t = np.asarray(self.time_full, dtype=float)
        f = np.asarray(self.flux_full, dtype=float)
        if len(t) != len(f) or len(t) == 0:
            messagebox.showinfo("Valider", "Données LC invalides.")
            return
        total_added = 0
        for d in self.periods_data:
            P, T0 = d['period'], d['t0']
            if not np.isfinite(P) or not np.isfinite(T0):
                continue
            # Recalcule tous les points pour éviter de conserver d'anciens mid-times incohérents
            # (ex: issus d'un ancien mode de validation à l'épheméride).
            mid_times = {}
            # Fenêtre de recherche autour du transit théorique: 5% de P, bornée.
            half_window = max(0.01, min(P * 0.05, 0.4))
            for epoch, tc_calc in _transit_windows(self.time_full, P, T0):
                tc_obs = _pick_transit_mid_time(t, f, tc_calc, half_window)
                if tc_obs is not None:
                    mid_times[epoch] = (tc_obs, 0.0, "auto")
                    total_added += 1
            d['mid_times'] = mid_times
        self._refresh_listbox()
        self._redraw()
        messagebox.showinfo(
            "Valider",
            f"Mid-times auto-recalculés (fenêtres avec données uniquement). {total_added} point(s). "
            "Les époques sans mesure de flux dans la fenêtre sont ignorées. "
            "Vous pouvez encore ajuster un point en cliquant sur la courbe."
        )

    # Nombre max de bandes de transit dessinées par période (évite fond rouge avec ~1000 transits)
    MAX_DRAWN_TRANSITS = 200

    def _redraw(self):
        self.ax.clear()
        self.ax.set_facecolor('white')
        rel = self.display_relative.get()
        t0 = self.t0_display
        # Courbe en premier (toujours visible)
        if self.time is not None and self.flux is not None:
            x = (self.time - t0) if rel else self.time
            self.ax.plot(x, self.flux, 'k.', markersize=1.2, alpha=0.7)
        # Puis bandes et lignes (limité pour ne pas remplir tout l'écran)
        for i, d in enumerate(self.periods_data):
            P, T0 = d['period'], d['t0']
            if not np.isfinite(P) or not np.isfinite(T0):
                continue
            color = d['color']
            mid_times = d.get('mid_times') or {}
            n_drawn = 0
            # Uniquement les époques réellement validées (données de flux dans la fenêtre ou clic manuel)
            for epoch in sorted(mid_times.keys()):
                if n_drawn >= self.MAX_DRAWN_TRANSITS:
                    break
                tc_obs, _, _ = self._unpack_mid_time_value(mid_times[epoch])
                x_tc = (tc_obs - t0) if rel else tc_obs
                if not np.isfinite(x_tc):
                    continue
                half = max(P * 0.02, 0.01)
                self.ax.axvspan(x_tc - half, x_tc + half, alpha=0.25, color=color)
                self.ax.axvline(x_tc, color=color, lw=1.2, alpha=0.9)
                n_drawn += 1
        self.ax.set_ylabel("Flux (rel.)")
        if rel:
            self.ax.set_xlabel("Jours depuis le début")
            title = f"LC concaténée – t₀ = {t0:.2f} (une couleur = une période)"
        else:
            self.ax.set_xlabel("Temps (BJD-TDB ou JD)")
            title = "LC concaténée – une couleur = une période"
        self.ax.set_title(title)
        self.ax.grid(True, alpha=0.3)
        if self.time is not None:
            x_min = (self.time.min() - t0) if rel else self.time.min()
            x_max = (self.time.max() - t0) if rel else self.time.max()
            self.ax.set_xlim(x_min, x_max)
        self.ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=False))
        fmt = ScalarFormatter()
        fmt.set_scientific(False)
        fmt.set_powerlimits((-10, 10))
        self.ax.xaxis.set_major_formatter(fmt)
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def _on_click(self, event):
        # Si l'utilisateur est en mode zoom/pan, ne pas déclencher
        # l'ajout/modification de transit par clic.
        mode = getattr(self.nav_toolbar, 'mode', None) or getattr(self.nav_toolbar, '_active', None)
        if mode is not None and str(mode).strip().upper() not in ("", "NONE"):
            return

        if event.inaxes != self.ax or event.xdata is None or self.time is None:
            return
        rel = self.display_relative.get()
        t0 = self.t0_display
        t_click = event.xdata + t0 if rel else event.xdata
        idx = self.selected_index.get()
        if idx < 0 or idx >= len(self.periods_data):
            return
        d = self.periods_data[idx]
        P, T0 = d['period'], d['t0']
        epoch = int(round((t_click - T0) / P))
        d.setdefault('mid_times', {})[epoch] = (t_click, 0.0, "manual")
        self.listbox.delete(idx)
        n = len(d.get('mid_times') or {})
        depth = d.get('depth') or 0
        self.listbox.insert(idx, f"  P={d['period']:.5f} j  prof.={depth:.4f}  ({n} TTV)")
        self.listbox.selection_set(idx)
        self._redraw()

    def _open_transit_fit_dialog(self):
        """Ajustement local trapézoïdal pour affiner Tc (époque par époque)."""
        idx = self.selected_index.get()
        if idx < 0 or idx >= len(self.periods_data):
            messagebox.showinfo("Fit transit", "Sélectionnez une période dans la liste.")
            return
        d = self.periods_data[idx]
        mid_times = d.get("mid_times") or {}
        if not mid_times:
            messagebox.showinfo(
                "Fit transit",
                "Aucun mid-time pour cette période. Utilisez « Valider tous les transits » ou cliquez sur la courbe.",
            )
            return

        def on_apply(epoch, t0, err):
            d.setdefault("mid_times", {})[int(epoch)] = (float(t0), float(err), "fit")
            self._refresh_listbox()
            self._redraw()

        TransitFitDialog(
            self,
            self.time_full,
            self.flux_full,
            None,
            d,
            on_apply,
        )

    def _export(self):
        out_dir = filedialog.askdirectory(title="Dossier pour les fichiers mid-time.csv")
        if not out_dir:
            return
        out_path = Path(out_dir)
        exported = []
        for i, d in enumerate(self.periods_data):
            mid_times = d.get('mid_times') or {}
            if not mid_times:
                continue
            P = d['period']
            depth = d.get('depth') or 0
            label = f"P{P:.4f}j".replace('.', '_')
            rows = []
            for ep, raw in sorted(mid_times.items()):
                tc, err, _ = self._unpack_mid_time_value(raw)
                rows.append({"Epoch": ep, "Mid_time": tc, "Uncertainty": err, "Depth": depth})
            df = pd.DataFrame(rows)
            # O-C demandé: T0 = premier mid-time exporté, P = période retenue (marqueur périodogramme, sans ajustement).
            ep = df["Epoch"].values.astype(float)
            mt = df["Mid_time"].values.astype(float)
            P_ref = float(d["period"])
            if np.isfinite(P_ref) and len(df) >= 1 and np.all(np.isfinite(ep)) and np.all(np.isfinite(mt)):
                ep0 = ep[0]
                T0_ref = mt[0]
                tc_linear = T0_ref + (ep - ep0) * P_ref
                oc_days = mt - tc_linear
                df["O-C_minutes"] = oc_days * (24 * 60)
            else:
                df["O-C_minutes"] = 0.0
            df = df[["Epoch", "Mid_time", "O-C_minutes", "Uncertainty", "Depth"]]
            fname = f"mid-time_{label}.csv"
            filepath = out_path / fname
            df.to_csv(filepath, index=False)
            exported.append(fname)
        if exported:
            messagebox.showinfo("Export", "Fichiers créés :\n" + "\n".join(exported))
        else:
            messagebox.showinfo("Export", "Aucun mid-time ajusté. Cliquez sur la courbe (période sélectionnée) pour en ajouter.")

    def _export_selected(self):
        """Exporte uniquement les mid-times MANUELS de la période sélectionnée."""
        idx = self.selected_index.get()
        if idx < 0 or idx >= len(self.periods_data):
            messagebox.showinfo("Export sélection", "Aucune période sélectionnée.")
            return

        d = self.periods_data[idx]
        mid_times = d.get('mid_times') or {}
        if not mid_times:
            messagebox.showinfo(
                "Export sélection",
                "La période sélectionnée ne contient aucun mid-time.\n"
                "Ajoutez des points par clic ou utilisez « Valider tous les transits »."
            )
            return
        manual_rows = []
        for ep, raw in sorted(mid_times.items()):
            tc, err, source = self._unpack_mid_time_value(raw)
            if str(source).lower() == "manual":
                manual_rows.append((ep, tc, err))
        if not manual_rows:
            messagebox.showinfo(
                "Export sélection",
                "Aucun mid-time manuel dans cette période.\n"
                "Cliquez sur la courbe pour créer des points manuels."
            )
            return

        out_dir = filedialog.askdirectory(title="Dossier pour le mid-time.csv de la période sélectionnée")
        if not out_dir:
            return
        out_path = Path(out_dir)

        P = d['period']
        depth = d.get('depth') or 0
        label = f"P{P:.4f}j".replace('.', '_')
        rows = [{"Epoch": ep, "Mid_time": tc, "Uncertainty": err, "Depth": depth} for ep, tc, err in manual_rows]
        df = pd.DataFrame(rows)

        ep = df["Epoch"].values.astype(float)
        mt = df["Mid_time"].values.astype(float)
        P_ref = float(d["period"])
        if np.isfinite(P_ref) and len(df) >= 1 and np.all(np.isfinite(ep)) and np.all(np.isfinite(mt)):
            ep0 = ep[0]
            T0_ref = mt[0]
            tc_linear = T0_ref + (ep - ep0) * P_ref
            oc_days = mt - tc_linear
            df["O-C_minutes"] = oc_days * (24 * 60)
        else:
            df["O-C_minutes"] = 0.0

        df = df[["Epoch", "Mid_time", "O-C_minutes", "Uncertainty", "Depth"]]
        fname = f"mid-time_{label}.csv"
        filepath = out_path / fname
        df.to_csv(filepath, index=False)
        messagebox.showinfo(
            "Export sélection",
            f"Fichier créé (mid-times manuels uniquement) :\n{filepath}\n"
            f"Lignes exportées : {len(df)}"
        )
