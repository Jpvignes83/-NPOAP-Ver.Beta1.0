# gui/transit_fit_dialog.py
"""Fenêtre d'ajustement trapézoïdal local pour affiner le mid-time d'un transit.

Option Wotan : détrendage du segment avec masque transit (Tc, T14) avant le fit.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import logging

import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from core.transit_fitter import fit_trapezoid_transit, model_flux, slice_window
from core.wotan_detrend import (
    default_window_length,
    flatten_segment_with_transit_mask,
    is_wotan_available,
    transit_boolean_mask,
)

logger = logging.getLogger(__name__)

# Méthodes Wotan sans dépendances optionnelles lourdes (voir doc wotan).
WOTAN_METHOD_CHOICES = ("biweight", "median", "mean", "huber", "lowess")


class TransitFitDialog(tk.Toplevel):
    """
    Sélection d'une époque parmi les mid-times déjà présents, découpe locale,
    fit trapèze + baseline, application du Tc raffiné.
    """

    def __init__(self, parent, time_full, flux_full, flux_err_full, period_dict, on_apply):
        super().__init__(parent)
        self.title("Affiner le mid-time (fit trapèze)")
        self.geometry("920x700")
        self._on_apply = on_apply
        self._d = period_dict
        self._time = np.asarray(time_full, dtype=float)
        self._flux = np.asarray(flux_full, dtype=float)
        self._ferr = None
        if flux_err_full is not None:
            fe = np.asarray(flux_err_full, dtype=float)
            if len(fe) == len(self._time):
                self._ferr = fe

        self._last_fit = None
        self._last_used_wotan = False
        self._wotan_ok = is_wotan_available()

        main = ttk.Frame(self, padding=8)
        main.pack(fill=tk.BOTH, expand=True)

        top = ttk.Frame(main)
        top.pack(fill=tk.X, pady=(0, 6))

        ttk.Label(top, text="Époque :").pack(side=tk.LEFT, padx=(0, 4))
        mid = period_dict.get("mid_times") or {}
        self._epoch_keys = sorted(mid.keys())
        labels = [f"Époque {ep}  (Tc ≈ {self._unpack_tc(mid[ep]):.6f})" for ep in self._epoch_keys]
        self._epoch_var = tk.StringVar()
        self._combo = ttk.Combobox(top, textvariable=self._epoch_var, values=labels, width=52, state="readonly")
        if labels:
            self._combo.current(0)
        self._combo.pack(side=tk.LEFT, padx=4)
        self._combo.bind("<<ComboboxSelected>>", lambda e: self._refresh_plot_data_only())

        ttk.Label(top, text="Demi-fenêtre (j)").pack(side=tk.LEFT, padx=(12, 2))
        self._half_win_var = tk.StringVar(value="")
        ttk.Entry(top, textvariable=self._half_win_var, width=8).pack(side=tk.LEFT, padx=2)

        ttk.Label(top, text="Largeur transit T14 (j)").pack(side=tk.LEFT, padx=(8, 2))
        P = float(period_dict.get("period") or 0.1)
        dur = float(period_dict.get("duration") or max(0.02, P * 0.05))
        self._width_var = tk.StringVar(value=f"{dur:.6f}")
        ttk.Entry(top, textvariable=self._width_var, width=10).pack(side=tk.LEFT, padx=2)

        mid_row = ttk.Frame(main)
        mid_row.pack(fill=tk.X, pady=4)
        ttk.Label(mid_row, text="Ingress (fraction demi-largeur)").pack(side=tk.LEFT, padx=(0, 4))
        self._ing_var = tk.StringVar(value="0.25")
        ttk.Entry(mid_row, textvariable=self._ing_var, width=6).pack(side=tk.LEFT, padx=2)
        self._slope_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(mid_row, text="Ajuster pente baseline", variable=self._slope_var).pack(
            side=tk.LEFT, padx=16
        )

        wotan_fr = ttk.LabelFrame(main, text="Détrendage Wotan (optionnel)", padding=6)
        wotan_fr.pack(fill=tk.X, pady=(2, 4))
        self._wotan_var = tk.BooleanVar(value=False)
        self._wotan_cb = ttk.Checkbutton(
            wotan_fr,
            text="Appliquer Wotan avant le fit (masque transit = Tc ± 1.15 × T14/2)",
            variable=self._wotan_var,
        )
        self._wotan_cb.pack(anchor=tk.W)
        if not self._wotan_ok:
            self._wotan_cb.state(["disabled"])
            ttk.Label(
                wotan_fr,
                text="Module « wotan » absent : pip install wotan  (voir requirements.txt)",
                foreground="gray",
                font=("", 8),
            ).pack(anchor=tk.W, pady=(0, 2))

        wotan_params = ttk.Frame(wotan_fr)
        wotan_params.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(wotan_params, text="Fenêtre Wotan (j), vide = auto").pack(side=tk.LEFT, padx=(0, 4))
        self._wotan_wl_var = tk.StringVar(value="")
        self._wotan_wl_entry = ttk.Entry(wotan_params, textvariable=self._wotan_wl_var, width=10)
        self._wotan_wl_entry.pack(side=tk.LEFT, padx=2)
        if not self._wotan_ok:
            self._wotan_wl_entry.state(["disabled"])

        ttk.Label(wotan_params, text="Méthode").pack(side=tk.LEFT, padx=(14, 4))
        self._wotan_method_var = tk.StringVar(value="biweight")
        self._wotan_method_combo = ttk.Combobox(
            wotan_params,
            textvariable=self._wotan_method_var,
            values=WOTAN_METHOD_CHOICES,
            width=12,
            state="readonly",
        )
        self._wotan_method_combo.pack(side=tk.LEFT, padx=2)
        if not self._wotan_ok:
            self._wotan_method_combo.state(["disabled"])

        btn_row = ttk.Frame(main)
        btn_row.pack(fill=tk.X, pady=4)
        ttk.Button(btn_row, text="Lancer le fit", command=self._run_fit).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row, text="Appliquer Tc à cette époque", command=self._apply).pack(side=tk.LEFT, padx=8)
        ttk.Button(btn_row, text="Fermer", command=self.destroy).pack(side=tk.RIGHT, padx=2)

        self._status = tk.StringVar(value="Choisissez une époque puis « Lancer le fit ».")
        ttk.Label(main, textvariable=self._status, font=("", 9)).pack(anchor=tk.W, pady=2)

        fig_fr = ttk.Frame(main)
        fig_fr.pack(fill=tk.BOTH, expand=True)
        self.fig = plt.Figure(figsize=(8.5, 4.2), facecolor="white")
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_fr)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._nav = NavigationToolbar2Tk(self.canvas, fig_fr)
        self._nav.update()

        self._init_defaults_half_window()
        self._refresh_plot_data_only()

    @staticmethod
    def _unpack_tc(raw):
        if isinstance(raw, (tuple, list)):
            return float(raw[0])
        return float(raw)

    def _selected_epoch(self):
        i = self._combo.current()
        if i < 0 or i >= len(self._epoch_keys):
            return None, None
        ep = self._epoch_keys[i]
        tc = self._unpack_tc((self._d.get("mid_times") or {}).get(ep))
        return ep, tc

    def _init_defaults_half_window(self):
        P = float(self._d.get("period") or 0.1)
        hw = max(0.02, min(P * 0.08, 0.6))
        self._half_win_var.set(f"{hw:.6f}")

    def _refresh_plot_data_only(self):
        ep, tc = self._selected_epoch()
        self.ax.clear()
        if tc is None:
            self.ax.text(0.5, 0.5, "Aucun mid-time pour cette période.", ha="center", va="center", transform=self.ax.transAxes)
            self.canvas.draw_idle()
            return
        try:
            hw = float(self._half_win_var.get().replace(",", "."))
        except ValueError:
            hw = max(0.02, float(self._d.get("period") or 0.1) * 0.08)
        t, f, e = slice_window(self._time, self._flux, self._ferr, tc, hw)
        self.ax.plot(t, f, "k.", markersize=3, alpha=0.75, label="Données")
        self.ax.axvline(tc, color="C1", ls="--", lw=1, label="Tc initial")
        self.ax.set_xlabel("Temps")
        self.ax.set_ylabel("Flux")
        self.ax.set_title(f"Fenêtre locale — époque {ep}")
        self.ax.legend(loc="best", fontsize=8)
        self.ax.grid(True, alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self._last_fit = None

    def _run_fit(self):
        ep, tc = self._selected_epoch()
        if tc is None:
            messagebox.showinfo("Fit", "Aucune époque valide.")
            return
        try:
            hw = float(self._half_win_var.get().replace(",", "."))
            width = float(self._width_var.get().replace(",", "."))
            ing = float(self._ing_var.get().replace(",", "."))
        except ValueError:
            messagebox.showerror("Fit", "Valeurs numériques invalides (fenêtre, largeur, ingress).")
            return

        t, f, e = slice_window(self._time, self._flux, self._ferr, tc, hw)
        if len(t) < 8:
            messagebox.showwarning("Fit", "Pas assez de points dans la fenêtre (élargir demi-fenêtre).")
            return

        self._last_used_wotan = False
        f_fit = f
        e_fit = e
        wotan_note = ""

        if self._wotan_var.get():
            if not self._wotan_ok:
                messagebox.showerror(
                    "Wotan",
                    "Le module wotan n'est pas installé.\n"
                    "Installez-le avec : pip install wotan",
                )
                return
            wl_user = None
            s = (self._wotan_wl_var.get() or "").strip().replace(",", ".")
            if s:
                try:
                    wl_user = float(s)
                    if wl_user <= 0:
                        raise ValueError("window_length doit être > 0")
                except ValueError:
                    messagebox.showerror("Wotan", "Fenêtre Wotan (j) invalide (nombre > 0 ou laisser vide).")
                    return
            method = (self._wotan_method_var.get() or "biweight").strip().lower()
            if method not in WOTAN_METHOD_CHOICES:
                method = "biweight"

            try:
                f_fit, _trend, _mask, scale = flatten_segment_with_transit_mask(
                    t,
                    f,
                    tc,
                    width,
                    window_length=wl_user,
                    method=method,
                    break_tolerance=0.0,
                    edge_cutoff=0.0,
                    mask_edge_factor=1.15,
                )
            except Exception as ex:
                logger.exception("Wotan flatten")
                messagebox.showerror("Wotan", f"Échec du détrendage Wotan :\n{ex}")
                return

            self._last_used_wotan = True
            if e is not None:
                e_fit = e / max(scale, 1e-15)
            wl_eff = wl_user if wl_user is not None else default_window_length(t)
            wotan_note = f" | Wotan ({method}, fen.≈{wl_eff:.4f} j)"

        out = fit_trapezoid_transit(
            t,
            f_fit,
            e_fit,
            t0_init=tc,
            half_window=hw,
            width_fixed=width,
            ingress_ratio=ing,
            fit_slope=self._slope_var.get(),
        )
        self._last_fit = out
        if not out.get("success"):
            self._status.set(out.get("message", "Échec"))
            messagebox.showwarning("Fit", out.get("message", "Échec du fit"))
            self._refresh_plot_data_only()
            return

        t_ref = out["t_ref"]
        tm = model_flux(
            t,
            out["t0"],
            out["depth"],
            out["a0"],
            out["a1"],
            out["width"],
            out["ingress_ratio"],
            t_ref,
        )
        self.ax.clear()
        if self._last_used_wotan:
            mvis = transit_boolean_mask(t, tc, width, edge_factor=1.15)
            if np.any(mvis):
                t_lo = float(np.min(t[mvis]))
                t_hi = float(np.max(t[mvis]))
                self.ax.axvspan(t_lo, t_hi, alpha=0.12, color="0.5", label="Masque transit (Wotan)")
        self.ax.plot(t, f_fit, "k.", markersize=3, alpha=0.75, label="Données (fit)" + (" — Wotan+norm." if self._last_used_wotan else ""))
        ti = np.linspace(float(t.min()), float(t.max()), max(200, len(t) * 3))
        tm_fine = model_flux(
            ti,
            out["t0"],
            out["depth"],
            out["a0"],
            out["a1"],
            out["width"],
            out["ingress_ratio"],
            t_ref,
        )
        self.ax.plot(ti, tm_fine, "r-", lw=1.5, label="Modèle")
        self.ax.axvline(tc, color="C1", ls="--", lw=1, alpha=0.7, label="Tc initial")
        self.ax.axvline(out["t0"], color="C2", ls="-", lw=1.2, label="Tc fit")
        self.ax.set_xlabel("Temps")
        self.ax.set_ylabel("Flux")
        self.ax.set_title(f"Époque {ep} — Tc = {out['t0']:.6f} ± {out['t0_err']:.6f} j")
        self.ax.legend(loc="best", fontsize=8)
        self.ax.grid(True, alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self._status.set(
            f"Tc = {out['t0']:.6f} ± {out['t0_err']:.6f} j | RMS rés. = {out['rms']:.6f} | prof. ≈ {out['depth']:.5f}"
            + wotan_note
        )

    def _apply(self):
        if not self._last_fit or not self._last_fit.get("success"):
            messagebox.showinfo("Appliquer", "Lancez d'abord un fit réussi.")
            return
        ep, _ = self._selected_epoch()
        if ep is None:
            return
        t0 = self._last_fit["t0"]
        err = float(self._last_fit.get("t0_err") or 0.0)
        try:
            self._on_apply(ep, t0, err)
        except Exception as e:
            logger.exception("on_apply transit fit")
            messagebox.showerror("Erreur", str(e))
            return
        messagebox.showinfo("Appliquer", f"Époque {ep} : Mid_time mis à jour ({t0:.6f} j, σ≈{err:.6f}).")
        self.destroy()
