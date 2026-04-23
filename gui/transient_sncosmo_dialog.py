# gui/transient_sncosmo_dialog.py
"""Boîte de dialogue : ajustement SN Ia (SALT2 / SALT3) avec sncosmo, filtres Gaia G / G_Bp / G_Rp."""

from __future__ import annotations

import logging
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from threading import Thread

logger = logging.getLogger(__name__)


class TransientSncosmoDialog(tk.Toplevel):
    """
    Charge un CSV de courbe de lumière (MJD, filtre Gaia, mag, erreur) et lance ``sncosmo.fit_lc``.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.title("SN Ia — sncosmo (SALT2 / SALT3, filtres Gaia)")
        self.geometry("780x620")
        self.minsize(640, 480)

        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        info = (
            "CSV avec en-tête : colonnes mjd, band, mag, mag_err (noms modifiables ci-dessous).\n"
            "Colonne « band » : G, G_Bp ou G_Rp (équivalent Gaia DR3 ; mappées vers gaia::g, gaia::gbp, gaia::grp dans sncosmo).\n"
            "Magnitudes et zéro-point doivent être cohérents avec le système photométrique choisi (AB ou Vega)."
        )
        ttk.Label(main, text=info, wraplength=740, justify=tk.LEFT).pack(anchor=tk.W, pady=(0, 8))

        file_fr = ttk.LabelFrame(main, text="Fichier", padding=6)
        file_fr.pack(fill=tk.X, pady=4)
        self._csv_var = tk.StringVar()
        ttk.Entry(file_fr, textvariable=self._csv_var, width=70).pack(side=tk.LEFT, padx=(0, 6), fill=tk.X, expand=True)
        ttk.Button(file_fr, text="Parcourir…", command=self._browse_csv).pack(side=tk.LEFT)

        cols = ttk.LabelFrame(main, text="Colonnes CSV", padding=6)
        cols.pack(fill=tk.X, pady=4)
        self._col_mjd = tk.StringVar(value="mjd")
        self._col_band = tk.StringVar(value="band")
        self._col_mag = tk.StringVar(value="mag")
        self._col_merr = tk.StringVar(value="mag_err")
        self._delimiter = tk.StringVar(value=",")
        for i, (lab, var) in enumerate(
            [
                ("Temps (MJD)", self._col_mjd),
                ("Filtre (G / G_Bp / G_Rp)", self._col_band),
                ("Magnitude", self._col_mag),
                ("Erreur mag.", self._col_merr),
            ]
        ):
            ttk.Label(cols, text=lab + " :").grid(row=i // 2, column=(i % 2) * 2, sticky=tk.E, padx=4, pady=2)
            ttk.Entry(cols, textvariable=var, width=14).grid(row=i // 2, column=(i % 2) * 2 + 1, sticky=tk.W, padx=4, pady=2)
        ttk.Label(cols, text="Séparateur :").grid(row=2, column=0, sticky=tk.E, padx=4, pady=2)
        ttk.Entry(cols, textvariable=self._delimiter, width=4).grid(row=2, column=1, sticky=tk.W, padx=4, pady=2)

        phot = ttk.LabelFrame(main, text="Photométrie (sncosmo)", padding=6)
        phot.pack(fill=tk.X, pady=4)
        ttk.Label(phot, text="Zéro-point zp :").grid(row=0, column=0, sticky=tk.E, padx=4)
        self._zp = tk.StringVar(value="25.0")
        ttk.Entry(phot, textvariable=self._zp, width=10).grid(row=0, column=1, sticky=tk.W)
        ttk.Label(phot, text="Système :").grid(row=0, column=2, sticky=tk.E, padx=(16, 4))
        self._zpsys = tk.StringVar(value="ab")
        ttk.Combobox(phot, textvariable=self._zpsys, values=("ab", "vega"), width=8, state="readonly").grid(
            row=0, column=3, sticky=tk.W
        )

        model_fr = ttk.LabelFrame(main, text="Modèle et redshift", padding=6)
        model_fr.pack(fill=tk.X, pady=4)
        ttk.Label(model_fr, text="Source :").grid(row=0, column=0, sticky=tk.E, padx=4)
        self._source = tk.StringVar(value="salt2")
        ttk.Combobox(
            model_fr,
            textvariable=self._source,
            values=("salt2", "salt3"),
            width=12,
            state="readonly",
        ).grid(row=0, column=1, sticky=tk.W)
        self._fit_z = tk.BooleanVar(value=False)
        ttk.Checkbutton(model_fr, text="Ajuster le redshift z", variable=self._fit_z).grid(
            row=0, column=2, padx=(20, 4)
        )
        ttk.Label(model_fr, text="z fixe :").grid(row=1, column=0, sticky=tk.E, padx=4, pady=4)
        self._z_fixed = tk.StringVar(value="0.05")
        ttk.Entry(model_fr, textvariable=self._z_fixed, width=10).grid(row=1, column=1, sticky=tk.W, pady=4)
        ttk.Label(model_fr, text="Bornes z (min, max) si ajustement :").grid(row=1, column=2, sticky=tk.E, padx=(12, 4))
        self._zmin = tk.StringVar(value="0.01")
        self._zmax = tk.StringVar(value="0.5")
        zf = ttk.Frame(model_fr)
        zf.grid(row=1, column=3, sticky=tk.W)
        ttk.Entry(zf, textvariable=self._zmin, width=8).pack(side=tk.LEFT)
        ttk.Label(zf, text=" — ").pack(side=tk.LEFT)
        ttk.Entry(zf, textvariable=self._zmax, width=8).pack(side=tk.LEFT)

        mw_fr = ttk.LabelFrame(main, text="Extinction MW (ligne de visée)", padding=6)
        mw_fr.pack(fill=tk.X, pady=4)
        self._apply_mw = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            mw_fr,
            text="Corriger extinction Galactique avant fit",
            variable=self._apply_mw,
        ).grid(row=0, column=0, columnspan=4, sticky=tk.W, padx=4, pady=(0, 4))
        ttk.Label(mw_fr, text="RA (deg) :").grid(row=1, column=0, sticky=tk.E, padx=4, pady=2)
        self._ra_deg = tk.StringVar(value="")
        ttk.Entry(mw_fr, textvariable=self._ra_deg, width=12).grid(row=1, column=1, sticky=tk.W, padx=4, pady=2)
        ttk.Label(mw_fr, text="Dec (deg) :").grid(row=1, column=2, sticky=tk.E, padx=4, pady=2)
        self._dec_deg = tk.StringVar(value="")
        ttk.Entry(mw_fr, textvariable=self._dec_deg, width=12).grid(row=1, column=3, sticky=tk.W, padx=4, pady=2)
        ttk.Label(mw_fr, text="E(B-V) MW :").grid(row=2, column=0, sticky=tk.E, padx=4, pady=2)
        self._ebv_mw = tk.StringVar(value="")
        ttk.Entry(mw_fr, textvariable=self._ebv_mw, width=12).grid(row=2, column=1, sticky=tk.W, padx=4, pady=2)
        ttk.Button(mw_fr, text="Récupérer E(B-V) (SFD)", command=self._fetch_mw_ebv).grid(
            row=2, column=2, columnspan=2, sticky=tk.W, padx=4, pady=2
        )

        btn_fr = ttk.Frame(main)
        btn_fr.pack(fill=tk.X, pady=8)
        self._run_btn = ttk.Button(btn_fr, text="Lancer l’ajustement", command=self._run_fit)
        self._run_btn.pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(btn_fr, text="Fermer", command=self.destroy).pack(side=tk.LEFT)

        out_lab = ttk.Label(main, text="Résultat")
        out_lab.pack(anchor=tk.W)
        text_fr = ttk.Frame(main)
        text_fr.pack(fill=tk.BOTH, expand=True, pady=(4, 0))
        scroll = ttk.Scrollbar(text_fr)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self._text = tk.Text(
            text_fr, height=16, wrap=tk.WORD, font=("Consolas", 9), yscrollcommand=scroll.set
        )
        self._text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.config(command=self._text.yview)

        self._append_help()

    def _append_help(self) -> None:
        self._text.insert(
            tk.END,
            "sncosmo est installé via requirements.txt (installation standard NPOAP).\n\n",
        )
        self._text.config(state=tk.DISABLED)

    def _browse_csv(self) -> None:
        path = filedialog.askopenfilename(
            parent=self,
            title="CSV courbe de lumière",
            filetypes=[("CSV", "*.csv"), ("Texte", "*.txt"), ("Tous", "*.*")],
        )
        if path:
            self._csv_var.set(path)

    def _run_fit(self) -> None:
        from core.sncosmo_sn_ia import (
            SNCOSMO_AVAILABLE,
            SNCOSMO_IMPORT_ERROR,
            apply_mw_extinction_correction_gaia,
            build_photometric_table,
            load_lc_csv,
            run_salt_fit,
        )

        if not SNCOSMO_AVAILABLE:
            messagebox.showerror(
                "sncosmo absent",
                "Le module sncosmo n’est pas installé.\n\n"
                f"Détail : {SNCOSMO_IMPORT_ERROR or 'import impossible'}\n\n"
                "Installez le profil optionnel :\n"
                "  pip install -r requirements-cosmology-sne.txt",
                parent=self,
            )
            return

        path = self._csv_var.get().strip()
        if not path:
            messagebox.showwarning("Fichier", "Choisissez un fichier CSV.", parent=self)
            return

        try:
            zp = float(self._zp.get().replace(",", "."))
        except ValueError:
            messagebox.showerror("Paramètre", "Zéro-point zp invalide.", parent=self)
            return

        delim = self._delimiter.get().strip() or ","
        try:
            t, bands, mag, merr = load_lc_csv(
                path,
                col_time=self._col_mjd.get().strip(),
                col_band=self._col_band.get().strip(),
                col_mag=self._col_mag.get().strip(),
                col_mag_err=self._col_merr.get().strip(),
                delimiter=delim,
            )
        except Exception as e:
            messagebox.showerror("CSV", str(e), parent=self)
            return

        try:
            phot = build_photometric_table(t, bands, mag, merr, zp=zp, zpsys=self._zpsys.get().strip())
        except Exception as e:
            messagebox.showerror("Données", str(e), parent=self)
            return

        correction_note = None
        if self._apply_mw.get():
            ebv_txt = self._ebv_mw.get().strip()
            if not ebv_txt:
                messagebox.showerror(
                    "Extinction MW",
                    "Renseignez E(B-V) MW ou utilisez le bouton « Récupérer E(B-V) (SFD) ».",
                    parent=self,
                )
                return
            try:
                ebv_val = float(ebv_txt.replace(",", "."))
            except ValueError:
                messagebox.showerror("Extinction MW", "Valeur E(B-V) invalide.", parent=self)
                return
            mag_corr = apply_mw_extinction_correction_gaia(mag, bands, ebv_val)
            try:
                phot = build_photometric_table(t, bands, mag_corr, merr, zp=zp, zpsys=self._zpsys.get().strip())
            except Exception as e:
                messagebox.showerror("Extinction MW", f"Correction impossible : {e}", parent=self)
                return
            correction_note = f"Correction MW appliquée: E(B-V)={ebv_val:.4f} (Gaia G/BP/RP)"

        fit_z = self._fit_z.get()
        z_fixed: float | None = None
        z_bounds = (0.001, 1.2)
        if fit_z:
            try:
                zlo = float(self._zmin.get().replace(",", "."))
                zhi = float(self._zmax.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Redshift", "Bornes z invalides.", parent=self)
                return
            if zlo >= zhi:
                messagebox.showerror("Redshift", "z_min doit être < z_max.", parent=self)
                return
            z_bounds = (zlo, zhi)
            try:
                z_fixed = float(self._z_fixed.get().replace(",", "."))
            except ValueError:
                z_fixed = None
        else:
            try:
                z_fixed = float(self._z_fixed.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Redshift", "Valeur de z fixe invalide.", parent=self)
                return

        source = self._source.get().strip()
        self._run_btn.state(["disabled"])
        self._text.config(state=tk.NORMAL)
        self._text.delete("1.0", tk.END)
        self._text.insert(tk.END, "Ajustement en cours…\n")
        self._text.config(state=tk.DISABLED)

        def task():
            try:
                res = run_salt_fit(
                    phot,
                    source=source,
                    z_fixed=z_fixed,
                    fit_z=fit_z,
                    z_bounds=z_bounds,
                )
                if correction_note:
                    res["mw_correction"] = correction_note
            except Exception as e:
                logger.exception("Fit sncosmo")
                res = {"ok": False, "error": str(e), "source": source}

            def done():
                self._run_btn.state(["!disabled"])
                self._text.config(state=tk.NORMAL)
                self._text.delete("1.0", tk.END)
                self._format_result(res)
                self._text.config(state=tk.DISABLED)

            self.after(0, done)

        Thread(target=task, daemon=True).start()

    def _format_result(self, res: dict) -> None:
        lines: list[str] = []
        if res.get("error"):
            lines.append("Erreur : " + str(res["error"]))
            lines.append("")
        lines.append(f"Modèle : {res.get('source', '')}")
        lines.append(f"Succès optimisation : {res.get('ok', False)}")
        if res.get("message"):
            lines.append(f"Message : {res['message']}")
        if res.get("mw_correction"):
            lines.append(str(res["mw_correction"]))
        dropped = res.get("dropped_bands") or []
        if dropped:
            lines.append("Bandes rejetées (hors domaine spectral du modèle) : " + ", ".join(dropped))
        chisq = res.get("chisq")
        ndof_val = res.get("ndof", -1)
        try:
            cq = float(chisq)
            if cq == cq:
                lines.append(f"chi² = {cq:.4f}   ndof = {ndof_val}")
        except (TypeError, ValueError):
            pass
        parms = res.get("parameters") or {}

        # Drapeau de fiabilité du fit (heuristiques pragmatiques)
        reliability_issues: list[str] = []
        try:
            if int(ndof_val) < 5:
                reliability_issues.append(f"ndof faible ({ndof_val} < 5)")
        except (TypeError, ValueError):
            reliability_issues.append("ndof invalide")
        c_val = parms.get("c")
        if c_val is not None:
            try:
                if abs(float(c_val)) > 1.0:
                    reliability_issues.append(f"|c| élevé ({float(c_val):.3g} > 1)")
            except (TypeError, ValueError):
                pass
        msg = str(res.get("message", "") or "").lower()
        if "covariance" in msg and ("not positive definite" in msg or "may not be accurate" in msg):
            reliability_issues.append("covariance potentiellement non fiable")
        if dropped:
            reliability_issues.append(f"{len(dropped)} bande(s) rejetée(s)")

        if reliability_issues:
            lines.append("")
            lines.append("[ALERTE] FIT NON FIABLE")
            lines.append("Raisons : " + "; ".join(reliability_issues))

        if parms:
            lines.append("")
            lines.append("Paramètres :")
            for k, v in parms.items():
                e = (res.get("errors") or {}).get(k)
                if e is not None:
                    lines.append(f"  {k} = {float(v):.6g} ± {float(e):.6g}")
                else:
                    lines.append(f"  {k} = {float(v):.6g}")
        self._text.insert(tk.END, "\n".join(lines) + "\n")

    def _fetch_mw_ebv(self) -> None:
        from core.sncosmo_sn_ia import get_mw_ebv_sfd

        try:
            ra = float(self._ra_deg.get().strip().replace(",", "."))
            dec = float(self._dec_deg.get().strip().replace(",", "."))
        except ValueError:
            messagebox.showerror("Coordonnées", "RA/Dec invalides (en degrés).", parent=self)
            return
        try:
            ebv = get_mw_ebv_sfd(ra, dec)
            self._ebv_mw.set(f"{ebv:.4f}")
            messagebox.showinfo("Extinction MW", f"E(B-V) (SFD) = {ebv:.4f}", parent=self)
        except Exception as e:
            logger.exception("Récupération E(B-V) MW")
            messagebox.showerror("Extinction MW", f"Échec récupération E(B-V):\n{e}", parent=self)


def open_transient_sncosmo_dialog(parent) -> None:
    """Ouvre la boîte de dialogue (Toplevel liée à la fenêtre parente)."""
    dlg = TransientSncosmoDialog(parent)
    try:
        dlg.transient(parent)
    except tk.TclError:
        pass
