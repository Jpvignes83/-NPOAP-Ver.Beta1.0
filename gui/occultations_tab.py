# gui/occultations_tab.py
"""Onglet Occultations : SORA (prédictions, analyse), PyMovie."""

import json
import logging
import subprocess
import threading
import tkinter as tk
from datetime import timedelta, timezone
from pathlib import Path
from typing import Any, Optional
from tkinter import filedialog, messagebox, ttk
import webbrowser

try:
    from zoneinfo import ZoneInfo
except ImportError:
    ZoneInfo = None  # type: ignore[misc, assignment]

import config

logger = logging.getLogger(__name__)

from core.sora_ephemeris import (
    format_prediction_table,
    run_occultation_prediction,
    try_import_sora,
)


def _open_url(url: str) -> None:
    try:
        webbrowser.open(url)
    except Exception:
        pass


def _default_pymovie_conda_root() -> str:
    """Racine type de l’environnement dédié PyMovie (conda user, Windows)."""
    return str(Path.home() / ".conda" / "envs" / "pymovie")


def _python_exe_in_conda_env(root: Path) -> Optional[Path]:
    if not root.is_dir():
        return None
    win_py = root / "python.exe"
    if win_py.is_file():
        return win_py
    nix_py = root / "bin" / "python"
    if nix_py.is_file():
        return nix_py
    return None


def _observatory_zoneinfo():
    """Fuseau IANA pour conversion UTC → heure locale (repli UTC)."""
    if ZoneInfo is None:
        return timezone.utc, "UTC (zoneinfo indisponible)"
    o = getattr(config, "OBSERVATORY", None) or {}
    raw = str(o.get("timezone", "") or "").strip()
    if not raw:
        return ZoneInfo("UTC"), "UTC"
    try:
        return ZoneInfo(raw), raw
    except Exception:
        pass
    aliases = {
        "santiago, chili": "America/Santiago",
        "chili": "America/Santiago",
        "chile": "America/Santiago",
    }
    key = raw.lower()
    if key in aliases:
        try:
            return ZoneInfo(aliases[key]), aliases[key]
        except Exception:
            pass
    return ZoneInfo("UTC"), f"UTC (fuseau inconnu : {raw!r})"


def _sharpcap_event_params(row, event_label: str) -> dict:
    """Paramètres communs aux exports SharpCap (.scs) depuis une ligne de prédiction SORA."""
    ep = row["Epoch"]
    if hasattr(ep, "isscalar") and not ep.isscalar:
        ep = ep[0]
    dt_utc = ep.utc.to_datetime(timezone=timezone.utc)
    tz, tz_note = _observatory_zoneinfo()
    dt_loc = dt_utc.astimezone(tz)
    dt_loc_start = dt_loc - timedelta(minutes=2)

    tbl = getattr(row, "table", None)
    meta = getattr(tbl, "meta", None) or {}
    object_name = str(meta.get("name", "") or "").strip()
    if not object_name and " \u2014 " in event_label:
        object_name = event_label.split(" \u2014 ", 1)[0].strip()
    if not object_name:
        object_name = "Occultation"

    duration_s = 300
    try:
        vel_km_s = abs(float(row["Vel"]))
        rad_km = float(meta.get("radius", 500))
        if vel_km_s > 1e-6:
            maxdur_s = 2.0 * rad_km / vel_km_s
            duration_s = int(max(90, min(3600, maxdur_s * 2.0 + 120.0)))
    except (TypeError, ValueError, KeyError):
        duration_s = 300

    return {
        "object_name": object_name,
        "duration_s": duration_s,
        "iso_utc": dt_utc.strftime("%Y-%m-%d %H:%M:%S UTC"),
        "iso_loc": dt_loc.strftime("%Y-%m-%d %H:%M:%S"),
        "start_local_hms": dt_loc_start.strftime("%H:%M:%S"),
        "tz_note": tz_note,
    }


def _ephemeris_error_mas_for_map(row) -> Optional[float]:
    """
    Incertitude éphéméride (mas) pour ``plot_occ_map(..., error=...)`` : SORA trace des lignes
    pointillées à rayon ± erreur projetée. Métadonnées ``error_ra`` / ``error_dec`` (mas).
    """
    tbl = getattr(row, "table", None)
    if tbl is None:
        return None
    meta = getattr(tbl, "meta", None) or {}
    ra = meta.get("error_ra")
    dec = meta.get("error_dec")
    if ra is None and dec is None:
        return None
    try:
        r = float(ra) if ra is not None else 0.0
        d = float(dec) if dec is not None else 0.0
    except (TypeError, ValueError):
        return None
    m = max(abs(r), abs(d))
    return m if m > 0.0 else None


class OccultationsTab(ttk.Frame):
    """Occultations : éphémérides SORA, PyMovie, analyse SORA."""

    _MAX_PREDICTION_BODIES = 20
    _PRED_LEFT_COL_WEIGHT = 1
    _PRED_RIGHT_COL_WEIGHT = 2
    # Taille minimale d’affichage si Pillow est indisponible (PhotoImage.subsample)
    _MAP_VIEW_FALLBACK_W = 960
    _MAP_VIEW_FALLBACK_H = 720

    def __init__(self, parent, night_observation_tab=None):
        super().__init__(parent, padding=8)
        from gui.manual_help import add_manual_help_header

        add_manual_help_header(self, "13-occultations-sora-pymovie-analyse")
        self._project_root = Path(__file__).resolve().parent.parent
        self._npoap_dir = Path.home() / ".npoap"
        try:
            self._npoap_dir.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass
        self._pymovie_settings_path = self._npoap_dir / "occultations_pymovie.json"
        self._pymovie_conda_root_var = tk.StringVar(master=self, value=_default_pymovie_conda_root())
        self._load_pymovie_settings()
        self._pymovie_proc: Optional[subprocess.Popen] = None
        self._pyote_proc: Optional[subprocess.Popen] = None
        self._night_observation_tab = night_observation_tab
        self._pred_occ_labels = []
        self._pred_occ_rows = []
        self._pred_map_tk_photo = None
        self._pred_map_prefs_path = self._npoap_dir / "sora_prediction_map_prefs.json"
        self._pred_map_resize_job: Optional[Any] = None
        self._pred_last_map_png_path: Optional[Path] = None

        nb = ttk.Notebook(self)
        nb.pack(fill=tk.BOTH, expand=True)

        f_sora_pred = ttk.Frame(nb, padding=10)
        f_pym = ttk.Frame(nb, padding=10)
        f_sora = ttk.Frame(nb, padding=10)

        nb.add(f_sora_pred, text="1. SORA \u2014 Prédictions")
        nb.add(f_pym, text="2. Photométrie PyMovie")
        nb.add(f_sora, text="3. Analyse SORA")

        self._build_sora_prediction_panel(f_sora_pred)
        self._build_pymovie_panel(f_pym)
        self._build_sora_analysis_panel(f_sora)

    def _link_row(self, parent, label: str, url: str):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=2)
        ttk.Label(row, text=label, width=28, anchor="w").pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(row, text="Ouvrir", width=10, command=lambda u=url: _open_url(u)).pack(
            side=tk.LEFT
        )

    def _readonly_text(self, parent, text: str, height: int = 10):
        wrap = ttk.Frame(parent)
        wrap.pack(fill=tk.BOTH, expand=True, pady=(0, 4))
        t = tk.Text(
            wrap,
            height=height,
            wrap=tk.WORD,
            font=("Consolas", 9) if tk.TkVersion >= 8.6 else ("Courier", 9),
            relief=tk.FLAT,
            padx=6,
            pady=6,
        )
        t.insert("1.0", text)
        t.configure(state=tk.DISABLED, background=self._text_bg())
        sb = ttk.Scrollbar(wrap, command=t.yview)
        t.configure(yscrollcommand=sb.set)
        t.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

    def _text_bg(self):
        try:
            bg = ttk.Style().lookup("TFrame", "background")
            return bg if bg else "#f5f5f5"
        except tk.TclError:
            return "#f5f5f5"

    @staticmethod
    def _observatory_snapshot_from_config():
        """
        Lit ``config.OBSERVATORY`` (onglet Accueil / sauvegarde : latitude, longitude, elevation ;
        repli sur lat, lon, elev du fichier ``config.py``).
        """
        o = getattr(config, "OBSERVATORY", None) or {}
        name = str(o.get("name", "") or "").strip() or "Site"
        lat = o.get("latitude", o.get("lat", 0.0))
        lon = o.get("longitude", o.get("lon", 0.0))
        elev = o.get("elevation", o.get("elev", 0.0))
        try:
            lat_f = float(lat)
        except (TypeError, ValueError):
            lat_f = 0.0
        try:
            lon_f = float(lon)
        except (TypeError, ValueError):
            lon_f = 0.0
        try:
            elev_f = float(elev)
        except (TypeError, ValueError):
            elev_f = 0.0
        return name, lon_f, lat_f, elev_f

    def _refresh_pred_observatory_labels(self):
        if not hasattr(self, "_pred_site_name_lbl"):
            return
        name, lon, lat, h = self._observatory_snapshot_from_config()
        self._pred_site_name_lbl.config(text=name)
        self._pred_site_lon_lbl.config(text=f"{lon:.6f}")
        self._pred_site_lat_lbl.config(text=f"{lat:.6f}")
        self._pred_site_h_lbl.config(text=f"{h:.1f} m")

    def _build_sora_prediction_panel(self, parent):
        ttk.Label(
            parent,
            text="SORA \u2014 prédictions / éphémérides",
            font=("Helvetica", 12, "bold"),
        ).pack(anchor="w", pady=(0, 6))

        body = ttk.Frame(parent)
        body.pack(fill=tk.BOTH, expand=True)
        body.columnconfigure(
            0,
            weight=self._PRED_LEFT_COL_WEIGHT,
            uniform="sora_pred",
            minsize=280,
        )
        body.columnconfigure(
            1,
            weight=self._PRED_RIGHT_COL_WEIGHT,
            uniform="sora_pred",
            minsize=400,
        )
        body.rowconfigure(0, weight=1)

        left = ttk.Frame(body, padding=(0, 0, 8, 0))
        right = ttk.Frame(body, padding=(4, 0, 0, 0))
        left.grid(row=0, column=0, sticky="nsew")
        right.grid(row=0, column=1, sticky="nsew")

        if self._night_observation_tab is not None:
            import_row = ttk.Frame(left)
            import_row.pack(fill=tk.X, pady=(0, 6))
            ttk.Button(
                import_row,
                text="Importer depuis Observation de la Nuit",
                command=self._pred_import_from_night,
            ).pack(side=tk.LEFT)

        ff = ttk.LabelFrame(left, text="Paramètres de prédiction", padding=8)
        ff.pack(fill=tk.X, pady=(0, 0))

        self._pred_body_var = tk.StringVar(value="136199")
        self._pred_t0_var = tk.StringVar(value="2026-04-01T00:00:00")
        self._pred_t1_var = tk.StringVar(value="2026-04-02T00:00:00")
        self._pred_mag_var = tk.StringVar(value="")
        self._pred_step_var = tk.StringVar(value="60")
        self._pred_divs_var = tk.StringVar(value="1")
        self._pred_cat_var = tk.StringVar(value="gaiadr3")

        r = 0
        ttk.Label(ff, text="Corps unique (si liste vide)").grid(row=r, column=0, sticky="nw", padx=(0, 8), pady=2)
        ttk.Entry(ff, textvariable=self._pred_body_var, width=40).grid(row=r, column=1, sticky="ew", pady=2)
        r += 1
        ttk.Label(ff, text="Liste de corps").grid(row=r, column=0, sticky="nw", padx=(0, 8), pady=2)
        list_wrap = ttk.Frame(ff)
        list_wrap.grid(row=r, column=1, sticky="ew", pady=2)
        self._pred_body_list = tk.Text(list_wrap, height=5, width=40, font=("Consolas", 9), wrap=tk.NONE)
        lsb = ttk.Scrollbar(list_wrap, command=self._pred_body_list.yview)
        self._pred_body_list.configure(yscrollcommand=lsb.set)
        self._pred_body_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lsb.pack(side=tk.RIGHT, fill=tk.Y)
        self._pred_body_list.insert(
            "1.0",
            "136199\nEuropa\nPluto\n",
        )
        ttk.Label(
            ff,
            text="(une ligne = un objet ; si au moins une ligne non vide, ce champ remplace \u00ab corps unique \u00bb ; max 20)",
            font=("TkDefaultFont", 8),
        ).grid(row=r + 1, column=1, sticky="w", pady=(0, 4))
        r += 2
        ttk.Label(ff, text="Début (UTC, ISO)").grid(row=r, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(ff, textvariable=self._pred_t0_var, width=40).grid(row=r, column=1, sticky="ew", pady=2)
        r += 1
        ttk.Label(ff, text="Fin (UTC, ISO)").grid(row=r, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(ff, textvariable=self._pred_t1_var, width=40).grid(row=r, column=1, sticky="ew", pady=2)
        r += 1
        ttk.Label(ff, text="Mag. limite (optionnel)").grid(row=r, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(ff, textvariable=self._pred_mag_var, width=40).grid(row=r, column=1, sticky="ew", pady=2)
        r += 1
        ttk.Label(ff, text="Pas (s)").grid(row=r, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(ff, textvariable=self._pred_step_var, width=12).grid(row=r, column=1, sticky="w", pady=2)
        r += 1
        ttk.Label(ff, text="Découpages éphéméride (divs)").grid(row=r, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(ff, textvariable=self._pred_divs_var, width=12).grid(row=r, column=1, sticky="w", pady=2)
        r += 1
        ttk.Label(ff, text="Catalogue").grid(row=r, column=0, sticky="w", padx=(0, 8), pady=2)
        cat_cb = ttk.Combobox(
            ff,
            textvariable=self._pred_cat_var,
            values=("gaiadr3", "gaiaedr3", "gaiadr2", "ucac4"),
            width=18,
            state="readonly",
        )
        cat_cb.grid(row=r, column=1, sticky="w", pady=2)
        r += 1
        ttk.Label(
            ff,
            text=(
                "UCAC4 : VizieR I/322A (type OccultWatcher) ; requêtes plus lourdes. "
                "Mag. limite numérique \u2192 bande V."
            ),
            font=("TkDefaultFont", 8),
            wraplength=360,
        ).grid(row=r, column=0, columnspan=2, sticky="w", padx=(0, 0), pady=(0, 4))
        r += 1

        self._pred_obs_frame = ttk.LabelFrame(
            ff,
            text="Site d\u2019observation (config \u2014 onglet Accueil)",
            padding=4,
        )
        self._pred_obs_frame.grid(row=r, column=0, columnspan=2, sticky="ew", pady=(6, 4))
        r += 1
        ff.columnconfigure(1, weight=1)

        ofg = self._pred_obs_frame
        ttk.Label(ofg, text="Nom").grid(row=0, column=0, sticky="w", padx=4)
        self._pred_site_name_lbl = ttk.Label(ofg, text="\u2014", anchor="w")
        self._pred_site_name_lbl.grid(row=0, column=1, sticky="ew", padx=4)
        ttk.Label(ofg, text="Longitude (\u00b0)").grid(row=1, column=0, sticky="w", padx=4)
        self._pred_site_lon_lbl = ttk.Label(ofg, text="\u2014", anchor="w", font=("Consolas", 9))
        self._pred_site_lon_lbl.grid(row=1, column=1, sticky="w", padx=4)
        ttk.Label(ofg, text="Latitude (\u00b0)").grid(row=2, column=0, sticky="w", padx=4)
        self._pred_site_lat_lbl = ttk.Label(ofg, text="\u2014", anchor="w", font=("Consolas", 9))
        self._pred_site_lat_lbl.grid(row=2, column=1, sticky="w", padx=4)
        ttk.Label(ofg, text="Altitude").grid(row=3, column=0, sticky="w", padx=4)
        self._pred_site_h_lbl = ttk.Label(ofg, text="\u2014", anchor="w", font=("Consolas", 9))
        self._pred_site_h_lbl.grid(row=3, column=1, sticky="w", padx=4)
        ttk.Label(
            ofg,
            text="Modifiable dans l\u2019onglet Accueil ; rechargé à chaque lancement de prédiction.",
            font=("TkDefaultFont", 8),
        ).grid(row=4, column=0, columnspan=2, sticky="w", padx=4, pady=(2, 0))
        ofg.columnconfigure(1, weight=1)

        btn_row = ttk.Frame(left)
        btn_row.pack(fill=tk.X, pady=(8, 0))
        ttk.Button(btn_row, text="Lancer la prédiction", command=self._pred_run).pack(side=tk.LEFT, padx=(0, 8))
        self._pred_map_btn = ttk.Button(btn_row, text="Voir la carte", command=self._pred_show_occ_map)
        self._pred_map_btn.pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(btn_row, text="Exporter la table (TXT)", command=self._pred_export_table_txt).pack(
            side=tk.LEFT, padx=(0, 8)
        )
        ttk.Button(btn_row, text="Exporter SharpCap (.scs)", command=self._pred_export_sharpcap_scs).pack(
            side=tk.LEFT
        )

        self._sora_pred_status = ttk.Label(left, text="", wraplength=360)
        self._sora_pred_status.pack(anchor="w", pady=(6, 0))

        rf = ttk.LabelFrame(left, text="Résultat (table de prédiction)", padding=4)
        rf.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        rw = ttk.Frame(rf)
        rw.pack(fill=tk.BOTH, expand=True)
        self._pred_result = tk.Text(
            rw,
            height=16,
            wrap=tk.NONE,
            font=("Consolas", 8),
            relief=tk.FLAT,
            padx=4,
            pady=4,
        )
        rsb = ttk.Scrollbar(rw, command=self._pred_result.yview)
        self._pred_result.configure(yscrollcommand=rsb.set)
        self._pred_result.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        rsb.pack(side=tk.RIGHT, fill=tk.Y)
        self._pred_result.insert("1.0", "Résultat affiché ici après \u00ab Lancer la prédiction \u00bb.")

        map_lf = ttk.LabelFrame(right, text="Carte d\u2019occultation", padding=6)
        map_lf.pack(fill=tk.BOTH, expand=True)
        ev_row = ttk.Frame(map_lf)
        ev_row.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(ev_row, text="Occultation :").pack(side=tk.LEFT, padx=(0, 6))
        self._pred_occ_combo = ttk.Combobox(ev_row, width=32, state="readonly")
        self._pred_occ_combo.pack(side=tk.LEFT, fill=tk.X, expand=True)

        self._pred_map_res_var = tk.IntVar(master=self, value=2)
        self._pred_map_zoom_var = tk.StringVar(master=self, value="4")
        self._pred_map_states_var = tk.BooleanVar(master=self, value=True)
        self._pred_map_labels_var = tk.BooleanVar(master=self, value=True)
        self._pred_map_meridian_var = tk.StringVar(master=self, value="30")
        self._pred_map_parallels_var = tk.StringVar(master=self, value="30")
        self._pred_map_style_var = tk.IntVar(master=self, value=1)
        self._pred_map_dpi_var = tk.StringVar(master=self, value="120")
        self._pred_map_fmt_var = tk.StringVar(master=self, value="png")
        self._pred_map_size_w_var = tk.StringVar(master=self, value="42")
        self._pred_map_size_h_var = tk.StringVar(master=self, value="44")
        self._pred_map_alpha_var = tk.StringVar(master=self, value="0.2")
        self._pred_map_pref_vars = {
            "resolution": self._pred_map_res_var,
            "zoom": self._pred_map_zoom_var,
            "states": self._pred_map_states_var,
            "labels": self._pred_map_labels_var,
            "meridian": self._pred_map_meridian_var,
            "parallels": self._pred_map_parallels_var,
            "mapstyle": self._pred_map_style_var,
            "dpi": self._pred_map_dpi_var,
            "fmt": self._pred_map_fmt_var,
            "mapsize_w": self._pred_map_size_w_var,
            "mapsize_h": self._pred_map_size_h_var,
            "alpha": self._pred_map_alpha_var,
        }
        self._pred_load_map_prefs()

        opt_lf = ttk.LabelFrame(map_lf, text="Options carte (SORA, plot_occ_map)", padding=6)
        opt_lf.pack(fill=tk.X, pady=(0, 6))
        og = opt_lf
        r0 = 0
        ttk.Label(og, text="Résolution Cartopy (1=10 m, 2=50 m, 3=100 m)").grid(
            row=r0, column=0, sticky="w", padx=(0, 6), pady=1
        )
        tk.Spinbox(og, from_=1, to=3, width=4, textvariable=self._pred_map_res_var).grid(
            row=r0, column=1, sticky="w", pady=1
        )
        ttk.Label(og, text="Zoom").grid(row=r0, column=2, sticky="w", padx=(12, 6), pady=1)
        ttk.Entry(og, textvariable=self._pred_map_zoom_var, width=8).grid(
            row=r0, column=3, sticky="w", pady=1
        )
        r0 += 1
        ttk.Checkbutton(
            og,
            text="Frontières (states)",
            variable=self._pred_map_states_var,
            command=self._pred_save_map_prefs,
        ).grid(row=r0, column=0, columnspan=2, sticky="w", pady=1)
        ttk.Checkbutton(
            og,
            text="Légendes (labels)",
            variable=self._pred_map_labels_var,
            command=self._pred_save_map_prefs,
        ).grid(row=r0, column=2, columnspan=2, sticky="w", pady=1)
        r0 += 1
        ttk.Label(og, text="Méridien (°)").grid(row=r0, column=0, sticky="w", padx=(0, 6), pady=1)
        ttk.Entry(og, textvariable=self._pred_map_meridian_var, width=6).grid(
            row=r0, column=1, sticky="w", pady=1
        )
        ttk.Label(og, text="Parallèles (°)").grid(row=r0, column=2, sticky="w", padx=(12, 6), pady=1)
        ttk.Entry(og, textvariable=self._pred_map_parallels_var, width=6).grid(
            row=r0, column=3, sticky="w", pady=1
        )
        r0 += 1
        ttk.Label(og, text="Style (1=N&B, 2=couleur)").grid(
            row=r0, column=0, sticky="w", padx=(0, 6), pady=1
        )
        tk.Spinbox(og, from_=1, to=2, width=4, textvariable=self._pred_map_style_var).grid(
            row=r0, column=1, sticky="w", pady=1
        )
        ttk.Label(og, text="DPI").grid(row=r0, column=2, sticky="w", padx=(12, 6), pady=1)
        ttk.Entry(og, textvariable=self._pred_map_dpi_var, width=6).grid(
            row=r0, column=3, sticky="w", pady=1
        )
        r0 += 1
        ttk.Label(og, text="Format fichier").grid(row=r0, column=0, sticky="w", padx=(0, 6), pady=1)
        fmt_cb = ttk.Combobox(
            og,
            textvariable=self._pred_map_fmt_var,
            values=("png", "pdf"),
            width=8,
            state="readonly",
        )
        fmt_cb.grid(row=r0, column=1, sticky="w", pady=1)
        fmt_cb.bind("<<ComboboxSelected>>", lambda _e: self._pred_save_map_prefs())
        ttk.Label(
            og,
            text="(aperçu intégré : PNG uniquement)",
            font=("TkDefaultFont", 8),
        ).grid(row=r0, column=2, columnspan=2, sticky="w", padx=(8, 0), pady=1)
        r0 += 1
        ttk.Label(og, text="Taille figure (cm) — l × h").grid(
            row=r0, column=0, sticky="w", padx=(0, 6), pady=1
        )
        szr = ttk.Frame(og)
        szr.grid(row=r0, column=1, columnspan=3, sticky="w", pady=1)
        ttk.Entry(szr, textvariable=self._pred_map_size_w_var, width=6).pack(side=tk.LEFT)
        ttk.Label(szr, text="\u00d7").pack(side=tk.LEFT, padx=4)
        ttk.Entry(szr, textvariable=self._pred_map_size_h_var, width=6).pack(side=tk.LEFT)
        r0 += 1
        ttk.Label(og, text="Alpha ombre nuit (0–1)").grid(
            row=r0, column=0, sticky="w", padx=(0, 6), pady=1
        )
        ttk.Entry(og, textvariable=self._pred_map_alpha_var, width=8).grid(
            row=r0, column=1, sticky="w", pady=1
        )
        ttk.Label(
            og,
            text="Préréglages enregistrés dans .npoap/sora_prediction_map_prefs.json",
            font=("TkDefaultFont", 8),
        ).grid(row=r0 + 1, column=0, columnspan=4, sticky="w", pady=(4, 0))

        def _bind_pref_save(widget):
            if isinstance(widget, (ttk.Entry, ttk.Combobox)):
                widget.bind("<FocusOut>", lambda _e: self._pred_save_map_prefs())
            elif isinstance(widget, tk.Spinbox):
                widget.bind("<FocusOut>", lambda _e: self._pred_save_map_prefs())
                widget.bind("<ButtonRelease-1>", lambda _e: self._pred_save_map_prefs())

        for w in og.winfo_children():
            _bind_pref_save(w)
            if isinstance(w, ttk.Frame):
                for c in w.winfo_children():
                    _bind_pref_save(c)

        map_canvas_wrap = tk.Frame(map_lf, bg=self._text_bg())
        map_canvas_wrap.pack(fill=tk.BOTH, expand=True, pady=(4, 0))
        self._pred_map_canvas_wrap = map_canvas_wrap
        map_canvas_wrap.bind("<Configure>", self._pred_on_map_canvas_configure)
        self._pred_map_placeholder = ttk.Label(
            map_canvas_wrap,
            text=(
                "Lancez une prédiction, sélectionnez une occultation ci-dessus, "
                "puis \u00ab Voir la carte \u00bb. La figure est générée par SORA (matplotlib / cartopy)."
            ),
            wraplength=400,
            justify=tk.LEFT,
        )
        self._pred_map_placeholder.pack(anchor="nw", pady=(0, 6), padx=4)
        self._pred_map_image = tk.Label(map_canvas_wrap, text="", bg=self._text_bg())
        self._pred_map_image.pack(expand=True, fill=tk.BOTH)

        self._refresh_pred_observatory_labels()

    @staticmethod
    def _pred_row_epoch_label(tab, j):
        try:
            ep = tab["Epoch"][j]
            if hasattr(ep, "iso"):
                s = ep.iso
                if hasattr(s, "value"):
                    s = s.value
                return str(s).replace("T", " ")
            return str(ep)
        except Exception:
            pass
        try:
            row = tab[j]
            occs = getattr(row, "occs", None)
            if occs is not None:
                d = getattr(occs, "datas", None)
                if d is not None and hasattr(d, "isot"):
                    return str(d.isot)
        except Exception:
            pass
        return f"événement {j + 1}"

    def _clear_pred_map_display(self):
        self._pred_map_tk_photo = None
        self._pred_last_map_png_path = None
        if not hasattr(self, "_pred_map_image"):
            return
        self._pred_map_image.configure(image="")
        if self._pred_map_placeholder.winfo_manager() == "":
            self._pred_map_placeholder.pack(anchor="nw", pady=(0, 6), before=self._pred_map_image)

    def _pred_load_map_prefs(self) -> None:
        try:
            if not self._pred_map_prefs_path.is_file():
                return
            with self._pred_map_prefs_path.open("r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            return
        if not isinstance(data, dict):
            return
        for key, var in getattr(self, "_pred_map_pref_vars", {}).items():
            if key not in data:
                continue
            val = data[key]
            try:
                if isinstance(var, tk.BooleanVar):
                    var.set(bool(val))
                elif isinstance(var, tk.IntVar):
                    var.set(int(round(float(val))))
                elif isinstance(var, tk.StringVar):
                    var.set("" if val is None else str(val))
                elif isinstance(var, tk.DoubleVar):
                    var.set(float(val))
            except (tk.TclError, ValueError, TypeError):
                pass

    def _pred_save_map_prefs(self) -> None:
        d = {}
        for key, var in getattr(self, "_pred_map_pref_vars", {}).items():
            try:
                if isinstance(var, tk.BooleanVar):
                    d[key] = bool(var.get())
                elif isinstance(var, (tk.IntVar, tk.StringVar, tk.DoubleVar)):
                    d[key] = var.get()
            except tk.TclError:
                continue
        try:
            with self._pred_map_prefs_path.open("w", encoding="utf-8") as f:
                json.dump(d, f, indent=2, ensure_ascii=False)
        except Exception:
            pass

    def _pred_parse_float(self, raw: str, default: float) -> float:
        try:
            return float((raw or "").strip().replace(",", "."))
        except (TypeError, ValueError):
            return default

    def _pred_parse_int(self, raw: str, default: int) -> int:
        try:
            return int(float((raw or "").strip().replace(",", ".")))
        except (TypeError, ValueError):
            return default

    def _pred_collect_plot_occ_map_kw(self) -> dict:
        """Paramètres SORA plot_occ_map (doc : cartes d’occultation)."""
        zoom = self._pred_parse_float(self._pred_map_zoom_var.get(), 4.0)
        if zoom <= 0:
            zoom = 4.0
        res = self._pred_map_res_var.get()
        if res not in (1, 2, 3):
            res = 2
        mer = max(1, self._pred_parse_int(self._pred_map_meridian_var.get(), 30))
        par = max(1, self._pred_parse_int(self._pred_map_parallels_var.get(), 30))
        dpi = max(30, min(600, self._pred_parse_int(self._pred_map_dpi_var.get(), 120)))
        ms_w = max(8.0, min(120.0, self._pred_parse_float(self._pred_map_size_w_var.get(), 42.0)))
        ms_h = max(8.0, min(120.0, self._pred_parse_float(self._pred_map_size_h_var.get(), 44.0)))
        fmt = (self._pred_map_fmt_var.get() or "png").strip().lower()
        if fmt not in ("png", "pdf", "svg"):
            fmt = "png"
        mstyle = self._pred_map_style_var.get()
        if mstyle not in (1, 2):
            mstyle = 1
        alpha = self._pred_parse_float(self._pred_map_alpha_var.get(), 0.2)
        alpha = max(0.0, min(1.0, alpha))
        return dict(
            fmt=fmt,
            dpi=int(dpi),
            mapsize=[float(ms_w), float(ms_h)],
            resolution=int(res),
            zoom=float(zoom),
            states=bool(self._pred_map_states_var.get()),
            labels=bool(self._pred_map_labels_var.get()),
            meridian=int(mer),
            parallels=int(par),
            mapstyle=int(mstyle),
            alpha=float(alpha),
        )

    def _pred_on_map_canvas_configure(self, _event=None) -> None:
        if self._pred_map_resize_job is not None:
            try:
                self.after_cancel(self._pred_map_resize_job)
            except Exception:
                pass
        self._pred_map_resize_job = self.after(150, self._pred_redraw_map_to_canvas)

    def _pred_redraw_map_to_canvas(self) -> None:
        self._pred_map_resize_job = None
        png_path = self._pred_last_map_png_path
        if png_path is None or not png_path.is_file():
            return
        if not hasattr(self, "_pred_map_canvas_wrap"):
            return
        try:
            cw = int(self._pred_map_canvas_wrap.winfo_width())
            ch = int(self._pred_map_canvas_wrap.winfo_height())
        except tk.TclError:
            return
        if cw < 32 or ch < 32:
            return
        bg = self._text_bg()
        try:
            from PIL import Image, ImageChops, ImageOps, ImageTk

            im = Image.open(str(png_path)).convert("RGBA")
            # Retire les bordures uniformes de la figure SORA pour mieux remplir le conteneur.
            try:
                bg_ref = Image.new("RGBA", im.size, im.getpixel((0, 0)))
                diff = ImageChops.difference(im, bg_ref)
                bbox = diff.getbbox()
                if bbox:
                    im = im.crop(bbox)
            except Exception:
                pass
            try:
                _resample = Image.Resampling.LANCZOS
            except AttributeError:
                _resample = Image.LANCZOS  # type: ignore[attr-defined]
            # "cover" : remplit tout l'espace disponible (peut rogner légèrement les bords).
            im = ImageOps.fit(im, (cw, ch), method=_resample, centering=(0.5, 0.5))
            photo = ImageTk.PhotoImage(im)
        except Exception:
            try:
                photo = tk.PhotoImage(file=str(png_path))
                mw, mh = max(cw, self._MAP_VIEW_FALLBACK_W), max(ch, self._MAP_VIEW_FALLBACK_H)
                while photo.width() > mw or photo.height() > mh:
                    photo = photo.subsample(2, 2)
            except tk.TclError as ex2:
                logger.error("Affichage carte PNG : %s", ex2)
                return
        self._pred_map_tk_photo = photo
        self._pred_map_image.configure(image=photo)

    def _pred_export_table_txt(self):
        """Exporte le contenu de la zone résultat (table formatée) vers un fichier texte UTF-8."""
        txt = self._pred_result.get("1.0", tk.END).strip()
        placeholder = "Résultat affiché ici après"
        if not txt or txt.startswith(placeholder) or txt == "Calcul\u2026":
            logger.warning("Export TXT : aucune table de prédiction à exporter.")
            self._sora_pred_status.configure(
                text="Export : lancez d'abord une prédiction, puis exportez."
            )
            return
        out_base = self._project_root / "output"
        out_base.mkdir(parents=True, exist_ok=True)
        path = filedialog.asksaveasfilename(
            title="Exporter la table de prédiction",
            initialdir=str(out_base),
            initialfile="prediction_sora.txt",
            defaultextension=".txt",
            filetypes=[("Texte", "*.txt"), ("Tous les fichiers", "*.*")],
        )
        if not path:
            return
        try:
            Path(path).write_text(txt + "\n", encoding="utf-8")
        except OSError as ex:
            logger.error("Export TXT impossible : %s", ex, exc_info=True)
            self._sora_pred_status.configure(text=f"Export : échec écriture — {ex}")
            return
        logger.info("Table de prédiction exportée : %s", path)
        self._sora_pred_status.configure(text=f"Table exportée : {path}")

    def _pred_get_export_prediction_row(self):
        """Occultation sélectionnée pour export ; met à jour le statut si indisponible."""
        rows = getattr(self, "_pred_occ_rows", None) or []
        if not rows:
            logger.warning("Export SCS : aucune prédiction en mémoire.")
            self._sora_pred_status.configure(
                text="Export SCS : lancez une prédiction et sélectionnez une occultation."
            )
            return None
        idx = self._pred_occ_combo.current()
        if idx < 0 or idx >= len(rows):
            logger.warning("Export SCS : aucune occultation sélectionnée.")
            self._sora_pred_status.configure(text="Export SCS : choisissez une occultation dans la liste.")
            return None
        row = rows[idx]
        labels = getattr(self, "_pred_occ_labels", None) or []
        label = labels[idx] if idx < len(labels) else f"occultation #{idx + 1}"
        return row, label

    def _pred_build_sharpcap_scs(self, row, event_label: str) -> str:
        """
        Séquence minimale type « SharpCap Just Record » (Occultation Manager / labstercam).
        Pour l'automatisation complète UTC (recommandée), copiez les modèles du dépôt officiel.
        """
        p = _sharpcap_event_params(row, event_label)
        iso_utc = p["iso_utc"]
        iso_loc = p["iso_loc"]
        start_local_hms = p["start_local_hms"]
        tz_note = p["tz_note"]
        object_name = p["object_name"]
        duration_s = p["duration_s"]

        doc_repo = "https://github.com/labstercam/occultation-tools/tree/main"
        doc_readme = (
            "https://github.com/labstercam/occultation-tools/blob/main/occultation-manager/ReadMe.md"
        )
        tpl_utc = (
            "https://github.com/labstercam/occultation-tools/blob/main/"
            "occultation-manager/python/SharpCap%20Sequence%20UTC%20template.txt"
        )
        tpl_just = (
            "https://github.com/labstercam/occultation-tools/blob/main/"
            "occultation-manager/python/SharpCap%20Just%20Record%20template.txt"
        )
        countdown_scs = (
            "https://github.com/labstercam/occultation-tools/blob/main/"
            "occultation-manager/python/countdown%20python%20for%20sequencer.scs"
        )

        lines = [
            "# NPOAP — export séquence SharpCap (.scs)",
            "#",
            "# Référence : Occultation Manager (Michael Camilleri) — modèle proche de",
            "# « SharpCap Just Record template » (WAIT UNTIL LOCALTIME + CAPTURE … SECONDS).",
            f"# Documentation : {doc_readme}",
            f"# Dépôt : {doc_repo}",
            f"# Modèle UTC recommandé (compte à rebours sûr, minuit) : {tpl_utc}",
            f"# Modèle « Just Record » source : {tpl_just}",
            f"# Snippets countdown UTC (à coller pour séquences avancées) : {countdown_scs}",
            "#",
            f"# Événement : {event_label}",
            f"# C/A SORA (UTC) : {iso_utc}",
            f"# Début enregistrement proposé (local, config Accueil) : {iso_loc} (−2 min avant C/A)",
            f"# Fuseau pour conversion : {tz_note}",
            "#",
            "# ATTENTION (doc Occultation Manager) : WAIT UNTIL LOCALTIME ne connaît pas la date ;",
            "# l'heure est celle de l'horloge Windows sous SharpCap. Pour nuit / minuit, préférez",
            "# le modèle UTC + countdown du dépôt ci-dessus.",
            "#",
            "# Doc éditeur de séquences SharpCap :",
            "# https://docs.sharpcap.co.uk/4.0/27_TheSharpCapSequenceEditor.htm",
            "#",
            "SEQUENCE",
            f' TARGETNAME "{object_name}"',
            (
                f' SHOW NOTIFICATION "NPOAP : attente début enregistrement ({start_local_hms} local)" '
                f"COLOUR Green DURATION 58"
            ),
            f' WAIT UNTIL LOCALTIME "{start_local_hms}"',
            f" CAPTURE {duration_s} SECONDS LIVE FRAMES",
            ' SHOW NOTIFICATION "NPOAP : enregistrement terminé." COLOUR GREEN DURATION 30',
            "END SEQUENCE",
            "",
        ]
        return "\n".join(lines)

    def _pred_export_sharpcap_scs(self):
        """Exporte une séquence .scs pour SharpCap à partir de l'occultation sélectionnée."""
        sel = self._pred_get_export_prediction_row()
        if sel is None:
            return
        row, label = sel
        body = self._pred_build_sharpcap_scs(row, label)
        out_base = self._project_root / "output"
        out_base.mkdir(parents=True, exist_ok=True)
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in label[:40]).strip("_") or "occultation"
        path = filedialog.asksaveasfilename(
            title="Exporter une séquence SharpCap (.scs)",
            initialdir=str(out_base),
            initialfile=f"npoap_sharpcap_{safe}.scs",
            defaultextension=".scs",
            filetypes=[("SharpCap sequence", "*.scs"), ("Tous les fichiers", "*.*")],
        )
        if not path:
            return
        try:
            Path(path).write_text(body, encoding="utf-8", newline="\r\n")
        except OSError as ex:
            logger.error("Export SCS impossible : %s", ex, exc_info=True)
            self._sora_pred_status.configure(text=f"Export SCS : échec — {ex}")
            return
        logger.info("Séquence SharpCap exportée : %s", path)
        self._sora_pred_status.configure(text=f"Séquence SharpCap exportée : {path}")

    def _pred_collect_bodies(self):
        """Une ligne du texte = un corps ; sinon le champ \u00ab corps unique \u00bb."""
        raw = self._pred_body_list.get("1.0", tk.END)
        lines = [
            ln.strip()
            for ln in raw.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        if lines:
            return lines[: self._MAX_PREDICTION_BODIES]
        one = self._pred_body_var.get().strip()
        return [one] if one else []

    def _pred_import_from_night(self):
        """Remplit liste / dates / site depuis l'onglet Observation de la Nuit (sélection)."""
        tab = self._night_observation_tab
        if tab is None:
            logger.warning("Import nuit impossible : onglet Observation de la Nuit non disponible.")
            self._sora_pred_status.configure(text="Import nuit : onglet non disponible.")
            return
        try:
            bodies, time_beg, time_end, _, _, _, _ = tab.get_export_bundle_for_sora()
        except Exception as ex:
            logger.error("Impossible de lire l'onglet nuit : %s", ex, exc_info=True)
            self._sora_pred_status.configure(text=f"Import nuit : erreur \u2014 {ex}")
            return
        if not bodies:
            logger.warning("Import nuit : aucun asteroide/comete coche dans Objets observables.")
            self._sora_pred_status.configure(
                text="Import nuit : aucun astéroïde ou comète coché dans Objets observables."
            )
            return
        if not time_beg or not time_end:
            logger.warning("Import nuit : date d'observation invalide ou vide.")
            self._sora_pred_status.configure(
                text="Import nuit : date d'observation invalide ou vide (attendu : AAAA-MM-JJ)."
            )
            return
        if len(bodies) > self._MAX_PREDICTION_BODIES:
            bodies = bodies[: self._MAX_PREDICTION_BODIES]
            logger.info("Import nuit : liste tronquee aux %d premiers corps.", self._MAX_PREDICTION_BODIES)
        self._pred_body_list.delete("1.0", tk.END)
        self._pred_body_list.insert("1.0", "\n".join(bodies) + "\n")
        self._pred_t0_var.set(time_beg)
        self._pred_t1_var.set(time_end)
        if not self._pred_mag_var.get().strip():
            self._pred_mag_var.set("16")
        self._refresh_pred_observatory_labels()
        self._sora_pred_status.configure(
            text=f"Import nuit : {len(bodies)} corps, {time_beg} \u2192 {time_end}. "
            "Site : config (onglet Accueil)."
        )

    def _pred_run(self):
        ok, msg = try_import_sora()
        if not ok:
            logger.error("SORA non disponible : %s", msg)
            self._sora_pred_status.configure(
                text=f"SORA non disponible : {msg}. "
                "Réinstallez via astroenv (pip install -r requirements.txt)."
            )
            return

        mag_raw = self._pred_mag_var.get().strip()
        mag_lim = None
        if mag_raw:
            try:
                mag_lim = float(mag_raw.replace(",", "."))
            except ValueError:
                logger.warning("Mag. limite invalide : %r", mag_raw)
                self._sora_pred_status.configure(
                    text="Mag. limite : nombre invalide (laisser vide pour aucune limite)."
                )
                return
        try:
            step = float(self._pred_step_var.get().replace(",", "."))
            divs = int(self._pred_divs_var.get().strip())
        except ValueError:
            logger.warning("Pas (s) ou divs invalide.")
            self._sora_pred_status.configure(text="Pas (s) et divs : valeurs numériques requises.")
            return

        self._refresh_pred_observatory_labels()
        obs_name, lon, lat, h_m = self._observatory_snapshot_from_config()

        if lon == 0.0 and lat == 0.0:
            logger.warning("Coordonnees observatoire a (0, 0) — verifiez config / onglet Accueil.")
            self._sora_pred_status.configure(
                text="Coordonnées observatoire (0, 0). Vérifiez l'onglet Accueil (latitude, longitude)."
            )
            return

        logger.info(
            "Prediction SORA : site=%s lon=%.6f lat=%.6f alt=%.1f m",
            obs_name, lon, lat, h_m,
        )
        self._sora_pred_status.configure(
            text=f"Prédiction en cours \u2014 site : {obs_name} (lon {lon:.4f}\u00b0, lat {lat:.4f}\u00b0, alt {h_m:.0f} m)\u2026"
        )
        self._pred_result.configure(state=tk.NORMAL)
        self._pred_result.delete("1.0", tk.END)
        self._pred_result.insert("1.0", "Calcul\u2026")
        self._pred_occ_labels = []
        self._pred_occ_rows = []
        if hasattr(self, "_pred_occ_combo"):
            self._pred_occ_combo.configure(values=[])
            self._pred_occ_combo.set("")
        self._clear_pred_map_display()
        self.update_idletasks()

        bodies = self._pred_collect_bodies()
        if not bodies:
            logger.warning("Aucun corps indique pour la prediction.")
            self._sora_pred_status.configure(text="Indiquez au moins un corps (liste ou champ unique).")
            self._pred_result.delete("1.0", tk.END)
            return
        raw_lines = [
            ln.strip()
            for ln in self._pred_body_list.get("1.0", tk.END).splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        if len(raw_lines) > self._MAX_PREDICTION_BODIES:
            logger.info("Liste tronquee aux %d premiers corps.", self._MAX_PREDICTION_BODIES)
            self._sora_pred_status.configure(
                text=f"Seuls les {self._MAX_PREDICTION_BODIES} premiers corps seront traités."
            )

        t0 = self._pred_t0_var.get().strip()
        t1 = self._pred_t1_var.get().strip()
        cat = self._pred_cat_var.get().strip() or "gaiadr3"

        nb = len(bodies)

        def work():
            parts = [
                f"Site : {obs_name} \u2014 lon {lon:.6f}\u00b0, lat {lat:.6f}\u00b0, alt {h_m:.1f} m\n"
                f"Intervalle : {t0}  \u2192  {t1}  |  Catalogue : {cat}\n"
            ]
            stored_rows = []
            try:
                for i, body in enumerate(bodies):

                    def upd(ii=i, b=body):
                        self._sora_pred_status.configure(
                            text=f"Prédiction {ii + 1}/{nb} : {b} (réseau + catalogue)\u2026"
                        )

                    self.after(0, upd)
                    try:
                        tab = run_occultation_prediction(
                            body,
                            t0,
                            t1,
                            mag_lim=mag_lim,
                            step=step,
                            divs=divs,
                            catalogue=cat,
                            use_geocenter=False,
                            observer_name=obs_name,
                            lon_deg=lon,
                            lat_deg=lat,
                            height_m=h_m,
                            verbose=False,
                        )
                        parts.append(f"\n{'=' * 72}\n  Corps : {body}\n{'=' * 72}\n")
                        parts.append(format_prediction_table(tab))
                        try:
                            nrows = len(tab)
                        except Exception:
                            nrows = 0
                        for j in range(nrows):
                            ep_str = self._pred_row_epoch_label(tab, j)
                            stored_rows.append((f"{body} \u2014 {ep_str}", tab[j]))
                    except Exception as ex:
                        logger.error("Prediction corps=%s : %s", body, ex, exc_info=True)
                        parts.append(
                            f"\n{'=' * 72}\n  Corps : {body} \u2014 ERREUR\n{'=' * 72}\n{ex}\n"
                        )

                txt = "".join(parts).lstrip("\n")

                def done_ok():
                    self._pred_occ_labels = [lbl for lbl, _ in stored_rows]
                    self._pred_occ_rows = [r for _, r in stored_rows]
                    if hasattr(self, "_pred_occ_combo"):
                        self._pred_occ_combo.configure(values=self._pred_occ_labels)
                        if self._pred_occ_labels:
                            self._pred_occ_combo.current(0)
                        else:
                            self._pred_occ_combo.set("")
                    self._sora_pred_status.configure(
                        text=(
                            f"Terminé : {nb} prédiction(s) \u2014 site : {obs_name} "
                            f"(lon {lon:.4f}\u00b0, lat {lat:.4f}\u00b0)."
                        )
                    )
                    self._pred_result.delete("1.0", tk.END)
                    self._pred_result.insert("1.0", txt or "(aucun résultat)")
                    self._pred_result.see("1.0")
                    logger.info("Prediction terminee : %d corps, %d occultations.", nb, len(stored_rows))

                self.after(0, done_ok)
            except Exception as e:
                err = str(e)
                logger.error("Prediction globale : %s", err, exc_info=True)

                def done_err():
                    self._sora_pred_status.configure(text=f"Erreur : {err[:200]}")
                    self._pred_result.delete("1.0", tk.END)
                    self._pred_result.insert("1.0", err)

                self.after(0, done_err)

        threading.Thread(target=work, daemon=True).start()

    def _pred_show_occ_map(self):
        ok, msg = try_import_sora()
        if not ok:
            logger.error("SORA non disponible pour la carte : %s", msg)
            self._sora_pred_status.configure(
                text=f"Carte : SORA non disponible ({msg}). Réinstallez via astroenv."
            )
            return
        rows = getattr(self, "_pred_occ_rows", None) or []
        if not rows:
            logger.warning("Carte : aucune ligne de prediction en memoire.")
            self._sora_pred_status.configure(
                text="Carte : aucune prédiction en mémoire. Lancez d'abord une prédiction."
            )
            return
        idx = self._pred_occ_combo.current()
        if idx < 0 or idx >= len(rows):
            logger.warning("Carte : aucune occultation selectionnee (idx=%d, rows=%d).", idx, len(rows))
            self._sora_pred_status.configure(
                text="Carte : sélectionnez une occultation dans la liste (droite)."
            )
            return
        row = rows[idx]
        self._pred_map_btn.configure(state=tk.DISABLED)
        self._sora_pred_status.configure(text="Génération de la carte d'occultation (SORA)\u2026")
        self.after(0, lambda r=row: self._pred_plot_occ_map_on_main_thread(r))

    def _pred_plot_occ_map_on_main_thread(self, row):
        """
        SORA plot_occ_map + matplotlib/cartopy doivent s'exécuter sur le thread Tk principal
        (évite RuntimeError sur plt.close et backends TkAgg). Avec zoom, SORA trace une flèche
        via plt.annotate avec des scalaires numpy 0-dim : matplotlib échoue ; ``arrow=False`` contourne ce bug.
        """
        err = None
        png_path = None
        try:
            import matplotlib.pyplot as plt

            out_dir = self._project_root / "output" / "sora_occ_maps"
            out_dir.mkdir(parents=True, exist_ok=True)
            obs_name, obs_lon, obs_lat, _ = self._observatory_snapshot_from_config()
            obs_lon = float(obs_lon)
            obs_lat = float(obs_lat)
            site_dict = {
                obs_name: [
                    obs_lon,
                    obs_lat,
                    6,
                    10,
                    "#c62828",
                    "^",
                ],
            }
            err_mas = _ephemeris_error_mas_for_map(row)
            if err_mas is not None:
                logger.info("plot_occ_map : bande d'incertitude éphéméride error=%.3f mas", err_mas)
            else:
                logger.info(
                    "plot_occ_map : pas d'error_ra/error_dec dans les métadonnées ; "
                    "pas de lignes d'incertitude (voir doc SORA, paramètre error=)."
                )
            opts = self._pred_collect_plot_occ_map_kw()
            map_kw = dict(
                path=str(out_dir),
                nameimg="npoap_occ_map",
                sites=site_dict,
                centerproj=[obs_lon, obs_lat],
                centermap_geo=[obs_lon, obs_lat],
                arrow=False,
                sscale=2.4,
                site_box_alpha=0.88,
                **opts,
            )
            if err_mas is not None:
                map_kw["error"] = err_mas
                map_kw["lncolor"] = "darkorange"
            ext = str(map_kw.get("fmt", "png") or "png").lower()
            if ext not in ("png", "pdf"):
                ext = "png"
                map_kw["fmt"] = "png"
            logger.info(
                "plot_occ_map : site=%s lon=%.6f lat=%.6f, opts=%s",
                obs_name, obs_lon, obs_lat, {k: map_kw[k] for k in sorted(map_kw) if k != "sites"},
            )

            def _attempt_plot(kw: dict, site_nm: bool) -> None:
                row.plot_occ_map(site_name=site_nm, **kw)

            def _plot_with_fallbacks(kw: dict) -> None:
                fallbacks = []
                seen = set()

                def _add_variant(d: dict) -> None:
                    sig = json.dumps(d, sort_keys=True, default=str)
                    if sig not in seen:
                        seen.add(sig)
                        fallbacks.append(dict(d))

                _add_variant(kw)
                _add_variant({k: v for k, v in kw.items() if k not in ("meridian", "parallels")})
                _add_variant(
                    {
                        k: v
                        for k, v in kw.items()
                        if k not in ("meridian", "parallels", "resolution", "mapstyle", "alpha")
                    }
                )
                _add_variant(
                    {
                        k: v
                        for k, v in kw.items()
                        if k
                        not in (
                            "meridian",
                            "parallels",
                            "resolution",
                            "mapstyle",
                            "alpha",
                            "site_box_alpha",
                        )
                    }
                )
                last_exc: Optional[BaseException] = None
                for kw_try in fallbacks:
                    for site_nm in (True, False):
                        try:
                            _attempt_plot(kw_try, site_nm)
                            return
                        except TypeError as te:
                            last_exc = te
                            try:
                                plt.close("all")
                            except Exception:
                                pass
                if last_exc is not None:
                    raise last_exc

            _plot_with_fallbacks(map_kw)
            png_path = out_dir / f"npoap_occ_map.{ext}"
            if not png_path.is_file():
                err = f"Fichier carte introuvable : {png_path}"
        except Exception as ex:
            logger.error("plot_occ_map : %s", ex, exc_info=True)
            err = str(ex)
        finally:
            self._pred_map_btn.configure(state=tk.NORMAL)

        if err:
            short = err if len(err) <= 220 else err[:220] + "\u2026"
            self._sora_pred_status.configure(text=f"Carte : erreur \u2014 {short}")
            logger.error("Carte : %s", err)
            return
        self._pred_save_map_prefs()
        self._pred_last_map_png_path = png_path
        self._sora_pred_status.configure(text=f"Carte générée : {png_path}")
        self._pred_map_placeholder.pack_forget()
        if png_path.suffix.lower() != ".png":
            self._pred_map_tk_photo = None
            self._pred_map_image.configure(image="")
            self._sora_pred_status.configure(
                text=f"Carte enregistrée ({png_path.suffix}) : {png_path} — aperçu intégré : PNG."
            )
            return
        self._pred_redraw_map_to_canvas()

    def _load_pymovie_settings(self) -> None:
        try:
            if not self._pymovie_settings_path.is_file():
                return
            with self._pymovie_settings_path.open("r", encoding="utf-8") as f:
                data = json.load(f)
            p = str(data.get("conda_env_root", "") or "").strip()
            if p:
                self._pymovie_conda_root_var.set(p)
        except Exception:
            pass

    def _save_pymovie_settings(self) -> None:
        try:
            payload = {"conda_env_root": (self._pymovie_conda_root_var.get() or "").strip()}
            with self._pymovie_settings_path.open("w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, ensure_ascii=False)
        except Exception:
            pass

    def _browse_pymovie_conda_root(self) -> None:
        cur = (self._pymovie_conda_root_var.get() or "").strip()
        initial = cur if cur and Path(cur).is_dir() else str(Path.home())
        d = filedialog.askdirectory(
            parent=self.winfo_toplevel(),
            title="Racine de l'environnement Conda (ex. …\\.conda\\envs\\pymovie)",
            initialdir=initial,
        )
        if d:
            self._pymovie_conda_root_var.set(d)
            self._save_pymovie_settings()
            if hasattr(self, "_pymovie_launch_status"):
                self._pymovie_launch_status.configure(text="Chemin enregistré.")

    def _pymovie_conda_python(self) -> Optional[Path]:
        root = Path((self._pymovie_conda_root_var.get() or "").strip())
        return _python_exe_in_conda_env(root)

    def _launch_pymovie_external(self) -> None:
        py = self._pymovie_conda_python()
        root = Path((self._pymovie_conda_root_var.get() or "").strip())
        if not root.is_dir():
            messagebox.showerror(
                "PyMovie",
                "Le répertoire racine de l'environnement Conda est invalide ou introuvable.\n"
                "Indiquez le dossier qui contient python.exe (ex. …\\.conda\\envs\\pymovie).",
                parent=self.winfo_toplevel(),
            )
            return
        if py is None:
            messagebox.showerror(
                "PyMovie",
                f"Aucun interpréteur Python trouvé dans :\n{root}\n\n"
                "Attendu : python.exe (Windows) ou bin/python (Linux, macOS).",
                parent=self.winfo_toplevel(),
            )
            return
        if self._pymovie_proc is not None and self._pymovie_proc.poll() is None:
            messagebox.showinfo("PyMovie", "PyMovie semble déjà lancé.", parent=self.winfo_toplevel())
            return
        try:
            self._pymovie_proc = subprocess.Popen(
                [str(py), "-c", "from pymovie import main; main.main()"],
                cwd=str(py.parent),
            )
            self._save_pymovie_settings()
            if hasattr(self, "_pymovie_launch_status"):
                self._pymovie_launch_status.configure(
                    text=f"PyMovie démarré (PID {self._pymovie_proc.pid})."
                )
        except Exception as ex:
            self._pymovie_proc = None
            logger.exception("Lancement PyMovie")
            messagebox.showerror(
                "PyMovie",
                f"Impossible de lancer PyMovie :\n{ex}",
                parent=self.winfo_toplevel(),
            )

    def _launch_pyote_external(self) -> None:
        py = self._pymovie_conda_python()
        root = Path((self._pymovie_conda_root_var.get() or "").strip())
        if not root.is_dir():
            messagebox.showerror(
                "PyOTE",
                "Le répertoire racine de l'environnement Conda est invalide ou introuvable.",
                parent=self.winfo_toplevel(),
            )
            return
        if py is None:
            messagebox.showerror(
                "PyOTE",
                f"Aucun interpréteur Python trouvé dans :\n{root}",
                parent=self.winfo_toplevel(),
            )
            return
        if self._pyote_proc is not None and self._pyote_proc.poll() is None:
            messagebox.showinfo("PyOTE", "PyOTE semble déjà lancé.", parent=self.winfo_toplevel())
            return
        try:
            self._pyote_proc = subprocess.Popen(
                [str(py), "-c", "from pyoteapp import pyote; pyote.main()"],
                cwd=str(py.parent),
            )
            self._save_pymovie_settings()
            if hasattr(self, "_pymovie_launch_status"):
                self._pymovie_launch_status.configure(
                    text=f"PyOTE démarré (PID {self._pyote_proc.pid})."
                )
        except Exception as ex:
            self._pyote_proc = None
            logger.exception("Lancement PyOTE")
            messagebox.showerror(
                "PyOTE",
                f"Impossible de lancer PyOTE :\n{ex}",
                parent=self.winfo_toplevel(),
            )

    def _run_pymovie_installer(self) -> None:
        """Ouvre une console qui exécute install_pymovie_env.bat (création env conda pymovie)."""
        bat = self._project_root / "install_pymovie_env.bat"
        if not bat.is_file():
            messagebox.showerror(
                "Installation PyMovie",
                f"Script introuvable :\n{bat}",
                parent=self.winfo_toplevel(),
            )
            return
        bat_s = str(bat.resolve())
        root_dir = str(self._project_root.resolve())
        try:
            # start : nouvelle fenêtre ; /k : laisser la console ouverte (logs, pause).
            subprocess.Popen(
                f'start "NPOAP - Installation PyMovie" cmd.exe /k call "{bat_s}"',
                cwd=root_dir,
                shell=True,
            )
        except Exception as ex:
            logger.exception("Lancement install_pymovie_env.bat")
            messagebox.showerror(
                "Installation PyMovie",
                f"Impossible de lancer l'installateur :\n{ex}",
                parent=self.winfo_toplevel(),
            )

    def _build_pymovie_panel(self, parent):
        section = ttk.LabelFrame(parent, text="PyMovie et PyOTE (IOTA)", padding=10)
        section.pack(fill=tk.BOTH, expand=True, anchor=tk.N)

        install_row = ttk.Frame(section)
        install_row.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(
            install_row,
            text="Installer l'environnement Conda (PyMovie + PyOTE)…",
            command=self._run_pymovie_installer,
        ).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Label(
            install_row,
            text="Ouvre une fenêtre console (conda create / pip). Journal : install_pymovie_env.log",
            wraplength=420,
        ).pack(side=tk.LEFT, fill=tk.X, expand=True)

        env_lf = ttk.LabelFrame(
            section,
            text="Environnement Conda",
            padding=8,
        )
        env_lf.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(
            env_lf,
            text="Racine du préfixe Conda (dossier contenant python.exe) :",
        ).pack(anchor="w")
        path_row = ttk.Frame(env_lf)
        path_row.pack(fill=tk.X, pady=(4, 6))
        entry = ttk.Entry(path_row, textvariable=self._pymovie_conda_root_var)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 6))
        entry.bind("<FocusOut>", lambda _e: self._save_pymovie_settings())
        ttk.Button(path_row, text="Parcourir…", command=self._browse_pymovie_conda_root).pack(
            side=tk.LEFT
        )
        btn_row = ttk.Frame(env_lf)
        btn_row.pack(fill=tk.X, pady=(0, 4))
        ttk.Button(btn_row, text="Lancer PyMovie", command=self._launch_pymovie_external).pack(
            side=tk.LEFT, padx=(0, 8)
        )
        ttk.Button(btn_row, text="Lancer PyOTE", command=self._launch_pyote_external).pack(
            side=tk.LEFT
        )
        self._pymovie_launch_status = ttk.Label(
            env_lf,
            text="Indiquez le chemin (ex. C:\\Users\\…\\.conda\\envs\\pymovie) puis lancez.",
            wraplength=640,
        )
        self._pymovie_launch_status.pack(anchor="w", pady=(4, 0))

        lf = ttk.LabelFrame(section, text="Liens utiles", padding=8)
        lf.pack(fill=tk.X, pady=(0, 0))
        self._link_row(
            lf,
            "Page IOTA \u2014 PyMovie",
            "https://occultations.org/observing/software/pymovie/",
        )
        self._link_row(
            lf,
            "Manuel PyMovie (PDF)",
            "https://occultations.org/sw/pymovie/PyMovie-doc.pdf",
        )
        self._link_row(
            lf,
            "Tutoriel rolling shutter (LiMovie + PYOTE, PDF)",
            "https://occultations.org/documents/LimovieTutorialRollingShutter.pdf",
        )

    def _run_sora_analysis_installer(self) -> None:
        """Ouvre une console qui exécute install_sora_analysis.bat (sora-astro dans astroenv)."""
        bat = self._project_root / "install_sora_analysis.bat"
        if not bat.is_file():
            messagebox.showerror(
                "Installation SORA",
                f"Script introuvable :\n{bat}",
                parent=self.winfo_toplevel(),
            )
            return
        bat_s = str(bat.resolve())
        root_dir = str(self._project_root.resolve())
        try:
            subprocess.Popen(
                f'start "NPOAP - Installation SORA (analyse)" cmd.exe /k call "{bat_s}"',
                cwd=root_dir,
                shell=True,
            )
        except Exception as ex:
            logger.exception("Lancement install_sora_analysis.bat")
            messagebox.showerror(
                "Installation SORA",
                f"Impossible de lancer l'installateur :\n{ex}",
                parent=self.winfo_toplevel(),
            )

    def _sora_analysis_set_status(self, text: str) -> None:
        if hasattr(self, "_sora_analysis_status"):
            self._sora_analysis_status.configure(text=text)

    def _sora_lc_pick_file(self) -> None:
        base = self._project_root / "output"
        path = filedialog.askopenfilename(
            title="Choisir un fichier de courbe de lumière (temps, flux)",
            initialdir=str(base if base.is_dir() else self._project_root),
            filetypes=[("Données", "*.dat *.txt *.csv"), ("Tous les fichiers", "*.*")],
        )
        if path:
            self._sora_lc_file_var.set(path)

    def _sora_lc_load_from_file(self) -> None:
        path = Path((self._sora_lc_file_var.get() or "").strip())
        if not path.is_file():
            self._sora_analysis_set_status("Fichier LC introuvable.")
            return
        exptime = self._pred_parse_float(self._sora_lc_exptime_var.get(), -1.0)
        if exptime <= 0:
            self._sora_analysis_set_status("Exptime invalide (doit être > 0).")
            return
        raw_cols = (self._sora_lc_usecols_var.get() or "0,1").replace(";", ",")
        try:
            usecols = [int(p.strip()) for p in raw_cols.split(",") if p.strip()]
        except ValueError:
            self._sora_analysis_set_status("usecols invalide (ex: 0,1).")
            return
        if len(usecols) != 2:
            self._sora_analysis_set_status("usecols doit contenir 2 colonnes: temps, flux.")
            return
        try:
            from sora import LightCurve

            lc = LightCurve(name=path.stem, file=str(path), exptime=float(exptime), usecols=usecols)
        except Exception as ex:
            logger.error("SORA LightCurve load: %s", ex, exc_info=True)
            self._sora_analysis_set_status(f"Chargement LC impossible: {ex}")
            return
        self._sora_lc_obj = lc
        self._sora_lc_fit = None
        self._sora_analysis_result.delete("1.0", tk.END)
        self._sora_analysis_result.insert(
            "1.0",
            (
                f"LC chargée: {path.name}\n"
                f"Initial: {getattr(lc, 'initial_time', '?')}\n"
                f"Mean:    {getattr(lc, 'time_mean', '?')}\n"
                f"End:     {getattr(lc, 'end_time', '?')}\n"
                f"Exptime: {getattr(lc, 'exptime', '?')} s\n"
                f"Cycle:   {getattr(lc, 'cycle', '?')} s\n"
            ),
        )
        self._sora_analysis_set_status("LC chargée. Lancez Auto-détection ou Ajustement.")

    def _sora_lc_autodetect(self) -> None:
        lc = getattr(self, "_sora_lc_obj", None)
        if lc is None:
            self._sora_analysis_set_status("Chargez d'abord une courbe de lumière.")
            return
        self._sora_analysis_set_status("Auto-détection en cours (occ_detect)…")

        def work():
            try:
                d = lc.occ_detect()
            except Exception as ex:
                logger.error("SORA occ_detect: %s", ex, exc_info=True)
                self.after(0, lambda: self._sora_analysis_set_status(f"Auto-détection échouée: {ex}"))
                return

            def done():
                self._sora_lc_detect = d
                imm = float(d.get("immersion_time"))
                eme = float(d.get("emersion_time"))
                terr = float(d.get("time_err", 0.05))
                self._sora_lc_imm_var.set(f"{imm:.6f}")
                self._sora_lc_eme_var.set(f"{eme:.6f}")
                self._sora_lc_delta_var.set(f"{max(terr * 5.0, 0.001):.6f}")
                self._sora_lc_tmin_var.set(f"{imm - 5.0:.6f}")
                self._sora_lc_tmax_var.set(f"{eme + 5.0:.6f}")
                self._sora_analysis_result.insert(
                    tk.END,
                    (
                        "\nAuto-détection:\n"
                        + "\n".join(f"  {k}: {v}" for k, v in d.items() if k != "occ_mask")
                        + "\n"
                    ),
                )
                self._sora_analysis_set_status("Auto-détection OK. Vérifiez les bornes puis lancez l'ajustement.")

            self.after(0, done)

        threading.Thread(target=work, daemon=True).start()

    def _sora_lc_run_fit(self) -> None:
        lc = getattr(self, "_sora_lc_obj", None)
        if lc is None:
            self._sora_analysis_set_status("Chargez d'abord une courbe de lumière.")
            return
        vel = self._pred_parse_float(self._sora_lc_vel_var.get(), 22.0)
        dist = self._pred_parse_float(self._sora_lc_dist_var.get(), 15.0)
        d_star = self._pred_parse_float(self._sora_lc_dstar_var.get(), 0.2)
        imm = self._pred_parse_float(self._sora_lc_imm_var.get(), float("nan"))
        eme = self._pred_parse_float(self._sora_lc_eme_var.get(), float("nan"))
        delta_t = max(self._pred_parse_float(self._sora_lc_delta_var.get(), 0.1), 1e-6)
        tmin = self._pred_parse_float(self._sora_lc_tmin_var.get(), imm - 5.0)
        tmax = self._pred_parse_float(self._sora_lc_tmax_var.get(), eme + 5.0)
        flux_min = self._pred_parse_float(self._sora_lc_fmin_var.get(), 0.0)
        flux_max = self._pred_parse_float(self._sora_lc_fmax_var.get(), 1.0)
        loop = max(100, self._pred_parse_int(self._sora_lc_loop_var.get(), 2000))
        threads = max(1, self._pred_parse_int(self._sora_lc_threads_var.get(), 4))
        method = (self._sora_lc_method_var.get() or "ls").strip().lower()
        if method not in ("ls", "de", "chisqr", "fastchi"):
            method = "ls"
        self._sora_analysis_set_status("Ajustement SORA en cours (occ_lcfit)…")

        def work():
            try:
                lc.set_vel(vel=float(vel))
                lc.set_dist(dist=float(dist))
                lc.set_star_diam(d_star=float(d_star))
                fit = lc.occ_lcfit(
                    immersion_time=float(imm),
                    emersion_time=float(eme),
                    delta_t=float(delta_t),
                    flux_min=float(flux_min),
                    flux_max=float(flux_max),
                    tmin=float(tmin),
                    tmax=float(tmax),
                    sigma="auto",
                    sigma_result=1,
                    loop=int(loop),
                    method=method,
                    threads=int(threads),
                )
            except Exception as ex:
                logger.error("SORA occ_lcfit: %s", ex, exc_info=True)
                self.after(0, lambda: self._sora_analysis_set_status(f"Ajustement échoué: {ex}"))
                return

            def done():
                self._sora_lc_fit = fit
                self._sora_analysis_result.insert(tk.END, "\nAjustement:\n" + str(fit) + "\n")
                self._sora_analysis_result.see(tk.END)
                self._sora_analysis_set_status("Ajustement terminé. Bouton χ² disponible.")

            self.after(0, done)

        threading.Thread(target=work, daemon=True).start()

    def _sora_lc_plot_chi2(self) -> None:
        fit = getattr(self, "_sora_lc_fit", None)
        if fit is None:
            self._sora_analysis_set_status("Aucun fit disponible. Lancez d'abord occ_lcfit.")
            return
        try:
            fit.plot_chi2()
            self._sora_analysis_set_status("Graphiques χ² générés.")
        except Exception as ex:
            logger.error("SORA plot_chi2: %s", ex, exc_info=True)
            self._sora_analysis_set_status(f"Affichage χ² échoué: {ex}")

    def _build_sora_analysis_panel(self, parent):
        section = ttk.LabelFrame(parent, text="Analyse SORA (IOTA)", padding=10)
        section.pack(fill=tk.BOTH, expand=True, anchor=tk.N)

        install_row = ttk.Frame(section)
        install_row.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(
            install_row,
            text="Installer / mettre à jour sora-astro (environnement astroenv)…",
            command=self._run_sora_analysis_installer,
        ).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Label(
            install_row,
            text="Ouvre une console (pip dans conda astroenv). Journal : install_sora_analysis.log",
            wraplength=420,
        ).pack(side=tk.LEFT, fill=tk.X, expand=True)

        note = (
            "Périmètre implémenté : analyse de courbe photométrique (LightCurve fichier, paramètres physiques, "
            "occ_detect, occ_lcfit, méthodes/threads, résultats et χ²).\n"
            "Hors périmètre pour l’instant : cordes multiples, modèle 3D OBJ, texture, fit_shape, "
            "boucles d’orientation."
        )
        self._readonly_text(section, note, height=4)

        lc_lf = ttk.LabelFrame(section, text="Flux LightCurve (SORA)", padding=8)
        lc_lf.pack(fill=tk.X, pady=(8, 6))

        self._sora_lc_obj = None
        self._sora_lc_fit = None
        self._sora_lc_detect = None
        self._sora_lc_file_var = tk.StringVar(value="")
        self._sora_lc_exptime_var = tk.StringVar(value="0.1")
        self._sora_lc_usecols_var = tk.StringVar(value="0,1")
        self._sora_lc_vel_var = tk.StringVar(value="22.0")
        self._sora_lc_dist_var = tk.StringVar(value="15.0")
        self._sora_lc_dstar_var = tk.StringVar(value="0.2")
        self._sora_lc_imm_var = tk.StringVar(value="")
        self._sora_lc_eme_var = tk.StringVar(value="")
        self._sora_lc_delta_var = tk.StringVar(value="0.1")
        self._sora_lc_tmin_var = tk.StringVar(value="")
        self._sora_lc_tmax_var = tk.StringVar(value="")
        self._sora_lc_fmin_var = tk.StringVar(value="0.0")
        self._sora_lc_fmax_var = tk.StringVar(value="1.0")
        self._sora_lc_loop_var = tk.StringVar(value="2000")
        self._sora_lc_threads_var = tk.StringVar(value="4")
        self._sora_lc_method_var = tk.StringVar(value="ls")

        r = 0
        ttk.Label(lc_lf, text="Fichier LC (temps, flux)").grid(row=r, column=0, sticky="w", padx=(0, 6), pady=2)
        ttk.Entry(lc_lf, textvariable=self._sora_lc_file_var).grid(row=r, column=1, sticky="ew", pady=2)
        ttk.Button(lc_lf, text="Parcourir…", command=self._sora_lc_pick_file).grid(
            row=r, column=2, padx=(6, 0), pady=2
        )
        r += 1
        ttk.Label(lc_lf, text="Exptime (s)").grid(row=r, column=0, sticky="w", padx=(0, 6), pady=2)
        ttk.Entry(lc_lf, textvariable=self._sora_lc_exptime_var, width=10).grid(row=r, column=1, sticky="w", pady=2)
        ttk.Label(lc_lf, text="usecols (temps,flux)").grid(row=r, column=1, sticky="e", padx=(0, 84), pady=2)
        ttk.Entry(lc_lf, textvariable=self._sora_lc_usecols_var, width=10).grid(row=r, column=2, sticky="w", pady=2)
        r += 1
        ttk.Button(lc_lf, text="Charger LightCurve", command=self._sora_lc_load_from_file).grid(
            row=r, column=0, sticky="w", pady=(4, 2)
        )
        ttk.Button(lc_lf, text="Auto-détection (occ_detect)", command=self._sora_lc_autodetect).grid(
            row=r, column=1, sticky="w", pady=(4, 2)
        )
        r += 1

        phys = ttk.LabelFrame(lc_lf, text="Physique / Fit", padding=6)
        phys.grid(row=r, column=0, columnspan=3, sticky="ew", pady=(6, 2))
        rr = 0
        ttk.Label(phys, text="vel (km/s)").grid(row=rr, column=0, sticky="w", padx=(0, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_vel_var, width=8).grid(row=rr, column=1, sticky="w", pady=1)
        ttk.Label(phys, text="dist (AU)").grid(row=rr, column=2, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_dist_var, width=8).grid(row=rr, column=3, sticky="w", pady=1)
        ttk.Label(phys, text="d_star (km)").grid(row=rr, column=4, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_dstar_var, width=8).grid(row=rr, column=5, sticky="w", pady=1)
        rr += 1
        ttk.Label(phys, text="immersion").grid(row=rr, column=0, sticky="w", padx=(0, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_imm_var, width=12).grid(row=rr, column=1, sticky="w", pady=1)
        ttk.Label(phys, text="emersion").grid(row=rr, column=2, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_eme_var, width=12).grid(row=rr, column=3, sticky="w", pady=1)
        ttk.Label(phys, text="delta_t").grid(row=rr, column=4, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_delta_var, width=8).grid(row=rr, column=5, sticky="w", pady=1)
        rr += 1
        ttk.Label(phys, text="tmin").grid(row=rr, column=0, sticky="w", padx=(0, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_tmin_var, width=12).grid(row=rr, column=1, sticky="w", pady=1)
        ttk.Label(phys, text="tmax").grid(row=rr, column=2, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_tmax_var, width=12).grid(row=rr, column=3, sticky="w", pady=1)
        ttk.Label(phys, text="flux min/max").grid(row=rr, column=4, sticky="w", padx=(12, 4), pady=1)
        ff = ttk.Frame(phys)
        ff.grid(row=rr, column=5, sticky="w", pady=1)
        ttk.Entry(ff, textvariable=self._sora_lc_fmin_var, width=5).pack(side=tk.LEFT)
        ttk.Label(ff, text="/").pack(side=tk.LEFT, padx=2)
        ttk.Entry(ff, textvariable=self._sora_lc_fmax_var, width=5).pack(side=tk.LEFT)
        rr += 1
        ttk.Label(phys, text="Méthode").grid(row=rr, column=0, sticky="w", padx=(0, 4), pady=1)
        ttk.Combobox(
            phys,
            textvariable=self._sora_lc_method_var,
            values=("ls", "de", "chisqr", "fastchi"),
            width=10,
            state="readonly",
        ).grid(row=rr, column=1, sticky="w", pady=1)
        ttk.Label(phys, text="loops").grid(row=rr, column=2, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_loop_var, width=10).grid(row=rr, column=3, sticky="w", pady=1)
        ttk.Label(phys, text="threads").grid(row=rr, column=4, sticky="w", padx=(12, 4), pady=1)
        ttk.Entry(phys, textvariable=self._sora_lc_threads_var, width=10).grid(row=rr, column=5, sticky="w", pady=1)
        rr += 1
        ttk.Button(phys, text="Lancer occ_lcfit", command=self._sora_lc_run_fit).grid(
            row=rr, column=0, sticky="w", pady=(6, 0)
        )
        ttk.Button(phys, text="Afficher χ²", command=self._sora_lc_plot_chi2).grid(
            row=rr, column=1, sticky="w", pady=(6, 0)
        )
        ttk.Label(
            phys,
            text="Conseil SORA: threads <= coeurs CPU - 1 pour éviter de figer la machine.",
            font=("TkDefaultFont", 8),
        ).grid(row=rr, column=2, columnspan=4, sticky="w", pady=(6, 0))
        lc_lf.columnconfigure(1, weight=1)

        self._sora_analysis_status = ttk.Label(section, text="", wraplength=980)
        self._sora_analysis_status.pack(anchor="w", pady=(2, 0))

        rf = ttk.LabelFrame(section, text="Résultats (Auto-détection + Fit)", padding=4)
        rf.pack(fill=tk.BOTH, expand=True, pady=(6, 0))
        rw = ttk.Frame(rf)
        rw.pack(fill=tk.BOTH, expand=True)
        self._sora_analysis_result = tk.Text(
            rw, height=10, wrap=tk.NONE, font=("Consolas", 8), relief=tk.FLAT, padx=4, pady=4
        )
        rsb = ttk.Scrollbar(rw, command=self._sora_analysis_result.yview)
        self._sora_analysis_result.configure(yscrollcommand=rsb.set)
        self._sora_analysis_result.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        rsb.pack(side=tk.RIGHT, fill=tk.Y)
        self._sora_analysis_result.insert(
            "1.0",
            "Workflow: 1) Choisir fichier LC  2) Charger  3) Auto-détection  4) Ajuster occ_lcfit  5) χ².",
        )
        self._sora_analysis_set_status("Prêt pour l'analyse LightCurve SORA (sans couche 3D).")

        lf = ttk.LabelFrame(section, text="Liens utiles", padding=8)
        lf.pack(fill=tk.X, pady=(8, 0))
        self._link_row(lf, "Dépôt GitHub SORA", "https://github.com/riogroup/SORA")
        self._link_row(lf, "Documentation ReadTheDocs", "https://sora.readthedocs.io/")
        self._link_row(
            lf,
            "Article MNRAS (DOI)",
            "https://doi.org/10.1093/mnras/stac032",
        )
