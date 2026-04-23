# gui/occultations_tab.py
"""Onglet Occultations : SORA (prédictions, analyse), PyMovie."""

import logging
import threading
import tkinter as tk
from datetime import timedelta, timezone
from pathlib import Path
from typing import Optional
from tkinter import filedialog, ttk
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
    # Zone fixe pour la carte (tk + subsampling du PNG) ; ratio ~ 4:3
    _MAP_VIEW_W = 960
    _MAP_VIEW_H = 720

    def __init__(self, parent, night_observation_tab=None):
        super().__init__(parent, padding=8)
        self._project_root = Path(__file__).resolve().parent.parent
        self._night_observation_tab = night_observation_tab
        self._pred_occ_labels = []
        self._pred_occ_rows = []
        self._pred_map_tk_photo = None

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
            minsize=self._MAP_VIEW_W + 24,
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

        map_canvas_wrap = tk.Frame(
            map_lf,
            width=self._MAP_VIEW_W,
            height=self._MAP_VIEW_H,
            bg=self._text_bg(),
        )
        map_canvas_wrap.pack(anchor="center", pady=(4, 0))
        map_canvas_wrap.pack_propagate(False)
        self._pred_map_placeholder = ttk.Label(
            map_canvas_wrap,
            text=(
                "Lancez une prédiction, sélectionnez une occultation ci-dessus, "
                "puis \u00ab Voir la carte \u00bb. La figure est générée par SORA (matplotlib / cartopy)."
            ),
            wraplength=max(320, self._MAP_VIEW_W - 80),
            justify=tk.LEFT,
        )
        self._pred_map_placeholder.pack(anchor="nw", pady=(0, 6), padx=4)
        self._pred_map_image = tk.Label(map_canvas_wrap, text="", bg=self._text_bg())
        self._pred_map_image.pack(anchor="center")

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
        if not hasattr(self, "_pred_map_image"):
            return
        self._pred_map_image.configure(image="")
        if self._pred_map_placeholder.winfo_manager() == "":
            self._pred_map_placeholder.pack(anchor="nw", pady=(0, 6), before=self._pred_map_image)

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
            mapsize_cm = [42.0, 44.0]
            logger.info(
                "plot_occ_map : site=%s lon=%.6f lat=%.6f, zoom=4, mapsize=%s cm, arrow=False",
                obs_name, obs_lon, obs_lat, mapsize_cm,
            )
            map_kw = dict(
                path=str(out_dir),
                nameimg="npoap_occ_map",
                fmt="png",
                dpi=120,
                mapsize=mapsize_cm,
                sites=site_dict,
                centerproj=[obs_lon, obs_lat],
                centermap_geo=[obs_lon, obs_lat],
                zoom=4,
                arrow=False,
                sscale=2.4,
                site_box_alpha=0.88,
            )
            if err_mas is not None:
                map_kw["error"] = err_mas
                map_kw["ercolor"] = "darkorange"
            try:
                row.plot_occ_map(site_name=True, **map_kw)
            except TypeError:
                logger.info(
                    "site_name=True : TypeError (SORA/cartopy) ; nouvelle tentative sans libellé de site."
                )
                plt.close()
                row.plot_occ_map(site_name=False, **map_kw)
            png_path = out_dir / "npoap_occ_map.png"
            if not png_path.is_file():
                err = f"Fichier PNG introuvable : {png_path}"
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
        self._sora_pred_status.configure(text=f"Carte générée : {png_path}")
        self._pred_map_placeholder.pack_forget()
        try:
            photo = tk.PhotoImage(file=str(png_path))
            mw, mh = self._MAP_VIEW_W, self._MAP_VIEW_H
            while photo.width() > mw or photo.height() > mh:
                photo = photo.subsample(2, 2)
        except tk.TclError as ex2:
            logger.error("Lecture PNG pour affichage Tkinter : %s", ex2)
            self._sora_pred_status.configure(text=f"Carte : impossible de charger le PNG ({ex2}).")
            return
        self._pred_map_tk_photo = photo
        self._pred_map_image.configure(image=photo)

    def _build_pymovie_panel(self, parent):
        ttk.Label(
            parent,
            text="PyMovie et PyOTE (IOTA)",
            font=("Helvetica", 12, "bold"),
        ).pack(anchor="w", pady=(0, 6))

        intro = (
            "PyMovie : photométrie par ouverture sur vidéos d'occultations stellaires.\n"
            "PyOTE : analyse des temps / diffraction en Python, conçu pour compléter PyMovie.\n\n"
            "Installez les paquets fournis par l'IOTA (Windows / Mac) et consultez le manuel PDF "
            "et les tutoriels vidéo sur le site IOTA."
        )
        self._readonly_text(parent, intro, height=7)

        lf = ttk.LabelFrame(parent, text="Liens utiles", padding=8)
        lf.pack(fill=tk.X, pady=(10, 0))
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

    def _build_sora_analysis_panel(self, parent):
        ttk.Label(
            parent,
            text="SORA \u2014 analyse (courbes, cordes, occultation)",
            font=("Helvetica", 12, "bold"),
        ).pack(anchor="w", pady=(0, 6))

        intro = (
            "Après acquisition et photométrie, SORA permet la modélisation : LightCurve, cordes, "
            "objet Occultation, etc.\n\n"
            "SORA est installé avec le reste des paquets (requirements.txt) dans astroenv. "
            "En dépannage Cartopy / GEOS, voir la doc SORA ou conda-forge.\n\n"
            "Citation : Gomes-Júnior et al. (2022), MNRAS 511, 1167\u20131181, DOI 10.1093/mnras/stac032"
        )
        self._readonly_text(parent, intro, height=9)

        lf = ttk.LabelFrame(parent, text="Liens utiles", padding=8)
        lf.pack(fill=tk.X, pady=(10, 0))
        self._link_row(lf, "Dépôt GitHub SORA", "https://github.com/riogroup/SORA")
        self._link_row(lf, "Documentation ReadTheDocs", "https://sora.readthedocs.io/")
        self._link_row(
            lf,
            "Article MNRAS (DOI)",
            "https://doi.org/10.1093/mnras/stac032",
        )
