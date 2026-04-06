import os
import shutil
import subprocess
import socket
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import webbrowser
from datetime import datetime
from typing import List, Literal, Optional, Tuple
from zoneinfo import ZoneInfo

import config

# Séquence type « exemple officiel » C2A (aide TCP). Champ 1°×1°. Date/heure : chaîne locale
# calculée depuis le fuseau de l’observatoire dans config.json (onglet Accueil), pas l’horloge Windows seule.
# https://c2a.astrosurf.com/english/support/help/Html/Comment_communiquer_avec_C2A_au_travers_d_un_socket_TCP.htm
C2A_TCP_DEFAULT_TEMPLATE = (
    "SetRA={ra_hours:.6f};SetDE={dec_deg:.6f};"
    "SetFieldX=1;SetFieldY=1;"
    "SetDateTime={obs_c2a_local};"
    "SetLatitude={lat_obs:.6f};SetLongitude={lon_obs:.6f};"
    "SetMapType=field;"
)

# Anciens modèles : migration vers le défaut actuel.
_C2A_TCP_LEGACY_TEMPLATES = frozenset(
    {
        "SetRa={ra_hours:.6f};SetDe={dec_deg:.6f};SetMapType=field;",
        "SetRA={ra_hours:.6f};SetDE={dec_deg:.6f};SetMapType=field;",
        "CENTER {ra_deg:.6f} {dec_deg:.6f}",
        "SetRA={ra_hours:.6f};SetDE={dec_deg:.6f};SetFieldX=2;SetFieldY=2;SetDateTime=currentlocal;SetLatitude={lat_obs:.6f};SetLongitude={lon_obs:.6f};SetMapType=field;",
    }
)


def _config_json_path() -> Path:
    return Path(getattr(config, "BASE_DIR", Path(__file__).resolve().parents[1])) / "config.json"


def _effective_observatory() -> dict:
    """
    Observatoire effectif : config.json (bloc « observatory », sauvé par l’onglet Accueil)
    prime sur les clés présentes ; complété par config.OBSERVATORY sinon.
    """
    out: dict = {}
    base = getattr(config, "OBSERVATORY", None) or {}
    out.update(dict(base))
    try:
        p = _config_json_path()
        if p.is_file():
            with p.open("r", encoding="utf-8") as f:
                data = json.load(f)
            job = data.get("observatory")
            if isinstance(job, dict):
                out.update(job)
    except Exception:
        pass
    return out


def _obs_latitude(obs: dict) -> float:
    for k in ("latitude", "lat"):
        if k not in obs:
            continue
        v = obs.get(k)
        if v is None or (isinstance(v, str) and not str(v).strip()):
            continue
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return 0.0


def _obs_longitude(obs: dict) -> float:
    for k in ("longitude", "lon"):
        if k not in obs:
            continue
        v = obs.get(k)
        if v is None or (isinstance(v, str) and not str(v).strip()):
            continue
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return 0.0


def _obs_elevation(obs: dict) -> float:
    for k in ("elevation", "elev"):
        if k not in obs:
            continue
        v = obs.get(k)
        if v is None or (isinstance(v, str) and not str(v).strip()):
            continue
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return 0.0


def _observatory_iana_tz_name() -> str:
    """Dérive un fuseau IANA depuis l’observatoire effectif (config.json + config.py)."""
    obs = _effective_observatory()
    tz_name = str(obs.get("timezone", "") or "").strip()
    if not tz_name:
        return "UTC"
    low = tz_name.lower()
    if "santiago" in low or "chili" in low or "chile" in low:
        return "America/Santiago"
    if "/" in tz_name:
        return tz_name
    return "UTC"


def observatory_local_now_c2a_string() -> str:
    """
    Date et heure « locales observatoire » au format fixe C2A : DD/MM/AAAA HH:MM:SS
    (heure civile dans le fuseau défini dans config.json → observatory.timezone).
    """
    tzid = _observatory_iana_tz_name()
    try:
        tz = ZoneInfo(tzid)
    except Exception:
        tz = ZoneInfo("UTC")
    return datetime.now(tz).strftime("%d/%m/%Y %H:%M:%S")


class PlanetariumTab(ttk.Frame):
    """
    Onglet Planétarium.

    Objectif principal : fournir un accès direct au logiciel C2A
    (voir documentation officielle : https://c2a.astrosurf.com/).
    """

    def __init__(self, parent):
        super().__init__(parent, padding=10)

        self.c2a_process = None
        self.c2a_path_var = tk.StringVar()
        # Liste locale des cibles envoyées par l’onglet « Observation de la Nuit »
        self.targets = []
        # Évite un double envoi TCP quand on programme la sélection après remplissage de la liste
        self._tcp_select_suppress = False
        # Fichier de configuration locale pour les paramètres TCP/IP (~/.npoap, comme l’onglet Nuit)
        self._npoap_dir = Path.home() / ".npoap"
        try:
            self._npoap_dir.mkdir(exist_ok=True)
        except Exception:
            pass
        # TCP + chemin exécutable C2A (C2A.exe, C2aw.exe, etc.)
        self._tcp_settings_path = self._npoap_dir / "planetarium_c2a_tcp.json"

        # Essai d'auto‑détection de C2A ; le chemin enregistré (.npoap) prime si présent
        self._auto_detect_c2a()

        self._build_ui()
        self._load_tcp_settings()

    def _refresh_observatory_info_label(self) -> None:
        """Met à jour le résumé lieu/fuseau depuis config.json (prioritaire) et config.py."""
        if not getattr(self, "_observatory_info_label", None):
            return
        obs = _effective_observatory()
        lat = _obs_latitude(obs)
        lon = _obs_longitude(obs)
        elev = _obs_elevation(obs)
        tz = str(obs.get("timezone", "") or "UTC")
        name = str(obs.get("name", "") or "N/A")
        cfg_path = _config_json_path()
        src = f"Fichier : {cfg_path}" if cfg_path.is_file() else f"{cfg_path} (absent — valeurs config.py)"
        self._observatory_info_label.configure(
            text=(
                f"Nom : {name}\n"
                f"Latitude : {lat:.4f}°   Longitude : {lon:.4f}°   Altitude : {elev:.0f} m\n"
                f"Fuseau horaire : {tz}\n"
                f"{src}\n"
                "Modifiez ces valeurs dans l’onglet Accueil (enregistrement dans config.json), "
                "puis reportez-les si besoin dans C2A (Options → Location)."
            )
        )

    # ------------------------------------------------------------------
    # Construction de l'interface
    # ------------------------------------------------------------------
    def _build_ui(self):
        title = ttk.Label(
            self,
            text="Planétarium – Intégration C2A",
            font=("Helvetica", 14, "bold"),
        )
        title.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 10))

        desc = ttk.Label(
            self,
            text=(
                "Cet onglet permet de lancer le logiciel de planétarium C2A.\n"
                "Assurez‑vous que C2A est installé sur votre machine "
                "(voir la page de support / version history du site officiel)."
            ),
            justify="left",
            wraplength=900,
        )
        desc.grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 15))

        # Ligne sélection du chemin vers C2A
        ttk.Label(self, text="Chemin de l'exécutable C2A :").grid(
            row=2, column=0, sticky="w"
        )

        entry = ttk.Entry(self, textvariable=self.c2a_path_var, width=90)
        entry.grid(row=2, column=1, sticky="we", padx=(5, 5))

        browse_btn = ttk.Button(self, text="Parcourir…", command=self._browse_c2a)
        browse_btn.grid(row=2, column=2, sticky="w")

        ttk.Label(
            self,
            text=(
                "Le chemin (C2A.exe, C2aw.exe, …) est mémorisé avec les autres réglages C2A dans "
                f"{self._npoap_dir}."
            ),
            foreground="gray45",
            font=("TkDefaultFont", 8),
            wraplength=880,
            justify="left",
        ).grid(row=3, column=0, columnspan=3, sticky="w", pady=(2, 0))

        # Ligne boutons Installer / Lancer / Quitter
        btn_row = ttk.Frame(self)
        btn_row.grid(row=4, column=0, columnspan=3, sticky="w", pady=(10, 0))

        ttk.Button(
            btn_row,
            text="Installer / Réinstaller C2A",
            command=self._install_c2a,
        ).pack(side=tk.LEFT, padx=(0, 8))

        ttk.Button(
            btn_row,
            text="Lancer C2A",
            command=self._launch_c2a,
        ).pack(side=tk.LEFT, padx=(0, 8))

        ttk.Button(
            btn_row,
            text="Quitter C2A",
            command=self._quit_c2a,
        ).pack(side=tk.LEFT)

        # Informations observatoire (config.json + repli config.py)
        obs_frame = ttk.LabelFrame(self, text="Localisation de l'observatoire (config.json)", padding=6)
        obs_frame.grid(row=5, column=0, columnspan=3, sticky="we", pady=(15, 0))

        self._observatory_info_label = ttk.Label(
            obs_frame,
            text="",
            justify="left",
            wraplength=900,
        )
        self._observatory_info_label.pack(anchor="w")
        self._refresh_observatory_info_label()

        # Cadre liste des cibles reçues depuis l’onglet Nuit
        targets_frame = ttk.LabelFrame(self, text="Cibles envoyées depuis « Observation de la Nuit »", padding=6)
        targets_frame.grid(row=6, column=0, columnspan=3, sticky="nsew", pady=(15, 0))

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=0)
        self.rowconfigure(6, weight=1)

        cols = ("Nom", "Type", "RA", "Dec", "Mag vis")
        self.targets_tree = ttk.Treeview(
            targets_frame,
            columns=cols,
            show="headings",
            height=8,
        )
        for col in cols:
            self.targets_tree.heading(col, text=col)
        self.targets_tree.column("Nom", width=160, anchor="w")
        self.targets_tree.column("Type", width=90, anchor="center")
        self.targets_tree.column("RA", width=120, anchor="center")
        self.targets_tree.column("Dec", width=120, anchor="center")
        self.targets_tree.column("Mag vis", width=70, anchor="center")
        self.targets_tree.pack(fill=tk.BOTH, expand=True)

        tgt_btn_row = ttk.Frame(targets_frame)
        tgt_btn_row.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(
            tgt_btn_row,
            text="Pointer C2A sur la cible sélectionnée",
            command=self._on_point_c2a_button,
        ).pack(side=tk.LEFT)

        # Mapping item Treeview → objet d'éphémérides
        self._targets_by_iid = {}

        # Contrôle TCP/IP C2A (protocole officiel : port 5876, commandes terminées par « ; »)
        tcp_frame = ttk.LabelFrame(self, text="Contrôle TCP/IP C2A (socket, port 5876)", padding=6)
        tcp_frame.grid(row=7, column=0, columnspan=3, sticky="we", pady=(8, 0))

        self.tcp_enabled_var = tk.BooleanVar(value=False)
        self.tcp_host_var = tk.StringVar(value="127.0.0.1")
        self.tcp_port_var = tk.StringVar(value="5876")
        # Doc : SetRA en heures décimales, SetDE en degrés ; plusieurs commandes séparées par « ; »
        self.tcp_command_template_var = tk.StringVar(value=C2A_TCP_DEFAULT_TEMPLATE)

        ttk.Checkbutton(
            tcp_frame,
            text=(
                "Activer le contrôle TCP vers C2A (sélection dans la liste, envoi depuis la Nuit, "
                "bouton « Pointer C2A… »)"
            ),
            variable=self.tcp_enabled_var,
        ).grid(row=0, column=0, columnspan=3, sticky="w")

        ttk.Label(
            tcp_frame,
            text=(
                "C2A écoute sur le port 5876. Chaque commande se termine par « ; » ; "
                "les nombres utilisent le point décimal. Pour que la carte se mette à jour, "
                "mimiquez l’exemple du manuel : champ (SetFieldX/Y), date (SetDateTime), "
                "lieu (SetLatitude/Longitude), type field, puis SetRA / SetDE."
            ),
            foreground="gray35",
            wraplength=880,
            justify="left",
        ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(0, 4))

        ttk.Label(tcp_frame, text="Hôte:").grid(row=2, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(tcp_frame, textvariable=self.tcp_host_var, width=16).grid(row=2, column=1, sticky="w", pady=(4, 0))

        ttk.Label(tcp_frame, text="Port:").grid(row=2, column=2, sticky="w", padx=(8, 0), pady=(4, 0))
        ttk.Entry(tcp_frame, textvariable=self.tcp_port_var, width=8).grid(row=2, column=3, sticky="w", pady=(4, 0))

        ttk.Label(
            tcp_frame,
            text=(
                "Séquence (format Python) : {name}, {ra_deg}, {dec_deg}, {ra_hours}, {lat_obs}, {lon_obs}, "
                "{obs_c2a_local} (fuseau = observatory.timezone dans config.json) :"
            ),
        ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(6, 0))
        ttk.Entry(
            tcp_frame,
            textvariable=self.tcp_command_template_var,
            width=80,
        ).grid(row=4, column=0, columnspan=4, sticky="we", pady=(2, 0))

        for i in range(4):
            tcp_frame.columnconfigure(i, weight=1)

        # Sélection dans la liste → éventuelle commande TCP
        self.targets_tree.bind("<<TreeviewSelect>>", self._on_target_selected)

        # Sauvegarde automatique des paramètres TCP à chaque modification
        def _on_var_changed(*_):
            self._save_tcp_settings_safe()

        self.tcp_enabled_var.trace_add("write", _on_var_changed)
        self.tcp_host_var.trace_add("write", _on_var_changed)
        self.tcp_port_var.trace_add("write", _on_var_changed)
        self.tcp_command_template_var.trace_add("write", _on_var_changed)
        self.c2a_path_var.trace_add("write", _on_var_changed)

        help_lbl = ttk.Label(
            self,
            text=(
                "Les cibles de l’onglet « Observation de la Nuit » apparaissent ici. "
                "Après envoi, la première cible est sélectionnée et une commande TCP est envoyée "
                "si le contrôle TCP est activé (C2A doit être ouvert, port 5876).\n"
                "Réf. protocole : "
                "https://c2a.astrosurf.com/english/support/help/Html/"
                "Comment_communiquer_avec_C2A_au_travers_d_un_socket_TCP.htm"
            ),
            foreground="gray40",
            wraplength=900,
            justify="left",
        )
        help_lbl.grid(row=8, column=0, columnspan=3, sticky="w", pady=(4, 0))

    # ------------------------------------------------------------------
    # Utilitaires
    # ------------------------------------------------------------------
    def _auto_detect_c2a(self):
        """
        Tente de localiser C2A dans quelques emplacements classiques.
        L'utilisateur peut toujours surcharger via le bouton 'Parcourir…'.
        """
        candidates = [
            r"C:\Program Files\C2A\C2A.exe",
            r"C:\Program Files (x86)\C2A\C2A.exe",
        ]

        for path in candidates:
            if os.path.isfile(path):
                self.c2a_path_var.set(path)
                return

        # Dernier recours : chercher 'C2A.exe' dans le PATH
        exe_name = "C2A.exe"
        found = shutil.which(exe_name)
        if found:
            self.c2a_path_var.set(found)

    def _browse_c2a(self):
        """
        Ouvre un sélecteur de fichier pour choisir C2A.exe.
        """
        initial_dir = None
        current = self.c2a_path_var.get().strip()
        if current:
            p = Path(current)
            if p.exists():
                initial_dir = str(p.parent)

        filepath = filedialog.askopenfilename(
            title="Sélectionner l'exécutable C2A",
            filetypes=[("Fichiers exécutables", "*.exe"), ("Tous les fichiers", "*.*")],
            initialdir=initial_dir or "C:\\",
        )
        if filepath:
            self.c2a_path_var.set(filepath)

    def _install_c2a(self):
        """
        Installation / configuration de C2A.

        1) Si une distribution locale est trouvée sous external_apps, on propose
           de lancer directement l’installateur C2A.
        2) À défaut, on ouvre la page officielle dans le navigateur.
        """
        base_dir = Path(config.BASE_DIR) if hasattr(config, "BASE_DIR") else Path(__file__).resolve().parents[1]
        external_apps = base_dir / "external_apps"
        local_installer = None
        if external_apps.exists():
            # Cherche un exécutable C2A dans external_apps (distribution déposée par l’utilisateur).
            for root, _, files in os.walk(external_apps):
                for f in files:
                    if f.lower().startswith("c2a") and f.lower().endswith(".exe"):
                        local_installer = Path(root) / f
                        break
                if local_installer is not None:
                    break

        if local_installer is not None:
            answer = messagebox.askyesno(
                "C2A",
                f"Une distribution C2A a été trouvée :\n{local_installer}\n\n"
                "Voulez‑vous lancer cet installateur maintenant ?",
            )
            if answer:
                try:
                    subprocess.Popen([str(local_installer)], cwd=str(local_installer.parent))
                except Exception as e:
                    messagebox.showerror(
                        "C2A",
                        f"Impossible de lancer l’installateur C2A local : {e}",
                    )
            return

        # Fallback : page officielle
        url = "https://c2a.astrosurf.com/english/support.htm#VersionHistory"
        try:
            webbrowser.open(url)
        except Exception as e:
            messagebox.showerror(
                "C2A",
                f"Impossible d'ouvrir la page C2A dans le navigateur.\nURL : {url}\n\nErreur : {e}",
            )

    # ------------------------------------------------------------------
    # Gestion persistance : TCP/IP + chemin exécutable C2A
    # ------------------------------------------------------------------
    def _load_tcp_settings(self) -> None:
        """Charge TCP/IP et le chemin C2A depuis ~/.npoap/planetarium_c2a_tcp.json."""
        try:
            if not self._tcp_settings_path.exists():
                # Valeurs par défaut TCP ; exécutable : garde l'auto-détection
                self.tcp_host_var.set("127.0.0.1")
                self.tcp_port_var.set("5876")
                self.tcp_command_template_var.set(C2A_TCP_DEFAULT_TEMPLATE)
                self.tcp_enabled_var.set(False)
                return
            with self._tcp_settings_path.open("r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            return

        try:
            host = data.get("host", "127.0.0.1")
            port = data.get("port", "5876")
            enabled = bool(data.get("enabled", False))
            template = data.get("command_template", C2A_TCP_DEFAULT_TEMPLATE)
            tnorm = str(template).strip().replace("\r", "").replace("\n", "")
            upgraded = tnorm in _C2A_TCP_LEGACY_TEMPLATES
            if upgraded:
                template = C2A_TCP_DEFAULT_TEMPLATE
            self.tcp_host_var.set(str(host))
            self.tcp_port_var.set(str(port))
            self.tcp_enabled_var.set(enabled)
            self.tcp_command_template_var.set(str(template))
            if upgraded:
                self.after_idle(self._save_tcp_settings_safe)

            exe_path = (data.get("c2a_exe_path") or "").strip()
            if exe_path and os.path.isfile(exe_path):
                self.c2a_path_var.set(exe_path)
        except Exception:
            pass

    def _save_tcp_settings_safe(self) -> None:
        """Enveloppe de sauvegarde silencieuse des paramètres TCP/IP."""
        try:
            self._save_tcp_settings()
        except Exception:
            # On ne bloque jamais l'interface pour une erreur de sauvegarde
            pass

    def _save_tcp_settings(self) -> None:
        """Sauvegarde TCP/IP et le chemin vers l'exécutable C2A (C2A.exe, C2aw.exe, …)."""
        data = {
            "host": (self.tcp_host_var.get() or "").strip() or "127.0.0.1",
            "port": (self.tcp_port_var.get() or "").strip() or "5876",
            "enabled": bool(self.tcp_enabled_var.get()),
            "command_template": self.tcp_command_template_var.get()
            or C2A_TCP_DEFAULT_TEMPLATE,
            "c2a_exe_path": (self.c2a_path_var.get() or "").strip(),
        }
        try:
            if not self._npoap_dir.exists():
                self._npoap_dir.mkdir(exist_ok=True)
        except Exception:
            pass
        try:
            with self._tcp_settings_path.open("w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        except Exception:
            pass

    def _launch_c2a(self):
        """
        Lance C2A via subprocess.
        """
        path = self.c2a_path_var.get().strip()
        if not path:
            messagebox.showerror(
                "C2A",
                "Merci d'indiquer le chemin vers l'exécutable C2A (C2A.exe).",
            )
            return

        exe_path = Path(path)
        if not exe_path.is_file():
            messagebox.showerror(
                "C2A",
                f"Le fichier indiqué n'existe pas :\n{exe_path}",
            )
            return

        # Si l'utilisateur a sélectionné par erreur un installateur dans external_apps,
        # on prévient au lieu de lancer l'installation.
        try:
            lower_parts = [p.lower() for p in exe_path.parts]
        except Exception:
            lower_parts = []
        is_external_apps = any("external_apps" in p for p in lower_parts)
        is_installer_name = any(
            kw in exe_path.name.lower()
            for kw in ("setup", "install", "installe", "installer")
        )
        if is_external_apps and is_installer_name:
            messagebox.showwarning(
                "C2A",
                "Le chemin sélectionné semble être un installateur C2A (dans external_apps).\n\n"
                "Merci de sélectionner l'exécutable installé, par exemple :\n"
                "  C:\\Program Files\\C2A\\C2A.exe\n"
                "via le bouton « Parcourir… ».",
            )
            return

        # Si un processus C2A existe déjà et est encore vivant, ne pas relancer.
        if self.c2a_process is not None and self.c2a_process.poll() is None:
            messagebox.showinfo("C2A", "C2A semble déjà lancé.")
            return

        try:
            self.c2a_process = subprocess.Popen([str(exe_path)], cwd=str(exe_path.parent))
        except Exception as e:
            self.c2a_process = None
            messagebox.showerror(
                "C2A",
                f"Impossible de lancer C2A : {e}",
            )

    def _quit_c2a(self):
        """
        Tente de fermer C2A si celui‑ci a été lancé depuis cet onglet.
        """
        if self.c2a_process is None or self.c2a_process.poll() is not None:
            messagebox.showinfo("C2A", "Aucun processus C2A géré par NPOAP n'est en cours d'exécution.")
            return

        try:
            self.c2a_process.terminate()
        except Exception as e:
            messagebox.showerror("C2A", f"Impossible de fermer C2A proprement : {e}")
        finally:
            self.c2a_process = None

    # ------------------------------------------------------------------
    # Intégration avec l'onglet « Observation de la Nuit »
    # ------------------------------------------------------------------
    def update_targets_from_night_tab(self, objects) -> None:
        """
        Reçoit une liste d'objets (EphemerisObject ou équivalent) depuis
        l'onglet « Observation de la Nuit » et les affiche dans la liste locale.
        Si le contrôle TCP est activé, tente tout de suite de pointer C2A sur la première cible.
        """
        self.targets = list(objects or [])

        self._tcp_select_suppress = True
        first_iid = None
        first_obj = None
        try:
            # Rafraîchir la Treeview
            self._targets_by_iid.clear()
            for item in self.targets_tree.get_children():
                self.targets_tree.delete(item)

            for obj in self.targets:
                name = getattr(obj, "name", "Cible")
                obj_type = getattr(obj, "obj_type", "")
                if hasattr(obj, "ra_sexagesimal") and hasattr(obj, "dec_sexagesimal"):
                    ra_str = obj.ra_sexagesimal()
                    dec_str = obj.dec_sexagesimal()
                else:
                    ra_deg = getattr(obj, "ra", None)
                    dec_deg = getattr(obj, "dec", None)
                    ra_str = f"{ra_deg:.5f}°" if ra_deg is not None else "N/A"
                    dec_str = f"{dec_deg:.5f}°" if dec_deg is not None else "N/A"
                mag = getattr(obj, "magnitude", None)
                mag_str = f"{mag:.2f}" if isinstance(mag, (int, float)) else ""

                iid = self.targets_tree.insert(
                    "",
                    tk.END,
                    values=(name, obj_type, ra_str, dec_str, mag_str),
                )
                self._targets_by_iid[iid] = obj
                if first_iid is None:
                    first_iid = iid
                    first_obj = obj

            if first_iid is not None:
                self.targets_tree.selection_set(first_iid)
                self.targets_tree.focus(first_iid)
                self.targets_tree.see(first_iid)

            if first_obj is not None:
                self._refresh_observatory_info_label()
                self._point_c2a_to_object(first_obj, feedback="errors")
        finally:
            self._tcp_select_suppress = False

    @staticmethod
    def _finalize_c2a_tcp_payload(cmd: str) -> str:
        """
        Protocole C2A : commandes ASCII, séparées par « ; », chaque commande se termine par « ; ».
        Pas de retour ligne requis côté client (doc officielle).
        """
        s = (cmd or "").strip().replace("\r", "").replace("\n", "")
        if not s:
            return ""
        if not s.endswith(";"):
            s += ";"
        return s

    def _send_c2a_tcp(self, host: str, port: int, payload: str) -> Tuple[bool, Optional[str]]:
        """
        Envoie une séquence de commandes et lit la réponse C2A (… puis OK\\r\\n).
        Retourne (True, None) si la réponse contient OK sans Error en tête de réponse.
        """
        data = self._finalize_c2a_tcp_payload(payload)
        if not data:
            return False, "Séquence de commandes vide."
        raw = data.encode("ascii", errors="replace")
        # Plusieurs stacks Windows lisent une « ligne » : CRLF après la séquence aide parfois.
        raw_line = raw + b"\r\n"
        try:
            with socket.create_connection((host, port), timeout=5.0) as sock:
                sock.sendall(raw_line)
                sock.settimeout(3.0)
                chunks: List[bytes] = []
                for _ in range(8):
                    try:
                        part = sock.recv(8192)
                    except socket.timeout:
                        break
                    if not part:
                        break
                    chunks.append(part)
                    if b"OK" in part or b"Error" in part:
                        break
                reply = b"".join(chunks).decode("ascii", errors="replace").strip()
                if not reply:
                    return (
                        False,
                        "Aucune réponse de C2A. Vérifiez qu’il s’agit bien de C2A sur le port 5876.",
                    )
                if reply.startswith("Error") or "Error;" in reply:
                    return False, reply
                if "OK" not in reply:
                    return False, f"Réponse inattendue de C2A : {reply!r}"
        except OSError as e:
            return False, str(e)
        return True, None

    def _point_c2a_to_object(
        self,
        obj,
        feedback: Literal["none", "errors", "all"] = "none",
    ) -> bool:
        """
        Construit la commande TCP (SetRa / SetDe / …) et l'envoie à C2A.
        feedback : none = silencieux ; errors = boîtes seulement si problème ;
        all = message aussi en cas de succès (bouton manuel).
        """
        if not self.tcp_enabled_var.get():
            if feedback != "none":
                messagebox.showinfo(
                    "Pointer C2A",
                    "Activez la case « Activer le contrôle TCP vers C2A » dans cet onglet.\n\n"
                    "C2A doit être ouvert. Le port d’écoute par défaut est 5876 (voir l’aide C2A, "
                    "section TCP socket).",
                )
            return False

        ra_deg = getattr(obj, "ra", None)
        dec_deg = getattr(obj, "dec", None)
        if not isinstance(ra_deg, (int, float)) or not isinstance(dec_deg, (int, float)):
            if feedback != "none":
                messagebox.showwarning(
                    "Pointer C2A",
                    "Cette cible n’a pas de coordonnées RA/DEC numériques exploitables.",
                )
            return False

        host = (self.tcp_host_var.get() or "127.0.0.1").strip()
        port_str = (self.tcp_port_var.get() or "5876").strip()
        try:
            port = int(port_str)
        except ValueError:
            if feedback != "none":
                messagebox.showwarning("Pointer C2A", "Le numéro de port TCP est invalide.")
            return False

        name = getattr(obj, "name", "")
        tpl = (self.tcp_command_template_var.get() or "").strip()
        if not tpl:
            if feedback != "none":
                messagebox.showwarning("Pointer C2A", "Le modèle de commande TCP est vide.")
            return False

        try:
            ra_hours = ra_deg / 15.0
        except Exception:
            ra_hours = 0.0

        self._refresh_observatory_info_label()
        obs_eff = _effective_observatory()
        lat_obs = _obs_latitude(obs_eff)
        lon_obs = _obs_longitude(obs_eff)

        try:
            cmd = tpl.format(
                name=name,
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                ra_hours=ra_hours,
                lat_obs=lat_obs,
                lon_obs=lon_obs,
                obs_c2a_local=observatory_local_now_c2a_string(),
            )
        except Exception as e:
            if feedback != "none":
                messagebox.showwarning(
                    "Pointer C2A",
                    f"Erreur dans le modèle de commande (accolades / champs) :\n{e}",
                )
            return False

        ok, err = self._send_c2a_tcp(host, port, cmd)
        if not ok:
            if feedback != "none":
                messagebox.showwarning(
                    "Pointer C2A",
                    f"Connexion vers {host}:{port} impossible ou erreur C2A.\n\n"
                    f"Vérifiez que C2A est lancé et que le pare-feu autorise le port 5876 en local.\n\n"
                    f"Détail : {err}",
                )
            return False

        if feedback == "all":
            messagebox.showinfo(
                "Pointer C2A",
                f"Commande envoyée à C2A pour « {name} » (RA={ra_hours:.6f} h, Dec={dec_deg:.6f}°).",
            )
        return True

    def _on_point_c2a_button(self) -> None:
        """Bouton : renvoie la commande pour la ligne sélectionnée, avec retour utilisateur."""
        selection = self.targets_tree.selection()
        if not selection:
            messagebox.showinfo(
                "Pointer C2A",
                "Sélectionnez d’abord une cible dans la liste ci-dessus.",
            )
            return
        iid = selection[0]
        obj = self._targets_by_iid.get(iid)
        if obj is None:
            return
        self._point_c2a_to_object(obj, feedback="all")

    def _on_target_selected(self, event=None) -> None:
        """
        Handler de sélection dans la liste des cibles.

        Si le contrôle TCP/IP est activé, envoie une commande formatée à C2A.
        """
        if self._tcp_select_suppress:
            return

        selection = self.targets_tree.selection()
        if not selection:
            return

        iid = selection[0]
        obj = self._targets_by_iid.get(iid)
        if obj is None:
            return

        self._point_c2a_to_object(obj, feedback="none")


