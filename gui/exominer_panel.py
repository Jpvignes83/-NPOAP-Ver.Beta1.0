# gui/exominer_panel.py
"""
Panneau « ExoMiner — validation de TCE » (appel subprocess au pipeline NASA).

Voir core.ExoMiner_analysis pour le pont subprocess.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _exominer_clone_path(install_base: Path) -> Path:
    """Chemin du clone : ``<racine NPOAP>/external_apps/Exominer``.

    Si la racine choisie est déjà le dossier ``external_apps``, on n'ajoute pas un second segment.
    """
    base = install_base.expanduser().resolve()
    if base.name.lower() == "external_apps":
        return base / "Exominer"
    return base / "external_apps" / "Exominer"


def _parse_tic_id_entry(raw: str) -> int:
    """Extrait l'entier TIC depuis « 147660201 », « TIC 147660201 », etc."""
    s = raw.strip()
    if not s:
        raise ValueError("Indiquez un ID TIC.")
    m = re.search(r"(?:TIC\s*)?(\d{6,})", s, re.I)
    if not m:
        raise ValueError("ID TIC invalide (ex. 147660201 ou TIC 147660201).")
    return int(m.group(1))


def _npoap_root_from_exominer_clone(exominer_dir: Path) -> Path:
    """Déduit la racine NPOAP depuis le dossier clone ``…/Exominer`` (corrige un ``external_apps`` en trop)."""
    ex = exominer_dir.expanduser().resolve()
    if ex.name.lower() != "exominer":
        return ex.parent.parent
    ep = ex.parent
    if ep.name.lower() != "external_apps":
        return ep.parent
    gg = ep.parent
    if gg.name.lower() == "external_apps":
        return gg.parent
    return gg


def _discover_bundled_exominer() -> Optional[Path]:
    """Clone présent avec les sources NPOAP : ``…/external_apps/Exominer`` (sans variable d'environnement)."""
    candidate = Path(__file__).resolve().parent.parent / "external_apps" / "Exominer"
    script = candidate / "exominer_pipeline" / "run_pipeline.py"
    if script.is_file():
        return candidate
    return None


class ExominerPanel(ttk.Frame):
    """Installation sous ``<racine>/external_apps/Exominer``, lancement via ``run_pipeline.py``."""

    def __init__(self, parent):
        super().__init__(parent, padding=8)
        self._run_thread: Optional[threading.Thread] = None
        self._clone_thread: Optional[threading.Thread] = None
        self._npoap_install_root: Optional[Path] = None
        _inst = os.environ.get("NPOAP_INSTALL_ROOT", "").strip()
        if _inst:
            self._npoap_install_root = Path(_inst).expanduser().resolve()
        else:
            _bundled = _discover_bundled_exominer()
            if _bundled is not None:
                self._npoap_install_root = _npoap_root_from_exominer_clone(_bundled)

        row_install = ttk.Frame(self)
        row_install.pack(fill=tk.X, pady=(0, 4))
        ttk.Button(
            row_install,
            text="Installer / MAJ ExoMiner (git)",
            command=self._install_exominer_clicked,
        ).pack(side=tk.LEFT)
        ttk.Button(
            row_install,
            text="Copier poids ExoMiner (Docker/Podman)",
            command=self._copy_exominer_data_clicked,
        ).pack(side=tk.LEFT, padx=(8, 0))

        inp = ttk.LabelFrame(self, text="Entrées pipeline", padding=8)
        inp.pack(fill=tk.X, pady=4)

        row_tic = ttk.Frame(inp)
        row_tic.pack(fill=tk.X)
        ttk.Label(row_tic, text="ID TIC:").pack(side=tk.LEFT)
        self.var_tic_id = tk.StringVar(value="")
        self.entry_tic_id = ttk.Entry(row_tic, textvariable=self.var_tic_id, width=44)
        self.entry_tic_id.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)

        row_src = ttk.Frame(inp)
        row_src.pack(fill=tk.X, pady=6)
        self.var_use_local = tk.BooleanVar(value=False)
        ttk.Radiobutton(
            row_src,
            text="Télécharger depuis MAST (réseau)",
            variable=self.var_use_local,
            value=False,
        ).pack(side=tk.LEFT)
        ttk.Radiobutton(
            row_src,
            text="Dossier local SPOC (FITS + DV XML)",
            variable=self.var_use_local,
            value=True,
        ).pack(side=tk.LEFT, padx=(12, 0))

        self.row_csv_local = ttk.Frame(inp)
        ttk.Label(self.row_csv_local, text="CSV tic_id / sector_run:").pack(side=tk.LEFT)
        self.var_tics_csv = tk.StringVar(value="")
        ttk.Entry(self.row_csv_local, textvariable=self.var_tics_csv, width=40).pack(
            side=tk.LEFT, padx=4, fill=tk.X, expand=True
        )
        ttk.Button(self.row_csv_local, text="Parcourir…", command=self._browse_tics).pack(side=tk.LEFT)
        ttk.Button(self.row_csv_local, text="Modèle vide…", command=self._save_tics_template).pack(
            side=tk.LEFT, padx=(4, 0)
        )

        self.row_ext = ttk.Frame(inp)
        self.row_ext.pack(fill=tk.X)
        ttk.Label(self.row_ext, text="Dossier données locales:").pack(side=tk.LEFT)
        self.var_external = tk.StringVar(value="")
        self.entry_external = ttk.Entry(self.row_ext, textvariable=self.var_external, width=48)
        self.entry_external.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        ttk.Button(self.row_ext, text="Parcourir…", command=self._browse_external).pack(side=tk.LEFT)

        opt = ttk.LabelFrame(self, text="Options & sortie", padding=8)
        opt.pack(fill=tk.X, pady=4)

        row_o = ttk.Frame(opt)
        row_o.pack(fill=tk.X)
        ttk.Label(row_o, text="Répertoire sortie run:").pack(side=tk.LEFT)
        self.var_output = tk.StringVar(value=str(Path.home() / ".npoap" / "Exominer"))
        ttk.Entry(row_o, textvariable=self.var_output, width=46).pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        ttk.Button(row_o, text="Parcourir…", command=self._browse_output).pack(side=tk.LEFT)

        row_m = ttk.Frame(opt)
        row_m.pack(fill=tk.X, pady=4)
        ttk.Label(row_m, text="Mode données:").pack(side=tk.LEFT)
        self.var_mode = tk.StringVar(value="2min")
        ttk.Combobox(
            row_m,
            textvariable=self.var_mode,
            values=("2min", "FFI"),
            width=8,
            state="readonly",
        ).pack(side=tk.LEFT, padx=4)
        ttk.Label(row_m, text="Modèle:").pack(side=tk.LEFT, padx=(16, 0))
        self.var_model = tk.StringVar(value="exominer++_single")
        ttk.Combobox(
            row_m,
            textvariable=self.var_model,
            values=(
                "exominer++_single",
                "exominer++_cviter-mean-ensemble",
                "exominer++_cv-super-mean-ensemble",
            ),
            width=28,
            state="readonly",
        ).pack(side=tk.LEFT, padx=4)

        row_p = ttk.Frame(opt)
        row_p.pack(fill=tk.X)
        ttk.Label(row_p, text="Processes:").pack(side=tk.LEFT)
        self.var_nproc = tk.StringVar(value="1")
        ttk.Entry(row_p, textvariable=self.var_nproc, width=6).pack(side=tk.LEFT, padx=4)
        ttk.Label(row_p, text="Jobs:").pack(side=tk.LEFT)
        self.var_njobs = tk.StringVar(value="1")
        ttk.Entry(row_p, textvariable=self.var_njobs, width=6).pack(side=tk.LEFT, padx=4)
        self.var_dl_dv_urls = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            row_p,
            text="CSV URLs rapports DV MAST",
            variable=self.var_dl_dv_urls,
        ).pack(side=tk.LEFT, padx=(16, 0))

        row_cat = ttk.Frame(opt)
        row_cat.pack(fill=tk.X, pady=4)
        ttk.Label(row_cat, text="Paramètres stellaires:").pack(side=tk.LEFT)
        self.var_stellar = tk.StringVar(value="ticv8")
        ttk.Entry(row_cat, textvariable=self.var_stellar, width=18).pack(side=tk.LEFT, padx=4)
        ttk.Label(row_cat, text="RUWE:").pack(side=tk.LEFT, padx=(8, 0))
        self.var_ruwe = tk.StringVar(value="gaiadr2")
        ttk.Entry(row_cat, textvariable=self.var_ruwe, width=14).pack(side=tk.LEFT, padx=4)

        btn_row = ttk.Frame(self)
        btn_row.pack(fill=tk.X, pady=6)
        ttk.Button(btn_row, text="Lancer Exominer", command=self._run_clicked).pack(side=tk.LEFT)

        log_fr = ttk.LabelFrame(self, text="Journal", padding=4)
        log_fr.pack(fill=tk.BOTH, expand=True, pady=4)
        self.log_text = tk.Text(log_fr, wrap=tk.WORD, height=12, font=("Consolas", 8))
        sb = ttk.Scrollbar(log_fr, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=sb.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        self.var_use_local.trace_add("write", lambda *_: self._toggle_external())

        self._toggle_external()

    def _toggle_external(self):
        local = self.var_use_local.get()
        self.entry_external.configure(state=tk.NORMAL if local else tk.DISABLED)
        self.entry_tic_id.configure(state=tk.DISABLED if local else tk.NORMAL)
        if local:
            self.row_csv_local.pack(fill=tk.X, pady=(0, 4), before=self.row_ext)
        else:
            self.row_csv_local.pack_forget()

    def _log(self, s: str):
        self.log_text.insert(tk.END, s + "\n")
        self.log_text.see(tk.END)

    def _resolve_install_base(self) -> Optional[Path]:
        env_ib = os.environ.get("NPOAP_INSTALL_ROOT", "").strip()
        if env_ib:
            return Path(env_ib).expanduser().resolve()
        if self._npoap_install_root is not None:
            return self._npoap_install_root
        er = os.environ.get("NPOAP_EXOMINER_ROOT", "").strip()
        if er:
            p = Path(er).expanduser().resolve()
            script = p / "exominer_pipeline" / "run_pipeline.py"
            if p.is_dir() and script.is_file():
                return _npoap_root_from_exominer_clone(p)
        return None

    def _exominer_repo_root(self) -> Optional[Path]:
        """Racine du clone (dossier contenant ``exominer_pipeline/run_pipeline.py``)."""
        base = self._resolve_install_base()
        if base is not None:
            candidate = _exominer_clone_path(base)
            script = candidate / "exominer_pipeline" / "run_pipeline.py"
            if script.is_file():
                return candidate.resolve()
        er = os.environ.get("NPOAP_EXOMINER_ROOT", "").strip()
        if er:
            p = Path(er).expanduser().resolve()
            script = p / "exominer_pipeline" / "run_pipeline.py"
            if script.is_file():
                return p.resolve()
        bundled = _discover_bundled_exominer()
        return bundled.resolve() if bundled is not None else None

    def _sync_exominer_env(self, repo: Path) -> None:
        """Persistance de session : même après redémarrage de NPOAP sans réinstaller."""
        repo = repo.resolve()
        self._npoap_install_root = _npoap_root_from_exominer_clone(repo)
        os.environ["NPOAP_INSTALL_ROOT"] = str(self._npoap_install_root)
        os.environ["NPOAP_EXOMINER_ROOT"] = str(repo)

    def _install_exominer_clicked(self) -> None:
        if self._clone_thread is not None and self._clone_thread.is_alive():
            messagebox.showwarning("ExoMiner", "Une installation ou mise à jour est déjà en cours.")
            return
        try:
            from shutil import which

            from core.ExoMiner_analysis import EXOMINER_REPO_URL

            if not which("git"):
                messagebox.showerror(
                    "ExoMiner",
                    "Git est introuvable dans le PATH.\nInstallez Git for Windows puis rouvrez NPOAP.",
                )
                return
            install_base = filedialog.askdirectory(title="")
            if not install_base:
                return
            base_p = Path(install_base).expanduser().resolve()
            dest = _exominer_clone_path(base_p)
            parent_p = dest.parent  # …/external_apps
            try:
                parent_p.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                messagebox.showerror("ExoMiner", f"Impossible de créer le dossier :\n{parent_p}\n\n{e}")
                return
            action = "clone"
            if dest.exists():
                if (dest / ".git").is_dir():
                    if not messagebox.askyesno(
                        "ExoMiner",
                        f"Le dossier existe déjà :\n{dest}\n\nMettre à jour avec git pull ?",
                    ):
                        return
                    action = "pull"
                else:
                    messagebox.showerror(
                        "ExoMiner",
                        f"Un dossier « Exominer » existe déjà sans être un dépôt Git :\n{dest}",
                    )
                    return

            repo_git = f"{EXOMINER_REPO_URL}.git"
            self._log(f"=== ExoMiner ({action}) → {dest} ===")
            self._clone_thread = threading.Thread(
                target=self._worker_clone_install,
                args=(parent_p, dest, action, repo_git),
                daemon=True,
            )
            self._clone_thread.start()
        except Exception as e:
            logger.exception("install ExoMiner")
            messagebox.showerror("ExoMiner", str(e))

    def _copy_exominer_data_clicked(self) -> None:
        """Lance le script PowerShell qui extrait ``exominer_pipeline/data`` depuis ghcr.io/nasa/exominer."""
        npoap = Path(__file__).resolve().parent.parent
        ps1 = npoap / "scripts" / "exominer_copy_pipeline_data.ps1"
        if not ps1.is_file():
            messagebox.showerror("ExoMiner", f"Script introuvable :\n{ps1}")
            return
        if sys.platform != "win32":
            messagebox.showinfo(
                "ExoMiner",
                "Sous Linux/macOS, exécutez dans un terminal :\n"
                f"  pwsh -File {ps1}\n\n"
                "Ou copiez manuellement depuis Podman/Docker (voir la doc NASA).",
            )
            return
        dest_data = npoap / "external_apps" / "Exominer" / "exominer_pipeline" / "data"
        if not (dest_data.parent).is_dir():
            messagebox.showwarning(
                "ExoMiner",
                "Le dossier exominer_pipeline du clone est absent.\n"
                "Utilisez d’abord « Installer / MAJ ExoMiner (git) ».",
            )
            return
        try:
            win = os.environ.get("WINDIR") or os.environ.get("SystemRoot") or r"C:\Windows"
            ps_exe = Path(win) / "System32" / "WindowsPowerShell" / "v1.0" / "powershell.exe"
            if not ps_exe.is_file():
                p = shutil.which("powershell.exe")
                ps_exe = Path(p) if p else Path("powershell.exe")
            # cmd /k : la console reste ouverte même si PowerShell quitte en erreur avant Read-Host
            subprocess.Popen(
                [
                    str(Path(win) / "System32" / "cmd.exe"),
                    "/k",
                    str(ps_exe),
                    "-NoLogo",
                    "-NoProfile",
                    "-ExecutionPolicy",
                    "Bypass",
                    "-WindowStyle",
                    "Normal",
                    "-File",
                    str(ps1),
                ],
                cwd=str(npoap),
                creationflags=subprocess.CREATE_NEW_CONSOLE,  # type: ignore[attr-defined]
            )
            self._log(f"Copie des données ExoMiner (nouvelle console PowerShell) : {ps1}")
            messagebox.showinfo(
                "ExoMiner",
                "Une fenêtre PowerShell va télécharger l’image ghcr.io/nasa/exominer et copier "
                "exominer_pipeline/data vers votre clone.\n\n"
                "Prérequis : Docker Desktop ou Podman installé.\n"
                "À la fin, refermez la fenêtre et relancez le pipeline dans NPOAP.",
            )
        except OSError as e:
            logger.exception("copy ExoMiner data")
            messagebox.showerror("ExoMiner", str(e))

    def _subprocess_kw_hide_console(self) -> dict:
        if sys.platform == "win32":
            return {"creationflags": subprocess.CREATE_NO_WINDOW}
        return {}

    def _worker_clone_install(self, parent_p: Path, dest: Path, action: str, repo_git: str) -> None:
        kw = dict(capture_output=True, text=True, timeout=3600, encoding="utf-8", errors="replace")
        kw.update(self._subprocess_kw_hide_console())
        proc: subprocess.CompletedProcess[str]
        try:
            if action == "clone":
                proc = subprocess.run(
                    ["git", "clone", "--depth", "1", repo_git, str(dest)],
                    cwd=str(parent_p),
                    **kw,
                )
            else:
                proc = subprocess.run(["git", "pull"], cwd=str(dest), **kw)
        except subprocess.TimeoutExpired:
            err = "Dépassement du délai (réseau lent ou dépôt volumineux)."
            logger.error(err)

            def tout():
                self._log(err)
                messagebox.showerror("ExoMiner", err)

            self.after(0, tout)
            return
        except Exception as e:
            logger.exception("clone/pull ExoMiner")
            err_msg = str(e)

            def ex():
                self._log(err_msg)
                messagebox.showerror("ExoMiner", err_msg)

            self.after(0, ex)
            return

        out = (proc.stdout or "") + (proc.stderr or "")

        def done():
            self._log(out[-16000:] if len(out) > 16000 else out)
            if proc.returncode != 0:
                messagebox.showerror(
                    "ExoMiner",
                    f"git a échoué (code {proc.returncode}).\nVoir le journal.",
                )
                return
            script = dest / "exominer_pipeline" / "run_pipeline.py"
            if not script.is_file():
                messagebox.showwarning(
                    "ExoMiner",
                    f"Clone terminé mais script introuvable :\n{script}\n"
                    "Vérifiez la racine dépôt ou réessayez sans --depth si besoin.",
                )
            install_base = _npoap_root_from_exominer_clone(dest)
            self._npoap_install_root = install_base
            os.environ["NPOAP_INSTALL_ROOT"] = str(install_base)
            root_str = str(dest.resolve())
            os.environ["NPOAP_EXOMINER_ROOT"] = root_str
            messagebox.showinfo(
                "ExoMiner",
                f"{'Mise à jour' if action == 'pull' else 'Installation'} terminée.\n\nRacine :\n{root_str}",
            )

        self.after(0, done)

    def _browse_tics(self):
        f = filedialog.askopenfilename(
            title="Table TIC (CSV)",
            filetypes=[("CSV", "*.csv"), ("Tous", "*.*")],
        )
        if f:
            self.var_tics_csv.set(f)

    def _browse_external(self):
        d = filedialog.askdirectory(title="Dossier MAST / SPOC (FITS + DV XML)")
        if d:
            self.var_external.set(d)

    def _browse_output(self):
        d = filedialog.askdirectory(title="Répertoire de sortie du run")
        if d:
            self.var_output.set(d)

    def _save_tics_template(self):
        fp = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
            title="Enregistrer modèle tic_id, sector_run",
        )
        if not fp:
            return
        Path(fp).write_text("tic_id,sector_run\n167526485,6-6\n", encoding="utf-8")
        self.var_tics_csv.set(fp)
        messagebox.showinfo("Modèle", f"Fichier créé :\n{fp}")

    def _validate_tics_csv(self, fp: str) -> bool:
        try:
            import pandas as pd

            df = pd.read_csv(fp)
            for col in ("tic_id", "sector_run"):
                if col not in df.columns:
                    messagebox.showerror(
                        "CSV invalide",
                        f"Colonne obligatoire manquante : {col}\nColonnes trouvées : {list(df.columns)}",
                    )
                    return False
            return True
        except Exception as e:
            messagebox.showerror("CSV", str(e))
            return False

    def _run_clicked(self):
        if self._run_thread is not None and self._run_thread.is_alive():
            messagebox.showwarning("ExoMiner", "Un calcul est déjà en cours.")
            return

        root_path = self._exominer_repo_root()
        script = (root_path / "exominer_pipeline" / "run_pipeline.py") if root_path else None
        if root_path is None or not root_path.is_dir() or script is None or not script.is_file():
            hint = Path(__file__).resolve().parent.parent / "external_apps" / "Exominer"
            messagebox.showerror(
                "ExoMiner",
                "Installation introuvable ou pipeline incomplet.\n\n"
                "• Utilisez « Installer / MAJ ExoMiner (git) » et choisissez la racine NPOAP, ou\n"
                f"• Vérifiez que le clone existe (ex. :\n{hint})",
            )
            return

        self._sync_exominer_env(root_path)
        self._log(f"ExoMiner détecté : {root_path.resolve()}")

        local = self.var_use_local.get()
        if local:
            tics = self.var_tics_csv.get().strip()
            if not tics or not Path(tics).is_file():
                messagebox.showerror("ExoMiner", "Choisissez un fichier CSV tic_id / sector_run (mode local).")
                return
            if not self._validate_tics_csv(tics):
                return
            ext_dir = self.var_external.get().strip()
            if not ext_dir or not Path(ext_dir).is_dir():
                messagebox.showerror("ExoMiner", "Dossier données locales invalide.")
                return
        else:
            try:
                _parse_tic_id_entry(self.var_tic_id.get())
            except ValueError as e:
                messagebox.showerror("ExoMiner", str(e))
                return

        try:
            max(1, int(self.var_nproc.get()))
            max(1, int(self.var_njobs.get()))
        except ValueError:
            messagebox.showerror("ExoMiner", "Processes et Jobs doivent être des entiers.")
            return

        self._run_thread = threading.Thread(target=self._worker_run, daemon=True)
        self._run_thread.start()

    def _worker_run(self):
        try:
            from core.ExoMiner_analysis import (
                check_exominer_dependencies,
                check_exominer_pipeline_data,
                mast_sector_run_for_tic,
                resolve_exominer_python,
                run_exominer_pipeline,
                write_tics_table_csv,
            )

            root_path = self._exominer_repo_root()
            if root_path is None:
                raise RuntimeError("Racine ExoMiner introuvable.")
            script_chk = root_path / "exominer_pipeline" / "run_pipeline.py"
            if not script_chk.is_file():
                raise RuntimeError(f"Pipeline introuvable : {script_chk}")
            self._sync_exominer_env(root_path)
            root = str(root_path.resolve())
            output_dir = Path(self.var_output.get().strip()).expanduser().resolve()
            output_dir.mkdir(parents=True, exist_ok=True)

            self.after(0, lambda: self._log("Vérification dépendances ExoMiner (TensorFlow, lightkurve, pydl)…"))
            dep_err = check_exominer_dependencies(resolve_exominer_python(None))
            if dep_err:
                raise RuntimeError(dep_err)
            self.after(0, lambda: self._log("Dépendances principales OK pour l'interpréteur ExoMiner."))

            self.after(0, lambda: self._log("Vérification poids / norm_stats ExoMiner (.keras, norm_stats)…"))
            pdata_err = check_exominer_pipeline_data(
                root_path,
                exominer_model=self.var_model.get().strip(),
            )
            if pdata_err:
                raise RuntimeError(pdata_err)
            self.after(0, lambda: self._log("Fichiers modèle / normalisation présents (ou chemin .keras valide)."))

            local = self.var_use_local.get()
            mode = self.var_mode.get().strip().lower()

            if not local:
                tic_int = _parse_tic_id_entry(self.var_tic_id.get())
                self.after(0, lambda t=tic_int: self._log(f"Requête MAST pour TIC {t} (sector_run)…"))
                sector_run = mast_sector_run_for_tic(tic_int, data_collection_mode=mode)
                self.after(0, lambda sr=sector_run: self._log(f"sector_run déduit : {sr}"))
                tics_fp = output_dir / "npoap_tics_input.csv"
                write_tics_table_csv(tics_fp, [{"tic_id": tic_int, "sector_run": sector_run}])
                tic_ids_fp = str(tics_fp)
            else:
                tic_ids_fp = self.var_tics_csv.get().strip()

            nproc = max(1, int(self.var_nproc.get()))
            njobs = max(1, int(self.var_njobs.get()))
            ext_dir = self.var_external.get().strip()

            self.after(0, lambda: self._log("=== Début pipeline ExoMiner ==="))
            res = run_exominer_pipeline(
                output_dir=output_dir,
                tic_ids_fp=tic_ids_fp,
                tic_ids=None,
                data_collection_mode=mode,
                num_processes=nproc,
                num_jobs=njobs,
                download_spoc_data_products="true" if self.var_dl_dv_urls.get() else "false",
                external_data_repository=None if not local else ext_dir,
                stellar_parameters_source=self.var_stellar.get().strip(),
                ruwe_source=self.var_ruwe.get().strip(),
                exominer_model=self.var_model.get().strip(),
                exominer_root=root,
                python_executable=None,
                timeout_sec=None,
                skip_tf_check=True,
                skip_pipeline_data_check=True,
            )

            def done():
                self._log((res.stdout or "")[-12000:])
                if res.stderr:
                    self._log("--- stderr ---")
                    self._log((res.stderr or "")[-8000:])
                if res.success and res.predictions_csv:
                    self._log(f"Prédictions : {res.predictions_csv}")
                    messagebox.showinfo("ExoMiner", f"Terminé.\n{res.predictions_csv}")
                else:
                    messagebox.showerror("ExoMiner", f"Code retour {res.returncode}\nVoir journal.")

            self.after(0, done)
        except Exception as e:
            logger.exception("run ExoMiner")
            err_msg = str(e)

            def err():
                self._log(err_msg)
                messagebox.showerror("ExoMiner", err_msg)

            self.after(0, err)
