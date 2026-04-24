# gui/manual_help.py
"""
Aide contextuelle : ouverture du manuel utilisateur (MANUEL_UTILISATEUR.md) sur le chapitre correspondant.
"""

from __future__ import annotations

import logging
from pathlib import Path
import tkinter as tk
from tkinter import messagebox, ttk
from tkinter.scrolledtext import ScrolledText

logger = logging.getLogger(__name__)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def user_manual_path() -> Path:
    return _repo_root() / "docs" / "MANUEL_UTILISATEUR.md"


# Fragments d’ancrage = cibles du sommaire du manuel (titres ## ou balise <a id="...">)
ANCHOR_LINE_MARKERS: dict[str, str] = {
    "1-vue-densemble": "## 1. Vue d'ensemble",
    "2-accueil": "## 2. Accueil",
    "3-réduction-de-données": "## 3. Réduction de Données",
    "coef-transformation-reduction": '<a id="coef-transformation-reduction"></a>',
    "4-photométrie-exoplanètes": "## 4. Photométrie Exoplanètes",
    "5-photométrie-astéroïdes": "## 5. Photométrie Astéroïdes",
    "6-photométrie-transitoires": "## 6. Photométrie Transitoires",
    "7-analyse-de-données": "## 7. Analyse de Données",
    "8-étoiles-binaires": "## 8. Étoiles Binaires",
    "9-easy-lucky-imaging": "## 9. Easy Lucky Imaging",
    "10-analyse-des-outils-et-fichiers-compatibles-lctools-tess": "## 10. Analyse des outils et fichiers compatibles LcTools (TESS)",
    "11-spectroscopie": "## 11. Spectroscopie",
    "12-observation-de-la-nuit": "## 12. Observation de la nuit",
    "13-occultations-sora-pymovie-analyse": "## 13. Occultations (SORA, PyMovie, analyse)",
    "14-conseils-généraux": "## 14. Conseils généraux",
    # Ancres supplémentaires (balises <a id="...">) dans la section « Vue d'ensemble »
    "planétarium-c2a": '<a id="planétarium-c2a"></a>',
    "onglet-catalogues": '<a id="onglet-catalogues"></a>',
    "analyse-amas-open-clusters": '<a id="analyse-amas-open-clusters"></a>',
}


def _line_index_for_anchor(lines: list[str], anchor: str) -> int:
    marker = ANCHOR_LINE_MARKERS.get(anchor)
    if not marker:
        return 0
    for i, line in enumerate(lines):
        if marker in line:
            return i
    return 0


def open_user_manual_chapter(anchor: str, parent: tk.Misc | None = None) -> None:
    """
    Affiche une fenêtre avec le texte du manuel à partir du chapitre identifié par ``anchor``
    (fragment type ``3-réduction-de-données``).
    """
    path = user_manual_path()
    if not path.is_file():
        messagebox.showerror(
            "Manuel",
            f"Fichier introuvable :\n{path}\n\nRéinstallez ou restaurez docs/MANUEL_UTILISATEUR.md.",
            parent=parent,
        )
        return
    try:
        raw = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        messagebox.showerror("Manuel", f"Lecture impossible :\n{e}", parent=parent)
        return

    lines = raw.splitlines()
    start = _line_index_for_anchor(lines, anchor)
    if anchor not in ANCHOR_LINE_MARKERS:
        logger.warning("Ancre manuel inconnue : %s — affichage depuis le début.", anchor)

    title = ANCHOR_LINE_MARKERS.get(anchor, "Manuel utilisateur")
    win = tk.Toplevel(parent)
    win.title(f"Manuel NPOAP — {title.replace('## ', '')[:72]}")
    win.geometry("920x640")
    outer = ttk.Frame(win, padding=8)
    outer.pack(fill=tk.BOTH, expand=True)
    ttk.Label(
        outer,
        text="docs/MANUEL_UTILISATEUR.md",
        font=("Arial", 8),
        foreground="gray",
    ).pack(anchor=tk.W)
    st = ScrolledText(
        outer,
        wrap=tk.WORD,
        font=("Consolas", 10) if _is_windows() else ("monospace", 10),
        height=28,
        width=100,
    )
    st.pack(fill=tk.BOTH, expand=True, pady=(4, 0))
    body = "\n".join(lines[start:])
    st.insert(tk.END, body)
    st.configure(state=tk.DISABLED)
    try:
        st.focus_set()
    except tk.TclError:
        pass


def _is_windows() -> bool:
    import sys

    return sys.platform == "win32"


def add_manual_help_header(parent: tk.Misc, anchor: str) -> ttk.Frame:
    """
    Insère une barre en tête du widget ``parent`` avec un bouton « Aide (manuel)… ».
    À appeler **avant** les autres ``pack`` / ``grid`` du contenu principal, ou adapter le placement.
    """
    bar = ttk.Frame(parent)
    bar.pack(side=tk.TOP, fill=tk.X, padx=4, pady=(2, 4))
    ttk.Label(bar, text="").pack(side=tk.LEFT, fill=tk.X, expand=True)
    ttk.Button(
        bar,
        text="Aide (manuel)…",
        width=18,
        command=lambda a=anchor: open_user_manual_chapter(
            a, parent.winfo_toplevel()
        ),
    ).pack(side=tk.RIGHT, padx=2, pady=0)
    return bar


def add_manual_help_header_grid(
    parent: tk.Misc,
    anchor: str,
    *,
    row: int = 0,
    column: int = 0,
    columnspan: int = 1,
    padx: int = 4,
    pady: tuple[int, int] = (2, 4),
) -> ttk.Frame:
    """
    Même rôle que ``add_manual_help_header`` mais placement en **grid** (ex. onglet entièrement en grid).
    """
    bar = ttk.Frame(parent)
    bar.grid(
        row=row,
        column=column,
        columnspan=columnspan,
        sticky="ew",
        padx=padx,
        pady=pady,
    )
    ttk.Label(bar, text="").pack(side=tk.LEFT, fill=tk.X, expand=True)
    ttk.Button(
        bar,
        text="Aide (manuel)…",
        width=18,
        command=lambda a=anchor: open_user_manual_chapter(
            a, parent.winfo_toplevel()
        ),
    ).pack(side=tk.RIGHT, padx=2, pady=0)
    return bar
