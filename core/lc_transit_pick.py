"""
Éditeur matplotlib des fenêtres « en transit » pour normalisation OOT (TESS → LcTools).

Raccourcis (voir aussi docstring du script tess_lc_fits_to_txt) :
  - glisser (bouton gauche) : ajouter une bande temporelle [tmin, tmax]
  - clic sur une bande : la sélectionner
  - d : supprimer la bande sélectionnée
  - + / - : élargir / rétrécir la largeur du masque de 5 % (centré sur le minimum de flux dans la bande)
  - Échap : désélectionner
  - u : annuler la dernière modification
  - q ou Entrée : valider et fermer
"""

from __future__ import annotations

from copy import deepcopy
from typing import List, Optional, Tuple

import numpy as np

import matplotlib.pyplot as plt


def _mask_to_intervals(time: np.ndarray, mask: np.ndarray) -> List[Tuple[float, float]]:
    m = np.asarray(mask, dtype=bool)
    if not np.any(m):
        return []
    idx = np.where(m)[0]
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.concatenate([[0], breaks + 1])
    ends = np.concatenate([breaks, [len(idx) - 1]])
    out: List[Tuple[float, float]] = []
    for s, e in zip(starts, ends):
        ia, ib = int(idx[s]), int(idx[e])
        out.append((float(time[ia]), float(time[ib])))
    return out


def _intervals_to_mask(time: np.ndarray, intervals: List[Tuple[float, float]]) -> np.ndarray:
    t = np.asarray(time, dtype=float)
    mask = np.zeros(len(t), dtype=bool)
    for a, b in intervals:
        lo, hi = (a, b) if a <= b else (b, a)
        mask |= (t >= lo) & (t <= hi)
    return mask


def pick_in_transit_mask(
    time_btjd: np.ndarray,
    flux_arr: np.ndarray,
    base_in_transit: Optional[np.ndarray] = None,
    title: str = "Transits (baseline OOT)",
) -> np.ndarray:
    """
    Affiche la LC ; l’utilisateur définit des intervalles « en transit ».
    Retourne un masque booléen de même longueur que *time_btjd*.
    """
    time_btjd = np.asarray(time_btjd, dtype=float).reshape(-1)
    flux_arr = np.asarray(flux_arr, dtype=float).reshape(-1)
    if time_btjd.shape != flux_arr.shape:
        raise ValueError("time_btjd et flux_arr doivent avoir la même longueur.")

    n = len(time_btjd)
    if base_in_transit is not None:
        bi = np.asarray(base_in_transit, dtype=bool).reshape(-1)
        if bi.shape[0] != n:
            raise ValueError("base_in_transit : longueur incompatible.")
        intervals: List[List[float]] = [list(x) for x in _mask_to_intervals(time_btjd, bi)]
    else:
        intervals = []

    history: List[List[List[float]]] = []
    _span_raw = float(np.nanmax(time_btjd) - np.nanmin(time_btjd))
    span = _span_raw if np.isfinite(_span_raw) and _span_raw > 0 else 1.0
    t_plot_min = float(np.nanmin(time_btjd))
    t_plot_max = float(np.nanmax(time_btjd))

    def _time_at_transit_bottom(a: float, b: float) -> float:
        """Temps du minimum de flux dans [a,b] (bas du transit) ; sinon milieu géométrique."""
        lo_i, hi_i = (a, b) if a <= b else (b, a)
        m = (time_btjd >= lo_i) & (time_btjd <= hi_i)
        if not np.any(m):
            return 0.5 * (lo_i + hi_i)
        idx = np.where(m)[0]
        flux_sel = np.asarray(flux_arr[idx], dtype=np.float64)
        if not np.any(np.isfinite(flux_sel)):
            return 0.5 * (lo_i + hi_i)
        k_rel = int(np.nanargmin(flux_sel))
        return float(time_btjd[idx[k_rel]])

    def _clamp_interval_to_plot(tc: float, half: float) -> tuple[float, float] | None:
        """[tc-half, tc+half] ramené dans [t_plot_min, t_plot_max] ; None si largeur nulle."""
        nlo, nhi = tc - half, tc + half
        if nlo < t_plot_min:
            sh = t_plot_min - nlo
            nlo += sh
            nhi += sh
        if nhi > t_plot_max:
            sh = nhi - t_plot_max
            nlo -= sh
            nhi -= sh
        nlo = max(nlo, t_plot_min)
        nhi = min(nhi, t_plot_max)
        if nhi - nlo <= span * 1e-12:
            return None
        return nlo, nhi

    fig, (ax, ax_help) = plt.subplots(
        2,
        1,
        figsize=(11, 6.4),
        gridspec_kw={"height_ratios": [4.0, 1.25], "hspace": 0.32},
    )
    try:
        mgr = fig.canvas.manager
        if mgr is not None and hasattr(mgr, "set_window_title"):
            mgr.set_window_title("NPOAP — Éditeur masques transit (LcTools)")
    except Exception:
        pass

    ax.plot(time_btjd, flux_arr, "k.", markersize=2, alpha=0.65, rasterized=True)
    ax.set_xlabel("Temps (BTJD)")
    ax.set_ylabel("Flux")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    ax_help.set_axis_off()
    _instructions = (
        "Outils d’édition — Souris : glisser horizontalement (bouton gauche) = ajouter une bande « en transit » ; "
        "clic sans glisser sur une bande = la sélectionner (bord plus foncé). "
        "Si vous utilisez le zoom ou le déplacement de la barre matplotlib, repassez sur l’outil « curseur » (flèche) "
        "avant de dessiner des bandes — sinon le glissé de zoom serait pris pour un masque.\n"
        "Touches :  d  = supprimer la bande sélectionnée  ·  + / -  = élargir / rétrécir la largeur de **5 %** "
        "(centré sur le **minimum de flux** dans la bande)  ·  Échap  = désélectionner  ·  "
        "u  = annuler la dernière modification  ·  q  ou  Entrée  = valider et fermer la fenêtre."
    )
    ax_help.text(
        0.5,
        0.95,
        _instructions,
        ha="center",
        va="top",
        fontsize=9,
        transform=ax_help.transAxes,
        linespacing=1.4,
        color="#222222",
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "#f5f5f5", "edgecolor": "#cccccc", "linewidth": 1},
    )

    span_artists: List[object] = []
    selected: List[int] = []

    def intervals_snapshot() -> List[List[float]]:
        return deepcopy(intervals)

    def push_history() -> None:
        history.append(intervals_snapshot())
        if len(history) > 40:
            history.pop(0)

    def draw_intervals() -> None:
        for art in span_artists:
            try:
                if hasattr(art, "remove"):
                    art.remove()
            except Exception:
                pass
        span_artists.clear()
        for i, (lo, hi) in enumerate(intervals):
            lo, hi = (lo, hi) if lo <= hi else (hi, lo)
            face = (1.0, 0.3, 0.3, 0.35) if selected and i == selected[0] else (0.9, 0.4, 0.4, 0.22)
            edg = "darkred" if selected and i == selected[0] else "red"
            poly = ax.axvspan(
                lo,
                hi,
                facecolor=face,
                edgecolor=edg,
                linewidth=1.0,
                zorder=2,
            )
            span_artists.append(poly)
        fig.canvas.draw_idle()

    press_x: List[Optional[float]] = [None]

    def _toolbar_or_lock_uses_mouse() -> bool:
        """True si zoom/pan matplotlib ou verrou canvas : ne pas traiter le glissé comme bande transit."""
        canvas = fig.canvas
        tb = getattr(canvas, "toolbar", None)
        if tb is not None:
            mode = getattr(tb, "mode", None)
            mode_s = ("" if mode is None else str(mode)).strip()
            if mode_s != "":
                return True
        try:
            if canvas.widgetlock.locked():
                return True
        except Exception:
            pass
        return False

    def on_press(event) -> None:
        if event.inaxes != ax or event.button != 1:
            return
        if _toolbar_or_lock_uses_mouse():
            press_x[0] = None
            return
        press_x[0] = float(event.xdata)

    def on_release(event) -> None:
        if event.inaxes != ax or event.button != 1:
            return
        if _toolbar_or_lock_uses_mouse():
            press_x[0] = None
            return
        if press_x[0] is None or event.xdata is None:
            press_x[0] = None
            return
        x0 = press_x[0]
        x1 = float(event.xdata)
        press_x[0] = None
        if abs(x1 - x0) < 1e-9:
            # sélection : quel intervalle contient x ?
            xc = x1
            selected.clear()
            for i, (lo, hi) in enumerate(intervals):
                a, b = (lo, hi) if lo <= hi else (hi, lo)
                if a <= xc <= b:
                    selected.append(i)
                    break
            draw_intervals()
            return
        push_history()
        intervals.append([min(x0, x1), max(x0, x1)])
        selected.clear()
        selected.append(len(intervals) - 1)
        draw_intervals()

    def on_key(event) -> None:
        key = (event.key or "").lower()
        if key in ("q", "enter"):
            plt.close(fig)
            return
        if key == "escape":
            selected.clear()
            draw_intervals()
            return
        if key == "u" and history:
            prev = history.pop()
            intervals.clear()
            intervals.extend(deepcopy(prev))
            selected.clear()
            draw_intervals()
            return
        if key == "d" and selected:
            push_history()
            i = selected[0]
            if 0 <= i < len(intervals):
                intervals.pop(i)
            selected.clear()
            draw_intervals()
            return
        if key in ("+", "=", "plus") and selected:
            push_history()
            i = selected[0]
            if 0 <= i < len(intervals):
                lo, hi = intervals[i]
                a, b = (lo, hi) if lo <= hi else (hi, lo)
                w = b - a
                if w > span * 1e-15:
                    tc = _time_at_transit_bottom(a, b)
                    w_new = w * 1.05
                    half = 0.5 * w_new
                    clamped = _clamp_interval_to_plot(tc, half)
                    if clamped is not None:
                        intervals[i] = [clamped[0], clamped[1]]
            draw_intervals()
            return
        if key in ("-", "minus", "underscore") and selected:
            push_history()
            i = selected[0]
            if 0 <= i < len(intervals):
                lo, hi = intervals[i]
                a, b = (lo, hi) if lo <= hi else (hi, lo)
                w = b - a
                w_min = max(span * 1e-8, 1e-9)
                if w > w_min * 1.001:
                    tc = _time_at_transit_bottom(a, b)
                    w_new = max(w * 0.95, w_min)
                    half = 0.5 * w_new
                    clamped = _clamp_interval_to_plot(tc, half)
                    if clamped is not None:
                        intervals[i] = [clamped[0], clamped[1]]
            draw_intervals()
            return

    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("button_release_event", on_release)
    fig.canvas.mpl_connect("key_press_event", on_key)

    draw_intervals()
    # Pas de tight_layout() : incompatible avec l’axe d’aide (set_axis_off + texte encadré) → UserWarning.
    fig.subplots_adjust(left=0.07, right=0.98, top=0.93, bottom=0.05, hspace=0.35)
    plt.show()

    return _intervals_to_mask(time_btjd, [tuple(x) for x in intervals])
