# core/alcdef_export.py
"""
Lecture de light_curve.txt (NPOAP) et export vers blocs S-ALCDEF.

Le format S-ALCDEF (Simple ALCDEF) est décrit sur https://alcdef.org/
(cf. « Brief History II », blocs STARTBLOCK / STARTDATA / ENDDATA / ENDBLOCK).
Les métadonnées complètes se souvent via le formulaire du site après génération des données.

Les magnitudes sont obtenues par m = -2.5 * log10(f_rel) où f_rel est le flux relatif
(cible / étoiles de comparaison, issues de la photométrie NPOAP). Il s'agit de magnitudes
différentielles à zéro arbitraire ; une transformation vers un système standard doit être
indiquée dans le formulaire de soumission / le rapport de publication.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

# Première ligne de light_curve.txt lorsque les données sont en magnitude (3 déc.) et non en flux.
LIGHT_CURVE_MAGNITUDE_HEADER = (
    "# NPOAP LIGHT_CURVE_MAGNITUDE — JD-UTC, mag (diff. -2.5*log10(fn)), mag_err (3 dec.)\n"
)

def _time_column_looks_like_rotation_phase(t: np.ndarray) -> bool:
    """
    Heuristique pour les exports « phase de rotation » (ex. DAMIT) :
    abscisse dans ~[0, 1] (sans en-tête magnitude NPOAP).
    """
    arr = np.asarray(t, dtype=float)
    if arr.size < 5:
        return False
    lo = float(np.nanmin(arr))
    hi = float(np.nanmax(arr))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return False
    if lo < -0.02 or hi > 1.02:
        return False
    return True


def damit_model_flux_at_phase(
    phase_obs: np.ndarray,
    damit_phase: np.ndarray,
    damit_flux: np.ndarray,
    shift: float,
) -> np.ndarray:
    """
    Évalue le modèle DAMIT sur les phases d'observation ``phase_obs`` (0–1),
    avec déphasage ``shift`` : valeur = D(mod(φ + shift, 1)).
    """
    phi = np.mod(np.asarray(phase_obs, dtype=float) + float(shift), 1.0)
    dph = np.asarray(damit_phase, dtype=float).ravel()
    dfl = np.asarray(damit_flux, dtype=float).ravel()
    if dph.size < 2:
        return np.full_like(phi, np.nan, dtype=float)
    o = np.argsort(dph)
    dph_s = dph[o]
    dfl_s = dfl[o]
    ph_ext = np.concatenate([dph_s - 1.0, dph_s, dph_s + 1.0])
    fl_ext = np.concatenate([dfl_s, dfl_s, dfl_s])
    return np.interp(phi, ph_ext, fl_ext, left=np.nan, right=np.nan)


def load_damit_light_curve_model(path: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Charge un ``light_curve.txt`` exporté par DAMIT : phase (0–1), flux relatif.

    Deux colonnes numériques attendues (pas le format magnitude NPOAP).
    Retourne ``(phase, flux)`` triés par phase.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with open(path, encoding="utf-8", errors="replace") as fp:
        first_line = fp.readline()
    if "LIGHT_CURVE_MAGNITUDE" in first_line.upper():
        raise ValueError("Fichier magnitude NPOAP : utilisez une courbe d'observations ou un export DAMIT brut.")
    try:
        data = np.loadtxt(path, comments="#")
    except Exception as e:
        raise ValueError(f"Lecture impossible : {path}") from e
    data = np.atleast_2d(data)
    if data.shape[1] != 2:
        raise ValueError(
            "Le modèle DAMIT attendu a exactement 2 colonnes : phase (0–1), flux relatif."
        )
    phase = np.asarray(data[:, 0], dtype=float)
    flux = np.asarray(data[:, 1], dtype=float)
    if not _time_column_looks_like_rotation_phase(phase):
        raise ValueError(
            "La première colonne ne ressemble pas à une phase de rotation dans [0, 1]."
        )
    o = np.argsort(phase)
    phase = phase[o]
    flux = flux[o]
    return phase, flux


__all__ = [
    "load_light_curve_txt",
    "load_damit_light_curve_model",
    "damit_model_flux_at_phase",
    "LIGHT_CURVE_MAGNITUDE_HEADER",
    "relative_flux_to_differential_magnitude",
    "build_simple_alcdef_text",
    "build_submission_report_lines",
]


def load_light_curve_txt(path: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Charge un fichier light_curve.txt NPOAP (temps en JD-UTC).

    Formats supportés :
    - Flux : Time (JD), Relative_flux_fn, relative_flux_fn_err (éventuelle ligne d'en-tête texte).
    - Magnitude : première ligne ``# ... LIGHT_CURVE_MAGNITUDE ...`` puis JD-UTC, mag, mag_err ;
      conversion vers flux linéaire : f = 10**(-0.4 * m).

    Les exports « phase + flux » du site DAMIT ne sont pas des courbes d'observation :
    utilisez ``load_damit_light_curve_model`` et la surimpression dans l'interface.

    Retourne ``(jd, flux, flux_err)`` ; ``flux_err`` est None si absent.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with open(path, encoding="utf-8", errors="replace") as fp:
        first_line = fp.readline()
    is_mag_format = "LIGHT_CURVE_MAGNITUDE" in first_line.upper()

    try:
        data = np.loadtxt(path, comments="#")
    except Exception:
        data = np.loadtxt(path, skiprows=1, comments="#")
    data = np.atleast_2d(data)
    if data.shape[1] < 2:
        raise ValueError(f"Fichier courbe invalide (colonnes insuffisantes) : {path}")
    jd = np.asarray(data[:, 0], dtype=float)

    if is_mag_format:
        mag = np.asarray(data[:, 1], dtype=float)
        flux = np.full_like(mag, np.nan, dtype=float)
        ok = np.isfinite(mag)
        flux[ok] = np.power(10.0, -0.4 * mag[ok])
        flux_err: Optional[np.ndarray]
        if data.shape[1] >= 3:
            smag = np.asarray(data[:, 2], dtype=float)
            flux_err = np.full_like(flux, np.nan, dtype=float)
            o2 = ok & np.isfinite(smag) & np.isfinite(flux) & (flux > 0)
            # σ_f = 0,4 · ln(10) · f · σ_m
            flux_err[o2] = flux[o2] * (0.4 * math.log(10.0)) * np.abs(smag[o2])
        else:
            flux_err = None
        return jd, flux, flux_err

    flux = np.asarray(data[:, 1], dtype=float)
    if data.shape[1] >= 3:
        flux_err = np.asarray(data[:, 2], dtype=float)
    else:
        flux_err = None
    if data.shape[1] == 2 and _time_column_looks_like_rotation_phase(jd):
        raise ValueError(
            "Ce fichier ressemble à un export DAMIT (phase 0–1 + flux relatif), pas à une courbe "
            "NPOAP en jours julien. Dans l’onglet courbe, utilisez « Charger modèle DAMIT… » pour le "
            "superposer à vos observations."
        )
    return jd, flux, flux_err


def relative_flux_to_differential_magnitude(
    flux: np.ndarray,
    flux_err: Optional[np.ndarray],
    *,
    normalize_median: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Passe du flux relatif à une magnitude différentielle : m = -2.5 * log10(f).

    Parameters
    ----------
    normalize_median : bool
        Si True, divise les flux par la médiane (après filtrage des valeurs > 0) pour
        centrer la courbe autour de 0 mag en moyenne.
    """
    f = np.asarray(flux, dtype=float)
    mask = np.isfinite(f) & (f > 0)
    if normalize_median and np.any(mask):
        m0 = float(np.median(f[mask]))
        if m0 > 0:
            f = f / m0
            if flux_err is not None:
                fe = np.asarray(flux_err, dtype=float) / m0
            else:
                fe = None
        else:
            fe = np.asarray(flux_err, dtype=float) if flux_err is not None else None
    else:
        fe = np.asarray(flux_err, dtype=float) if flux_err is not None else None

    mag = np.full_like(f, np.nan, dtype=float)
    safe = np.isfinite(f) & (f > 0)
    mag[safe] = -2.5 * np.log10(f[safe])

    mag_err = np.full_like(f, np.nan, dtype=float)
    if fe is not None:
        ok = safe & np.isfinite(fe) & (f > 0)
        mag_err[ok] = (2.5 / math.log(10.0)) * np.abs(fe[ok] / f[ok])
    else:
        mag_err = np.full_like(f, np.nan, dtype=float)

    return mag, mag_err


def build_simple_alcdef_text(
    jd: np.ndarray,
    mag: np.ndarray,
    mag_err: np.ndarray,
    *,
    delimiter: str = "|",
    comment_lines: Optional[List[str]] = None,
    alcdef_metadata_lines: Optional[List[str]] = None,
) -> str:
    """
    Assemble un ou plusieurs blocs S-ALCDEF (un seul bloc ici : une session fichier).
    Les lignes commençant par # sont ignorées par le parseur S-ALCDEF.

    ``alcdef_metadata_lines`` : lignes ``CLÉ=valeur`` (ex. EQUINOX, COMPNAME1, …)
    insérées **après** ``STARTBLOCK`` et **avant** ``STARTDATA``, conformément au format ALCDEF.
    """
    out: List[str] = []
    out.append("# Fichier généré par NPOAP — soumission : https://alcdef.org/ (Simple ALCDEF)")
    out.append("# Vérifier avec ALCDEFVerify avant envoi. Métadonnées : formulaire web ALCDEF.")
    if comment_lines:
        for c in comment_lines:
            for part in str(c).splitlines():
                out.append("# " + part.strip())
    out.append("STARTBLOCK")
    if alcdef_metadata_lines:
        for raw in alcdef_metadata_lines:
            line = str(raw).strip()
            if not line or line.startswith("#"):
                continue
            out.append(line)
    out.append("STARTDATA")
    for t, m, me in zip(jd, mag, mag_err):
        if not np.isfinite(t) or not np.isfinite(m):
            continue
        if np.isfinite(me):
            out.append(f"{t:.6f}{delimiter}{m:+.6f}{delimiter}{me:.6f}")
        else:
            out.append(f"{t:.6f}{delimiter}{m:+.6f}")
    out.append("ENDDATA")
    out.append("ENDBLOCK")
    out.append("")
    return "\n".join(out)


@dataclass
class ReportMeta:
    """Informations rédigées par l'utilisateur (publication / LCDB / formulaire ALCDEF)."""

    asteroid_number: str = ""
    asteroid_name: str = ""
    observers: str = ""
    telescope: str = ""
    instrument: str = ""
    filter_name: str = ""
    mag_band: str = ""
    mpc_code: str = ""
    comment: str = ""


def build_submission_report_lines(
    jd: np.ndarray,
    flux: np.ndarray,
    flux_err: Optional[np.ndarray],
    mag: np.ndarray,
    mag_err: np.ndarray,
    meta: ReportMeta,
    *,
    fit: Optional[Dict[str, Any]] = None,
    alcdef_comparator_lines: Optional[List[str]] = None,
) -> str:
    """
    Texte d'aide pour rédaction (Minor Planet Bulletin, notes LCDB, commentaire ALCDEF).
    Ce n'est pas un fichier caméra-prête : l'auteur doit suivre le guide des auteurs MPB.

    ``alcdef_comparator_lines`` : copie des métadonnées COMP* / EQUINOX (Gaia DR3) figées à la validation.
    """
    jd = np.asarray(jd, dtype=float)
    mag = np.asarray(mag, dtype=float)
    lines: List[str] = []
    lines.append("=== Notes de soumission (NPOAP) — astéroïde / photométrie ===")
    lines.append("")
    oid = " ".join(x for x in (meta.asteroid_number.strip(), meta.asteroid_name.strip()) if x).strip()
    lines.append(f"Objet : {oid or '(à compléter)'}")
    lines.append(f"Observateurs : {meta.observers or '(à compléter)'}")
    lines.append(f"Lieu / code MPC : {meta.mpc_code or '(à compléter)'}")
    lines.append(f"Télescope : {meta.telescope or '(à compléter)'}")
    lines.append(f"Détecteur / caméra : {meta.instrument or '(à compléter)'}")
    lines.append(f"Filtre (instrument) : {meta.filter_name or '(à compléter)'}")
    lines.append(f"MagBand ALCDEF / système : {meta.mag_band or '(à compléter — voir standard ALCDEF 2024+ ; Gaia G = GG)'}")
    lines.append("")
    if alcdef_comparator_lines is not None:
        if alcdef_comparator_lines:
            lines.append("Métadonnées ALCDEF — comparateurs (Gaia DR3, résolus à la validation) :")
            for ln in alcdef_comparator_lines:
                lines.append(f"  {ln}")
        else:
            lines.append(
                "Comparateurs ALCDEF : non enregistrés (aucun comparateur C* dans la session "
                "photométrie au moment de « Valider pour publication » — pas de lignes COMP / EQUINOX)."
            )
        lines.append("")
    if len(jd):
        lines.append(
            f"Période couverte (JD UTC) : {float(np.min(jd)):.6f} – {float(np.max(jd)):.6f} "
            f"({len(jd)} points)"
        )
    lines.append(
        f"Amplitude indicative (mag diff.) : {float(np.nanmax(mag) - np.nanmin(mag)):.3f} mag "
        "(non corrigée phase angle ; à affiner pour publication)"
    )
    if fit and fit.get("success") and fit.get("P"):
        lines.append(
            f"Ajustement modèle cosinus (outil NPOAP) : P = {fit['P']:.6f} j, "
            f"A = {fit.get('A', float('nan')):.4f}, t0 = {fit.get('t0', float('nan')):.6f} JD ; "
            f"chi²={fit.get('chi2', float('nan')):.2f}"
        )
    lines.append("")
    lines.append("Données réduites : magnitudes différentielles m = -2.5 log10(flux relatif NPOAP).")
    lines.append("Archivage recommandé : dépôt ALCDEF (https://alcdef.org/) + citation dans l'article.")
    lines.append("MPB : modèles et guide auteurs sur https://mpbulletin.org/")
    if meta.comment.strip():
        lines.append("")
        lines.append("Commentaire :")
        lines.append(meta.comment.strip())
    lines.append("")
    return "\n".join(lines)
