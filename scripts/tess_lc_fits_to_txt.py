#!/usr/bin/env python3
"""
Convertit une courbe de lumière TESS au format FITS (produit *_lc.fits) en fichier texte.

Usage:
  python scripts/tess_lc_fits_to_txt.py "chemin/vers/tess..._lc.fits"
  python scripts/tess_lc_fits_to_txt.py fichier.fits -o sortie.csv --delimiter comma

Sortie par défaut : même nom que le FITS avec extension .txt (CSV : virgules).
Colonnes : TIME (BTJD), flux PDCSAP (ou SAP si pas de PDCSAP), erreurs, QUALITY si présent.

Mode --normalized-lctools : export type LcTools avec :
  - temps en BJD-TDB (TESS SPOC : TIME en BTJD → BJD-TDB = TIME + 2457000.0) ;
  - normalisation par la moyenne des points hors transit (--t0, --period, --duration-*) ;
  - ou --norm-baseline all pour moyenne sur toute la série (sans masque transit).

Éphemeride **fusionnée** avec ``--planet-name`` (alias ``--exoplanet-eu-name``) :
  **exoplanet.eu** (TAP) + **NASA Exoplanet Archive** (astroquery) ; priorité CLI > EU > NASA
  pour P et T₀ ; **T₁₄** : CLI > NASA (souvent ``pl_trandur`` en heures). Voir ``core/transit_catalog_merge.py``.

Sélection à la souris : ``--pick-transit-windows`` (matplotlib). Peut suivre ``--planet-name`` :
  bandes initiales = éphemeride ; puis **édition** : clic = sélection, **d** supprimer, **+**/**-** élargir/rétrécir,
  glisser = ajouter une bande, **Échap** désélectionner, **u** annuler, **q** valider (voir ``core/lc_transit_pick.py``).

Mode ``--edit-save-fits CHEMIN.fits`` : copie le FITS d’entrée et met à **NaN** les colonnes flux / fond (PDCSAP/SAP/…)
sur les cadences « en transit » (masque catalogue si T₀,P,T₁₄ connus, puis éditeur si ``--pick-transit-windows`` ;
sans éditeur : masque épéméride seule si T₀, P et T₁₄ sont fournis par CLI ou ``--planet-name``).

**Export --normalized-lctools sur un \*\_lc_edited.fits** : les transits sont souvent à **NaN** dans ce FITS ; sans les
valeurs d’origine, ils disparaissent du fichier .txt. Le script tente alors automatiquement de relire PDCSAP/SAP (et SAP_BKG)
depuis le \*\_lc.fits* non édité du même dossier (remplacement ``_lc_edited`` → ``_lc``). Sinon utilisez
``--flux-source-fits``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

try:
    import numpy as np
    from astropy.table import Table
    from astropy.io import fits
except ImportError:
    print("Installez astropy : pip install astropy", file=sys.stderr)
    sys.exit(1)


def _pick_flux_columns(table) -> tuple[str | None, str | None]:
    """Choisit colonnes flux / erreur (priorité PDCSAP)."""
    names = [c.lower() for c in table.colnames]
    colmap = {c.lower(): c for c in table.colnames}

    for flux_k, err_k in [
        ("pdcsap_flux", "pdcsap_flux_err"),
        ("sap_flux", "sap_flux_err"),
        ("flux", "flux_err"),
    ]:
        if flux_k in colmap:
            return colmap[flux_k], colmap.get(err_k)

    return None, None


def _pick_time_column(table) -> str | None:
    for key in ("time", "tim", "t"):
        for c in table.colnames:
            if c.lower() == key:
                return c
    for c in table.colnames:
        if "time" in c.lower() and "corr" not in c.lower():
            return c
    return table.colnames[0] if table.colnames else None


def fits_lc_to_table(fits_path: Path) -> Table:
    """Charge l'extension LIGHTCURVE (TESS / Kepler-like)."""
    with fits.open(fits_path, memmap=False) as hdul:
        # Extension 1 = LIGHTCURVE pour TESS SPOC standard
        for i, hdu in enumerate(hdul):
            if hdu.data is None:
                continue
            if getattr(hdu, "name", "").upper() == "LIGHTCURVE":
                return Table.read(hdul, hdu=i, format="fits")
        # Fallback : première table non vide après le primaire
        for i in range(1, len(hdul)):
            if hdul[i].data is not None and hdul[i].is_image is False:
                return Table.read(hdul, hdu=i, format="fits")
    raise ValueError("Aucune extension LIGHTCURVE / table trouvée dans le FITS.")


def _col_by_name(table: Table, *candidates: str) -> str | None:
    lowmap = {c.lower(): c for c in table.colnames}
    for cand in candidates:
        if cand.lower() in lowmap:
            return lowmap[cand.lower()]
    return None


def _as_float_array(table: Table, col: str) -> np.ndarray:
    x = np.asarray(table[col].data if hasattr(table[col], "data") else table[col], dtype=np.float64)
    if np.ma.isMaskedArray(x):
        x = np.ma.filled(x, np.nan)
    return np.asarray(x, dtype=np.float64)


def _nanmean_positive_denom(a: np.ndarray) -> float:
    m = np.nanmean(a[np.isfinite(a)])
    if m is None or not np.isfinite(m) or m == 0.0:
        raise ValueError("Impossible de normaliser : moyenne du flux non finie ou nulle.")
    return float(m)


# Décalage standard TESS SPOC / HLSP : TIME = BJD_TDB - 2457000.0 (BTJD).
TESS_BTJD_OFFSET = 2457000.0


def btjd_to_bjd_tdb(time_btjd: np.ndarray, offset: float = TESS_BTJD_OFFSET) -> np.ndarray:
    """Convertit le temps catalogue TESS (BTJD) en BJD au temps TDB."""
    return np.asarray(time_btjd, dtype=np.float64) + float(offset)


def _resolve_raw_lc_fits_path(edited_path: Path) -> Path | None:
    """
    Pour une copie « éditée » (flux en transit → NaN), retrouve le FITS d'origine dans le même dossier.

    Gère notamment :
    - SPOC ``..._lc_edited.fits`` → ``..._lc.fits``
    - HLSP / fast-lc ``..._fast-lc_edited.fits`` → ``..._fast-lc.fits``
    - Tout nom se terminant par ``_edited.fits`` / ``_edited.fits.gz`` → même base sans ``_edited``.
    """
    n = edited_path.name
    candidates: list[Path] = []

    if n.endswith("_edited.fits"):
        candidates.append(edited_path.with_name(n[: -len("_edited.fits")] + ".fits"))
    elif n.endswith("_edited.fits.gz"):
        candidates.append(edited_path.with_name(n[: -len("_edited.fits.gz")] + ".fits.gz"))

    # Anciens motifs explicites (si le nom ne finit pas exactement par _edited, ex. doublon rare)
    for old, new in (
        ("_lc_edited.fits.gz", "_lc.fits.gz"),
        ("_lc_edited.fits", "_lc.fits"),
        ("_fast-lc_edited.fits.gz", "_fast-lc.fits.gz"),
        ("_fast-lc_edited.fits", "_fast-lc.fits"),
    ):
        if old in n:
            c = edited_path.with_name(n.replace(old, new))
            if c not in candidates:
                candidates.append(c)

    seen: set[str] = set()
    for c in candidates:
        key = str(c.resolve())
        if key in seen:
            continue
        seen.add(key)
        if c.is_file():
            return c.resolve()
    return None


def try_restore_flux_from_raw_lc(
    input_fits: Path,
    table: Table,
    *,
    flux_source: Path | None = None,
) -> str | None:
    """
    Recopie les colonnes flux (et fond) depuis un FITS non édité aligné sur *table* (même nombre de cadences, TIME identique).
    Permet d'exporter les transits en --normalized-lctools alors que l'entrée est un *_lc_edited.fits* avec NaN en transit.

    Retourne un message pour stderr, ou None si aucune action.
    """
    if flux_source is not None:
        raw_path = Path(flux_source).resolve()
        if not raw_path.is_file():
            return f"Erreur : --flux-source-fits introuvable : {raw_path}"
    else:
        raw_path = _resolve_raw_lc_fits_path(input_fits)
        if raw_path is None or not raw_path.is_file():
            time_col = _pick_time_column(table)
            flux_col, _ = _pick_flux_columns(table)
            if time_col and flux_col:
                flux_cur = _as_float_array(table, flux_col)
                n_bad = int(np.sum(~np.isfinite(flux_cur)))
                if n_bad > 0:
                    return (
                        f"Avertissement : {n_bad} cadences avec flux non fini dans {input_fits.name} ; "
                        f"placez le FITS d'origine (ex. même nom sans _edited) à côté ou utilisez --flux-source-fits."
                    )
            return None

    time_col = _pick_time_column(table)
    flux_col, _ = _pick_flux_columns(table)
    if not time_col or not flux_col:
        return None
    bkg_col = _col_by_name(table, "sap_bkg", "SAP_BKG", "sap_bg", "background")

    try:
        tr = fits_lc_to_table(raw_path)
    except Exception as e:
        return f"Avertissement : impossible de lire {raw_path.name} ({e}), flux non restauré."

    time_col_r = _pick_time_column(tr)
    flux_col_r, _ = _pick_flux_columns(tr)
    if not time_col_r or not flux_col_r:
        return f"Avertissement : {raw_path.name} sans colonnes TIME/flux utilisables."

    time_btjd = _as_float_array(table, time_col)
    time_r = _as_float_array(tr, time_col_r)
    if time_r.shape != time_btjd.shape:
        return (
            f"Avertissement : {raw_path.name} a {time_r.shape[0]} cadences, "
            f"{input_fits.name} en a {time_btjd.shape[0]} — flux non restauré."
        )
    dt = float(np.nanmax(np.abs(time_btjd - time_r)))
    if dt > 1e-5:
        return f"Avertissement : TIME différent entre FITS (max |Δ|={dt}) — flux non restauré."

    flux_r = _as_float_array(tr, flux_col_r)
    try:
        table[flux_col] = flux_r
    except Exception as e:
        return f"Avertissement : impossible d'injecter le flux depuis {raw_path.name} ({e})."

    if bkg_col:
        bkg_col_r = _col_by_name(tr, "sap_bkg", "SAP_BKG", "sap_bg", "background")
        if bkg_col_r:
            try:
                table[bkg_col] = _as_float_array(tr, bkg_col_r)
            except Exception:
                pass

    return (
        f"Flux relu depuis {raw_path.name} ({input_fits.name} avait des NaN en transit — transits conservés pour LcTools)."
    )


def in_transit_mask(time_btjd: np.ndarray, t0_btjd: float, period_days: float, duration_days: float) -> np.ndarray:
    """
    True pour les cadences dont le milieu tombe à l'intérieur d'un transit
    (fenêtre de largeur duration_days centrée sur chaque T0 + n*P).
    t0_btjd doit être dans la même échelle que la colonne TIME du FITS (BTJD).
    """
    if period_days <= 0 or duration_days <= 0:
        raise ValueError("period et duration doivent être strictement positives.")
    half = duration_days * 0.5
    cycles = (np.asarray(time_btjd, dtype=np.float64) - float(t0_btjd)) / float(period_days)
    dist_days = np.abs(cycles - np.round(cycles)) * float(period_days)
    # Inclusif sur les bords : fenêtre [T0 - T14/2, T0 + T14/2] si duration = T14 (jours).
    return dist_days <= half


def write_normalized_lctools_format(
    table: Table,
    out_path: Path,
    *,
    good_quality_only: bool,
    norm_baseline: str,
    t0_btjd: float | None,
    period_days: float | None,
    duration_days: float | None,
    btjd_offset: float,
    in_transit: np.ndarray | None = None,
) -> None:
    """
    Export type LcTools : temps en BJD-TDB, flux / moyenne des points hors transit (par défaut).

    norm_baseline : 'oot' (hors transit) ou 'all' (toute la série valide).
    Si *in_transit* (booléen, une valeur par cadence) est fourni, il définit les transits
    et remplace le masque basé sur T₀ / P / T₁₄.
    """
    time_col = _pick_time_column(table)
    flux_col, _ = _pick_flux_columns(table)
    if not time_col or not flux_col:
        raise ValueError("Colonnes TIME ou PDCSAP/SAP introuvables.")

    bkg_col = _col_by_name(table, "sap_bkg", "SAP_BKG", "sap_bg", "background")
    qual_col = _col_by_name(table, "quality", "QUALITY")

    time_btjd = _as_float_array(table, time_col)
    flux = _as_float_array(table, flux_col)
    bkg = _as_float_array(table, bkg_col) if bkg_col else None

    valid = np.isfinite(time_btjd) & np.isfinite(flux)
    if good_quality_only and qual_col:
        q = np.asarray(table[qual_col].data if hasattr(table[qual_col], "data") else table[qual_col])
        if np.ma.isMaskedArray(q):
            q = np.ma.filled(q, 0)
        q = np.asarray(q, dtype=np.int64)
        valid &= q == 0

    if not np.any(valid):
        raise ValueError("Aucun point valide après filtrage (quality / NaN).")

    if in_transit is not None and norm_baseline == "all":
        raise ValueError("Masque in_transit (souris) incompatible avec --norm-baseline all.")

    if in_transit is not None:
        it = np.asarray(in_transit, dtype=bool)
        if it.shape != (len(time_btjd),):
            raise ValueError(f"in_transit : longueur {it.shape} != {len(time_btjd)} cadences.")
        baseline = valid & ~it
        if not np.any(baseline):
            raise ValueError(
                "Aucun point hors transit pour la normalisation (masque souris / catalogue trop large ?)."
            )
        n_oot = int(np.count_nonzero(baseline))
        if n_oot < 3:
            raise ValueError(f"Trop peu de points hors transit pour une moyenne fiable (n={n_oot}).")
        norm_mask = baseline
    elif norm_baseline == "oot":
        if t0_btjd is None or period_days is None or duration_days is None:
            raise ValueError("Pour --norm-baseline oot, fournissez --t0-btjd, --period-days et --duration-days ou --duration-hours.")
        it = in_transit_mask(time_btjd, t0_btjd, period_days, duration_days)
        baseline = valid & ~it
        if not np.any(baseline):
            raise ValueError(
                "Aucun point hors transit pour la normalisation : vérifiez T0, P, durée ou élargissez la durée."
            )
        n_oot = int(np.count_nonzero(baseline))
        if n_oot < 3:
            raise ValueError(f"Trop peu de points hors transit pour une moyenne fiable (n={n_oot}).")
        norm_mask = baseline
    elif norm_baseline == "all":
        norm_mask = valid
    else:
        raise ValueError(f"norm_baseline inconnu : {norm_baseline!r}")

    mean_flux = _nanmean_positive_denom(flux[norm_mask])
    flux_n = flux / mean_flux

    if bkg is not None:
        mean_bkg = _nanmean_positive_denom(bkg[norm_mask])
        bkg_n = bkg / mean_bkg
        header = "#Time (BJD-TDB),Normalized PDCSAP_FLUX,Normalized SAP_BKG"
        ncols = 3
    else:
        bkg_n = None
        header = "#Time (BJD-TDB),Normalized PDCSAP_FLUX"
        ncols = 2

    time_out = btjd_to_bjd_tdb(time_btjd, offset=btjd_offset)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(header + "\n")
        for i in range(len(time_out)):
            if ncols == 3:
                if not np.isfinite(time_out[i]) or not np.isfinite(flux_n[i]) or not np.isfinite(bkg_n[i]):
                    continue
                f.write(
                    f"{time_out[i]:.17f},{flux_n[i]:.14f},{bkg_n[i]:.14f}\n"
                )
            else:
                if not np.isfinite(time_out[i]) or not np.isfinite(flux_n[i]):
                    continue
                f.write(f"{time_out[i]:.17f},{flux_n[i]:.14f}\n")


def apply_in_transit_nan_to_fits(fits_in: Path, fits_out: Path, in_transit: np.ndarray) -> None:
    """
    Copie *fits_in* vers *fits_out* et remplace par NaN les colonnes de flux (et SAP_BKG) sur les cadences
    où *in_transit* est True.
    """
    import shutil

    fits_in = fits_in.resolve()
    fits_out = fits_out.resolve()
    if fits_in == fits_out:
        raise ValueError("Le fichier FITS de sortie doit être différent de l’entrée (copie + modification).")
    fits_out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(fits_in, fits_out)

    with fits.open(fits_out, mode="update", memmap=False) as hdul:
        lc_hdu = None
        for hdu in hdul:
            if hdu.data is None:
                continue
            if getattr(hdu, "name", "").upper() == "LIGHTCURVE":
                lc_hdu = hdu
                break
        if lc_hdu is None:
            raise ValueError("Extension LIGHTCURVE introuvable dans le FITS.")

        data = lc_hdu.data
        n = len(data)
        m = np.asarray(in_transit, dtype=bool).reshape(-1)
        if m.shape[0] != n:
            raise ValueError(f"Longueur du masque transit ({m.shape[0]}) != nombre de cadences ({n}).")

        for name in data.dtype.names:
            u = name.upper()
            if "FLUX" in u or "SAP_BKG" in u:
                if not np.issubdtype(data.dtype[name], np.floating):
                    continue
                data[name][m] = np.nan

        lc_hdu.header.add_history("NPOAP: cadences en transit -> NaN (flux / SAP_BKG).")
    # context manager closes and flushes


def main() -> int:
    p = argparse.ArgumentParser(description="TESS/Kepler *_lc.fits → fichier texte (série temporelle)")
    p.add_argument("fits_file", type=Path, help="Fichier *_lc.fits")
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Fichier de sortie (.txt, .csv, .tsv). Défaut : même nom que le FITS + .txt",
    )
    p.add_argument(
        "--delimiter",
        choices=("comma", "tab", "space"),
        default="comma",
        help="comma = CSV, tab = TSV, space = colonnes séparées par espaces",
    )
    p.add_argument(
        "--all-columns",
        action="store_true",
        help="Exporter toutes les colonnes numériques de l'extension LIGHTCURVE",
    )
    p.add_argument(
        "--normalized-lctools",
        action="store_true",
        help=(
            "Export type LcTools : BJD-TDB + PDCSAP normalisé (moyenne hors transit par défaut) + SAP_BKG si présent"
        ),
    )
    p.add_argument(
        "--good-quality-only",
        action="store_true",
        help="Pour --normalized-lctools : ne garder que QUALITY==0 pour la baseline et l'export",
    )
    p.add_argument(
        "--norm-baseline",
        choices=("oot", "all"),
        default="oot",
        help=(
            "oot = normaliser sur la moyenne des points hors transit (défaut, exige --t0-btjd, --period-days, durée) ; "
            "all = moyenne sur tous les points valides"
        ),
    )
    p.add_argument(
        "--t0-btjd",
        type=float,
        default=None,
        help="Temps de conjonction inférieure T0 en BTJD (même échelle que la colonne TIME du FITS TESS)",
    )
    p.add_argument(
        "--period-days",
        type=float,
        default=None,
        help="Période orbitale (jours), pour exclure les transits de la baseline",
    )
    dur = p.add_mutually_exclusive_group()
    dur.add_argument(
        "--duration-days",
        type=float,
        default=None,
        help="Durée totale du transit (jours), pour le masque hors transit",
    )
    dur.add_argument(
        "--duration-hours",
        type=float,
        default=None,
        help="Durée totale du transit (heures), convertie en jours",
    )
    p.add_argument(
        "--btjd-offset",
        type=float,
        default=TESS_BTJD_OFFSET,
        help=f"Décalage BJD-TDB = TIME + offset (défaut TESS SPOC : {TESS_BTJD_OFFSET})",
    )
    p.add_argument(
        "--planet-name",
        "--exoplanet-eu-name",
        dest="catalog_planet_name",
        type=str,
        default=None,
        metavar="NOM",
        help=(
            "Avec --normalized-lctools : fusion exoplanet.eu (TAP) + NASA Archive (P, T₀, T₁₄ si dispo). "
            "Priorité : CLI > EU > NASA. Combinable avec --pick-transit-windows (catalogue puis affinage souris)."
        ),
    )
    p.add_argument(
        "--exoplanet-eu-timeout",
        type=float,
        default=120.0,
        help="Délai HTTP max (s) pour la requête TAP exoplanet.eu",
    )
    p.add_argument(
        "--skip-exoplanet-eu",
        action="store_true",
        help="Avec --planet-name : ne pas interroger exoplanet.eu (NASA / CLI seuls).",
    )
    p.add_argument(
        "--skip-nasa-ephemeris",
        action="store_true",
        help="Avec --planet-name : ne pas interroger la NASA (exoplanet.eu / CLI seuls).",
    )
    p.add_argument(
        "--pick-transit-windows",
        action="store_true",
        help=(
            "Éditeur matplotlib des masques transit : clic = sélection, d = supprimer, +/- = élargir/rétrécir, "
            "glisser = nouvelle bande, Échap = désélection, u = annuler, q/Entrée = valider. "
            "Avec --planet-name : bandes initiales depuis l’éphemeride puis édition. Nécessite matplotlib."
        ),
    )
    p.add_argument(
        "--edit-save-fits",
        type=Path,
        default=None,
        metavar="OUT.fits",
        help=(
            "Copie le FITS d’entrée vers OUT.fits et met à NaN flux / SAP_BKG sur les cadences « en transit ». "
            "Utiliser --pick-transit-windows pour l’éditeur ; sinon T₀, P et T₁₄ requis (CLI ou --planet-name) "
            "pour le masque épéméride seul."
        ),
    )
    p.add_argument(
        "--flux-source-fits",
        type=Path,
        default=None,
        metavar="CHEMIN_lc.fits",
        help=(
            "Avec --normalized-lctools : relire PDCSAP/SAP (et SAP_BKG) depuis ce FITS au lieu de l’entrée. "
            "Utile si l’entrée est *_lc_edited.fits* (NaN en transit). Sinon le script cherche automatiquement *_lc.fits*."
        ),
    )
    args = p.parse_args()

    if args.normalized_lctools and args.edit_save_fits is not None:
        print("--normalized-lctools et --edit-save-fits sont incompatibles.", file=sys.stderr)
        return 1

    if args.catalog_planet_name and not args.normalized_lctools and args.edit_save_fits is None:
        print("--planet-name nécessite --normalized-lctools ou --edit-save-fits.", file=sys.stderr)
        return 1
    if args.pick_transit_windows and not args.normalized_lctools and args.edit_save_fits is None:
        print("--pick-transit-windows nécessite --normalized-lctools ou --edit-save-fits.", file=sys.stderr)
        return 1

    fits_path = args.fits_file.resolve()
    if not fits_path.is_file():
        print(f"Fichier introuvable : {fits_path}", file=sys.stderr)
        return 1

    try:
        t = fits_lc_to_table(fits_path)
    except Exception as e:
        print(f"Lecture FITS impossible : {e}", file=sys.stderr)
        return 1

    out = args.output
    if out is None:
        out = fits_path.with_suffix(".txt")

    # --- Enregistrer un FITS édité (masque transit -> NaN sur les flux) ---
    if args.edit_save_fits is not None:
        if args.skip_exoplanet_eu and args.skip_nasa_ephemeris and args.catalog_planet_name:
            print(
                "--planet-name avec --skip-exoplanet-eu et --skip-nasa-ephemeris : aucune source catalogue.",
                file=sys.stderr,
            )
            return 1

        duration_days = None
        if args.duration_hours is not None:
            duration_days = args.duration_hours / 24.0
        elif args.duration_days is not None:
            duration_days = args.duration_days

        period_days = args.period_days
        t0_btjd = args.t0_btjd

        if args.catalog_planet_name:
            try:
                from core.transit_catalog_merge import merge_transit_ephemeris
            except ImportError as e:
                print(f"Module core.transit_catalog_merge introuvable : {e}", file=sys.stderr)
                return 1
            m = merge_transit_ephemeris(
                args.catalog_planet_name.strip(),
                cli_period_days=args.period_days,
                cli_t0_btjd=args.t0_btjd,
                cli_duration_days=duration_days,
                btjd_offset=args.btjd_offset,
                eu_timeout_s=float(args.exoplanet_eu_timeout),
                query_exoplanet_eu=not args.skip_exoplanet_eu,
                query_nasa=not args.skip_nasa_ephemeris,
            )
            period_days = m.period_days
            t0_btjd = m.t0_btjd
            duration_days = m.duration_days
            print(
                f"Catalogues fusionnés : P={period_days} (src={m.sources.get('period')}), "
                f"T0_BTJD={t0_btjd} (src={m.sources.get('t0')}), "
                f"T14_jours={duration_days} (src={m.sources.get('duration')}).",
                file=sys.stderr,
            )

        time_col = _pick_time_column(t)
        flux_col, _ = _pick_flux_columns(t)
        if not time_col or not flux_col:
            print("Colonnes TIME ou flux introuvables.", file=sys.stderr)
            return 1
        time_btjd = _as_float_array(t, time_col)
        flux_arr = _as_float_array(t, flux_col)

        base_it = None
        if period_days is not None and t0_btjd is not None and duration_days is not None:
            try:
                base_it = in_transit_mask(time_btjd, t0_btjd, period_days, duration_days)
            except Exception as e:
                print(f"Masque épéméride impossible : {e}", file=sys.stderr)
                return 1

        in_transit_mask_arr: np.ndarray | None = None
        if args.pick_transit_windows:
            try:
                from core.lc_transit_pick import pick_in_transit_mask
            except ImportError as e:
                print(f"Sélection souris impossible (lc_transit_pick) : {e}", file=sys.stderr)
                return 1
            try:
                in_transit_mask_arr = pick_in_transit_mask(
                    time_btjd,
                    flux_arr,
                    base_in_transit=base_it,
                    title=f"Édition LC (FITS) — {fits_path.name}",
                )
            except Exception as e:
                print(f"Sélection interactive échouée : {e}", file=sys.stderr)
                return 1
            n_tot = int(np.count_nonzero(in_transit_mask_arr))
            if base_it is not None:
                n_cat = int(np.count_nonzero(base_it))
                print(
                    f"Masque édité : {n_tot} cadences « en transit » "
                    f"(éphemeride initiale ~{n_cat} cadences, puis modifications interactives).",
                    file=sys.stderr,
                )
            else:
                print(f"Masque édité : {n_tot} cadences « en transit ».", file=sys.stderr)
        else:
            if base_it is None:
                print(
                    "Sans --pick-transit-windows : fournissez T₀, P et T₁₄ (CLI ou --planet-name) "
                    "pour appliquer le masque épéméride seul.",
                    file=sys.stderr,
                )
                return 1
            in_transit_mask_arr = base_it
            print(
                f"Masque épéméride seul : {int(np.count_nonzero(in_transit_mask_arr))} cadences « en transit ».",
                file=sys.stderr,
            )

        try:
            apply_in_transit_nan_to_fits(fits_path, args.edit_save_fits, in_transit_mask_arr)
        except Exception as e:
            print(f"Écriture FITS édité impossible : {e}", file=sys.stderr)
            return 1
        print(f"FITS édité enregistré : {args.edit_save_fits.resolve()}", file=sys.stderr)
        return 0

    if args.normalized_lctools:
        if args.all_columns:
            print("--normalized-lctools est incompatible avec --all-columns", file=sys.stderr)
            return 1
        if args.skip_exoplanet_eu and args.skip_nasa_ephemeris and args.catalog_planet_name:
            print("--planet-name avec --skip-exoplanet-eu et --skip-nasa-ephemeris : aucune source catalogue.", file=sys.stderr)
            return 1

        restore_msg = try_restore_flux_from_raw_lc(
            fits_path, t, flux_source=args.flux_source_fits
        )
        if restore_msg:
            print(restore_msg, file=sys.stderr)
            if restore_msg.startswith("Erreur :"):
                return 1

        duration_days = None
        if args.duration_hours is not None:
            duration_days = args.duration_hours / 24.0
        elif args.duration_days is not None:
            duration_days = args.duration_days

        period_days = args.period_days
        t0_btjd = args.t0_btjd
        in_transit_mask_arr = None

        if args.catalog_planet_name:
            if args.norm_baseline == "all":
                print(
                    "Note : --planet-name ignoré avec --norm-baseline all (aucune éphemeride requise).",
                    file=sys.stderr,
                )
            else:
                try:
                    from core.transit_catalog_merge import merge_transit_ephemeris
                except ImportError as e:
                    print(f"Module core.transit_catalog_merge introuvable : {e}", file=sys.stderr)
                    return 1
                m = merge_transit_ephemeris(
                    args.catalog_planet_name.strip(),
                    cli_period_days=args.period_days,
                    cli_t0_btjd=args.t0_btjd,
                    cli_duration_days=duration_days,
                    btjd_offset=args.btjd_offset,
                    eu_timeout_s=float(args.exoplanet_eu_timeout),
                    query_exoplanet_eu=not args.skip_exoplanet_eu,
                    query_nasa=not args.skip_nasa_ephemeris,
                )
                period_days = m.period_days
                t0_btjd = m.t0_btjd
                duration_days = m.duration_days
                print(
                    f"Catalogues fusionnés : P={period_days} (src={m.sources.get('period')}), "
                    f"T0_BTJD={t0_btjd} (src={m.sources.get('t0')}), "
                    f"T14_jours={duration_days} (src={m.sources.get('duration')}).",
                    file=sys.stderr,
                )

        if args.pick_transit_windows:
            if args.norm_baseline == "all":
                print("--pick-transit-windows nécessite --norm-baseline oot (défaut).", file=sys.stderr)
                return 1
            try:
                from core.lc_transit_pick import pick_in_transit_mask
            except ImportError as e:
                print(f"Sélection souris impossible (lc_transit_pick) : {e}", file=sys.stderr)
                return 1
            time_col = _pick_time_column(t)
            flux_col, _ = _pick_flux_columns(t)
            if not time_col or not flux_col:
                print("Colonnes TIME ou flux introuvables pour l'affichage.", file=sys.stderr)
                return 1
            time_btjd = _as_float_array(t, time_col)
            flux_arr = _as_float_array(t, flux_col)

            base_it = None
            if args.norm_baseline != "all":
                if t0_btjd is not None and period_days is not None and duration_days is not None:
                    base_it = in_transit_mask(time_btjd, t0_btjd, period_days, duration_days)
                elif args.catalog_planet_name:
                    print(
                        "Avec --planet-name et --pick-transit-windows, il faut T₀, P et T₁₄ complets "
                        "(catalogues + CLI) pour tracer le masque catalogue avant affinage.",
                        file=sys.stderr,
                    )
                    return 1

            try:
                in_transit_mask_arr = pick_in_transit_mask(
                    time_btjd,
                    flux_arr,
                    base_in_transit=base_it,
                    title=f"Transits (baseline OOT) — {fits_path.name}",
                )
            except Exception as e:
                print(f"Sélection interactive échouée : {e}", file=sys.stderr)
                return 1
            n_tot = int(np.count_nonzero(in_transit_mask_arr))
            if base_it is not None:
                n_cat = int(np.count_nonzero(base_it))
                print(
                    f"Masque édité : {n_tot} cadences « en transit » "
                    f"(éphemeride initiale ~{n_cat} cadences, puis modifications interactives).",
                    file=sys.stderr,
                )
            else:
                print(f"Masque édité : {n_tot} cadences « en transit ».", file=sys.stderr)

        if args.norm_baseline == "oot" and in_transit_mask_arr is None:
            if t0_btjd is None or period_days is None or duration_days is None:
                print(
                    "Avec --norm-baseline oot : fournissez T₀, P et T₁₄ (CLI), ou --planet-name (EU+NASA), "
                    "ou --pick-transit-windows. Sinon --norm-baseline all.",
                    file=sys.stderr,
                )
                return 1
        try:
            write_normalized_lctools_format(
                t,
                out,
                good_quality_only=args.good_quality_only,
                norm_baseline=args.norm_baseline,
                t0_btjd=t0_btjd,
                period_days=period_days,
                duration_days=duration_days,
                btjd_offset=args.btjd_offset,
                in_transit=in_transit_mask_arr,
            )
        except Exception as e:
            print(f"Export normalisé impossible : {e}", file=sys.stderr)
            print(f"Colonnes disponibles : {t.colnames}", file=sys.stderr)
            return 1
        print(f"Écrit (format LcTools normalisé) : {out}")
        return 0

    if args.all_columns:
        sub = t
    else:
        time_col = _pick_time_column(t)
        flux_col, err_col = _pick_flux_columns(t)
        cols = []
        if time_col:
            cols.append(time_col)
        if flux_col:
            cols.append(flux_col)
        if err_col:
            cols.append(err_col)
        if "quality" in [c.lower() for c in t.colnames]:
            q = next(c for c in t.colnames if c.lower() == "quality")
            cols.append(q)
        # Retirer doublons en gardant l'ordre
        seen = set()
        cols = [c for c in cols if c and c not in seen and not seen.add(c)]
        missing = [c for c in cols if c not in t.colnames]
        if missing:
            print(f"Colonnes manquantes : {missing}. Colonnes disponibles : {t.colnames}", file=sys.stderr)
            return 1
        sub = t[cols]

    fmt_map = {"comma": "ascii.csv", "tab": "ascii.tab", "space": "ascii.basic"}
    fmt = fmt_map[args.delimiter]

    try:
        sub.write(str(out), format=fmt, overwrite=True)
    except Exception as e:
        print(f"Écriture impossible : {e}", file=sys.stderr)
        return 1

    print(f"Écrit : {out} ({len(sub)} lignes, {len(sub.colnames)} colonnes)")
    print(f"Colonnes : {', '.join(sub.colnames)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
