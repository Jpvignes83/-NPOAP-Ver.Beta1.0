"""
Outils de concaténation de courbes (ex. exports TESS / LcTools) pour l’analyse (périodogrammes, TTV).
"""

import pandas as pd
import numpy as np
import os
import glob
import re

TIME_SYNONYMS = (
    "TIME",
    "BJD",
    "JD",
    "BJD_TDB",
    "BTJD",
    "HJD",
    "DATE",
    "MID_TIME",
)
FLUX_SYNONYMS = (
    "DETRENDED_FLUX",
    "PDCSAP_FLUX",
    "CORR_FLUX",
    "FLUX",
    "SAP_FLUX",
    "RAW_FLUX",
    "NORMALIZED",  # ex. Normalized PDCSAP_FLUX (export LcTools)
)


def _normalize_colname(c: str) -> str:
    return (
        str(c)
        .strip()
        .upper()
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("-", "_")
    )


def _pick_time_flux_columns(df: pd.DataFrame):
    cols = [_normalize_colname(c) for c in df.columns]
    df = df.copy()
    df.columns = cols
    col_t = next((c for c in df.columns if any(s in c for s in TIME_SYNONYMS)), None)
    col_f = next((c for c in df.columns if any(s in c for s in FLUX_SYNONYMS)), None)
    return df, col_t, col_f


def _try_read_lctools_style(filepath: str, filename: str) -> pd.DataFrame | None:
    """Fichiers ``tess_lc_fits_to_txt.py --normalized-lctools`` : première ligne #Time (BJD-TDB),..."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            first = f.readline()
    except OSError:
        return None
    s = first.strip()
    if not s.startswith("#") or "time" not in s.lower():
        return None
    header_line = s[1:].strip()
    names = [_normalize_colname(x) for x in header_line.split(",")]
    if not names or not names[0]:
        return None
    try:
        df = pd.read_csv(filepath, skiprows=1, names=names, engine="python", comment="#")
    except Exception:
        return None
    if df.empty or len(df.columns) < 2:
        return None
    df, col_t, col_f = _pick_time_flux_columns(df)
    if not col_t or not col_f:
        print(f"Avertissement : {filename} (entête LcTools) : colonnes temps/flux non reconnues.")
        return None
    out = pd.DataFrame(
        {
            "Time": pd.to_numeric(df[col_t], errors="coerce"),
            "Flux": pd.to_numeric(df[col_f], errors="coerce"),
        }
    )
    out = out.dropna(subset=["Time", "Flux"])
    return out if not out.empty else None


def _try_read_standard_csv(filepath: str, filename: str) -> pd.DataFrame | None:
    """CSV/TXT avec ligne d’en-tête (Time, Flux, …)."""
    try:
        df0 = pd.read_csv(filepath, engine="python", nrows=2)
    except Exception:
        return None
    if df0.empty or len(df0.columns) < 2:
        return None
    df, col_t, col_f = _pick_time_flux_columns(df0)
    if not col_t or not col_f:
        return None
    try:
        df = pd.read_csv(filepath, engine="python")
    except Exception:
        return None
    df, col_t, col_f = _pick_time_flux_columns(df)
    if not col_t or not col_f:
        return None
    out = pd.DataFrame(
        {
            "Time": pd.to_numeric(df[col_t], errors="coerce"),
            "Flux": pd.to_numeric(df[col_f], errors="coerce"),
        }
    )
    out = out.dropna(subset=["Time", "Flux"])
    return out if not out.empty else None


def _skip_lc_basename(filename: str) -> bool:
    return filename in ("concatenated_lightcurve.csv", "mid-time.csv")


def _process_one_lc_file(filepath: str, all_data: list) -> None:
    """Lit un fichier .txt/.csv et ajoute un bloc Time/Flux à *all_data* (même logique que l’ancienne boucle)."""
    filename = os.path.basename(filepath)
    if _skip_lc_basename(filename):
        return
    try:
        df_quick = _try_read_lctools_style(filepath, filename)
        if df_quick is not None:
            mean_flux = np.nanmedian(df_quick["Flux"])
            if np.isfinite(mean_flux) and mean_flux != 0:
                df_quick = df_quick.copy()
                df_quick["Flux"] = df_quick["Flux"] / mean_flux
            all_data.append(df_quick)
            return

        df_quick = _try_read_standard_csv(filepath, filename)
        if df_quick is not None:
            mean_flux = np.nanmedian(df_quick["Flux"])
            if np.isfinite(mean_flux) and mean_flux != 0:
                df_quick = df_quick.copy()
                df_quick["Flux"] = df_quick["Flux"] / mean_flux
            all_data.append(df_quick)
            return

        custom_headers = {}
        data_start_line = 0

        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if line.startswith("#"):
                    match = re.search(r"Column\s*(\d+):\s*(.+)", line)
                    if match:
                        col_index = int(match.group(1)) - 1
                        col_name = (
                            match.group(2)
                            .strip()
                            .upper()
                            .replace("(", "")
                            .replace(")", "")
                            .replace("-", "_")
                        )
                        custom_headers[col_index] = col_name
                elif not line.strip() or line.strip().startswith("2"):
                    data_start_line = i
                    break

        if not custom_headers:
            print(
                f"Avertissement : Fichier {filename} ignoré. "
                f"Pas d’entête LcTools, pas de CSV standard, pas de « Column N: » dans les #."
            )
            return

        df = pd.read_csv(
            filepath,
            sep=r"[,\s]+",
            comment="#",
            header=None,
            skiprows=data_start_line,
            engine="python",
        )

        rename_map = {k: v for k, v in custom_headers.items() if k < df.shape[1]}
        if not rename_map:
            return

        df.rename(columns=rename_map, inplace=True)

        normalized_cols = [
            c.strip().upper().replace("(", "").replace(")", "").replace("-", "_") for c in df.columns
        ]
        df.columns = normalized_cols

        col_t = next((c for c in df.columns if any(syn in c for syn in TIME_SYNONYMS)), None)
        col_f = next((c for c in df.columns if any(syn in c for syn in FLUX_SYNONYMS)), None)

        if not col_t or not col_f:
            print(
                f"Avertissement : Fichier {filename} ignoré. "
                f"Headers TIME ({col_t}) ou FLUX ({col_f}) non trouvés après l'extraction."
            )
            return

        df["Time"] = pd.to_numeric(df[col_t], errors="coerce")
        df["Flux"] = pd.to_numeric(df[col_f], errors="coerce")

        df = df.dropna(subset=["Time", "Flux"])

        if df.empty:
            print(f"Avertissement : Fichier {filename} ignoré (données non valides).")
            return

        mean_flux = np.nanmedian(df["Flux"])
        if np.isfinite(mean_flux) and mean_flux != 0:
            df["Flux"] = df["Flux"] / mean_flux

        all_data.append(df[["Time", "Flux"]].copy())

    except Exception as e:
        print(f"Erreur CRITIQUE lors du traitement de {filename}: {e}. Fichier ignoré.")


def concatenate_lightcurve_paths(file_paths, output_directory=None):
    """
    Concatène une liste explicite de fichiers .txt / .csv (mêmes formats que ``concatenate_lightcurves``).

    Parameters
    ----------
    file_paths : iterable de chemins (str ou Path)
    output_directory : dossier pour ``concatenated_lightcurve.csv`` (défaut : dossier du premier fichier)
    """
    paths = [os.path.abspath(os.path.normpath(p)) for p in file_paths if p and str(p).strip()]
    paths = [p for p in paths if os.path.isfile(p)]
    paths = list(dict.fromkeys(paths))

    if not paths:
        raise ValueError("Aucun fichier .txt ou .csv valide dans la liste.")

    out_dir = output_directory
    if not out_dir:
        out_dir = os.path.dirname(paths[0])
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    all_data = []
    for filepath in paths:
        _process_one_lc_file(filepath, all_data)

    if not all_data:
        raise ValueError("Aucun fichier de light curve valide lu dans la sélection.")

    final_df = pd.concat(all_data, ignore_index=True)
    final_df = final_df.sort_values("Time").reset_index(drop=True)

    output_filepath = os.path.join(out_dir, "concatenated_lightcurve.csv")
    final_df.to_csv(output_filepath, index=False)
    print(f"SUCCESS: Fichier concaténé sauvegardé sous: {output_filepath}")

    return final_df["Time"].values, final_df["Flux"].values


def concatenate_lightcurves(directory):
    """
    Concatène tous les fichiers TXT/CSV du répertoire (headers commentés #, LcTools, etc.).
    """
    directory = os.path.abspath(directory)
    all_files = glob.glob(os.path.join(directory, "*.txt")) + glob.glob(os.path.join(directory, "*.csv"))

    if not all_files:
        raise ValueError("Aucun fichier .txt ou .csv trouvé dans le répertoire.")

    return concatenate_lightcurve_paths(all_files, output_directory=directory)