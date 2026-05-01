# core/ExoMiner_analysis.py
"""
Pont subprocess vers le pipeline ExoMiner (NASA) : exécute ``exominer_pipeline/run_pipeline.py``
depuis la racine du dépôt, avec un interpréteur Python dédié.

Variables d'environnement :
    NPOAP_INSTALL_ROOT    Racine NPOAP (GUI) ; clone attendu sous ``external_apps/Exominer``.
    NPOAP_EXOMINER_ROOT   Racine du dépôt ExoMiner (contient ``exominer_pipeline/run_pipeline.py``).
    NPOAP_EXOMINER_PYTHON Interpréteur Python pour ce pipeline (défaut : ``sys.executable``).

Référence amont : https://github.com/nasa/Exominer
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

logger = logging.getLogger(__name__)

EXOMINER_REPO_URL = "https://github.com/nasa/Exominer"
# Chemins relatifs à la racine du clone (pipeline_run_config.yaml / config_normalize_data.yaml).
_EXOMINER_MODEL_RELPATHS = {
    "exominer++_single": "exominer_pipeline/data/exominer-plusplus_cv-iter0-model0_tess-spoc-2min-s1s67_tess-kepler.keras",
    "exominer++_cviter-mean-ensemble": "exominer_pipeline/data/exominer-plusplus_cv-iter0-ensemble_tess-spoc-2min-s1s67_tess.keras",
    "exominer++_cv-super-mean-ensemble": "exominer_pipeline/data/exominer-plusplus_cv-mean-ensemble_tess-spoc-2min-s1s67_tess-kepler.keras",
}
_EXOMINER_NORM_STATS_RELPATHS = (
    "exominer_pipeline/data/norm_stats/train_scalarparam_norm_stats.npy",
    "exominer_pipeline/data/norm_stats/train_centroid_norm_stats.npy",
    "exominer_pipeline/data/norm_stats/train_diffimg_norm_stats.npy",
)


def resolve_exominer_root(root: Optional[Union[str, Path]] = None) -> Optional[Path]:
    """Renvoie le chemin racine du clone ExoMiner, ou None si non défini."""
    raw = root if root is not None else os.environ.get("NPOAP_EXOMINER_ROOT", "").strip()
    if not raw:
        return None
    p = Path(raw).expanduser().resolve()
    return p


def resolve_exominer_python(explicit: Optional[Union[str, Path]] = None) -> str:
    """Python utilisé pour lancer ``run_pipeline.py`` (TensorFlow / deps NASA)."""
    if explicit is not None:
        return str(Path(explicit).expanduser())
    env_py = os.environ.get("NPOAP_EXOMINER_PYTHON", "").strip()
    if env_py:
        return str(Path(env_py).expanduser())
    return sys.executable


def pipeline_script_path(root: Path) -> Path:
    return root / "exominer_pipeline" / "run_pipeline.py"


def check_exominer_dependencies(python_executable: str) -> Optional[str]:
    """
    Vérifie que l'interpréteur peut importer les modules lourds utilisés dès le prétraitement
    Exominer (TensorFlow, lightkurve, etc.).

    Returns:
        None si OK, sinon un message d'erreur à afficher à l'utilisateur.
    """
    kw: dict = dict(
        capture_output=True,
        text=True,
        timeout=240,
        encoding="utf-8",
        errors="replace",
    )
    if sys.platform == "win32":
        kw["creationflags"] = subprocess.CREATE_NO_WINDOW  # type: ignore[attr-defined]
    probe = (
        "import tensorflow as tf\n"
        "import lightkurve as lk\n"
        "from pydl.pydlutils import bspline\n"
        "assert tf.__version__\n"
    )
    try:
        proc = subprocess.run(
            [python_executable, "-c", probe],
            **kw,
        )
    except subprocess.TimeoutExpired:
        return (
            "Le test des imports ExoMiner a dépassé le délai (premier chargement TensorFlow souvent long).\n\n"
            f"Interpréteur : {python_executable}"
        )
    except OSError as e:
        return f"Impossible d'exécuter l'interpréteur ExoMiner :\n{python_executable}\n\n{e}"

    if proc.returncode == 0:
        return None

    err = (proc.stderr or proc.stdout or "").strip()
    low = err.lower()
    if "no module named" in low and "tensorflow" in low:
        return (
            "TensorFlow n'est pas installé pour l'interpréteur qui lance ExoMiner.\n\n"
            f"Interpréteur actuel : {python_executable}\n\n"
            "Installez TensorFlow dans cet environnement, par exemple :\n"
            "  pip install tensorflow\n\n"
            "Vous pouvez aussi définir la variable d'environnement NPOAP_EXOMINER_PYTHON "
            "vers un autre python.exe (conda) où les dépendances Exominer sont installées ; "
            "voir le dépôt NASA Exominer (conda_env_exoplnt_dl_*.yml)."
        )
    if "no module named" in low and "lightkurve" in low:
        return (
            "Le paquet lightkurve est requis par le prétraitement des courbes Exominer.\n\n"
            f"Interpréteur : {python_executable}\n\n"
            "Installez-le par exemple avec :\n"
            "  pip install lightkurve\n\n"
            "(Astropy est une dépendance de lightkurve ; pip l'installera si besoin.)"
        )
    if "no module named" in low and "pydl" in low:
        return (
            "Le paquet pydl est requis par le module kepler_spline (lissage spline) d'Exominer.\n\n"
            f"Interpréteur : {python_executable}\n\n"
            "Installez-le avec :\n"
            "  pip install pydl\n"
        )
    return (
        f"Import des dépendances ExoMiner impossible pour :\n{python_executable}\n\n{err[:1600]}"
    )


def check_exominer_pipeline_data(
    exominer_root: Union[str, Path],
    *,
    exominer_model: str = "exominer++_single",
) -> Optional[str]:
    """
    Vérifie la présence des poids (.keras) et des statistiques de normalisation attendus par le pipeline
    lorsqu'on lance ExoMiner hors conteneur (ces fichiers ne sont pas dans le dépôt Git public).

    Returns:
        None si tout est présent, sinon un message d'erreur pour l'utilisateur.
    """
    root = Path(exominer_root).expanduser().resolve()
    model_arg = (exominer_model or "").strip()
    candidates: List[Path] = []

    if model_arg in _EXOMINER_MODEL_RELPATHS:
        candidates.append(root / _EXOMINER_MODEL_RELPATHS[model_arg])
    else:
        mp = Path(model_arg)
        candidates.append(mp if mp.is_absolute() else (root / mp))

    for rel in _EXOMINER_NORM_STATS_RELPATHS:
        candidates.append(root / rel)

    missing = [p for p in candidates if not p.is_file()]
    if not missing:
        return None

    lines = "\n".join(f"  - {p}" for p in missing[:12])
    if len(missing) > 12:
        lines += f"\n  … et {len(missing) - 12} autre(s) fichier(s)."

    script_ps1 = Path(__file__).resolve().parent.parent / "scripts" / "exominer_copy_pipeline_data.ps1"
    script_hint = ""
    if script_ps1.is_file():
        script_hint = (
            "\n\nScript NPOAP (Docker ou Podman ; télécharge l’image puis copie ``data``) :\n"
            f"  {script_ps1}\n"
            f'  Double-clic : scripts\\exominer_copy_pipeline_data.bat\n'
            "  Prérequis : Docker Desktop ou Podman installé et dans le PATH.\n"
        )

    return (
        "Fichiers attendus par le pipeline ExoMiner introuvables sous le clone (run hors Podman) :\n"
        f"{lines}\n\n"
        "Le dépôt GitHub public ne contient pas les poids TensorFlow ni les ``norm_stats`` ; ils sont livrés dans "
        "l’image officielle ``ghcr.io/nasa/exominer`` (voir la doc NASA « Running the pipeline without Podman »)."
        f"{script_hint}\n"
        "Exemple manuel (Podman ou Docker) depuis la racine du clone ExoMiner :\n"
        "  podman create --name em-extract ghcr.io/nasa/exominer:latest\n"
        "  podman cp em-extract:/app/exominer_pipeline/data ./exominer_pipeline/data\n"
        "  podman rm em-extract\n\n"
        "Documentation : https://github.com/nasa/Exominer/blob/main/docs/running-exominer-pipeline.md"
    )


@dataclass
class ExoMinerProbeResult:
    ok: bool
    root: Optional[Path]
    message: str


def probe_exominer_installation(root: Optional[Union[str, Path]] = None) -> ExoMinerProbeResult:
    """
    Vérifie la présence du clone et du script principal du pipeline TESS (sans importer TensorFlow).
    """
    r = resolve_exominer_root(root)
    if r is None:
        return ExoMinerProbeResult(
            False,
            None,
            "NPOAP_EXOMINER_ROOT non défini : indiquez la racine du dépôt ExoMiner "
            f"({EXOMINER_REPO_URL}).",
        )
    if not r.is_dir():
        return ExoMinerProbeResult(False, r, f"Répertoire inexistant : {r}")
    script = pipeline_script_path(r)
    if not script.is_file():
        return ExoMinerProbeResult(
            False,
            r,
            f"Script pipeline introuvable : {script}. Vérifiez que la racine est bien le clone Git.",
        )
    return ExoMinerProbeResult(True, r, f"Pipeline détecté : {script}")


@dataclass
class ExoMinerPipelineResult:
    returncode: int
    output_dir: Path
    stdout: str
    stderr: str
    predictions_csv: Optional[Path]

    @property
    def success(self) -> bool:
        return self.returncode == 0


def _normalize_external_repo(val: Optional[str]) -> Optional[str]:
    """None si telechargement MAST ; sinon chemin expandi."""
    if val is None or str(val).strip().lower() in ("", "null", "none"):
        return None
    return str(Path(val).expanduser().resolve())


def write_tics_table_csv(
    path: Union[str, Path],
    rows: Sequence[Mapping[str, Any]],
) -> Path:
    """
    Écrit le CSV attendu par ExoMiner : colonnes ``tic_id``, ``sector_run`` (ex. ``6-6``, ``1-39``).

    Args:
        path: Chemin du fichier à créer.
        rows: Séquence de dicts avec au moins ``tic_id`` (int) et ``sector_run`` (str).
    """
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas est requis pour write_tics_table_csv ; installez pandas dans l'environnement NPOAP."
        ) from e

    df = pd.DataFrame(list(rows))
    for col in ("tic_id", "sector_run"):
        if col not in df.columns:
            raise ValueError(f"Colonne manquante : {col}")
    df = df[["tic_id", "sector_run"]].copy()
    df.to_csv(path, index=False)
    logger.info("[ExoMiner] Table TIC écrite : %s (%s lignes)", path, len(df))
    return path


def mast_sector_run_for_tic(tic_id: int, *, data_collection_mode: str = "2min") -> str:
    """
    Interroge MAST (astroquery) et déduit un ``sector_run`` ``début-fin`` pour ExoMiner,
    à partir des noms de fichiers des rapports DV SPOC ``dvr.xml`` disponibles pour ce TIC.
    """
    try:
        from astroquery.mast import Observations
    except ImportError as e:
        raise ImportError(
            "astroquery est requis pour résoudre les secteurs TESS depuis MAST. "
            "Installez astroquery dans l'environnement NPOAP."
        ) from e

    # Optionnel : nécessite boto3 ; les requêtes catalogue MAST fonctionnent sans.
    try:
        Observations.enable_cloud_dataset()
    except ImportError:
        logger.info("[ExoMiner] boto3 absent — poursuite sans enable_cloud_dataset().")
    mode = (data_collection_mode or "2min").strip().lower()
    obs_collection = "TESS" if mode == "2min" else "HLSP"

    obs_table = Observations.query_criteria(
        target_name=str(int(tic_id)),
        obs_collection=obs_collection,
    )
    if len(obs_table) == 0:
        raise RuntimeError(f"Aucune observation MAST TESS pour le TIC {tic_id}.")

    products = Observations.get_product_list(obs_table)
    if len(products) == 0:
        raise RuntimeError(f"Aucun produit MAST pour le TIC {tic_id}.")

    span_pat = re.compile(r"-s(\d{4})-s(\d{4})")
    pairs: List[Tuple[int, int]] = []
    for fn in products["productFilename"]:
        if not str(fn).lower().endswith("dvr.xml"):
            continue
        m = span_pat.search(str(fn))
        if m:
            pairs.append((int(m.group(1)), int(m.group(2))))

    if not pairs:
        raise RuntimeError(
            f"Aucun fichier dvr.xml exploitable pour déterminer sector_run (TIC {tic_id}). "
            "Utilisez le mode dossier local avec un CSV tic_id / sector_run."
        )

    best = max(pairs, key=lambda t: (t[1], t[0]))
    return f"{best[0]}-{best[1]}"


def find_predictions_csv(output_dir: Union[str, Path]) -> Optional[Path]:
    """
    Localise le CSV agrégé des prédictions après un run réussi.

    Le script amont écrit ``predictions_<basename(output_dir)>.csv`` et des
    ``ranked_predictions_predictset.csv`` dans chaque ``job_*``.
    """
    out = Path(output_dir).expanduser().resolve()
    # Fichier agrégé principal
    agg = out / f"predictions_{out.name}.csv"
    if agg.is_file():
        return agg
    # Repli : premier ranked_predictions_predictset.csv
    ranked = sorted(out.rglob("ranked_predictions_predictset.csv"))
    if ranked:
        return ranked[0]
    globs = sorted(out.glob("predictions*.csv"))
    return globs[0] if globs else None


def build_run_pipeline_command(
    *,
    exominer_root: Path,
    python_executable: str,
    output_dir: Union[str, Path],
    tic_ids_fp: Optional[Union[str, Path]] = None,
    tic_ids: Optional[str] = None,
    data_collection_mode: str = "2min",
    num_processes: int = 1,
    num_jobs: int = 1,
    download_spoc_data_products: str = "false",
    external_data_repository: Optional[str] = None,
    stellar_parameters_source: str = "ticv8",
    ruwe_source: str = "gaiadr2",
    exominer_model: str = "exominer++_single",
) -> List[str]:
    """Construit la liste d'arguments pour subprocess ([exe, script] + args)."""
    root = Path(exominer_root).resolve()
    script = pipeline_script_path(root)
    out = Path(output_dir).expanduser()

    cmd: List[str] = [
        python_executable,
        str(script),
        "--output_dir",
        str(out),
        "--data_collection_mode",
        data_collection_mode,
        "--num_processes",
        str(num_processes),
        "--num_jobs",
        str(num_jobs),
        "--download_spoc_data_products",
        download_spoc_data_products,
    ]
    # Ne pas passer --external_data_repository en mode MAST : le pipeline NASA attend
    # l'absence d'argument (argparse default None), pas la chaine "null" (chemin invalide).
    ext_resolved = _normalize_external_repo(external_data_repository)
    if ext_resolved is not None:
        cmd.extend(["--external_data_repository", ext_resolved])
    cmd.extend(
        [
            "--stellar_parameters_source",
            stellar_parameters_source,
            "--ruwe_source",
            ruwe_source,
            "--exominer_model",
            exominer_model,
        ]
    )
    if tic_ids_fp is not None:
        cmd.extend(["--tic_ids_fp", str(Path(tic_ids_fp).expanduser().resolve())])
    elif tic_ids:
        cmd.extend(["--tic_ids", tic_ids])
    else:
        raise ValueError("Spécifiez tic_ids_fp ou tic_ids (voir argparse du pipeline NASA).")

    return cmd


def run_exominer_pipeline(
    output_dir: Union[str, Path],
    *,
    tic_ids_fp: Optional[Union[str, Path]] = None,
    tic_ids: Optional[str] = None,
    data_collection_mode: str = "2min",
    num_processes: int = 1,
    num_jobs: int = 1,
    download_spoc_data_products: str = "false",
    external_data_repository: Optional[str] = None,
    stellar_parameters_source: str = "ticv8",
    ruwe_source: str = "gaiadr2",
    exominer_model: str = "exominer++_single",
    exominer_root: Optional[Union[str, Path]] = None,
    python_executable: Optional[Union[str, Path]] = None,
    timeout_sec: Optional[float] = None,
    env: Optional[Mapping[str, str]] = None,
    skip_tf_check: bool = False,
    skip_pipeline_data_check: bool = False,
) -> ExoMinerPipelineResult:
    """
    Lance ``exominer_pipeline/run_pipeline.py`` depuis la racine du clone via subprocess.

    Args:
        output_dir: Répertoire de sortie du run (créé si besoin).
        tic_ids_fp: CSV avec colonnes tic_id, sector_run.
        tic_ids: Liste séparée par virgules (syntaxe NASA), ignorée si tic_ids_fp est fourni.
        data_collection_mode: ``2min`` ou ``ffi``.
        download_spoc_data_products: ``\"true\"`` ou ``\"false\"``.
        external_data_repository: Chemin local vers FITS + DV XML, ou None pour téléchargement MAST (aucun flag NASA).
        exominer_model: ``exominer++_single``, ``exominer++_cviter-mean-ensemble``, ``exominer++_cv-super-mean-ensemble``, ou chemin ``.keras``.
        exominer_root: Surcharge de NPOAP_EXOMINER_ROOT.
        python_executable: Surcharge de NPOAP_EXOMINER_PYTHON / sys.executable.
        timeout_sec: Optionnel ; None = pas de limite (les runs peuvent être longs).
        env: Variables d'environnement additionnelles (fusionnées avec os.environ).
        skip_tf_check: Si True, ne pas vérifier TensorFlow (déjà fait par l'appelant).
        skip_pipeline_data_check: Si True, ne pas vérifier les fichiers .keras / norm_stats (avancé).

    Returns:
        ExoMinerPipelineResult avec chemins et code retour.
    """
    probe = probe_exominer_installation(exominer_root)
    if not probe.ok or probe.root is None:
        raise RuntimeError(probe.message)

    root = probe.root
    py = resolve_exominer_python(python_executable)
    if not skip_tf_check:
        dep_msg = check_exominer_dependencies(py)
        if dep_msg:
            raise RuntimeError(dep_msg)
    if not skip_pipeline_data_check:
        data_msg = check_exominer_pipeline_data(root, exominer_model=exominer_model)
        if data_msg:
            raise RuntimeError(data_msg)
    out_path = Path(output_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    cmd = build_run_pipeline_command(
        exominer_root=root,
        python_executable=py,
        output_dir=out_path,
        tic_ids_fp=tic_ids_fp,
        tic_ids=tic_ids,
        data_collection_mode=data_collection_mode,
        num_processes=num_processes,
        num_jobs=num_jobs,
        download_spoc_data_products=download_spoc_data_products,
        external_data_repository=external_data_repository,
        stellar_parameters_source=stellar_parameters_source,
        ruwe_source=ruwe_source,
        exominer_model=exominer_model,
    )

    logger.info("[ExoMiner] cwd=%s", root)
    logger.info("[ExoMiner] cmd=%s", " ".join(cmd))

    run_env = {**os.environ, **dict(env or {})}
    # run_pipeline.py vit sous exominer_pipeline/ : sans ceci, seul ce dossier est sur sys.path
    # et les imports ``exominer_pipeline``, ``src_preprocessing``, ``query_dv_reports`` échouent.
    root_str = str(root)
    _prev_pp = run_env.get("PYTHONPATH", "").strip()
    run_env["PYTHONPATH"] = (
        root_str if not _prev_pp else f"{root_str}{os.pathsep}{_prev_pp}"
    )
    # Journaux TensorFlow / oneDNN plus sobres (surcharge possible avant le lancement NPOAP)
    run_env.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    run_env.setdefault("TF_ENABLE_ONEDNN_OPTS", "0")
    # lightkurve (oktopus), numpy/pandas FutureWarning dans les deps NASA
    run_env.setdefault("PYTHONWARNINGS", "ignore::UserWarning,ignore::FutureWarning")

    proc = subprocess.run(
        cmd,
        cwd=str(root),
        env=run_env,
        capture_output=True,
        text=True,
        timeout=timeout_sec,
    )

    pred = find_predictions_csv(out_path) if proc.returncode == 0 else None
    if proc.returncode != 0:
        logger.error("[ExoMiner] échec rc=%s\nstderr=%s", proc.returncode, proc.stderr[:4000])

    return ExoMinerPipelineResult(
        returncode=proc.returncode,
        output_dir=out_path,
        stdout=proc.stdout or "",
        stderr=proc.stderr or "",
        predictions_csv=pred,
    )
