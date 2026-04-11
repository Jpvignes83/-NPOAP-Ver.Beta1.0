import json
import os
import logging
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Callable, Tuple
from dataclasses import dataclass
import shutil
import time
import requests
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import subprocess
import numpy as np
import tempfile

logger = logging.getLogger(__name__)


def run_plate_solve_parallel(
    items: list[Path],
    work_one: Callable[[Path], None],
    *,
    max_workers: int,
    progress_callback: Optional[Callable[[int], None]] = None,
    pipeline_control=None,
    log_prefix: str = "Plate-solve",
) -> None:
    """
    Lance work_one(path) sur plusieurs FITS en parallèle (pause/stop respectés au début de chaque tâche).
    La progression est mise à jour après chaque fichier réellement traité (pas les tâches annulées par Stop).
    """
    total = len(items)
    if total == 0:
        return
    mw = max(1, min(max(1, max_workers), total))
    completed = 0
    lock = threading.Lock()

    def work(p: Path) -> None:
        nonlocal completed
        if pipeline_control is not None:
            pipeline_control.wait_if_paused()
            if pipeline_control.should_stop():
                return
        try:
            work_one(p)
        except Exception as e:
            logger.error("%s — échec sur %s : %s", log_prefix, p.name, e)
        with lock:
            completed += 1
            if progress_callback is not None:
                progress_callback(min(100, int(completed / total * 100)))

    with ThreadPoolExecutor(max_workers=mw) as ex:
        futures = [ex.submit(work, p) for p in items]
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception:
                pass


class AstrometrySolverNova:
    def __init__(
        self,
        api_key_file: str | Path | None = None,
        output_dir: str | Path | None = None,
        downsample_factor: int = 1,
        parallel_workers: int = 4,
    ):
        self.api_url = "https://nova.astrometry.net/api"
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.api_key_file = Path(api_key_file) if api_key_file else Path.home() / ".astrometry_api_key"
        self.api_key = self._load_api_key()
        self.downsample_factor = downsample_factor
        self.parallel_workers = max(1, int(parallel_workers))

    # -------------------------------------------------------------------
    # HELPERS
    # -------------------------------------------------------------------
    def _safe_ascii(self, x) -> str:
        return str(x).encode("ascii", "ignore").decode("ascii")

    def _log(self, msg: str):
        print(f"[NOVA] {msg}")

    def _estimate(self, start, phase, total=5):
        elapsed = time.time() - start
        remaining = (elapsed / (phase + 1)) * (total - phase - 1)
        return f"{remaining:.1f}s restantes"

    def _load_api_key(self) -> str:
        if self.api_key_file.exists():
            return self.api_key_file.read_text().strip()
        raise FileNotFoundError(f"Clé API manquante : {self.api_key_file}")

    # -------------------------------------------------------------------
    # LOGIN
    # -------------------------------------------------------------------
    def _login(self, session: requests.Session) -> str:
        login_url = f"{self.api_url}/login"
        data = {"request-json": json.dumps({"apikey": self.api_key})}

        resp = session.post(login_url, data=data, timeout=60)
        try:
            response = resp.json()
        except Exception:
            raise RuntimeError(f"Réponse non-JSON Nova (login) : {self._safe_ascii(resp.text)[:500]}")

        if response.get("status") != "success":
            raise RuntimeError(f"Échec login Astrometry.net : {self._safe_ascii(response)}")

        return response["session"]

    # -------------------------------------------------------------------
    # ANTI-429 : SAFE JSON
    # -------------------------------------------------------------------
    def _safe_get_json(self, session, url, context, sleep=2, retries=10):
        for attempt in range(retries):
            resp = session.get(url)

            if resp.status_code == 429:
                time.sleep(sleep * (attempt + 1))
                continue

            try:
                return resp.json()
            except Exception:
                raise RuntimeError(
                    f"Réponse non-JSON Nova ({context}) : "
                    f"{self._safe_ascii(resp.text)[:500]}"
                )

        raise RuntimeError(f"Abandon après {retries} tentatives sur {context} (429 répétés)")

    # -------------------------------------------------------------------
    # ANTI-429 : SAFE WCS DOWNLOAD
    # -------------------------------------------------------------------
    def _safe_get_wcs(self, session, url, sleep=2, retries=10):
        for attempt in range(retries):
            resp = session.get(url)

            if resp.status_code == 429:
                time.sleep(sleep * (attempt + 1))
                continue

            if resp.status_code == 200:
                return resp.content

            raise RuntimeError(
                f"Échec téléchargement WCS ({resp.status_code}) : "
                f"{self._safe_ascii(resp.text)[:500]}"
            )

        raise RuntimeError(f"Abandon après {retries} tentatives de téléchargement WCS (429 répétés)")

    # -------------------------------------------------------------------
    # MAIN SOLVER
    # -------------------------------------------------------------------
    def solve_file(self, fits_path: Path, progress_callback=None):

        start_time = time.time()
        phase = 0
        self._log("Démarrage solveur NOVA…")

        # ----------------------------------------------------------
        # 1) Downsampling local si FITS > 30 Mo
        # ----------------------------------------------------------
        self._log("Préparation du fichier FITS…")
        fits_to_upload = fits_path

        ds = self.downsample_factor
        if ds > 1:
            
            with fits.open(fits_path, memmap=False) as hdul:
                data = np.array(hdul[0].data, copy=True)
                header = hdul[0].header.copy()

            data_small = data[::ds, ::ds]
            tmp_path = self.output_dir / f"{fits_path.stem}_ds{ds}.fits"
            fits.writeto(tmp_path, data_small, header, overwrite=True)
            fits_to_upload = tmp_path

        # ----------------------------------------------------------
        # 2) Upload
        # ----------------------------------------------------------
        self._log("PHASE 1/5 : Upload…")
        self._log("   Estimation : " + self._estimate(start_time, phase))

        with requests.Session() as session:

            session_key = self._login(session)
            upload_url = f"{self.api_url}/upload"

            with open(fits_to_upload, "rb") as f:
                files = {"file": (fits_to_upload.name, f, "application/octet-stream")}
                data = {"request-json": json.dumps({"session": session_key})}
                r = session.post(upload_url, data=data, files=files, timeout=300)

            try:
                upload_resp = r.json()
            except Exception:
                raise RuntimeError(f"Réponse non-JSON Nova (upload) : {self._safe_ascii(r.text)[:500]}")

            subid = upload_resp.get("subid")
            if not subid:
                raise RuntimeError(f"Échec upload : {self._safe_ascii(upload_resp)}")

        if progress_callback: progress_callback(20)
        phase += 1

        # ----------------------------------------------------------
        # 3) Polling submissions -> job_id
        # ----------------------------------------------------------
        self._log("PHASE 2/5 : Polling des submissions…")
        self._log("   Estimation : " + self._estimate(start_time, phase))

        with requests.Session() as session:
            session_key = self._login(session)

            job_id = None
            while job_id is None:
                url = f"{self.api_url}/submissions/{subid}"
                status = self._safe_get_json(session, url, "submissions")

                jobs = status.get("jobs", [])
                if jobs and jobs[0] is not None:
                    job_id = jobs[0]
                    break

                time.sleep(2)
                if progress_callback: progress_callback(30)

        phase += 1

        # ----------------------------------------------------------
        # 4) Polling job
        # ----------------------------------------------------------
        self._log("PHASE 3/5 : Polling du job…")
        self._log("   Estimation : " + self._estimate(start_time, phase))

        with requests.Session() as session:
            session_key = self._login(session)

            while True:
                url = f"{self.api_url}/jobs/{job_id}"
                job_status = self._safe_get_json(session, url, "jobs")

                if job_status.get("status") == "success":
                    break
                if job_status.get("status") == "failure":
                    raise RuntimeError("Astrométrie échouée")

                time.sleep(5)
                if progress_callback:
                    progress_callback(50)

        phase += 1

        # ----------------------------------------------------------
        # 5) Télécharger WCS
        # ----------------------------------------------------------
        self._log("PHASE 4/5 : Téléchargement WCS…")
        self._log("   Estimation : " + self._estimate(start_time, phase))

        with requests.Session() as session:
            session_key = self._login(session)

            base_url = self.api_url.rstrip("/")
            download_base = base_url[:-4] if base_url.endswith("/api") else base_url
            wcs_url = f"{download_base}/wcs_file/{job_id}"

            wcs_path = self.output_dir / f"{fits_path.stem}.wcs"
            wcs_content = self._safe_get_wcs(session, wcs_url)

            with open(wcs_path, "wb") as f:
                f.write(wcs_content)

        if progress_callback: progress_callback(70)
        phase += 1

        # ----------------------------------------------------------
        # 6) Injection WCS minimale
        # ----------------------------------------------------------
        self._log("PHASE 5/5 : Injection du WCS…")
        self._log("   Estimation : " + self._estimate(start_time, phase))

        w = WCS(str(wcs_path))
        wcs_header = w.to_header()

        with fits.open(fits_path, mode="update", memmap=False) as hdul:
            hdr = hdul[0].header
            for key, val in wcs_header.items():
                try:
                    if isinstance(val, str):
                        val = val.encode("ascii", "ignore").decode()
                    hdr[key] = val
                except Exception:
                    pass
            hdul.flush()

        if progress_callback: progress_callback(100)

        # ----------------------------------------------------------
        # 7) Copier dans dossier science
        # ----------------------------------------------------------
        science_dir = self.output_dir.parent / "science"
        science_dir.mkdir(exist_ok=True)

        solved_path = science_dir / f"{fits_path.stem}-platesolved.fits"
        shutil.copy2(fits_path, solved_path)
        self._log(f"✔ Plate-solved FITS sauvegardé : {solved_path}")

        # ----------------------------------------------------------
        # 8) Nettoyage des fichiers temporaires
        # ----------------------------------------------------------
        if wcs_path.exists():
            os.remove(wcs_path)
            self._log(f"Supprimé : {wcs_path.name}")

        ds_files = list(self.output_dir.glob(f"{fits_path.stem}_ds*.fits"))
        for d in ds_files:
            try:
                os.remove(d)
                self._log(f"Supprimé : {d.name}")
            except:
                pass

        # FIN
        self._log("✔ Astrométrie NOVA terminée.")
        return True

    # -------------------------------------------------------------------
    # Solve full directory
    # -------------------------------------------------------------------
    def solve_directory(self, directory: Path, progress_callback=None, pipeline_control=None):
        """
        Résout les FITS du répertoire en parallèle (jusqu'à ``parallel_workers`` à la fois).
        Progression : une mise à jour par fichier terminé (compatible barre monotone côté UI).
        """
        files = sorted(directory.glob("*.fits"))
        if not files:
            return
        logger.info(
            "NOVA : %d fichier(s), jusqu'à %d résolutions en parallèle.",
            len(files),
            self.parallel_workers,
        )
        run_plate_solve_parallel(
            files,
            lambda p: self.solve_file(p, None),
            max_workers=self.parallel_workers,
            progress_callback=progress_callback,
            pipeline_control=pipeline_control,
            log_prefix="NOVA",
        )

@dataclass
class SolverConfig:
    """Paramètres de résolution astrométrique."""
    scale_low: float = 0.3      # Echelle min (arcsec/pix)
    scale_high: float = 2.0     # Echelle max (arcsec/pix)
    downsample: int = 2         # Réduction de l'image pour accélérer
    timeout: int = 300          # Temps max en secondes (augmenté à 5 min)
    cpulimit: int = 120         # Temps CPU max (augmenté à 2 min)
    omp_threads: int = 4        # Threads OpenMP par solve-field (WSL), voir parallel_workers
    parallel_workers: int = 4   # Nombre d'images traitées en parallèle (LOCAL / NOVA via ImageProcessor)
    retry_on_failure: bool = True  # Réessayer en cas d'échec
    max_retries: int = 1        # Nombre max de tentatives

class AstrometrySolverLocal:
    """
    Solveur utilisant le moteur Astrometry.net installé dans WSL (Ubuntu).
    """

    def __init__(self, science_dir: Path, config: SolverConfig = SolverConfig()):
        self.science_dir = Path(science_dir)
        self.config = config
        self.bash_cmd = "wsl"  # On appelle directement WSL

    def _to_wsl_path(self, windows_path: Path) -> str:
        r"""
        Convertit un chemin Windows (C:\...) en chemin WSL (/mnt/c/...).
        """
        resolved = windows_path.resolve()
        drive = resolved.drive.replace(':', '').lower()
        posix_path = resolved.as_posix().split(':', 1)[1]
        return f"/mnt/{drive}{posix_path}"

    def _update_fits_header(self, fits_path: Path, wcs: WCS, time_keys: dict) -> None: # NOUVEL ARGUMENT
        """Injecte le WCS et les clés de temps dans le fichier FITS final."""
        with fits.open(fits_path, mode="update") as hdul:
            hdr = hdul[0].header
            wcs_header = wcs.to_header()
            
            # 1. Mise à jour du WCS
            hdr.update(wcs_header)
            
            # 2. Rétablissement des clés de temps critiques
            if "JD-UTC" in time_keys:
                hdr["JD-UTC"] = (time_keys["JD-UTC"], "Julian Date (UTC) at mid-exposure")
            if "BJD-TDB" in time_keys:
                hdr["BJD-TDB"] = (time_keys["BJD-TDB"], "Barycentric Julian Date (TDB) mid-exposure")
            
            # 3. Ajout de traçabilité
            hdr.add_history("Plate-solved via Local WSL Astrometry")
            hdr["PLTSOLVD"] = (True, "Solved using local Astrometry.net")
            
            hdul.flush()

    def _validate_fits_file(self, fits_path: Path) -> Tuple[bool, str]:
        """
        Valide un fichier FITS avant de l'envoyer à WSL.
        
        Returns
        -------
        Tuple[bool, str]
            (is_valid, error_message)
        """
        try:
            # 1. Vérifier que le fichier existe et est lisible
            if not fits_path.exists():
                return False, "Fichier introuvable"
            
            # 2. Vérifier la taille du fichier (doit être > 0)
            if fits_path.stat().st_size == 0:
                return False, "Fichier vide"
            
            # 3. Essayer d'ouvrir le fichier FITS
            try:
                with fits.open(fits_path, mode='readonly') as hdul:
                    # 4. Vérifier qu'il y a au moins un HDU
                    if len(hdul) == 0:
                        return False, "Aucun HDU dans le fichier FITS"
                    
                    # 5. Vérifier que le premier HDU contient des données
                    hdu0 = hdul[0]
                    if hdu0.data is None:
                        return False, "HDU primaire sans données"
                    
                    # 6. Vérifier que les données sont valides (pas toutes NaN/inf)
                    data = hdu0.data
                    if data.size == 0:
                        return False, "Données vides"
                    
                    # 7. Vérifier le type de données (doit être numérique)
                    if not np.issubdtype(data.dtype, np.number):
                        return False, f"Type de données non numérique : {data.dtype}"
                    
                    # 8. Vérifier les dimensions (doit être 2D pour une image)
                    if data.ndim != 2:
                        return False, f"Dimensions invalides : {data.ndim}D (attendu 2D)"
                    
                    # 9. Vérifier que les dimensions sont raisonnables
                    h, w = data.shape
                    if h < 10 or w < 10:
                        return False, f"Image trop petite : {w}x{h} pixels"
                    if h > 50000 or w > 50000:
                        return False, f"Image trop grande : {w}x{h} pixels (peut causer des problèmes WSL)"
                    
                    # 10. Vérifier qu'il y a des valeurs finies
                    finite_count = np.isfinite(data).sum()
                    if finite_count == 0:
                        return False, "Aucune valeur finie dans l'image"
                    
                    # 11. Vérifier le format BITPIX (doit être standard)
                    bitpix = hdu0.header.get('BITPIX', None)
                    if bitpix is not None:
                        # BITPIX standard : 8, 16, 32, -32, -64
                        valid_bitpix = [8, 16, 32, -32, -64]
                        if bitpix not in valid_bitpix:
                            logger.warning(f"BITPIX non standard : {bitpix} (peut causer des problèmes)")
                    
            except (TypeError, OSError, ValueError) as e:
                error_msg = str(e)
                if "buffer is too small" in error_msg or "too small for requested array" in error_msg:
                    return False, "Fichier FITS corrompu ou incomplet (buffer trop petit)"
                return False, f"Erreur lecture FITS : {error_msg}"
            
            return True, ""
            
        except Exception as e:
            return False, f"Erreur validation : {str(e)}"
    
    def solve_file(self, fits_path: Path, progress_callback: Optional[Callable[[int], None]] = None) -> bool:
        """
        Exécute la résolution via WSL avec retry automatique.
        """
        fits_path = Path(fits_path).resolve()
        # --- Étape 1: EXTRAIRE LES TEMPS AVANT TOUTE MODIFICATION ---
        try:
            original_hdr = fits.getheader(fits_path)
            time_keys = {}
            if "JD-UTC" in original_hdr:
                time_keys["JD-UTC"] = original_hdr["JD-UTC"]
            if "BJD-TDB" in original_hdr:
                time_keys["BJD-TDB"] = original_hdr["BJD-TDB"]
        except Exception:
            time_keys = {} # Aucune donnée temps à conserver
        # -------------------------------------------------------------
        
        if not fits_path.exists():
            logger.error(f"Fichier introuvable : {fits_path}")
            return False
        
        # --- VALIDATION DU FICHIER FITS AVANT ENVOI À WSL ---
        is_valid, error_msg = self._validate_fits_file(fits_path)
        if not is_valid:
            logger.error(f"❌ Fichier FITS invalide : {fits_path.name}")
            logger.error(f"   Raison : {error_msg}")
            logger.error("   Ce fichier ne peut pas être traité par WSL (risque de Bus error)")
            return False

        self.science_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Solving: {fits_path.name}")
        if progress_callback: progress_callback(10)
        
        # Tentatives avec retry
        max_attempts = (self.config.max_retries + 1) if self.config.retry_on_failure else 1
        
        for attempt in range(max_attempts):
            if attempt > 0:
                logger.info(f"Tentative {attempt + 1}/{max_attempts} pour {fits_path.name}")
                # Ajuster les paramètres pour la retry (downsample plus agressif)
                current_downsample = self.config.downsample * (2 ** attempt)
                current_timeout = int(self.config.timeout * 1.5)  # Timeout plus long
            else:
                current_downsample = self.config.downsample
                current_timeout = self.config.timeout

            with tempfile.TemporaryDirectory() as temp_dir:
                work_dir = Path(temp_dir)
                temp_fits = work_dir / fits_path.name
                shutil.copy2(fits_path, temp_fits)
                
                wsl_input_path = self._to_wsl_path(temp_fits)
                
                # Nom de fichier de solution prévu côté Linux /tmp
                linux_solution_name = f"{temp_fits.stem}.new"
                linux_solution_path = f"/tmp/{linux_solution_name}"
                
                # 1. Création de la commande Linux (OMP_NUM_THREADS : parallélisme interne solve-field)
                nt = max(1, int(self.config.omp_threads))
                linux_cmd_str = (
                    f" unset LD_LIBRARY_PATH; "
                    f"export OMP_NUM_THREADS={nt}; "
                    f"/usr/bin/solve-field --overwrite --no-plots --no-verify "
                    f"--temp-dir /tmp --dir /tmp "
                    f"--scale-units arcsecperpix --scale-low {self.config.scale_low} "
                    f"--scale-high {self.config.scale_high} "
                    f"--downsample {current_downsample} "
                    f"--cpulimit {self.config.cpulimit} "
                    f"'{wsl_input_path}'"
                )
                
                cmd = ["wsl", "bash", "-l", "-c", linux_cmd_str]
                logger.info(
                    f"Exécution solve-field via WSL (timeout={current_timeout}s, "
                    f"downsample={current_downsample}, OMP_NUM_THREADS={nt})…"
                )

                if progress_callback: progress_callback(30)

                try:
                    result = subprocess.run(
                        cmd,
                        check=True,
                        capture_output=True,
                        text=True,
                        timeout=current_timeout
                    )
                except FileNotFoundError:
                    logger.error("❌ WSL introuvable. Installez WSL (wsl --install) ou utilisez l'astrométrie 'Via Astrometry.net (NOVA)'.")
                    return False
                except subprocess.CalledProcessError as e:
                    error_code = e.returncode
                    stderr_text = (e.stderr or "")[:2000]
                    stdout_text = (e.stdout or "")[:2000]
                    logger.error(f"solve-field a échoué (code {error_code}) pour {fits_path.name}")
                    if stdout_text:
                        logger.error(f"STDOUT solve-field:\n{stdout_text}")
                    if stderr_text:
                        logger.error(f"STDERR solve-field:\n{stderr_text}")
                    
                    # Diagnostic du code d'erreur
                    if error_code == 255:
                        logger.warning(f"❌ Échec WSL (Code 255) pour {fits_path.name}")
                        logger.warning("   Code 255 = Erreur système (processus tué ou crash)")
                        
                        # Détection spécifique du Bus error
                        if "Bus error" in stderr_text or "core dumped" in stderr_text:
                            logger.error("   ⚠️ Bus error détecté : problème de compatibilité WSL/astrometry.net")
                            logger.error("   Causes possibles :")
                            logger.error("   - Fichier FITS corrompu ou format non standard")
                            logger.error("   - Problème de mémoire/alignement dans WSL")
                            logger.error("   - Incompatibilité entre version astrometry.net et WSL")
                            logger.error("   💡 Suggestion : essayez l'astrométrie via Astrometry.net (NOVA).")
                        
                        if "an-fitstopnm" in stderr_text:
                            logger.error("   ⚠️ Échec lors de la conversion FITS→PNM")
                            logger.error("   Le fichier FITS pourrait être corrompu ou dans un format non supporté")
                    else:
                        logger.warning(f"❌ Échec WSL (Code {error_code}) pour {fits_path.name}")
                    
                    if attempt < max_attempts - 1:
                        logger.info(f"   → Nouvelle tentative avec downsample={current_downsample*2}")
                        continue  # Réessayer
                    if error_code == 255 and ("Bus error" in stderr_text or "an-fitstopnm" in stderr_text):
                        logger.error("   💡 Conseil : essayez l'astrométrie via Astrometry.net (NOVA).")
                    return False
                except subprocess.TimeoutExpired:
                    logger.warning(f"Timeout ({current_timeout}s) dépassé pour {fits_path.name}")
                    if attempt < max_attempts - 1:
                        logger.info(f"   → Nouvelle tentative avec timeout={int(current_timeout*1.5)}s")
                        continue  # Réessayer
                    else:
                        logger.error(f"Timeout ({current_timeout}s) dépassé après {max_attempts} tentatives.")
                        return False

                if progress_callback: progress_callback(80)

                # --- PHASE 2 : Récupération et Copie du Fichier de Solution ---
                
                # Chemin d'arrivée souhaité (Windows)
                wcs_file = temp_fits.with_suffix(".new")
                
                # Commande pour copier le fichier de solution de /tmp Linux vers le dossier temp Windows
                copy_cmd = ["wsl", "cp", linux_solution_path, self._to_wsl_path(wcs_file)]

                try:
                    # Exécution de la copie
                    subprocess.run(copy_cmd, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError:
                    logger.warning("Échec : Le solveur n'a pas créé de fichier de solution valide.")
                    if 'result' in locals():
                         # Affiche la dernière partie du log pour le diagnostic
                        log_end = result.stdout[-300:] 
                        logger.warning(f"Log solveur : {log_end}")
                    if attempt < max_attempts - 1:
                        continue  # Réessayer
                    return False
                    
                # 5. Lecture du WCS et validation (le fichier est maintenant sur le disque Windows)
                if not wcs_file.exists():
                    logger.warning("Fichier de solution non trouvé après la copie.")
                    if attempt < max_attempts - 1:
                        continue  # Réessayer
                    return False
                    
                try:
                    with fits.open(wcs_file) as hdul:
                        # Si la solution est valide, elle doit avoir les champs WCS
                        if "CRVAL1" not in hdul[0].header:
                            logger.warning("Fichier .new trouvé, mais header WCS manquant (échec validation astrométrique).")
                            if attempt < max_attempts - 1:
                                continue  # Réessayer
                            return False

                        wcs_solution = WCS(hdul[0].header)
                except Exception as e:
                    logger.error(f"Impossible de lire le WCS généré : {e}")
                    if attempt < max_attempts - 1:
                        continue  # Réessayer
                    return False
                    
                # --- PHASE 3 : Création du fichier Science ---
                final_name = f"{fits_path.stem}_solved{fits_path.suffix}"
                science_file = self.science_dir / final_name
                shutil.copy2(fits_path, science_file)
                self._update_fits_header(science_file, wcs_solution, time_keys=time_keys)
                
                logger.info(f"Succès ! -> {science_file.name}")
                if progress_callback: progress_callback(100)
                
                # ──────── NETTOYAGE EXPLICITE DES FICHIERS ──────────────────────
                try:
                    # wcs_file est le chemin Windows du fichier de solution copié
                    if wcs_file.exists():
                        os.remove(wcs_file)
                    # temp_fits est la copie FITS de travail (input)
                    if temp_fits.exists():
                        os.remove(temp_fits)
                except Exception as e:
                    logger.debug(f"Erreur lors du nettoyage: {e}")
                
                return True  # Succès, sortir de la boucle
        
        # Si on arrive ici, toutes les tentatives ont échoué
        logger.error(f"Échec après {max_attempts} tentatives pour {fits_path.name}")
        return False

