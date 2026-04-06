import json
import os
import logging
from typing import Optional, Callable, Tuple, List, Dict
from dataclasses import dataclass, asdict
import shutil
import time
import requests
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import subprocess
import numpy as np
import queue
import threading
import tempfile
import gzip
from utils.progress_manager import ProgressManager
import astropy.units as u
from astropy.coordinates import SkyCoord, solar_system_ephemeris
from astropy.table import Table
from astropy.wcs.utils import fit_wcs_from_points
from scipy.spatial import cKDTree
from scipy.optimize import least_squares
from photutils.detection import DAOStarFinder
from photutils.centroids import centroid_2dg
import pandas as pd

logger = logging.getLogger(__name__)

class AstrometrySolverNova:
    def __init__(self, api_key_file: str | Path | None = None, output_dir: str | Path | None = None, downsample_factor: int = 1):
        self.api_url = "https://nova.astrometry.net/api"
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.api_key_file = Path(api_key_file) if api_key_file else Path.home() / ".astrometry_api_key"
        self.api_key = self._load_api_key()
        self.downsample_factor = downsample_factor

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
    def solve_directory(self, directory: Path, progress_callback=None):
        files = sorted(directory.glob("*.fits"))
        total = len(files)
        for i, f in enumerate(files):
            cb = (
                lambda p: progress_callback(int((i + p / 100.0) / total * 100))
                if progress_callback else None
            )
            self.solve_file(f, cb)

logger = logging.getLogger(__name__)

@dataclass
class SolverConfig:
    """Paramètres de résolution astrométrique."""
    scale_low: float = 0.3      # Echelle min (arcsec/pix)
    scale_high: float = 2.0     # Echelle max (arcsec/pix)
    downsample: int = 2         # Réduction de l'image pour accélérer
    timeout: int = 300          # Temps max en secondes (augmenté à 5 min)
    cpulimit: int = 120         # Temps CPU max (augmenté à 2 min)
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
                
                # 1. Création de la commande Linux
                linux_cmd_str = (
                    f" unset LD_LIBRARY_PATH; /usr/bin/solve-field --overwrite --no-plots --no-verify " 
                    f"--temp-dir /tmp --dir /tmp " # Force l'usage de /tmp Linux pour tous les fichiers intermédiaires
                    f"--scale-units arcsecperpix --scale-low {self.config.scale_low} "
                    f"--scale-high {self.config.scale_high} "
                    f"--downsample {current_downsample} "
                    f"--cpulimit {self.config.cpulimit} "
                    f"'{wsl_input_path}'" # Chemin d'entrée Windows
                )
                
                cmd = ["wsl", "bash", "-l", "-c", linux_cmd_str]
                logger.info(f"Exécution solve-field via WSL (timeout={current_timeout}s, downsample={current_downsample})…")

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
                            logger.error("   💡 Suggestion : Utilisez la méthode 'LocalG' (catalogues locaux) à la place")
                        
                        if "an-fitstopnm" in stderr_text:
                            logger.error("   ⚠️ Échec lors de la conversion FITS→PNM")
                            logger.error("   Le fichier FITS pourrait être corrompu ou dans un format non supporté")
                    else:
                        logger.warning(f"❌ Échec WSL (Code {error_code}) pour {fits_path.name}")
                    
                    if attempt < max_attempts - 1:
                        logger.info(f"   → Nouvelle tentative avec downsample={current_downsample*2}")
                        continue  # Réessayer
                    if error_code == 255 and ("Bus error" in stderr_text or "an-fitstopnm" in stderr_text):
                        logger.error("   💡 Conseil : Essayez la méthode 'LocalG' (catalogues locaux) ou NOVA.")
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
            
            
class AstrometrySolverLocalG:
    """
    Solveur utilisant les catalogues Gaia DR3 locaux (CSV/CSV.GZIP ou format 1476-6).
    Analogue à ASTAP mais utilise directement la magnitude G de Gaia (pas de transformation vers V).
    """
    
    def __init__(self, catalog_dir: Path, science_dir: Path):
        """
        Initialise le solveur avec le répertoire contenant les catalogues.
        
        Parameters
        ----------
        catalog_dir : Path
            Répertoire contenant les fichiers de catalogues (CSV, CSV.GZIP, ou .1476)
        science_dir : Path
            Répertoire de sortie pour les images résolues
        """
        self.catalog_dir = Path(catalog_dir)
        self.science_dir = Path(science_dir)
        
        if not self.catalog_dir.exists():
            raise ValueError(f"Répertoire de catalogues introuvable : {catalog_dir}")
        
        self.science_dir.mkdir(parents=True, exist_ok=True)
        
        # Index des catalogues pour recherche rapide
        self._catalog_index: Optional[Dict] = None
        self._index_file = self.catalog_dir / ".catalog_index.json"
        self._index_lock = threading.Lock()
        
        logger.info(f"AstrometrySolverLocalG initialisé avec catalogues : {catalog_dir}")
    
    def _parse_catalog_filename(self, file_path: Path) -> Optional[Dict]:
        """
        Parse le nom d'un fichier de catalogue pour extraire les plages RA/DEC couvertes.
        
        Formats supportés :
        - gaia_dr3_nord_ra00h-05h_mag18.0.1476 (heures seulement)
        - gaia_dr3_sud_ra00h00m-00h20m_mag18.0.1476 (heures et minutes)
        - gaia_1476_gaia_dr3_sud_ra00h00m-00h20m_mag18.0.1476 (avec préfixe)
        - gaia_dr3_ra00h-05h_mag18.0.1476 (hémisphère non spécifié)
        
        Returns
        -------
        dict avec keys: 'ra_min', 'ra_max', 'dec_min', 'dec_max', 'hemisphere', 'file_path'
        ou None si le parsing échoue
        """
        name = file_path.stem.replace('.csv', '').replace('.gz', '')
        
        try:
            # Chercher le pattern RA dans le nom
            if '_ra' in name:
                ra_part = name.split('_ra')[1].split('_')[0]
                # Formats supportés :
                # - 00h-05h (heures seulement)
                # - 00h00m-00h20m (heures et minutes)
                # - ra00h-05h (avec préfixe "ra")
                # - ra00h00m-00h20m (avec préfixe "ra" et minutes)
                
                # Nettoyer le préfixe "ra" si présent
                if ra_part.startswith('ra'):
                    ra_part = ra_part[2:]
                
                # Détecter le séparateur : chercher "-" entre deux parties avec "h"
                if '-' in ra_part:
                    # Trouver la position du "-" qui sépare les deux parties
                    dash_pos = ra_part.find('-')
                    ra_min_str = ra_part[:dash_pos]
                    ra_max_str = ra_part[dash_pos+1:]
                    
                    # Parser ra_min
                    if 'm' in ra_min_str:
                        # Format avec minutes : 00h00m
                        parts = ra_min_str.split('h')
                        ra_min_h = int(parts[0]) if parts[0] else 0
                        ra_min_m = int(parts[1].replace('m', '')) if len(parts) > 1 and parts[1] else 0
                        ra_min = (ra_min_h + ra_min_m / 60.0) * 15.0
                    else:
                        # Format simple : 00h
                        ra_min_h = int(ra_min_str) if ra_min_str else 0
                        ra_min = ra_min_h * 15.0
                    
                    # Parser ra_max
                    if 'm' in ra_max_str:
                        # Format avec minutes : 00h20m
                        parts = ra_max_str.split('h')
                        ra_max_h = int(parts[0]) if parts[0] else 0
                        ra_max_m = int(parts[1].replace('m', '')) if len(parts) > 1 and parts[1] else 0
                        ra_max = (ra_max_h + ra_max_m / 60.0) * 15.0
                    else:
                        # Format simple : 05h
                        ra_max_h = int(ra_max_str.split('h')[0]) if ra_max_str else 0
                        ra_max = ra_max_h * 15.0
                    
                    # Déterminer DEC selon l'hémisphère
                    if 'nord' in name or 'north' in name:
                        dec_min = -5.0
                        dec_max = 90.0
                        hemisphere = 'north'
                    elif 'sud' in name or 'south' in name:
                        dec_min = -90.0
                        dec_max = 5.0
                        hemisphere = 'south'
                    else:
                        # Par défaut, supposer couverture complète en DEC
                        dec_min = -90.0
                        dec_max = 90.0
                        hemisphere = 'both'
                    
                    return {
                        'ra_min': ra_min,
                        'ra_max': ra_max,
                        'dec_min': dec_min,
                        'dec_max': dec_max,
                        'hemisphere': hemisphere,
                        'file_path': str(file_path),
                        'file_name': file_path.name
                    }
            else:
                # Si le format n'est pas reconnu, essayer de lire le fichier pour déterminer la couverture
                # Pour l'instant, on retourne None et on log un avertissement
                logger.debug(f"Format de nom non reconnu pour {file_path.name} (attendu: *_ra*_*h-*h_*)")
        except Exception as e:
            logger.warning(f"Erreur parsing nom fichier {file_path.name}: {e}")
        
        return None
    
    def _build_catalog_index(self, force_rebuild: bool = False) -> Dict:
        """
        Construit un index de tous les fichiers de catalogue avec leurs plages RA/DEC.
        
        Parameters
        ----------
        force_rebuild : bool
            Si True, reconstruit l'index même s'il existe déjà
        
        Returns
        -------
        dict
            Index des catalogues avec structure :
            {
                'files': [
                    {
                        'ra_min': float,
                        'ra_max': float,
                        'dec_min': float,
                        'dec_max': float,
                        'hemisphere': str,
                        'file_path': str,
                        'file_name': str
                    },
                    ...
                ],
                'last_updated': str (timestamp)
            }
        """
        with self._index_lock:
            # Vérifier si l'index existe et est à jour
            should_rebuild = force_rebuild
            if not should_rebuild and self._index_file.exists():
                try:
                    with open(self._index_file, 'r', encoding='utf-8') as f:
                        cached_index = json.load(f)
                    
                    # Si l'index en cache est vide, forcer la reconstruction
                    cached_files_count = len(cached_index.get('files', []))
                    if cached_files_count == 0:
                        logger.info("Index en cache vide, reconstruction nécessaire...")
                        should_rebuild = True
                    else:
                        # Vérifier si les fichiers existent toujours
                        all_exist = True
                        for entry in cached_index.get('files', []):
                            file_path = Path(entry['file_path'])
                            if not file_path.exists():
                                all_exist = False
                                break
                        
                        if all_exist:
                            # Vérifier s'il y a de nouveaux fichiers .1476 non indexés
                            current_files = set(self.catalog_dir.glob("*.1476"))
                            indexed_files = {Path(entry['file_path']) for entry in cached_index.get('files', [])}
                            if current_files != indexed_files:
                                logger.info(f"Nouveaux fichiers .1476 détectés ({len(current_files)} vs {len(indexed_files)} indexés), reconstruction...")
                                should_rebuild = True
                            else:
                                logger.debug(f"Index de catalogues chargé depuis {self._index_file} ({cached_files_count} fichiers)")
                                self._catalog_index = cached_index
                                return cached_index
                except Exception as e:
                    logger.warning(f"Erreur chargement index : {e}, reconstruction...")
                    should_rebuild = True
            
            # Construire l'index si nécessaire
            if should_rebuild:
                # Construire l'index
                logger.info("Construction de l'index des catalogues...")
                catalog_files = list(self.catalog_dir.glob("*.1476"))
                logger.info(f"Trouvé {len(catalog_files)} fichier(s) .1476 dans {self.catalog_dir}")
                
                if len(catalog_files) == 0:
                    logger.warning(f"Aucun fichier .1476 trouvé dans {self.catalog_dir}")
                    # Lister tous les fichiers pour diagnostic
                    all_files = list(self.catalog_dir.glob("*"))
                    logger.debug(f"Fichiers présents dans le répertoire : {[f.name for f in all_files[:10]]}")
                
                index = {
                    'files': [],
                    'last_updated': time.strftime('%Y-%m-%d %H:%M:%S')
                }
                
                parsed_count = 0
                for cat_file in catalog_files:
                    parsed = self._parse_catalog_filename(cat_file)
                    if parsed:
                        index['files'].append(parsed)
                        parsed_count += 1
                    else:
                        logger.warning(f"Impossible de parser le nom du fichier : {cat_file.name}")
                
                if parsed_count < len(catalog_files):
                    logger.warning(f"Seulement {parsed_count}/{len(catalog_files)} fichier(s) .1476 ont pu être indexés")
                
                # Sauvegarder l'index
                try:
                    with open(self._index_file, 'w', encoding='utf-8') as f:
                        json.dump(index, f, indent=2, ensure_ascii=False)
                    logger.info(f"Index de {len(index['files'])} catalogues sauvegardé dans {self._index_file}")
                except Exception as e:
                    logger.warning(f"Impossible de sauvegarder l'index : {e}")
                
                self._catalog_index = index
                return index
            else:
                # Retourner l'index en cache
                return self._catalog_index if self._catalog_index else {}
    
    def _find_catalog_files(self, center_coord: SkyCoord, search_radius_deg: float, max_files: Optional[int] = None) -> List[Path]:
        """
        Trouve les fichiers de catalogues .1476 couvrant la région du ciel.
        Utilise un index pour une recherche rapide.
        
        Parameters
        ----------
        center_coord : SkyCoord
            Centre de la région à couvrir
        search_radius_deg : float
            Rayon de recherche en degrés
        max_files : int, optional
            Nombre maximum de fichiers à retourner (pour nearby-solve, utiliser 1-2)
        
        Returns
        -------
        List[Path]
            Liste des fichiers de catalogue couvrant la région, triés par proximité
        """
        # Charger/construire l'index
        if self._catalog_index is None:
            self._build_catalog_index()
        
        center_ra = center_coord.ra.deg
        center_dec = center_coord.dec.deg
        
        # Calculer les limites de recherche
        ra_min_search = center_ra - search_radius_deg
        ra_max_search = center_ra + search_radius_deg
        dec_min_search = center_dec - search_radius_deg
        dec_max_search = center_dec + search_radius_deg
        
        # Rechercher dans l'index - calculer d'abord la distance pour filtrer
        candidate_files = []
        
        logger.debug(f"Recherche catalogues : RA={center_ra:.2f}°, DEC={center_dec:.2f}°, rayon={search_radius_deg:.2f}°")
        logger.debug(f"Limites recherche : RA [{ra_min_search:.2f}°, {ra_max_search:.2f}°], DEC [{dec_min_search:.2f}°, {dec_max_search:.2f}°]")
        
        for entry in self._catalog_index.get('files', []):
            file_ra_min = entry['ra_min']
            file_ra_max = entry['ra_max']
            file_dec_min = entry['dec_min']
            file_dec_max = entry['dec_max']
            
            # Calculer la distance au centre AVANT de vérifier le chevauchement
            # pour éliminer rapidement les fichiers trop éloignés
            file_ra_center = (file_ra_min + file_ra_max) / 2.0
            if file_ra_max < file_ra_min:  # Cas 360°
                file_ra_center = 180.0  # Approximatif
            
            file_dec_center = (file_dec_min + file_dec_max) / 2.0
            
            # Distance angulaire approximative (en degrés)
            # Utiliser une distance qui tient compte de la latitude (DEC)
            # Pour les petites distances, on peut approximer avec une distance euclidienne
            # mais en pondérant RA par cos(DEC) pour tenir compte de la réduction de l'échelle en RA
            ra_diff = abs(center_ra - file_ra_center)
            if ra_diff > 180:
                ra_diff = 360 - ra_diff
            dec_diff = abs(center_dec - file_dec_center)
            # Distance pondérée : RA est réduit par cos(DEC) pour tenir compte de la sphère
            cos_dec = np.cos(np.radians(center_dec))
            distance = np.sqrt((ra_diff * cos_dec)**2 + dec_diff**2)
            
            # Filtrer les fichiers trop éloignés (distance > 5x le rayon de recherche)
            # Augmenter le seuil car la distance est maintenant plus précise
            if distance > search_radius_deg * 5.0:
                continue
            
            # Vérifier chevauchement RA (gérer le cas 360°)
            if file_ra_max < file_ra_min:  # Cas où le fichier couvre 360°
                ra_overlap = (ra_min_search < file_ra_max) or (ra_max_search > file_ra_min) or (ra_min_search < 0) or (ra_max_search > 360)
            else:
                ra_overlap = (ra_min_search < file_ra_max) and (ra_max_search > file_ra_min)
            
            # Vérifier chevauchement DEC
            dec_overlap = (dec_min_search < file_dec_max) and (dec_max_search > file_dec_min)
            
            if ra_overlap and dec_overlap:
                file_path = Path(entry['file_path'])
                if file_path.exists():
                    candidate_files.append((file_path, distance))
                    logger.debug(f"Fichier candidat: {file_path.name} (RA [{file_ra_min:.1f}°, {file_ra_max:.1f}°], DEC [{file_dec_min:.1f}°, {file_dec_max:.1f}°], distance={distance:.2f}°)")
        
        # Trier par distance (les plus proches en premier)
        candidate_files.sort(key=lambda x: x[1])
        
        # Limiter le nombre de fichiers si spécifié (pour nearby-solve)
        if max_files is not None:
            candidate_files = candidate_files[:max_files]
        
        catalog_files = [f for f, _ in candidate_files]
        
        if len(catalog_files) > 0:
            max_dist = max(d for _, d in candidate_files) if candidate_files else 0
            logger.debug(f"Trouvé {len(catalog_files)} fichier(s) de catalogue pertinent(s) (distance max: {max_dist:.2f}°)")
        else:
            logger.warning(f"Aucun fichier trouvé dans l'index ({len(self._catalog_index.get('files', []))} fichiers indexés)")
            # Log les premiers fichiers de l'index pour diagnostic
            if len(self._catalog_index.get('files', [])) > 0:
                logger.debug(f"Premiers fichiers indexés:")
                for i, entry in enumerate(self._catalog_index.get('files', [])[:5]):
                    logger.debug(f"  {i+1}. {entry.get('file_name', 'unknown')}: RA [{entry.get('ra_min', 0):.1f}°, {entry.get('ra_max', 0):.1f}°], DEC [{entry.get('dec_min', 0):.1f}°, {entry.get('dec_max', 0):.1f}°]")
        
        return catalog_files
    
    def _file_covers_region(self, file_path: Path, center_coord: SkyCoord, radius_deg: float) -> bool:
        """Vérifie si un fichier de catalogue couvre la région spécifiée."""
        name = file_path.stem.replace('.csv', '').replace('.gz', '')
        
        # Parser le nom : gaia_dr3_{hem}_ra{min}h-{max}h_mag{mag}[_suffix]
        try:
            if '_ra' in name:
                ra_part = name.split('_ra')[1].split('_')[0]
                # Formats supportés :
                # - 00h-05h (heures seulement)
                # - 00h00m-00h20m (heures et minutes)
                # - ra00h-05h (avec préfixe "ra")
                # - ra00h00m-00h20m (avec préfixe "ra" et minutes)
                
                # Nettoyer le préfixe "ra" si présent
                if ra_part.startswith('ra'):
                    ra_part = ra_part[2:]
                
                # Détecter le séparateur : chercher "-" entre deux parties avec "h"
                if '-' in ra_part:
                    # Trouver la position du "-" qui sépare les deux parties
                    dash_pos = ra_part.find('-')
                    ra_min_str = ra_part[:dash_pos]
                    ra_max_str = ra_part[dash_pos+1:]
                    
                    # Parser ra_min
                    if 'm' in ra_min_str:
                        # Format avec minutes : 00h00m
                        parts = ra_min_str.split('h')
                        ra_min_h = int(parts[0]) if parts[0] else 0
                        ra_min_m = int(parts[1].replace('m', '')) if len(parts) > 1 and parts[1] else 0
                        file_ra_min = (ra_min_h + ra_min_m / 60.0) * 15.0
                    else:
                        # Format simple : 00h
                        ra_min_h = int(ra_min_str) if ra_min_str else 0
                        file_ra_min = ra_min_h * 15.0
                    
                    # Parser ra_max
                    if 'm' in ra_max_str:
                        # Format avec minutes : 00h20m
                        parts = ra_max_str.split('h')
                        ra_max_h = int(parts[0]) if parts[0] else 0
                        ra_max_m = int(parts[1].replace('m', '')) if len(parts) > 1 and parts[1] else 0
                        file_ra_max = (ra_max_h + ra_max_m / 60.0) * 15.0
                    else:
                        # Format simple : 05h
                        ra_max_h = int(ra_max_str.split('h')[0]) if ra_max_str else 0
                        file_ra_max = ra_max_h * 15.0
                    
                    # Vérifier chevauchement RA
                    center_ra = center_coord.ra.deg
                    if file_ra_max < file_ra_min:  # Cas 360°
                        ra_overlap = (center_ra - radius_deg < file_ra_max) or (center_ra + radius_deg > file_ra_min)
                    else:
                        ra_overlap = (center_ra - radius_deg < file_ra_max) and (center_ra + radius_deg > file_ra_min)
                    
                    # Vérifier DEC
                    center_dec = center_coord.dec.deg
                    if 'nord' in name or 'north' in name:
                        dec_overlap = (center_dec - radius_deg < 90) and (center_dec + radius_deg > -5)
                    elif 'sud' in name or 'south' in name:
                        dec_overlap = (center_dec - radius_deg < 5) and (center_dec + radius_deg > -90)
                    else:
                        dec_overlap = True
                    
                    return ra_overlap and dec_overlap
        except Exception:
            pass
        
        return True  # Si parsing échoue, inclure par sécurité
    
    def _read_1476_file(self, file_path: Path) -> Optional[Table]:
        """
        Lit un fichier de catalogue au format binaire 1476-6 (ASTAP).
        
        Format 1476-6 (ASTAP) :
        Le format peut varier, mais généralement :
        - En-tête optionnel avec métadonnées
        - Enregistrements d'étoiles : chaque étoile = 6 bytes :
          * RA : 2 bytes (unsigned short, 0-65535 représente 0-360 degrés)
          * DEC : 2 bytes (signed short, -32768 à 32767 représente -90 à +90 degrés)
          * Magnitude : 1 byte (0-255, magnitude * 10, donc 0.0 à 25.5)
          * Flags : 1 byte (réservé ou flags additionnels)
        
        Note: Cette implémentation tente de détecter automatiquement la structure.
        """
        try:
            file_size = file_path.stat().st_size
            
            with open(file_path, 'rb') as f:
                # Essayer différentes structures d'en-tête
                # Structure 1: En-tête de 16 bytes avec nombre d'étoiles
                header = f.read(16)
                if len(header) < 16:
                    logger.error(f"Fichier trop petit : {file_path.name}")
                    return None
                
                # Essayer de détecter le format
                header_size = 16  # Par défaut
                try:
                    n_stars_header_little = int.from_bytes(header[8:12], byteorder='little', signed=False)
                    n_stars_header_big = int.from_bytes(header[8:12], byteorder='big', signed=False)
                    n_stars_calculated = (file_size - 16) // 6
                    
                    # Vérifier si le nombre d'étoiles dans l'en-tête est cohérent
                    if (n_stars_header_little > 0 and n_stars_header_little < 100000000 and 
                        (file_size - 16) % 6 == 0 and 
                        abs(n_stars_header_little - n_stars_calculated) < 10):
                        n_stars = n_stars_header_little
                        header_size = 16
                        logger.debug(f"Format détecté (little-endian) : en-tête 16 bytes, n_stars={n_stars}")
                    elif (n_stars_header_big > 0 and n_stars_header_big < 100000000 and 
                          (file_size - 16) % 6 == 0 and 
                          abs(n_stars_header_big - n_stars_calculated) < 10):
                        n_stars = n_stars_header_big
                        header_size = 16
                        logger.debug(f"Format détecté (big-endian) : en-tête 16 bytes, n_stars={n_stars}")
                    elif file_size % 6 == 0:
                        # Pas d'en-tête, directement les données
                        f.seek(0)
                        n_stars = file_size // 6
                        header_size = 0
                        logger.debug(f"Format détecté : pas d'en-tête, n_stars={n_stars}")
                    else:
                        # Format inconnu, essayer avec en-tête de 16 bytes
                        n_stars = n_stars_calculated
                        header_size = 16
                        logger.warning(f"Format non standard, tentative avec en-tête 16 bytes, n_stars={n_stars}")
                except Exception as e:
                    # En cas d'erreur, essayer sans en-tête
                    if file_size % 6 == 0:
                        f.seek(0)
                        n_stars = file_size // 6
                        header_size = 0
                        logger.debug(f"Format détecté (fallback) : pas d'en-tête, n_stars={n_stars}")
                    else:
                        n_stars = (file_size - 16) // 6
                        header_size = 16
                        logger.warning(f"Format non standard (fallback), tentative avec en-tête 16 bytes, n_stars={n_stars}")
                
                # Lire les enregistrements d'étoiles (6 bytes chacun)
                stars = []
                byteorder_options = ['little', 'big']
                
                for byteorder in byteorder_options:
                    try:
                        f.seek(header_size)  # Position de départ après l'en-tête
                        stars = []
                        
                        for i in range(min(n_stars, 10000000)):  # Limite de sécurité
                            record = f.read(6)
                            if len(record) < 6:
                                break
                            
                            # Parser l'enregistrement
                            ra_raw = int.from_bytes(record[0:2], byteorder=byteorder, signed=False)
                            dec_raw = int.from_bytes(record[2:4], byteorder=byteorder, signed=True)
                            mag_raw = record[4]
                            
                            # Convertir en coordonnées réelles
                            ra_deg = (ra_raw / 65536.0) * 360.0
                            dec_deg = (dec_raw / 32768.0) * 90.0
                            
                            # Correction pour les catalogues "sud" : si le nom du fichier indique "sud" mais que DEC est positif,
                            # inverser le signe (certains formats 1476-6 stockent les DEC sud comme valeurs positives)
                            if ('sud' in file_path.name.lower() or 'south' in file_path.name.lower()) and dec_deg > 0:
                                dec_deg = -dec_deg
                            
                            mag = mag_raw / 10.0
                            
                            # Validation des valeurs
                            if 0 <= ra_deg <= 360 and -90 <= dec_deg <= 90 and 0 <= mag <= 30:
                                stars.append({
                                    'RA_ICRS': ra_deg,
                                    'DE_ICRS': dec_deg,
                                    'Gmag': mag
                                })
                            # Si valeurs invalides, continuer quand même (peut être un problème de byteorder)
                        
                        if len(stars) > 100:  # Si on a trouvé assez d'étoiles valides
                            logger.debug(f"Format détecté avec byteorder={byteorder}, {len(stars)} étoiles valides")
                            break
                    except Exception as e:
                        logger.debug(f"Erreur avec byteorder={byteorder}: {e}")
                        continue
                
                if not stars:
                    logger.error(f"Aucune étoile valide trouvée dans {file_path.name}")
                    return None
                
                catalog_table = Table(rows=stars)
                logger.info(f"Fichier {file_path.name} chargé : {len(catalog_table)} étoiles")
                return catalog_table
                
        except Exception as e:
            logger.error(f"Erreur lecture fichier 1476-6 {file_path.name}: {e}", exc_info=True)
            return None
    
    def _load_catalog_stars(self, center_coord: SkyCoord, search_radius_deg: float, nearby_mode: bool = False, min_stars: int = 100) -> Optional[Table]:
        """
        Charge les étoiles du catalogue .1476 dans la région spécifiée.
        
        Parameters
        ----------
        center_coord : SkyCoord
            Centre de la région
        search_radius_deg : float
            Rayon de recherche en degrés
        nearby_mode : bool
            Si True, charge d'abord un seul fichier, puis élargit si nécessaire
        min_stars : int
            Nombre minimum d'étoiles requis (pour nearby_mode)
        
        Returns
        -------
        Optional[Table]
            Table des étoiles chargées, ou None si aucune trouvée
        """
        if nearby_mode:
            # Nearby-solve : charger d'abord un seul fichier
            catalog_files = self._find_catalog_files(center_coord, search_radius_deg, max_files=1)
            
            if not catalog_files:
                logger.warning(f"Aucun fichier de catalogue .1476 trouvé pour nearby-solve")
                return None
            
            # Essayer d'abord avec un seul fichier
            all_stars = []
            for cat_file in catalog_files:
                catalog_table = self._read_1476_file(cat_file)
                if catalog_table is None or len(catalog_table) == 0:
                    continue
                
                # Filtrer par région
                try:
                    # Log des statistiques avant filtrage
                    if len(catalog_table) > 0:
                        ra_min = np.min(catalog_table['RA_ICRS'])
                        ra_max = np.max(catalog_table['RA_ICRS'])
                        dec_min = np.min(catalog_table['DE_ICRS'])
                        dec_max = np.max(catalog_table['DE_ICRS'])
                        logger.info(f"Avant filtrage {cat_file.name}: {len(catalog_table)} étoiles, RA [{ra_min:.2f}°, {ra_max:.2f}°], DEC [{dec_min:.2f}°, {dec_max:.2f}°]")
                        logger.info(f"Recherche: centre RA={center_coord.ra.deg:.2f}°, DEC={center_coord.dec.deg:.2f}°, rayon={search_radius_deg:.2f}°")
                    
                    catalog_coords = SkyCoord(
                        ra=catalog_table['RA_ICRS'] * u.deg,
                        dec=catalog_table['DE_ICRS'] * u.deg
                    )
                    separations = center_coord.separation(catalog_coords).deg
                    mask = separations <= search_radius_deg
                    n_matches = np.sum(mask)
                    
                    logger.info(f"Après filtrage {cat_file.name}: {n_matches} étoiles dans le rayon de {search_radius_deg:.2f}°")
                    
                    if n_matches > 0:
                        filtered_table = catalog_table[mask]
                        for row in filtered_table:
                            all_stars.append({
                                'RA_ICRS': row['RA_ICRS'],
                                'DE_ICRS': row['DE_ICRS'],
                                'Gmag': row['Gmag']
                            })
                    else:
                        # Log quelques exemples d'étoiles pour diagnostic
                        if len(catalog_table) > 0:
                            sample_indices = np.linspace(0, len(catalog_table)-1, min(5, len(catalog_table)), dtype=int)
                            logger.warning(f"Aucune étoile trouvée après filtrage pour {cat_file.name}. Exemples d'étoiles:")
                            for idx in sample_indices:
                                row = catalog_table[idx]
                                sep = center_coord.separation(SkyCoord(ra=row['RA_ICRS']*u.deg, dec=row['DE_ICRS']*u.deg)).deg
                                logger.warning(f"  Étoile {idx}: RA={row['RA_ICRS']:.2f}°, DEC={row['DE_ICRS']:.2f}°, séparation={sep:.2f}°")
                except Exception as e:
                    logger.error(f"Erreur traitement {cat_file.name}: {e}", exc_info=True)
                    continue
            
            # Si pas assez d'étoiles, élargir la recherche en ajoutant les fichiers les plus proches
            processed_files = set(catalog_files)  # Fichiers déjà traités
            if len(all_stars) < min_stars:
                logger.debug(f"Pas assez d'étoiles ({len(all_stars)} < {min_stars}), ajout de fichiers supplémentaires")
                # Chercher les 2 fichiers les plus proches suivants (sans élargir le rayon de recherche)
                # On utilise l'index pour trouver les fichiers les plus proches sans les charger
                additional_files = self._find_catalog_files(center_coord, search_radius_deg * 3.0, max_files=5)
                
                # Filtrer pour ne garder que les nouveaux fichiers
                new_files = [f for f in additional_files if f not in processed_files]
                
                # Limiter à 2 fichiers supplémentaires pour éviter de charger trop
                for cat_file in new_files[:2]:
                    catalog_table = self._read_1476_file(cat_file)
                    if catalog_table is None or len(catalog_table) == 0:
                        continue
                    
                    try:
                        catalog_coords = SkyCoord(
                            ra=catalog_table['RA_ICRS'] * u.deg,
                            dec=catalog_table['DE_ICRS'] * u.deg
                        )
                        separations = center_coord.separation(catalog_coords).deg
                        # Utiliser le rayon de recherche initial, pas élargi
                        mask = separations <= search_radius_deg
                        
                        if np.sum(mask) > 0:
                            filtered_table = catalog_table[mask]
                            for row in filtered_table:
                                all_stars.append({
                                    'RA_ICRS': row['RA_ICRS'],
                                    'DE_ICRS': row['DE_ICRS'],
                                    'Gmag': row['Gmag']
                                })
                    except Exception as e:
                        logger.error(f"Erreur traitement {cat_file.name}: {e}")
                        continue
                
                catalog_files = list(processed_files) + new_files[:2]  # Mettre à jour pour le log
            else:
                # Mode normal : charger tous les fichiers couvrant la région
                # En blind-solve (rayon >= 90°), charger tous les catalogues disponibles
                if search_radius_deg >= 90.0:
                    # Blind-solve : charger tous les catalogues sans filtrage par distance
                    catalog_files = list(self.catalog_dir.glob("*.1476"))
                    logger.info(f"Mode blind-solve : chargement de tous les catalogues ({len(catalog_files)} fichiers)")
                else:
                    catalog_files = self._find_catalog_files(center_coord, search_radius_deg)
            
            if not catalog_files:
                logger.warning(f"Aucun fichier de catalogue .1476 trouvé pour RA={center_coord.ra.deg:.2f}°, DEC={center_coord.dec.deg:.2f}°")
                return None
            
            all_stars = []
            
            for cat_file in catalog_files:
                catalog_table = self._read_1476_file(cat_file)
                if catalog_table is None or len(catalog_table) == 0:
                    continue
                
                # Filtrer par région (sauf en blind-solve avec rayon très large)
                try:
                    if search_radius_deg >= 90.0:
                        # Blind-solve : pas de filtrage, charger toutes les étoiles du catalogue
                        for row in catalog_table:
                            all_stars.append({
                                'RA_ICRS': row['RA_ICRS'],
                                'DE_ICRS': row['DE_ICRS'],
                                'Gmag': row['Gmag']
                            })
                    else:
                        catalog_coords = SkyCoord(
                            ra=catalog_table['RA_ICRS'] * u.deg,
                            dec=catalog_table['DE_ICRS'] * u.deg
                        )
                        separations = center_coord.separation(catalog_coords).deg
                        mask = separations <= search_radius_deg
                        
                        if np.sum(mask) > 0:
                            filtered_table = catalog_table[mask]
                            for row in filtered_table:
                                all_stars.append({
                                    'RA_ICRS': row['RA_ICRS'],
                                    'DE_ICRS': row['DE_ICRS'],
                                    'Gmag': row['Gmag']
                                })
                except Exception as e:
                    logger.error(f"Erreur traitement {cat_file.name}: {e}")
                    continue
        
        if not all_stars:
            return None
        
        catalog_table = Table(rows=all_stars)
        logger.info(f"Chargé {len(catalog_table)} étoiles depuis {len(catalog_files)} fichier(s) .1476")
        return catalog_table
    
    def _detect_stars(self, data: np.ndarray, fwhm: float = 5.0, threshold_sigma: float = 3.0, max_sources: int = 500) -> Optional[Table]:
        """Détecte les étoiles dans l'image."""
        try:
            # Calculer le seuil
            mean, median, std = np.nanmean(data), np.nanmedian(data), np.nanstd(data)
            threshold = median + threshold_sigma * std
            
            # Détection avec DAOStarFinder
            finder = DAOStarFinder(fwhm=fwhm, threshold=threshold, roundlo=-1.0, roundhi=1.0)
            sources = finder(data)
            
            if sources is None or len(sources) == 0:
                return None
            
            # Trier par flux et limiter
            sources.sort('flux', reverse=True)
            sources = sources[:max_sources]
            
            # Affinage des centroïdes
            refined_x = []
            refined_y = []
            for src in sources:
                try:
                    y_int = int(src['ycentroid'])
                    x_int = int(src['xcentroid'])
                    y_slice = slice(max(0, y_int-5), min(data.shape[0], y_int+6))
                    x_slice = slice(max(0, x_int-5), min(data.shape[1], x_int+6))
                    cutout = data[y_slice, x_slice]
                    
                    if cutout.size < 6:
                        refined_x.append(src['xcentroid'])
                        refined_y.append(src['ycentroid'])
                        continue
                    
                    xc, yc = centroid_2dg(cutout)
                    refined_x.append(src['xcentroid'] - 5 + xc)
                    refined_y.append(src['ycentroid'] - 5 + yc)
                except Exception:
                    refined_x.append(src['xcentroid'])
                    refined_y.append(src['ycentroid'])
            
            sources['xcentroid'] = refined_x
            sources['ycentroid'] = refined_y
            
            return sources
        except Exception as e:
            logger.error(f"Erreur détection étoiles: {e}")
            return None
    
    def _match_stars(self, detected_coords: SkyCoord, catalog_table: Table, match_radius_arcsec: float = 5.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Fait le matching entre étoiles détectées et catalogue."""
        catalog_coords = SkyCoord(ra=catalog_table['RA_ICRS'] * u.deg, dec=catalog_table['DE_ICRS'] * u.deg)
        
        # KD-tree pour matching rapide
        catalog_xyz = catalog_coords.cartesian.xyz.value.T
        kdtree = cKDTree(catalog_xyz)
        
        detected_xyz = detected_coords.cartesian.xyz.value.T
        match_radius_rad = np.radians(match_radius_arcsec / 3600.0)
        
        distances, indices = kdtree.query(detected_xyz, distance_upper_bound=match_radius_rad)
        
        # Convertir distances en arcsec
        valid_mask = indices < len(catalog_table)
        distances_arcsec = np.full(len(detected_coords), np.inf)
        distances_arcsec[valid_mask] = np.degrees(distances[valid_mask]) * 3600.0
        
        pixel_indices = np.where(valid_mask)[0]
        catalog_indices = indices[valid_mask]
        
        return pixel_indices, catalog_indices, distances_arcsec[valid_mask]
    
    def _calculate_wcs(self, pixel_coords: np.ndarray, catalog_coords: SkyCoord, wcs_init: Optional[WCS] = None) -> WCS:
        """Calcule le WCS à partir des matches."""
        if len(pixel_coords) < 3:
            raise ValueError(f"Pas assez de matches ({len(pixel_coords)} < 3) pour calculer WCS")
        
        # Utiliser fit_wcs_from_points
        wcs_fitted = fit_wcs_from_points(
            pixel_coords.T,
            catalog_coords,
            proj_point='center',
            sip_degree=None
        )
        
        return wcs_fitted
    
    def _add_sip_correction(self, wcs: WCS, pixel_coords: np.ndarray, catalog_coords: SkyCoord, order: int = 2) -> WCS:
        """
        Ajoute la correction SIP (Simple Imaging Polynomial) au WCS.
        
        La correction SIP modélise les distorsions optiques non-linéaires avec des polynômes :
        - u' = u + Σ A_i_j * u^i * v^j
        - v' = v + Σ B_i_j * u^i * v^j
        
        où u et v sont les coordonnées pixel normalisées par rapport à CRPIX.
        
        Parameters
        ----------
        wcs : WCS
            WCS linéaire initial
        pixel_coords : np.ndarray
            Coordonnées pixel des matches (N, 2)
        catalog_coords : SkyCoord
            Coordonnées célestes du catalogue
        order : int
            Ordre du polynôme SIP (1 ou 2)
        
        Returns
        -------
        WCS
            WCS avec correction SIP ajoutée
        """
        # Vérifier le nombre minimum de points
        min_points = {1: 6, 2: 10}.get(order, 10)
        if len(pixel_coords) < min_points:
            logger.debug(f"Pas assez de points ({len(pixel_coords)} < {min_points}) pour SIP order {order}")
            return wcs
        
        try:
            # Obtenir CRPIX depuis le WCS
            crpix1 = wcs.wcs.crpix[0]
            crpix2 = wcs.wcs.crpix[1]
            
            # Coordonnées pixel normalisées (u, v) en pixels par rapport à CRPIX
            # Dans le standard SIP, u et v sont en pixels, pas en degrés
            u = pixel_coords[:, 0] - crpix1
            v = pixel_coords[:, 1] - crpix2
            
            # Calculer les résidus en coordonnées tangent planes (en degrés)
            predicted = wcs.pixel_to_world(pixel_coords[:, 0], pixel_coords[:, 1])
            
            # Convertir en coordonnées tangent planes pour le calcul SIP
            # Utiliser la projection tangent plane locale
            center_ra = catalog_coords.ra.deg.mean()
            center_dec = catalog_coords.dec.deg.mean()
            
            # Calculer les résidus en coordonnées tangent planes
            # Delta RA et DEC en degrés
            delta_ra = (predicted.ra.deg - catalog_coords.ra.deg) * np.cos(np.radians(catalog_coords.dec.deg))
            delta_dec = predicted.dec.deg - catalog_coords.dec.deg
            
            # Convertir en coordonnées tangent planes (x, y) en degrés
            # Approximation locale : x ≈ delta_ra, y ≈ delta_dec
            x_residuals = delta_ra
            y_residuals = delta_dec
            
            # Générer les termes polynomiaux pour l'ordre spécifié
            def generate_polynomial_terms(u_vals, v_vals, max_order):
                """Génère les termes polynomiaux u^i * v^j pour i+j <= max_order."""
                terms = []
                for i in range(max_order + 1):
                    for j in range(max_order + 1 - i):
                        if i == 0 and j == 0:
                            continue  # Pas de terme constant (déjà dans le WCS linéaire)
                        terms.append((u_vals ** i) * (v_vals ** j))
                return np.column_stack(terms) if terms else np.zeros((len(u_vals), 0))
            
            # Générer la matrice de design pour les polynômes
            poly_terms = generate_polynomial_terms(u, v, order)
            
            if poly_terms.shape[1] == 0:
                logger.warning("Aucun terme polynomial généré pour SIP")
                return wcs
            
            # Ajuster les coefficients A_i_j pour RA (x)
            # x_corrected = x_linear + Σ A_i_j * u^i * v^j
            try:
                # Utiliser least_squares pour l'ajustement robuste
                def residuals_ra(params):
                    return x_residuals - (poly_terms @ params)
                
                result_ra = least_squares(residuals_ra, np.zeros(poly_terms.shape[1]), method='lm')
                a_coeffs = result_ra.x
                
                # Ajuster les coefficients B_i_j pour DEC (y)
                def residuals_dec(params):
                    return y_residuals - (poly_terms @ params)
                
                result_dec = least_squares(residuals_dec, np.zeros(poly_terms.shape[1]), method='lm')
                b_coeffs = result_dec.x
                
                # Obtenir le header FITS depuis le WCS
                header = wcs.to_header()
                
                # Ajouter les métadonnées SIP au header
                header['A_ORDER'] = order
                header['B_ORDER'] = order
                header['A_DMAX'] = (order, 'Maximum polynomial degree for A coefficients')
                header['B_DMAX'] = (order, 'Maximum polynomial degree for B coefficients')
                
                # Ajouter les coefficients A_i_j et B_i_j au header
                coeff_idx = 0
                for i in range(order + 1):
                    for j in range(order + 1 - i):
                        if i == 0 and j == 0:
                            continue
                        # Format FITS standard : A_i_j et B_i_j
                        header[f'A_{i}_{j}'] = (a_coeffs[coeff_idx], f'SIP distortion coefficient A_{i}_{j}')
                        header[f'B_{i}_{j}'] = (b_coeffs[coeff_idx], f'SIP distortion coefficient B_{i}_{j}')
                        coeff_idx += 1
                
                # Mettre à jour CTYPE pour indiquer SIP
                ctype1 = header.get('CTYPE1', 'RA---TAN')
                ctype2 = header.get('CTYPE2', 'DEC--TAN')
                
                if not ctype1.endswith('-SIP'):
                    if '-TAN' in ctype1:
                        header['CTYPE1'] = ctype1.replace('-TAN', '-TAN-SIP')
                    elif '-SIN' in ctype1:
                        header['CTYPE1'] = ctype1.replace('-SIN', '-SIN-SIP')
                    else:
                        header['CTYPE1'] = ctype1 + '-SIP'
                
                if not ctype2.endswith('-SIP'):
                    if '-TAN' in ctype2:
                        header['CTYPE2'] = ctype2.replace('-TAN', '-TAN-SIP')
                    elif '-SIN' in ctype2:
                        header['CTYPE2'] = ctype2.replace('-SIN', '-SIN-SIP')
                    else:
                        header['CTYPE2'] = ctype2 + '-SIP'
                
                # Créer un nouveau WCS depuis le header avec SIP
                wcs_sip = WCS(header)
                
                # Calculer les résidus après correction SIP pour validation
                predicted_sip = wcs_sip.pixel_to_world(pixel_coords[:, 0], pixel_coords[:, 1])
                delta_ra_sip = (predicted_sip.ra.deg - catalog_coords.ra.deg) * np.cos(np.radians(catalog_coords.dec.deg)) * 3600.0
                delta_dec_sip = (predicted_sip.dec.deg - catalog_coords.dec.deg) * 3600.0
                
                rms_ra_before = np.sqrt(np.mean(x_residuals**2)) * 3600.0
                rms_dec_before = np.sqrt(np.mean(y_residuals**2)) * 3600.0
                rms_ra_after = np.sqrt(np.mean(delta_ra_sip**2))
                rms_dec_after = np.sqrt(np.mean(delta_dec_sip**2))
                
                logger.info(f"Correction SIP order {order} appliquée : "
                          f"RMS RA {rms_ra_before:.2f}→{rms_ra_after:.2f} arcsec, "
                          f"RMS DEC {rms_dec_before:.2f}→{rms_dec_after:.2f} arcsec")
                
                return wcs_sip
                
            except Exception as e:
                logger.warning(f"Erreur ajustement coefficients SIP: {e}")
                return wcs
                
        except Exception as e:
            logger.warning(f"Erreur correction SIP: {e}", exc_info=True)
            return wcs
    
    def solve_file(
        self, 
        fits_path: Path, 
        science_dir: Optional[Path] = None,
        progress_callback: Optional[Callable[[int], None]] = None,
        ra_nearby: Optional[float] = None,
        dec_nearby: Optional[float] = None
    ) -> bool:
        """
        Résout l'astrométrie d'un fichier FITS en utilisant les catalogues locaux.
        
        Parameters
        ----------
        fits_path : Path
            Chemin vers le fichier FITS à résoudre
        science_dir : Path, optional
            Répertoire de sortie pour les images résolues (utilise self.science_dir si None)
        progress_callback : Callable, optional
            Fonction de callback pour la progression (0-100)
        ra_nearby : float, optional
            RA approximative en degrés (pour nearby solve)
        dec_nearby : float, optional
            DEC approximative en degrés (pour nearby solve)
        
        Returns
        -------
        bool
            True si la résolution a réussi, False sinon
        """
        if progress_callback:
            progress_callback(5)
        
        if science_dir is None:
            science_dir = self.science_dir
        else:
            science_dir = Path(science_dir)
            science_dir.mkdir(parents=True, exist_ok=True)
        
        fits_path = Path(fits_path).resolve()
        if not fits_path.exists():
            logger.error(f"Fichier introuvable : {fits_path}")
            return False
        
        # Ignorer les fichiers temporaires
        if '_temp_' in fits_path.name or fits_path.name.endswith('.bak.fits'):
            logger.debug(f"Fichier temporaire ignoré : {fits_path.name}")
            return False
        
        logger.info(f"Résolution astrométrie locale G : {fits_path.name}")
        
        try:
            # 1. Charger l'image FITS
            try:
                with fits.open(fits_path, mode='readonly') as hdul:
                    data = hdul[0].data.astype(float)
                    header = hdul[0].header.copy()
            except (TypeError, OSError, ValueError) as e:
                if "buffer is too small" in str(e) or "too small for requested array" in str(e):
                    logger.error(f"Fichier FITS corrompu ou incomplet : {fits_path.name}")
                    return False
                raise
            
            if progress_callback:
                progress_callback(10)
            
            # 2. Obtenir WCS initial ou position approximative
            wcs_init = None
            center_coord = None
            
            try:
                wcs_init = WCS(header)
                if wcs_init.is_celestial:
                    h, w = data.shape
                    center_coord = wcs_init.pixel_to_world(w/2, h/2)
                    logger.info(f"WCS initial trouvé : RA={center_coord.ra.deg:.6f}°, DEC={center_coord.dec.deg:.6f}°")
            except Exception:
                pass
            
            # Si pas de WCS, utiliser ra_nearby/dec_nearby ou OBJCTRA/OBJCTDEC
            if center_coord is None:
                # Priorité 1 : ra_nearby/dec_nearby fournis explicitement
                if ra_nearby is not None and dec_nearby is not None:
                    center_coord = SkyCoord(ra=ra_nearby * u.deg, dec=dec_nearby * u.deg)
                    logger.info(f"Position approximative fournie (nearby-solve) : RA={ra_nearby:.6f}°, DEC={dec_nearby:.6f}°")
                else:
                    # Priorité 2 : OBJCTRA/OBJCTDEC depuis le header
                    ra_h = header.get("OBJCTRA")
                    dec_h = header.get("OBJCTDEC")
                    if ra_h and dec_h:
                        try:
                            center_coord = SkyCoord(ra=ra_h, dec=dec_h, unit=(u.hourangle, u.deg))
                            logger.info(f"Position depuis OBJCTRA/OBJCTDEC (nearby-solve auto) : RA={center_coord.ra.deg:.6f}°, DEC={center_coord.dec.deg:.6f}°")
                        except Exception as e:
                            logger.debug(f"Impossible de parser OBJCTRA/OBJCTDEC : {e}")
            
            # 3. Déterminer si on est en mode nearby-solve ou blind-solve
            # Nearby-solve si :
            # - ra_nearby/dec_nearby sont fournis explicitement, OU
            # - OBJCTRA/OBJCTDEC sont disponibles (même sans WCS), OU
            # - WCS valide est disponible
            # Blind-solve si aucune position n'est disponible
            has_nearby_hint = (ra_nearby is not None and dec_nearby is not None) or (center_coord is not None)
            nearby_mode = has_nearby_hint and (wcs_init is None or not wcs_init.is_celestial)
            
            # Pour blind-solve sans position, utiliser une position par défaut (centre du ciel visible)
            # et charger tous les catalogues disponibles
            if center_coord is None:
                # Blind-solve : utiliser une position par défaut (équateur, méridien 0h)
                # Le solveur chargera tous les catalogues disponibles
                center_coord = SkyCoord(ra=0.0 * u.deg, dec=0.0 * u.deg)
                logger.info(f"Mode blind-solve : aucune position initiale, utilisation position par défaut (RA=0°, DEC=0°)")
                nearby_mode = False  # En blind-solve, on charge tous les catalogues
            
            # 4. Calculer le rayon de recherche
            h, w = data.shape
            if wcs_init and wcs_init.is_celestial:
                # Calculer FOV depuis WCS
                corners = wcs_init.pixel_to_world([0, w, 0, w], [0, 0, h, h])
                fov_radius = center_coord.separation(corners).max()
                # Pour nearby-solve, utiliser un rayon plus petit (FOV + marge)
                if nearby_mode:
                    search_radius = fov_radius * 1.1  # Petite marge pour nearby-solve
                else:
                    search_radius = fov_radius * 1.2  # Marge normale
            else:
                # Estimation FOV par défaut
                if nearby_mode:
                    # Nearby-solve : rayon plus petit (0.5-1 degré)
                    search_radius = 0.75 * u.deg
                else:
                    # Blind solve : rayon très large pour charger tous les catalogues (180 degrés)
                    search_radius = 180.0 * u.deg
            
            search_radius_deg = search_radius.to(u.deg).value
            
            if nearby_mode:
                logger.info(f"Mode nearby-solve : rayon de recherche = {search_radius_deg:.2f}°")
            else:
                logger.info(f"Mode blind-solve : rayon de recherche = {search_radius_deg:.2f}°")
            
            if progress_callback:
                progress_callback(20)
            
            # 5. Charger les étoiles du catalogue (avec nearby_mode pour optimisation)
            catalog_table = self._load_catalog_stars(center_coord, search_radius_deg, nearby_mode=nearby_mode, min_stars=100)
            if catalog_table is None or len(catalog_table) == 0:
                raise ValueError("Aucune étoile trouvée dans les catalogues pour cette région")
            
            if progress_callback:
                progress_callback(30)
            
            # 6. Détecter les étoiles dans l'image
            sources = self._detect_stars(data, fwhm=5.0, threshold_sigma=3.0, max_sources=500)
            if sources is None or len(sources) == 0:
                raise ValueError("Aucune étoile détectée dans l'image")
            
            logger.info(f"Étoiles détectées : {len(sources)}")
            
            if progress_callback:
                progress_callback(40)
            
            # 7. Créer WCS initial si nécessaire
            if wcs_init is None or not wcs_init.is_celestial:
                # Créer WCS minimal
                cdelt = 2.0 / 3600.0  # 2 arcsec/pixel par défaut
                wcs_init = WCS(naxis=2)
                wcs_init.wcs.crpix = [w/2, h/2]
                wcs_init.wcs.crval = [center_coord.ra.deg, center_coord.dec.deg]
                wcs_init.wcs.cdelt = [-cdelt, cdelt]
                wcs_init.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            
            # 8. Matching initial
            detected_coords = wcs_init.pixel_to_world(sources['xcentroid'], sources['ycentroid'])
            match_radius_arcsec = 10.0
            
            pixel_indices, catalog_indices, distances = self._match_stars(
                detected_coords, catalog_table, match_radius_arcsec
            )
            
            n_matches = len(pixel_indices)
            logger.info(f"Matches initiaux : {n_matches}")
            
            # Si pas assez de matches, augmenter le rayon
            if n_matches < 3:
                match_radius_arcsec = 30.0
                pixel_indices, catalog_indices, distances = self._match_stars(
                    detected_coords, catalog_table, match_radius_arcsec
                )
                n_matches = len(pixel_indices)
                logger.info(f"Matches avec rayon élargi : {n_matches}")
            
            if n_matches < 3:
                raise ValueError(f"Pas assez de matches ({n_matches} < 3) pour calculer WCS")
            
            if progress_callback:
                progress_callback(60)
            
            # 9. Calculer le WCS avec filtrage itératif des outliers
            # Sauvegarder les coordonnées originales pour SIP
            pixel_coords_original = np.array([sources['xcentroid'][pixel_indices], sources['ycentroid'][pixel_indices]]).T
            catalog_coords_original = SkyCoord(
                ra=catalog_table['RA_ICRS'][catalog_indices] * u.deg,
                dec=catalog_table['DE_ICRS'][catalog_indices] * u.deg
            )
            
            # Copies pour l'itération
            pixel_coords = pixel_coords_original.copy()
            catalog_coords = catalog_coords_original
            
            # Fit itératif avec rejet d'outliers pour améliorer la précision
            max_iterations = 3
            outlier_threshold_arcsec = 5.0  # Seuil initial pour rejet d'outliers
            
            best_wcs = None
            best_n_matches = 0
            best_rms = np.inf
            best_outlier_mask = None
            
            for iteration in range(max_iterations):
                if len(pixel_coords) < 3:
                    break
                
                # Calculer WCS
                wcs_fitted = self._calculate_wcs(pixel_coords, catalog_coords, wcs_init)
                
                # Calculer les résidus
                predicted = wcs_fitted.pixel_to_world(pixel_coords[:, 0], pixel_coords[:, 1])
                separations = catalog_coords.separation(predicted)
                residuals_arcsec = separations.arcsec
                
                # Calculer RMS
                rms = np.sqrt(np.mean(residuals_arcsec**2))
                
                # Garder le meilleur fit
                if rms < best_rms and len(pixel_coords) >= 3:
                    best_wcs = wcs_fitted
                    best_n_matches = len(pixel_coords)
                    best_rms = rms
                    # Sauvegarder le masque d'outliers pour utiliser les bons points avec SIP
                    best_outlier_mask = np.ones(len(pixel_coords_original), dtype=bool)
                    # Reconstruire le masque complet depuis les indices filtrés
                    if iteration > 0:
                        # Si on a filtré, on doit reconstruire le masque complet
                        # Pour simplifier, on utilisera tous les points originaux pour SIP
                        best_outlier_mask = np.ones(len(pixel_coords_original), dtype=bool)
                
                # Si c'est la dernière itération, ne pas filtrer
                if iteration == max_iterations - 1:
                    break
                
                # Filtrer les outliers (garder seulement ceux avec résidu < seuil)
                outlier_mask = residuals_arcsec < outlier_threshold_arcsec
                n_outliers = np.sum(~outlier_mask)
                
                if n_outliers == 0:
                    # Pas d'outliers, on peut arrêter
                    break
                
                if np.sum(outlier_mask) < 3:
                    # Trop peu de points après filtrage, garder tous les points
                    break
                
                # Filtrer les coordonnées
                pixel_coords = pixel_coords[outlier_mask]
                catalog_coords = catalog_coords[outlier_mask]
                
                # Réduire le seuil pour la prochaine itération
                outlier_threshold_arcsec = max(2.0, outlier_threshold_arcsec * 0.7)
                
                logger.debug(f"Iteration {iteration + 1}: {n_outliers} outliers rejetés, {len(pixel_coords)} matches restants, RMS={rms:.2f}\"")
            
            if best_wcs is None:
                # Fallback : utiliser le premier fit
                wcs_fitted = self._calculate_wcs(pixel_coords_original, catalog_coords_original, wcs_init)
                pixel_coords_for_sip = pixel_coords_original
                catalog_coords_for_sip = catalog_coords_original
            else:
                wcs_fitted = best_wcs
                logger.info(f"Fit optimal : {best_n_matches} matches, RMS={best_rms:.2f}\"")
                # Utiliser tous les points originaux pour SIP (le fit est déjà optimisé)
                pixel_coords_for_sip = pixel_coords_original
                catalog_coords_for_sip = catalog_coords_original
            
            # Calculer les résidus finaux pour logging
            predicted_final = wcs_fitted.pixel_to_world(pixel_coords_for_sip[:, 0], pixel_coords_for_sip[:, 1])
            separations_final = catalog_coords_for_sip.separation(predicted_final)
            rms_final = np.sqrt(np.mean(separations_final.arcsec**2))
            logger.info(f"WCS calculé : {len(pixel_coords_for_sip)} matches, RMS={rms_final:.2f}\"")
            
            # 10. Ajouter correction SIP pour améliorer la précision
            # Utiliser order=1 si moins de 10 matches, order=2 sinon
            if len(pixel_coords_for_sip) >= 6:
                sip_order = 1 if len(pixel_coords_for_sip) < 10 else 2
                wcs_fitted = self._add_sip_correction(wcs_fitted, pixel_coords_for_sip, catalog_coords_for_sip, order=sip_order)
            
            if progress_callback:
                progress_callback(80)
            
            # 11. Sauvegarder l'image résolue
            output_path = science_dir / f"{fits_path.stem}_solved.fits"
            
            # Lire l'image originale en readonly
            with fits.open(fits_path, mode='readonly') as hdul_in:
                hdul_out = fits.HDUList()
                for i, hdu in enumerate(hdul_in):
                    if i == 0:
                        # Mettre à jour le header avec le nouveau WCS
                        new_header = hdu.header.copy()
                        wcs_header = wcs_fitted.to_header()
                        new_header.update(wcs_header)
                        
                        # Extraire directement la matrice CD du WCS
                        # Le WCS peut utiliser PC+CDELT ou CD directement
                        # On calcule toujours la matrice CD complète pour garantir la précision
                        try:
                            # Méthode 1 : Si le WCS a déjà une matrice CD, l'utiliser
                            if hasattr(wcs_fitted.wcs, 'cd') and wcs_fitted.wcs.cd is not None:
                                wcs_cd = wcs_fitted.wcs.cd
                                if wcs_cd.shape == (2, 2):
                                    cd_matrix = wcs_cd
                                else:
                                    cd_matrix = None
                            else:
                                cd_matrix = None
                            
                            # Méthode 2 : Calculer depuis PC + CDELT si CD n'existe pas
                            if cd_matrix is None:
                                if hasattr(wcs_fitted.wcs, 'pc') and wcs_fitted.wcs.pc is not None:
                                    pc = wcs_fitted.wcs.pc
                                    cdelt = wcs_fitted.wcs.cdelt
                                    if pc.shape == (2, 2) and cdelt is not None and len(cdelt) == 2:
                                        # CD = PC * diag(CDELT)
                                        cd_matrix = np.array([
                                            [pc[0, 0] * cdelt[0], pc[0, 1] * cdelt[0]],
                                            [pc[1, 0] * cdelt[1], pc[1, 1] * cdelt[1]]
                                        ])
                            
                            # Méthode 3 : Calculer depuis pixel_scale_matrix si disponible
                            if cd_matrix is None:
                                try:
                                    # pixel_scale_matrix donne la transformation pixel -> world en degrés
                                    # On doit la convertir en matrice CD
                                    pixel_scale = wcs_fitted.pixel_scale_matrix
                                    if pixel_scale is not None and pixel_scale.shape == (2, 2):
                                        # pixel_scale est en degrés/pixel, c'est exactement ce qu'on veut pour CD
                                        cd_matrix = pixel_scale
                                except:
                                    pass
                            
                            # Si on a une matrice CD, l'utiliser
                            if cd_matrix is not None:
                                new_header['CD1_1'] = (cd_matrix[0, 0], 'Matrix element 1,1')
                                new_header['CD1_2'] = (cd_matrix[0, 1], 'Matrix element 1,2')
                                new_header['CD2_1'] = (cd_matrix[1, 0], 'Matrix element 2,1')
                                new_header['CD2_2'] = (cd_matrix[1, 1], 'Matrix element 2,2')
                                
                                # Supprimer les anciennes clés PC et CDELT
                                for key in ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2']:
                                    if key in new_header:
                                        del new_header[key]
                                
                                # Ajouter CUNIT si nécessaire
                                if 'CUNIT1' not in new_header:
                                    new_header['CUNIT1'] = ('deg', 'Units of coordinate increment and value')
                                if 'CUNIT2' not in new_header:
                                    new_header['CUNIT2'] = ('deg', 'Units of coordinate increment and value')
                            else:
                                # Fallback : utiliser to_header() mais corriger CDELT si nécessaire
                                if 'CDELT1' in new_header and abs(new_header['CDELT1']) > 0.1:
                                    # CDELT semble correct, utiliser PC*CDELT
                                    if 'PC1_1' in new_header:
                                        pc1_1 = new_header['PC1_1']
                                        pc1_2 = new_header.get('PC1_2', 0.0)
                                        pc2_1 = new_header.get('PC2_1', 0.0)
                                        pc2_2 = new_header.get('PC2_2', 0.0)
                                        cdelt1 = new_header['CDELT1']
                                        cdelt2 = new_header.get('CDELT2', cdelt1)
                                        
                                        new_header['CD1_1'] = (pc1_1 * cdelt1, 'Matrix element 1,1')
                                        new_header['CD1_2'] = (pc1_2 * cdelt1, 'Matrix element 1,2')
                                        new_header['CD2_1'] = (pc2_1 * cdelt2, 'Matrix element 2,1')
                                        new_header['CD2_2'] = (pc2_2 * cdelt2, 'Matrix element 2,2')
                                        
                                        for key in ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1', 'CDELT2']:
                                            if key in new_header:
                                                del new_header[key]
                        except Exception as e:
                            logger.warning(f"Erreur lors de l'extraction de la matrice CD : {e}, utilisation du header par défaut")
                        
                        new_header['PLTSOLVD'] = (True, 'Solved using local Gaia DR3 catalogues')
                        new_header.add_history('Plate-solved via LocalAstrometrySolver (Gaia DR3 G-magnitude)')
                        hdul_out.append(fits.PrimaryHDU(data=hdu.data, header=new_header))
                    else:
                        hdul_out.append(hdu)
            
            # Sauvegarder avec gestion robuste des verrous de fichiers
            max_retries = 5
            for attempt in range(max_retries):
                try:
                    temp_path = output_path.parent / f"{output_path.stem}_temp_{int(time.time() * 1000)}.fits"
                    hdul_out.writeto(temp_path, overwrite=True, output_verify='ignore')
                    hdul_out.close()
                    
                    # Remplacer le fichier final
                    if output_path.exists():
                        try:
                            output_path.unlink()
                        except PermissionError:
                            time.sleep(0.1)
                            output_path.unlink()
                    
                    temp_path.rename(output_path)
                    break
                except (PermissionError, OSError) as e:
                    if attempt < max_retries - 1:
                        time.sleep(0.2)
                        continue
                    else:
                        raise
            
            if progress_callback:
                progress_callback(100)
            
            logger.info(f"✅ Image résolue sauvegardée : {output_path.name}")
            return True
            
        except Exception as e:
            logger.error(f"Erreur résolution {fits_path.name}: {e}", exc_info=True)
            return False
            
            
















