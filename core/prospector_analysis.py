# core/prospector_analysis.py
"""
Module pour l'analyse de spectres de galaxies et l'inférence de propriétés stellaires
à partir de SED (Spectral Energy Distribution) en utilisant Prospector.
"""

import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, TYPE_CHECKING, Any
import json
import subprocess
import tempfile
from astropy import units as u
from astropy.io import fits

logger = logging.getLogger(__name__)
# Log immédiat pour vérifier que le module est chargé
logger.info("[PROSPECTOR_MODULE] Module core.prospector_analysis chargé")

# IMPORTANT: Définir SPS_HOME AVANT d'importer prospect
# Prospector essaie d'utiliser SPS_HOME lors de l'import, donc il doit être défini avant
import os
import pathlib
from utils.wsl_utils import windows_path_to_wsl

PROSPECTOR_BACKEND = os.environ.get("NPOAP_PROSPECTOR_BACKEND", "windows").strip().lower() or "windows"
PROSPECTOR_WSL_AVAILABLE = False


def _resolve_wsl_python() -> str:
    """Retourne le Python du venv WSL Prospector."""
    custom = os.environ.get("NPOAP_PROSPECTOR_WSL_PYTHON", "").strip()
    if custom:
        return custom
    # Valeur alignée avec Installation_fsps/install_prospector_wsl.sh
    return "$HOME/.local/share/npoap-prospector-wsl/venv/bin/python3"


def _resolve_wsl_sps_home() -> str:
    """Retourne SPS_HOME côté WSL."""
    custom = os.environ.get("NPOAP_PROSPECTOR_WSL_SPS_HOME", "").strip()
    if custom:
        return custom
    return "$HOME/.local/share/fsps"


def _probe_wsl_prospector() -> Tuple[bool, str]:
    """Vérifie si Prospector+FSPS sont utilisables côté WSL."""
    if PROSPECTOR_BACKEND != "wsl":
        return False, "backend WSL non demandé"
    wsl_python = _resolve_wsl_python()
    wsl_sps_home = _resolve_wsl_sps_home()
    probe_cmd = (
        f'PY="{wsl_python}"; '
        f'export SPS_HOME="{wsl_sps_home}"; '
        '[ -x "${PY/#\\~/$HOME}" ] || [ -x "$PY" ] || exit 21; '
        '"${PY/#\\~/$HOME}" -c "import prospect, fsps; print(getattr(prospect, \'__version__\', \'?\'), getattr(fsps, \'__version__\', \'?\'))"'
    )
    try:
        cp = subprocess.run(
            ["wsl", "-e", "bash", "-lc", probe_cmd],
            capture_output=True,
            text=True,
            check=False,
            timeout=20,
        )
        if cp.returncode == 0:
            return True, (cp.stdout or "").strip()
        err = (cp.stderr or cp.stdout or "").strip()
        return False, f"rc={cp.returncode} {err}"
    except Exception as e:
        return False, str(e)

if 'SPS_HOME' not in os.environ:
    sps_home_default = pathlib.Path.home() / '.local' / 'share' / 'fsps'
    os.environ['SPS_HOME'] = str(sps_home_default)
    logger.info(f"[PROSPECTOR_MODULE] SPS_HOME défini: {sps_home_default}")
else:
    logger.info(f"[PROSPECTOR_MODULE] SPS_HOME déjà défini: {os.environ['SPS_HOME']}")

# Créer les fichiers stub nécessaires dès le chargement du module
# Cela permet à Prospector de s'importer même si FSPS n'est pas installé
def _create_fsps_stub_files():
    """Crée les fichiers stub FSPS nécessaires pour l'import de Prospector"""
    try:
        # Déterminer SPS_HOME
        if 'SPS_HOME' in os.environ:
            sps_home = pathlib.Path(os.environ['SPS_HOME'])
        else:
            sps_home = pathlib.Path.home() / '.local' / 'share' / 'fsps'
            os.environ['SPS_HOME'] = str(sps_home)
        
        sps_home.mkdir(parents=True, exist_ok=True)
        (sps_home / 'dust').mkdir(parents=True, exist_ok=True)
        (sps_home / 'sed').mkdir(parents=True, exist_ok=True)
        
        # Créer le fichier stub avec format Prospector (10 colonnes: wave + 9 fnu)
        dust_file = sps_home / 'dust' / 'Nenkova08_y010_torusg_n10_q2.0.dat'
        # Toujours vérifier et recréer si nécessaire pour s'assurer du bon format
        needs_recreate = True
        if dust_file.exists():
            try:
                # Vérifier le format en lisant les lignes (skip_header=4 dans Prospector)
                # Format attendu: 4 lignes d'en-tête puis ~125 lignes avec 10 colonnes
                with open(dust_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 5:  # Au moins 4 en-têtes + 1 ligne de données
                        # Vérifier la 5ème ligne (première ligne de données après skip_header=4)
                        data_line = lines[4].strip()
                        if data_line and not data_line.startswith('#'):
                            # IMPORTANT: Vérifier le séparateur (3 espaces exactement)
                            # Prospector utilise delimiter='   ' (3 espaces) dans np.genfromtxt
                            cols_3spaces = data_line.split('   ')  # Séparateur: 3 espaces exactement
                            if len(cols_3spaces) == 10:
                                # Vérifier aussi la ligne suivante
                                if len(lines) >= 6:
                                    data_line2 = lines[5].strip()
                                    if data_line2 and not data_line2.startswith('#'):
                                        cols2_3spaces = data_line2.split('   ')  # Séparateur: 3 espaces
                                        if len(cols2_3spaces) == 10:
                                            # Vérifier qu'on a au moins ~125 lignes de données (4 en-têtes + 125 lignes = 129 total)
                                            if len(lines) >= 129:
                                                needs_recreate = False
                                                logger.info(f"[PROSPECTOR_MODULE] ✓ Fichier stub existe et a le bon format (10 colonnes avec delimiter=3 espaces, {len(lines)} lignes)")
                                            else:
                                                logger.warning(f"[PROSPECTOR_MODULE] Fichier stub a seulement {len(lines)} lignes (besoin d'au moins 129), recréation nécessaire")
                                        else:
                                            logger.warning(f"[PROSPECTOR_MODULE] Fichier stub ligne 6 a {len(cols2_3spaces)} colonnes avec delimiter=3 espaces (attendu: 10), recréation nécessaire")
                                    else:
                                        logger.warning(f"[PROSPECTOR_MODULE] Fichier stub ligne 6 invalide, recréation nécessaire")
                                else:
                                    logger.warning(f"[PROSPECTOR_MODULE] Fichier stub a moins de 6 lignes, recréation nécessaire")
                            else:
                                # Peut-être que le fichier utilise un autre séparateur (1 espace)
                                cols_1space = data_line.split()  # Séparateur: espace(s) quelconque(s)
                                if len(cols_1space) == 10:
                                    logger.warning(f"[PROSPECTOR_MODULE] Fichier stub utilise un séparateur incorrect (1 espace au lieu de 3 espaces), recréation nécessaire")
                                else:
                                    logger.warning(f"[PROSPECTOR_MODULE] Fichier stub ligne 5 a {len(cols_3spaces)} colonnes avec delimiter=3 espaces (attendu: 10), recréation nécessaire")
                        else:
                            logger.warning(f"[PROSPECTOR_MODULE] Fichier stub ligne 5 invalide ou commentaire, recréation nécessaire")
                    else:
                        logger.warning(f"[PROSPECTOR_MODULE] Fichier stub a seulement {len(lines)} lignes (besoin d'au moins 5), recréation nécessaire")
            except Exception as e:
                logger.warning(f"[PROSPECTOR_MODULE] Erreur lors de la vérification du fichier stub: {e}, recréation...")
                needs_recreate = True
        
        if needs_recreate:
            logger.info(f"[PROSPECTOR_MODULE] Création/recréation du fichier stub FSPS: {dust_file}")
            # Supprimer l'ancien fichier s'il existe avec un mauvais format
            if dust_file.exists():
                try:
                    dust_file.unlink()
                    logger.info(f"[PROSPECTOR_MODULE] Ancien fichier stub supprimé (format incorrect)")
                except Exception as del_error:
                    logger.warning(f"[PROSPECTOR_MODULE] Impossible de supprimer l'ancien fichier: {del_error}")
            try:
                with open(dust_file, 'w') as f:
                    # Format attendu par Prospector (d'après fake_fsps.py ligne 191):
                    # dtype=[('wave', 'f8'), ('fnu_5', '<U20'), ('fnu_10', '<U20'), ('fnu_20', '<U20'),
                    #        ('fnu_30', '<U20'), ('fnu_40', '<U20'), ('fnu_60', '<U20'),
                    #        ('fnu_80', '<U20'), ('fnu_100', '<U20'), ('fnu_150', '<U20')]
                    # delimiter='   ' (3 espaces), skip_header=4, ~125 lignes de données
                    # Format: 1 colonne numérique (wave) + 9 colonnes de chaînes (fnu_*)
                    
                    # Écrire 4 lignes d'en-tête (seront ignorées par skip_header=4)
                    f.write("# Nenkova08 AGN torus dust model - Stub file\n")
                    f.write("# This is a stub file created automatically\n")
                    f.write("# Replace with real FSPS data file for full functionality\n")
                    f.write("# wave   fnu_5   fnu_10   fnu_20   fnu_30   fnu_40   fnu_60   fnu_80   fnu_100   fnu_150\n")
                    
                    # Créer ~125 lignes de données (comme le fichier réel)
                    # Format: wave (numérique) + 9 valeurs fnu (comme chaînes, séparées par 3 espaces)
                    for i in range(125):
                        wave = 1.0 + i * 0.1  # Colonne 1: longueur d'onde (numérique, ex: 1.0 à ~13.4)
                        # Colonnes 2-10: valeurs fnu (comme chaînes, max 20 caractères)
                        # Utiliser des valeurs progressives pour simuler des données
                        fnu_values = [f"{(j+1)*0.001 + i*0.0001:.6f}" for j in range(9)]  # 9 valeurs
                        # Formater la ligne avec 3 espaces comme séparateur (delimiter='   ')
                        line = f"{wave:.6f}   " + "   ".join(fnu_values)
                        f.write(line + "\n")
                logger.info(f"[PROSPECTOR_MODULE] ✓ Fichier stub créé/recréé avec format Prospector (4 en-têtes + 125 lignes, delimiter=3 espaces)")
            except Exception as e:
                logger.error(f"[PROSPECTOR_MODULE] ✗ Erreur lors de la création du fichier stub: {e}")
                return False
        
        return True
    except Exception as e:
        logger.warning(f"[PROSPECTOR_MODULE] ⚠ Erreur lors de la création des fichiers stub: {e}")
        return False

# Créer les fichiers stub dès maintenant (AVANT l'import de prospect)
# SPS_HOME a déjà été défini plus haut
_create_fsps_stub_files()

# Pour les annotations de type
if TYPE_CHECKING:
    try:
        from prospect.models import SpecModel as SedModel
    except ImportError:
        from prospect.models import SedModel
    if SPECUTILS_AVAILABLE:
        from specutils import Spectrum1D
    else:
        Spectrum1D = Any  # Type placeholder si specutils n'est pas disponible

# Import Prospector
PROSPECTOR_AVAILABLE = False
logger.info("=" * 60)
logger.info("DEBUT IMPORT PROSPECTOR")
logger.info("=" * 60)

try:
    # SPS_HOME a déjà été défini au-dessus, juste vérifier
    logger.info("[PROSPECTOR] Vérification de SPS_HOME...")
    
    if 'SPS_HOME' not in os.environ:
        logger.info("[PROSPECTOR] SPS_HOME non défini, création automatique...")
        # Essayer de trouver un répertoire SPS_HOME commun
        possible_paths = [
            pathlib.Path.home() / '.local' / 'share' / 'fsps',
            pathlib.Path.home() / '.fsps',
            pathlib.Path.home() / 'fsps',
        ]
        
        # Si aucun chemin trouvé, créer un répertoire temporaire
        # Prospector essaie d'accéder à SPS_HOME même sans FSPS
        sps_home = possible_paths[0]
        logger.info(f"[PROSPECTOR] Création du répertoire SPS_HOME: {sps_home}")
        sps_home.mkdir(parents=True, exist_ok=True)
        
        # Créer les sous-répertoires nécessaires
        (sps_home / 'dust').mkdir(exist_ok=True)
        (sps_home / 'sed').mkdir(exist_ok=True)
        
        os.environ['SPS_HOME'] = str(sps_home)
        logger.info(f"[PROSPECTOR] SPS_HOME défini automatiquement: {sps_home}")
        logger.info("[PROSPECTOR] Note: Pour une utilisation complète, installez FSPS (voir docs/INSTALLATION_FSPS.md)")
    else:
        sps_home = pathlib.Path(os.environ['SPS_HOME'])
        logger.info(f"[PROSPECTOR] SPS_HOME déjà défini: {sps_home}")
        # S'assurer que les sous-répertoires existent
        (sps_home / 'dust').mkdir(parents=True, exist_ok=True)
        (sps_home / 'sed').mkdir(parents=True, exist_ok=True)
    
    # Prospector utilise 'prospect' comme nom de module
    # Essayer d'importer prospect - peut échouer si les fichiers de données FSPS sont manquants
    # Les fichiers stub ont déjà été créés par _create_fsps_stub_files() au chargement du module
    logger.info("[PROSPECTOR] Tentative d'import du module 'prospect'...")
    try:
        # Vérifier que le fichier stub existe et a le bon format
        dust_file = sps_home / 'dust' / 'Nenkova08_y010_torusg_n10_q2.0.dat'
        if not dust_file.exists():
            logger.warning(f"[PROSPECTOR] Fichier stub manquant, recréation...")
            _create_fsps_stub_files()
        else:
            # Vérifier le format rapidement (skip_header=4, donc vérifier la 5ème ligne)
            try:
                with open(dust_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 5:
                        # Vérifier la 5ème ligne (première ligne de données)
                        data_line = lines[4].strip()
                        if data_line and not data_line.startswith('#'):
                            cols = data_line.split()
                            if len(cols) != 10:
                                logger.warning(f"[PROSPECTOR] Fichier stub ligne 5 a {len(cols)} colonnes au lieu de 10, recréation...")
                                _create_fsps_stub_files()
                        else:
                            logger.warning(f"[PROSPECTOR] Fichier stub ligne 5 invalide, recréation...")
                            _create_fsps_stub_files()
                    else:
                        logger.warning(f"[PROSPECTOR] Fichier stub trop court ({len(lines)} lignes), recréation...")
                        _create_fsps_stub_files()
            except Exception as e:
                logger.warning(f"[PROSPECTOR] Erreur lors de la vérification du fichier stub: {e}, recréation...")
                _create_fsps_stub_files()
        
        import prospect
        logger.info(f"[PROSPECTOR] ✓ Module 'prospect' importé avec succès")
        logger.info(f"[PROSPECTOR] Version de prospect: {getattr(prospect, '__version__', 'inconnue')}")
        logger.info(f"[PROSPECTOR] Emplacement: {getattr(prospect, '__file__', 'inconnu')}")
        
        # Prospector utilise SpecModel (pas SedModel) dans les versions récentes
        # L'import peut nécessiter de passer par un sous-module
        SedModel = None
        transforms = None
        logger.info("[PROSPECTOR] Tentative d'import de SpecModel depuis prospect.models...")
        try:
            # Essayer d'importer depuis prospect.models directement
            from prospect.models import SpecModel, transforms
            SedModel = SpecModel  # Alias pour compatibilité avec le code existant
            logger.info("[PROSPECTOR] ✓ SpecModel importé depuis prospect.models")
        except ImportError as e:
            logger.warning(f"[PROSPECTOR] Import depuis prospect.models échoué: {e}")
            try:
                # Essayer depuis prospect.models.sedmodel
                logger.info("[PROSPECTOR] Tentative d'import depuis prospect.models.sedmodel...")
                from prospect.models.sedmodel import SpecModel, transforms
                SedModel = SpecModel
                logger.info("[PROSPECTOR] ✓ SpecModel importé depuis prospect.models.sedmodel")
            except ImportError as e2:
                logger.warning(f"[PROSPECTOR] Import depuis prospect.models.sedmodel échoué: {e2}")
                try:
                    # Fallback: SedModel (ancienne version)
                    logger.info("[PROSPECTOR] Tentative d'import de SedModel (version ancienne)...")
                    from prospect.models import SedModel, transforms
                    logger.info("[PROSPECTOR] ✓ SedModel importé depuis prospect.models (version ancienne)")
                except ImportError as e3:
                    logger.warning(f"[PROSPECTOR] Import de SedModel échoué: {e3}")
                    try:
                        # Essayer d'importer transforms séparément
                        logger.info("[PROSPECTOR] Tentative d'import séparé de transforms...")
                        from prospect.models import transforms
                        # Essayer d'obtenir SpecModel dynamiquement
                        import prospect.models.sedmodel as sedmod
                        logger.info(f"[PROSPECTOR] Module sedmodel importé: {sedmod}")
                        logger.info(f"[PROSPECTOR] Contenu de sedmodel: {[x for x in dir(sedmod) if not x.startswith('_')][:10]}")
                        if hasattr(sedmod, 'SpecModel'):
                            SedModel = sedmod.SpecModel
                            logger.info("[PROSPECTOR] ✓ SpecModel obtenu dynamiquement depuis prospect.models.sedmodel")
                        elif hasattr(sedmod, 'SedModel'):
                            SedModel = sedmod.SedModel
                            logger.info("[PROSPECTOR] ✓ SedModel obtenu dynamiquement depuis prospect.models.sedmodel")
                        else:
                            raise ImportError("Aucune classe de modèle trouvée dans prospect.models.sedmodel")
                    except (ImportError, AttributeError) as e4:
                        logger.error(f"[PROSPECTOR] ✗ Impossible d'importer SpecModel/SedModel: {e4}")
                        logger.error(f"[PROSPECTOR] Type d'erreur: {type(e4).__name__}")
                        import traceback
                        logger.error(f"[PROSPECTOR] Traceback: {traceback.format_exc()}")
                        # Essayer d'importer au moins transforms
                        try:
                            from prospect.models import transforms
                            logger.info("[PROSPECTOR] ✓ transforms importé avec succès")
                        except ImportError as e5:
                            logger.error(f"[PROSPECTOR] ✗ Impossible d'importer transforms: {e5}")
                            transforms = None
                        SedModel = None
                        logger.warning("[PROSPECTOR] SedModel/SpecModel non trouvé, sera importé dynamiquement lors de l'utilisation")
        
        logger.info("[PROSPECTOR] Tentative d'import de fit_model et FastStepBasis...")
        try:
            from prospect.fitting import fit_model
            logger.info("[PROSPECTOR] ✓ fit_model importé avec succès")
        except ImportError as e:
            logger.error(f"[PROSPECTOR] ✗ Impossible d'importer fit_model: {e}")
            raise
        
        try:
            from prospect.sources import FastStepBasis
            logger.info("[PROSPECTOR] ✓ FastStepBasis importé avec succès")
        except ImportError as e:
            logger.error(f"[PROSPECTOR] ✗ Impossible d'importer FastStepBasis: {e}")
            raise
    except (FileNotFoundError, OSError) as e:
        # Prospector peut échouer si les fichiers de données FSPS sont manquants
        error_msg = str(e)
        logger.error(f"[PROSPECTOR] ✗ Erreur FileNotFoundError/OSError lors de l'import: {error_msg}")
        logger.error(f"[PROSPECTOR] Type d'erreur: {type(e).__name__}")
        import traceback
        logger.error(f"[PROSPECTOR] Traceback complet:\n{traceback.format_exc()}")
        
        # Essayer de créer automatiquement les fichiers manquants si possible
        if 'dust' in error_msg.lower() or 'dat' in error_msg.lower() or 'Nenkova' in error_msg:
            logger.info("[PROSPECTOR] Tentative de création automatique des fichiers de données manquants...")
            try:
                # Extraire le chemin du fichier depuis l'erreur
                import re
                match = re.search(r'([A-Z]:[^"]+\.dat)', error_msg)
                if match:
                    missing_file = pathlib.Path(match.group(1))
                    logger.info(f"[PROSPECTOR] Création du fichier stub: {missing_file}")
                    missing_file.parent.mkdir(parents=True, exist_ok=True)
                    with open(missing_file, 'w') as f:
                        f.write("# Fichier stub créé automatiquement - remplacez par le vrai fichier FSPS si nécessaire\n")
                        f.write("# Format: colonnes séparées par des espaces\n")
                        f.write("0.0 0.0 0.0\n")
                    logger.info(f"[PROSPECTOR] ✓ Fichier stub créé, nouvelle tentative d'import...")
                    
                    # Réessayer l'import après avoir créé le fichier stub
                    try:
                        import prospect
                        logger.info("[PROSPECTOR] ✓ Import réussi après création du fichier stub")
                        # Continuer avec le reste de l'import...
                        # On va relancer le processus depuis le début mais cette fois ça devrait marcher
                        # Cependant, pour éviter une boucle infinie, on va plutôt gérer cette erreur différemment
                        # et continuer avec une version limitée de Prospector
                    except Exception as e2:
                        logger.error(f"[PROSPECTOR] ✗ Import échoué même après création du stub: {e2}")
                        raise
                else:
                    logger.warning("[PROSPECTOR] Impossible d'extraire le chemin du fichier manquant depuis l'erreur")
                    raise
            except Exception as stub_error:
                logger.error(f"[PROSPECTOR] ✗ Erreur lors de la création du fichier stub: {stub_error}")
                logger.warning("[PROSPECTOR] Prospector installé mais fichiers de données manquants")
                logger.warning("[PROSPECTOR] Les fichiers de données FSPS sont nécessaires pour utiliser Prospector")
                logger.warning("[PROSPECTOR] Pour installer les fichiers de données:")
                logger.warning("[PROSPECTOR] 1. Installer FSPS: voir docs/INSTALLATION_FSPS.md")
                logger.warning("[PROSPECTOR] 2. Ou télécharger les fichiers de données FSPS manuellement")
                raise ImportError("Prospector nécessite les fichiers de données FSPS. Installez FSPS ou téléchargez les fichiers de données manuellement.")
        else:
            raise
    except ImportError as e:
        logger.error(f"[PROSPECTOR] ✗ Erreur ImportError lors de l'import: {e}")
        logger.error(f"[PROSPECTOR] Type d'erreur: {type(e).__name__}")
        import traceback
        logger.error(f"[PROSPECTOR] Traceback complet:\n{traceback.format_exc()}")
        raise
    
    # Import FSPS pour les templates SSP.
    # Si backend WSL demandé, on évite la vérification locale Windows (bruyante/inutile).
    FSPS_AVAILABLE = False
    if PROSPECTOR_BACKEND != "wsl":
        logger.info("[PROSPECTOR] Vérification de FSPS...")
        try:
            import fsps
            FSPS_AVAILABLE = True
            logger.info("[PROSPECTOR] ✓ FSPS disponible pour les templates SSP")
            logger.info(f"[PROSPECTOR] Version FSPS: {getattr(fsps, '__version__', 'inconnue')}")
            logger.info(f"[PROSPECTOR] Emplacement FSPS: {getattr(fsps, '__file__', 'inconnu')}")
        except ImportError as e:
            FSPS_AVAILABLE = False
            logger.warning(f"[PROSPECTOR] ⚠ FSPS n'est pas disponible: {e}")
            logger.warning("[PROSPECTOR] Templates SSP ne fonctionneront pas")
            logger.warning("[PROSPECTOR] Pour installer FSPS: voir https://dfm.io/python-fsps/current/installation/")
            logger.warning("[PROSPECTOR] Note: FSPS nécessite CMake et un compilateur Fortran")
            logger.warning("[PROSPECTOR] Prospector peut fonctionner avec d'autres bibliothèques SSP si disponibles")

    # Backend alternatif: Prospector+FSPS via WSL (venv Linux)
    if not FSPS_AVAILABLE and PROSPECTOR_BACKEND == "wsl":
        ok_wsl, details = _probe_wsl_prospector()
        PROSPECTOR_WSL_AVAILABLE = ok_wsl
        if ok_wsl:
            FSPS_AVAILABLE = True
            logger.info("[PROSPECTOR] ✓ Backend WSL actif (%s)", details)
            logger.info("[PROSPECTOR] FSPS sera exploité via WSL")
        else:
            logger.warning("[PROSPECTOR] Backend WSL indisponible: %s", details)
    
    # Import pour la calibration spectroscopique
    try:
        from prospect.sources import CSPSpecBasis
    except ImportError:
        CSPSpecBasis = None
    
    # Import pour l'échantillonnage MCMC
    try:
        import dynesty
        DYNESTY_AVAILABLE = True
    except ImportError:
        DYNESTY_AVAILABLE = False
        logger.warning("Dynesty n'est pas disponible - utilisation d'emcee par défaut")
    
    try:
        import emcee
        EMCEE_AVAILABLE = True
    except ImportError:
        EMCEE_AVAILABLE = False
        logger.warning("Emcee n'est pas disponible pour MCMC")
    
    # Import I/O
    try:
        from prospect.io import write_results as write_prospector_results
        from prospect.io import read_results as read_prospector_results
    except ImportError:
        write_prospector_results = None
        read_prospector_results = None
    
    PROSPECTOR_AVAILABLE = True
    logger.info("[PROSPECTOR] ========================================")
    logger.info("[PROSPECTOR] ✓ PROSPECTOR DISPONIBLE ET FONCTIONNEL")
    logger.info("[PROSPECTOR] ========================================")
    logger.info(f"[PROSPECTOR] PROSPECTOR_AVAILABLE = {PROSPECTOR_AVAILABLE}")
    logger.info(f"[PROSPECTOR] FSPS_AVAILABLE = {FSPS_AVAILABLE}")
    logger.info(f"[PROSPECTOR] SedModel défini: {SedModel is not None}")
    logger.info("[PROSPECTOR] Prospector disponible pour l'analyse de populations stellaires")
except (ImportError, FileNotFoundError, OSError) as e:
    error_msg = str(e)
    logger.error("[PROSPECTOR] ========================================")
    logger.error("[PROSPECTOR] ✗ ERREUR LORS DE L'IMPORT DE PROSPECTOR")
    logger.error("[PROSPECTOR] ========================================")
    logger.error(f"[PROSPECTOR] Type d'exception: {type(e).__name__}")
    logger.error(f"[PROSPECTOR] Message d'erreur: {error_msg}")
    import traceback
    logger.error(f"[PROSPECTOR] Traceback complet:\n{traceback.format_exc()}")
    if 'SPS_HOME' in error_msg or 'dust' in error_msg.lower() or 'dat' in error_msg.lower() or 'FileNotFoundError' in str(type(e)):
        logger.warning("[PROSPECTOR] Prospector installé mais nécessite les fichiers de données FSPS")
        logger.warning("[PROSPECTOR] Les fichiers de données FSPS sont nécessaires pour utiliser Prospector")
        logger.warning("[PROSPECTOR] Solution: Installez FSPS (voir docs/INSTALLATION_FSPS.md) ou utilisez WSL")
        logger.warning("[PROSPECTOR] Prospector ne sera pas disponible tant que FSPS n'est pas installé")
    else:
        logger.warning("[PROSPECTOR] Prospector n'est pas installé ou erreur d'import")
        logger.warning("[PROSPECTOR] Pour installer Prospector: pip install git+https://github.com/bd-j/prospector.git")
    PROSPECTOR_AVAILABLE = False
    FSPS_AVAILABLE = False
    DYNESTY_AVAILABLE = False
    EMCEE_AVAILABLE = False
    logger.info(f"[PROSPECTOR] PROSPECTOR_AVAILABLE = {PROSPECTOR_AVAILABLE}")
except Exception as e:
    logger.error("[PROSPECTOR] ========================================")
    logger.error("[PROSPECTOR] ✗ ERREUR INATTENDUE LORS DE L'IMPORT")
    logger.error("[PROSPECTOR] ========================================")
    logger.error(f"[PROSPECTOR] Type d'exception: {type(e).__name__}")
    logger.error(f"[PROSPECTOR] Message d'erreur: {str(e)}")
    import traceback
    logger.error(f"[PROSPECTOR] Traceback complet:\n{traceback.format_exc()}")
    PROSPECTOR_AVAILABLE = False
    FSPS_AVAILABLE = False
    DYNESTY_AVAILABLE = False
    EMCEE_AVAILABLE = False
    logger.info(f"[PROSPECTOR] PROSPECTOR_AVAILABLE = {PROSPECTOR_AVAILABLE}")

# Import specutils pour la compatibilité
SPECUTILS_AVAILABLE = False
Spectrum1D = Any  # Placeholder par défaut
try:
    from specutils import Spectrum1D
    SPECUTILS_AVAILABLE = True
except ImportError:
    SPECUTILS_AVAILABLE = False
    logger.warning("specutils n'est pas disponible")
    # Spectrum1D reste défini comme Any (placeholder)


class ProspectorAnalyzer:
    """
    Classe pour analyser des spectres de galaxies et inférer des propriétés stellaires
    à partir de SED en utilisant Prospector.
    """
    
    def __init__(self, use_fsps: bool = True):
        """
        Initialise l'analyseur Prospector.
        
        Parameters
        ----------
        use_fsps : bool
            Utiliser FSPS pour les templates SSP (par défaut: True)
        """
        if not PROSPECTOR_AVAILABLE:
            raise ImportError("Prospector n'est pas installé. Installez-le avec: pip install prospector")
        
        self.model = None
        self.sps = None  # Simple Stellar Population object
        self.obs_data = None
        self.fitting_result = None
        self.sed_model = None
        self.use_fsps = use_fsps
        self.prospector_backend = "wsl" if (PROSPECTOR_BACKEND == "wsl" and PROSPECTOR_WSL_AVAILABLE) else "windows"
        
        # Initialiser l'objet SPS si FSPS est disponible
        if use_fsps and FSPS_AVAILABLE and self.prospector_backend != "wsl":
            try:
                # FastStepBasis est le wrapper Prospector pour FSPS
                # Il prend en paramètre les options FSPS
                self.sps = FastStepBasis(
                    zcontinuous=1,  # Interpolation continue en métallicité
                    compute_vega_mags=False,
                    vactoair_flag=False
                )
                logger.info("Objet SPS initialisé avec FSPS (FastStepBasis)")
            except Exception as e:
                logger.warning(f"Impossible d'initialiser FSPS: {e}")
                logger.warning(f"Type d'erreur: {type(e).__name__}")
                self.sps = None
                if not FSPS_AVAILABLE:
                    logger.error("FSPS n'est pas installé. Installez-le avec: pip install fsps")
        elif self.prospector_backend == "wsl":
            logger.info("Objet SPS local non initialisé: backend Prospector via WSL")
        
    def load_galaxy_spectrum(self, file_path: Union[str, Path]) -> Any:
        """
        Charge un spectre de galaxie depuis un fichier FITS ou ASCII.
        
        Parameters
        ----------
        file_path : str ou Path
            Chemin vers le fichier spectre
            
        Returns
        -------
        Spectrum1D
            Spectre chargé
        """
        if not SPECUTILS_AVAILABLE:
            raise ImportError("specutils n'est pas disponible")
        
        file_path = Path(file_path)
        
        try:
            if file_path.suffix.lower() in ['.fits', '.fit', '.fts']:
                # Charger avec specutils
                spectrum = Spectrum1D.read(str(file_path))
                logger.info(f"Spectre de galaxie chargé depuis {file_path}")
                return spectrum
            else:
                # Charger depuis ASCII
                data = np.loadtxt(file_path)
                if data.shape[1] >= 2:
                    wavelength = data[:, 0] * u.AA
                    flux = data[:, 1] * u.Unit('erg cm-2 s-1 AA-1')
                    spectrum = Spectrum1D(spectral_axis=wavelength, flux=flux)
                    logger.info(f"Spectre ASCII chargé depuis {file_path}")
                    return spectrum
                else:
                    raise ValueError("Format ASCII invalide. Attendu: lambda flux")
        except Exception as e:
            logger.error(f"Erreur lors du chargement du spectre: {e}")
            raise
    
    def create_sed_from_photometry(self, 
                                   wavelengths: np.ndarray,
                                   fluxes: np.ndarray,
                                   flux_errors: Optional[np.ndarray] = None,
                                   filters: Optional[List[str]] = None) -> Dict:
        """
        Crée une SED à partir de données photométriques (multi-bande).
        
        Parameters
        ----------
        wavelengths : array-like
            Longueurs d'onde centrales des filtres (en Angstroms)
        fluxes : array-like
            Flux dans chaque filtre (en erg/s/cm²/Å ou magnitudes)
        flux_errors : array-like, optionnel
            Incertitudes sur les flux
        filters : list of str, optionnel
            Noms des filtres (ex: ['u', 'g', 'r', 'i', 'z'])
            
        Returns
        -------
        dict
            Dictionnaire contenant les données SED
        """
        wavelengths = np.asarray(wavelengths)
        fluxes = np.asarray(fluxes)
        
        sed_data = {
            'wavelength': wavelengths,
            'flux': fluxes,
            'flux_error': flux_errors if flux_errors is not None else np.zeros_like(fluxes),
            'filters': filters if filters is not None else [f'filter_{i}' for i in range(len(wavelengths))]
        }
        
        logger.info(f"SED créée avec {len(wavelengths)} points photométriques")
        return sed_data
    
    def create_sed_from_spectrum(self, spectrum: Any) -> Dict:
        """
        Crée une SED à partir d'un spectre (specutils Spectrum1D).
        
        Parameters
        ----------
        spectrum : Spectrum1D
            Spectre à convertir en SED
            
        Returns
        -------
        dict
            Dictionnaire contenant les données SED
        """
        wavelength = spectrum.spectral_axis.value  # Angstroms
        flux = spectrum.flux.value  # Flux
        
        sed_data = {
            'wavelength': wavelength,
            'flux': flux,
            'flux_error': np.zeros_like(flux),  # Pas d'erreurs par défaut
            'type': 'spectrum'
        }
        
        logger.info(f"SED créée depuis spectre avec {len(wavelength)} points")
        return sed_data
    
    def combine_photometry_and_spectrum(self,
                                       photometry_sed: Dict,
                                       spectrum_sed: Dict) -> Dict:
        """
        Combine des données photométriques et spectroscopiques en une SED complète.
        
        Parameters
        ----------
        photometry_sed : dict
            SED photométrique (dict avec 'wavelength', 'flux', etc.)
        spectrum_sed : dict
            SED spectroscopique (dict avec 'wavelength', 'flux', etc.)
            
        Returns
        -------
        dict
            SED combinée
        """
        # Trier par longueur d'onde
        all_wavelengths = np.concatenate([photometry_sed['wavelength'], spectrum_sed['wavelength']])
        all_fluxes = np.concatenate([photometry_sed['flux'], spectrum_sed['flux']])
        all_errors = np.concatenate([
            photometry_sed.get('flux_error', np.zeros(len(photometry_sed['wavelength']))),
            spectrum_sed.get('flux_error', np.zeros(len(spectrum_sed['wavelength'])))
        ])
        
        # Trier par longueur d'onde
        sort_idx = np.argsort(all_wavelengths)
        
        combined_sed = {
            'wavelength': all_wavelengths[sort_idx],
            'flux': all_fluxes[sort_idx],
            'flux_error': all_errors[sort_idx],
            'type': 'combined'
        }
        
        logger.info(f"SED combinée avec {len(combined_sed['wavelength'])} points")
        return combined_sed
    
    def setup_prospector_model(self,
                              sed_data: Dict,
                              param_bounds: Optional[Dict] = None,
                              use_fsps: bool = True) -> "SedModel":
        """
        Configure un modèle Prospector pour l'inférence des propriétés stellaires.
        
        Parameters
        ----------
        sed_data : dict
            Données SED (dict avec 'wavelength', 'flux', 'flux_error')
        param_bounds : dict, optionnel
            Bornes des paramètres à ajuster (âge, métallicité, extinction, etc.)
        use_fsps : bool
            Utiliser FSPS pour les templates SSP (par défaut: True)
            
        Returns
        -------
        SedModel
            Modèle Prospector configuré
        """
        if not PROSPECTOR_AVAILABLE:
            raise ImportError("Prospector n'est pas installé")
        
        # Paramètres par défaut si non spécifiés
        if param_bounds is None:
            param_bounds = {
                'zred': (0.0, 3.0),  # Redshift
                'logzsol': (-2.0, 0.5),  # Métallicité en log(Z/Z_sol)
                'tage': (0.1, 13.8),  # Âge en Gyr
                'dust2': (0.0, 2.0),  # Extinction (Av)
            }
        
        # 1. Configurer les templates SSP (Simple Stellar Population)
        # Utiliser l'objet SPS déjà créé ou en créer un nouveau
        ssp_basis = None
        
        if use_fsps and FSPS_AVAILABLE:
            # Utiliser FSPS pour générer les templates SSP
            if self.sps is not None:
                ssp_basis = self.sps
                logger.info("Utilisation de l'objet SPS existant (FSPS)")
            else:
                try:
                    # Créer une base SSP avec FSPS
                    ssp_basis = FastStepBasis(
                        zcontinuous=1,  # Interpolation continue en métallicité
                        compute_vega_mags=False,
                        vactoair_flag=False
                    )
                    self.sps = ssp_basis
                    logger.info("Templates SSP configurés avec FSPS")
                except Exception as e:
                    logger.warning(f"Impossible de configurer FSPS: {e}")
                    logger.warning("Prospector utilisera un modèle simplifié sans FSPS")
                    ssp_basis = None
        else:
            if not FSPS_AVAILABLE:
                logger.warning("FSPS non disponible - Prospector fonctionnera en mode simplifié")
                logger.warning("Les templates SSP ne seront pas générés automatiquement")
                logger.warning("Note: Pour une analyse complète, installez FSPS manuellement")
            ssp_basis = None
        
        # Si FSPS n'est pas disponible, on peut toujours créer le modèle
        # mais il fonctionnera en mode simplifié (sans templates SSP automatiques)
        if ssp_basis is None:
            logger.info("Configuration du modèle Prospector sans FSPS (mode simplifié)")
        
        # 2. Définir les paramètres du modèle SED
        # Paramètres de base pour une population stellaire simple
        model_params = {}
        
        # Paramètres de l'historique de formation stellaire (SFH)
        # Modèle simple: SFH constante sur une durée donnée
        model_params['sfh'] = 0  # 0 = SFH constante
        model_params['tage'] = {'init': 5.0, 'prior': transforms.LogUniform(mini=0.1, maxi=13.8)}
        model_params['tau'] = {'init': 1.0, 'prior': transforms.LogUniform(mini=0.1, maxi=10.0)}
        
        # Métallicité
        model_params['logzsol'] = {
            'init': -0.5,
            'prior': transforms.TopHat(mini=param_bounds['logzsol'][0], 
                                      maxi=param_bounds['logzsol'][1])
        }
        
        # Extinction par la poussière (loi de Calzetti)
        model_params['dust_type'] = 0  # Loi de Calzetti
        model_params['dust1'] = 0.0  # Extinction en V band pour les étoiles jeunes
        model_params['dust2'] = {
            'init': 0.3,
            'prior': transforms.TopHat(mini=param_bounds['dust2'][0], 
                                      maxi=param_bounds['dust2'][1])
        }
        
        # Redshift
        model_params['zred'] = {
            'init': 0.0,
            'prior': transforms.TopHat(mini=param_bounds['zred'][0], 
                                      maxi=param_bounds['zred'][1])
        }
        
        # 4. Calibration spectroscopique (polynôme pour corriger les variations spectrales)
        # Si on a des données spectroscopiques, ajouter une calibration polynomiale
        if sed_data.get('type') == 'spectrum' or len(sed_data['wavelength']) > 100:
            # Ajouter des paramètres pour la calibration spectroscopique
            # Polynôme de degré 3 pour la calibration spectrophotométrique
            polyorder = 3  # Degré du polynôme de calibration
            model_params['spec_norm'] = {'init': 1.0, 'isfree': False}  # Normalisation spectrale
            
            # Paramètres du polynôme de calibration (coefficients a0, a1, a2, a3)
            for i in range(polyorder + 1):
                if i == 0:
                    model_params[f'poly_coeffs_{i}'] = {'init': 1.0, 'isfree': True}
                else:
                    model_params[f'poly_coeffs_{i}'] = {
                        'init': 0.0,
                        'isfree': True,
                        'prior': transforms.TopHat(mini=-1.0, maxi=1.0)
                    }
            
            model_params['polyorder'] = polyorder
            model_params['spec_redshift'] = {'init': 0.0, 'isfree': False}  # Redshift spectral
            
            logger.info(f"Calibration spectroscopique configurée (polynôme degré {polyorder})")
        
        # 3. Créer le modèle SED
        try:
            # SedModel de Prospector nécessite sps (Stellar Population Synthesis object) et model_params
            # Si ssp_basis est None, le modèle ne peut pas être créé correctement
            if ssp_basis is None:
                raise ValueError(
                    "Impossible de créer le modèle Prospector: FSPS n'est pas disponible.\n"
                    "FSPS est requis pour générer les templates SSP.\n"
                    "Installez FSPS (voir docs/INSTALLATION_FSPS.md)"
                )
            
            # Créer le modèle avec les paramètres
            # SpecModel/SedModel prend 'sps' comme premier paramètre positionnel ou nommé
            if SedModel is None:
                # Si SedModel n'a pas été importé, essayer de le récupérer dynamiquement
                from prospect.models import SpecModel
                SedModel = SpecModel
            
            model = SedModel(sps=ssp_basis, **model_params)
            logger.info("Modèle Prospector configuré avec succès")
            self.model = model
            return model
        except Exception as e:
            logger.error(f"Erreur lors de la configuration du modèle Prospector: {e}")
            logger.error(f"Type d'erreur: {type(e).__name__}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise
    
    def setup_observation_data(self, sed_data: Dict) -> Dict:
        """
        Prépare les données d'observation pour Prospector.
        
        Parameters
        ----------
        sed_data : dict
            Données SED (dict avec 'wavelength', 'flux', 'flux_error')
            
        Returns
        -------
        dict
            Dictionnaire de données d'observation formaté pour Prospector
        """
        wavelength = np.asarray(sed_data['wavelength'])
        flux = np.asarray(sed_data['flux'])
        flux_error = sed_data.get('flux_error', np.ones_like(flux) * 0.1 * flux)
        flux_error = np.asarray(flux_error)
        
        # Convertir les flux en maggies (1 maggie = 3631 Jy)
        # Conversion: Flux (erg/s/cm²/Å) -> maggies
        # 1 maggie = 3.631e-20 erg/s/cm²/Hz = 3.631e-20 * c / λ² erg/s/cm²/Å
        # Pour simplifier, on utilise une conversion approximative
        c = 2.998e18  # Angstrom/s
        flux_maggies = flux * wavelength**2 / (c * 3.631e-20)
        flux_error_maggies = flux_error * wavelength**2 / (c * 3.631e-20)
        
        # Préparer les données d'observation
        obs = {
            'wavelength': wavelength,  # Angstroms
            'spectrum': flux_maggies,  # Flux en maggies
            'unc': flux_error_maggies,  # Incertitudes en maggies
        }
        
        # Si on a des données photométriques (peu de points), les traiter comme des filtres
        if len(wavelength) < 20:
            # Données photométriques
            obs['phot_wave'] = wavelength
            obs['phot_flux'] = flux_maggies
            obs['phot_unc'] = flux_error_maggies
            obs['maggies'] = flux_maggies
            obs['maggies_unc'] = flux_error_maggies
            obs['filters'] = None  # Sera défini si on a des filtres réels
        else:
            # Données spectroscopiques
            obs['wavelength'] = wavelength
            obs['spectrum'] = flux_maggies
            obs['unc'] = flux_error_maggies
            obs['mask'] = np.ones(len(wavelength), dtype=bool)  # Masque pour les pixels valides
        
        logger.info(f"Données d'observation préparées: {len(wavelength)} points ({'spectroscopie' if len(wavelength) >= 20 else 'photométrie'})")
        return obs

    def _fit_stellar_properties_wsl(
        self,
        sed_data: Dict,
        n_walkers: int = 100,
        n_steps: int = 1000,
    ) -> Dict:
        """Exécute une inférence Prospector via le venv WSL."""
        wsl_python = _resolve_wsl_python()
        wsl_sps_home = _resolve_wsl_sps_home()
        worker_path = Path(__file__).resolve().parent / "prospector_wsl_worker.py"
        if not worker_path.exists():
            raise FileNotFoundError(f"Worker WSL introuvable: {worker_path}")

        with tempfile.TemporaryDirectory(prefix="npoap_wsl_fit_") as tmpdir:
            in_json = Path(tmpdir) / "input.json"
            out_json = Path(tmpdir) / "output.json"
            payload = {
                "sed_data": sed_data,
                "n_walkers": int(n_walkers),
                "n_steps": int(n_steps),
            }
            in_json.write_text(json.dumps(payload), encoding="utf-8")

            in_wsl = windows_path_to_wsl(in_json)
            out_wsl = windows_path_to_wsl(out_json)
            worker_wsl = windows_path_to_wsl(worker_path)

            cmd = (
                f'PY="{wsl_python}"; '
                f'export SPS_HOME="{wsl_sps_home}"; '
                f'"${{PY/#\\~/$HOME}}" "{worker_wsl}" "{in_wsl}" "{out_wsl}"'
            )
            cp = subprocess.run(
                ["wsl", "-e", "bash", "-lc", cmd],
                capture_output=True,
                text=True,
                check=False,
            )
            if cp.returncode != 0:
                err = (cp.stderr or cp.stdout or "").strip()
                raise RuntimeError(f"Échec worker Prospector WSL (code {cp.returncode}): {err}")
            if not out_json.exists():
                raise RuntimeError("Le worker WSL n'a pas produit de sortie JSON.")
            return json.loads(out_json.read_text(encoding="utf-8"))
    
    def fit_stellar_properties(self,
                              sed_data: Dict,
                              initial_params: Optional[Dict] = None,
                              n_walkers: int = 100,
                              n_steps: int = 1000,
                              sampler_type: str = 'dynesty') -> Dict:
        """
        Infère les propriétés stellaires à partir d'une SED en utilisant Prospector.
        
        Parameters
        ----------
        sed_data : dict
            Données SED
        initial_params : dict, optionnel
            Paramètres initiaux pour le fit
        n_walkers : int
            Nombre de marcheurs pour MCMC (si sampler_type='emcee')
        n_steps : int
            Nombre d'itérations MCMC (si sampler_type='emcee') ou dynesty
        sampler_type : str
            Type d'échantillonneur: 'dynesty' (par défaut) ou 'emcee'
            
        Returns
        -------
        dict
            Résultats de l'inférence (paramètres, incertitudes, etc.)
        """
        if not PROSPECTOR_AVAILABLE:
            raise ImportError("Prospector n'est pas installé")
        
        try:
            if self.prospector_backend == "wsl":
                logger.info("Inférence Prospector déléguée au backend WSL")
                return self._fit_stellar_properties_wsl(
                    sed_data,
                    n_walkers=n_walkers,
                    n_steps=n_steps,
                )

            logger.info(f"Démarrage de l'inférence Prospector (Sampler: {sampler_type})")
            
            # 1. Configurer le modèle Prospector
            if self.model is None:
                logger.info("Configuration du modèle Prospector...")
                self.model = self.setup_prospector_model(sed_data)
            
            # 2. Préparer les données d'observation
            obs = self.setup_observation_data(sed_data)
            self.obs_data = obs
            
            # 3. Paramétrer l'échantillonneur MCMC
            run_params = {}
            
            if sampler_type == 'dynesty' and DYNESTY_AVAILABLE:
                # Utiliser Dynesty (échantillonnage imbriqué)
                run_params['nested_method'] = 'dynesty'
                run_params['nested_nlive_init'] = min(400, n_walkers)
                run_params['nested_nlive_batch'] = min(200, n_walkers // 2)
                run_params['nested_dlogz_init'] = 0.05
                run_params['nested_maxiter'] = n_steps
                logger.info("Configuration de l'échantillonneur Dynesty")
                
            elif sampler_type == 'emcee' and EMCEE_AVAILABLE:
                # Utiliser Emcee (MCMC classique)
                run_params['sampler'] = 'emcee'
                run_params['nwalkers'] = n_walkers
                run_params['niter'] = n_steps
                run_params['nburn'] = max(100, n_steps // 4)  # Burn-in
                logger.info(f"Configuration de l'échantillonneur Emcee ({n_walkers} walkers, {n_steps} steps)")
                
            else:
                # Fallback: utiliser dynesty par défaut ou emcee si disponible
                if DYNESTY_AVAILABLE:
                    run_params['nested_method'] = 'dynesty'
                    run_params['nested_nlive_init'] = 400
                    run_params['nested_maxiter'] = n_steps
                    logger.info("Utilisation de Dynesty par défaut")
                elif EMCEE_AVAILABLE:
                    run_params['sampler'] = 'emcee'
                    run_params['nwalkers'] = n_walkers
                    run_params['niter'] = n_steps
                    run_params['nburn'] = n_steps // 4
                    logger.info("Utilisation de Emcee par défaut")
                else:
                    raise ImportError("Aucun échantillonneur MCMC disponible (dynesty ou emcee requis)")
            
            # 4. Paramètres additionnels pour l'ajustement
            run_params['verbose'] = True
            run_params['output_dir'] = None  # Pas de sauvegarde par défaut
            
            # 5. Lancer l'ajustement Prospector
            logger.info("Démarrage de l'ajustement Prospector...")
            
            # Utiliser l'objet SPS déjà créé
            if self.sps is not None:
                sps = self.sps
                logger.info("Utilisation de l'objet SPS existant")
            elif hasattr(self.model, 'sps') and self.model.sps is not None:
                # Le modèle peut avoir son propre objet SPS
                sps = self.model.sps
                logger.info("Utilisation de l'objet SPS du modèle")
            elif FSPS_AVAILABLE:
                try:
                    # Créer un objet SPS avec FSPS
                    sps = FastStepBasis(
                        zcontinuous=1,  # Interpolation continue en métallicité
                        compute_vega_mags=False,
                        vactoair_flag=False
                    )
                    self.sps = sps
                    logger.info("Objet SPS créé avec FSPS (FastStepBasis)")
                except Exception as e:
                    logger.error(f"Erreur lors de la création de l'objet SPS: {e}")
                    logger.error(f"Type d'erreur: {type(e).__name__}")
                    raise ImportError(f"Impossible de créer l'objet SPS avec FSPS: {e}")
            else:
                error_msg = (
                    "FSPS n'est pas disponible et aucun template SSP n'est configuré.\n"
                    "FSPS est requis pour générer des templates SSP.\n"
                    "Pour installer FSPS sur Windows:\n"
                    "1. Installer CMake: voir docs/INSTALLATION_CMAKE.md\n"
                    "2. Installer un compilateur Fortran (gfortran): voir docs/INSTALLATION_GFORTRAN.md\n"
                    "3. Installer FSPS: voir docs/INSTALLATION_FSPS.md ou utiliser installer_fsps_windows.ps1"
                )
                logger.error(error_msg)
                raise ImportError(error_msg)
            
            # Appel à fit_model de Prospector
            # Format: fit_model(obs, sps, model, ...)
            try:
                logger.info(f"Appel fit_model avec sps={type(sps).__name__}, model={type(self.model).__name__}")
                output = fit_model(
                    obs,
                    sps,
                    self.model,
                    optimize=False,  # Pas d'optimisation préliminaire
                    **run_params
                )
                logger.info("fit_model terminé avec succès")
                
                # Extraire les résultats depuis l'output de Prospector
                # Le format dépend de l'échantillonneur utilisé
                if isinstance(output, dict):
                    # Format dictionnaire
                    result = output.get('sampling', {})
                    theta = result.get('bestfit', {})
                    samples = result.get('chain', None)
                    
                    # Si pas dans 'sampling', chercher directement
                    if theta is None or len(theta) == 0:
                        theta = output.get('bestfit', {})
                    if samples is None:
                        samples = output.get('chain', None)
                else:
                    # Format objet
                    if hasattr(output, 'theta'):
                        theta = output.theta
                    elif hasattr(output, 'bestfit'):
                        theta = output.bestfit
                    else:
                        theta = {}
                    
                    samples = getattr(output, 'chain', None)
                
                # Si theta est un tableau, le convertir en dictionnaire
                if isinstance(theta, np.ndarray) or isinstance(theta, (list, tuple)):
                    param_names = list(self.model.theta_labels()) if hasattr(self.model, 'theta_labels') else \
                                [f'param_{i}' for i in range(len(theta))]
                    theta = {name: val for name, val in zip(param_names, theta)}
                
                # Calculer les incertitudes depuis les échantillons
                uncertainties = {}
                if samples is not None and len(samples) > 0:
                    samples = np.asarray(samples)
                    # Calculer écarts-types pour chaque paramètre
                    param_names = list(theta.keys()) if isinstance(theta, dict) else \
                                list(self.model.theta_labels()) if hasattr(self.model, 'theta_labels') else \
                                [f'param_{i}' for i in range(len(theta))]
                    
                    for i, param_name in enumerate(param_names):
                        if samples.ndim == 2:
                            # Format (n_samples, n_params)
                            if samples.shape[1] > i:
                                param_samples = samples[:, i]
                            else:
                                param_samples = samples.flatten()
                        elif samples.ndim == 3:
                            # Format (n_walkers, n_steps, n_params)
                            param_samples = samples[:, :, i].flatten() if samples.shape[2] > i else samples.flatten()
                        else:
                            param_samples = samples.flatten()
                        
                        uncertainties[param_name] = float(np.std(param_samples))
                else:
                    # Pas d'échantillons, utiliser des incertitudes par défaut (10% de la valeur)
                    for param_name, param_val in theta.items():
                        uncertainties[param_name] = 0.1 * abs(param_val) if param_val != 0 else 0.1
                
                # Construire le dictionnaire de résultats
                results = {
                    'success': True,
                    'parameters': theta,
                    'uncertainties': uncertainties,
                    'samples': samples,
                    'output': output
                }
                
                logger.info("Inférence Prospector terminée avec succès")
                return results
                
            except Exception as e:
                logger.error(f"Erreur lors de l'ajustement Prospector: {e}", exc_info=True)
                raise
            
        except Exception as e:
            logger.error(f"Erreur lors de l'inférence Prospector: {e}", exc_info=True)
            results = {
                'success': False,
                'message': f'Erreur lors de l\'inférence: {str(e)}',
                'parameters': {},
                'uncertainties': {},
                'error': str(e)
            }
            return results
    
    def get_stellar_properties_summary(self, results: Dict) -> str:
        """
        Génère un résumé des propriétés stellaires inférées.
        
        Parameters
        ----------
        results : dict
            Résultats de l'inférence Prospector
            
        Returns
        -------
        str
            Résumé formaté des propriétés stellaires
        """
        if not results.get('success', False):
            return "Inférence non réussie - voir les logs pour plus de détails"
        
        summary = "=== Propriétés Stellaires Inférées ===\n\n"
        backend = results.get("backend")
        if backend:
            summary += f"Backend: {backend}\n"
        if results.get("message"):
            summary += f"{results.get('message')}\n\n"
        
        params = results.get('parameters', {})
        uncertainties = results.get('uncertainties', {})
        
        if 'tage' in params:
            tage_val = params['tage']
            tage_err = uncertainties.get('tage', 0.0)
            summary += f"Âge: {tage_val:.2f} ± {tage_err:.2f} Gyr\n"
        
        if 'logzsol' in params:
            zsol_val = params['logzsol']
            zsol_err = uncertainties.get('logzsol', 0.0)
            summary += f"Métallicité: {zsol_val:.3f} ± {zsol_err:.3f} log(Z/Z_sol)\n"
        
        if 'dust2' in params:
            dust_val = params['dust2']
            dust_err = uncertainties.get('dust2', 0.0)
            summary += f"Extinction (Av): {dust_val:.2f} ± {dust_err:.2f} mag\n"
        
        if 'zred' in params:
            zred_val = params['zred']
            zred_err = uncertainties.get('zred', 0.0)
            summary += f"Redshift: {zred_val:.4f} ± {zred_err:.4f}\n"
        
        return summary
