#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
NPOAP - Script de vérification d'installation
Vérifie toutes les dépendances et modules requis pour NPOAP
"""

import sys
import os
import importlib.metadata as _importlib_metadata
from pathlib import Path


def _patch_distributions_skip_null_metadata():
    """
    Contourne un bug CuPy / importlib.metadata : certaines entrees de site-packages
    ont distribution.metadata is None, ce qui fait echouer cupy._detect_duplicate_installation
    puis tout import transitif (dask -> specutils).
    """
    _orig = _importlib_metadata.distributions

    def _filtered():
        for dist in _orig():
            if getattr(dist, "metadata", None) is not None:
                yield dist

    _importlib_metadata.distributions = _filtered  # type: ignore[assignment]


_patch_distributions_skip_null_metadata()

# Codes de couleur ANSI pour le terminal
class Colors:
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    UNDERLINE = '\033[4m'

def print_header(text):
    """Affiche un en-tête formaté"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text:^70}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.RESET}\n")

def print_success(text):
    """Affiche un message de succès"""
    print(f"{Colors.GREEN}✓{Colors.RESET} {text}")

def print_error(text):
    """Affiche un message d'erreur"""
    print(f"{Colors.RED}✗{Colors.RESET} {text}")

def print_warning(text):
    """Affiche un message d'avertissement"""
    print(f"{Colors.YELLOW}⚠{Colors.RESET} {text}")

def print_info(text):
    """Affiche un message d'information"""
    print(f"{Colors.CYAN}ℹ{Colors.RESET} {text}")

def check_module(module_name, required=True, version_attr='__version__', optional_reason=""):
    """
    Vérifie si un module est disponible
    
    Parameters:
    -----------
    module_name : str
        Nom du module à vérifier
    required : bool
        Si True, module est requis (erreur si absent). Si False, optionnel (avertissement si absent)
    version_attr : str
        Attribut contenant la version (par défaut '__version__')
    optional_reason : str
        Raison pour laquelle le module est optionnel (si required=False)
    
    Returns:
    --------
    tuple : (bool, str, str)
        (disponible, version, message)
    """
    try:
        module = __import__(module_name)
        version = getattr(module, version_attr, "version inconnue")
        return True, version, f"{module_name} {version}"
    except ImportError as e:
        if required:
            return False, None, f"{module_name} non installé (REQUIS)"
        else:
            reason = f" ({optional_reason})" if optional_reason else ""
            return False, None, f"{module_name} non installé{reason} (optionnel)"
    except Exception as e:
        hint = ""
        msg = str(e)
        if "'NoneType' object has no attribute 'get'" in msg or "_detect_duplicate_installation" in msg:
            hint = (
                " — indice: CuPy ou métadonnées pip incohérentes (pip check ; "
                "pip uninstall cupy cupy-cuda11x cupy-cuda12x … puis une seule variante si GPU requis)."
            )
        if required:
            return False, None, f"{module_name} import échoue (REQUIS): {type(e).__name__}: {e}{hint}"
        reason = f" ({optional_reason})" if optional_reason else ""
        return False, None, f"{module_name} import échoue{reason} (optionnel): {type(e).__name__}: {e}{hint}"

def check_prospector():
    """Vérifie l'installation complète de Prospector"""
    results = {
        'available': False,
        'version': None,
        'components': {},
        'errors': [],
        'warnings': []
    }
    
    # Test 1: prospect module (ImportError ou erreur si SPS_HOME / FSPS incomplet)
    try:
        import prospect
        results['available'] = True
        results['version'] = getattr(prospect, '__version__', 'version inconnue')
        results['components']['prospect'] = True
    except ImportError as e:
        results['components']['prospect'] = False
        results['errors'].append(f"prospect non disponible: {e}")
        return results
    except Exception as e:
        results['available'] = False
        results['components']['prospect'] = False
        results['errors'].append(
            f"prospect installe mais import echoue (souvent SPS_HOME ou FSPS manquant): {e}"
        )
        results['warnings'].append(
            "Definir SPS_HOME et donnees FSPS, ou reinstaller avec INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat"
        )
        return results
    
    # Test 2: sedpy.observate
    try:
        from sedpy import observate
        results['components']['sedpy.observate'] = True
    except ImportError as e:
        results['components']['sedpy.observate'] = False
        results['errors'].append(f"sedpy.observate non disponible: {e}")
        results['warnings'].append("sedpy doit être installé depuis GitHub, pas PyPI")
    
    # Test 3: SpecModel
    try:
        from prospect.models import SpecModel
        results['components']['SpecModel'] = True
    except ImportError:
        try:
            from prospect.models import SedModel
            results['components']['SedModel'] = True
            results['warnings'].append("SedModel disponible (version ancienne), utilisez SpecModel dans les versions récentes")
        except ImportError as e:
            results['components']['SpecModel'] = False
            results['components']['SedModel'] = False
            results['warnings'].append(f"SpecModel/SedModel non disponible: {e}")
    
    # Test 4: FastStepBasis
    try:
        from prospect.sources import FastStepBasis
        results['components']['FastStepBasis'] = True
    except ImportError as e:
        results['components']['FastStepBasis'] = False
        results['warnings'].append(f"FastStepBasis non disponible: {e} (FSPS peut être requis)")
    
    # Test 5: fit_model
    try:
        from prospect.fitting import fit_model
        results['components']['fit_model'] = True
    except ImportError as e:
        results['components']['fit_model'] = False
        results['warnings'].append(f"fit_model non disponible: {e}")
    
    # Test 6: FSPS (optionnel)
    try:
        import fsps
        results['components']['FSPS'] = True
        results['fsps_version'] = getattr(fsps, '__version__', 'version inconnue')
    except ImportError:
        results['components']['FSPS'] = False
        results['warnings'].append("FSPS non installé (Prospector utilisera des fichiers stub)")
    
    # Test 7: SPS_HOME
    sps_home = os.environ.get('SPS_HOME')
    if sps_home:
        results['SPS_HOME'] = sps_home
        # Vérifier que le fichier stub existe
        stub_file = Path(sps_home) / 'dust' / 'Nenkova08_y010_torusg_n10_q2.0.dat'
        if stub_file.exists():
            results['components']['FSPS_stub'] = True
            # Vérifier le format du fichier stub
            try:
                with open(stub_file, 'r') as f:
                    lines = f.readlines()
                if len(lines) >= 129:
                    # Vérifier le séparateur (3 espaces)
                    if len(lines) >= 5:
                        data_line = lines[4].strip()
                        if data_line and not data_line.startswith('#'):
                            cols = data_line.split('   ')  # 3 espaces
                            if len(cols) == 10:
                                results['components']['FSPS_stub_format'] = True
                            else:
                                results['warnings'].append("Fichier stub FSPS: format incorrect (nombre de colonnes ou séparateur)")
                                results['components']['FSPS_stub_format'] = False
                else:
                    results['warnings'].append(f"Fichier stub FSPS: nombre de lignes insuffisant ({len(lines)} < 129)")
                    results['components']['FSPS_stub_format'] = False
            except Exception as e:
                results['warnings'].append(f"Erreur lors de la vérification du fichier stub: {e}")
                results['components']['FSPS_stub_format'] = False
        else:
            results['components']['FSPS_stub'] = False
            results['warnings'].append(f"Fichier stub FSPS manquant: {stub_file}")
    else:
        results['SPS_HOME'] = None
        results['warnings'].append("SPS_HOME non défini (sera créé automatiquement si nécessaire)")
    
    return results

def check_local_modules():
    """Vérifie que les modules locaux de NPOAP sont disponibles"""
    results = {
        'core': False,
        'gui': False,
        'utils': False,
        'config': False,
        'errors': []
    }
    
    # Vérifier que nous sommes dans le répertoire NPOAP
    npoap_dir = Path(__file__).parent
    
    # Test 1: core
    try:
        import core
        results['core'] = True
    except ImportError as e:
        results['errors'].append(f"Module 'core' non disponible: {e}")
        if not (npoap_dir / 'core').exists():
            results['errors'].append(f"Répertoire 'core' introuvable dans {npoap_dir}")
    
    # Test 2: gui
    try:
        import gui
        results['gui'] = True
    except ImportError as e:
        results['errors'].append(f"Module 'gui' non disponible: {e}")
        if not (npoap_dir / 'gui').exists():
            results['errors'].append(f"Répertoire 'gui' introuvable dans {npoap_dir}")
    
    # Test 3: utils
    try:
        import utils
        results['utils'] = True
    except ImportError as e:
        results['errors'].append(f"Module 'utils' non disponible: {e}")
        if not (npoap_dir / 'utils').exists():
            results['errors'].append(f"Répertoire 'utils' introuvable dans {npoap_dir}")
    
    # Test 4: config.py
    config_file = npoap_dir / 'config.py'
    if config_file.exists():
        results['config'] = True
    else:
        results['errors'].append(f"Fichier 'config.py' introuvable dans {npoap_dir}")
    
    # Test 5: main.py
    main_file = npoap_dir / 'main.py'
    if not main_file.exists():
        results['errors'].append(f"Fichier 'main.py' introuvable dans {npoap_dir}")
    
    return results

def main():
    """Fonction principale de vérification"""
    
    # En-tête
    print_header("VÉRIFICATION DE L'INSTALLATION NPOAP")
    
    # Compteurs
    required_ok = 0
    required_total = 0
    optional_ok = 0
    optional_total = 0
    errors = []
    warnings = []
    
    # === VÉRIFICATION DES DÉPENDANCES CORE (REQUISES) ===
    print_header("1. Dépendances Core (Requis)")
    
    core_modules = [
        ('numpy', True),
        ('scipy', True),
        ('pandas', True),
        ('astropy', True),
        ('photutils', True),
        ('ccdproc', True),
        ('astroquery', True),
        ('matplotlib', True),
        ('requests', True),
        # Paquet PyPI "Pillow", module Python "PIL" — ne pas ajouter de 3e element (sinon import "Pillow" echoue)
        ('PIL', True),
        ('emcee', True),
        ('specutils', True),
        ('statsmodels', True),
    ]
    
    for module_info in core_modules:
        if len(module_info) == 2:
            module_name, required = module_info
            import_name = module_name
        else:
            module_name, required, import_name = module_info
        
        available, version, message = check_module(import_name, required=required)
        
        if required:
            required_total += 1
            if available:
                print_success(message)
                required_ok += 1
            else:
                print_error(message)
                errors.append(message)
        else:
            optional_total += 1
            if available:
                print_success(f"{message} (optionnel)")
                optional_ok += 1
            else:
                print_warning(message)
                warnings.append(message)
    
    # === VÉRIFICATION DES DÉPENDANCES OPTIONNELLES ===
    print_header("2. Dépendances Optionnelles")
    
    optional_modules = [
        ('phoebe', False, "Étoiles binaires (PHOEBE2)"),
        ('rebound', False, "Simulation N-body"),
        ('ultranest', False, "Nested sampling bayésien"),
        ('pylightcurve', False, "Modélisation de courbes de lumière d'exoplanètes"),
    ]
    
    for module_info in optional_modules:
        if len(module_info) == 2:
            module_name, optional_reason = module_info
            import_name = module_name
        else:
            module_name, _, optional_reason = module_info
            import_name = module_name
        
        available, version, message = check_module(import_name, required=False, optional_reason=optional_reason)
        optional_total += 1
        if available:
            print_success(f"{message} (optionnel)")
            optional_ok += 1
        else:
            print_warning(message)
            warnings.append(message)
    
    # CuPy (GPU) - vérification spéciale
    print_info("Vérification de CuPy (accélération GPU)...")
    try:
        import cupy as cp
        gpu_available = cp.cuda.is_available()
        if gpu_available:
            device = cp.cuda.Device(0)
            compute_capability = device.compute_capability
            print_success(f"CuPy {cp.__version__} installé - GPU disponible (Compute Capability: {compute_capability})")
            optional_ok += 1
        else:
            print_warning(f"CuPy {cp.__version__} installé - GPU non disponible (CUDA requis)")
            warnings.append("CuPy installé mais GPU non disponible")
        optional_total += 1
    except Exception as e:
        if isinstance(e, ImportError):
            print_warning("CuPy non installé (optionnel - accélération GPU)")
            warnings.append("CuPy non installé (optionnel)")
        else:
            msg = str(e)
            hint = ""
            if "'NoneType' object has no attribute 'get'" in msg or "_detect_duplicate_installation" in msg:
                hint = (
                    " — indice: métadonnées pip / plusieurs variantes CuPy ; "
                    "voir docs/MANUEL_INSTALLATION.md (Dépannage CuPy)."
                )
            print_warning(f"CuPy import ou GPU : {type(e).__name__}: {e}{hint}")
            warnings.append(f"CuPy: {e}")
        optional_total += 1
    
    # === VÉRIFICATION DE PROSPECTOR ===
    print_header("3. Prospector (Optionnel - Analyse Spectroscopique)")
    
    prospector_results = check_prospector()
    
    if prospector_results['available']:
        print_success(f"prospect {prospector_results['version']} disponible")
        optional_ok += 1
        
        # Afficher les composants
        print(f"  {Colors.CYAN}Composants Prospector:{Colors.RESET}")
        for component, available in prospector_results['components'].items():
            if available:
                print(f"    {Colors.GREEN}✓{Colors.RESET} {component}")
            else:
                print(f"    {Colors.YELLOW}⚠{Colors.RESET} {component} non disponible")
        
        # Afficher SPS_HOME
        if prospector_results['SPS_HOME']:
            print(f"  {Colors.CYAN}SPS_HOME:{Colors.RESET} {prospector_results['SPS_HOME']}")
        else:
            print_warning("  SPS_HOME non défini")
        
        # Afficher FSPS si installé
        if prospector_results['components'].get('FSPS'):
            print_success(f"  FSPS {prospector_results.get('fsps_version', 'version inconnue')} installé")
        else:
            print_warning("  FSPS non installé (utilisation de fichiers stub)")
    else:
        print_warning("Prospector indisponible ou incomplet (optionnel)")
        if prospector_results['errors']:
            for error in prospector_results['errors']:
                print_error(f"  {error}")
                errors.append(f"Prospector: {error}")
    
    optional_total += 1
    
    # Ajouter les avertissements de Prospector
    warnings.extend([f"Prospector: {w}" for w in prospector_results.get('warnings', [])])
    
    # === VÉRIFICATION DES MODULES LOCAUX ===
    print_header("4. Modules Locaux NPOAP")
    
    local_results = check_local_modules()
    
    if local_results['core']:
        print_success("Module 'core' disponible")
    else:
        print_error("Module 'core' non disponible")
        errors.extend([f"Modules locaux: {e}" for e in local_results['errors'] if 'core' in e.lower()])
    
    if local_results['gui']:
        print_success("Module 'gui' disponible")
    else:
        print_error("Module 'gui' non disponible")
        errors.extend([f"Modules locaux: {e}" for e in local_results['errors'] if 'gui' in e.lower()])
    
    if local_results['utils']:
        print_success("Module 'utils' disponible")
    else:
        print_error("Module 'utils' non disponible")
        errors.extend([f"Modules locaux: {e}" for e in local_results['errors'] if 'utils' in e.lower()])
    
    if local_results['config']:
        print_success("Fichier 'config.py' trouvé")
    else:
        print_error("Fichier 'config.py' introuvable")
        errors.extend([f"Modules locaux: {e}" for e in local_results['errors'] if 'config' in e.lower()])
    
    # === VÉRIFICATION DE L'ENVIRONNEMENT PYTHON ===
    print_header("5. Environnement Python")
    
    print_info(f"Version Python: {sys.version}")
    print_info(f"Exécutable Python: {sys.executable}")
    print_info(f"Plateforme: {sys.platform}")
    
    # Vérifier si nous sommes dans un environnement virtuel ou Conda
    _conda = bool(os.environ.get("CONDA_DEFAULT_ENV") or os.environ.get("CONDA_PREFIX"))
    _venv = hasattr(sys, "real_prefix") or (
        hasattr(sys, "base_prefix") and sys.base_prefix != sys.prefix
    )
    if _conda or _venv:
        if _conda and os.environ.get("CONDA_DEFAULT_ENV"):
            print_success(f"Conda detecte (env: {os.environ.get('CONDA_DEFAULT_ENV')})")
        else:
            print_success("Environnement virtuel/Conda detecte")
    else:
        print_warning("Aucun environnement virtuel/Conda detecte (utilisation de Python systeme)")
    
    # === RÉSUMÉ FINAL ===
    print_header("RÉSUMÉ")
    
    print(f"{Colors.BOLD}Dépendances Requises:{Colors.RESET}")
    print(f"  {Colors.GREEN if required_ok == required_total else Colors.RED}"
          f"{required_ok}/{required_total} installées{Colors.RESET}")
    
    print(f"\n{Colors.BOLD}Dépendances Optionnelles:{Colors.RESET}")
    print(f"  {Colors.CYAN}{optional_ok}/{optional_total} installées{Colors.RESET}")
    
    # Afficher les erreurs
    if errors:
        print(f"\n{Colors.BOLD}{Colors.RED}Erreurs ({len(errors)}):{Colors.RESET}")
        for i, error in enumerate(errors, 1):
            print(f"  {i}. {Colors.RED}{error}{Colors.RESET}")
    
    # Afficher les avertissements
    if warnings:
        print(f"\n{Colors.BOLD}{Colors.YELLOW}Avertissements ({len(warnings)}):{Colors.RESET}")
        for i, warning in enumerate(warnings, 1):
            print(f"  {i}. {Colors.YELLOW}{warning}{Colors.RESET}")
    
    # Conclusion
    print()
    if required_ok == required_total and not local_results['errors']:
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.RESET}")
        print(f"{Colors.BOLD}{Colors.GREEN}✓ Installation NPOAP vérifiée avec succès!{Colors.RESET}")
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.RESET}")
        if warnings:
            print(f"\n{Colors.YELLOW}Note: Certaines fonctionnalités optionnelles ne sont pas disponibles,{Colors.RESET}")
            print(f"{Colors.YELLOW}mais les fonctionnalités de base de NPOAP fonctionneront correctement.{Colors.RESET}")
        return 0
    else:
        print(f"{Colors.BOLD}{Colors.RED}{'='*70}{Colors.RESET}")
        print(f"{Colors.BOLD}{Colors.RED}✗ Installation incomplète - Veuillez corriger les erreurs ci-dessus{Colors.RESET}")
        print(f"{Colors.BOLD}{Colors.RED}{'='*70}{Colors.RESET}")
        if errors:
            print(f"\n{Colors.YELLOW}Pour installer les dépendances manquantes:{Colors.RESET}")
            print(f"  conda activate astroenv")
            print(f"  pip install -r requirements.txt")
        return 1

if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print(f"\n\n{Colors.YELLOW}Vérification interrompue par l'utilisateur{Colors.RESET}")
        sys.exit(130)
    except Exception as e:
        print(f"\n\n{Colors.RED}Erreur inattendue lors de la vérification: {e}{Colors.RESET}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
