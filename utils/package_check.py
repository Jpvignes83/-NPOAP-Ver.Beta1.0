# utils/package_check.py
"""
Utilitaires pour vérifier les mises à jour des packages requis.
"""
import logging
import subprocess
import sys
import json
import urllib.request
import urllib.error
from typing import Optional, Tuple

logger = logging.getLogger(__name__)

def get_installed_version(package_name: str) -> Optional[str]:
    """
    Récupère la version installée d'un package.
    
    Args:
        package_name: Nom du package (ex: 'pylightcurve')
    
    Returns:
        Version installée ou None si le package n'est pas installé
    """
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "show", package_name],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if line.startswith('Version:'):
                    return line.split(':', 1)[1].strip()
    except Exception as e:
        logger.debug(f"Erreur lors de la récupération de la version installée de {package_name}: {e}")
    return None

def get_latest_version_pypi(package_name: str) -> Optional[str]:
    """
    Récupère la dernière version disponible sur PyPI via l'API JSON.
    
    Args:
        package_name: Nom du package (ex: 'pylightcurve')
    
    Returns:
        Dernière version disponible ou None en cas d'erreur
    """
    try:
        url = f"https://pypi.org/pypi/{package_name}/json"
        with urllib.request.urlopen(url, timeout=5) as response:
            data = json.loads(response.read())
            return data.get('info', {}).get('version')
    except urllib.error.URLError as e:
        logger.debug(f"Erreur de connexion lors de la vérification de {package_name} sur PyPI: {e}")
    except Exception as e:
        logger.debug(f"Erreur lors de la récupération de la dernière version de {package_name}: {e}")
    return None

def compare_versions(v1: str, v2: str) -> int:
    """
    Compare deux versions (format semver ou similaire).
    Retourne: -1 si v1 < v2, 0 si v1 == v2, 1 si v1 > v2
    
    Args:
        v1: Première version
        v2: Deuxième version
    
    Returns:
        Résultat de la comparaison
    """
    try:
        # Essayer d'utiliser packaging si disponible
        from packaging import version
        v1_parsed = version.parse(v1)
        v2_parsed = version.parse(v2)
        if v1_parsed < v2_parsed:
            return -1
        elif v1_parsed > v2_parsed:
            return 1
        else:
            return 0
    except ImportError:
        # Fallback: comparaison simple par tuples de nombres
        def parse_version(v):
            parts = []
            for part in v.split('.'):
                try:
                    parts.append(int(part))
                except ValueError:
                    # Si ce n'est pas un nombre, ignorer
                    pass
            return tuple(parts)
        
        v1_parts = parse_version(v1)
        v2_parts = parse_version(v2)
        if v1_parts < v2_parts:
            return -1
        elif v1_parts > v2_parts:
            return 1
        else:
            return 0

def check_pylightcurve_update() -> Tuple[bool, str, Optional[str]]:
    """
    Vérifie si une mise à jour de pylightcurve est disponible.
    
    Returns:
        Tuple (update_available, installed_version, latest_version)
        - update_available: True si une mise à jour est disponible
        - installed_version: Version installée
        - latest_version: Dernière version disponible (None si erreur)
    """
    package_name = "pylightcurve"
    
    installed = get_installed_version(package_name)
    if installed is None:
        logger.warning(f"{package_name} n'est pas installé")
        return False, "non installé", None
    
    latest = get_latest_version_pypi(package_name)
    if latest is None:
        logger.debug(f"Impossible de récupérer la dernière version de {package_name} (pas de connexion internet ?)")
        return False, installed, None
    
    try:
        # Comparer les versions
        comparison = compare_versions(installed, latest)
        if comparison < 0:
            logger.info(f"Mise à jour disponible pour {package_name}: {installed} → {latest}")
            return True, installed, latest
        else:
            logger.debug(f"{package_name} est à jour (version {installed})")
            return False, installed, latest
    except Exception as e:
        logger.debug(f"Erreur lors de la comparaison des versions de {package_name}: {e}")
        return False, installed, latest

