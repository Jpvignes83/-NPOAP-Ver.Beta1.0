# utils/wsl_utils.py
"""
Utilitaires pour l'interaction avec WSL (Windows Subsystem for Linux).
Utilise par la detection KBMOD via WSL depuis NPOAP sous Windows.
"""
from pathlib import Path
from typing import Union


def windows_path_to_wsl(windows_path: Union[str, Path]) -> str:
    """
    Convertit un chemin Windows (ex. C:\\Users\\...) en chemin WSL (/mnt/c/Users/...).

    Args:
        windows_path: Chemin absolu ou relatif sous Windows.

    Returns:
        Chemin equivalent sous WSL, utilisable dans une commande `wsl ...`.
    """
    resolved = Path(windows_path).resolve()
    drive = resolved.drive.replace(":", "").lower()
    # as_posix() donne des / ; enlever le "C:" ou "c:" du debut
    posix_path = resolved.as_posix()
    if ":" in posix_path:
        posix_path = posix_path.split(":", 1)[1]
    if posix_path.startswith("/"):
        return f"/mnt/{drive}{posix_path}"
    return f"/mnt/{drive}/{posix_path}"
