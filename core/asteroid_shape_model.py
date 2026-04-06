# core/asteroid_shape_model.py
"""
Chargement de mod�les de forme 3D d'ast�ro�des pour visualisation.

Formats support�s :
- OBJ (.obj) : lignes "v x y z" et "f i j k" (triangles, indices 1-based).
- DAMIT shape.txt : nombre de sommets, puis liste (x y z), puis nombre de facettes,
  puis liste (i j k) en indices 1-based (style Kaasalainen / DAMIT).
"""
import logging
import numpy as np
from pathlib import Path

logger = logging.getLogger(__name__)


def load_shape_obj(path):
    """
    Charge un maillage depuis un fichier OBJ (vertices + faces triangulaires).

    Returns
    -------
    vertices : (N, 3) float
    faces : (M, 3) int, indices 0-based
    """
    path = Path(path)
    vertices = []
    faces = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "v" and len(parts) >= 4:
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif parts[0] == "f" and len(parts) >= 4:
                # f i j k ou f i// j// k// ou f i/j/n k/l/n ...
                idx = []
                for p in parts[1:4]:
                    i = p.split("/")[0].split("//")[0]
                    idx.append(int(i) - 1)  # 1-based -> 0-based
                faces.append(idx)
    if not vertices or not faces:
        raise ValueError(f"OBJ invalide ou vide : {path}")
    return np.array(vertices), np.array(faces, dtype=np.int32)


def load_shape_txt(path):
    """
    Charge un maillage au format DAMIT / Kaasalainen shape.txt :
    - Premi�re ligne : nombre de sommets (nv)
    - nv lignes : x y z (s�par�s par espaces)
    - Ligne suivante : nombre de facettes (nf)
    - nf lignes : i j k (indices 1-based)

    Returns
    -------
    vertices : (N, 3) float
    faces : (M, 3) int, indices 0-based
    """
    path = Path(path)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        raise ValueError(f"Fichier shape.txt vide : {path}")

    # Premi�re ligne = nombre de sommets
    try:
        nv = int(lines[0].split()[0])
    except (ValueError, IndexError):
        raise ValueError(f"shape.txt : premi�re ligne doit �tre le nombre de sommets : {path}")

    vertices = []
    for i in range(1, 1 + nv):
        if i >= len(lines):
            raise ValueError(f"shape.txt : attendu {nv} sommets, ligne {i} manquante")
        parts = lines[i].split()
        if len(parts) < 3:
            continue
        vertices.append([float(parts[0]), float(parts[1]), float(parts[2])])

    start_facets = 1 + nv
    if start_facets >= len(lines):
        raise ValueError(f"shape.txt : nombre de facettes manquant apr�s les sommets")
    try:
        nf = int(lines[start_facets].split()[0])
    except (ValueError, IndexError):
        raise ValueError(f"shape.txt : ligne {start_facets} doit �tre le nombre de facettes")

    faces = []
    for i in range(start_facets + 1, start_facets + 1 + nf):
        if i >= len(lines):
            break
        parts = lines[i].split()
        if len(parts) < 3:
            continue
        # 1-based -> 0-based
        faces.append([int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2]) - 1])

    if not vertices or not faces:
        raise ValueError(f"shape.txt : aucun sommet ou facette valide : {path}")
    return np.array(vertices), np.array(faces, dtype=np.int32)


def load_shape(path):
    """
    Charge un mod�le de forme depuis un fichier .obj ou shape.txt (DAMIT).

    Parameters
    ----------
    path : str or Path
        Chemin vers le fichier.

    Returns
    -------
    vertices : (N, 3) np.ndarray
    faces : (M, 3) np.ndarray, indices 0-based
    """
    path = Path(path)
    suf = path.suffix.lower()
    if suf == ".obj":
        return load_shape_obj(path)
    if suf == ".txt" or "shape" in path.name.lower():
        return load_shape_txt(path)
    # Essai par contenu
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline().strip()
    if first.startswith("v ") or first.startswith("v\t"):
        return load_shape_obj(path)
    try:
        int(first.split()[0])
        return load_shape_txt(path)
    except (ValueError, IndexError):
        raise ValueError(f"Format non reconnu pour {path} (attendu .obj ou shape.txt DAMIT)")
