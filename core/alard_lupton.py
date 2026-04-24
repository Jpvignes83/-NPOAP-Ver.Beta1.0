# core/alard_lupton.py
"""
Soustraction d'images Alard-Lupton (1998) en pur Python.

Référence: Alard & Lupton, ApJ 503, 325 (1998), eq. (1)-(3).
On trouve un noyau de convolution K et un fond optionnel b tels que
  Ref(x,y) ⊗ Kernel(u,v) = I(x,y)  [eq. 1]
  avec fond différentiel: Ref ⊗ K = I + bg  [eq. 3]
donc en pratique on ajuste  I_science ≈ I_ref ⊗ K + b  en moindres carrés.
L'image soustraite est I_science - (I_ref ⊗ K + b).

Simplifications par rapport à l'article:
- Noyau constant (pas de variation spatiale du noyau, §3 de l'article).
- Base: Gaussiennes 2D pures (pas de termes polynomiaux u^dx v^dy de l'eq. 2).
- Fond: constant (pas de polynôme bg(x,y) complet).
- Pondération: optionnelle Poisson σ² ∝ I (article §2.3).

Dépendances: numpy, scipy (fftconvolve).
"""

import logging
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

try:
    from scipy.signal import fftconvolve
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    fftconvolve = None


def _gaussian_kernel_2d(half_size: int, sigma: float) -> np.ndarray:
    """Noyau 2D Gaussien centré, taille (2*half_size+1) x (2*half_size+1)."""
    x = np.arange(-half_size, half_size + 1, dtype=float)
    xx, yy = np.meshgrid(x, x)
    k = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
    k /= k.sum()
    return k


def alard_lupton_subtract(
    science: np.ndarray,
    reference: np.ndarray,
    kernel_half_size: int = 10,
    kernel_sigmas: Tuple[float, ...] = (0.5, 1.0, 2.0, 3.0),
    fit_background: bool = True,
    regularize: float = 1e-8,
    use_poisson_weights: bool = False,
    gain: float = 1.0,
) -> np.ndarray:
    """
    Soustraction optimale Alard-Lupton (noyau constant).

    Correspond à l'article Alard & Lupton (1998): eq. (1) Ref ⊗ K = I,
    avec colonnes C_i = Ref ⊗ B_i et système normal M a = V.
    Option Poisson: M_ij = ∫ C_i C_j/σ², V_i = ∫ I C_i/σ², σ² = gain * I (article §2.3).

    Parameters
    ----------
    science : np.ndarray
        Image science (2D).
    reference : np.ndarray
        Image de référence, même forme que science (déjà alignée).
    kernel_half_size : int
        Demi-taille du noyau en pixels (noyau (2*h+1) x (2*h+1)). Défaut 10.
    kernel_sigmas : tuple of float
        Sigmas des Gaussiennes de la base. Défaut (0.5, 1.0, 2.0, 3.0).
    fit_background : bool
        Inclure un terme de fond constant (équivalent eq. 3 avec bg constant).
    regularize : float
        Régularisation (diagonale de AtA) pour stabilité. Défaut 1e-8.
    use_poisson_weights : bool
        Si True, pondération 1/σ² avec σ² = gain * max(I, eps) comme dans l'article.
    gain : float
        Facteur gain (ADU/photon) pour les poids Poisson; utilisé si use_poisson_weights=True.

    Returns
    -------
    np.ndarray
        Image soustraite (science - modèle), même forme que science.
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("Alard-Lupton nécessite scipy (scipy.signal.fftconvolve).")

    if science.shape != reference.shape:
        raise ValueError("science et reference doivent avoir la même forme.")

    science = np.asarray(science, dtype=float)
    reference = np.asarray(reference, dtype=float)

    # Masque des pixels valides (on évite NaN/Inf)
    valid = np.isfinite(science) & np.isfinite(reference)
    if not np.any(valid):
        raise ValueError("Aucun pixel fini commun dans science et reference.")

    # Pour les convolutions, remplir les non-finies par 0 pour ne pas propager NaN
    R = np.where(np.isfinite(reference), reference, 0.0)
    I = np.where(np.isfinite(science), science, 0.0)

    # Base de noyaux Gaussiennes
    basis_kernels = [
        _gaussian_kernel_2d(kernel_half_size, s) for s in kernel_sigmas
    ]

    # Convolutions C_i = R ⊗ K_i (mode='same' => même taille que R)
    convolved = []
    for Ki in basis_kernels:
        Ci = fftconvolve(R, Ki, mode="same")
        convolved.append(Ci)

    n_basis = len(convolved)
    n_params = n_basis + (1 if fit_background else 0)

    # Poids Poisson (article §2.3): σ(x,y) = k * √I, poids w = 1/σ²
    flat_valid = valid.ravel()
    if use_poisson_weights:
        eps = 1e-6 * np.nanmax(I) if np.nanmax(I) > 0 else 1e-6
        var = gain * np.where(I.ravel() > 0, I.ravel(), eps)
        weight = np.where(var > 0, 1.0 / var, 0.0)
        weight = np.where(flat_valid, weight, 0.0)
    else:
        weight = np.where(flat_valid, 1.0, 0.0)

    # Équations normales (article: M_ij = ∫ C_i C_j/σ², V_i = ∫ I C_i/σ²)
    # Colonnes C_i = R ⊗ K_i, cible y = I. AtA[p,q] = sum_p (col_p * col_q * w_p), Aty[p] = sum_p (col_p * y_p * w_p)
    AtA = np.zeros((n_params, n_params), dtype=float)
    Aty = np.zeros(n_params, dtype=float)
    I_flat = I.ravel()

    for i in range(n_basis):
        ci = convolved[i].ravel()
        for j in range(i + 1):
            cj = convolved[j].ravel()
            AtA[i, j] = AtA[j, i] = np.dot(ci, cj * weight)
        Aty[i] = np.dot(ci, I_flat * weight)

    if fit_background:
        ones_flat = np.ones_like(I).ravel()
        for i in range(n_basis):
            AtA[i, n_basis] = AtA[n_basis, i] = np.dot(convolved[i].ravel(), ones_flat * weight)
        AtA[n_basis, n_basis] = np.sum(weight)
        Aty[n_basis] = np.dot(I_flat, weight)

    # Régularisation
    AtA += regularize * np.eye(n_params)

    try:
        c = np.linalg.solve(AtA, Aty)
    except np.linalg.LinAlgError as e:
        logger.warning(f"Alard-Lupton: solve échoué ({e}), retour soustraction simple.")
        out = science - reference
        out -= np.nanmedian(out)
        return out.astype(np.float64)

    # Modèle = sum c_i * C_i + c_bg
    model = np.zeros_like(I)
    for i in range(n_basis):
        model += c[i] * convolved[i]
    if fit_background:
        model += c[n_basis]

    subtracted = science - model
    # Recentrer la médiane à 0 (convention courante)
    med = np.nanmedian(subtracted)
    subtracted = subtracted - med

    logger.debug(
        f"Alard-Lupton: noyau {len(basis_kernels)} Gaussiennes (sigmas={kernel_sigmas}), "
        f"coeffs={c[:4].round(4).tolist()}{(' ...' if len(c) > 4 else '')}"
    )
    return subtracted.astype(np.float64)
