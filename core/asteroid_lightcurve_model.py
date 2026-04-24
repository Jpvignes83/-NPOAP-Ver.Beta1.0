# core/asteroid_lightcurve_model.py
"""
Mod�le de courbe de lumi�re pour ast�ro�des (rotation) et inversion (ajustement des param�tres).

Mod�le : flux = F0 * (1 + A * cos(2*pi*(t - t0) / P))
  - P : p�riode de rotation (jours)
  - t0 : �poque du maximum de lumi�re (JD)
  - A : amplitude (sans dimension, 0 < A < 1)
  - F0 : flux moyen / niveau de r�f�rence
"""
import logging
import numpy as np
from scipy.optimize import minimize

logger = logging.getLogger(__name__)


def fit_harmonic_series_wls(
    time,
    flux,
    flux_err,
    P,
    n_harmonics=2,
    linear_drift=False,
    t_ref=None,
    ridge_lam=None,
):
    """
    Moindres carrés (pondérés) pour P fixée avec plusieurs harmoniques :

        flux ≈ b0 + Σ_k [ a_k cos(k ω (t-t_ref)) + b_k sin(k ω (t-t_ref)) ]  (+ dérive linéaire optionnelle)

    Une seule harmonique (k=1) avec une fenêtre courte devant une longue P est souvent **illisible** :
    cos/sin varient presque linéairement → modèle trop rigide. Augmenter ``n_harmonics`` (souvent 2–3)
    améliore l’ajustement quand la forme n’est pas une sinusoïde pure.

    Parameters
    ----------
    n_harmonics : int
        Nombre d'harmoniques (≥ 1). Réduit automatiquement si n_points insuffisants.
    ridge_lam : float, optional
        Régularisation léger sur les coefficients non constants (stabilise si mal conditionné).
    """
    t = np.asarray(time, dtype=float)
    y = np.asarray(flux, dtype=float)
    n = len(t)
    if n < 3:
        return {
            "y_model": None,
            "beta": None,
            "success": False,
            "message": "Au moins 3 points requis.",
            "chi2": np.nan,
            "n_dof": 0,
            "phase_frac": np.nan,
            "n_harmonics": 0,
        }
    if P is None or not np.isfinite(P) or P <= 0:
        return {
            "y_model": None,
            "beta": None,
            "success": False,
            "message": "P invalide.",
            "chi2": np.nan,
            "n_dof": 0,
            "phase_frac": np.nan,
            "n_harmonics": 0,
        }

    span = float(np.ptp(t))
    phase_frac = span / float(P) if span > 0 else 0.0

    if t_ref is None:
        t_ref = float(np.median(t))
    omega = 2.0 * np.pi / float(P)
    phi = omega * (t - t_ref)

    n_h = int(max(1, n_harmonics))
    drift_cols = 1 if linear_drift else 0
    max_h_from_n = max(1, (n - 1 - drift_cols) // 2)
    n_h = min(n_h, max_h_from_n, 8)

    cols = [np.ones(n, dtype=float)]
    for k in range(1, n_h + 1):
        cols.append(np.cos(k * phi))
        cols.append(np.sin(k * phi))
    if linear_drift:
        td = t - float(np.mean(t))
        scale = float(np.std(td))
        if scale < 1e-12:
            scale = 1.0
        cols.append(td / scale)
    X = np.column_stack(cols)
    p = X.shape[1]

    if flux_err is not None:
        sig = np.asarray(flux_err, dtype=float)
        good = (sig > 0) & np.isfinite(sig)
        if not np.any(good):
            sig = np.ones(n, dtype=float)
        else:
            med = float(np.nanmedian(sig[good]))
            if med <= 0:
                med = 1.0
            sig = np.where(good, sig, med)
        Xw = X / sig[:, np.newaxis]
        yw = y / sig
        G = Xw.T @ Xw
        rhs = Xw.T @ yw
    else:
        Xw = X
        yw = y
        G = X.T @ X
        rhs = X.T @ y
        sig = np.ones(n, dtype=float)

    if ridge_lam is None:
        cond = float(np.linalg.cond(G)) if p > 1 else 1.0
        ridge_lam = 0.0
        if cond > 120.0 or (phase_frac < 0.22 and n_h == 1):
            ridge_lam = 1e-5 * float(np.trace(G)) / max(p, 1)
    if ridge_lam and ridge_lam > 0:
        D = np.eye(p, dtype=float)
        D[0, 0] = 0.0
        G = G + ridge_lam * D
    try:
        beta = np.linalg.solve(G, rhs)
    except np.linalg.LinAlgError:
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)

    y_hat = X @ beta
    chi2 = float(np.sum(((y - y_hat) / sig) ** 2))
    n_dof = n - p

    b0 = float(beta[0])
    b1 = float(beta[1]) if p > 2 else 0.0
    b2 = float(beta[2]) if p > 2 else 0.0
    amplitude_R = float(np.hypot(b1, b2))
    phase_shift_alpha = float(np.arctan2(b2, b1))
    t0_equiv = t_ref + phase_shift_alpha / omega

    msg_parts = [f"LS Fourier, {n_h} harmonique(s)"]
    if phase_frac < 0.12:
        msg_parts.append(
            f"Attention : seulement {100.0 * phase_frac:.1f}% d'une période couverte — modèle fragile."
        )
    elif phase_frac < 0.25:
        msg_parts.append(
            f"Couverture courte ({100.0 * phase_frac:.1f}% de P) — préférer ≥2 harmoniques ou vérifier P."
        )

    return {
        "y_model": y_hat,
        "beta": beta,
        "t_ref": t_ref,
        "omega": omega,
        "P": float(P),
        "amplitude_R": amplitude_R,
        "offset_b0": b0,
        "t0_equiv": t0_equiv,
        "phase_shift_alpha": phase_shift_alpha,
        "chi2": chi2,
        "n_dof": n_dof,
        "success": True,
        "message": " · ".join(msg_parts),
        "linear_drift": bool(linear_drift),
        "phase_frac": float(phase_frac),
        "n_harmonics": n_h,
        "ridge_lam": float(ridge_lam) if ridge_lam else 0.0,
    }


def fit_first_harmonic_wls(time, flux, flux_err, P, linear_drift=False, t_ref=None):
    """Compatibilité : une seule harmonique (préférez ``fit_harmonic_series_wls`` avec n_harmonics≥2)."""
    return fit_harmonic_series_wls(
        time, flux, flux_err, P, n_harmonics=1, linear_drift=linear_drift, t_ref=t_ref
    )


def light_curve_model(t, P, t0, A, F0):
    """
    Mod�le de courbe de lumi�re (rotation, premi�re harmonique).

    Parameters
    ----------
    t : array-like
        Temps (JD).
    P : float
        P�riode (jours).
    t0 : float
        �poque du maximum (JD).
    A : float
        Amplitude (0 < A < 1).
    F0 : float
        Flux de r�f�rence (niveau moyen).

    Returns
    -------
    np.ndarray
        Flux mod�lis�.
    """
    t = np.asarray(t, dtype=float)
    if P <= 0:
        return np.full_like(t, F0)
    phase = 2.0 * np.pi * (t - t0) / P
    return F0 * (1.0 + A * np.cos(phase))


def fit_light_curve(
    time,
    flux,
    flux_err=None,
    P_init=None,
    t0_init=None,
    A_init=0.1,
    F0_init=None,
    bounds_P=(0.05, 100.0),
):
    """
    Ajuste les param�tres du mod�le (inversion) par minimisation du chi�.

    Parameters
    ----------
    time : array-like
        Temps (JD).
    flux : array-like
        Flux observ� (normalis� ou non).
    flux_err : array-like, optional
        Incertitudes sur le flux. Si None, chi� = somme des r�sidus�.
    P_init : float, optional
        P�riode initiale. Si None, utilise la m�diane des �carts entre points.
    t0_init : float, optional
        �poque initiale. Si None, utilise le temps du premier maximum approch�.
    A_init : float
        Amplitude initiale (d�faut 0.1).
    F0_init : float, optional
        Flux moyen initial. Si None, utilise np.median(flux).
    bounds_P : tuple
        (P_min, P_max) en jours.

    Returns
    -------
    dict
        Cl�s: P, t0, A, F0, chi2, n_dof, success, message.
    """
    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)
    if flux_err is not None:
        flux_err = np.asarray(flux_err, dtype=float)
        if np.any(flux_err <= 0):
            flux_err = None
    n = len(time)
    if n < 4:
        return {
            "P": np.nan,
            "t0": np.nan,
            "A": np.nan,
            "F0": np.nan,
            "chi2": np.nan,
            "n_dof": 0,
            "success": False,
            "message": "Moins de 4 points.",
        }

    if F0_init is None:
        F0_init = float(np.median(flux))
    if P_init is None or P_init <= 0:
        dt = np.diff(np.sort(time))
        dt = dt[dt > 0]
        P_init = float(np.median(dt)) * 2.0 if len(dt) else 0.1
        P_init = np.clip(P_init, bounds_P[0], bounds_P[1])
    if t0_init is None:
        t0_init = float(np.median(time))

    def chi2(x):
        P, t0, A, F0 = x[0], x[1], x[2], x[3]
        if P <= 0 or A < 0 or A > 1 or F0 <= 0:
            return 1e30
        fmod = light_curve_model(time, P, t0, A, F0)
        res = flux - fmod
        if flux_err is not None:
            res = res / flux_err
        return np.sum(res ** 2)

    x0 = [P_init, t0_init, A_init, F0_init]
    bnds = [
        (bounds_P[0], bounds_P[1]),
        (time.min() - 10 * P_init, time.max() + 10 * P_init),
        (0.001, 1.0),
        (1e-6, np.inf),
    ]
    try:
        res = minimize(
            chi2,
            x0,
            method="L-BFGS-B",
            bounds=bnds,
            options=dict(maxiter=500),
        )
        P, t0, A, F0 = res.x
        chi2_val = float(res.fun)
        n_dof = n - 4
        return {
            "P": P,
            "t0": t0,
            "A": A,
            "F0": F0,
            "chi2": chi2_val,
            "n_dof": n_dof,
            "success": res.success,
            "message": res.message,
        }
    except Exception as e:
        logger.exception("Erreur fit_light_curve: %s", e)
        return {
            "P": np.nan,
            "t0": np.nan,
            "A": np.nan,
            "F0": np.nan,
            "chi2": np.nan,
            "n_dof": n - 4,
            "success": False,
            "message": str(e),
        }
