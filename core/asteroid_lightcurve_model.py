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
        (1e-6, None),
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
