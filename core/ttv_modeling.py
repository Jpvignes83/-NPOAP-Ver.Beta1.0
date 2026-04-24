import numpy as np
import emcee


def _sanitize_yerr(yerr: np.ndarray, y: np.ndarray, y_amp: float, y_std: float) -> np.ndarray:
    """
    Évite sigma² = 0 (erreurs O-C nulles / manquantes dans le CSV) qui provoquent NaN dans la log-vraisemblance.
    """
    yerr = np.asarray(yerr, dtype=float)
    y = np.asarray(y, dtype=float)
    scale = float(np.nanstd(y)) if y.size > 1 else 0.0
    if not np.isfinite(scale) or scale <= 0:
        scale = max(float(y_amp), float(y_std), 1e-12)
    floor = max(1e-15, 1e-10 * scale)
    yerr = np.where(np.isfinite(yerr), yerr, floor)
    return np.maximum(yerr, floor)


def multi_sine_model(t, *params):
    """
    Modèle multi-fréquences : Somme de A_i * sin(2*pi*t/P_i + phi_i) + Offset
    """
    offset = params[-1]
    sine_params = params[:-1]
    n_sines = len(sine_params) // 3
    
    if len(sine_params) % 3 != 0 or n_sines == 0:
        return np.full_like(t, np.nan) 

    model = np.zeros_like(t, dtype=float)
    
    for i in range(n_sines):
        A = sine_params[i * 3]
        P = sine_params[i * 3 + 1]
        phi = sine_params[i * 3 + 2]
        
        # Sécurité division par zéro
        if P != 0:
            model += A * np.sin(2 * np.pi * t / P + phi)
        else:
            return np.full_like(t, np.nan) 

    return model + offset

# AJOUT de y_std dans la signature
def log_prior(theta, t_min, t_max, y_amp, y_std):
    """
    Définit les limites de recherche pour les paramètres. (Priors)
    [A1, P1, phi1, A2, P2, phi2, ..., Offset]
    """
    offset = theta[-1]
    sine_params = theta[:-1]
    n_sines = len(sine_params) // 3
    duration = float(t_max - t_min)
    if duration <= 0 or not np.isfinite(duration):
        duration = 100.0

    # Bornes sur P (mêmes unités que l'axe x, souvent numéro d'époque) : évite intervalle vide si duration < 1
    p_lo = max(1e-6, min(0.5, duration * 1e-3))
    p_hi = max(p_lo * 2.0, duration * 20.0, 1.0)

    # CALCUL DE LA LIMITE ROBUSTE : 5x l'écart-type des données.
    robust_limit = max(y_amp, y_std, 1e-12) * 5.0

    # 1. Prior sur l'Offset
    if not (-robust_limit < offset < robust_limit):
        return -np.inf

    # 2. Prior sur chaque sinusoïde (A, P, phi)
    for i in range(n_sines):
        A = sine_params[i * 3]
        P = sine_params[i * 3 + 1]
        phi = sine_params[i * 3 + 2]

        # A : Positif, max 2x la limite robuste (Amplitude max autorisée)
        if not (0 < A < 2 * robust_limit):
            return -np.inf

        # P : dans une plage cohérente avec l'étendue des données
        if not (p_lo < P < p_hi):
            return -np.inf

        # phi : Entre -pi et +pi (inchangé)
        if not (-np.pi < phi < np.pi):
            return -np.inf

    return 0.0

# AJOUT de y_std dans la signature
def log_probability(theta, x, y, yerr, t_min, t_max, y_amp, y_std):
    lp = log_prior(theta, t_min, t_max, y_amp, y_std)  # AJOUT de y_std
    if not np.isfinite(lp):
        return -np.inf

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)
    if x.shape != y.shape or y.shape != yerr.shape:
        return -np.inf
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return -np.inf

    model = multi_sine_model(x, *theta)

    if np.any(np.isnan(model)) or not np.all(np.isfinite(model)):
        return -np.inf

    yerr_safe = _sanitize_yerr(yerr, y, y_amp, y_std)
    sigma2 = yerr_safe**2
    diff = y - model
    terms = diff**2 / sigma2 + np.log(sigma2)
    if not np.all(np.isfinite(terms)):
        return -np.inf
    out = lp - 0.5 * float(np.sum(terms))
    return out if np.isfinite(out) else -np.inf

# MISE À JOUR de la signature
def fit_sine_model(x, y, yerr, nwalkers=32, nsteps=10000, n_frequences=1):
    """
    Lance le MCMC pour trouver les N_frequences sinusoïdes.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)
    if x.shape != y.shape or y.shape != yerr.shape:
        raise ValueError("x, y et yerr doivent avoir la même longueur.")

    t_min, t_max = float(np.min(x)), float(np.max(x))
    y_amp = (np.max(y) - np.min(y)) / 2
    y_std = float(np.std(y)) if len(y) > 1 else 1.0  # CALCUL DE L'ÉCART-TYPE
    if y_amp == 0:
        y_amp = 1.0

    yerr = _sanitize_yerr(yerr, y, y_amp, y_std)

    # Nombre de paramètres: 3 * n_frequences (A, P, phi) + 1 (Offset)
    ndim = 3 * n_frequences + 1

    # Initialisation des walkers (Guess)
    initial_guess = []
    duration = max(float(t_max - t_min), 1e-9)
    p_lo = max(1e-6, min(0.5, duration * 1e-3))
    p_hi = max(p_lo * 2.0, duration * 20.0, 1.0)

    for i in range(n_frequences):
        A_guess = y_amp / n_frequences
        P_guess = duration / (i + 2)
        P_guess = float(np.clip(P_guess, p_lo * 1.001, p_hi * 0.999))
        phi_guess = 0.0
        initial_guess.extend([A_guess, P_guess, phi_guess])

    initial_guess.append(float(np.mean(y)))  # Ajout de l'Offset à la fin

    # Initialisation MCMC
    pos = np.array(initial_guess) + 1e-4 * np.random.randn(nwalkers, ndim)

    # S'assurer que P reste dans le prior
    for i in range(n_frequences):
        col = i * 3 + 1
        pos[:, col] = np.abs(pos[:, col])
        pos[:, col] = np.clip(pos[:, col], p_lo * 1.001, p_hi * 0.999)

    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability, 
        # AJOUT de y_std dans les arguments passés à log_probability
        args=(x, y, yerr, t_min, t_max, y_amp, y_std) 
    )
    
    sampler.run_mcmc(pos, nsteps, progress=True)
    
    # --- CORRECTION MAJEURE ICI ---
    # 1. On ignore les 25% premiers pas (burn-in) pour la stabilité
    discard = int(nsteps * 0.25)
    
    # 2. On extrait la chaîne et on l'aplatit pour avoir un tableau 2D [samples, parameters]
    # flat=True est essentiel pour que fit_ttv puisse lire les données
    flat_samples = sampler.get_chain(discard=discard, thin=1, flat=True)
    
    return flat_samples