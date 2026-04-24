# core/nbody_simulation.py
# Simulation N-body pour analyse de systèmes planétaires multiples
# Adapté depuis EXOTIC (https://github.com/pearsonkyle/Nbody-ai)

import numpy as np
import logging
import copy
from astropy import units as u
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
from scipy.interpolate import interp1d

logger = logging.getLogger(__name__)

# Tentative d'import de rebound (dépendance optionnelle)
try:
    import rebound
    REBOUND_AVAILABLE = True
except ImportError:
    REBOUND_AVAILABLE = False
    logger.warning("rebound n'est pas installé. Les simulations N-body ne sont pas disponibles.")
    logger.warning("Installez avec: pip install rebound")

# Tentative d'import de ultranest (dépendance optionnelle pour le fitting)
try:
    from ultranest import ReactiveNestedSampler
    ULTRANEST_AVAILABLE = True
except ImportError:
    ULTRANEST_AVAILABLE = False
    logger.warning("ultranest n'est pas installé. Le fitting N-body n'est pas disponible.")
    logger.warning("Installez avec: pip install ultranest")

# Constantes
mearth = u.M_earth.to(u.kg)
msun = u.M_sun.to(u.kg)
mjup = u.M_jup.to(u.kg)

def empty_data(N):
    """Crée un dictionnaire vide pour stocker les paramètres orbitaux d'une simulation N-body."""
    return {
        'x': np.zeros(N),
        'y': np.zeros(N),
        'z': np.zeros(N),
        'P': np.zeros(N),
        'a': np.zeros(N),
        'e': np.zeros(N),
        'inc': np.zeros(N),
        'Omega': np.zeros(N),
        'omega': np.zeros(N),
        'M': np.zeros(N),
    }

def maxavg(x):
    """Moyenne robuste aux valeurs aberrantes."""
    return (np.percentile(np.abs(x), 75) + np.max(np.abs(x))) * 0.5

def generate_simulation(objects, Ndays=None, Noutputs=None):
    """
    Génère une simulation REBOUND pour un système multi-planètes.
    
    Parameters
    ----------
    objects : list
        Liste de dictionnaires avec les paramètres des objets (étoile + planètes).
        Format: [{'m': masse_étoile}, {'m': masse_planète1, 'P': période1, ...}, ...]
    Ndays : float, optional
        Durée de la simulation en jours.
    Noutputs : int, optional
        Nombre de points de sortie.
    
    Returns
    -------
    sim_data : dict ou rebound.Simulation
        Résultats de la simulation si Ndays et Noutputs sont fournis,
        sinon retourne l'objet Simulation REBOUND.
    """
    if not REBOUND_AVAILABLE:
        raise ImportError("rebound n'est pas installé. Installez-le avec: pip install rebound")
    
    sim = rebound.Simulation()
    sim.units = ('day', 'AU', 'Msun')
    
    for i in range(len(objects)):
        sim.add(**objects[i])
    
    sim.move_to_com()
    
    if Ndays and Noutputs:
        return integrate_simulation(sim, objects, Ndays, Noutputs)
    else:
        return sim

def integrate_simulation(sim, objects, Ndays, Noutputs):
    """
    Intègre une simulation REBOUND et enregistre les positions orbitales.
    
    Parameters
    ----------
    sim : rebound.Simulation
        Simulation REBOUND initialisée.
    objects : list
        Liste des objets du système.
    Ndays : float
        Durée de la simulation en jours.
    Noutputs : int
        Nombre de points de sortie.
    
    Returns
    -------
    sim_data : dict
        Dictionnaire contenant les données de la simulation:
        - 'pdata': liste des données orbitales pour chaque planète
        - 'star': données de position de l'étoile
        - 'times': tableau des temps
        - 'objects': liste des objets initiaux
        - 'dt': facteur de conversion pour les vitesses radiales
    """
    if not REBOUND_AVAILABLE:
        raise ImportError("rebound n'est pas installé.")
    
    ps = sim.particles
    times = np.linspace(0., Ndays, Noutputs)
    
    pdata = [empty_data(Noutputs) for i in range(len(objects) - 1)]
    star = {'x': np.zeros(Noutputs), 'y': np.zeros(Noutputs), 'z': np.zeros(Noutputs)}
    
    for i, time in enumerate(times):
        sim.integrate(time)
        
        # Enregistrer les positions de l'étoile
        for k in star.keys():
            star[k][i] = getattr(ps[0], k)
        
        # Enregistrer les données orbitales des planètes
        for j in range(1, len(objects)):
            for k in pdata[j-1].keys():
                pdata[j-1][k][i] = getattr(ps[j], k)
    
    sim_data = {
        'pdata': pdata,
        'star': star,
        'times': times,
        'objects': objects,
        'dt': Noutputs / (Ndays * 24 * 60 * 60)  # Conversion pour RV
    }
    
    return sim_data

def find_zero(t1, x1, t2, x2):
    """
    Trouve le zéro d'une ligne entre deux points (pour trouver les temps de transit).
    
    Parameters
    ----------
    t1, x1 : float
        Temps et position du premier point.
    t2, x2 : float
        Temps et position du second point.
    
    Returns
    -------
    T0 : float
        Temps du zéro (transit).
    """
    m = (x2 - x1) / (t2 - t1)
    T0 = -x1 / m + t1
    return T0

def transit_times(xp, xs, times):
    """
    Trouve les temps de transit depuis les données de position.
    
    Parameters
    ----------
    xp : array
        Positions x de la planète.
    xs : array
        Positions x de l'étoile.
    times : array
        Temps de la simulation.
    
    Returns
    -------
    tt : array
        Temps de transit.
    """
    dx = xp - xs
    tt = []
    for i in range(1, len(dx)):
        if dx[i-1] >= 0 and dx[i] <= 0:
            tt.append(find_zero(times[i-1], dx[i-1], times[i], dx[i]))
    return np.array(tt)

def predict_future_transits(objects, days, outputs, planet_index=1, t0=0.0):
    """
    Simule un système N-body et renvoie les temps de transit futurs.
    
    Parameters
    ----------
    objects : list
        Liste des objets du système (étoile + planètes).
    days : float
        Durée de la simulation en jours.
    outputs : int
        Nombre de points de sortie pour l'intégration.
    planet_index : int
        Index de la planète (1 = première planète).
    t0 : float
        Référence temporelle (BJD/JD) correspondant au temps 0 de la simulation.
    
    Returns
    -------
    result : dict
        {'epochs': array, 'Tc_rel': array, 'Tc_abs': array}
    """
    if not REBOUND_AVAILABLE:
        raise ImportError("rebound n'est pas installé. Installez-le avec: pip install rebound")
    if days <= 0:
        raise ValueError("La durée de simulation doit être > 0")
    if outputs < 10:
        raise ValueError("Le nombre de sorties doit être >= 10")
    if planet_index < 1 or planet_index >= len(objects):
        raise ValueError("planet_index invalide (doit référencer une planète)")
    
    sim_data = generate_simulation(objects, days, outputs)
    planet_data = sim_data['pdata'][planet_index - 1]
    
    Tc_rel = transit_times(planet_data['x'], sim_data['star']['x'], sim_data['times'])
    Tc_abs = Tc_rel + float(t0)
    epochs = np.arange(len(Tc_rel), dtype=int)
    
    return {'epochs': epochs, 'Tc_rel': Tc_rel, 'Tc_abs': Tc_abs}

def TTV(epochs, tt):
    """
    Calcule les TTV (O-C) en soustrayant une tendance linéaire.
    
    Parameters
    ----------
    epochs : array
        Époques orbitales.
    tt : array
        Temps de transit observés.
    
    Returns
    -------
    ttv : array
        Variations temporelles de transit (résidus O-C).
    m : float
        Période moyenne (pente).
    b : float
        Temps de référence (ordonnée à l'origine).
    """
    N = len(epochs)
    A = np.vstack([np.ones(N), epochs]).T
    b, m = np.linalg.lstsq(A, tt, rcond=None)[0]
    ttv = tt - m * np.array(epochs) - b
    return ttv, m, b

def lomb_scargle_analysis(t, y, dy=None, min_freq=None, max_freq=None, peaktol=0.05):
    """
    Analyse de périodogramme de Lomb-Scargle.
    
    Parameters
    ----------
    t : array
        Temps.
    y : array
        Valeurs observées.
    dy : array, optional
        Erreurs sur les valeurs.
    min_freq : float, optional
        Fréquence minimale.
    max_freq : float, optional
        Fréquence maximale.
    peaktol : float
        Seuil pour la détection de pics.
    
    Returns
    -------
    freq : array
        Fréquences.
    power : array
        Puissances.
    peak_periods : array
        Périodes des pics détectés.
    """
    if dy is not None:
        ls = LombScargle(t, y, dy)
    else:
        ls = LombScargle(t, y)
    
    if min_freq is None:
        min_freq = 1. / (1.5 * (max(t) - min(t)))
    
    if max_freq is None:
        max_freq = 1. / 2.5
    
    freq, power = ls.autopower(
        maximum_frequency=max_freq,
        minimum_frequency=min_freq,
        nyquist_factor=2,
        samples_per_peak=5,
        method='cython'
    )
    
    peaks, amps = find_peaks(power, height=peaktol)
    peaks = peaks[np.argsort(amps['peak_heights'])[::-1]]
    peak_periods = 1. / freq[peaks]
    
    return freq, power, peak_periods

def analyze_simulation(m, ttvfast=False):
    """
    Analyse les résultats d'une simulation N-body.
    
    Parameters
    ----------
    m : dict
        Résultats de la simulation (retour de integrate_simulation).
    ttvfast : bool
        Si True, retourne uniquement les TTV (mode rapide).
    
    Returns
    -------
    data : dict
        Dictionnaire contenant les analyses (TTV, RV, périodogrammes, etc.).
    """
    if ttvfast:
        tt = transit_times(m['pdata'][0]['x'], m['star']['x'], m['times'])
        ttv, per, t0 = TTV(np.arange(len(tt)), tt)
        return np.arange(len(ttv)), ttv, tt
    
    # Vitesses radiales
    RV = np.diff(m['star']['x']) * 1.496e11 * m['dt']  # Conversion en m/s
    
    freq, power, peak_periods = lomb_scargle_analysis(
        m['times'][1:], RV, min_freq=1./365, max_freq=1.
    )
    
    data = {
        'times': m['times'],
        'RV': {
            'freq': freq,
            'power': power,
            'signal': RV,
            'max': maxavg(RV),
            'peak_periods': peak_periods
        },
        'mstar': m['objects'][0]['m'],
        'planets': [],
        'objects': m['objects']
    }
    
    # Analyser chaque planète
    for j in range(1, len(m['objects'])):
        pdata = {}
        
        # Temps de transit
        tt = transit_times(m['pdata'][j-1]['x'], m['star']['x'], m['times'])
        if len(tt) >= 3:
            ttv, per, t0 = TTV(np.arange(len(tt)), tt)
        else:
            per = m['objects'][j]['P']
            t0 = 0
            ttv = np.array([0])
        
        # Périodogramme des TTV
        freq, power, peak_periods = lomb_scargle_analysis(
            np.arange(len(ttv)), ttv
        )
        
        # Sauvegarder les données
        for k in ['e', 'inc', 'a', 'omega']:
            pdata[k] = np.mean(m['pdata'][j-1][k])
        
        pdata['mass'] = m['objects'][j]['m']
        pdata['P'] = per
        pdata['t0'] = t0
        pdata['tt'] = tt
        pdata['ttv'] = ttv
        pdata['max'] = maxavg(ttv)
        pdata['freq'] = freq
        pdata['power'] = power
        pdata['peak_periods'] = peak_periods
        
        # Sous-échantillonnage pour le tracé
        pdata['x'] = m['pdata'][j-1]['x'][::4]
        pdata['y'] = m['pdata'][j-1]['z'][::4]
        
        data['planets'].append(pdata)
    
    return data


def interp_distribution(values, nbins=50):
    """
    Crée une fonction d'interpolation pour une distribution de probabilité.
    
    Parameters
    ----------
    values : array
        Valeurs à utiliser pour créer la distribution.
    nbins : int
        Nombre de bins pour l'histogramme.
    
    Returns
    -------
    interp_func : function
        Fonction d'interpolation de la PDF.
    """
    if len(values) == 0:
        # Retourner une fonction constante si pas de valeurs
        return lambda x: np.ones_like(x) * 0.0
    
    value_grid = np.linspace(np.min(values), np.max(values), nbins)
    heights, edges = np.histogram(values, bins=value_grid, density=True)
    edge_center = (edges[:-1] + edges[1:]) / 2
    
    # Normaliser l'histogramme
    bin_widths = np.diff(edges)
    total_area = np.sum(bin_widths * heights)
    if total_area > 0:
        normalized_heights = heights / total_area
    else:
        normalized_heights = heights
    
    # Interpoler pour créer une PDF continue
    return interp1d(edge_center, normalized_heights, kind='linear', bounds_error=False, fill_value=0)


class NBodyFitter:
    """
    Classe pour ajuster les paramètres d'un système planétaire aux observations TTV.
    
    Utilise ultranest (ReactiveNestedSampler) pour l'ajustement bayésien.
    """
    
    def __init__(self, data, prior=None, bounds=None, verbose=True):
        """
        Initialise le fitteur N-body.
        
        Parameters
        ----------
        data : list
            Liste de dictionnaires contenant les données observées.
            Format: [{}, {'Tc': array, 'Tc_err': array}, ...]
            - data[0]: données de l'étoile (vide pour l'instant)
            - data[1]: données de la planète intérieure (temps de transit observés)
            - data[2]: données de la planète extérieure (optionnel)
        prior : list
            Liste de dictionnaires avec les paramètres initiaux des objets.
            Format: [{'m': masse_étoile}, {'m': masse_p1, 'P': période_p1, ...}, ...]
        bounds : list
            Liste de dictionnaires définissant les bornes et priors pour chaque objet.
            Format: [{}, {'P': [min, max]}, {'P': [min, max], 'm': [min, max], 'P_logl': func, 'm_logl': func}, ...]
        verbose : bool
            Afficher les informations de progression.
        """
        if not REBOUND_AVAILABLE:
            raise ImportError("rebound n'est pas installé. Installez-le avec: pip install rebound")
        if not ULTRANEST_AVAILABLE:
            raise ImportError("ultranest n'est pas installé. Installez-le avec: pip install ultranest")
        
        self.data = data
        self.bounds = bounds
        self.prior = copy.deepcopy(prior)
        self.verbose = verbose
        self.results = None
        self.parameters = None
        self.errors = None
        
        if self.prior and self.bounds:
            self.fit_nested()
    
    def fit_nested(self):
        """Lance l'ajustement bayésien avec ultranest."""
        # Configuration des paramètres libres
        freekeys = []
        boundarray = []
        
        for i, planet_bounds in enumerate(self.bounds):
            for bound_key, bound_value in planet_bounds.items():
                if '_logl' in bound_key or '_fn' in bound_key:
                    continue
                freekeys.append(f"{i}_{bound_key}")
                if isinstance(bound_value, (list, tuple)) and len(bound_value) == 2:
                    boundarray.append([bound_value[0], bound_value[1]])
                else:
                    boundarray.append([bound_value, bound_value])
        
        # Trouver la durée de la simulation
        min_time = np.min(self.data[1]['Tc'])
        max_time = np.max(self.data[1]['Tc'])
        sim_time = (max_time - min_time) * 1.05
        
        # Normaliser les temps de transit
        Tc_norm = self.data[1]['Tc'] - min_time
        if self.prior[1]['P'] > 0:
            self.orbit = np.rint(Tc_norm / self.prior[1]['P']).astype(int)
        else:
            raise ValueError("Période de la planète intérieure doit être > 0")
        
        # Convertir en array numpy
        boundarray = np.array(boundarray)
        bounddiff = np.diff(boundarray, axis=1).reshape(-1)
        
        # Fonction de vraisemblance
        def loglike(pars):
            chi2 = 0.0
            
            # Mettre à jour les paramètres
            for i, par in enumerate(pars):
                idx, key = freekeys[i].split('_')
                idx = int(idx)
                if key == 'tmid':
                    continue
                self.prior[idx][key] = par
            
            # Vérifier les priors de likelihood (pour la planète extérieure typiquement)
            if len(self.bounds) > 2:
                if 'm_logl' in self.bounds[2]:
                    likelihood_m = self.bounds[2]['m_logl'](self.prior[2]['m'])
                    if likelihood_m <= 0:
                        return -1e6
                if 'P_logl' in self.bounds[2]:
                    likelihood_P = self.bounds[2]['P_logl'](self.prior[2]['P'])
                    if likelihood_P <= 0:
                        return -1e6
            
            # Lancer la simulation N-body
            try:
                sim_data = generate_simulation(self.prior, sim_time, int(sim_time * 24))
            except:
                return -1e6
            
            sim_shift = 0
            
            # Comparer avec les observations
            for i, planet in enumerate(self.prior):
                if i > 0 and self.data[i]:  # Skip étoile
                    try:
                        # Calculer les temps de transit depuis la simulation
                        Tc_sim = transit_times(
                            sim_data['pdata'][i-1]['x'],
                            sim_data['star']['x'],
                            sim_data['times']
                        )
                        
                        # Décalage temporel
                        if i-1 == 0:
                            sim_shift = Tc_sim.min()
                        
                        Tc_sim -= sim_shift
                        
                        # Comparer avec les observations
                        # Vérifier que les indices sont valides
                        if len(self.orbit) > len(Tc_sim):
                            return -1e6
                        
                        # Aligner les temps de transit simulés aux observations
                        residual = self.data[i]['Tc'] - Tc_sim[self.orbit]
                        Tc_sim_aligned = Tc_sim[self.orbit] + residual.mean()
                        
                        # Calcul du chi2
                        chi2 += -0.5 * np.sum(
                            ((self.data[i]['Tc'] - Tc_sim_aligned) / self.data[i]['Tc_err'])**2
                        )
                    except:
                        return -1e6
            
            return chi2
        
        # Transformation des priors
        def prior_transform(upars):
            return boundarray[:, 0] + bounddiff * upars
        
        # Lancer ultranest
        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e5)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(
                max_ncalls=1e5, show_status=False, viz_callback=False
            )
        
        # Extraire les résultats
        self.parameters = copy.deepcopy(self.prior)
        self.errors = {}
        self.quantiles = {}
        
        # Mettre à jour avec les résultats (si disponibles)
        if hasattr(self.results, 'maximum_likelihood') and hasattr(self.results, 'posterior'):
            for i, key in enumerate(freekeys):
                idx, param_key = key.split('_')
                idx = int(idx)
                if param_key != 'tmid':
                    if hasattr(self.results.maximum_likelihood, 'point'):
                        self.parameters[idx][param_key] = self.results.maximum_likelihood.point[i]
                    if hasattr(self.results.posterior, 'stdev'):
                        self.errors[key] = self.results.posterior.stdev[i]
                    if hasattr(self.results.posterior, 'errlo') and hasattr(self.results.posterior, 'errup'):
                        self.quantiles[key] = [
                            self.results.posterior.errlo[i],
                            self.results.posterior.errup[i]
                        ]

