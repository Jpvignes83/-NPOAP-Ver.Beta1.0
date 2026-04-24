# core/parsec_isochrones.py
"""
Isochrones PARSEC (Padova) pour l'estimation d'age d'amas.
Utilise le package ezpadova pour interroger le serveur CMD (stev.oapd.inaf.it).
Si ezpadova n'est pas installe: pip install git+https://github.com/mfouesneau/ezpadova
"""
import logging
import numpy as np

logger = logging.getLogger(__name__)

PARSEC_AVAILABLE = False
try:
    import ezpadova
    PARSEC_AVAILABLE = True
except ImportError:
    ezpadova = None


def _column_names(df):
    """Retourne les noms de colonnes pour G, G_BP, G_RP (variants possibles)."""
    cols = list(df.columns) if hasattr(df, 'columns') else []
    g = None
    gbp = None
    grp = None
    for c in cols:
        c_lower = str(c).lower()
        if c_lower in ('gmag', 'g_mag', 'phot_g_mean_mag'):
            g = c
        elif c_lower in ('g_bpmag', 'g_bp', 'bp mag', 'phot_bp_mean_mag'):
            gbp = c
        elif c_lower in ('g_rpmag', 'g_rp', 'rp mag', 'phot_rp_mean_mag'):
            grp = c
    if g is None and 'Gmag' in cols:
        g = 'Gmag'
    if gbp is None and 'G_BPmag' in cols:
        gbp = 'G_BPmag'
    if grp is None and 'G_RPmag' in cols:
        grp = 'G_RPmag'
    return g, gbp, grp


def _turnoff_from_isochrone(iso, g_col, gbp_col, grp_col):
    """
    Estime le tour de courbe (turn-off) d'une isochrone PARSEC.
    Le turn-off est le point le plus brillant (min Gmag) sur la sequence principale.
    On restreint a la partie haute (couleur pas trop rouge) pour eviter les geantes.
    """
    g = np.asarray(iso[g_col])
    gbp = np.asarray(iso[gbp_col])
    grp = np.asarray(iso[grp_col])
    color = gbp - grp
    valid = np.isfinite(g) & np.isfinite(color)
    if np.sum(valid) < 5:
        return None, None
    g = g[valid]
    color = color[valid]
    # Sequence principale: couleur typiquement 0.2 a 1.2; le turn-off est vers le bleu
    ms = (color >= 0.15) & (color <= 1.5)
    if np.sum(ms) < 3:
        return np.nan, np.nan
    idx_min = np.nanargmin(g[ms])
    color_to = float(color[ms][idx_min])
    m_g_to = float(g[ms][idx_min])
    return color_to, m_g_to


def fetch_parsec_turnoff_grid(logage_min=8.0, logage_max=10.2, step=0.2, MH=0.0):
    """
    Recupere une grille de positions tour de courbe (couleur, M_G) pour des isochrones PARSEC.
    Metallicitie solaire par defaut (MH=0).

    Returns:
        list of (log_age, color_to, M_G_to) ou None si ezpadova indisponible / erreur.
    """
    if not PARSEC_AVAILABLE or ezpadova is None:
        return None
    try:
        # logage = (min, max, step) en log10(age en annees)
        r = ezpadova.get_isochrones(
            logage=(logage_min, logage_max, step),
            MH=(MH, MH, 1),
            photsys_file='gaiaEDR3'
        )
    except Exception as e:
        logger.warning("ezpadova get_isochrones failed: %s", e)
        return None
    if r is None:
        return None
    # ezpadova peut retourner un DataFrame unique avec colonne logAge, ou un dict de DataFrames
    if hasattr(r, 'columns'):
        df = r
        if 'logAge' not in df.columns:
            # Une seule isochrone
            g_col, gbp_col, grp_col = _column_names(df)
            if g_col is None or gbp_col is None or grp_col is None:
                logger.warning("PARSEC: colonnes Gaia non trouvees: %s", list(df.columns))
                return None
            color_to, m_g_to = _turnoff_from_isochrone(df, g_col, gbp_col, grp_col)
            if color_to is None:
                return None
            log_age = logage_min
            return [(log_age, color_to, m_g_to)]
        log_ages = np.unique(df['logAge'])
        result = []
        for log_age in log_ages:
            sub = df[df['logAge'] == log_age]
            g_col, gbp_col, grp_col = _column_names(sub)
            if g_col is None or gbp_col is None or grp_col is None:
                continue
            color_to, m_g_to = _turnoff_from_isochrone(sub, g_col, gbp_col, grp_col)
            if color_to is not None and np.isfinite(color_to):
                result.append((float(log_age), color_to, m_g_to))
        return result if result else None
    if isinstance(r, dict):
        result = []
        for (log_age, mh), df in r.items():
            if hasattr(df, 'columns'):
                g_col, gbp_col, grp_col = _column_names(df)
                if g_col and gbp_col and grp_col:
                    color_to, m_g_to = _turnoff_from_isochrone(df, g_col, gbp_col, grp_col)
                    if color_to is not None and np.isfinite(color_to):
                        result.append((float(log_age), color_to, m_g_to))
        return result if result else None
    if isinstance(r, (list, tuple)):
        result = []
        for i, df in enumerate(r):
            log_age = logage_min + i * step
            if hasattr(df, 'columns'):
                g_col, gbp_col, grp_col = _column_names(df)
                if g_col and gbp_col and grp_col:
                    color_to, m_g_to = _turnoff_from_isochrone(df, g_col, gbp_col, grp_col)
                    if color_to is not None and np.isfinite(color_to):
                        result.append((float(log_age), color_to, m_g_to))
        return result if result else None
    return None


def get_parsec_turnoff_for_age(log_age, MH=0.0):
    """
    Retourne (color_to, M_G_to) pour une isochrone PARSEC a un age donne.
    log_age = log10(age en annees), ex. 9.0 pour 1 Gyr.
    """
    grid = fetch_parsec_turnoff_grid(logage_min=log_age, logage_max=log_age, step=0.1, MH=MH)
    if not grid or len(grid) == 0:
        return None, None
    return grid[0][1], grid[0][2]
