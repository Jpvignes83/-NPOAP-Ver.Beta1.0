import os
from pathlib import Path

# ──────────────────────────────────────────────────────────────
# Observatoire (valeurs par défaut – modifiables via HomeTab)
# ──────────────────────────────────────────────────────────────
OBSERVATORY = {
    "name": "",
    "lat": ,   
    "lon": ,
    "elev": ,     # mètres
    "timezone": "",
}

# ──────────────────────────────────────────────────────────────
# Chemins de sortie
# ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent.resolve()
OUTPUT_DIR = BASE_DIR / "output"
PHOTOMETRY_DIR = OUTPUT_DIR / "photometry"
CALIBRATED_DIR = OUTPUT_DIR / "calibrated"
ASTROMETRY_DIR = OUTPUT_DIR / "astrometry"

# ──────────────────────────────────────────────────────────────
# Paramètres photométrie (aperture par défaut)
# ──────────────────────────────────────────────────────────────
PHOTOMETRY_DEFAULTS = {
    "aperture_radius_fwhm": 1.3,
    "inner_annulus_fwhm": 1.8,
    "outer_annulus_fwhm": 2.5,
    "sigma_clip": 3.0,
    "fwhm_kernel_size": 21,
    "gaia_max_magnitude": 14.0,
    "gaia_cone_radius_arcmin": 15.0,
}

# ──────────────────────────────────────────────────────────────
# KBMOD via WSL (détection Synthetic Tracking)
# ──────────────────────────────────────────────────────────────
# Commande Python sous WSL pour exécuter scripts/kbmod_wsl_detect.py.
# Par défaut "python3" (Python par défaut dans le PATH WSL).
# Si vous avez installé KBMOD dans un env conda/venv, indiquez le chemin complet
# dans WSL, ex: "/home/user/miniconda3/bin/python3" ou "/home/user/miniconda3/envs/kbmod/bin/python3"
KBMOD_WSL_PYTHON = "/home/user/miniconda3/envs/astroenv/bin/python3"

# ──────────────────────────────────────────────────────────────
# Astrometry.net (clé API locale)
# ──────────────────────────────────────────────────────────────
ASTROMETRY_API_KEY_FILE = Path.home() / ".astrometry_api_key"


EQUIPMENT_OBSERVATION = {
    "obs_code": "",           # Code observateur AAVSO (5 caractères max)
    "camera": "CCD",             # Nom de la caméra
    "binning": "1x1",         # Binning (1x1, 2x2, 3x3, 4x4)
    "delim": "tab",             # Délimiteur pour rapports (, ; | : ! / ? ou tab)
    "telescope_diameter_mm": 500.0,  # Diamètre du télescope (mm) – configurable via Accueil
    "focal_length_mm": 1939,  # Focale en millimètres
    "sensor_width_mm": 0.0,   # Largeur capteur (mm) – configurable via Accueil
    "sensor_height_mm": 0.0,  # Hauteur capteur (mm) – configurable via Accueil
    "pixel_size_um": 3.76,  # Taille du pixel en micromètres
    "pixel_scale_arcsec": 0.3800,  # Échelle pixel calculée (arcsec/pixel)
}

# ──────────────────────────────────────────────────────────────
# Affiliations (modifiables via HomeTab)
# ──────────────────────────────────────────────────────────────
AFFILIATIONS = []  # Liste de dicts: [{"text": "...", "selected": bool}]

# ──────────────────────────────────────────────────────────────
# Astro-COLIBRI (recherche transitoires, cone search)
# UID optionnel : inscription gratuite sur astro-colibri.com, 100 req/jour
# ──────────────────────────────────────────────────────────────
ASTRO_COLIBRI_UID = ''

# ──────────────────────────────────────────────────────────────
# TNS (Transient Name Server) API Configuration
# User-Agent: "user" (compte utilisateur) ou "bot" (compte bot)
# ──────────────────────────────────────────────────────────────
# Chaque utilisateur doit saisir et sauvegarder son TNS ID / nom ; laisser vides dans la distribution.
TNS_CONFIG = {
    "tns_marker_type": "user",   # "user" | "bot"
    "tns_id": "",                # TNS ID utilisateur (à remplir par l'utilisateur)
    "tns_name": "",              # Nom utilisateur TNS (à remplir par l'utilisateur)
    "bot_id": "",                # ID du bot (pour type "bot")
    "api_key": "",               # Clé API TNS (à configurer via GUI)
    "bot_name": "NPOAP",
    "use_sandbox": False,
}
