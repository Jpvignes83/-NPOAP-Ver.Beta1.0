import logging
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u

import config  # <-- récupération des paramètres d'observatoire

logger = logging.getLogger(__name__)


class JD_BJD_Calculator:
    def __init__(self, obs_lat=None, obs_long=None, obs_elev=None):
        """
        Initialise l'objet avec les coordonnées de l'observatoire.

        Si aucun argument n'est fourni, les valeurs sont récupérées depuis
        config.OBSERVATORY (mis à jour par l'onglet d'accueil).

        Parameters
        ----------
        obs_lat : float | None
            Latitude en degrés (Nord positif, Sud négatif).
        obs_long : float | None
            Longitude en degrés (Est positif, Ouest négatif).
        obs_elev : float | None
            Altitude en mètres.
        """
        # Récupération depuis config si non fournis
        obs_cfg = getattr(config, "OBSERVATORY", {})

        if obs_lat is None:
            obs_lat = obs_cfg.get("latitude", obs_cfg.get("lat", 0.0))
        if obs_long is None:
            obs_long = obs_cfg.get("longitude", obs_cfg.get("lon", 0.0))
        if obs_elev is None:
            obs_elev = obs_cfg.get("elevation", obs_cfg.get("elev", 0.0))

        self.location = EarthLocation(
            lat=float(obs_lat) * u.deg,
            lon=float(obs_long) * u.deg,
            height=float(obs_elev) * u.m,
        )

    # ------------------------------------------------------------------
    # Parsing DATE-OBS
    # ------------------------------------------------------------------
    def _parse_date_obs(self, date_str: str) -> Time | None:
        date_str = str(date_str).strip()
        if not date_str:
            return None
        try:
            return Time(date_str, format="isot", scale="utc")
        except Exception:
            pass
        try:
            return Time(date_str, format="fits", scale="utc")
        except Exception:
            logger.error(f"DATE-OBS illisible : {date_str}")
            return None

    # ------------------------------------------------------------------
    # Parsing RA/DEC (cible) depuis le header FITS
    # ------------------------------------------------------------------
    def _parse_target_coord(self, header) -> SkyCoord | None:
        """
        Lit les coordonnées équatoriales de la cible depuis le header FITS.

        Ordre de priorité :
          1. OBJCTRA / OBJCTDEC (RA en heures, DEC en degrés, format FITS classique)
          2. CRVAL1 / CRVAL2 (référence WCS, en degrés)

        Les coordonnées sont retournées en J2000 (FK5, équinoxe J2000.0),
        conforme à l'usage recommandé avec les éphémérides JPL (DE405, HORIZONS).

        Returns
        -------
        SkyCoord | None
            Coordonnées J2000 (FK5) de la cible, ou None si absentes / invalides.
        """
        # Cadre J2000 : Earth mean equator and equinox of J2000.0 (compatible DE405)
        j2000 = "fk5"
        equinox_j2000 = "J2000.0"

        # 1. OBJCTRA / OBJCTDEC (format typique : "HH MM SS.S" et "±DD MM SS.S")
        ra_str = header.get("OBJCTRA")
        dec_str = header.get("OBJCTDEC")
        if ra_str is not None and dec_str is not None:
            try:
                return SkyCoord(
                    ra_str,
                    dec_str,
                    unit=(u.hourangle, u.deg),
                    frame=j2000,
                    equinox=equinox_j2000,
                )
            except Exception as e:
                logger.debug(f"OBJCTRA/OBJCTDEC invalides : {e}")

        # 2. CRVAL1 / CRVAL2 (WCS, en degrés)
        crval1 = header.get("CRVAL1")
        crval2 = header.get("CRVAL2")
        if crval1 is not None and crval2 is not None:
            try:
                return SkyCoord(
                    ra=float(crval1) * u.deg,
                    dec=float(crval2) * u.deg,
                    frame=j2000,
                    equinox=equinox_j2000,
                )
            except (TypeError, ValueError) as e:
                logger.debug(f"CRVAL1/CRVAL2 invalides : {e}")

        return None

    def compute_jd_bjd(self, file: str) -> bool:
        """
        Calcule et écrit dans le header FITS :
          - JD-UTC  : Julian Date (UTC) au milieu de l'exposition
          - BJD-TDB : Barycentric JD (TDB) au milieu de l'exposition

        La position de l'observatoire est prise en priorité dans le header
        (SITELAT/SITELONG/SITEELEV), sinon on utilise self.location (config).
        """
        try:
            with fits.open(file, mode="update") as hdul:
                header = hdul[0].header

                # --- DATE-OBS ---
                if "DATE-OBS" not in header:
                    logger.error(f"DATE-OBS manquant dans {file}. Calcul JD/BJD annulé.")
                    return False

                date_obs = self._parse_date_obs(header["DATE-OBS"])
                if date_obs is None:
                    logger.error(f"DATE-OBS invalide dans {file}. Calcul JD/BJD annulé.")
                    return False

                # --- EXPTIME ---
                try:
                    exptime = float(header.get("EXPTIME", 0))
                except Exception as e:
                    logger.error(f"EXPTIME invalide dans {file} : {e}")
                    return False

                if exptime <= 0:
                    logger.error(f"EXPTIME non positif dans {file}. Calcul JD/BJD annulé.")
                    return False

                # Milieu de pose (UTC)
                midpoint_obs = date_obs + (exptime / 2.0) * u.s

                # --- Coordonnées de la cible ---
                target = self._parse_target_coord(header)
                if target is None:
                    logger.warning(
                        f"RA/DEC manquants ou invalides dans {file}. "
                        f"JD-UTC écrit, BJD-TDB non calculé."
                    )
                    header["JD-UTC"] = (
                        float(midpoint_obs.jd),
                        "Julian Date at mid-exposure (UTC)",
                    )
                    return True

                # --- Coordonnées de l'observatoire ---
                lat = header.get("SITELAT")
                lon = header.get("SITELONG")
                elev = header.get("SITEELEV")

                obs_location = self.location  # valeur par défaut (config)
                if lat is not None and lon is not None:
                    try:
                        obs_location = EarthLocation(
                            lat=float(lat) * u.deg,
                            lon=float(lon) * u.deg,
                            height=float(elev) * u.m if elev is not None else 0.0 * u.m,
                        )
                    except Exception as e:
                        logger.warning(
                            f"SITELAT/SITELONG invalides dans {file}, "
                            f"utilisation de config.OBSERVATORY : {e}"
                        )

                # --- Correction barycentrique ---
                try:
                    ltt_bary = midpoint_obs.light_travel_time(target, location=obs_location)
                    bjd_tdb = midpoint_obs.tdb + ltt_bary
                except Exception as e:
                    logger.error(f"Erreur light_travel_time pour {file} : {e}")
                    # Au minimum, on écrit JD-UTC
                    header["JD-UTC"] = (
                        float(midpoint_obs.jd),
                        "Julian Date at mid-exposure (UTC)",
                    )
                    return False

                # --- Écriture dans le header ---
                header["JD-UTC"] = (
                    float(midpoint_obs.jd),
                    "Julian Date at mid-exposure (UTC)",
                )
                header["BJD-TDB"] = (
                    float(bjd_tdb.jd),
                    "Barycentric JD (TDB) at mid-exposure",
                )
                hdul.flush() 
                logger.info(f"✅ JD-UTC et BJD-TDB enregistrés pour {file}.")
                return True

        except Exception as e:
            logger.error(f"Erreur compute_jd_bjd({file}) : {e}")
            return False



