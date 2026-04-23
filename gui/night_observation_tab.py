# gui/night_observation_tab.py
"""
Onglet pour les éphémérides de nuit.
Utilise l'AAVSO Target Tool pour les étoiles binaires à éclipses (index.csv).
Utilise le NASA Exoplanet Archive (TAP / pscomppars, tran_flag=1) pour les exoplanètes en transit,
avec requête restreinte à la bande de déclinaison observable depuis la latitude (config.OBSERVATORY)
et filtrage : nuit astronomique (Soleil ≤ −18°), transit avec pl_tranmid / pl_trandur / pl_orbper,
altitude > 25° à ingress−1 h et egress+1 h (dates couvertes par la section Paramètres),
séparation ≥ 50° par rapport à la Lune à ces instants (exoplanètes ; astéroïdes, comètes et transitoires à la culmination sur la nuit).
Utilise le système MPCORB du Minor Planet Center pour les astéroïdes et comètes.
Exoplanètes, astéroïdes, comètes et transitoires : exclusion si la cible est à moins de 50°
de la Lune (instant de culmination pour MPC ; ingress−1 h, milieu de transit, egress+1 h pour les transits).
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import logging
from pathlib import Path
from datetime import datetime, timedelta, time as dt_time
import pytz
import json
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import requests
import csv
import io
import re
import math
from collections import Counter
from typing import List, Optional, Tuple
from threading import Thread

import config
from config import OBSERVATORY
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, get_sun, get_body
from astropy.utils.iers import conf as iers_conf

try:
    from core.astro_colibri_client import AstroColibriClient
    ASTRO_COLIBRI_AVAILABLE = True
except Exception:
    AstroColibriClient = None
    ASTRO_COLIBRI_AVAILABLE = False

# Désactiver les avertissements IERS pour éviter les délais de téléchargement
iers_conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'

logger = logging.getLogger(__name__)

try:
    import exoclock
    EXOCLOCK_AVAILABLE = True
except Exception:
    exoclock = None
    EXOCLOCK_AVAILABLE = False

# NASA Exoplanet Archive — service TAP (données « transiting planets », voir TransitView transits) :
# https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TransitView/nph-visibletbls?dataset=transits
NASA_EXOARCHIVE_TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
NASA_EXOPLANET_LOCAL_CSV = "nasa_exoplanet_transits.csv"
ETD_EXOPLANETS_URL = "https://var.astro.cz/fr/Exoplanets"
ETD_SEARCH_API_URL = "https://var.astro.cz/api/Search/Exoplanets"
EXOPLANET_PROVIDER_NASA = "nasa"
EXOPLANET_PROVIDER_EXOCLOCK = "exoclock"
EXOPLANET_PROVIDER_ETD = "etd"


def moon_altitude_deg_at_utc(t_utc: datetime, location: EarthLocation) -> float:
    """Altitude apparente de la Lune (°) au lieu et instant UTC donnés."""
    if t_utc.tzinfo is None:
        t_utc = pytz.UTC.localize(t_utc)
    else:
        t_utc = t_utc.astimezone(pytz.UTC)
    t_ast = Time(t_utc)
    moon = get_body("moon", t_ast, location=location)
    frame = AltAz(obstime=t_ast, location=location)
    return float(moon.transform_to(frame).alt.deg)


def moon_object_separation_deg_icrs(
    ra_deg: float, dec_deg: float, obstime: Time, location: EarthLocation
) -> float:
    """Séparation angulaire apparente (°) entre la cible ICRS et la Lune au lieu et instant donnés."""
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    moon = get_body("moon", obstime, location=location)
    return float(coord.separation(moon).deg)


def moon_phase_french_name_and_illumination(t_utc: datetime) -> Tuple[str, float, float]:
    """
    Phase de la Lune à l’instant donné (UTC).

    Returns
    -------
    name_fr : str
    illumination : float
        Fraction illuminée 0–1 (0 = nouvelle, 1 = pleine).
    elong_deg : float
        Élongation géocentrique Soleil–Lune (°).
    """
    if t_utc.tzinfo is None:
        t_utc = pytz.UTC.localize(t_utc)
    else:
        t_utc = t_utc.astimezone(pytz.UTC)
    t_ast = Time(t_utc)
    sun = get_sun(t_ast)
    moon = get_body("moon", t_ast, location=None)
    elong = float(sun.separation(moon).deg)
    elong_rad = math.radians(elong)
    illumination = 0.5 * (1.0 - math.cos(elong_rad))
    e = elong % 360.0
    if e < 22.5 or e >= 337.5:
        name = "Nouvelle lune"
    elif e < 67.5:
        name = "Premier croissant"
    elif e < 112.5:
        name = "Premier quartier"
    elif e < 157.5:
        name = "Gibbeuse croissante"
    elif e < 202.5:
        name = "Pleine lune"
    elif e < 247.5:
        name = "Gibbeuse décroissante"
    elif e < 292.5:
        name = "Dernier quartier"
    else:
        name = "Dernier croissant"
    return name, illumination, elong


def declination_bounds_for_latitude(lat_deg: float) -> Tuple[float, float]:
    """
    Bande de déclinaison nécessaire pour qu'un objet puisse dépasser l'horizon géométrique
    au moins une fois depuis la latitude donnée (φ) : δ ∈ [max(-90, φ−90), min(90, φ+90)].
    Utilisé pour restreindre la requête TAP (config.OBSERVATORY lat).
    """
    dec_min = max(-90.0, lat_deg - 90.0)
    dec_max = min(90.0, lat_deg + 90.0)
    return dec_min, dec_max


def build_nasa_exoplanet_tap_adql(dec_min: float, dec_max: float) -> str:
    """ADQL : planètes en transit dont la déclinaison peut être observable depuis la latitude du site."""
    return (
        "SELECT pl_name, ra, dec, sy_vmag, sy_gaiamag, sy_jmag, sy_kmag, pl_orbper, pl_tranmid, pl_trandur, "
        "pl_ratror, pl_rade, pl_radj, st_rad FROM pscomppars "
        "WHERE tran_flag=1 "
        f"AND dec >= {dec_min:.6f} AND dec <= {dec_max:.6f}"
    )


def host_star_magnitude_from_nasa_row(row: dict) -> Optional[float]:
    """
    Magnitude de l’étoile hôte depuis une ligne TAP NASA (pscomppars).
    Essaie sy_vmag puis bandes de repli ; plage réaliste pour affichage / filtrage.
    """
    for key in ("sy_vmag", "sy_gaiamag", "sy_jmag", "sy_kmag"):
        raw = row.get(key)
        if raw is None or str(raw).strip() == "":
            continue
        try:
            x = float(raw)
            if math.isfinite(x) and -4.0 <= x <= 25.0:
                return x
        except (TypeError, ValueError):
            continue
    return None


# Rayons solaires (pour recréer « transitdepthcalc » depuis pscomppars, cf. NASA Transit Service)
R_EARTH_OVER_R_SUN = 6.3781e6 / 6.957e8
R_JUP_OVER_R_SUN = 7.1492e7 / 6.957e8


def transitdepthcalc_percent_from_row(row: dict) -> float:
    """
    Profondeur « Transit Depth - Calculated [%] » comme `transitdepthcalc` du service
    Transit & Ephemeris NASA : (Rp/R*)^2 × 100 à partir des rayons dans pscomppars ;
    R* = 1 R☉ si st_rad absent ; 0 si aucun rayon planétaire (comme l’archive).
    https://exoplanetarchive.ipac.caltech.edu/docs/transit/transit_parameters.html
    """
    def _pos_float(key: str) -> Optional[float]:
        v = row.get(key)
        if v is None or str(v).strip() == "":
            return None
        try:
            x = float(v)
            return x if x > 0 else None
        except (TypeError, ValueError):
            return None

    st_rad_solar = _pos_float("st_rad")
    if st_rad_solar is None:
        st_rad_solar = 1.0

    rat = _pos_float("pl_ratror")
    if rat is not None:
        rp_over_rs = rat
    else:
        pl_rade = _pos_float("pl_rade")
        pl_radj = _pos_float("pl_radj")
        if pl_rade is not None:
            rp_over_rs = pl_rade * R_EARTH_OVER_R_SUN / st_rad_solar
        elif pl_radj is not None:
            rp_over_rs = pl_radj * R_JUP_OVER_R_SUN / st_rad_solar
        else:
            return 0.0

    return 100.0 * rp_over_rs * rp_over_rs


def _jd_to_utc_datetime(jd: float) -> datetime:
    """Convertit un jour julien (NASA pl_tranmid) en datetime UTC aware."""
    dt = Time(jd, format="jd").to_datetime()
    if dt.tzinfo is None:
        return pytz.UTC.localize(dt)
    return dt.astimezone(pytz.UTC)


def exo_transit_k_for_night_utc_interval(
    pl_tranmid_jd: float,
    pl_orbper: float,
    pl_trandur_h: float,
    night_start_utc: datetime,
    night_end_utc: datetime,
) -> Optional[int]:
    """
    Premier indice d’époque k tel que [ingress−1 h, egress+1 h] est inclus dans
    [night_start_utc, night_end_utc] (même critère que le filtre exoplanète).
    """
    if pl_orbper <= 0 or pl_trandur_h <= 0:
        return None
    half_dur_d = (pl_trandur_h / 2.0) / 24.0
    margin_d = 1.0 / 24.0
    jd_start = float(Time(night_start_utc).jd)
    jd_end = float(Time(night_end_utc).jd)
    jd_lo = jd_start + half_dur_d + margin_d
    jd_hi = jd_end - half_dur_d - margin_d
    if jd_lo > jd_hi:
        return None
    t0 = pl_tranmid_jd
    p = pl_orbper
    k_min = math.ceil((jd_lo - t0) / p)
    k_max = math.floor((jd_hi - t0) / p)
    if k_min > k_max:
        return None
    return k_min


# Filtrage exoplanètes : fenêtre [ingress−1 h, egress+1 h] dans la nuit astronomique, alt > 25°.
EXO_TRANSIT_MIN_ALT_DEG = 25.0
EXO_ASTRO_SUN_ALT_DEG = -18.0
# Exoplanètes, astéroïdes, comètes, transitoires : écart minimal à la Lune (critère de recherche nuit).
NIGHT_MOON_MIN_SEPARATION_DEG = 50.0


def comet_mpes_designation_from_mpc_line(line: str) -> Optional[str]:
    """
    Désignation utilisable par le MPES à partir d'une ligne de AllCometEls.txt.
    Le fichier mélange formats ; ce n'est en principe pas « C/… » ou « P/… » dans les 12
    premiers caractères (souvent « 0123P », « CD04Y010 », etc.).
    """
    s = line.strip()
    if len(s) < 8:
        return None
    tokens = s.split()
    if not tokens:
        return None

    # 1) Lignes récentes : après une date AAAA/MM/JJ au format YYYYMMDD, deux colonnes
    #    numériques (souvent magnitudes), puis la désignation jusqu'à MPC / MPEC / réf. « n, ».
    for i, t in enumerate(tokens):
        if not re.fullmatch(r"(19|20)\d{6}", t):
            continue
        j = i + 1
        for _ in range(2):
            if j < len(tokens) and re.match(
                r"^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$", tokens[j]
            ):
                j += 1
        parts: List[str] = []
        while j < len(tokens):
            tj = tokens[j]
            if tj == "MPEC":
                break
            if tj == "MPC" or tj.startswith("MPC"):
                break
            if re.match(r"^\d+,$", tj):
                break
            parts.append(tj)
            j += 1
        if parts:
            des = " ".join(parts)
            if "/" in des:
                return des

    # 2) Comète périodique compacte en tête : 0109P, 0141P, fragment éventuel « a »… sur le token suivant.
    m = re.match(r"^0*(\d+)P$", tokens[0], re.I)
    if m:
        n = m.group(1)
        base = f"{n}P"
        if len(tokens) > 1 and len(tokens[1]) == 1 and tokens[1].isalpha():
            return f"{base}-{tokens[1].upper()}"
        return base

    # 3) Désignation « normale » dès le premier champ ou plus loin (C/, D/, X/, P/…),
    #    y compris « C/ 240 V1 » (tokens séparés après le slash).
    for i, tok in enumerate(tokens):
        if tok.startswith(("C/", "D/", "X/")) or (
            tok.startswith("P/") and len(tok) > 2
        ):
            parts = [tok]
            j = i + 1
            while j < len(tokens):
                tj = tokens[j]
                if tj == "MPEC" or tj == "MPC" or tj.startswith("MPC"):
                    break
                if re.match(r"^\d+,$", tj):
                    break
                if re.fullmatch(r"[-+]?\d+\.\d+", tj):
                    break
                # Année / numéro de comète non périodique (ex. 240 entre C/ et V1).
                if re.fullmatch(r"-?\d{1,5}", tj):
                    parts.append(tj)
                    j += 1
                    continue
                # Fragment type A1, Y1, P1.
                if re.fullmatch(r"[A-Za-z]\d*", tj) and len(tj) <= 6:
                    parts.append(tj)
                    j += 1
                    continue
                break
            cand = " ".join(parts).strip()
            if "/" in cand and len(cand) > 3:
                return cand

    return None


# Classes et fonctions utilitaires
class EphemerisObject:
    """Représente un objet avec ses éphémérides."""
    
    def __init__(self, name: str, ra: float, dec: float, magnitude: Optional[float] = None,
                 obj_type: str = "unknown", period: Optional[float] = None,
                 exo_pl_tranmid_jd: Optional[float] = None,
                 exo_pl_trandur_h: Optional[float] = None,
                 exo_transitdepthcalc_pct: Optional[float] = None,
                 exo_source: Optional[str] = None):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.magnitude = magnitude
        self.obj_type = obj_type
        self.period = period
        # Données TAP NASA (transit), pour schéma sur le graphique d’altitude
        self.exo_pl_tranmid_jd = exo_pl_tranmid_jd
        self.exo_pl_trandur_h = exo_pl_trandur_h
        # Équivalent colonne transitdepthcalc (profondeur calculée depuis Rp/R*), [%]
        self.exo_transitdepthcalc_pct = exo_transitdepthcalc_pct
        # Source exoplanète (ex. "nasa", "exoclock", "etd")
        self.exo_source = (exo_source or "").strip().lower() or None
        
        # Éphémérides calculées
        self.alt_max = 0.0
        self.alt_max_time = None
        self.rise_time = None
        self.set_time = None
        self.transit_time = None
        
    def ra_sexagesimal(self) -> str:
        """Retourne l'ascension droite au format hh:mm:ss.ss."""
        ra_hours = self.ra / 15.0
        hours = int(ra_hours)
        minutes_float = (ra_hours - hours) * 60.0
        minutes = int(minutes_float)
        seconds = (minutes_float - minutes) * 60.0
        return f"{hours:02d}:{minutes:02d}:{seconds:05.2f}"
    
    def dec_sexagesimal(self) -> str:
        """Retourne la déclinaison au format deg:min:sec.ss."""
        dec_abs = abs(self.dec)
        sign = '+' if self.dec >= 0 else '-'
        degrees = int(dec_abs)
        minutes_float = (dec_abs - degrees) * 60.0
        minutes = int(minutes_float)
        seconds = (minutes_float - minutes) * 60.0
        return f"{sign}{degrees:02d}:{minutes:02d}:{seconds:05.2f}"
    
    def coordinates_sexagesimal(self) -> str:
        """Retourne les coordonnées au format RA=hh:mm:ss.ss DEC=deg:min:sec.ss."""
        return f"RA={self.ra_sexagesimal()} DEC={self.dec_sexagesimal()}"


def parse_ra_dec(ra_str: str, dec_str: str) -> Tuple[Optional[float], Optional[float]]:
    """Parse RA et DEC depuis des chaînes au format AAVSO."""
    ra_deg = None
    try:
        if not ra_str or not isinstance(ra_str, str):
            return None, None
        
        ra_str = ra_str.strip()
        # Remplacer les séparateurs par des espaces
        ra_clean = ra_str.replace('h', ' ').replace('m', ' ').replace('s', ' ').strip()
        parts = ra_clean.split()
        
        if len(parts) >= 3:
            h = float(parts[0])
            m = float(parts[1])
            s = float(parts[2])
            hours = h + m/60.0 + s/3600.0
            ra_deg = hours * 15.0
        elif len(parts) == 2:
            h = float(parts[0])
            m = float(parts[1])
            hours = h + m/60.0
            ra_deg = hours * 15.0
        elif len(parts) == 1:
            # Format décimal en heures
            ra_deg = float(parts[0]) * 15.0
    except (ValueError, IndexError):
        pass
    
    dec_deg = None
    try:
        if not dec_str or not isinstance(dec_str, str):
            return ra_deg, None
        
        dec_str = dec_str.strip()
        # Remplacer les séparateurs par des espaces
        dec_clean = dec_str.replace('d', ' ').replace('°', ' ').replace("'", ' ').replace('"', ' ').replace('m', ' ').replace('s', ' ').strip()
        parts = dec_clean.split()
        
        if len(parts) >= 3:
            deg = float(parts[0])
            m = float(parts[1])
            s = float(parts[2])
            dec_deg = deg + (m/60.0 if deg >= 0 else -m/60.0) + (s/3600.0 if deg >= 0 else -s/3600.0)
        elif len(parts) == 2:
            deg = float(parts[0])
            m = float(parts[1])
            dec_deg = deg + (m/60.0 if deg >= 0 else -m/60.0)
        elif len(parts) == 1:
            # Format décimal en degrés
            dec_deg = float(parts[0])
    except (ValueError, IndexError):
        pass
    
    return ra_deg, dec_deg


def normalize_exoplanet_name(name: str) -> str:
    """
    Normalise le format des noms exoplanètes pour éviter les variantes
    du type "WASP-12 b" vs "WASP-12b".
    Convention retenue: suffixe planète collé au numéro (ex: WASP-12b).
    """
    txt = " ".join(str(name or "").strip().split())
    if not txt:
        return txt
    # Ex: "WASP-12 b" -> "WASP-12b"
    txt = re.sub(r"^(.+?\d)\s+([a-z])$", r"\1\2", txt, flags=re.IGNORECASE)
    return txt


def parse_icrs_sexagesimal_to_deg(
    ra_h: str,
    ra_m: str,
    ra_s: str,
    dec_d: str,
    dec_m: str,
    dec_s: str,
    dec_south: bool,
) -> Optional[Tuple[float, float]]:
    """
    Convertit RA (h m s) et DEC (d m s) en degrés ICRS.
    DEC : case « Sud » ou degrés négatifs → déclinaison négative.
    """
    def _f(x: str) -> Optional[float]:
        t = (x or "").strip().replace(",", ".")
        if t == "":
            return None
        try:
            return float(t)
        except ValueError:
            return None

    h = _f(ra_h)
    if h is None:
        return None
    m_ra = _f(ra_m) or 0.0
    s_ra = _f(ra_s) or 0.0
    ra_hours = h + m_ra / 60.0 + s_ra / 3600.0
    if ra_hours < 0 or ra_hours >= 24.0:
        return None
    ra_deg = ra_hours * 15.0

    d0 = _f(dec_d)
    if d0 is None:
        return None
    m_dec = _f(dec_m) or 0.0
    s_dec = _f(dec_s) or 0.0
    sign = -1.0 if (dec_south or d0 < 0) else 1.0
    dec_abs = abs(d0) + m_dec / 60.0 + s_dec / 3600.0
    if dec_abs > 90.0:
        return None
    dec_deg = sign * dec_abs
    return ra_deg, dec_deg


class NightEphemerides:
    """Calcule les éphémérides pour une nuit donnée."""
    
    def __init__(self, latitude: float, longitude: float, elevation: float):
        self.location = EarthLocation(
            lat=latitude * u.deg,
            lon=longitude * u.deg,
            height=elevation * u.m
        )

    def _sun_altitude_deg_local(self, t_local: datetime) -> float:
        """Altitude du Soleil (°) au lieu / instant donnés, t_local timezone-aware."""
        t_utc = t_local.astimezone(pytz.UTC)
        t_ast = Time(t_utc)
        sun = get_sun(t_ast)
        frame = AltAz(obstime=t_ast, location=self.location)
        return float(sun.transform_to(frame).alt.deg)
    
    def get_sunrise_sunset(self, obs_date: datetime, timezone_str: str = "America/Santiago"):
        """
        Calcule le lever et coucher du soleil (altitude = 0°) pour une date donnée.
        
        Parameters
        ----------
        obs_date : datetime
            Date d'observation (datetime naive, sera interprétée comme UTC)
        timezone_str : str
            Nom du fuseau horaire (ex: "America/Santiago")
        
        Returns
        -------
        sunset_time : datetime
            Heure de coucher du soleil (heure locale)
        sunrise_time : datetime
            Heure de lever du soleil (heure locale)
        local_timezone : pytz.timezone
            Fuseau horaire local
        """
        try:
            # Obtenir le fuseau horaire local (gère automatiquement heure d'été/hiver)
            local_tz = pytz.timezone(timezone_str)
        except pytz.exceptions.UnknownTimeZoneError:
            # Fallback : calculer depuis la longitude
            logger.warning(f"Timezone {timezone_str} inconnu, utilisation de la longitude")
            lon = self.location.lon.deg
            tz_offset = int(lon / 15.0)
            local_tz = pytz.FixedOffset(tz_offset * 60)
        
        # Convertir la date d'observation en UTC si nécessaire
        if obs_date.tzinfo is None:
            obs_date_utc = pytz.UTC.localize(obs_date)
        else:
            obs_date_utc = obs_date.astimezone(pytz.UTC)
        
        # Convertir en heure locale
        obs_date_local = obs_date_utc.astimezone(local_tz)
        
        # Date locale à minuit pour les calculs
        local_date = obs_date_local.date()
        local_midnight = local_tz.localize(datetime.combine(local_date, datetime.min.time()))
        
        # Coucher / lever géométriques : balayage minute par minute (midi → midi+1 jour)
        # (l’ancienne boucle par heure entière ratrait dans la mauvaise heure si le passage
        #  horizon / −18° tombait avant HH:00, d’où des défauts 18 h / 6 h.)
        next_date = local_date + timedelta(days=1)
        next_midnight = local_tz.localize(datetime.combine(next_date, datetime.min.time()))
        noon = local_tz.localize(datetime.combine(local_date, dt_time(12, 0, 0)))
        noon_next = noon + timedelta(days=1)

        sunset_time = None
        sunrise_time = None
        cur = noon
        prev_alt = None
        step = timedelta(minutes=1)
        while cur <= noon_next:
            alt = self._sun_altitude_deg_local(cur)
            if prev_alt is not None:
                if sunset_time is None and prev_alt > 0.0 and alt <= 0.0:
                    sunset_time = cur
                if sunset_time is not None and sunrise_time is None and prev_alt < 0.0 and alt >= 0.0:
                    sunrise_time = cur
                    break
            prev_alt = alt
            cur += step

        # Valeurs par défaut si non trouvées
        if sunset_time is None:
            sunset_time = local_midnight.replace(hour=18, minute=0, second=0, microsecond=0)
        if sunrise_time is None:
            sunrise_time = next_midnight.replace(hour=6, minute=0, second=0, microsecond=0)
        
        return sunset_time, sunrise_time, local_tz

    def get_astronomical_twilight_bounds(
        self, obs_date: datetime, timezone_str: str = "America/Santiago"
    ):
        """
        Début et fin de la nuit astronomique (Soleil à −18°) : fin du crépuscule astronomique le soir,
        début de l’aube astronomique le lendemain matin. Même convention de fuseau que get_sunrise_sunset.
        """
        thresh = EXO_ASTRO_SUN_ALT_DEG
        try:
            local_tz = pytz.timezone(timezone_str)
        except pytz.exceptions.UnknownTimeZoneError:
            logger.warning(f"Timezone {timezone_str} inconnu, utilisation de la longitude")
            lon = self.location.lon.deg
            tz_offset = int(lon / 15.0)
            local_tz = pytz.FixedOffset(tz_offset * 60)

        if obs_date.tzinfo is None:
            obs_date_utc = pytz.UTC.localize(obs_date)
        else:
            obs_date_utc = obs_date.astimezone(pytz.UTC)

        obs_date_local = obs_date_utc.astimezone(local_tz)
        local_date = obs_date_local.date()
        local_midnight = local_tz.localize(datetime.combine(local_date, datetime.min.time()))
        next_date = local_date + timedelta(days=1)
        next_midnight = local_tz.localize(datetime.combine(next_date, datetime.min.time()))

        # Midi local du jour d’observation → midi du lendemain : une rotation diurne complète,
        # pas de raffinement piégé sur une heure entière incorrecte.
        noon = local_tz.localize(datetime.combine(local_date, dt_time(12, 0, 0)))
        scan_end = noon + timedelta(hours=30)
        astro_dusk = None
        astro_dawn = None
        cur = noon
        prev_alt = None
        step = timedelta(minutes=1)
        while cur <= scan_end:
            alt = self._sun_altitude_deg_local(cur)
            if prev_alt is not None:
                if astro_dusk is None and prev_alt > thresh and alt <= thresh:
                    astro_dusk = cur
                if astro_dusk is not None and astro_dawn is None and prev_alt <= thresh and alt > thresh:
                    astro_dawn = cur
                    break
            prev_alt = alt
            cur += step

        if astro_dusk is None:
            logger.warning(
                "Crépuscule astronomique (−18°) non détecté sur la fenêtre scannée — repli ~21 h locale"
            )
            astro_dusk = local_tz.localize(datetime.combine(local_date, dt_time(21, 0, 0)))
        if astro_dawn is None:
            logger.warning(
                "Aube astronomique (−18°) non détectée — repli ~5 h locale (jour suivant)"
            )
            astro_dawn = local_tz.localize(datetime.combine(next_date, dt_time(5, 0, 0)))

        return astro_dusk, astro_dawn, local_tz
    
    def get_night_period(self, obs_date: datetime):
        """Calcule les crépuscules astronomiques, nautiques et civils."""
        t_obs = Time(obs_date)
        
        # Calculer le lever et coucher du soleil
        sun = get_sun(t_obs)
        frame = AltAz(obstime=t_obs, location=self.location)
        sun_altaz = sun.transform_to(frame)
        
        # Crépuscule astronomique: altitude du soleil < -18°
        # Crépuscule nautique: altitude du soleil < -12°
        # Crépuscule civil: altitude du soleil < -6°
        
        # Approche: calculer pour plusieurs heures autour de minuit
        sunset_astro = None
        sunrise_astro = None
        sunset_naut = None
        sunrise_naut = None
        sunset_civil = None
        sunrise_civil = None
        
        # Chercher le coucher (soir)
        for hour in range(12, 24):
            t_test = Time(obs_date.replace(hour=hour, minute=0, second=0))
            sun_test = get_sun(t_test)
            frame_test = AltAz(obstime=t_test, location=self.location)
            altaz_test = sun_test.transform_to(frame_test)
            
            if altaz_test.alt.deg < -18 and sunset_astro is None:
                # Interpoler pour trouver l'heure exacte
                for minute in range(0, 60, 5):
                    t_interp = Time(obs_date.replace(hour=hour, minute=minute, second=0))
                    sun_interp = get_sun(t_interp)
                    frame_interp = AltAz(obstime=t_interp, location=self.location)
                    altaz_interp = sun_interp.transform_to(frame_interp)
                    if altaz_interp.alt.deg < -18:
                        sunset_astro = t_interp.to_datetime()
                        break
            if altaz_test.alt.deg < -12 and sunset_naut is None:
                for minute in range(0, 60, 5):
                    t_interp = Time(obs_date.replace(hour=hour, minute=minute, second=0))
                    sun_interp = get_sun(t_interp)
                    frame_interp = AltAz(obstime=t_interp, location=self.location)
                    altaz_interp = sun_interp.transform_to(frame_interp)
                    if altaz_interp.alt.deg < -12:
                        sunset_naut = t_interp.to_datetime()
                        break
            if altaz_test.alt.deg < -6 and sunset_civil is None:
                for minute in range(0, 60, 5):
                    t_interp = Time(obs_date.replace(hour=hour, minute=minute, second=0))
                    sun_interp = get_sun(t_interp)
                    frame_interp = AltAz(obstime=t_interp, location=self.location)
                    altaz_interp = sun_interp.transform_to(frame_interp)
                    if altaz_interp.alt.deg < -6:
                        sunset_civil = t_interp.to_datetime()
                        break
        
        # Chercher le lever (matin)
        next_date = obs_date + timedelta(days=1)
        for hour in range(0, 12):
            t_test = Time(next_date.replace(hour=hour, minute=0, second=0))
            sun_test = get_sun(t_test)
            frame_test = AltAz(obstime=t_test, location=self.location)
            altaz_test = sun_test.transform_to(frame_test)
            
            if altaz_test.alt.deg > -18 and sunrise_astro is None:
                for minute in range(0, 60, 5):
                    t_interp = Time(next_date.replace(hour=hour, minute=minute, second=0))
                    sun_interp = get_sun(t_interp)
                    frame_interp = AltAz(obstime=t_interp, location=self.location)
                    altaz_interp = sun_interp.transform_to(frame_interp)
                    if altaz_interp.alt.deg > -18:
                        sunrise_astro = t_interp.to_datetime()
                        break
            if altaz_test.alt.deg > -12 and sunrise_naut is None:
                for minute in range(0, 60, 5):
                    t_interp = Time(next_date.replace(hour=hour, minute=minute, second=0))
                    sun_interp = get_sun(t_interp)
                    frame_interp = AltAz(obstime=t_interp, location=self.location)
                    altaz_interp = sun_interp.transform_to(frame_interp)
                    if altaz_interp.alt.deg > -12:
                        sunrise_naut = t_interp.to_datetime()
                        break
            if altaz_test.alt.deg > -6 and sunrise_civil is None:
                for minute in range(0, 60, 5):
                    t_interp = Time(next_date.replace(hour=hour, minute=minute, second=0))
                    sun_interp = get_sun(t_interp)
                    frame_interp = AltAz(obstime=t_interp, location=self.location)
                    altaz_interp = sun_interp.transform_to(frame_interp)
                    if altaz_interp.alt.deg > -6:
                        sunrise_civil = t_interp.to_datetime()
                        break
        
        # Valeurs par défaut si non trouvées
        if sunset_astro is None:
            sunset_astro = obs_date.replace(hour=18, minute=0, second=0, microsecond=0)
        if sunrise_astro is None:
            sunrise_astro = (obs_date + timedelta(days=1)).replace(hour=6, minute=0, second=0, microsecond=0)
        if sunset_naut is None:
            sunset_naut = sunset_astro - timedelta(hours=1)
        if sunrise_naut is None:
            sunrise_naut = sunrise_astro + timedelta(hours=1)
        if sunset_civil is None:
            sunset_civil = sunset_naut - timedelta(hours=1)
        if sunrise_civil is None:
            sunrise_civil = sunrise_naut + timedelta(hours=1)
        
        return sunset_naut, sunset_astro, sunrise_astro, sunrise_naut, sunset_civil, sunrise_civil
    
    def calculate_ephemeris(self, obj: EphemerisObject, obs_date: datetime,
                           start_time: datetime, end_time: datetime):
        """Calcule les éphémérides pour un objet."""
        coord = SkyCoord(ra=obj.ra * u.deg, dec=obj.dec * u.deg, frame='icrs')
        
        # S'assurer que start_time et end_time sont cohérents en termes de timezone
        if start_time.tzinfo is None:
            start_time_tz = pytz.UTC.localize(start_time)
        else:
            start_time_tz = start_time.astimezone(pytz.UTC) if start_time.tzinfo else pytz.UTC.localize(start_time)
        
        if end_time.tzinfo is None:
            end_time_tz = pytz.UTC.localize(end_time)
        else:
            end_time_tz = end_time.astimezone(pytz.UTC) if end_time.tzinfo else pytz.UTC.localize(end_time)
        
        times = []
        current = start_time_tz
        while current <= end_time_tz:
            times.append(current)
            current = current + timedelta(minutes=15)
        
        altitudes = []
        for t in times:
            # Convertir datetime en Time astropy (gère automatiquement les timezones)
            t_astro = Time(t)
            frame = AltAz(obstime=t_astro, location=self.location)
            altaz = coord.transform_to(frame)
            altitudes.append(altaz.alt.deg)
        
        if altitudes:
            max_idx = np.argmax(altitudes)
            obj.alt_max = altitudes[max_idx]
            obj.alt_max_time = times[max_idx]
            obj.transit_time = times[max_idx]
            
            for i, alt in enumerate(altitudes):
                if alt > 0:
                    obj.rise_time = times[i]
                    break
            for i in range(len(altitudes) - 1, -1, -1):
                if altitudes[i] > 0:
                    obj.set_time = times[i]
                    break


class NightObservationTab(ttk.Frame):
    """Onglet pour calculer les éphémérides des objets observables pendant la nuit."""
    
    def __init__(self, parent):
        super().__init__(parent)
        logger.info("Initialisation de l'onglet Observation de la Nuit")
        
        # Configuration observatoire
        self.location = None
        self.load_observatory_config()
        
        # Répertoire de sauvegarde
        self.npoap_dir = Path.home() / ".npoap"
        self.npoap_dir.mkdir(exist_ok=True)
        self.catalogues_dir = self.npoap_dir / "catalogues"
        self.catalogues_dir.mkdir(exist_ok=True)
        # Dossiers d'enregistrement (modifiables par l'utilisateur), défaut : .npoap/catalogues
        self.mpc_catalogues_dir_var = tk.StringVar(value=str(self.catalogues_dir))
        self.aavso_catalogues_dir_var = tk.StringVar(value=str(self.catalogues_dir))
        
        # Données
        self.objects = []
        self.filtered_objects = []
        self.selected_objects = {}
        self.last_exoplanet_provider_counts = {}
        self.etd_authenticated = False

        # Callback optionnel vers l'onglet Planétarium (C2A)
        # Renseigné par MainWindow via set_c2a_visualizer.
        self.c2a_visualizer = None
        
        # Variables Tkinter
        self.date_var = tk.StringVar(value=datetime.now().strftime("%Y-%m-%d"))
        self.date_end_var = tk.StringVar(value="")
        self.use_interval_var = tk.BooleanVar(value=False)
        
        # Cases à cocher types d'objets
        self.asteroids_var = tk.BooleanVar(value=False)
        self.exoplanets_var = tk.BooleanVar(value=True)
        self.comets_var = tk.BooleanVar(value=False)
        self.eclipsing_binaries_var = tk.BooleanVar(value=True)
        
        # Coordonnées manuelles ICRS (sexagésimal) pour le graphique d'altitude
        self.manual_ra_h_var = tk.StringVar(value="")
        self.manual_ra_m_var = tk.StringVar(value="")
        self.manual_ra_s_var = tk.StringVar(value="")
        self.manual_dec_d_var = tk.StringVar(value="")
        self.manual_dec_m_var = tk.StringVar(value="")
        self.manual_dec_s_var = tk.StringVar(value="")
        self.manual_dec_south_var = tk.BooleanVar(value=False)
        self.manual_export_nina_var = tk.BooleanVar(value=False)
        self.manual_coord_show_var = tk.BooleanVar(value=False)

        # Recherche transitoires Astro-COLIBRI (UID + plage de dates)
        self.colibri_uid_var = tk.StringVar(value=getattr(config, "ASTRO_COLIBRI_UID", "") or "")
        self.colibri_date_start_var = tk.StringVar(value="")
        self.colibri_date_end_var = tk.StringVar(value="")
        
        self.create_widgets()
        self.update_catalog_status()
    
    def load_observatory_config(self):
        """Charge la configuration de l'observatoire depuis config.py."""
        try:
            lat = float(OBSERVATORY.get('lat', 0.0))
            lon = float(OBSERVATORY.get('lon', 0.0))
            elev = float(OBSERVATORY.get('elev', 0.0))
            timezone_name = OBSERVATORY.get('timezone', 'America/Santiago')
            
            self.location = EarthLocation(
                lat=lat * u.deg,
                lon=lon * u.deg,
                height=elev * u.m
            )
            
            # Déterminer le fuseau horaire
            # Si c'est "Santiago, Chili", utiliser "America/Santiago"
            if 'Santiago' in timezone_name or 'Chili' in timezone_name:
                self.timezone_str = "America/Santiago"
            else:
                # Essayer d'extraire ou utiliser directement
                self.timezone_str = timezone_name if '/' in timezone_name else "America/Santiago"
            
            logger.info(f"Observatoire chargé: Lat={lat}°, Lon={lon}°, Alt={elev}m, Timezone={self.timezone_str}")
        except Exception as e:
            logger.error(f"Erreur chargement observatoire: {e}", exc_info=True)
            try:
                self.location = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)
                self.timezone_str = "UTC"
            except Exception:
                self.location = EarthLocation(lat=0, lon=0, height=0)
                self.timezone_str = "UTC"
    
    def create_widgets(self):
        """Crée les widgets de l'interface."""
        
        # Zone principale : gauche = catalogues + paramètres + liste ;
        # droite = graphique (right_col) sur toute la hauteur. L'espace horizontal restant après
        # left_col est partagé 90 % / 10 % : right_col plus étroit de 10 %, marge vide à l'extrême droite.
        content_frame = ttk.Frame(self)
        content_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        content_frame.grid_rowconfigure(0, weight=1)
        content_frame.grid_columnconfigure(0, weight=0)
        content_frame.grid_columnconfigure(1, weight=9)
        content_frame.grid_columnconfigure(2, weight=1)

        left_col = ttk.Frame(content_frame)
        right_col = ttk.Frame(content_frame)
        right_margin = ttk.Frame(content_frame)
        left_col.grid(row=0, column=0, sticky="nsew")
        right_col.grid(row=0, column=1, sticky="nsew", padx=(10, 0))
        right_margin.grid(row=0, column=2, sticky="nsew")
        
        # ===== CADRE 1 : TÉLÉCHARGEMENT DES CATALOGUES =====
        catalog_frame = ttk.LabelFrame(left_col, text="1. Téléchargement des Catalogues", padding=10)
        catalog_frame.pack(fill=tk.X, pady=(0, 5))
        
        # MPC : dossier choisi par l'utilisateur + téléchargement
        mpc_dir_row = ttk.Frame(catalog_frame)
        mpc_dir_row.pack(fill=tk.X, pady=(0, 2))
        ttk.Label(mpc_dir_row, text="Dossier MPC (NEA.txt, AllCometEls.txt) :").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(mpc_dir_row, textvariable=self.mpc_catalogues_dir_var, width=50).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=2
        )
        ttk.Button(mpc_dir_row, text="📁", command=self.browse_mpc_catalogues_dir, width=3).pack(side=tk.LEFT)

        mpc_row = ttk.Frame(catalog_frame)
        mpc_row.pack(fill=tk.X, pady=2)
        ttk.Button(
            mpc_row,
            text="Télécharger MPC (Astéroïdes/Comètes)",
            command=self.download_mpc_catalogs,
        ).pack(side=tk.LEFT, padx=5)
        self.mpc_status_label = ttk.Label(mpc_row, text="Non téléchargé", foreground="gray")
        self.mpc_status_label.pack(side=tk.LEFT, padx=10)
        
        # AAVSO : dossier pour index.csv + téléchargement
        aavso_dir_row = ttk.Frame(catalog_frame)
        aavso_dir_row.pack(fill=tk.X, pady=(4, 2))
        ttk.Label(
            aavso_dir_row,
            text="Dossier AAVSO / NASA (index.csv, nasa_exoplanet_transits.csv) :",
        ).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(aavso_dir_row, textvariable=self.aavso_catalogues_dir_var, width=50).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=2
        )
        ttk.Button(aavso_dir_row, text="📁", command=self.browse_aavso_catalogues_dir, width=3).pack(side=tk.LEFT)

        aavso_row = ttk.Frame(catalog_frame)
        aavso_row.pack(fill=tk.X, pady=2)
        ttk.Button(
            aavso_row,
            text="Télécharger AAVSO (binaires à éclipses)",
            command=self.download_aavso_index,
        ).pack(side=tk.LEFT, padx=5)
        self.aavso_status_label = ttk.Label(aavso_row, text="Non téléchargé", foreground="gray")
        self.aavso_status_label.pack(side=tk.LEFT, padx=10)
        
        # ===== CADRE 2 : PARAMÈTRES =====
        params_frame = ttk.LabelFrame(left_col, text="2. Paramètres", padding=10)
        params_frame.pack(fill=tk.X, pady=5)
        
        # Date / Intervalle
        date_row = ttk.Frame(params_frame)
        date_row.pack(fill=tk.X, pady=5)
        
        ttk.Label(date_row, text="Date d'observation:").pack(side=tk.LEFT, padx=5)
        date_entry = ttk.Entry(date_row, textvariable=self.date_var, width=12)
        date_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Checkbutton(date_row, text="Intervalle", variable=self.use_interval_var,
                       command=self.toggle_date_interval).pack(side=tk.LEFT, padx=5)
        
        self.date_end_entry = ttk.Entry(date_row, textvariable=self.date_end_var, width=12, state='disabled')
        self.date_end_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(date_row, text="Aujourd'hui", command=self.set_today).pack(side=tk.LEFT, padx=5)

        # UI NASA masquée temporairement (phrase + bouton).
        # Le backend de téléchargement est conservé pour réactivation ultérieure.
        self.nasa_exo_status_label = None
        
        # Types d'objets
        types_row = ttk.Frame(params_frame)
        types_row.pack(fill=tk.X, pady=5)
        
        ttk.Label(types_row, text="Types d'objets:").pack(side=tk.LEFT, padx=5)
        
        ttk.Checkbutton(types_row, text="Astéroïdes", variable=self.asteroids_var).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(types_row, text="Exoplanètes", variable=self.exoplanets_var).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(types_row, text="Comètes", variable=self.comets_var).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(types_row, text="Étoiles binaires à éclipses", 
                       variable=self.eclipsing_binaries_var).pack(side=tk.LEFT, padx=5)

        # Recherche de transitoires via Astro-COLIBRI (UID + plage de dates)
        colibri_frame = ttk.LabelFrame(params_frame, text="Transitoires (Astro-COLIBRI)", padding=6)
        colibri_frame.pack(fill=tk.X, pady=(8, 4))

        uid_row = ttk.Frame(colibri_frame)
        uid_row.pack(fill=tk.X, pady=2)
        ttk.Label(uid_row, text="UID Astro-COLIBRI:").pack(side=tk.LEFT, padx=(0, 4))
        ttk.Entry(uid_row, textvariable=self.colibri_uid_var, width=26).pack(side=tk.LEFT, padx=2)

        dates_row = ttk.Frame(colibri_frame)
        dates_row.pack(fill=tk.X, pady=2)
        ttk.Label(dates_row, text="Date début (YYYY-MM-DD):").pack(side=tk.LEFT, padx=(0, 4))
        ttk.Entry(dates_row, textvariable=self.colibri_date_start_var, width=12).pack(side=tk.LEFT, padx=2)
        ttk.Label(dates_row, text="Date fin:").pack(side=tk.LEFT, padx=(8, 4))
        ttk.Entry(dates_row, textvariable=self.colibri_date_end_var, width=12).pack(side=tk.LEFT, padx=2)

        ttk.Button(
            colibri_frame,
            text="🔍 Rechercher (Astro-COLIBRI → Objets observables)",
            command=self.search_colibri_for_observables,
        ).pack(fill=tk.X, pady=(6, 0))

        calc_row = ttk.Frame(params_frame)
        calc_row.pack(fill=tk.X, pady=2)
        ttk.Button(
            calc_row,
            text="🔭 Calculer les éphémérides",
            command=self.calculate_ephemerides,
            style="Accent.TButton",
        ).pack(side=tk.LEFT, padx=5)
        
        # Bouton exporter JSON NINA
        buttons_row = ttk.Frame(params_frame)
        buttons_row.pack(fill=tk.X, pady=5)
        
        ttk.Button(
            buttons_row,
            text="📄 Exporter JSON NINA (objets sélectionnés)",
            command=self.export_nina_json,
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            buttons_row,
            text="Réinitialiser",
            command=self.reset_observable_objects_view,
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            buttons_row,
            text="🪐 Envoyer vers l'onglet Planétarium (C2A)",
            command=self.on_visualize_in_c2a,
        ).pack(side=tk.LEFT, padx=5)
        
        # ===== LISTE DES OBJETS =====
        list_frame = ttk.LabelFrame(left_col, text="Objets observables", padding=5)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=(5, 0), ipadx=5)
        
        # Scrollbars pour la liste
        list_scrollbar_y = ttk.Scrollbar(list_frame, orient=tk.VERTICAL)
        list_scrollbar_x = ttk.Scrollbar(list_frame, orient=tk.HORIZONTAL)
        
        # Treeview avec case à cocher virtuelle (colonne 0 = checkbox)
        columns = ("Sélection", "Nom", "Type", "Mag vis", "Coordonnées", "Alt. max")
        self.results_tree = ttk.Treeview(list_frame, columns=columns, show="headings",
                                         yscrollcommand=list_scrollbar_y.set,
                                         xscrollcommand=list_scrollbar_x.set,
                                         height=25)
        
        # Configuration des colonnes — cliquer sur « Tout » coche/décoche toutes les lignes visibles
        self.results_tree.heading("Sélection", text="Tout")
        self.results_tree.column("Sélection", width=58, anchor='center')
        
        # Astéroïdes/comètes : désignation MPC (souvent « packed », ex. 00433 ou K07Tf8A), pas le nom IAU.
        self.results_tree.heading("Nom", text="Désignation / nom")
        self.results_tree.column("Nom", width=150, anchor='w')
        
        self.results_tree.heading("Type", text="Type")
        self.results_tree.column("Type", width=80, anchor='center')
        
        self.results_tree.heading("Mag vis", text="Mag vis")
        self.results_tree.column("Mag vis", width=58, anchor='center')
        
        self.results_tree.heading("Coordonnées", text="Coordonnées")
        self.results_tree.column("Coordonnées", width=180, anchor='center')
        
        self.results_tree.heading("Alt. max", text="Alt. max")
        self.results_tree.column("Alt. max", width=60, anchor='center')
        self.results_tree.tag_configure("source_nasa", foreground="blue")
        self.results_tree.tag_configure("source_exoclock", foreground="red")
        
        list_scrollbar_y.config(command=self.results_tree.yview)
        list_scrollbar_x.config(command=self.results_tree.xview)
        
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        list_scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        list_scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Bind pour sélection/désélection sur clic
        self.results_tree.bind('<Button-1>', self.on_tree_click)
        
        # Graphique : occupe toute la largeur de right_col (elle-même ~90 % de l'espace extensible).
        plot_frame = ttk.LabelFrame(right_col, text="Graphique d'altitude", padding=5)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        
        self.fig, self.ax = plt.subplots(figsize=(14, 9))
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def _toggle_select_all_observables(self) -> None:
        """Coche toutes les lignes visibles si au moins une est décochée ; sinon tout décoche."""
        children = self.results_tree.get_children()
        if not children:
            return
        rows: List[Tuple[str, str, List[str]]] = []
        for item in children:
            values = self.results_tree.item(item, "values")
            if len(values) > 1:
                obj_name = self._object_name_from_display_name(values[1])
                rows.append((item, obj_name, list(values)))
        if not rows:
            return
        n_selected = sum(1 for _, on, _ in rows if self.selected_objects.get(on, False))
        select_all = n_selected < len(rows)
        for item, obj_name, values_list in rows:
            self.selected_objects[obj_name] = select_all
            values_list[0] = "☑" if select_all else "☐"
            self.results_tree.item(item, values=tuple(values_list))
        self.update_plot()

    def on_tree_click(self, event):
        """Gère le clic sur le treeview pour sélectionner/désélectionner les objets
        et, lors de la sélection, affiche une fenêtre de détails + met en évidence la trajectoire."""
        region = self.results_tree.identify_region(event.x, event.y)
        if region == "heading":
            column = self.results_tree.identify_column(event.x)
            if column == "#1":
                self._toggle_select_all_observables()
            return
        if region == "cell":
            column = self.results_tree.identify_column(event.x)
            item = self.results_tree.identify_row(event.y)
            
            # Si on clique sur la colonne "Sélection" (colonne 0)
            if column == "#1" and item:
                # Récupérer l'objet associé
                values = self.results_tree.item(item, "values")
                if len(values) > 1:
                    obj_name = self._object_name_from_display_name(values[1])
                    # Toggle sélection
                    currently_selected = self.selected_objects.get(obj_name, False)
                    now_selected = not currently_selected
                    self.selected_objects[obj_name] = now_selected

                    # Mettre à jour l'affichage
                    new_check = "☑" if now_selected else "☐"
                    values_list = list(values)
                    values_list[0] = new_check
                    self.results_tree.item(item, values=values_list)

                    # Mettre à jour le graphique (trajetoire dans la fenêtre d'altitude)
                    self.update_plot()

                    # Si l'objet vient d'être sélectionné, afficher une fenêtre de détails
                    if now_selected:
                        obj = self._find_object_by_name(obj_name)
                        if obj is not None:
                            self._show_object_details_window(obj)

    @staticmethod
    def _object_name_from_display_name(display_name: str) -> str:
        """Retire les préfixes visuels de source dans la liste."""
        txt = str(display_name or "")
        for prefix in ("[NASA] ", "[ExoClock] "):
            if txt.startswith(prefix):
                return txt[len(prefix):]
        return txt

    @staticmethod
    def _display_name_for_object(obj: "EphemerisObject") -> str:
        source = (getattr(obj, "exo_source", None) or "").strip().lower()
        if source == EXOPLANET_PROVIDER_NASA:
            return f"[NASA] {obj.name}"
        if source == EXOPLANET_PROVIDER_EXOCLOCK:
            return f"[ExoClock] {obj.name}"
        return obj.name
    
    def _find_object_by_name(self, name: str) -> Optional["EphemerisObject"]:
        """Retrouve un objet par son nom dans la liste complète."""
        for obj in self.objects:
            if obj.name == name:
                return obj
        return None

    def _show_object_details_window(self, obj: "EphemerisObject") -> None:
        """Ouvre une petite fenêtre de détails pour l'objet sélectionné."""
        win = tk.Toplevel(self)
        win.title(f"Détails objet : {obj.name}")
        win.geometry("420x220")

        frame = ttk.Frame(win, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        mag_str = f"{obj.magnitude:.2f}" if obj.magnitude is not None else "N/A"

        id_label = (
            "Désignation MPC"
            if getattr(obj, "obj_type", "") in ("asteroid", "comet")
            else "Nom"
        )
        lw = 16
        lines = [
            f"{id_label:<{lw}}: {obj.name}",
            f"{'Type':<{lw}}: {obj.obj_type}",
            f"{'RA (hh:mm:ss)':<{lw}}: {obj.ra_sexagesimal()}",
            f"{'Dec (dd:mm:ss)':<{lw}}: {obj.dec_sexagesimal()}",
            f"{'Mag visuelle':<{lw}}: {mag_str}",
            (
                f"{'Alt. max':<{lw}}: {obj.alt_max:.1f}°"
                if hasattr(obj, "alt_max")
                else f"{'Alt. max':<{lw}}: N/A"
            ),
        ]

        text = tk.Text(frame, height=8, wrap=tk.WORD, font=("Consolas", 10))
        text.pack(fill=tk.BOTH, expand=True)
        text.insert("1.0", "\n".join(lines))
        text.config(state="disabled")

        btn_row = ttk.Frame(frame)
        btn_row.pack(fill=tk.X, pady=5)

        ttk.Button(btn_row, text="Fermer", command=win.destroy).pack(side=tk.LEFT)
        ttk.Button(
            btn_row,
            text="🪐 Visualiser cette cible dans C2A",
            command=lambda: self._visualize_single_object_in_c2a(obj),
        ).pack(side=tk.LEFT, padx=8)

    # ------------------------------------------------------------------
    # Intégration C2A
    # ------------------------------------------------------------------
    def set_c2a_visualizer(self, callback) -> None:
        """
        Enregistre une fonction pour pousser une liste d’objets vers l’onglet Planétarium.

        callback(objs: List[EphemerisObject]) -> None
        """
        self.c2a_visualizer = callback

    def _get_selected_objects(self) -> List["EphemerisObject"]:
        """
        Retourne la liste des objets actuellement cochés dans la liste.
        """
        selected_list = []
        for name, selected in self.selected_objects.items():
            if not selected:
                continue
            obj = self._find_object_by_name(name)
            if obj is not None:
                selected_list.append(obj)
        return selected_list

    def get_observation_date_interval_yyyy_mm_dd(self) -> Tuple[str, str]:
        """
        Lit l’intervalle calendaire des champs « Date d’observation » et « Intervalle »
        (date de fin), comme pour le calcul d’éphémérides / le TAP exoplanètes.

        Retourne (date_début, date_fin) au format ``YYYY-MM-DD`` avec ``date_fin`` >= ``date_début``.
        Si la date de début est invalide : ``("", "")``. Sans intervalle coché : fin = début.
        Si la date de fin est vide ou invalide alors que l’intervalle est coché : fin = début.
        """
        raw_start = self.date_var.get().strip()
        try:
            t_start = datetime.strptime(raw_start, "%Y-%m-%d")
        except ValueError:
            return "", ""

        d_start = t_start.strftime("%Y-%m-%d")
        d_end = d_start

        if self.use_interval_var.get():
            raw_end = (self.date_end_var.get() or "").strip()
            if raw_end:
                try:
                    t_end = datetime.strptime(raw_end, "%Y-%m-%d")
                    if t_end < t_start:
                        t_start, t_end = t_end, t_start
                    d_start = t_start.strftime("%Y-%m-%d")
                    d_end = t_end.strftime("%Y-%m-%d")
                except ValueError:
                    pass

        return d_start, d_end

    def get_export_bundle_for_sora(self) -> Tuple[List[str], str, str, float, float, float, str]:
        """
        Astéroïdes et comètes cochés + intervalle calendaire (UTC) + site observatoire pour SORA.

        Retourne
        --------
        bodies, time_beg, time_end, lon_deg, lat_deg, height_m, site_name

        Les dates reprennent ``get_observation_date_interval_yyyy_mm_dd()`` (champs Paramètres
        de l’onglet). ``time_beg`` / ``time_end`` sont ``AAAA-MM-JJT00:00:00`` et
        ``…T23:59:59`` du dernier jour de l’intervalle.
        Les coordonnées reprennent ``self.location`` (config observatoire).
        """
        bodies: List[str] = []
        for obj in self._get_selected_objects():
            if getattr(obj, "obj_type", "") in ("asteroid", "comet"):
                n = (obj.name or "").strip()
                if n:
                    bodies.append(n)
        seen = set()
        uniq: List[str] = []
        for b in bodies:
            if b not in seen:
                seen.add(b)
                uniq.append(b)
        bodies = uniq

        d0, d1 = self.get_observation_date_interval_yyyy_mm_dd()
        if d0:
            time_beg = f"{d0}T00:00:00"
            time_end = f"{d1}T23:59:59"
        else:
            time_beg = ""
            time_end = ""

        if self.location is not None:
            lat_deg = float(self.location.lat.to(u.deg).value)
            lon_deg = float(self.location.lon.to(u.deg).value)
            height_m = float(self.location.height.to(u.m).value)
        else:
            lat_deg = float(OBSERVATORY.get("lat", 0.0))
            lon_deg = float(OBSERVATORY.get("lon", 0.0))
            height_m = float(OBSERVATORY.get("elev", 0.0))
        site_name = str(OBSERVATORY.get("name", "Site") or "Site")

        return bodies, time_beg, time_end, lon_deg, lat_deg, height_m, site_name

    def on_visualize_in_c2a(self) -> None:
        """
        Handler du bouton 'Envoyer vers l'onglet Planétarium (C2A)'.
        """
        if self.c2a_visualizer is None:
            messagebox.showwarning(
                "C2A",
                "Le lien vers l'onglet Planétarium (C2A) n'est pas initialisé.",
            )
            return

        objs = self._get_selected_objects()
        if not objs:
            messagebox.showinfo(
                "C2A",
                "Aucune cible n'est sélectionnée dans la liste des objets observables.\n"
                "Cochez au moins un objet dans la colonne ☑ puis réessayez.",
            )
            return

        self.c2a_visualizer(objs)

    def _visualize_single_object_in_c2a(self, obj: "EphemerisObject") -> None:
        """
        Appelé depuis la fenêtre de détails pour envoyer directement cet objet
        dans la liste de l’onglet Planétarium.
        """
        if self.c2a_visualizer is None:
            messagebox.showwarning(
                "C2A",
                "Le lien vers l'onglet Planétarium (C2A) n'est pas initialisé.",
            )
            return
        self.c2a_visualizer([obj])
    
    def toggle_date_interval(self):
        """Active/désactive le champ de date de fin pour l'intervalle."""
        if self.use_interval_var.get():
            self.date_end_entry.config(state='normal')
            if not self.date_end_var.get():
                # Si vide, mettre la date de fin = date de début + 1 jour
                start_date = datetime.strptime(self.date_var.get(), "%Y-%m-%d")
                end_date = start_date + timedelta(days=1)
                self.date_end_var.set(end_date.strftime("%Y-%m-%d"))
        else:
            self.date_end_entry.config(state='disabled')
            self.date_end_var.set("")

    def search_colibri_for_observables(self):
        """Recherche des transitoires Astro-COLIBRI par plage de dates et les projette dans « Objets observables »."""
        if not ASTRO_COLIBRI_AVAILABLE or AstroColibriClient is None:
            messagebox.showerror("Erreur", "Module Astro-COLIBRI indisponible.")
            return

        uid = self.colibri_uid_var.get().strip()
        if not uid:
            messagebox.showwarning(
                "Attention",
                "Un UID Astro-COLIBRI est requis (inscription gratuite sur astro-colibri.com).",
            )
            return

        # Plage de dates : si non renseignée, utiliser la date d'observation / intervalle de l'onglet.
        date_start = (self.colibri_date_start_var.get() or "").strip()
        date_end = (self.colibri_date_end_var.get() or "").strip()

        try:
            if not date_start:
                date_start = self.date_var.get().strip()
            if not date_end:
                date_end = (self.date_end_var.get().strip() or date_start)

            datetime.strptime(date_start, "%Y-%m-%d")
            datetime.strptime(date_end, "%Y-%m-%d")
        except ValueError:
            messagebox.showerror("Erreur", "Format de date Astro-COLIBRI invalide. Utilisez YYYY-MM-DD.")
            return

        Thread(
            target=self._search_colibri_for_observables_task,
            args=(uid, date_start, date_end),
            daemon=True,
        ).start()

    def _search_colibri_for_observables_task(self, uid: str, date_start: str, date_end: str):
        """Thread de recherche Astro-COLIBRI puis projection vers la liste « Objets observables »."""
        try:
            client = AstroColibriClient(uid=uid)
            time_min = f"{date_start}T00:00:00"
            time_max = f"{date_end}T23:59:59"
            results = client.latest_transients(time_min=time_min, time_max=time_max) or []
        except Exception as e:
            logger.error("Erreur recherche Astro-COLIBRI (nuit) : %s", e, exc_info=True)
            self.after(
                0,
                lambda: messagebox.showerror(
                    "Erreur", f"Erreur recherche Astro-COLIBRI pour la nuit :\n{e}"
                ),
            )
            return

        if not results:
            self.after(0, lambda: messagebox.showinfo("Info", "Aucun transitoire Astro-COLIBRI trouvé."))
            return

        def _project_results():
            added = 0
            # S'assurer que la fenêtre de nuit et les éphémérides sont initialisées
            try:
                obs_date = datetime.strptime(self.date_var.get().strip(), "%Y-%m-%d")
            except ValueError:
                # Si la date de l'onglet est invalide, utiliser la date de début Astro-COLIBRI
                obs_date = datetime.strptime(date_start, "%Y-%m-%d")

            # Initialiser la nuit si besoin
            if getattr(self, "start_time", None) is None or getattr(self, "end_time", None) is None:
                self._setup_night_window_for_plot(obs_date)
            else:
                # Met à jour obs_date pour cohérence
                self.obs_date = obs_date

            for evt in results:
                name = (
                    evt.get("source_name")
                    or evt.get("trigger_id")
                    or "Transient"
                )
                ra = evt.get("ra")
                dec = evt.get("dec")
                try:
                    if ra is None or dec is None:
                        continue
                    ra_deg = float(ra)
                    dec_deg = float(dec)
                except (TypeError, ValueError):
                    continue

                # Magnitude approximative si disponible
                mag = None
                for key in ("mag", "magnitude", "mag_value", "mag_last"):
                    val = evt.get(key)
                    if val is not None and val not in ("", -1):
                        try:
                            mag = float(val)
                        except (TypeError, ValueError):
                            mag = None
                        break

                obj = EphemerisObject(
                    name=str(name),
                    ra=ra_deg,
                    dec=dec_deg,
                    magnitude=mag,
                    obj_type="transient",
                )

                # Calculer les éphémérides pour ce transitoire (alt_max, transit, etc.)
                try:
                    if hasattr(self, "night_eph") and self.start_time and self.end_time:
                        self.night_eph.calculate_ephemeris(obj, self.obs_date, self.start_time, self.end_time)
                except Exception as e:
                    logger.error("Erreur calcul éphémérides Astro-COLIBRI pour %s : %s", obj.name, e, exc_info=True)

                if not self._object_clear_of_moon(obj):
                    continue

                self.objects.append(obj)
                self.filtered_objects.append(obj)
                self.selected_objects.setdefault(obj.name, False)
                added += 1

            # Rafraîchir la fenêtre « Objets observables »
            if added > 0:
                self.display_results(
                    getattr(self, "sunset_local", None),
                    getattr(self, "astro_dusk_local", getattr(self, "sunset_local", None)),
                    getattr(self, "astro_dawn_local", getattr(self, "sunrise_local", None)),
                    getattr(self, "sunrise_local", None),
                )
                # Mettre à jour le graphique d'altitude avec les nouveaux objets
                if hasattr(self, "obs_date") and getattr(self, "start_time", None) and getattr(self, "end_time", None):
                    self.plot_ephemerides(self.obs_date, self.start_time, self.end_time)
                messagebox.showinfo(
                    "Astro-COLIBRI",
                    f"{added} transitoire(s) ajouté(s) à « Objets observables ».",
                )
            else:
                messagebox.showinfo(
                    "Astro-COLIBRI",
                    "Aucun transitoire exploitable (coordonnées manquantes ou invalides).",
                )

        self.after(0, _project_results)

    def reset_observable_objects_view(self):
        """
        Vide la zone résultats (liste + graphique) et supprime le CSV exoplanètes NASA
        enregistré dans le dossier catalogue (nasa_exoplanet_transits.csv).
        """
        nasa_path = self._nasa_exoplanet_csv_path()
        if nasa_path.exists():
            file_msg = f"\n\nSupprimer aussi le fichier :\n{nasa_path}"
        else:
            file_msg = "\n\n(Aucun fichier nasa_exoplanet_transits.csv dans ce dossier.)"
        if not messagebox.askyesno(
            "Réinitialiser",
            "Vider la liste « Objets observables », effacer le graphique d’altitude "
            "et les données de calcul en mémoire ?"
            + file_msg,
        ):
            return

        self.objects = []
        self.filtered_objects = []
        self.selected_objects = {}
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)

        for attr in (
            "obs_date",
            "start_time",
            "end_time",
            "sunset_local",
            "sunrise_local",
            "astro_dusk_local",
            "astro_dawn_local",
            "local_tz",
            "night_eph",
        ):
            if hasattr(self, attr):
                try:
                    delattr(self, attr)
                except AttributeError:
                    pass

        self.ax.clear()
        self.ax.text(
            0.5,
            0.5,
            "Aucun objet — utilisez « Calculer les éphémérides ».",
            ha="center",
            va="center",
            transform=self.ax.transAxes,
            fontsize=12,
        )
        self.ax.set_xlabel("Heure locale", fontsize=11)
        self.ax.set_ylabel("Altitude (°)", fontsize=11)
        self.ax.set_title("Graphique d’altitude", fontsize=13, fontweight="bold")
        self.canvas.draw()

        if nasa_path.exists():
            try:
                nasa_path.unlink()
                logger.info("Fichier supprimé : %s", nasa_path)
            except OSError as e:
                logger.warning("Impossible de supprimer %s : %s", nasa_path, e)
                messagebox.showwarning(
                    "Fichier",
                    f"La réinitialisation a réussi, mais le fichier n’a pas pu être supprimé :\n{nasa_path}\n{e}",
                )

        self.update_catalog_status()
    
    def update_plot(self):
        """Met à jour le graphique avec les objets sélectionnés (et la courbe coordonnées manuelles si cochée)."""
        if not hasattr(self, "obs_date") or not hasattr(self, "start_time"):
            return
        self.plot_ephemerides(self.obs_date, self.start_time, self.end_time)
    
    def _mpc_files_directory(self) -> Path:
        """Répertoire local des fichiers NEA.txt et AllCometEls.txt (champ ou défaut .npoap/catalogues)."""
        s = self.mpc_catalogues_dir_var.get().strip() if hasattr(self, "mpc_catalogues_dir_var") else ""
        return Path(s) if s else self.catalogues_dir

    def browse_mpc_catalogues_dir(self):
        """Ouvre un dialogue pour choisir le dossier d'enregistrement des catalogues MPC."""
        initial = self._mpc_files_directory()
        if not initial.is_dir():
            initial = self.catalogues_dir
        chosen = filedialog.askdirectory(
            title="Dossier pour NEA.txt et AllCometEls.txt",
            initialdir=str(initial),
        )
        if chosen:
            self.mpc_catalogues_dir_var.set(chosen)
            self.update_catalog_status()

    def _aavso_files_directory(self) -> Path:
        """Répertoire contenant index.csv (AAVSO)."""
        s = self.aavso_catalogues_dir_var.get().strip() if hasattr(self, "aavso_catalogues_dir_var") else ""
        return Path(s) if s else self.catalogues_dir

    def _aavso_index_csv_path(self) -> Path:
        return self._aavso_files_directory() / "index.csv"

    def _nasa_exoplanet_csv_path(self) -> Path:
        """Fichier CSV local issu du TAP NASA (pscomppars, planètes en transit)."""
        return self._aavso_files_directory() / NASA_EXOPLANET_LOCAL_CSV

    def _night_eph(self) -> NightEphemerides:
        """Éphémérides de nuit pour l’observatoire (config.OBSERVATORY → self.location)."""
        return NightEphemerides(
            latitude=self.location.lat.deg,
            longitude=self.location.lon.deg,
            elevation=self.location.height.to(u.m).value,
        )

    def _setup_night_window_for_plot(self, obs_date: datetime) -> None:
        """Définit la nuit astronomique (−18°), lever/coucher 0° et fuseau pour le graphique."""
        self.obs_date = obs_date
        self.night_eph = self._night_eph()
        sunset_local, sunrise_local, local_tz = self.night_eph.get_sunrise_sunset(
            obs_date, self.timezone_str
        )
        self.local_tz = local_tz
        self.sunset_local = sunset_local
        self.sunrise_local = sunrise_local
        astro_dusk_local, astro_dawn_local, _ = self.night_eph.get_astronomical_twilight_bounds(
            obs_date, self.timezone_str
        )
        self.astro_dusk_local = astro_dusk_local
        self.astro_dawn_local = astro_dawn_local
        self.start_time = astro_dusk_local.astimezone(pytz.UTC)
        self.end_time = astro_dawn_local.astimezone(pytz.UTC)

    def _manual_coord_deg_for_plot(self) -> Optional[Tuple[float, float]]:
        if not self.manual_coord_show_var.get():
            return None
        return parse_icrs_sexagesimal_to_deg(
            self.manual_ra_h_var.get(),
            self.manual_ra_m_var.get(),
            self.manual_ra_s_var.get(),
            self.manual_dec_d_var.get(),
            self.manual_dec_m_var.get(),
            self.manual_dec_s_var.get(),
            self.manual_dec_south_var.get(),
        )

    def refresh_altitude_plot_if_manual(self) -> None:
        """Recalcule le graphique lorsque la case « Afficher sur le graphique » change."""
        if self.manual_coord_show_var.get():
            parsed = parse_icrs_sexagesimal_to_deg(
                self.manual_ra_h_var.get(),
                self.manual_ra_m_var.get(),
                self.manual_ra_s_var.get(),
                self.manual_dec_d_var.get(),
                self.manual_dec_m_var.get(),
                self.manual_dec_s_var.get(),
                self.manual_dec_south_var.get(),
            )
            if parsed is None:
                messagebox.showwarning(
                    "Coordonnées",
                    "Complétez RA (h m s) et δ (° ′ ″) avec des valeurs valides, "
                    "ou décochez « Afficher sur le graphique ».",
                )
                self.manual_coord_show_var.set(False)
                return
        try:
            obs_date = datetime.strptime(self.date_var.get().strip(), "%Y-%m-%d")
        except ValueError:
            if self.manual_coord_show_var.get():
                messagebox.showerror("Date", "Date d'observation invalide (AAAA-MM-JJ).")
            return
        if getattr(self, "start_time", None) is None or getattr(self, "end_time", None) is None:
            self._setup_night_window_for_plot(obs_date)
        self.plot_ephemerides(obs_date, self.start_time, self.end_time)

    def refresh_altitude_plot_with_manual(self) -> None:
        """Bouton : met à jour le graphique (vérifie les coordonnées si la case est cochée)."""
        if self.manual_coord_show_var.get():
            parsed = parse_icrs_sexagesimal_to_deg(
                self.manual_ra_h_var.get(),
                self.manual_ra_m_var.get(),
                self.manual_ra_s_var.get(),
                self.manual_dec_d_var.get(),
                self.manual_dec_m_var.get(),
                self.manual_dec_s_var.get(),
                self.manual_dec_south_var.get(),
            )
            if parsed is None:
                messagebox.showerror(
                    "Coordonnées",
                    "RA (0–24 h en h, m, s) et δ (−90° à +90°) invalides ou incomplets.",
                )
                return
        try:
            obs_date = datetime.strptime(self.date_var.get().strip(), "%Y-%m-%d")
        except ValueError:
            messagebox.showerror("Date", "Date d'observation invalide (AAAA-MM-JJ).")
            return
        self._setup_night_window_for_plot(obs_date)
        self.plot_ephemerides(obs_date, self.start_time, self.end_time)

    def _utc_astronomical_night_segments(
        self, start_date: datetime, end_date: Optional[datetime]
    ) -> List[Tuple[datetime, datetime]]:
        """
        Pour chaque date locale entre start_date et end_date (inclus), retourne
        (fin crépuscule astro, fin nuit / aube astro) en UTC timezone-aware.
        """
        ne = self._night_eph()
        tz_name = self.timezone_str

        if end_date is None:
            days = [start_date.date()]
        else:
            d0, d1 = start_date.date(), end_date.date()
            if d1 < d0:
                d0, d1 = d1, d0
            days = []
            cur = d0
            while cur <= d1:
                days.append(cur)
                cur += timedelta(days=1)

        segments: List[Tuple[datetime, datetime]] = []
        for d in days:
            obs = datetime.combine(d, datetime.min.time())
            dusk_local, dawn_local, _ = ne.get_astronomical_twilight_bounds(obs, tz_name)
            segments.append(
                (
                    dusk_local.astimezone(pytz.UTC),
                    dawn_local.astimezone(pytz.UTC),
                )
            )
        return segments

    def _exo_observable_from_site_during_nights(
        self,
        ra_deg: float,
        dec_deg: float,
        astro_night_segments_utc: List[Tuple[datetime, datetime]],
        pl_orbper: Optional[float],
        pl_tranmid_jd: Optional[float],
        pl_trandur_h: Optional[float],
    ) -> bool:
        """
        True s’il existe au moins une nuit astronomique où un transit tombe tel que
        ingress−1 h et egress+1 h sont dans [crépuscule astro, aube astro] et que
        l’altitude de la cible est > 25° à ces deux instants (NASA : pl_tranmid en JD,
        pl_trandur durée totale du transit en heures).
        """
        if (
            pl_orbper is None
            or pl_tranmid_jd is None
            or pl_trandur_h is None
            or pl_orbper <= 0
            or pl_trandur_h <= 0
        ):
            return False

        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
        half_dur_d = (pl_trandur_h / 2.0) / 24.0
        margin_d = 1.0 / 24.0

        for t_dusk_utc, t_dawn_utc in astro_night_segments_utc:
            if t_dawn_utc <= t_dusk_utc:
                continue
            jd_start = float(Time(t_dusk_utc).jd)
            jd_end = float(Time(t_dawn_utc).jd)
            jd_lo = jd_start + half_dur_d + margin_d
            jd_hi = jd_end - half_dur_d - margin_d
            if jd_lo > jd_hi:
                continue

            t0 = pl_tranmid_jd
            p = pl_orbper
            k_min = math.ceil((jd_lo - t0) / p)
            k_max = math.floor((jd_hi - t0) / p)
            if k_min > k_max:
                continue

            for k in range(k_min, k_max + 1):
                t_mid = t0 + k * p
                t1_jd = t_mid - half_dur_d - margin_d
                t2_jd = t_mid + half_dur_d + margin_d
                if t1_jd < jd_start - 1e-9 or t2_jd > jd_end + 1e-9:
                    continue
                ta1 = Time(t1_jd, format="jd")
                ta2 = Time(t2_jd, format="jd")
                alt1 = coord.transform_to(AltAz(obstime=ta1, location=self.location)).alt.deg
                alt2 = coord.transform_to(AltAz(obstime=ta2, location=self.location)).alt.deg
                if alt1 > EXO_TRANSIT_MIN_ALT_DEG and alt2 > EXO_TRANSIT_MIN_ALT_DEG:
                    if self.location is not None:
                        t_mid_ast = Time(t_mid, format="jd")
                        seps = (
                            moon_object_separation_deg_icrs(ra_deg, dec_deg, ta1, self.location),
                            moon_object_separation_deg_icrs(ra_deg, dec_deg, t_mid_ast, self.location),
                            moon_object_separation_deg_icrs(ra_deg, dec_deg, ta2, self.location),
                        )
                        if min(seps) < NIGHT_MOON_MIN_SEPARATION_DEG:
                            continue
                    return True
        return False

    def _object_clear_of_moon(self, obj: EphemerisObject) -> bool:
        """
        True si la cible respecte l'écart minimal à la Lune (50°) pour les types
        exoplanète, astéroïde, comète et transitoire. Les binaires à éclipses ne sont pas filtrées ici.
        """
        if obj.obj_type not in ("exoplanet", "asteroid", "comet", "transient"):
            return True
        if self.location is None:
            return True
        min_sep = NIGHT_MOON_MIN_SEPARATION_DEG
        if obj.obj_type == "exoplanet":
            t0 = obj.exo_pl_tranmid_jd
            p = obj.period
            if (
                t0 is not None
                and p is not None
                and p > 0
                and getattr(self, "start_time", None) is not None
                and getattr(self, "end_time", None) is not None
            ):
                st = self.start_time
                et = self.end_time
                if st.tzinfo is None:
                    st = pytz.UTC.localize(st)
                else:
                    st = st.astimezone(pytz.UTC)
                if et.tzinfo is None:
                    et = pytz.UTC.localize(et)
                else:
                    et = et.astimezone(pytz.UTC)
                jd_start = float(Time(st).jd)
                jd_end = float(Time(et).jd)
                k_min = math.ceil((jd_start - t0) / p)
                k_max = math.floor((jd_end - t0) / p)
                if k_min <= k_max:
                    t_mid_jd = t0 + k_min * p
                    times_ast = [Time(t_mid_jd, format="jd")]
                    tdur_h = obj.exo_pl_trandur_h
                    if tdur_h is not None and tdur_h > 0:
                        half_d = (tdur_h / 2.0) / 24.0
                        margin_d = 1.0 / 24.0
                        times_ast = [
                            Time(t_mid_jd - half_d - margin_d, format="jd"),
                            Time(t_mid_jd, format="jd"),
                            Time(t_mid_jd + half_d + margin_d, format="jd"),
                        ]
                    for t_ast in times_ast:
                        if (
                            moon_object_separation_deg_icrs(obj.ra, obj.dec, t_ast, self.location)
                            < min_sep
                        ):
                            return False
                    return True
        t = obj.alt_max_time
        if t is None:
            return True
        if t.tzinfo is None:
            t = pytz.UTC.localize(t)
        else:
            t = t.astimezone(pytz.UTC)
        return (
            moon_object_separation_deg_icrs(obj.ra, obj.dec, Time(t), self.location) >= min_sep
        )

    def fetch_nasa_exoplanet_tap_rows(self, dec_min: float, dec_max: float) -> List[dict]:
        """Interroge le TAP NASA avec filtre déclinaison (site observatoire)."""
        adql = build_nasa_exoplanet_tap_adql(dec_min, dec_max)
        logger.info("TAP NASA (aperçu requête) : %s", adql[:200] + ("..." if len(adql) > 200 else ""))
        response = requests.get(
            NASA_EXOARCHIVE_TAP_SYNC,
            params={"query": adql, "format": "csv"},
            timeout=180,
        )
        response.raise_for_status()
        text = response.text.strip()
        if not text:
            return []
        reader = csv.DictReader(io.StringIO(text))
        return list(reader)

    def browse_aavso_catalogues_dir(self):
        """Ouvre un dialogue pour choisir le dossier d'enregistrement d'index.csv (AAVSO)."""
        initial = self._aavso_files_directory()
        if not initial.is_dir():
            initial = self.catalogues_dir
        chosen = filedialog.askdirectory(
            title="Dossier pour index.csv (AAVSO) et nasa_exoplanet_transits.csv (NASA)",
            initialdir=str(initial),
        )
        if chosen:
            self.aavso_catalogues_dir_var.set(chosen)
            self.update_catalog_status()

    def update_catalog_status(self):
        """Met à jour le statut des catalogues avec dates de mise à jour."""
        # Statut MPC
        mpc_dir = self._mpc_files_directory()
        nea_file = mpc_dir / "NEA.txt"
        comet_file = mpc_dir / "AllCometEls.txt"
        
        if nea_file.exists() or comet_file.exists():
            # Prendre la date la plus récente
            dates = []
            if nea_file.exists():
                dates.append(datetime.fromtimestamp(nea_file.stat().st_mtime))
            if comet_file.exists():
                dates.append(datetime.fromtimestamp(comet_file.stat().st_mtime))
            if dates:
                latest_date = max(dates)
                self.mpc_status_label.config(
                    text=f"Dernière mise à jour: {latest_date.strftime('%Y-%m-%d %H:%M')}", 
                    foreground="green")
            else:
                self.mpc_status_label.config(text="Non téléchargé", foreground="gray")
        else:
            self.mpc_status_label.config(text="Non téléchargé", foreground="gray")
        
        # Statut AAVSO (binaires à éclipses)
        index_csv = self._aavso_index_csv_path()
        if index_csv.exists():
            file_size = index_csv.stat().st_size / 1024
            mtime = datetime.fromtimestamp(index_csv.stat().st_mtime)
            self.aavso_status_label.config(
                text=f"Dernière mise à jour: {mtime.strftime('%Y-%m-%d %H:%M')} ({file_size:.1f} KB)",
                foreground="green",
            )
        else:
            self.aavso_status_label.config(text="Non téléchargé", foreground="gray")

        # Statut catalogue exoplanètes NASA (TAP)
        if self.nasa_exo_status_label is not None:
            nasa_csv = self._nasa_exoplanet_csv_path()
            if nasa_csv.exists():
                file_size = nasa_csv.stat().st_size / 1024
                mtime = datetime.fromtimestamp(nasa_csv.stat().st_mtime)
                self.nasa_exo_status_label.config(
                    text=f"Dernière mise à jour: {mtime.strftime('%Y-%m-%d %H:%M')} ({file_size:.1f} KB)",
                    foreground="green",
                )
            else:
                self.nasa_exo_status_label.config(text="Non téléchargé", foreground="gray")

    def download_mpc_catalogs(self):
        """Télécharge les catalogues MPC (NEA et comètes)."""
        try:
            self.mpc_status_label.config(text="Téléchargement en cours...", foreground="blue")
            self.update()
            
            # URLs des catalogues MPC
            nea_url = "https://www.minorplanetcenter.net/iau/MPCORB/NEA.txt"
            comet_url = "https://minorplanetcenter.net/iau/MPCORB/AllCometEls.txt"
            
            mpc_dir = self._mpc_files_directory()
            mpc_dir.mkdir(parents=True, exist_ok=True)

            # Télécharger NEA
            logger.info(f"Téléchargement du catalogue NEA depuis {nea_url}...")
            nea_file = mpc_dir / "NEA.txt"
            
            try:
                response = requests.get(nea_url, timeout=120)
                if response.status_code == 200:
                    with open(nea_file, 'wb') as f:
                        f.write(response.content)
                    logger.info(f"Catalogue NEA téléchargé: {len(response.content)} bytes")
                else:
                    logger.warning(f"Erreur téléchargement NEA: HTTP {response.status_code}")
            except Exception as e:
                logger.error(f"Erreur lors du téléchargement du catalogue NEA: {e}")
            
            # Télécharger comètes
            logger.info("Téléchargement du catalogue comètes...")
            response = requests.get(comet_url, timeout=120)
            if response.status_code == 200:
                comet_file = mpc_dir / "AllCometEls.txt"
                with open(comet_file, 'wb') as f:
                    f.write(response.content)
                logger.info(f"Catalogue comètes téléchargé: {len(response.content)} bytes")
            else:
                logger.warning(f"Erreur téléchargement comètes: HTTP {response.status_code}")
            
            # Mettre à jour le statut
            self.update_catalog_status()
            messagebox.showinfo("Succès", "Catalogues MPC téléchargés avec succès")
            
        except Exception as e:
            logger.error(f"Erreur téléchargement catalogues MPC: {e}", exc_info=True)
            self.mpc_status_label.config(text="Erreur de téléchargement", foreground="red")
            messagebox.showerror("Erreur", f"Erreur lors du téléchargement des catalogues MPC:\n{e}")
    
    def download_aavso_index(self):
        """Télécharge index.csv depuis l'AAVSO Target Tool (binaires à éclipses uniquement, pas d'exoplanètes)."""
        try:
            url = "https://targettool.aavso.org/TargetTool/default/index.csv?ev=on&settype=true"

            self.aavso_status_label.config(text="Téléchargement en cours...", foreground="blue")
            self.update()

            logger.info(f"Téléchargement depuis AAVSO (EV only): {url}")
            response = requests.get(url, timeout=60)
            response.raise_for_status()

            content = response.text
            if not content or len(content.strip()) == 0:
                raise ValueError("Réponse vide du serveur")

            index_csv = self._aavso_index_csv_path()
            index_csv.parent.mkdir(parents=True, exist_ok=True)
            with open(index_csv, "w", encoding="utf-8") as f:
                f.write(content)

            logger.info(f"Catalogue AAVSO sauvegardé: {index_csv}")
            self.update_catalog_status()
            messagebox.showinfo(
                "Succès",
                "Catalogue AAVSO (binaires à éclipses) téléchargé et sauvegardé dans :\n"
                f"{index_csv}",
            )

        except Exception as e:
            logger.error(f"Erreur téléchargement AAVSO: {e}", exc_info=True)
            self.aavso_status_label.config(text="Erreur de téléchargement", foreground="red")
            messagebox.showerror("Erreur", f"Erreur lors du téléchargement depuis AAVSO:\n{e}")

    def download_nasa_exoplanet_transits_catalog(self):
        """
        Interroge le TAP NASA pour les exoplanètes en transit filtrées par le site (δ) et par
        nuit astronomique : ingress−1 h et egress+1 h à altitude > 25° (dates « Paramètres »).
        """
        try:
            obs_start = datetime.strptime(self.date_var.get().strip(), "%Y-%m-%d")
            obs_end = None
            if self.use_interval_var.get() and self.date_end_var.get().strip():
                obs_end = datetime.strptime(self.date_end_var.get().strip(), "%Y-%m-%d")
        except ValueError:
            messagebox.showerror(
                "Date",
                "Indiquez une date d’observation valide (AAAA-MM-JJ). "
                "Le filtre NASA utilise cette date (et la date de fin si « Intervalle » est coché).",
            )
            return

        try:
            if self.nasa_exo_status_label is not None:
                self.nasa_exo_status_label.config(text="Téléchargement en cours...", foreground="blue")
            self.update()

            lat_deg = self.location.lat.deg
            dec_min, dec_max = declination_bounds_for_latitude(lat_deg)
            astro_segments = self._utc_astronomical_night_segments(obs_start, obs_end)

            rows = self.fetch_nasa_exoplanet_tap_rows(dec_min, dec_max)
            fieldnames = [
                "pl_name",
                "ra",
                "dec",
                "sy_vmag",
                "sy_gaiamag",
                "sy_jmag",
                "sy_kmag",
                "pl_orbper",
                "pl_tranmid",
                "pl_trandur",
                "transitdepthcalc",
            ]
            kept: List[dict] = []
            for row in rows:
                try:
                    ra = float(row["ra"])
                    dec = float(row["dec"])
                except (KeyError, TypeError, ValueError):
                    continue
                try:
                    p_orb = float(row["pl_orbper"]) if row.get("pl_orbper") not in (None, "") else None
                    tmid = float(row["pl_tranmid"]) if row.get("pl_tranmid") not in (None, "") else None
                    tdur = float(row["pl_trandur"]) if row.get("pl_trandur") not in (None, "") else None
                except (TypeError, ValueError):
                    continue
                if not self._exo_observable_from_site_during_nights(
                    ra, dec, astro_segments, p_orb, tmid, tdur
                ):
                    continue
                tdc = transitdepthcalc_percent_from_row(row)
                rec = {k: row.get(k, "") for k in fieldnames if k != "transitdepthcalc"}
                rec["transitdepthcalc"] = str(tdc)
                kept.append(rec)

            out_path = self._nasa_exoplanet_csv_path()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
                w.writeheader()
                w.writerows(kept)

            logger.info(
                "Catalogue NASA exoplanètes (site + nuit astro + transit) : %s lignes → %s",
                len(kept),
                out_path,
            )
            self.update_catalog_status()
            messagebox.showinfo(
                "Succès",
                f"{len(kept)} exoplanète(s) en transit retenues (TAP + nuit astronomique, "
                f"ingress−1 h / egress+1 h à alt > {EXO_TRANSIT_MIN_ALT_DEG:.0f}°, lat {lat_deg:.2f}°, "
                f"nuits du {obs_start.date()}"
                f"{' au ' + str(obs_end.date()) if obs_end else ''}).\nFichier :\n{out_path}",
            )
        except Exception as e:
            logger.error(f"Erreur téléchargement NASA Exoplanet Archive: {e}", exc_info=True)
            if self.nasa_exo_status_label is not None:
                self.nasa_exo_status_label.config(text="Erreur de téléchargement", foreground="red")
            messagebox.showerror(
                "Erreur",
                f"Erreur lors du téléchargement depuis le NASA Exoplanet Archive :\n{e}",
            )
    
    def set_today(self):
        """Définit la date d'observation à aujourd'hui."""
        self.date_var.set(datetime.now().strftime("%Y-%m-%d"))
    
    def load_nasa_exoplanet_objects(
        self, obs_date: datetime, date_end: Optional[datetime] = None
    ):
        """
        Interroge le TAP NASA avec filtre déclinaison (latitude observatoire / config),
        puis ne garde que les cibles avec transit (pl_tranmid, pl_trandur, pl_orbper) vérifiables :
        ingress−1 h et egress+1 h dans la nuit astronomique et altitude > 25° à ces instants.
        """
        objects: List[EphemerisObject] = []
        lat_deg = self.location.lat.deg
        dec_min, dec_max = declination_bounds_for_latitude(lat_deg)
        astro_segments = self._utc_astronomical_night_segments(obs_date, date_end)

        try:
            rows = self.fetch_nasa_exoplanet_tap_rows(dec_min, dec_max)
        except Exception as e:
            logger.error(f"Erreur TAP NASA (exoplanètes): {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Échec de la requête TAP NASA Exoplanet Archive :\n{e}")
            return objects

        for row in rows:
            try:
                name = normalize_exoplanet_name((row.get("pl_name") or "").strip())
                if not name:
                    continue
                ra = float(row["ra"])
                dec = float(row["dec"])
            except (KeyError, TypeError, ValueError):
                continue

            try:
                p_orb = float(row["pl_orbper"]) if row.get("pl_orbper") not in (None, "") else None
                tmid = float(row["pl_tranmid"]) if row.get("pl_tranmid") not in (None, "") else None
                tdur = float(row["pl_trandur"]) if row.get("pl_trandur") not in (None, "") else None
            except (TypeError, ValueError):
                continue

            if not self._exo_observable_from_site_during_nights(
                ra, dec, astro_segments, p_orb, tmid, tdur
            ):
                continue

            magnitude = host_star_magnitude_from_nasa_row(row)

            period = None
            p_raw = row.get("pl_orbper")
            if p_raw is not None and str(p_raw).strip() != "":
                try:
                    period = float(p_raw)
                except (TypeError, ValueError):
                    pass

            tdc_pct = transitdepthcalc_percent_from_row(row)

            objects.append(
                EphemerisObject(
                    name=name,
                    ra=ra,
                    dec=dec,
                    magnitude=magnitude,
                    obj_type="exoplanet",
                    period=period,
                    exo_pl_tranmid_jd=tmid,
                    exo_pl_trandur_h=tdur,
                    exo_transitdepthcalc_pct=tdc_pct,
                    exo_source=EXOPLANET_PROVIDER_NASA,
                )
            )

        logger.info(
            "%s exoplanètes retenues après TAP (δ) + nuit astro + fenêtre transit (alt > %.0f°) lat=%.2f°",
            len(objects),
            EXO_TRANSIT_MIN_ALT_DEG,
            lat_deg,
        )
        return objects

    def _iter_enabled_exoplanet_providers(self) -> List[str]:
        """
        Retourne les fournisseurs exoplanètes activés.
        NASA reste la source de référence ; ExoClock est activé en complément
        lorsqu'il est disponible localement.
        """
        providers = []
        if EXOCLOCK_AVAILABLE:
            providers.append(EXOPLANET_PROVIDER_EXOCLOCK)
        providers.append(EXOPLANET_PROVIDER_NASA)
        if self.etd_authenticated:
            providers.append(EXOPLANET_PROVIDER_ETD)
        else:
            logger.info("Provider ETD désactivé automatiquement (authentification absente).")
        logger.info(
            "Providers exoplanètes actifs: %s (EXOCLOCK_AVAILABLE=%s, ETD_AUTH=%s)",
            ", ".join(providers),
            EXOCLOCK_AVAILABLE,
            self.etd_authenticated,
        )
        return providers

    def _probe_etd_authentication(self) -> bool:
        """Vérifie si l'API ETD est accessible sans/auth session courante."""
        try:
            api_probe = requests.get(
                f"{ETD_SEARCH_API_URL}?pageId=1&pageSize=1",
                timeout=20,
            )
            if api_probe.status_code == 401:
                return False
            return api_probe.status_code == 200
        except Exception:
            return False

    @staticmethod
    def _nested_get(data: dict, *path: str):
        """Accès défensif à un champ potentiellement imbriqué."""
        cur = data
        for key in path:
            if not isinstance(cur, dict):
                return None
            cur = cur.get(key)
            if cur is None:
                return None
        return cur

    @staticmethod
    def _to_float_or_none(value) -> Optional[float]:
        if value is None:
            return None
        try:
            txt = str(value).strip()
            if txt == "":
                return None
            return float(txt)
        except (TypeError, ValueError):
            return None

    def _normalize_exoclock_planet_record(self, rec: dict) -> Optional[dict]:
        """
        Normalise un enregistrement ExoClock vers le schéma interne.
        Tolère plusieurs structures (top-level ou sous-clés star/planet).
        """
        if not isinstance(rec, dict):
            return None

        name = (
            rec.get("name")
            or rec.get("planet_name")
            or rec.get("pl_name")
            or self._nested_get(rec, "planet", "name")
        )
        if not name:
            return None

        ra = self._to_float_or_none(
            rec.get("ra")
            or rec.get("ra_deg")
            or self._nested_get(rec, "star", "ra_deg")
            or self._nested_get(rec, "star", "ra")
        )
        dec = self._to_float_or_none(
            rec.get("dec")
            or rec.get("dec_deg")
            or self._nested_get(rec, "star", "dec_deg")
            or self._nested_get(rec, "star", "dec")
        )
        if ra is None or dec is None:
            return None

        period = self._to_float_or_none(
            rec.get("pl_orbper")
            or rec.get("ephem_period")
            or self._nested_get(rec, "planet", "ephem_period")
            or self._nested_get(rec, "planet", "period")
        )
        tmid = self._to_float_or_none(
            rec.get("pl_tranmid")
            or rec.get("ephem_mid_time")
            or self._nested_get(rec, "planet", "ephem_mid_time")
            or self._nested_get(rec, "planet", "mid_time")
        )
        tdur = self._to_float_or_none(
            rec.get("pl_trandur")
            or rec.get("transit_duration_h")
            or rec.get("transit_duration")
            or self._nested_get(rec, "planet", "transit_duration")
            or self._nested_get(rec, "planet", "duration")
        )
        if tdur is not None and tdur > 24.0:
            # Certains catalogues fournissent la durée en minutes.
            tdur = tdur / 60.0

        rp_over_rs = self._to_float_or_none(
            rec.get("pl_ratror")
            or rec.get("rp_over_rs")
            or self._nested_get(rec, "planet", "rp_over_rs")
        )
        depth_pct = None
        if rp_over_rs is not None and rp_over_rs > 0:
            depth_pct = 100.0 * rp_over_rs * rp_over_rs

        mag = self._to_float_or_none(
            rec.get("sy_vmag")
            or rec.get("sy_gaiamag")
            or rec.get("sy_jmag")
            or rec.get("vmag")
            or rec.get("v_mag")
            or rec.get("Vmag")
            or rec.get("brightness")
            or rec.get("stellar_magnitude")
            or self._nested_get(rec, "star", "vmag")
            or self._nested_get(rec, "star", "mag")
            or self._nested_get(rec, "star", "sy_vmag")
            or self._nested_get(rec, "star", "gaia_mag")
        )
        if mag is not None and (not math.isfinite(mag) or mag < -4.0 or mag > 25.0):
            mag = None

        return {
            "name": str(name).strip(),
            "ra": ra,
            "dec": dec,
            "period": period,
            "tmid": tmid,
            "tdur_h": tdur,
            "depth_pct": depth_pct,
            "mag": mag,
        }

    def load_exoclock_exoplanet_objects(
        self, obs_date: datetime, date_end: Optional[datetime] = None
    ) -> List[EphemerisObject]:
        """
        Charge les exoplanètes depuis ExoClock (si le module est disponible),
        puis applique les mêmes filtres de nuit/altitude que pour NASA.
        """
        if not EXOCLOCK_AVAILABLE:
            logger.info("Module 'exoclock' indisponible : source ExoClock ignorée.")
            return []

        objects: List[EphemerisObject] = []
        astro_segments = self._utc_astronomical_night_segments(obs_date, date_end)

        try:
            raw = exoclock.get_all_planets()
        except Exception as e:
            logger.warning("ExoClock: échec get_all_planets(): %s", e)
            return []

        records = []
        if isinstance(raw, list):
            records = raw
        elif isinstance(raw, dict):
            # Selon versions: soit dict de records, soit payload avec clé planets.
            if isinstance(raw.get("planets"), list):
                records = raw["planets"]
            else:
                records = list(raw.values())
        else:
            logger.warning("ExoClock: format inattendu pour get_all_planets(): %s", type(raw))
            return []

        parsed = 0
        retained = 0
        for rec in records:
            if isinstance(rec, str):
                # Dans plusieurs versions d'exoclock, get_all_planets()
                # renvoie simplement une liste de noms de planètes.
                try:
                    rec = exoclock.get_planet(rec)
                except Exception:
                    continue

            norm = self._normalize_exoclock_planet_record(rec)
            if norm is None:
                continue
            parsed += 1

            if not self._exo_observable_from_site_during_nights(
                norm["ra"],
                norm["dec"],
                astro_segments,
                norm["period"],
                norm["tmid"],
                norm["tdur_h"],
            ):
                continue

            objects.append(
                EphemerisObject(
                    name=normalize_exoplanet_name(norm["name"]),
                    ra=norm["ra"],
                    dec=norm["dec"],
                    magnitude=norm["mag"],
                    obj_type="exoplanet",
                    period=norm["period"],
                    exo_pl_tranmid_jd=norm["tmid"],
                    exo_pl_trandur_h=norm["tdur_h"],
                    exo_transitdepthcalc_pct=norm["depth_pct"],
                    exo_source=EXOPLANET_PROVIDER_EXOCLOCK,
                )
            )
            retained += 1

        logger.info(
            "ExoClock: %d enregistrements exploitables, %d retenus après filtres nuit/transit",
            parsed,
            retained,
        )
        return objects

    def load_etd_exoplanet_objects(
        self, obs_date: datetime, date_end: Optional[datetime] = None
    ) -> List[EphemerisObject]:
        """
        Charge des cibles depuis le catalogue ETD/VarAstro si la page retourne
        un tableau de données exploitable côté serveur.
        """
        objects: List[EphemerisObject] = []
        astro_segments = self._utc_astronomical_night_segments(obs_date, date_end)

        if not self.etd_authenticated:
            logger.info("ETD: source ignorée (non authentifiée).")
            return objects

        try:
            response = requests.get(ETD_EXOPLANETS_URL, timeout=60)
            response.raise_for_status()
            html = response.text
        except Exception as e:
            logger.warning("ETD: échec d'accès à %s : %s", ETD_EXOPLANETS_URL, e)
            return objects

        try:
            tables = pd.read_html(io.StringIO(html))
        except Exception:
            tables = []
        if not tables:
            logger.info("ETD: aucun tableau HTML détecté sur %s", ETD_EXOPLANETS_URL)
            return objects

        etd_table = None
        for t in tables:
            cols = {str(c).strip().lower() for c in t.columns}
            if "nom" in cols and "ascension droite" in cols and "déclinaison" in cols:
                etd_table = t
                break
        if etd_table is None:
            logger.info("ETD: aucun tableau compatible trouvé (colonnes attendues absentes).")
            return objects

        parsed = 0
        retained = 0
        for _, row in etd_table.iterrows():
            try:
                name = normalize_exoplanet_name(str(row.get("Nom", "")).strip())
                if not name or name.lower() == "nan":
                    continue
                ra_raw = str(row.get("Ascension droite", "")).strip()
                dec_raw = str(row.get("Déclinaison", "")).strip()
                ra, dec = parse_ra_dec(ra_raw, dec_raw)
                if ra is None or dec is None:
                    continue

                mag = self._to_float_or_none(row.get("Luminosité (Mag)"))
                if mag is not None and (not math.isfinite(mag) or mag < -4.0 or mag > 25.0):
                    mag = None

                period = self._to_float_or_none(row.get("Période"))
                tmid = self._to_float_or_none(row.get("Époch"))
                dur_min = self._to_float_or_none(row.get("Durée (m)"))
                tdur_h = (dur_min / 60.0) if dur_min is not None else None
                depth_pct = self._to_float_or_none(row.get("Profondeur (%)"))
            except Exception:
                continue

            parsed += 1
            if not self._exo_observable_from_site_during_nights(
                ra, dec, astro_segments, period, tmid, tdur_h
            ):
                continue

            objects.append(
                EphemerisObject(
                    name=name,
                    ra=ra,
                    dec=dec,
                    magnitude=mag,
                    obj_type="exoplanet",
                    period=period,
                    exo_pl_tranmid_jd=tmid,
                    exo_pl_trandur_h=tdur_h,
                    exo_transitdepthcalc_pct=depth_pct,
                    exo_source=EXOPLANET_PROVIDER_ETD,
                )
            )
            retained += 1

        logger.info(
            "ETD: %d enregistrements exploitables, %d retenus après filtres nuit/transit",
            parsed,
            retained,
        )
        return objects

    def load_exoplanet_objects(
        self, obs_date: datetime, date_end: Optional[datetime] = None
    ) -> List[EphemerisObject]:
        """
        Agrège les exoplanètes depuis plusieurs fournisseurs et dédoublonne
        les cibles par nom (insensible à la casse).
        """
        providers = self._iter_enabled_exoplanet_providers()
        all_objects: List[EphemerisObject] = []
        provider_counts = {}

        for provider in providers:
            if provider == EXOPLANET_PROVIDER_NASA:
                provider_objects = self.load_nasa_exoplanet_objects(obs_date, date_end)
            elif provider == EXOPLANET_PROVIDER_EXOCLOCK:
                provider_objects = self.load_exoclock_exoplanet_objects(obs_date, date_end)
            elif provider == EXOPLANET_PROVIDER_ETD:
                provider_objects = self.load_etd_exoplanet_objects(obs_date, date_end)
            else:
                logger.warning("Fournisseur exoplanètes inconnu ignoré : %s", provider)
                continue

            logger.info(
                "Source exoplanètes '%s' : %d cible(s) retenue(s)",
                provider,
                len(provider_objects),
            )
            provider_counts[provider] = len(provider_objects)
            all_objects.extend(provider_objects)

        deduped: List[EphemerisObject] = []
        seen_names = {}

        def _canonical_exoplanet_key(name: str) -> str:
            """
            Clé de dédoublonnage robuste pour exoplanètes :
            ignore casse, espaces, tirets et ponctuation.
            Exemples équivalents: 'WASP-12 b', 'wasp 12b', 'WASP12B'.
            """
            raw = (name or "").strip().lower()
            return re.sub(r"[^a-z0-9]", "", raw)

        for obj in all_objects:
            key = _canonical_exoplanet_key(obj.name)
            if not key:
                continue

            if key not in seen_names:
                seen_names[key] = obj
                deduped.append(obj)
                continue

            # Fusion douce: la première source garde la priorité,
            # mais on complète les champs manquants avec les sources suivantes.
            # Magnitude : préférer NASA (sy_vmag / catalogue TAP) si disponible.
            kept = seen_names[key]
            obj_src = (getattr(obj, "exo_source", None) or "").strip().lower()
            kept_src = (getattr(kept, "exo_source", None) or "").strip().lower()
            if obj.magnitude is not None:
                if kept.magnitude is None:
                    kept.magnitude = obj.magnitude
                elif obj_src == EXOPLANET_PROVIDER_NASA and kept_src != EXOPLANET_PROVIDER_NASA:
                    kept.magnitude = obj.magnitude
            if kept.period is None and obj.period is not None:
                kept.period = obj.period
            if kept.exo_pl_tranmid_jd is None and obj.exo_pl_tranmid_jd is not None:
                kept.exo_pl_tranmid_jd = obj.exo_pl_tranmid_jd
            if kept.exo_pl_trandur_h is None and obj.exo_pl_trandur_h is not None:
                kept.exo_pl_trandur_h = obj.exo_pl_trandur_h
            if (
                kept.exo_transitdepthcalc_pct is None
                and obj.exo_transitdepthcalc_pct is not None
            ):
                kept.exo_transitdepthcalc_pct = obj.exo_transitdepthcalc_pct

        # Exiger une magnitude hôte connue (NASA et/ou ExoClock / ETD après fusion).
        with_mag: List[EphemerisObject] = []
        dropped_mag = 0
        for obj in deduped:
            if obj.magnitude is None:
                dropped_mag += 1
                logger.warning(
                    "Exoplanète '%s' exclue : magnitude de l’étoile hôte introuvable "
                    "(NASA : sy_vmag, sy_gaiamag, … ; ExoClock : sy_vmag, vmag, …).",
                    obj.name,
                )
                continue
            with_mag.append(obj)

        if dropped_mag:
            logger.info(
                "Exoplanètes sans magnitude après fusion NASA/ExoClock/ETD : %d exclue(s)",
                dropped_mag,
            )

        logger.info(
            "Exoplanètes agrégées : %d entrée(s) (%d après dédoublonnage, %d avec mag) via %s",
            len(all_objects),
            len(deduped),
            len(with_mag),
            ", ".join(providers) if providers else "aucune source",
        )
        self.last_exoplanet_provider_counts = {
            **provider_counts,
            "total_before_dedup": len(all_objects),
            "total_after_dedup": len(deduped),
            "total_with_magnitude": len(with_mag),
            "dropped_missing_magnitude": dropped_mag,
        }
        return with_mag

    def load_aavso_objects(self):
        """Charge les étoiles binaires à éclipses depuis index.csv AAVSO (types EB / EA)."""
        objects = []

        index_path = self._aavso_index_csv_path()
        # Ancien emplacement à la racine .npoap
        old_path = self.npoap_dir / "index.csv"
        if not index_path.exists() and old_path.exists():
            try:
                import shutil

                index_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(old_path), str(index_path))
                logger.info(f"Fichier déplacé de {old_path} vers {index_path}")
            except Exception as e:
                logger.warning(f"Impossible de déplacer le fichier: {e}")
                self.aavso_catalogues_dir_var.set(str(old_path.parent))
                index_path = self._aavso_index_csv_path()

        if not index_path.exists():
            logger.warning("Catalogue AAVSO (binaires) non trouvé. Utilisez le bouton de téléchargement.")
            return []

        try:
            with open(index_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                
                if reader.fieldnames:
                    logger.info(f"Colonnes AAVSO disponibles: {list(reader.fieldnames)}")
                
                for row_num, row in enumerate(reader, start=2):
                    try:
                        var_type = row.get('Var. Type', row.get('Var.Type', '')).strip()

                        if 'EB' not in var_type and 'EA' not in var_type:
                            continue
                        
                        name = row.get('Star Name', '').strip()
                        if not name:
                            continue
                        
                        # Chercher RA et DEC avec plusieurs noms possibles
                        ra_str = row.get('RA (J2000.0)', row.get('RA(J2000.0)', row.get('RA', '')))
                        dec_str = row.get('Dec (J2000.0)', row.get('Dec(J2000.0)', row.get('Dec', row.get('DEC', ''))))
                        
                        if not ra_str or not dec_str:
                            logger.debug(f"Ligne {row_num}: RA ou DEC manquant pour {name}")
                            continue
                        
                        ra, dec = parse_ra_dec(ra_str, dec_str)
                        if ra is None or dec is None:
                            logger.warning(f"Ligne {row_num}: Impossible de parser RA/DEC pour {name}: RA={ra_str}, DEC={dec_str}")
                            continue
                        
                        # Magnitude
                        mag_str = row.get('Max Mag', row.get('MaxMag', row.get('Magnitude', '')))
                        magnitude = None
                        if mag_str:
                            try:
                                # Extraire la valeur numérique (peut être "9.52 V" ou "None V")
                                mag_parts = mag_str.strip().split()
                                if mag_parts and mag_parts[0].lower() != 'none':
                                    mag_val = float(mag_parts[0])
                                    if mag_val <= 15.0:  # Filtrer par magnitude
                                        magnitude = mag_val
                                    else:
                                        logger.debug(f"{name} exclu: magnitude {mag_val} > 15.0")
                                        continue
                                else:
                                    # Magnitude "None", on continue quand même
                                    logger.debug(f"{name}: magnitude None, on continue")
                            except (ValueError, IndexError):
                                pass
                        
                        obj = EphemerisObject(
                            name=name,
                            ra=ra,
                            dec=dec,
                            magnitude=magnitude,
                            obj_type="eclipsing_binary",
                            period=None
                        )
                        objects.append(obj)
                        
                    except Exception as e:
                        logger.debug(f"Erreur ligne {row_num}: {e}")
                        continue
            
            logger.info(f"{len(objects)} binaires à éclipses chargés depuis l'AAVSO")
            
        except Exception as e:
            logger.error(f"Erreur lors de la lecture du catalogue: {e}", exc_info=True)
            messagebox.showerror("Erreur", 
                               f"Erreur lors de la lecture du catalogue:\n{e}")
        
        return objects
    
    def fetch_observable_mpc_objects(self, obs_date: datetime):
        """Parse les catalogues MPC locaux et interroge le service MPC pour obtenir les éphémérides."""
        objects = []
        
        try:
            mpc_dir = self._mpc_files_directory()
            nea_file = mpc_dir / "NEA.txt"
            comet_file = mpc_dir / "AllCometEls.txt"

            # Vérifier que les catalogues existent
            if not nea_file.exists() and not comet_file.exists():
                logger.warning("Catalogues MPC non téléchargés. Utilisez le bouton de téléchargement.")
                return objects
            
            # Parser les astéroïdes NEA si demandé
            if self.asteroids_var.get() and nea_file.exists():
                logger.info("Analyse des astéroïdes NEA...")
                nea_objects = self.parse_mpc_orbital_elements(nea_file, obs_date, obj_type="asteroid")
                objects.extend(nea_objects)
                logger.info("%s astéroïdes NEA retenus après éphémérides MPC", len(nea_objects))
            
            # Parser les comètes si demandé
            if self.comets_var.get() and comet_file.exists():
                logger.info("Analyse des comètes...")
                comet_objects = self.parse_mpc_orbital_elements(comet_file, obs_date, obj_type='comet')
                objects.extend(comet_objects)
                logger.info(f"{len(comet_objects)} comètes analysées")
            
        except Exception as e:
            logger.error(f"Erreur récupération catalogues MPC: {e}", exc_info=True)
            raise
        
        return objects
    
    def parse_mpc_orbital_elements(self, file_path: Path, obs_date: datetime, obj_type: str = 'asteroid'):
        """Parse un fichier MPCORB et calcule les éphémérides pour la date donnée."""
        objects = []
        nea_lines_scanned = 0
        nea_len_lt_202 = 0
        nea_h_unreadable = 0
        nea_h_gt_185 = 0
        nea_mpc_attempts = 0
        nea_mpc_fail = Counter()
        comet_line_candidates = 0
        comet_mpc_attempts = 0
        comet_mpc_fail = Counter()
        comet_no_designation = 0

        try:
            with open(file_path, 'r', encoding='latin-1') as f:
                lines = f.readlines()

            # NEA.txt : 3 lignes d'en-tête MPC. AllCometEls.txt : pas d'en-tête à sauter de la même façon.
            header_skip = 3 if obj_type == "asteroid" else 0

            for line_num, line in enumerate(lines, 1):
                if line_num <= header_skip or not line.strip():
                    continue

                try:
                    if obj_type == 'asteroid':
                        nea_lines_scanned += 1
                        # Format NEA/MPCORB: colonnes fixes
                        if len(line) < 202:
                            nea_len_lt_202 += 1
                            continue

                        # Colonnes MPCORB 1–7 : désignation permanente ou provisoire au format « packed »
                        # (voir https://www.minorplanetcenter.net/iau/info/PackedDes.html ), pas le nom propre.
                        designation = line[0:7].strip()
                        H = line[8:13].strip()  # Magnitude absolue

                        try:
                            H_val = float(H) if H else 99.0
                        except ValueError:
                            nea_h_unreadable += 1
                            continue

                        if H_val > 18.5:
                            nea_h_gt_185 += 1
                            continue

                        nea_mpc_attempts += 1
                        try:
                            mpc_obj, fail_reason = self.get_mpc_ephemeris_for_object(
                                designation, obs_date, obj_type="asteroid"
                            )
                            if mpc_obj:
                                objects.append(mpc_obj)
                                if len(objects) % 10 == 0:
                                    logger.info("Analyse NEA : %s objets avec éphéméride OK…", len(objects))
                            else:
                                nea_mpc_fail[fail_reason] += 1
                        except Exception as e:
                            logger.debug("Erreur éphémérides pour %s : %s", designation, e)
                            nea_mpc_fail[f"exception_{type(e).__name__}"] += 1

                        if len(objects) >= 200:
                            logger.info("Limite de 200 astéroïdes NEA atteinte (éphémérides MPC).")
                            break

                    elif obj_type == 'comet':
                        raw = line.rstrip("\n\r")
                        if len(raw) < 30:
                            continue

                        designation = comet_mpes_designation_from_mpc_line(raw)
                        if not designation:
                            comet_no_designation += 1
                            continue

                        comet_line_candidates += 1
                        comet_mpc_attempts += 1
                        try:
                            mpc_obj, fail_reason = self.get_mpc_ephemeris_for_object(
                                designation, obs_date, obj_type="comet"
                            )
                            if mpc_obj:
                                objects.append(mpc_obj)
                            else:
                                comet_mpc_fail[fail_reason] += 1
                        except Exception as e:
                            logger.debug("Erreur éphémérides comète %s : %s", designation, e)
                            comet_mpc_fail[f"exception_{type(e).__name__}"] += 1

                        if len(objects) >= 100:
                            logger.info("Limite de 100 comètes atteinte.")
                            break

                except Exception as e:
                    logger.debug("Erreur parsing ligne %s : %s", line_num, e)
                    continue

            if obj_type == "asteroid":
                nea_h_ok = nea_mpc_attempts
                logger.info(
                    "NEA diagnostic : lignes scannées (hors en-tête)=%s | len<202=%s | H illisible=%s | "
                    "H>18.5=%s | requêtes MPC (H OK)=%s | succès éphéméride=%s | échecs : %s",
                    nea_lines_scanned,
                    nea_len_lt_202,
                    nea_h_unreadable,
                    nea_h_gt_185,
                    nea_mpc_attempts,
                    len(objects),
                    dict(nea_mpc_fail) if nea_mpc_fail else "(aucun — toutes OK ou aucune requête)",
                )
            elif obj_type == "comet":
                logger.info(
                    "Comètes diagnostic : désignation introuvable (ligne ignorée)=%s | "
                    "candidates=%s | requêtes MPC=%s | succès=%s | échecs : %s",
                    comet_no_designation,
                    comet_line_candidates,
                    comet_mpc_attempts,
                    len(objects),
                    dict(comet_mpc_fail) if comet_mpc_fail else "(aucun)",
                )

        except Exception as e:
            logger.error("Erreur parsing fichier MPC : %s", e, exc_info=True)

        return objects
    
    def get_mpc_ephemeris_for_object(
        self, designation: str, obs_date: datetime, obj_type: str = "asteroid"
    ) -> Tuple[Optional[EphemerisObject], str]:
        """
        Récupère les éphémérides MPC. Retourne (objet, 'ok') ou (None, code) pour diagnostic
        (http_NNN, timeout, request_*, parse_no_coordinates, parse_magnitude_filter).
        """
        try:
            # MPES (même schéma que astroquery.mpc.MPCClass._args_to_ephemeris_payload) :
            # hôte cgi.minorplanetcenter.net, champs long/lat/alt (m), date « d » sans « : »,
            # intervalle i / u / l, raty/s/m pour type d'éphéméride et mouvement propre.
            base_url = "https://cgi.minorplanetcenter.net/cgi-bin/mpeph2.cgi"

            geod = self.location.geodetic
            lon_deg = geod[0].deg
            lat_deg = geod[1].deg
            alt_m = geod[2].to(u.m).value

            t0 = Time(obs_date, scale="utc", precision=0)
            d_str = t0.iso.replace(":", "")

            params = {
                "ty": "e",
                "TextArea": designation.strip(),
                "uto": "0",
                "igd": "n",
                "ibh": "n",
                "fp": "y",
                "adir": "N",
                "tit": "",
                "bu": "",
                "long": lon_deg,
                "lat": lat_deg,
                "alt": alt_m,
                "d": d_str,
                "i": "1",
                "u": "d",
                "l": "1",
                "raty": "a",
                "s": "t",
                "m": "h",
            }

            response = requests.post(base_url, data=params, timeout=10)

            if response.status_code != 200:
                return None, f"http_{response.status_code}"

            return self.parse_mpc_response(response.text, designation, obj_type)

        except requests.Timeout:
            return None, "timeout"
        except requests.RequestException as e:
            return None, f"request_{type(e).__name__}"
        except Exception as e:
            logger.debug("Erreur récupération éphémérides MPC pour %s : %s", designation, e)
            return None, f"request_{type(e).__name__}"

    def parse_mpc_response(
        self, content: str, designation: str, obj_type: str
    ) -> Tuple[Optional[EphemerisObject], str]:
        """Parse la réponse du service MPC. Retourne (EphemerisObject, 'ok') ou (None, raison)."""
        ra = None
        dec = None
        mag = None
        
        # Utiliser le parsing amélioré
        lines = content.strip().split('\n')
        
        for line in lines:
            line_stripped = line.strip()
            
            if len(line_stripped) > 50:
                date_pattern = r'\d{4}\s+\d{1,2}\s+\d{1,2}'
                time_pattern = r'\d{1,2}\s+\d{2}'
                
                if re.match(date_pattern, line_stripped) or re.match(time_pattern, line_stripped):
                    parts = line_stripped.split()
                    
                    if len(parts) >= 8:
                        try:
                            idx = 0
                            
                            if len(parts) >= 3 and parts[0].isdigit() and len(parts[0]) == 4:
                                idx = 3
                            # Heure UT MPES : souvent un seul champ HHMMSS (ex. 000000)
                            if (
                                idx < len(parts)
                                and parts[idx].isdigit()
                                and len(parts[idx]) == 6
                            ):
                                idx += 1
                            elif idx < len(parts) and parts[idx].isdigit() and len(parts[idx]) <= 2:
                                idx += 2
                            
                            if idx + 2 < len(parts):
                                ra_h = float(parts[idx])
                                ra_m = float(parts[idx + 1])
                                ra_s = float(parts[idx + 2])
                                ra = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0
                            
                            if idx + 5 < len(parts):
                                dec_str = parts[idx + 3]
                                dec_sign = 1
                                if dec_str.startswith('-') or dec_str.startswith('−'):
                                    dec_sign = -1
                                    dec_str = dec_str.lstrip('-').lstrip('−')
                                
                                dec_d = float(dec_str) * dec_sign
                                dec_m = float(parts[idx + 4])
                                dec_s = float(parts[idx + 5])
                                dec = dec_d + (dec_sign * dec_m/60.0) + (dec_sign * dec_s/3600.0)
                            
                            for j in range(min(idx + 6, len(parts)), len(parts)):
                                try:
                                    val = float(parts[j])
                                    if 8.0 <= val <= 25.0:
                                        mag = val
                                        break
                                except (ValueError, IndexError):
                                    continue
                            
                            if ra is not None and dec is not None:
                                break
                                
                        except (ValueError, IndexError):
                            continue
        
        # Utiliser regex si le parsing précédent n'a pas fonctionné
        if ra is None or dec is None:
            ra_patterns = [
                r'(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}(?:\.\d+)?)\s*[^\d]',
                r'(\d{1,2})h\s+(\d{1,2})m\s+(\d{1,2}(?:\.\d+)?)s',
            ]
            
            for pattern in ra_patterns:
                ra_match = re.search(pattern, content, re.IGNORECASE)
                if ra_match:
                    ra_h = float(ra_match.group(1))
                    ra_m = float(ra_match.group(2))
                    ra_s = float(ra_match.group(3))
                    ra = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0
                    break
            
            dec_patterns = [
                r'([+-]?)\s*(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}(?:\.\d+)?)\s*[^\d]',
                r'([+-]?)\s*(\d{1,2})°\s+(\d{1,2})[\'\s]\s*(\d{1,2}(?:\.\d+)?)',
            ]
            
            for pattern in dec_patterns:
                dec_match = re.search(pattern, content, re.IGNORECASE)
                if dec_match:
                    sign = -1 if dec_match.group(1) == '-' else 1
                    dec_d = float(dec_match.group(2))
                    dec_m = float(dec_match.group(3))
                    dec_s = float(dec_match.group(4))
                    dec = sign * (dec_d + dec_m/60.0 + dec_s/3600.0)
                    break
            
            mag_patterns = [
                r'V\s*[=:]\s*(\d+\.\d+)',
                r'Mag[^:]*[=:]\s*(\d+\.\d+)',
                r'(\d{1,2}\.\d)\s*(?:mag|V)',
            ]
            
            for pattern in mag_patterns:
                mag_match = re.search(pattern, content, re.IGNORECASE)
                if mag_match:
                    try:
                        mag_val = float(mag_match.group(1))
                        if 8.0 <= mag_val <= 25.0:
                            mag = mag_val
                            break
                    except (ValueError, IndexError):
                        continue
        
        if ra is None or dec is None:
            return None, "parse_no_coordinates"

        # Filtrer par magnitude (max 18 pour les astéroïdes, un peu plus pour les comètes)
        max_mag = 18.0 if obj_type == 'asteroid' else 20.0
        if mag and mag > max_mag:
            return None, "parse_magnitude_filter"

        return (
            EphemerisObject(
                name=designation,
                ra=ra,
                dec=dec,
                magnitude=mag,
                obj_type=obj_type,
                period=None,
            ),
            "ok",
        )
    
    def _validate_catalogues_for_ephemerides(self):
        """
        Vérifie que les fichiers catalogue nécessaires existent pour les types cochés.
        Retourne None si OK, sinon un message d'erreur détaillé.
        """
        parts = []
        mpc_dir = self._mpc_files_directory()
        if self.mpc_catalogues_dir_var.get().strip() and not mpc_dir.is_dir():
            parts.append(f"Dossier MPC inexistant ou invalide :\n{mpc_dir}")
        if self.aavso_catalogues_dir_var.get().strip() and not self._aavso_files_directory().is_dir():
            parts.append(f"Dossier AAVSO inexistant ou invalide :\n{self._aavso_files_directory()}")

        if self.eclipsing_binaries_var.get():
            idx = self._aavso_index_csv_path()
            if not idx.exists():
                parts.append(
                    "Binaires à éclipses : index.csv AAVSO introuvable.\n"
                    f"Attendu : {idx}\n"
                    "Indiquez le dossier ou téléchargez « AAVSO (binaires à éclipses) »."
                )

        if self.exoplanets_var.get():
            try:
                datetime.strptime(self.date_var.get().strip(), "%Y-%m-%d")
            except ValueError:
                parts.append(
                    "Exoplanètes : date d’observation invalide (AAAA-MM-JJ) — "
                    "requis pour interroger les sources exoplanètes avec les nuits correspondantes."
                )
            if self.use_interval_var.get() and self.date_end_var.get().strip():
                try:
                    datetime.strptime(self.date_end_var.get().strip(), "%Y-%m-%d")
                except ValueError:
                    parts.append(
                        "Exoplanètes : date de fin d’intervalle invalide (AAAA-MM-JJ)."
                    )

        if self.asteroids_var.get():
            nea = mpc_dir / "NEA.txt"
            if not nea.exists():
                parts.append(
                    "Astéroïdes : NEA.txt introuvable.\n"
                    f"Attendu : {nea}\n"
                    "Indiquez le dossier MPC ou téléchargez les fichiers MPC."
                )

        if self.comets_var.get():
            com = mpc_dir / "AllCometEls.txt"
            if not com.exists():
                parts.append(
                    "Comètes : AllCometEls.txt introuvable.\n"
                    f"Attendu : {com}\n"
                    "Indiquez le dossier MPC ou téléchargez les fichiers MPC."
                )

        if parts:
            return "\n\n".join(parts)
        return None

    def calculate_ephemerides(self):
        """Calcule les éphémérides pour la nuit sélectionnée."""
        try:
            # Parser la date
            obs_date_str = self.date_var.get()
            obs_date = datetime.strptime(obs_date_str, "%Y-%m-%d")
            self.obs_date = obs_date
            logger.info(
                "Calcul éphémérides lancé (types: ast=%s, exo=%s, com=%s, eb=%s)",
                self.asteroids_var.get(),
                self.exoplanets_var.get(),
                self.comets_var.get(),
                self.eclipsing_binaries_var.get(),
            )
            self.etd_authenticated = self._probe_etd_authentication()
            if not self.etd_authenticated:
                logger.info(
                    "ETD auto-désactivé: API non authentifiée (%s).",
                    ETD_SEARCH_API_URL,
                )
            
            # Vérifier les types d'objets
            if not (self.asteroids_var.get() or self.exoplanets_var.get() or 
                    self.comets_var.get() or self.eclipsing_binaries_var.get()):
                messagebox.showwarning("Attention", "Sélectionnez au moins un type d'objet")
                return

            missing = self._validate_catalogues_for_ephemerides()
            if missing:
                messagebox.showerror("Catalogues manquants ou dossiers invalides", missing)
                return
            
            # Initialiser la liste des objets
            self.objects = []
            self.selected_objects = {}
            
            # Charger les objets selon les types sélectionnés
            if self.eclipsing_binaries_var.get():
                self.objects.extend(self.load_aavso_objects())
            if self.exoplanets_var.get():
                date_end = None
                if self.use_interval_var.get() and self.date_end_var.get().strip():
                    date_end = datetime.strptime(self.date_end_var.get().strip(), "%Y-%m-%d")
                self.objects.extend(self.load_exoplanet_objects(obs_date, date_end))
            
            if self.asteroids_var.get() or self.comets_var.get():
                # Parser et interroger le service MPC pour les astéroïdes/comètes
                mpc_objects = self.fetch_observable_mpc_objects(obs_date)
                self.objects.extend(mpc_objects)
            
            if not self.objects:
                messagebox.showwarning("Aucune cible", "Aucune cible trouvée pour les critères sélectionnés")
                return
            
            logger.info(f"{len(self.objects)} objets au total")
            
            # Calculer les éphémérides
            logger.info("Calcul des éphémérides...")
            
            self._setup_night_window_for_plot(obs_date)

            logger.info(
                "Horizon 0° (local): coucher %s, lever %s",
                self.sunset_local.strftime("%H:%M:%S"),
                self.sunrise_local.strftime("%H:%M:%S"),
            )
            logger.info(
                "Nuit astronomique −18° (local): fin crépuscule %s, début aube %s",
                self.astro_dusk_local.strftime("%H:%M:%S"),
                self.astro_dawn_local.strftime("%H:%M:%S"),
            )

            for obj in self.objects:
                self.night_eph.calculate_ephemeris(obj, obs_date, self.start_time, self.end_time)
            
            # Filtrer selon les critères : altitude minimum de 25° pour tous les objets
            filtered_objects = []
            for obj in self.objects:
                # Vérifier que l'objet atteint au moins 25° d'altitude pendant la nuit
                if obj.alt_max < 25.0:
                    continue
                
                # Vérifier que l'objet est visible pendant au moins une partie de la nuit astro
                if obj.transit_time:
                    # S'assurer que transit_time et start_time/end_time sont comparables
                    transit_time = obj.transit_time
                    if transit_time.tzinfo is None:
                        transit_time = pytz.UTC.localize(transit_time)
                    else:
                        transit_time = transit_time.astimezone(pytz.UTC)
                    
                    if transit_time < self.start_time or transit_time > self.end_time:
                        continue
                
                if not self._object_clear_of_moon(obj):
                    continue

                filtered_objects.append(obj)
            
            self.filtered_objects = filtered_objects
            
            logger.info(f"{len(self.filtered_objects)} objets après filtrage")
            
            # Afficher les résultats (on passe les heures locales pour compatibilité)
            self.display_results(
                self.sunset_local, self.sunset_local, self.sunrise_local, self.sunrise_local
            )
            self.plot_ephemerides(obs_date, self.start_time, self.end_time)
            
            logger.info("Calcul terminé avec succès")
            details = ""
            if self.exoplanets_var.get() and self.last_exoplanet_provider_counts:
                c = self.last_exoplanet_provider_counts
                details = (
                    "\n\nDétail exoplanètes (sources): "
                    f"NASA={c.get(EXOPLANET_PROVIDER_NASA, 0)}, "
                    f"ExoClock={c.get(EXOPLANET_PROVIDER_EXOCLOCK, 0)}, "
                    f"ETD={c.get(EXOPLANET_PROVIDER_ETD, 0)}, "
                    f"avant dédoublonnage={c.get('total_before_dedup', 0)}, "
                    f"après dédoublonnage={c.get('total_after_dedup', 0)}"
                )
            messagebox.showinfo(
                "Succès",
                f"Éphémérides calculées pour {len(self.filtered_objects)} objets observables{details}",
            )
            
        except ValueError as e:
            messagebox.showerror("Erreur de date", f"Format de date invalide: {e}")
        except Exception as e:
            logger.error(f"Erreur lors du calcul: {e}", exc_info=True)
            messagebox.showerror("Erreur", f"Erreur lors du calcul: {e}")
    
    def ra_dec_to_nina_format(self, ra_deg, dec_deg):
        """Convertit les coordonnées RA/DEC en degrés vers le format NINA (heures/min/sec pour RA, deg/min/sec pour DEC)."""
        # Conversion RA : degrés -> heures (1h = 15°)
        ra_hours = ra_deg / 15.0
        ra_h = int(ra_hours)
        ra_minutes_float = (ra_hours - ra_h) * 60.0
        ra_m = int(ra_minutes_float)
        ra_s = (ra_minutes_float - ra_m) * 60.0
        
        # Conversion DEC : degrés -> degrés/minutes/secondes
        dec_abs = abs(dec_deg)
        dec_negative = dec_deg < 0
        dec_d = int(dec_abs)
        dec_minutes_float = (dec_abs - dec_d) * 60.0
        dec_m = int(dec_minutes_float)
        dec_s = (dec_minutes_float - dec_m) * 60.0
        
        return {
            'RAHours': ra_h,
            'RAMinutes': ra_m,
            'RASeconds': ra_s,
            'NegativeDec': dec_negative,
            'DecDegrees': dec_d,
            'DecMinutes': dec_m,
            'DecSeconds': dec_s
        }

    def _parse_manual_icrs_deg_for_export(self) -> Optional[Tuple[float, float]]:
        """RA/δ saisis (h m s, ° ′ ″) ; utilisé si la case « Exporter vers NINA » est cochée."""
        return parse_icrs_sexagesimal_to_deg(
            self.manual_ra_h_var.get(),
            self.manual_ra_m_var.get(),
            self.manual_ra_s_var.get(),
            self.manual_dec_d_var.get(),
            self.manual_dec_m_var.get(),
            self.manual_dec_s_var.get(),
            self.manual_dec_south_var.get(),
        )

    def _build_nina_deep_sky_json(self, target_name: str, ra_deg: float, dec_deg: float) -> dict:
        """Structure JSON NINA pour une cible (nom + coordonnées ICRS en degrés)."""
        coords = self.ra_dec_to_nina_format(ra_deg, dec_deg)
        return {
            "$id": "1",
            "$type": "NINA.Sequencer.Container.DeepSkyObjectContainer, NINA.Sequencer",
            "Target": {
                "$id": "2",
                "$type": "NINA.Astrometry.InputTarget, NINA.Astrometry",
                "Expanded": True,
                "TargetName": target_name,
                "PositionAngle": 0.0,
                "InputCoordinates": {
                    "$id": "3",
                    "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                    "RAHours": coords["RAHours"],
                    "RAMinutes": coords["RAMinutes"],
                    "RASeconds": coords["RASeconds"],
                    "NegativeDec": coords["NegativeDec"],
                    "DecDegrees": coords["DecDegrees"],
                    "DecMinutes": coords["DecMinutes"],
                    "DecSeconds": coords["DecSeconds"],
                },
            },
            "ExposureInfoListExpanded": False,
            "ExposureInfoList": {
                "$id": "4",
                "$type": "NINA.Core.Utility.AsyncObservableCollection`1[[NINA.Sequencer.Utility.ExposureInfo, NINA.Sequencer]], NINA.Core",
                "$values": [],
            },
            "Strategy": {
                "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
            },
            "Name": target_name,
            "Conditions": {
                "$id": "5",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System.Collections.ObjectModel",
                "$values": [],
            },
            "IsExpanded": True,
            "Items": {
                "$id": "6",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System.Collections.ObjectModel",
                "$values": [],
            },
            "Triggers": {
                "$id": "7",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System.Collections.ObjectModel",
                "$values": [],
            },
            "Parent": None,
            "ErrorBehavior": 0,
            "Attempts": 1,
        }

    def export_nina_json(self):
        """Exporte les objets sélectionnés ; exporte aussi les ICRS saisis si « Exporter vers NINA » est coché."""
        selected_objs = [
            obj
            for obj in self.filtered_objects
            if obj.name in self.selected_objects and self.selected_objects[obj.name]
        ]
        want_manual_nina = self.manual_export_nina_var.get()
        manual_deg = self._parse_manual_icrs_deg_for_export() if want_manual_nina else None

        if not selected_objs and manual_deg is None:
            if want_manual_nina:
                messagebox.showwarning(
                    "Export NINA",
                    "Coordonnées ICRS invalides ou incomplètes (h, m, s et °, ′, ″), "
                    "ou bien sélectionnez des objets dans la liste (case ☑).",
                )
            else:
                messagebox.showwarning(
                    "Aucun objet sélectionné",
                    "Sélectionnez au moins un objet dans la liste (case ☑), "
                    "ou cochez « Exporter vers NINA » avec des coordonnées ICRS valides.",
                )
            return

        default_dir = str(Path.home() / "Downloads")
        if not Path(default_dir).exists():
            default_dir = str(Path.home())

        dest_dir = filedialog.askdirectory(
            title="Choisir le dossier de destination pour les fichiers JSON NINA",
            initialdir=default_dir,
        )
        if not dest_dir:
            return

        dest_path = Path(dest_dir)
        exported_count = 0

        try:
            for obj in selected_objs:
                nina_json = self._build_nina_deep_sky_json(obj.name, obj.ra, obj.dec)
                safe_name = re.sub(r'[<>:"/\\|?*]', "_", obj.name)
                json_file = dest_path / f"{safe_name}.json"
                with open(json_file, "w", encoding="utf-8") as f:
                    json.dump(nina_json, f, indent=2, ensure_ascii=False)
                exported_count += 1
                logger.info("Fichier JSON exporté : %s", json_file)

            if want_manual_nina and manual_deg is None and selected_objs:
                messagebox.showwarning(
                    "Export NINA",
                    "Les coordonnées ICRS n'ont pas été exportées : valeurs invalides ou incomplètes "
                    "(vérifiez h m s et ° ′ ″), ou décochez « Exporter vers NINA ».",
                )

            if manual_deg is not None:
                ra_deg, dec_deg = manual_deg
                manual_target_name = "ICRS saisi"
                nina_json = self._build_nina_deep_sky_json(manual_target_name, ra_deg, dec_deg)
                # Nom de fichier distinct des cibles catalogue éventuellement nommées pareil
                json_file = dest_path / "ICRS_saisi_manuel.json"
                with open(json_file, "w", encoding="utf-8") as f:
                    json.dump(nina_json, f, indent=2, ensure_ascii=False)
                exported_count += 1
                logger.info("Fichier JSON NINA (ICRS saisi) exporté : %s", json_file)

            messagebox.showinfo(
                "Export réussi",
                f"{exported_count} fichier(s) JSON exporté(s) avec succès.\n\n"
                f"Dossier de destination :\n{dest_dir}\n\n"
                f"Les fichiers sont prêts à être utilisés dans NINA.",
            )

        except Exception as e:
            logger.error(f"Erreur lors de l'export JSON : {e}", exc_info=True)
            messagebox.showerror("Erreur d'export", f"Erreur lors de l'export JSON : {e}")
    
    def display_results(self, sunset_naut, sunset_astro, sunrise_astro, sunrise_naut):
        """Affiche les résultats dans le tableau."""
        # Effacer le tableau
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        
        # Afficher les objets filtrés avec case à cocher, triés par ascension droite (RA)
        # Ne pas afficher les objets dont l'altitude maximale est < 20°
        for obj in sorted(self.filtered_objects, key=lambda x: x.ra):
            if getattr(obj, "alt_max", 0.0) < 20.0:
                continue
            mag_str = f"{obj.magnitude:.2f}" if obj.magnitude else "N/A"
            source = (getattr(obj, "exo_source", None) or "").strip().lower()
            display_name = self._display_name_for_object(obj)
            
            # Case à cocher (non sélectionné par défaut)
            check = "☐"
            if obj.name in self.selected_objects and self.selected_objects[obj.name]:
                check = "☑"
            
            tags = ()
            if source == EXOPLANET_PROVIDER_NASA:
                tags = ("source_nasa",)
            elif source == EXOPLANET_PROVIDER_EXOCLOCK:
                tags = ("source_exoclock",)

            self.results_tree.insert(
                "",
                tk.END,
                values=(
                    check,  # Colonne Sélection
                    display_name,
                    obj.obj_type,
                    mag_str,
                    obj.coordinates_sexagesimal(),
                    f"{obj.alt_max:.1f}°",
                ),
                tags=tags,
            )
    
    def plot_ephemerides(self, obs_date, start_time, end_time):
        """Trace le graphique d'altitude sur la nuit astronomique (−18°), centré sur 00h."""
        self.ax.clear()
        manual_deg = self._manual_coord_deg_for_plot()

        if start_time is None or end_time is None:
            self._setup_night_window_for_plot(obs_date)
            start_time = self.start_time
            end_time = self.end_time

        if start_time is None or end_time is None:
            self.ax.text(
                0.5,
                0.5,
                "Impossible de définir la nuit (date / fuseau).",
                ha="center",
                va="center",
                transform=self.ax.transAxes,
            )
            self.canvas.draw()
            return

        local_tz = getattr(self, "local_tz", None)
        sunset_local = getattr(self, "sunset_local", None)
        sunrise_local = getattr(self, "sunrise_local", None)
        if local_tz is None or sunset_local is None or sunrise_local is None:
            self.ax.text(
                0.5,
                0.5,
                "Paramètres de nuit manquants — lancez un calcul d’éphémérides ou « Mettre à jour le graphique ».",
                ha="center",
                va="center",
                transform=self.ax.transAxes,
                fontsize=10,
            )
            self.canvas.draw()
            return

        astro_dusk_local = getattr(self, "astro_dusk_local", None)
        astro_dawn_local = getattr(self, "astro_dawn_local", None)
        if astro_dusk_local is None or astro_dawn_local is None:
            astro_dusk_local = sunset_local
            astro_dawn_local = sunrise_local

        if start_time.tzinfo is None:
            start_time_tz = pytz.UTC.localize(start_time)
        else:
            start_time_tz = start_time.astimezone(pytz.UTC)

        if end_time.tzinfo is None:
            end_time_tz = pytz.UTC.localize(end_time)
        else:
            end_time_tz = end_time.astimezone(pytz.UTC)

        times_utc = []
        current = start_time_tz
        while current <= end_time_tz:
            times_utc.append(current)
            current = current + timedelta(minutes=15)

        times_local = [t_utc.astimezone(local_tz) for t_utc in times_utc]

        obs_date_local = local_tz.localize(obs_date.replace(hour=0, minute=0, second=0, microsecond=0))
        midnight_local = obs_date_local

        if sunset_local < midnight_local:
            midnight_local = local_tz.localize(
                (obs_date - timedelta(days=1)).replace(hour=0, minute=0, second=0, microsecond=0)
            )

        def _hours_since_midnight(t_loc: datetime) -> float:
            dh = (t_loc - midnight_local).total_seconds() / 3600.0
            if dh > 12:
                dh -= 24
            return dh

        times_hours = []
        for t_local in times_local:
            delta = t_local - midnight_local
            hours_from_midnight = delta.total_seconds() / 3600.0
            if hours_from_midnight > 12:
                hours_from_midnight -= 24
            times_hours.append(hours_from_midnight)

        moon_alts = [moon_altitude_deg_at_utc(t, self.location) for t in times_utc]

        dusk_astro_h = _hours_since_midnight(astro_dusk_local)
        dawn_astro_h = _hours_since_midnight(astro_dawn_local)
        dusk_geo_h = _hours_since_midnight(sunset_local)
        dawn_geo_h = _hours_since_midnight(sunrise_local)

        self.ax.axvline(
            x=dusk_astro_h,
            color="darkred",
            linestyle="--",
            linewidth=2,
            alpha=0.9,
            label="Fin crépuscule astro (−18°)",
            zorder=2,
        )
        self.ax.axvline(
            x=dawn_astro_h,
            color="darkorange",
            linestyle="--",
            linewidth=2,
            alpha=0.9,
            label="Début aube astro (−18°)",
            zorder=2,
        )
        self.ax.axvline(
            x=dusk_geo_h,
            color="gray",
            linestyle=":",
            linewidth=3.0,
            alpha=0.85,
            label="Coucher 0°",
            zorder=2,
        )
        self.ax.axvline(
            x=dawn_geo_h,
            color="gray",
            linestyle=":",
            linewidth=3.0,
            alpha=0.85,
            label="Lever 0°",
            zorder=2,
        )

        if dusk_astro_h < dawn_astro_h:
            self.ax.axvspan(dusk_astro_h, dawn_astro_h, alpha=0.15, color="darkblue", zorder=0)
        else:
            self.ax.axvspan(dusk_astro_h, 0, alpha=0.15, color="darkblue", zorder=0)
            self.ax.axvspan(0, dawn_astro_h, alpha=0.15, color="darkblue", zorder=0)

        self.ax.plot(
            times_hours,
            moon_alts,
            color="goldenrod",
            linestyle="-",
            linewidth=2.6,
            alpha=0.95,
            label="Lune (trajectoire)",
            zorder=6,
        )

        mid_utc = start_time_tz + (end_time_tz - start_time_tz) / 2
        ph_name, ph_illum, ph_elong = moon_phase_french_name_and_illumination(mid_utc)
        phase_box = (
            f"Lune — {ph_name}\n"
            f"Illuminée ≈ {ph_illum * 100.0:.0f} %  |  élong. ≈ {ph_elong:.0f}°"
        )
        self.ax.text(
            0.02,
            0.98,
            phase_box,
            transform=self.ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            linespacing=1.2,
            bbox=dict(
                boxstyle="round,pad=0.4",
                facecolor="wheat",
                alpha=0.9,
                edgecolor="saddlebrown",
            ),
            zorder=12,
        )

        selected_obj_list = [
            obj
            for obj in self.filtered_objects
            if obj.name in self.selected_objects and self.selected_objects[obj.name]
        ]
        if not selected_obj_list and manual_deg is None:
            self.ax.text(
                0.5,
                0.06,
                "Aucune cible cochée — Lune et repères de nuit affichés.",
                ha="center",
                va="bottom",
                transform=self.ax.transAxes,
                fontsize=9,
                color="dimgray",
                style="italic",
                zorder=11,
            )

        altitude_labels_done = set()

        def _exo_source_style(obj: "EphemerisObject"):
            source = (getattr(obj, "exo_source", None) or "").strip().lower()
            if source == EXOPLANET_PROVIDER_NASA:
                return source, "NASA", "blue", "navy"
            if source == EXOPLANET_PROVIDER_EXOCLOCK:
                return source, "ExoClock", "red", "darkred"
            return source, "Autre", "purple", "indigo"

        for obj in selected_obj_list:
            coord = SkyCoord(ra=obj.ra * u.deg, dec=obj.dec * u.deg, frame='icrs')
            altitudes = []
            
            # Calculer l'altitude pour chaque instant de la nuit
            for t_utc in times_utc:
                # Time() d'astropy gère automatiquement les datetime avec timezone
                t_astro = Time(t_utc)
                frame = AltAz(obstime=t_astro, location=self.location)
                altaz = coord.transform_to(frame)
                altitudes.append(altaz.alt.deg)
            
            # Ne tracer que si l'objet atteint au moins 25° d'altitude
            if max(altitudes) >= 25.0:
                _, source_label, line_color, _ = _exo_source_style(obj)
                label = (
                    f"Altitude {source_label}"
                    if source_label not in altitude_labels_done
                    else None
                )
                self.ax.plot(
                    times_hours,
                    altitudes,
                    alpha=0.7,
                    linewidth=1.5,
                    color=line_color,
                    label=label,
                )
                altitude_labels_done.add(source_label)

        if manual_deg is not None:
            ra_m, dec_m = manual_deg
            coord_m = SkyCoord(ra=ra_m * u.deg, dec=dec_m * u.deg, frame="icrs")
            alt_manual = []
            for t_utc in times_utc:
                t_astro = Time(t_utc)
                frame = AltAz(obstime=t_astro, location=self.location)
                alt_manual.append(coord_m.transform_to(frame).alt.deg)
            self.ax.plot(
                times_hours,
                alt_manual,
                color="darkviolet",
                linewidth=2.2,
                alpha=0.92,
                linestyle="-",
                label="Coordonnées ICRS (saisie)",
                zorder=5,
            )

        # Lignes de référence d'altitude
        self.ax.axhline(y=30, color='green', linestyle=':', alpha=0.5, linewidth=1)
        if any(obj.obj_type in ['asteroid', 'comet'] for obj in self.filtered_objects):
            self.ax.axhline(y=20, color='orange', linestyle=':', alpha=0.5, linewidth=1)
        self.ax.axhline(y=80, color='red', linestyle=':', alpha=0.5, linewidth=1)
        
        # Configuration des axes
        self.ax.set_xlabel("Heure locale", fontsize=13)
        self.ax.set_ylabel("Altitude (°)", fontsize=13)
        self.ax.set_title(
            f"Altitude (Lune + objets) — {obs_date.strftime('%Y-%m-%d')} (nuit astronomique −18°)",
            fontsize=15,
            fontweight="bold",
        )
        y_lo = min(15.0, float(np.min(moon_alts)) - 5.0)
        y_lo = max(-12.0, y_lo)
        self.ax.set_ylim(y_lo, 90)
        
        # Limites X : englober nuit astro et repères 0°
        x_min = min(dusk_astro_h, dawn_astro_h, dusk_geo_h, dawn_geo_h)
        x_max = max(dusk_astro_h, dawn_astro_h, dusk_geo_h, dawn_geo_h)
        
        # S'assurer que 0h est au centre visuel : ajuster les limites symétriquement autour de 0h
        max_abs = max(abs(x_min), abs(x_max))
        x_min = -max_abs - 0.3
        x_max = max_abs + 0.3
        self.ax.set_xlim(x_min, x_max)
        
        # Ajouter une ligne verticale pour marquer minuit (00h) au centre
        self.ax.axvline(x=0, color='grey', linestyle='-', linewidth=1, alpha=0.5)
        
        # Créer des ticks d'heures (toutes les heures pour la nuit)
        import matplotlib.ticker as ticker
        self.ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        
        # Formater les heures pour afficher correctement autour de minuit (centré sur 0h)
        def format_hour(x, pos):
            h = round(x)
            if h < 0:
                # Heures négatives (avant minuit) : -3h devient 21h, -6h devient 18h
                h = 24 + h
            elif h >= 24:
                # Heures >= 24 (après minuit) : 25h devient 01h
                h = h - 24
            elif h == 0:
                # Minuit
                return "00h"
            return f"{h:02d}h"
        
        self.ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_hour))
        
        self.ax.grid(True, alpha=0.3)

        # Schéma de transit (exoplanètes sélectionnées), 1 h avant ingress puis trapèze.
        # Couleurs demandées : NASA en bleu, ExoClock en rouge.
        def _utc_to_plot_hours(t_utc: datetime) -> float:
            if t_utc.tzinfo is None:
                t_utc = pytz.UTC.localize(t_utc)
            else:
                t_utc = t_utc.astimezone(pytz.UTC)
            tl = t_utc.astimezone(local_tz)
            dh = (tl - midnight_local).total_seconds() / 3600.0
            if dh > 12:
                dh -= 24
            return dh

        schema_labels_done = set()
        exo_idx = 0
        for obj in selected_obj_list:
            if obj.obj_type != "exoplanet":
                continue
            tmid = obj.exo_pl_tranmid_jd
            tdur = obj.exo_pl_trandur_h
            porb = obj.period
            if tmid is None or tdur is None or porb is None or porb <= 0 or tdur <= 0:
                continue
            k = exo_transit_k_for_night_utc_interval(
                tmid, porb, tdur, start_time_tz, end_time_tz
            )
            if k is None:
                continue
            half_dur_d = (tdur / 2.0) / 24.0
            t_mid_jd = tmid + k * porb
            t_ing_jd = t_mid_jd - half_dur_d
            t_egr_jd = t_mid_jd + half_dur_d
            margin_d = 1.0 / 24.0
            t_im1_jd = t_ing_jd - margin_d
            t_ep1_jd = t_egr_jd + margin_d

            dep = obj.exo_transitdepthcalc_pct
            if dep is None or dep <= 0:
                dep = 1.0
            dip = max(0.35, min(22.0, dep * 1.15))
            y_base = 17.0 + exo_idx * 2.4
            exo_idx += 1
            if y_base + dip > 40:
                y_base = max(16.0, 40.0 - dip)
            dip = min(dip, max(0.25, y_base - 15.6))

            dur_jd = t_egr_jd - t_ing_jd
            ramp_jd = min(dur_jd * 0.15, dur_jd / 2.0 - 1e-9)
            if ramp_jd < 1e-10:
                ramp_jd = max(dur_jd / 3.0, 1e-8)
            t_lo_jd = t_ing_jd + ramp_jd
            t_hi_jd = t_egr_jd - ramp_jd
            if t_lo_jd >= t_hi_jd:
                t_c = 0.5 * (t_ing_jd + t_egr_jd)
                eps = 1e-7
                t_lo_jd = t_c - eps
                t_hi_jd = t_c + eps

            jd_pts = [t_im1_jd, t_ing_jd, t_lo_jd, t_hi_jd, t_egr_jd, t_ep1_jd]
            y_pts = [
                y_base,
                y_base,
                y_base - dip,
                y_base - dip,
                y_base,
                y_base,
            ]
            t_utc_pts = [_jd_to_utc_datetime(j) for j in jd_pts]
            x_pts = [_utc_to_plot_hours(t) for t in t_utc_pts]

            source, source_label, line_color, text_color = _exo_source_style(obj)
            if source_label == "Autre" and source:
                source_label = source.upper()

            legend_key = source_label
            lbl = (
                f"Transit {source_label} (schéma, transitdepthcalc)"
                if legend_key not in schema_labels_done
                else None
            )
            self.ax.plot(
                x_pts,
                y_pts,
                color=line_color,
                linewidth=2.2,
                alpha=0.92,
                zorder=8,
                label=lbl,
            )
            x_mid_plot = _utc_to_plot_hours(_jd_to_utc_datetime(t_mid_jd))
            raw_tdc = getattr(obj, "exo_transitdepthcalc_pct", None)
            if raw_tdc is None:
                tdc_label = "transitdepthcalc = —"
            else:
                tdc_ppt = raw_tdc * 10.0
                tdc_label = f"transitdepthcalc = {tdc_ppt:.4f} ‰"
            self.ax.text(
                x_mid_plot,
                y_base + 1.0,
                tdc_label,
                ha="center",
                va="bottom",
                fontsize=8,
                color=text_color,
                zorder=9,
            )
            schema_labels_done.add(legend_key)

        h, lab = self.ax.get_legend_handles_labels()
        if h:
            leg = self.ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
            for txt in leg.get_texts():
                if "0°" in txt.get_text():
                    txt.set_fontweight("bold")

        self.canvas.draw()
