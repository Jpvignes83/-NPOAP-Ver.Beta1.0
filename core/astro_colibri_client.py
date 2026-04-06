"""
Client pour l'API Astro-COLIBRI (transients, cone search, événements).
Remplace l'usage de TNS pour la recherche d'objets transitoires.

Documentation: https://astro-colibri.science/apidoc
Site: https://astro-colibri.com/
"""

import json
import logging
from typing import Optional, Dict, List, Any
from datetime import datetime
import urllib.request
import urllib.parse
import urllib.error

logger = logging.getLogger(__name__)

BASE_URL = "https://astro-colibri.science"


class AstroColibriClient:
    """
    Client pour l'API Astro-COLIBRI.
    Cone search, latest transients, event et source_details.
    Inscription gratuite pour obtenir un uid (100 requêtes/jour pour cone_search et latest_transients).
    """

    def __init__(self, uid: Optional[str] = None):
        """
        Parameters
        ----------
        uid : str, optional
            User ID Astro-COLIBRI (compte gratuit). Recommandé pour cone_search / latest_transients (quota 100 req/jour).
        """
        self.uid = (uid or "").strip() or None

    def _get(self, path: str, params: Optional[Dict[str, str]] = None) -> Optional[Dict[str, Any]]:
        """Requête GET."""
        url = f"{BASE_URL}{path}"
        if params:
            qs = urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})
            url = f"{url}?{qs}"
        req = urllib.request.Request(url, method="GET")
        req.add_header("Accept", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            logger.error("Astro-COLIBRI HTTP %s: %s", e.code, e.reason)
            try:
                logger.error("%s", e.read().decode("utf-8"))
            except Exception:
                pass
            return None
        except Exception as e:
            logger.error("Astro-COLIBRI request error: %s", e)
            return None

    def _post(self, path: str, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Requête POST avec body JSON."""
        url = f"{BASE_URL}{path}"
        body = json.dumps(data).encode("utf-8")
        req = urllib.request.Request(url, data=body, method="POST")
        req.add_header("Content-Type", "application/json")
        req.add_header("Accept", "application/json")
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            logger.error("Astro-COLIBRI HTTP %s: %s", e.code, e.reason)
            try:
                logger.error("%s", e.read().decode("utf-8"))
            except Exception:
                pass
            return None
        except Exception as e:
            logger.error("Astro-COLIBRI request error: %s", e)
            return None

    def cone_search(
        self,
        ra: float,
        dec: float,
        radius_deg: float,
        time_min: Optional[str] = None,
        time_max: Optional[str] = None,
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Recherche en cône (transients + catalogues).
        Limité à 100 req/jour si uid fourni.

        Parameters
        ----------
        ra, dec : float
            Centre en degrés (ICRS).
        radius_deg : float
            Rayon en degrés.
        time_min, time_max : str, optional
            Plage temporelle en ISO (ex. "2024-01-01T00:00:00").

        Returns
        -------
        Liste d'événements (voevents). Chaque entrée a trigger_id, source_name, ra, dec, time, type, etc.
        """
        body = {
            "properties": {"position": {"ra": ra, "dec": dec}, "radius": radius_deg},
            "time_range": {},
            "return_format": "json",
        }
        if self.uid:
            body["uid"] = self.uid
        if time_min:
            body["time_range"]["min"] = time_min if "T" in time_min else f"{time_min}T00:00:00"
        if time_max:
            body["time_range"]["max"] = time_max if "T" in time_max else f"{time_max}T23:59:59"
        if not body["time_range"]:
            body["time_range"] = {"min": "2020-01-01T00:00:00", "max": "2030-12-31T23:59:59"}

        out = self._post("/cone_search", body)
        if not out:
            return None
        voevents = out.get("voevents") or []
        return voevents

    def latest_transients(self, time_min: str, time_max: str) -> Optional[List[Dict[str, Any]]]:
        """
        Derniers transients dans une fenêtre temporelle.
        Nécessite uid (100 req/jour).

        Parameters
        ----------
        time_min, time_max : str
            ISO format (ex. "2024-01-01T00:00:00").
        """
        if not self.uid:
            logger.warning("latest_transients requiert un uid Astro-COLIBRI")
            return None
        body = {
            "time_range": {"min": time_min, "max": time_max},
            "uid": self.uid,
            "return_format": "json",
        }
        out = self._post("/latest_transients", body)
        if not out:
            return None
        return out.get("voevents") or []

    def get_event(
        self,
        trigger_id: Optional[str] = None,
        source_name: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Récupère un événement par trigger_id ou source_name (GET /event).
        """
        if trigger_id:
            out = self._get("/event", {"trigger_id": trigger_id})
        elif source_name:
            out = self._get("/event", {"source_name": source_name})
        else:
            return None
        if not out:
            return None
        if isinstance(out, list):
            return out[0] if out else None
        return out

    def get_source_details(self, name: str) -> Optional[Dict[str, Any]]:
        """
        Coordonnées et infos de base par nom (GET /source_details).
        """
        return self._get("/source_details", {"name": name})
