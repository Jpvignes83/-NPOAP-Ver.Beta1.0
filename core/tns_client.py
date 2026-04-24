"""
Client pour l'API Transient Name Server (TNS).
Permet de rechercher et récupérer des informations sur les objets transitoires.

Documentation: https://www.wis-tns.org/api
"""

import json
import logging
from typing import Optional, Dict, List, Any
from pathlib import Path
import urllib.request
import urllib.parse
import urllib.error

logger = logging.getLogger(__name__)


class TNSUnauthorizedError(Exception):
    """Levée lorsque l'API TNS retourne 401 Unauthorized (identifiants invalides ou mauvais environnement)."""
    pass


# URLs de l'API TNS
TNS_SANDBOX_URL = "https://sandbox.wis-tns.org/api"
TNS_PRODUCTION_URL = "https://www.wis-tns.org/api"


class TNSClient:
    """
    Client pour interagir avec l'API TNS (Transient Name Server).
    Supporte le user-agent en mode "user" (compte utilisateur) ou "bot" (compte bot).
    """

    def __init__(self, bot_id: Optional[str] = None, api_key: Optional[str] = None,
                 bot_name: Optional[str] = None, use_sandbox: bool = True,
                 tns_marker_type: str = "user",
                 tns_id: Optional[str] = None, tns_name: Optional[str] = None):
        """
        Initialise le client TNS.

        Parameters
        ----------
        bot_id : str, optional
            ID du bot TNS (pour type "bot")
        api_key : str, optional
            Clé API TNS (obligatoire pour les requêtes)
        bot_name : str, optional
            Nom du bot TNS (pour type "bot")
        use_sandbox : bool
            Si True, utilise l'environnement SANDBOX (recommandé pour les tests)
        tns_marker_type : str
            "user" ou "bot". Si "user", le user-agent utilise tns_id et tns_name.
        tns_id : str, optional
            TNS ID utilisateur (ex: "2661") pour type "user"
        tns_name : str, optional
            Nom utilisateur TNS (ex: "jpvignes") pour type "user"
        """
        self.bot_id = bot_id
        self.api_key = api_key
        self.bot_name = bot_name or "NPOAP"
        self.use_sandbox = use_sandbox
        self.base_url = TNS_SANDBOX_URL if use_sandbox else TNS_PRODUCTION_URL
        self.tns_marker_type = (tns_marker_type or "user").strip().lower()
        if self.tns_marker_type not in ("user", "bot"):
            self.tns_marker_type = "user"

        # Construire le user-agent selon la spécification TNS
        # Ex: tns_marker{"tns_id":2661,"type": "user", "name":"jpvignes"}
        # https://www.wis-tns.org/api
        if self.tns_marker_type == "user" and (tns_id or tns_name):
            tid = (tns_id or "").strip()
            name = (tns_name or "").strip() or "user"
            try:
                tid_int = int(tid) if tid.isdigit() else None
            except (ValueError, AttributeError):
                tid_int = None
            if tid_int is not None:
                self.user_agent = f'tns_marker{{"tns_id":{tid_int},"type": "user", "name":"{name}"}}'
            else:
                self.user_agent = f'tns_marker{{"tns_id":"{tid}","type": "user", "name":"{name}"}}'
        elif self.tns_marker_type == "bot" and bot_id:
            try:
                tid = int(bot_id) if str(bot_id).strip().isdigit() else bot_id
                self.user_agent = f'tns_marker{{"tns_id":{tid},"type": "bot", "name":"{self.bot_name}"}}'
            except (ValueError, AttributeError):
                self.user_agent = f'tns_marker{{"tns_id":"{bot_id}","type": "bot", "name":"{self.bot_name}"}}'
        else:
            # Config incomplète (ID/nom vides)
            if self.tns_marker_type == "user":
                self.user_agent = 'tns_marker{"tns_id":"","type": "user", "name":""}'
                logger.debug("TNS : TNS ID / nom non renseignés (requêtes refusées tant que non configuré).")
            else:
                self.user_agent = f'tns_marker{{"tns_id":"","type": "bot", "name":"{self.bot_name}"}}'
                logger.debug("TNS : bot_id non défini.")
    
    def _make_request(self, endpoint: str, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Effectue une requête POST vers l'API TNS.
        
        Parameters
        ----------
        endpoint : str
            Endpoint de l'API (ex: "get/search", "get/object")
        data : dict
            Données JSON à envoyer dans la requête
        
        Returns
        -------
        dict ou None
            Réponse JSON de l'API, ou None en cas d'erreur
        """
        # Pour Search/Get, le manuel TNS exige api_key + tns_marker. On envoie la requête si on a au moins le marker.
        if not self.api_key and self.tns_marker_type == "bot":
            logger.error("Clé API TNS obligatoire en mode bot")
            return None

        url = f"{self.base_url}/{endpoint}"

        params = {'data': json.dumps(data)}
        if self.api_key:
            params['api_key'] = self.api_key
        post_data = urllib.parse.urlencode(params).encode('utf-8')
        
        # Créer la requête
        req = urllib.request.Request(url, data=post_data, method='POST')
        req.add_header('User-Agent', self.user_agent)
        req.add_header('Content-Type', 'application/x-www-form-urlencoded')
        
        # Log pour déboguer
        logger.debug(f"Requête TNS: URL={url}, User-Agent={self.user_agent}")
        
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                result = json.loads(response.read().decode('utf-8'))
                
                # Vérifier les erreurs dans la réponse
                if result.get('id_code') == 200:
                    return result.get('data', {})
                else:
                    error_msg = result.get('id_message', 'Erreur inconnue')
                    logger.error(f"Erreur API TNS: {error_msg}")
                    return None
                    
        except urllib.error.HTTPError as e:
            logger.error(f"Erreur HTTP TNS: {e.code} - {e.reason}")
            try:
                error_body = e.read().decode('utf-8')
                logger.error(f"Corps de l'erreur: {error_body}")
            except Exception:
                pass
            if e.code == 401:
                logger.info(f"TNS User-Agent envoyé: {self.user_agent}")
                raise TNSUnauthorizedError(
                    "Authentification refusée (401).\n\n"
                    "Pour les APIs Search/Get, TNS exige api_key + tns_marker (manuel TNS).\n"
                    "→ Utilisez le mode Bot : Bot ID (ex. 197845), Nom (bot) et API Key.\n"
                    "→ Ou en mode User : TNS ID, Nom et une API Key de votre compte (My Account sur wis-tns.org).\n"
                    "Environnement : Production (pas Sandbox)."
                )
            return None
        except urllib.error.URLError as e:
            logger.error(f"Erreur URL TNS: {e.reason}")
            return None
        except json.JSONDecodeError as e:
            logger.error(f"Erreur décodage JSON TNS: {e}")
            return None
        except Exception as e:
            logger.error(f"Erreur requête TNS: {e}", exc_info=True)
            return None
    
    def search_objects(self, objname: Optional[str] = None, 
                      ra: Optional[float] = None, dec: Optional[float] = None,
                      radius: Optional[float] = None, radius_units: str = "arcsec",
                      discovery_date_start: Optional[str] = None,
                      discovery_date_end: Optional[str] = None,
                      discovery_mag_min: Optional[float] = None,
                      discovery_mag_max: Optional[float] = None,
                      internal_name: Optional[str] = None,
                      public_timestamp: Optional[str] = None,
                      objtype: Optional[int] = None,
                      redshift: Optional[float] = None,
                      hostname: Optional[str] = None,
                      groupid: Optional[int] = None,
                      classified_sne: Optional[int] = None,
                      include_frb: Optional[int] = None) -> Optional[List[Dict[str, Any]]]:
        """
        Recherche des objets transitoires selon des critères.
        
        Parameters
        ----------
        objname : str, optional
            Nom de l'objet (ex: "2021rf", "SN 2021abc")
        ra : float, optional
            Ascension droite (degrés)
        dec : float, optional
            Déclinaison (degrés)
        radius : float, optional
            Rayon de recherche (nécessite ra et dec)
        radius_units : str
            Unités du rayon ("arcsec", "arcmin", "deg")
        discovery_date_start : str, optional
            Date de découverte début (format: "YYYY-MM-DD")
        discovery_date_end : str, optional
            Date de découverte fin (format: "YYYY-MM-DD")
        discovery_mag_min : float, optional
            Magnitude de découverte minimale
        discovery_mag_max : float, optional
            Magnitude de découverte maximale
        internal_name : str, optional
            Nom interne de l'objet
        public_timestamp : str, optional
            Timestamp public
        objtype : int, optional
            Type d'objet (voir documentation TNS)
        redshift : float, optional
            Redshift
        hostname : str, optional
            Nom de l'hôte
        groupid : int, optional
            ID du groupe
        classified_sne : int, optional
            SNe classifiées (0 ou 1)
        include_frb : int, optional
            Inclure les FRB (0 ou 1)
        
        Returns
        -------
        list ou None
            Liste des objets trouvés, ou None en cas d'erreur
        """
        search_data = {}
        
        if objname:
            search_data['objname'] = objname
        if ra is not None and dec is not None:
            search_data['ra'] = ra
            search_data['dec'] = dec
            if radius is not None:
                search_data['radius'] = radius
                search_data['radius_units'] = radius_units
        if discovery_date_start:
            search_data['discovery_date_start'] = discovery_date_start
        if discovery_date_end:
            search_data['discovery_date_end'] = discovery_date_end
        if discovery_mag_min is not None:
            search_data['discovery_mag_min'] = discovery_mag_min
        if discovery_mag_max is not None:
            search_data['discovery_mag_max'] = discovery_mag_max
        if internal_name:
            search_data['internal_name'] = internal_name
        if public_timestamp:
            search_data['public_timestamp'] = public_timestamp
        if objtype is not None:
            search_data['objtype'] = objtype
        if redshift is not None:
            search_data['redshift'] = redshift
        if hostname:
            search_data['hostname'] = hostname
        if groupid is not None:
            search_data['groupid'] = groupid
        if classified_sne is not None:
            search_data['classified_sne'] = classified_sne
        if include_frb is not None:
            search_data['include_frb'] = include_frb
        
        if not search_data:
            logger.warning("Aucun critère de recherche fourni")
            return None
        
        result = self._make_request("get/search", search_data)
        
        if result and 'reply' in result:
            return result['reply']
        return None
    
    def get_object(self, objname: str, photometry: bool = False, 
                   spectra: bool = False) -> Optional[Dict[str, Any]]:
        """
        Récupère les détails d'un objet transitoire.
        
        Parameters
        ----------
        objname : str
            Nom de l'objet (ex: "2021rf", "SN 2021abc")
        photometry : bool
            Inclure les données photométriques
        spectra : bool
            Inclure les données spectrales
        
        Returns
        -------
        dict ou None
            Détails de l'objet, ou None en cas d'erreur
        """
        data = {'objname': objname}
        
        if photometry:
            data['photometry'] = '1'
        if spectra:
            data['spectra'] = '1'
        
        result = self._make_request("get/object", data)
        
        if result and 'reply' in result:
            return result['reply']
        return None
    
    def get_file(self, objname: str, file_type: str = "spectrum", 
                 file_id: Optional[int] = None) -> Optional[bytes]:
        """
        Récupère un fichier associé à un objet (spectre, etc.).
        
        Parameters
        ----------
        objname : str
            Nom de l'objet
        file_type : str
            Type de fichier ("spectrum", etc.)
        file_id : int, optional
            ID du fichier
        
        Returns
        -------
        bytes ou None
            Contenu du fichier, ou None en cas d'erreur
        """
        data = {
            'objname': objname,
            'file_type': file_type
        }
        if file_id is not None:
            data['file_id'] = file_id
        
        # Cette méthode nécessite une implémentation spécifique pour les fichiers
        # Pour l'instant, retourner None
        logger.warning("get_file n'est pas encore implémenté")
        return None
