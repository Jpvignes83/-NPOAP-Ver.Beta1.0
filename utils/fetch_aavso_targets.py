#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script pour récupérer les exoplanètes et étoiles binaires à éclipses
depuis l'AAVSO Target Tool et les sauvegarder en CSV.

Deux méthodes :
1. Télécharger le catalogue index.csv directement (sans authentification)
2. Se connecter à Target Tool pour une requête personnalisée (avec authentification)
"""

import requests
import csv
import sys
from datetime import datetime, timedelta
from pathlib import Path
import logging
import getpass

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def download_index_csv(start_date: datetime, end_date: datetime, 
                       output_file: str = None) -> str:
    """
    Télécharge le catalogue index.csv depuis l'AAVSO pour une période donnée.
    
    Parameters
    ----------
    start_date : datetime
        Date de début
    end_date : datetime
        Date de fin
    output_file : str, optional
        Fichier de sortie (par défaut: index_YYYYMMDD_YYYYMMDD.csv)
    
    Returns
    -------
    str
        Chemin du fichier téléchargé
    """
    # Format des dates pour l'URL
    start_str = start_date.strftime("%Y%m%d")
    end_str = end_date.strftime("%Y%m%d")
    
    # URL du catalogue index.csv
    # L'AAVSO publie des catalogues par période
    base_url = "https://www.aavso.org/vsx/index.php"
    
    # Plusieurs URLs possibles pour le catalogue
    possible_urls = [
        f"https://www.aavso.org/vsx/index.php?view=download.list&type=all&period1={start_str}&period2={end_str}",
        f"https://www.aavso.org/downloads/index.csv",
        f"https://www.aavso.org/vsx/download.php?view=index.csv&period1={start_str}&period2={end_str}",
    ]
    
    # Essayer aussi une URL directe du catalogue complet ou par type
    # Pour exoplanètes et étoiles binaires à éclipses
    exo_url = "https://www.aavso.org/vsx/download.php?view=index.csv&type=EP"
    ev_url = "https://www.aavso.org/vsx/download.php?view=index.csv&type=EB"
    
    logger.info(f"Tentative de téléchargement du catalogue index.csv pour la période {start_str} - {end_str}")
    
    if output_file is None:
        output_file = f"index_{start_str}_{end_str}.csv"
    
    output_path = Path(output_file)
    
    # Essayer de télécharger depuis différentes sources
    content = None
    for url in possible_urls + [exo_url, ev_url]:
        try:
            logger.info(f"Tentative: {url}")
            response = requests.get(url, timeout=60, stream=True)
            if response.status_code == 200:
                # Vérifier si c'est du CSV
                content_type = response.headers.get('Content-Type', '').lower()
                if 'csv' in content_type or response.text.startswith('#') or ',' in response.text[:200]:
                    content = response.text
                    logger.info(f"✓ Catalogue téléchargé depuis: {url}")
                    break
        except Exception as e:
            logger.debug(f"Erreur avec {url}: {e}")
            continue
    
    if content is None:
        logger.warning("Impossible de télécharger le catalogue index.csv automatiquement")
        logger.info("Vous pouvez le télécharger manuellement depuis:")
        logger.info("  https://www.aavso.org/vsx/index.php?view=download.list")
        logger.info("Ou utiliser l'authentification avec Target Tool")
        return None
    
    # Sauvegarder
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    logger.info(f"✓ Catalogue sauvegardé: {output_path}")
    
    # Afficher quelques statistiques
    try:
        lines = content.strip().split('\n')
        logger.info(f"Nombre de lignes: {len(lines)}")
        if len(lines) > 1:
            num_targets = len([l for l in lines if l.strip() and not l.strip().startswith('#')])
            logger.info(f"Nombre de cibles: {num_targets}")
    except Exception as e:
        logger.debug(f"Erreur lors de l'analyse: {e}")
    
    return str(output_path)


def fetch_aavso_targets_with_auth(output_file: str = None, exoplanets: bool = True, 
                                  eclipsing_binaries: bool = True,
                                  username: str = None, password: str = None):
    """
    Récupère les cibles depuis l'AAVSO Target Tool avec authentification.
    
    Parameters
    ----------
    output_file : str, optional
        Chemin du fichier CSV de sortie
    exoplanets : bool
        Inclure les exoplanètes
    eclipsing_binaries : bool
        Inclure les étoiles binaires à éclipses
    username : str, optional
        Nom d'utilisateur AAVSO (demandé interactivement si non fourni)
    password : str, optional
        Mot de passe AAVSO (demandé interactivement si non fourni)
    """
    base_url = "https://targettool.aavso.org/TargetTool/default/index"
    
    # Demander les identifiants si nécessaire
    if username is None:
        username = input("Nom d'utilisateur AAVSO: ")
    if password is None:
        password = getpass.getpass("Mot de passe AAVSO: ")
    
    # Construire l'URL avec les paramètres
    params = []
    if exoplanets:
        params.append('exo=on')
    if eclipsing_binaries:
        params.append('ev=on')
    if params:
        params.append('settype=true')
    
    url = f"{base_url}?{'&'.join(params)}"
    
    logger.info(f"Téléchargement depuis: {url}")
    
    # Créer une session pour gérer l'authentification
    session = requests.Session()
    
    try:
        # Page de connexion
        login_url = "https://targettool.aavso.org/TargetTool/default/login"
        
        # Obtenir la page de connexion pour récupérer les tokens CSRF si nécessaire
        login_page = session.get(login_url, timeout=30)
        login_page.raise_for_status()
        
        # Tenter de se connecter
        # Note: Les champs de formulaire peuvent varier, à ajuster selon le site réel
        login_data = {
            'username': username,
            'password': password,
        }
        
        # Chercher un token CSRF dans le HTML si présent
        import re
        csrf_match = re.search(r'name=["\']_token["\']\s+value=["\']([^"\']+)["\']', login_page.text)
        if csrf_match:
            login_data['_token'] = csrf_match.group(1)
        
        login_response = session.post(login_url, data=login_data, timeout=30)
        
        # Vérifier si la connexion a réussi
        if 'login' in login_response.url.lower() or 'error' in login_response.text.lower():
            logger.error("Échec de l'authentification. Vérifiez vos identifiants.")
            return None
        
        logger.info("✓ Authentification réussie")
        
        # Maintenant accéder à la page des cibles
        response = session.get(url, timeout=30)
        response.raise_for_status()
        
        # Tenter d'obtenir le CSV
        csv_urls = [
            url.replace('/index', '/exportCSV'),
            url + '&exportCSV=1',
            url + '&format=csv',
            url.replace('/index', '/export')
        ]
        
        content = None
        for csv_url in csv_urls:
            try:
                logger.info(f"Tentative export CSV: {csv_url}")
                csv_response = session.get(csv_url, timeout=30)
                if csv_response.status_code == 200:
                    content_type = csv_response.headers.get('Content-Type', '').lower()
                    if 'csv' in content_type or 'text/plain' in content_type:
                        content = csv_response.text
                        logger.info("✓ CSV récupéré avec succès")
                        break
            except Exception as e:
                logger.debug(f"Erreur avec {csv_url}: {e}")
                continue
        
        # Si pas de CSV direct, parser le HTML
        if content is None:
            logger.info("CSV direct non disponible, tentative de parsing HTML...")
            html_content = response.text
            
            try:
                from bs4 import BeautifulSoup
                soup = BeautifulSoup(html_content, 'html.parser')
                
                # Chercher un tableau
                table = soup.find('table') or soup.find('table', {'id': 'targets'}) or \
                       soup.find('table', {'class': 'dataTable'}) or \
                       soup.find('table', {'class': 'table'})
                
                if table:
                    rows = table.find_all('tr')
                    if len(rows) > 0:
                        headers = [th.get_text(strip=True) for th in rows[0].find_all(['th', 'td'])]
                        
                        import io
                        csv_buffer = io.StringIO()
                        writer = csv.writer(csv_buffer)
                        writer.writerow(headers)
                        
                        for row in rows[1:]:
                            cells = [td.get_text(strip=True) for td in row.find_all(['td', 'th'])]
                            if cells:
                                writer.writerow(cells)
                        
                        content = csv_buffer.getvalue()
                        logger.info("✓ Tableau HTML parsé avec succès")
            except ImportError:
                logger.warning("BeautifulSoup non disponible pour le parsing HTML")
            except Exception as e:
                logger.warning(f"Erreur lors du parsing HTML: {e}")
        
        if content is None:
            logger.error("Impossible de récupérer les données")
            return None
        
        # Sauvegarder
        if output_file is None:
            today = datetime.now().strftime("%Y%m%d")
            output_file = f"aavso_targets_{today}.csv"
        
        output_path = Path(output_file)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        logger.info(f"✓ Données sauvegardées: {output_path}")
        
        # Statistiques
        try:
            lines = content.strip().split('\n')
            logger.info(f"Nombre de lignes: {len(lines)}")
            if len(lines) > 1:
                num_targets = len(lines) - 1
                logger.info(f"Nombre de cibles: {num_targets}")
        except Exception as e:
            logger.debug(f"Erreur lors de l'analyse: {e}")
        
        return str(output_path)
        
    except requests.RequestException as e:
        logger.error(f"Erreur lors de la requête: {e}")
        return None
    except Exception as e:
        logger.error(f"Erreur inattendue: {e}", exc_info=True)
        return None
    finally:
        session.close()


def fetch_aavso_targets(output_file: str = None, exoplanets: bool = True, 
                        eclipsing_binaries: bool = True):
    """
    Récupère les cibles depuis l'AAVSO Target Tool.
    URL directe : https://targettool.aavso.org/TargetTool/default/index.csv?ev=on&exo=on&settype=true
    
    Parameters
    ----------
    output_file : str, optional
        Chemin du fichier CSV de sortie. Par défaut: aavso_targets_YYYYMMDD.csv
    exoplanets : bool
        Inclure les exoplanètes (True par défaut)
    eclipsing_binaries : bool
        Inclure les étoiles binaires à éclipses (True par défaut)
    """
    base_url = "https://targettool.aavso.org/TargetTool/default/index.csv"
    
    # Construire l'URL avec les paramètres
    params = []
    if exoplanets:
        params.append('exo=on')
    if eclipsing_binaries:
        params.append('ev=on')
    if params:
        params.append('settype=true')
    
    url = f"{base_url}?{'&'.join(params)}"
    
    logger.info(f"Téléchargement depuis: {url}")
    
    try:
        # L'URL retourne directement le CSV, pas besoin d'authentification
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        
        content = response.text
        
        # Vérifier que c'est bien du CSV
        if not content or len(content.strip()) == 0:
            logger.error("Réponse vide du serveur")
            return None
        
        logger.info("✓ CSV récupéré avec succès")
        
        # Vérifier que le contenu commence par un en-tête CSV
        if not content.strip().startswith('"') and ',' not in content[:100]:
            logger.warning("Le contenu ne semble pas être du CSV valide")
            logger.debug(f"Premières lignes: {content[:500]}")
        
        # Sauvegarder directement
        if output_file is None:
            today = datetime.now().strftime("%Y%m%d")
            output_file = f"aavso_targets_{today}.csv"
        
        output_path = Path(output_file)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        logger.info(f"✓ Données sauvegardées dans: {output_path}")
        
        # Afficher quelques statistiques
        try:
            lines = content.strip().split('\n')
            logger.info(f"Nombre de lignes: {len(lines)}")
            if len(lines) > 1:
                # Compter les données (ignorer les lignes vides et l'en-tête)
                num_targets = len([l for l in lines[1:] if l.strip()])
                logger.info(f"Nombre de cibles: {num_targets}")
                
                # Afficher les premières lignes pour vérification
                if len(lines) > 0:
                    logger.debug("Première ligne (en-têtes):")
                    logger.debug(lines[0][:200])
                    if len(lines) > 1:
                        logger.debug("Première cible:")
                        logger.debug(lines[1][:200])
        except Exception as e:
            logger.debug(f"Erreur lors de l'analyse du CSV: {e}")
        
        return str(output_path)
        
    except requests.RequestException as e:
        logger.error(f"Erreur lors de la requête: {e}")
        # Fallback: essayer avec BeautifulSoup si HTML retourné
        try:
            logger.warning("Réponse non CSV, vérification du format...")
            logger.debug(f"Content-Type: {response.headers.get('Content-Type')}")
            logger.debug(f"Premiers caractères: {response.text[:200]}")
            return None
        except Exception as e:
            logger.debug(f"Erreur lors de la vérification: {e}")
            return None
    except Exception as e:
        logger.error(f"Erreur inattendue: {e}", exc_info=True)
        return None


def main():
    """Point d'entrée principal."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Récupère les cibles AAVSO (exoplanètes et étoiles binaires à éclipses)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples:
  # Télécharger le catalogue index.csv pour une période
  python utils/fetch_aavso_targets.py --index --start-date 2025-01-01 --end-date 2025-12-31
  
  # Se connecter à Target Tool (identifiants demandés)
  python utils/fetch_aavso_targets.py --auth
  
  # Télécharger index.csv pour les 30 prochains jours
  python utils/fetch_aavso_targets.py --index --days 30
  
  # Ou depuis le répertoire utils/:
  cd utils
  python fetch_aavso_targets.py --index --days 30
        """
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help="Fichier CSV de sortie"
    )
    
    # Options pour index.csv
    parser.add_argument(
        '--index',
        action='store_true',
        help="Télécharger le catalogue index.csv (sans authentification)"
    )
    parser.add_argument(
        '--start-date',
        type=str,
        help="Date de début pour index.csv (format: YYYY-MM-DD)"
    )
    parser.add_argument(
        '--end-date',
        type=str,
        help="Date de fin pour index.csv (format: YYYY-MM-DD)"
    )
    parser.add_argument(
        '--days',
        type=int,
        help="Nombre de jours à partir d'aujourd'hui pour index.csv"
    )
    
    # Options pour authentification
    parser.add_argument(
        '--auth',
        action='store_true',
        help="Utiliser l'authentification pour Target Tool"
    )
    parser.add_argument(
        '--username',
        type=str,
        help="Nom d'utilisateur AAVSO (pour --auth)"
    )
    parser.add_argument(
        '--password',
        type=str,
        help="Mot de passe AAVSO (pour --auth, non recommandé en ligne de commande)"
    )
    
    # Filtres
    parser.add_argument(
        '--no-exoplanets',
        action='store_true',
        help="Exclure les exoplanètes"
    )
    parser.add_argument(
        '--no-eclipsing-binaries',
        action='store_true',
        help="Exclure les étoiles binaires à éclipses"
    )
    
    args = parser.parse_args()
    
    # Méthode index.csv
    if args.index:
        # Déterminer les dates
        if args.days:
            start_date = datetime.now()
            end_date = start_date + timedelta(days=args.days)
        elif args.start_date and args.end_date:
            try:
                start_date = datetime.strptime(args.start_date, "%Y-%m-%d")
                end_date = datetime.strptime(args.end_date, "%Y-%m-%d")
            except ValueError:
                logger.error("Format de date invalide. Utilisez YYYY-MM-DD")
                sys.exit(1)
        else:
            # Par défaut: aujourd'hui à aujourd'hui + 30 jours
            start_date = datetime.now()
            end_date = start_date + timedelta(days=30)
            logger.info(f"Utilisation de la période par défaut: {start_date.date()} - {end_date.date()}")
        
        output = download_index_csv(start_date, end_date, args.output)
    
    # Méthode avec authentification
    elif args.auth:
        output = fetch_aavso_targets_with_auth(
            output_file=args.output,
            exoplanets=not args.no_exoplanets,
            eclipsing_binaries=not args.no_eclipsing_binaries,
            username=args.username,
            password=args.password
        )
    
    # Par défaut: utiliser l'URL directe (sans authentification nécessaire)
    else:
        logger.info("Téléchargement direct depuis AAVSO Target Tool (pas d'authentification nécessaire)")
        output = fetch_aavso_targets(
            output_file=args.output,
            exoplanets=not args.no_exoplanets,
            eclipsing_binaries=not args.no_eclipsing_binaries
        )
    
    if output:
        print(f"\n✓ Fichier créé: {output}")
        sys.exit(0)
    else:
        print("\n✗ Échec de la récupération")
        sys.exit(1)


if __name__ == "__main__":
    main()

