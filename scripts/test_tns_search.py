#!/usr/bin/env python3
"""
Test rapide de l'outil de recherche TNS.
À lancer depuis la racine du projet : python scripts/test_tns_search.py
"""
import sys
from pathlib import Path

# Ajouter la racine du projet au path
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))

import config
from core.tns_client import TNSClient, TNSUnauthorizedError

def main():
    import os
    tns_config = getattr(config, "TNS_CONFIG", {})
    # Permettre override par variables d'environnement (pour test sans modifier config)
    api_key = (os.environ.get("TNS_API_KEY") or tns_config.get("api_key") or "").strip()
    bot_id = (os.environ.get("TNS_BOT_ID") or tns_config.get("bot_id") or "").strip()
    tns_id = (tns_config.get("tns_id") or "").strip()
    tns_name = (tns_config.get("tns_name") or "").strip()
    if os.environ.get("TNS_BOT_ID"):
        marker_type = "bot"
    else:
        marker_type = (tns_config.get("tns_marker_type") or "user").strip().lower()
    use_sandbox = tns_config.get("use_sandbox", False)

    if not api_key:
        print("Échec : API Key TNS non configurée.")
        print("Ouvrez NPOAP, onglet Photometrie transitoires, Configuration API :")
        print("  remplissez Type compte, Bot ID (ou TNS ID/Nom), API Key, puis Sauvegarder Config.")
        return 1

    if marker_type == "bot":
        if not bot_id:
            print("Échec : en mode 'bot', le Bot ID doit être renseigné.")
            return 1
    else:
        tns_id = (tns_config.get("tns_id") or "").strip()
        tns_name = (tns_config.get("tns_name") or "").strip()
        if not tns_id or not tns_name:
            print("Échec : en mode 'user', TNS ID et Nom doivent être renseignés.")
            return 1

    client = TNSClient(
        bot_id=bot_id or None,
        api_key=api_key,
        bot_name=tns_config.get("bot_name", "NPOAP"),
        use_sandbox=use_sandbox,
        tns_marker_type=marker_type,
        tns_id=tns_id or None,
        tns_name=tns_name or None,
    )

    print("Test recherche TNS (objet connu: 2021rf)...")
    print("User-Agent:", client.user_agent[:60] + "..." if len(client.user_agent) > 60 else client.user_agent)

    try:
        results = client.search_objects(objname="2021rf")
        if results is None:
            print("Réponse vide ou erreur API (vérifiez les logs).")
            return 1
        if not results:
            print("Aucun résultat (liste vide).")
            return 0
        print(f"OK — {len(results)} objet(s) trouvé(s).")
        for i, obj in enumerate(results[:3]):
            name = obj.get("objname") or obj.get("name", "?")
            print(f"  [{i+1}] {name}")
        if len(results) > 3:
            print(f"  ... et {len(results) - 3} autre(s).")
        return 0
    except TNSUnauthorizedError as e:
        print("Erreur 401 Unauthorized:", e)
        print("Vérifiez Bot ID (ou TNS ID/Nom) et API Key, et que vous utilisez Production (use_sandbox=False).")
        return 1
    except Exception as e:
        print("Erreur:", e)
        return 1

if __name__ == "__main__":
    sys.exit(main())
