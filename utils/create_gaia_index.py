#!/usr/bin/env python3
"""
Script utilitaire pour créer un index des fichiers Gaia déjà téléchargés.
Permet de reprendre un téléchargement en excluant les fichiers existants.
"""

import sys
from pathlib import Path

# Ajouter le répertoire parent au path pour importer les modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.gaia_dr3_extractor import GaiaDR3Extractor
import argparse
import json


def main():
    parser = argparse.ArgumentParser(
        description="Crée un index des fichiers Gaia déjà téléchargés"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Répertoire contenant les fichiers Gaia (défaut: demandé interactivement)"
    )
    parser.add_argument(
        "--hemisphere",
        type=str,
        choices=["north", "south"],
        default=None,
        help="Filtrer par hémisphère (north ou south)"
    )
    parser.add_argument(
        "--mag-limit",
        type=float,
        default=None,
        help="Filtrer par magnitude limite"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="gaia_index.json",
        help="Nom du fichier d'index (défaut: gaia_index.json)"
    )
    
    args = parser.parse_args()
    
    # Déterminer le répertoire
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        # Demander interactivement
        output_dir_str = input("Répertoire contenant les fichiers Gaia [C:\\Users\\docke\\OneDrive\\Documents\\catalogues]: ")
        if not output_dir_str.strip():
            output_dir = Path("C:/Users/docke/OneDrive/Documents/catalogues")
        else:
            output_dir = Path(output_dir_str)
    
    if not output_dir.exists():
        print(f"❌ Erreur : Le répertoire {output_dir} n'existe pas")
        return 1
    
    print(f"📁 Répertoire : {output_dir}")
    print(f"🔍 Création de l'index...")
    
    # Créer l'extracteur et l'index
    extractor = GaiaDR3Extractor(output_dir=output_dir)
    index = extractor.create_index(
        hemisphere=args.hemisphere,
        mag_limit=args.mag_limit
    )
    
    # Sauvegarder l'index
    index_path = extractor.save_index(index, filename=args.output)
    
    # Afficher un résumé
    print("\n" + "=" * 80)
    print("RÉSUMÉ DE L'INDEX")
    print("=" * 80)
    print(f"Total fichiers : {len(index['files'])}")
    print(f"Hémisphère nord : {len(index['by_hemisphere']['north'])} fichiers")
    print(f"Hémisphère sud : {len(index['by_hemisphere']['south'])} fichiers")
    print(f"Plages RA uniques : {len(index['ra_ranges'])}")
    
    if args.mag_limit:
        mag_key = f"{args.mag_limit:.1f}"
        if mag_key in index['by_mag']:
            print(f"Magnitude {args.mag_limit} : {len(index['by_mag'][mag_key])} fichiers")
    
    print(f"\n✅ Index sauvegardé : {index_path}")
    print("\n💡 Pour reprendre le téléchargement, utilisez l'onglet Réduction avec")
    print("   l'option 'Reprendre téléchargement' activée (elle le sera automatiquement).")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
