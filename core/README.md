# NPOAP

Application de traitement et d’analyse pour l’astrophotométrie et les transitoires (réduction CCD, astrométrie, photométrie, courbes de lumière, etc.).

## Configuration locale

Copiez **`config.example.json`** vers **`config.json`**, puis renseignez observatoire, matériel et clé API Astrometry.net. Le fichier **`config.json`** est ignoré par Git (ne pas le pousser sur GitHub).

## Installation

Voir **[README_INSTALLATION.md](README_INSTALLATION.md)** et la documentation dans le dossier `docs/`.

## Prérequis

- Windows recommandé pour le flux d’installation automatisé ; parties du code utilisent **WSL** pour certains outils (ex. astrometry.net local).
- Python **3.11.x** et dépendances listées dans `requirements.txt`.

## Lancement

Après installation, utiliser le script généré **`LANCER_NPOAP.bat`** ou exécuter `main.py` dans l’environnement configuré.

## Licence et crédits

Consulter `docs/ACKNOWLEDGEMENTS.md` et les fichiers `LICENSE` éventuels des sous-projets tiers (ex. `external_apps/`).  
Parmi les dépendances récentes pour la chaîne TTV / mid-times : **[Wōtan](https://github.com/hippke/wotan)** (détrendage optionnel dans le viewer d’affinage Tc), avec citation scientifique indiquée dans les remerciements.
