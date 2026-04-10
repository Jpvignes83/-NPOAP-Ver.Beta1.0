# Manuel Utilisateur NPOAP

**NPOAP - Nouvelle Plateforme d'Observation et d'Analyse Photométrique**

Version 1.0

**Responsable HOPS-modified** : J.P Vignes  
**Contact** : jeanpascal.vignes@gmail.com

---

## Table des matières

1. [Vue d'ensemble](#1-vue-densemble)
2. [Accueil](#2-accueil)
3. [Réduction de Données](#3-réduction-de-données)
4. [Photométrie Exoplanètes (HOPS intégré)](#4-photométrie-exoplanètes)
5. [Photométrie Astéroïdes](#5-photométrie-astéroïdes)
   - [5.1 Astrométrie Zero-Aperture](#51-astrométrie-zero-aperture)
   - [5.2 Courbe de lumière, modèle harmonique, détrendage et ALCDEF](#52-courbe-de-lumière-modèle-harmonique-détrendage-et-alcdef)
   - [5.2a Surimpression modèle DAMIT](#52a-surimpression-modèle-damit)
   - [5.3 Onglet Publication (jeu validé, S-ALCDEF)](#53-onglet-publication-jeu-validé-s-alcdef)
6. [Photométrie Transitoires](#6-photométrie-transitoires)
7. [Analyse de Données](#7-analyse-de-données)
   - [7.1 Sous-onglet A : Détermination Période](#71-sous-onglet-a--détermination-période)
   - [7.2 Sous-onglet B : Recherche & Analyse TTV](#72-sous-onglet-b--recherche--analyse-ttv)
   - [7.3 Sous-onglet C : Analyse Système Multiple](#73-sous-onglet-c--analyse-système-multiple)
   - [7.4 Sous-onglet D : Simulation N-body](#74-sous-onglet-d--simulation-n-body)
8. [Étoiles Binaires](#8-étoiles-binaires)
9. [Easy Lucky Imaging](#9-easy-lucky-imaging)
   - [9.1 Mesure de séparation](#91-mesure-de-séparation)
10. [Analyse des outils et fichiers compatibles LcTools (TESS)](#10-analyse-des-outils-et-fichiers-compatibles-lctools-tess)
11. [Spectroscopie](#11-spectroscopie)
12. [Observation de la nuit](#12-observation-de-la-nuit)
13. [Conseils généraux](#13-conseils-généraux)

---

## 1. Vue d'ensemble

NPOAP est une application complète de **réduction, photométrie et analyse** d'observations astronomiques. Elle permet de traiter des observations d'exoplanètes, d'astéroïdes, de transitoires, d'étoiles binaires et de spectroscopie, d'analyse de données.

L'interface principale contient plusieurs onglets dans un notebook :

- **🏠 Accueil** : Configuration de base et calculateur d'échelle de pixel
- **🛠️ Réduction de Données** : Traitement des images brutes = calibration et plate-solving
- **🔭 Photométrie Exoplanètes** : Chaîne **HOPS** intégrée (réduction → photométrie → ajustement transit dans HOPS ; analyse avancée possible dans *Analyse de Données*)
- **🛰️ Photométrie Astéroïdes** : Photométrie et astrométrie d'astéroïdes
- **💥 Photométrie Transitoires** : Analyse d'événements transitoires
- **📈 Analyse de Données** : Outils d'analyse avancés (périodogrammes, TTV, N-body)
- **⭐ Étoiles Binaires** : Modélisation de systèmes binaires
- **🌟 Easy Lucky Imaging** : Traitement d'images d'étoiles doubles et mesure de séparation
- **📚 Catalogues** : Extraction et gestion de catalogues astronomiques (étoiles, astéroïdes, comètes, binaires, exoplanètes)
- **🔬 Spectroscopie** : Analyse de spectres
- **🌙 Observation de la nuit** : Éphémérides pour la nuit (astéroïdes, exoplanètes, comètes, binaires à éclipses), graphique d’altitude et export NINA

### Workflow recommandé

1. **Réduction de données** → Calibrer les images brutes, astrométrie, alignement WCS (optionnel), empilement (optionnel)
2. **Photométrie** → Extraire les courbes de lumière
3. **Analyse de données** → Analyser les périodes, TTV, systèmes multiples



---

## 2. Accueil

L'onglet **Accueil** est le point d'entrée de l'application. Il permet de configurer les paramètres de base nécessaires au fonctionnement de NPOAP.

### Configuration de l'observatoire

Dans la section **Observatoire**, vous devez renseigner :

- **Nom de l'observatoire** : Le nom officiel de votre observatoire
- **Latitude** : Coordonnée géographique en degrés décimaux (exemple : 48.8566 pour Paris)
- **Longitude** : Coordonnée géographique en degrés décimaux (exemple : 2.3522 pour Paris)
- **Élévation** : Altitude de l'observatoire en mètres

Ces informations sont utilisées pour les calculs d'astrométrie et de photométrie.

### Clé API Astrometry.net

Dans la section **Clé API Astrometry.net**, vous pouvez :

- **Obtenir une clé** : 
  - Vous devez obtenir une clé API gratuite sur le site web d'Astrometry.net : https://nova.astrometry.net/
  - Créez un compte sur le site, puis demandez une clé API dans votre espace utilisateur
  - Une fois la clé obtenue, utilisez le bouton "Générer/Modifier" dans NPOAP pour l'entrer

- **Entrer une clé API** : 
  - Cliquez sur le bouton **"Générer/Modifier"** pour ouvrir une fenêtre de dialogue
  - Entrez votre clé API Astrometry.net dans le champ (la clé sera masquée pour la sécurité)
  - La clé sera sauvegardée automatiquement

### Calculateur d'échelle de pixel

Dans la section **Calculateur d'échelle de pixel**, vous pouvez calculer l'échelle de pixel de votre système d'observation :

1. **Taille du pixel (µm)** : Entrez la taille d'un pixel de votre capteur en micromètres (ex: 3.76)
2. **Focale (mm)** : Entrez la longueur focale de votre télescope en millimètres (ex: 1939)
3. **Calculer** : Cliquez sur le bouton pour calculer l'échelle de pixel en secondes d'arc par pixel
4. **Sauvegarder dans config.py** : Cliquez sur ce bouton pour sauvegarder les valeurs dans le fichier de configuration

L'échelle de pixel calculée sera automatiquement utilisée dans l'onglet **Easy Lucky Imaging** pour les mesures de séparation angulaire.

**Formule utilisée** : `Échelle de pixel (arcsec/pixel) = (Taille pixel en mm / Focale en mm) × 206265`

### Sauvegarde

Cliquez sur le bouton **Sauvegarder la configuration** pour enregistrer vos paramètres. Ces informations seront conservées pour les prochaines sessions.

---

## 3. Réduction de Données

L'onglet **Réduction de Données** permet de traiter les images brutes acquises par votre télescope pour obtenir des images scientifiques prêtes à l'analyse.

### Structure des données d'entrée

Organisez vos données dans des dossiers séparés :

- **Bias** : Images de biais (exposition minimale)
- **Dark** : Images de noir (obturateur fermé)
- **Flat** : Images de champ uniforme (images d'étalonnage)
- **Light** : Images scientifiques (objets d'intérêt)

### Structure des dossiers créés automatiquement

Lorsque vous définissez un répertoire de travail, NPOAP crée automatiquement la structure suivante :

```
répertoire_de_travail/
  ├── science/              # Images science finales (créé automatiquement)
  │   └── aligned/          # Images alignées (créé automatiquement)
  └── output/               # Dossier de sortie (créé automatiquement)
      ├── calibrated/       # Images calibrées (créé automatiquement)
      ├── astrometry/       # Résultats d'astrométrie (créé automatiquement mais utilisé que pour NOVA-astrometry.net)
      ├── master_bias.fits  # Master bias (créé lors de la calibration)
      ├── master_dark.fits  # Master dark (créé lors de la calibration)
      └── master_flat.fits  # Master flat (créé lors de la calibration)
```

### Processus de réduction détaillé

#### Étape 1 : Définir le répertoire de travail

1. Cliquez sur **"📁 Définir Répertoire"**
2. **Important** : Sélectionnez le répertoire **où se trouvent vos images brutes** 
   - Ce répertoire servira de base pour votre projet
   - Les images brutes sont dans ce répertoire 
3. NPOAP crée automatiquement les sous-dossiers dans ce répertoire :
  - `science/` → images résolues (astrométrie locale) et emplacement recommandé pour l’image empilée (Master)
  - `output/` → pour toutes les sorties
  - `output/calibrated/` → pour les images calibrées
  - `output/astrometry/` → pour les résultats d'astrométrie avec NOVA-astrometry.net seulement
  - `science/aligned/` → pour les images alignées (après « Aligner Images (WCS) »)

#### Étape 2 : Charger les fichiers

1. **Charger Lights** : Sélectionnez les fichiers contenant vos images scientifiques
2. **Charger Bias** : Sélectionnez les fichiers contenant vos images de biais
3. **Charger Darks** : Sélectionnez les fichiers contenant vos images de dark
4. **Charger Flats** : Sélectionnez les fichiers contenant vos images de flat
5. **Option — Scaler les darks** : Cochez **« Scaler les darks au temps d'exposition des lights (si différent) »** si le temps d'exposition de vos darks ne correspond pas exactement à celui des images science. NPOAP extrapolera alors le master dark en le multipliant par le rapport (temps d'exposition de la light / temps d'exposition des darks), comme dans AstroImageJ. Les en-têtes FITS doivent contenir **EXPTIME** (ou **EXPOSURE**) pour les lights et les darks.

#### Sécurités avant calibration

Avant de créer les masters et de calibrer, NPOAP vérifie automatiquement :

- **Binning** : Les images de calibration (bias, darks, flats) doivent avoir le **même binning** que les images science. Les en-têtes FITS lus sont **XBINNING** / **YBINNING**, **BINAXIS1** / **BINAXIS2** ou **CCDSUM** (ex. `2x2`). En cas d’incohérence, la calibration est refusée avec un message indiquant le fichier et les binning en cause.
- **Filtre** : Les **flats** doivent avoir un **filtre compatible** avec les lights : l’ensemble des filtres présents dans les flats doit être inclus dans l’ensemble des filtres des images science. Les clés utilisées sont **FILTER**, **FILT**, **FILTER1** ou **INSFLTNM**. En cas d’incohérence (ex. flat en filtre R alors que les lights sont en V), la calibration est refusée.

#### Étape 3 : Création des masters

Lorsque vous lancez la calibration, les masters sont créés automatiquement :

1. **Master Bias** (`output/master_bias.fits`) :
   - Combine toutes les images de biais par médiane

2. **Master Dark** (`output/master_dark.fits`) :
   - Pour chaque image dark, soustrait le master bias
   - Combine toutes les images dark corrigées par médiane
   - Si **« Scaler les darks »** est coché : le master dark n’est pas modifié ; le scaling est appliqué *par image science* à l’étape 4 (voir ci-dessous).

3. **Master Flat** (`output/master_flat.fits`) :
   - Pour chaque image flat :
     - Soustrait le master bias
     - Soustrait le master dark
   - Combine toutes les images flat corrigées par médiane
   - Normalise (divise par la moyenne)

#### Étape 4 : Calibration des images Light

Pour chaque image scientifique :

1. **Soustraction du master bias** : `image = image - master_bias`
2. **Soustraction du master dark** :
   - **Sans scaling** : `image = image - master_dark`
   - **Avec scaling** (case cochée) : `image = image - master_dark × (EXPTIME_light / EXPTIME_darks)`. Le temps de référence des darks est la médiane des **EXPTIME** des images dark. Si une image science n’a pas d’**EXPTIME** valide, le dark est appliqué sans scaling pour cette image (un avertissement est enregistré dans le journal).
3. **Division par le master flat** : `image = image / master_flat`
4. **Sauvegarde** : L'image calibrée est sauvegardée dans `output/calibrated/`

**Important** : Les images calibrées sont automatiquement sauvegardées dans le dossier `output/calibrated/`. Ce dossier est ensuite utilisé pour l'astrométrie et les étapes suivantes.

#### Étape 5 : Astrométrie 

L'astrométrie permet d'ajouter les informations de coordonnées célestes (WCS) dans les en-têtes FITS des images calibrées. 

**Prérequis** :
- Les images doivent être **calibrées** (étape 4 terminée)
- Les images calibrées se trouvent dans `output/calibrated/`

**Catalogues utilisés** :
- **NOVA** : Le serveur Astrometry.net utilise des index internes (dérivés de catalogues tels que **Tycho-2** ou **Tycho-2 + Gaia-DR2**, selon l’échelle du champ). Aucune configuration de catalogue côté NPOAP.
- **Astrométrie locale (WSL)** : Le catalogue effectif est celui des **fichiers d’index** que vous avez installés pour *solve-field* sous WSL (souvent index Tycho-2 ou Tycho-2+Gaia-DR2).

**Deux méthodes disponibles** :

1. **🌐 Astrométrie en ligne (NOVA)** :
   - **Quand l'utiliser** : Si vous avez une connexion internet et une clé API Astrometry.net (configurée dans l'onglet Accueil)
   - **Avantages** : Rapide, simple, pas d'installation supplémentaire
   - **Processus** :
     - Cliquez sur **"🌐 Via Astrometry.net (NOVA)"**
     - Les images calibrées du dossier `output/calibrated/` sont envoyées au serveur Astrometry.net
     - Les images astrométrées sont sauvegardées dans `output/astrometry/`

2. **🖥️ Astrométrie locale (WSL)** :
   - **Quand l'utiliser** : Si vous n'avez pas de connexion internet, ou si vous préférez un traitement local
   - **Avantages** : Fonctionne hors ligne, contrôle total
   - **Prérequis** : WSL (Windows Subsystem for Linux) et solve-field installé
   - **Processus** :
     - Cliquez sur **"🖥️ Astrométrie Locale (WSL)"**
     - Les images calibrées du dossier `output/calibrated/` sont traitées localement
     - Les images astrométrées sont sauvegardées dans `science/`

#### Étape 6 : Alignement (optionnel)

1. **Source** : Les images résolues (avec WCS) sont dans `science/` (après astrométrie locale) ou dans `output/astrometry/` (après NOVA).
2. **Alignement WCS** : Cliquez sur **« 📐 Aligner Images (WCS) »**. NPOAP utilise la bibliothèque **reproject** pour reprojeter chaque image sur le WCS de la première image (référence), ce qui supprime les décalages entre poses. **Prérequis** : le module `reproject` doit être installé (environnement **astroenv** recommandé).
3. **Sortie** : Les images alignées sont enregistrées dans **`science/aligned/`** (les images dans `science/` ne sont pas modifiées).

#### Étape 7 : Empilement (optionnel)

1. Cliquez sur **« 📚 Empiler Images (Stack) »**.
2. Sélectionnez les images à empiler (par défaut le dialogue ouvre **`science/aligned/`** si vous avez lancé l’alignement, sinon **`science/`**).
3. Lors de l’enregistrement du Master, le dossier proposé par défaut est **`science/`** : vous pouvez y enregistrer directement l’image empilée (ex. `Master_Stacked.fits`).

### Tutoriel pas à pas

1. **Définir le répertoire** :
   - Cliquez sur "📁 Définir Répertoire"
   - **Sélectionnez le répertoire où se trouvent vos images brutes** (ex: `C:/Observations/2024-01-15/`)
   - NPOAP crée automatiquement : `science/`, `output/`, `output/calibrated/`, etc.

2. **Charger les fichiers** :
   - Cliquez sur "📂 Charger Lights" et sélectionnez vos images scientifiques
   - Répétez pour Bias, Darks, Flats
   - Si les temps d'exposition des darks diffèrent de ceux des lights, cochez **« Scaler les darks au temps d'exposition des lights (si différent) »**

3. **Lancer la calibration** :
   - Cliquez sur **"🚀 Lancer Calibration"**
   - NPOAP vérifie le binning (lights vs bias/darks/flats) et le filtre (lights vs flats) ; en cas d’erreur, un message indique le problème
   - Les masters sont créés dans `output/` : `master_bias.fits`, `master_dark.fits`, `master_flat.fits`
   - Les images calibrées sont sauvegardées dans `output/calibrated/`

4. **Astrométrie** 
   - **Astrométrie en ligne (NOVA)** :
     - Assurez-vous d'avoir configuré votre clé API Astrometry.net dans l'onglet Accueil
     - Cliquez sur **"🌐 Via Astrometry.net (NOVA)"**
     - Les images astrométrées sont sauvegardées dans `output/astrometry/`
   - **Astrométrie locale** :
     - Cliquez sur **"🖥️ Astrométrie Locale (WSL)"** (nécessite WSL et solve-field installés)
     - Les images astrométrées sont sauvegardées dans `science/`

5. **Utiliser les images calibrées** :
   - Pour la photométrie, utilisez les images du dossier `astrometry/` ou `science/` (après astrométrie)
  

---

## 4. Photométrie Exoplanètes

Ce chapitre décrit l'onglet **Photométrie Exoplanètes** de NPOAP, qui intègre **HOPS** (HOlomon Photometric Software, v3.3.x) pour le traitement photométrique des transits. Il inclut un **résumé** du manuel utilisateur officiel HOPS (document *hops33_manual_en.pdf*, v3.3.3) puis les **spécificités NPOAP** (intégration, extensions et dépannage).

### 4.1 Résumé du manuel officiel HOPS (référence ExoWorlds Spies)

Le manuel PDF officiel présente HOPS comme une chaîne en **six étapes** pour analyser des séries d'images FITS d'exoplanètes. Points essentiels rappelés :

- **Prérequis** : Windows, macOS ou Linux ; au moins 4 Go de RAM et d'espace disque ; **Python 3.7+** (3.8+ recommandé), idéalement via Anaconda.
- **Installation** : télécharger *hops-master.zip* depuis le site ExoWorlds Spies, décompresser et lancer l'installateur adapté au système ; un raccourci « hops » peut être créé sur le bureau.
- **Avant de lancer** : calibrations **brutes** (pas de masters préfabriqués) — au moins ~5 bias, ~5 darks (temps d'exposition aligné sur la science), ~5 flats (comptes ~2/3 du plein puits si possible) ; **tous les fichiers** (science + calibrations) dans **un seul dossier** sans sous-dossiers ; **identifiants de noms distincts** pour bias, dark, flat et science (ex. préfixes différents).
- **Étape 1 — Data & Target** : choix du répertoire de travail, vérification des identifiants de fichiers et des mots-clés FITS (temps d'exposition, date/heure, filtre, estampille temporelle), choix de la cible (header, SIMBAD en ligne, ou saisie RA/Dec), optionnellement lieu et infos observateur.
- **Étape 2 — Reduction** : lancement automatique après l'étape 1 ; construction des masters et réduction des images science ; possibilité d'afficher toutes les images (plus lent) ; arrêt via le bouton dédié.
- **Étape 3 — Inspect Frames** (filtre recommandé) : graphiques SKY et PSF pour repérer nuages, crépuscule, mauvais guidage ; exclusion/inclusion d'images ; la **première image** doit être représentative (surexposition, tracking).  
  Dans HOPS-modified, un fichier `quality_criteria.txt` est généré dans le dossier de réduction, trié du **plus mauvais** au **meilleur** frame, et les ~10 % les plus dégradées sont mises en évidence dans l'inspection.
- **Étape 4 — Alignment** : alignement des images ; en cas d'échec de détection d'étoiles, HOPS propose d'exclure l'image si elle est défectueuse.
- **Étape 5 — Photometry** : sélection de la **cible** et d'au moins **deux étoiles de comparaison** dans le champ (carré rouge = FOV utile ; rectangles jaunes = flux proche de la cible) ; critères : proximité, magnitude et **couleur (Gbp−Grp)** comparables, stabilité (**non variables**) ; réglage des rayons d'ouverture et anneaux de ciel ; **RUN PHOTOMETRY** puis inspection des courbes (rapports cible/comparaisons normalisés) ; sauvegarde des résultats avant l'étape suivante.
- **Étape 6 — Exoplanet Fitting** : choix du fichier de photométrie (ouverture ou PSF gaussienne) ; paramètres planète (catalogue intégré ou saisie manuelle) ; **RUN TEST** (aperçu détrendu, modèle rouge vs attendu cyan, résidus) puis **RUN FITTING** après validation.

Pour le détail des captures d'écran et des cas limites, se reporter au **PDF officiel** ; la suite du présent chapitre documente **uniquement** le fonctionnement **dans NPOAP**.

### 4.2 Intégration de HOPS dans NPOAP

- **Emplacement** : après installation depuis le ZIP, HOPS est attendu sous `external_apps/hops/hops-master` (voir [Manuel d'installation](MANUEL_INSTALLATION.md)).
- **Interface** : la fenêtre principale HOPS est **embarquée** dans l'onglet (zone avec **défilement** si la grille ou les figures dépassent la hauteur du panneau).
- **Boutons NPOAP** :
  - **Installer / Réinstaller HOPS (ZIP)** : déploie ou met à jour les sources depuis **`HOPS-modified.zip`** (prioritaire s'il est présent dans `external_apps/hops/`) et applique les correctifs nécessaires à l'intégration.
  - **Lancer HOPS** : démarre HOPS dans le cadre réservé.
  - **Réinitialiser HOPS** : recrée le conteneur d'affichage et relance une instance propre (utile après incident ou après mise à jour).
- **Licence** : en bas du cadre, le texte **MIT License** et le copyright **Angelos Tsiaras and Konstantinos Karpouzas** (2017) sont affichés conformément à HOPS.

### 4.3 Comportement particulier sous NPOAP (fenêtres, isolation, journaux)

- **Fenêtres secondaires** : les sous-fenêtres HOPS (Data & Target, Reduction, Photometry, etc.) s'ouvrent en **Toplevel** ; NPOAP applique une politique pour les garder **au premier plan** (transient + topmost) afin qu'elles ne passent pas derrière la fenêtre principale. Les boîtes **showinfo** / **askyesno** suivent la même logique pour éviter un blocage apparent.
- **Isolation** : les callbacks HOPS restaurent le répertoire de travail global après exécution ; l'import HOPS ne laisse pas de modification permanente de `sys.path`. Les styles **ttk** de NPOAP ne doivent pas être modifiés par le thème Combobox de HOPS (style local **HOPS.TCombobox**).
- **Journal d'exécution** : les messages d'erreur et exceptions utiles au débogage sont écrits dans  
  `external_apps/hops/hops-master/runtime_logs/hops_runtime.log`  
  (niveau par défaut orienté **erreurs** ; voir le manuel d'installation).

### 4.4 Extensions et correctifs appliqués au code HOPS dans NPOAP

Les points suivants complètent ou adaptent le comportement du manuel PDF :

| Thème | Description |
|--------|-------------|
| **Réduction** | Résolution des chemins FITS en **chemins absolus** par rapport au répertoire de travail HOPS (évite un blocage à 0 % si le répertoire courant de NPOAP diffère du dossier des données). Message indiquant que la **barre « science »** n'avance qu'après la phase bias/dark/flat. |
| **Inspection — qualité des frames** | HOPS-modified calcule des critères de qualité (score composite) et écrit `quality_criteria.txt` dans le dossier de réduction. Le fichier est trié du pire au meilleur ; la fenêtre **Inspection** met en évidence les ~10 % de frames les plus problématiques pour accélérer le tri. |
| **Réduction — masters** | En fin de calcul, HOPS-modified enregistre les masters dans le dossier de réduction : `master_bias.fits`, `master_dark.fits`, `master_darkf.fits`, `master_flat.fits` (si les jeux correspondants sont fournis). |
| **Photometry — cible** | Bouton **Use OBJCTRA/OBJCTDEC** : lecture des coordonnées dans le header FITS, mise à jour de la cible dans le journal HOPS et positionnement de T1 si un **WCS** est disponible ; sinon invitation à lancer le **plate solve**. |
| **Photometry — Gaia DR3** | Après **PLATE SOLVE IMAGE**, HOPS interroge **astroquery** / table **gaiadr3.gaia_source** (positions et magnitudes). La colonne **Gbp−Grp** peut s'appuyer sur BP/RP ; une colonne **diagnostic** indique si le plate solve, le rapprochement ou les magnitudes BP/RP manquent. |
| **Photometry — variables** | Colonne **Var (Gaia DR3)** : affichage du **phot_variable_flag** (ex. **VAR** si VARIABLE, **const** si CONSTANT, **n/a** si non disponible). Pour une comparaison marquée variable, les **cercles** et le libellé **C*n*** **clignotent** (cyan / orange) pour signaler le risque photométrique. |
| **Filtres Gaia (passbands)** | Au lancement de HOPS depuis NPOAP, les passbands `GAIA_G`, `GAIA_BP`, `GAIA_RP` sont chargés automatiquement depuis `resources/filters/` (fichiers `.txt`) via `pylightcurve`, si ces fichiers sont présents. |
| **RUN PHOTOMETRY** | Meilleur retour utilisateur pendant le traitement ; garde-fous pour les appels réseau lents (ex. **timeout** sur services externes en fin de chaîne) afin de limiter les blocages ; boîtes de dialogue de fin de traitement **devant** la fenêtre de progression. |

### 4.5 Ancien volet NPOAP « Exoplanètes »

L'ancien panneau gauche NPOAP de cet onglet (image de référence locale, workflow T1 / ouvertures / batch / analyse / rapports dédiés) a été **retiré de l'interface**. Les modules Python sous-jacents peuvent rester présents pour d'autres usages dans le projet ; le parcours **recommandé** pour la photométrie exoplanète dans cet onglet est **exclusivement via HOPS**.

### 4.6 Rappels pratiques

1. Installer ou mettre à jour HOPS depuis le ZIP, puis **Lancer HOPS** (premier lancement : installation possible des dépendances Python, voir manuel d'installation).
2. Suivre les **six étapes** HOPS comme dans le manuel PDF ; en cas de problème réseau, vérifier la connexion pour le plate solve et les requêtes Gaia.
3. Consulter `hops_runtime.log` (dans le dossier ci-dessus) si une étape échoue sans message visible suffisant.
4. Pour l'**ajustement de transit** détaillé avec les graphiques Chi², Shapiro-Wilk, etc., utiliser l'onglet **[Analyse de Données](#7-analyse-de-données)** de NPOAP sur les produits exportés par votre chaîne ; ce n'est pas décrit ici comme prolongement direct du cadre HOPS intégré.

### 4.7 Choix du filtre (important)

Le filtre choisi dans **HOPS > Data & Target** impacte toute la chaîne :

- la cohérence de calibration (lights/flats),
- la comparaison avec les catalogues (couleur, diagnostics),
- le calcul des coefficients de **limb-darkening** au fitting.

#### Recommandations pratiques

- Utilisez le filtre réellement employé à l'acquisition (`R`, `V`, `I`, etc.).
- Si vous utilisez un filtre non standard, choisissez la bande la plus proche pour commencer, puis vérifiez les résidus au RUN TEST.
- Évitez de mélanger des séries prises avec des filtres différents dans une même courbe de lumière.

#### Filtres Gaia dans NPOAP/HOPS

Les filtres `Gaia G`, `Gaia BP`, `Gaia RP` sont disponibles dans la liste.

- Si les passbands Gaia sont installés dans l'environnement, le limb-darkening est calculé nativement (`GAIA_G/BP/RP`).
- Sinon, NPOAP applique automatiquement un **fallback** pour préserver le fitting :
  - `Gaia G` -> `JOHNSON_V`
  - `Gaia BP` -> `JOHNSON_B`
  - `Gaia RP` -> `COUSINS_R`

Ce fallback permet de terminer le traitement, mais la solution la plus fidèle reste d'installer les passbands Gaia personnalisés (voir `MANUEL_INSTALLATION.md`).

---

## 5. Photométrie Astéroïdes

L'onglet **Photométrie Astéroïdes** permet de mesurer la magnitude des astéroïdes, de réaliser de l'astrométrie précise et de générer des rapports au format ADES (Asteroid Data Exchange Standard).

### Chargement d'images

1. **Ouvrir un dossier** : Utilisez le bouton "📁 Charger images" pour charger un dossier contenant des images FITS calibrées.
2. **Navigation** : Utilisez le slider ou les boutons de navigation (⏮ ◀ ▶ ⏭) pour parcourir les images
3. **Défilement automatique** : Utilisez le bouton ▶ pour activer/désactiver le défilement automatique

### Récupération des éphémérides

1. **ID Cible** : Entrez le numéro ou le nom de l'astéroïde (ex: "433" pour Eros)
   - **Formats acceptés par JPL Horizons** :
     - Numéro d'astéroïde : `433` (Eros), `1` (Cérès), `2021` (numéros MPC)
     - Désignation provisoire : `2021 AB` (format MPC)
     - Nom officiel : `Eros`, `Ceres` (noms reconnus par JPL Horizons)
     - Comètes : `C/2020 F3 (NEOWISE)` ou numéro de comète
   - **Si l'objet n'a pas d'ID ou n'est pas trouvé** :
     - **ID vide** : Un message d'erreur s'affiche : "ID cible vide"
     - **Objet introuvable** : L'application tente plusieurs méthodes de récupération :
       * Tentative 1 : Requête standard avec période de 3 heures autour de l'observation
       * Tentative 2 : Requête avec paramètres minimaux (RA/Dec + magnitude uniquement)
       * Tentative 3 : Période réduite à 30 minutes
       * Tentative 4 : Requête pour une seule date précise
     - Si toutes les tentatives échouent, un message d'erreur détaillé est affiché
     - **Causes possibles d'échec** :
       * ID incorrect ou inconnu dans la base JPL Horizons
       * Nom d'objet mal orthographié
       * Objet trop récent (pas encore dans la base de données)
       * Problème de connexion internet (requête vers JPL Horizons)
   - **Solutions** :
     * Vérifiez l'ID sur le site [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) avant de l'utiliser
     * Pour les objets récents, utilisez la désignation provisoire complète du MPC
     * Assurez-vous d'avoir une connexion internet active
     * Vérifiez que la date d'observation (DATE-OBS dans le header FITS) est valide
2. **Code Obs.** : Entrez votre code observatoire MPC (ex: "500" pour Greenwich)
3. **Pas Éphémérides** : Choisissez le pas d'éphémérides :
   - **1m** = ~180 points pour 3h (très précis)
   - **2m** = ~90 points pour 3h (recommandé)
   - **5m** = ~36 points pour 3h (rapide)
   - **10m** = ~18 points pour 3h (moins précis)
4. **Récupérer** : Cliquez sur "🔭 Récupérer Éphémérides"

### Éphémérides et ancres : quelle position pour la collecte et les rapports ?

**Ordre de priorité dans NPOAP** pour remplir la position de la cible lors de la **collecte des observations** (messages de journal du type « Position T1 pour … », champs utiles aux rapports ADES, etc.) :

1. **Table d’éphémérides chargée** — Après « 🔭 Récupérer Éphémérides », si l’interpolation au **JD du milieu d’exposition** de chaque image réussit, la **RA/Dec** utilisée pour cette collecte est **uniquement celle issue de la table** (JPL Horizons). Les coordonnées **pixel** éventuellement affichées dans le journal sont obtenues en projetant cette position avec le **WCS** du fichier FITS.
2. **Ancres manuelles** — Si **aucune** position n’a pu être fixée par les éphémérides (pas de récupération, ou échec du calcul), mais que vous avez enregistré une position T1 sur la **première** et la **dernière** image de la série, la position sur les autres images est **interpolée** entre ces deux ancres (interpolation linéaire sur la sphère en fonction du JD), comme décrit plus bas pour les objets sans ID.

**Point important** : tant que les éphémérides sont **chargées et utilisées avec succès**, les **ancres ne remplacent pas** cette trajectoire dans la collecte. Si l’éphéméride est **fausse** (mauvais numéro MPC, mauvaise désignation, orbite inadaptée, confusion d’objet, etc.), les RA/Dec enregistrés pour les rapports suivront cette erreur, **même** si le pointage visuel ou les ancres sur la vraie tache seraient meilleurs.

**Recommandation** :

- Si vous constatez que la trajectoire modèle **ne colle pas** à l’objet sur les images : **n’utilisez pas** la table d’éphémérides pour cette session (ne pas lancer « Récupérer Éphémérides », ou enchaîner une procédure sans cette étape selon votre flux), puis placez **deux ancres** fiables sur la **vraie** cible (première et dernière image, ou deux poses bien séparées dans le temps).
- À l’inverse, une éphéméride **correcte** reste en général **préférable sur une longue série**, car deux ancres ne font qu’approximer le mouvement par une **ligne** sur la sphère entre deux instants ; sur quelques heures c’est souvent suffisant si les clics sont bons.

**Affichage à l’écran (cercle T1)** : à l’ouverture d’une image, le marqueur peut aussi provenir des mots-clés FITS **`OBJCTRA` / `OBJCTDEC`** (y compris sexagésimal avec **espaces**, ex. `06 57 57`, ou avec **deux-points**) ou d’un **clic manuel**. Ce marqueur **n’est pas** la même chose que la position « rapport » : un décalage entre cercle et log « Position T1 pour … » incite à vérifier séparément **éphéméride**, **header** et **pointage**.

### Workflow pour objets sans ID (nouvelle découverte)

Si vous observez un objet qui n'a **pas encore d'ID MPC** ou qui n'est pas dans la base JPL Horizons, vous pouvez suivre ce workflow :

1. **Sauter la récupération d'éphémérides** :
   - Ne remplissez pas le champ "ID Cible" ou laissez-le vide
   - Cliquez directement sur "🧭 Lancer l'Astrométrie" ou "🧭 Batch Astrométrie"

2. **Effectuer l'astrométrie** :
   - L'astrométrie fonctionne **indépendamment des éphémérides**
   - Elle nécessite seulement un WCS initial dans les headers FITS
   - L'astrométrie sera effectuée normalement sur toutes vos images

3. **Sélectionner T1 sur la première et la dernière image** :
   - Allez sur la **première image** de la série, puis **cliquez-gauche sur l'objet** pour désigner T1.
   - Allez sur la **dernière image**, puis **cliquez-gauche sur l'objet** pour désigner T1 à nouveau.
   - Ces deux positions servent d’« ancres » : le logiciel en déduit la position de T1 sur les images intermédiaires si besoin.
   - Les coordonnées (RA, Dec) sont calculées via le WCS.

   **Comment la position de T1 est-elle calculée sur chaque image sans éphémérides ?**  
   Pour chaque image, le logiciel utilise une **interpolation linéaire** entre la première et la dernière ancre, en fonction du **temps d’observation** (JD) de l’image :
   - Il calcule la fraction de temps entre la première et la dernière image : `frac = (JD_image − JD_première) / (JD_dernière − JD_première)`.
   - Les deux positions ancres (RA, Dec) sont converties en coordonnées **cartésiennes** sur la sphère céleste (pour éviter les discontinuités en RA à 0°/360°).
   - La position interpolée est : position = position_première + frac × (position_dernière − position_première), puis reconversion en RA/Dec.
   - Cette position (RA, Dec) est transformée en pixels via le WCS de l’image, puis un **centroïde** local affine la position pour chaque mesure photométrique.

   En pratique, on suppose que l’objet se déplace de façon quasi linéaire sur le ciel entre la première et la dernière prise ; c’est une approximation raisonnable sur une série courte (quelques heures).

4. **Configurer la photométrie** :
   - Sur une image de référence, cliquez sur **"⭕ SET-UP PHOTOMÉTRIE"** pour définir T1/comparateurs et les apertures.
   - En cas de comète diffuse (FWHM instable), privilégiez le mode image par image (voir ci-dessous).

5. **Résultats photométriques**
   - Mode **comète (image par image)** :
     - Bouton **« Photométrie image »**
     - CSV cumulatif : **`results/<objet>_photometrie_image_par_image.csv`**
     - Colonnes : `filename, JD-UTC, date_obs, filter_used, delta-to_G, mag_T1_G, rmsMag_T1`
   - Mode **astéroïdes (batch)** :
     - Bouton **"📊 PHOTOMÉTRIE BATCH (Astéroïdes)"**
     - Résultats batch standards dans `photometrie/results.csv` (+ `light_curve.txt`)
     - Copie compilée exportée dans `results`.

6. **Sans éphémérides** :
   - Le suivi T1 reste possible via priorité :
     1. positions ZA0/astrométriques déjà calculées,
     2. ancres manuelles (première/dernière image),
     3. position fixe de référence (fallback).
   - La magnitude peut être calculée via comparateurs Gaia (si comparateurs valides).

7. **Après obtention de l'ID** :
   - Une fois que l'objet a reçu un ID MPC, vous pouvez relancer "🔭 Récupérer Éphémérides"
   - Les éphémérides permettront le positionnement automatique de T1 et compléteront les rapports

### Sélection de la cible (T1)

1. **Position automatique** (si éphémérides disponibles) : Si les éphémérides sont chargées, la position de l'astéroïde est calculée automatiquement (voir la section **Éphémérides et ancres** ci-dessus pour la priorité dans les **rapports**).
2. **Clic manuel** : Cliquez sur l'astéroïde dans l'image pour le sélectionner manuellement.
3. **Header FITS** : lors du chargement d’une image, le logiciel peut placer T1 à partir de **`OBJCTRA` / `OBJCTDEC`** (ou `RA` / `DEC` en secours), en interprétant les formats **sexagésimaux** classiques (séparateurs **:** ou **espaces**) ainsi que les **valeurs numériques** en degrés ou heures.

### Sélection des étoiles de comparaison

1. **Ouvrir la fenêtre de sélection** : Cliquez sur "Sélectionner les comps"
2. **Critères de sélection** :
   - Étoiles dans un rayon de 15 arcminutes
   - Magnitude appropriée (ni trop faible, ni saturée)
   - Exclusion automatique des étoiles variables connues (Gaia DR3)
3. **Validation** : Cliquez sur "Valider" pour confirmer votre sélection

### Photométrie

La photométrie peut être effectuée sur une image unique ou en mode batch sur toute la série. La précision de la photométrie dépend de la qualité de l'astrométrie préalable.

#### Méthodes d'astrométrie et leur impact sur la photométrie

Avant de lancer la photométrie, vous devez choisir une méthode d'astrométrie qui influence la précision du positionnement de T1 :

1. **Astrométrie Zero-Aperture (extrapolation)** :
   - **Principe** : Utilise plusieurs ouvertures photométriques et extrapole les positions astrométriques à zéro ouverture pour minimiser les biais
   - **Avantages** : 
     - **Précision maximale** (typiquement RMS < 0.1 arcsec)
     - Minimise les erreurs systématiques liées au rayon d'ouverture
     - Idéal pour des observations de haute précision
   - **Inconvénients** : 
     - Plus lent (nécessite plusieurs mesures par étoile)
     - Optimisé en batch (4 images sur 5 utilisent une extrapolation rapide)
   - **Quand l'utiliser** : Observations nécessitant une précision astrométrique maximale (astéroïdes, objets faibles, mesures de haute précision)

2. **Astrométrie Classique (FWHM)** :
   - **Principe** : Utilise une seule ouverture optimisée selon le FWHM de l'image
   - **Avantages** :
     - **Rapide** (une seule mesure par étoile)
     - Suffisant pour la plupart des observations
     - Précision typique : RMS ~ 0.2-0.3 arcsec
   - **Inconvénients** : 
     - Légèrement moins précis que zero-aperture
     - Peut avoir des biais systématiques pour certaines configurations
   - **Quand l'utiliser** : Observations standards, traitement rapide de grandes séries

#### Modes de photométrie

1. **Photométrie image par image (comètes ou contrôle fin)** :
   - **Bouton** : **« Photométrie image »**
   - **Fonctionnement** :
     - Traite uniquement l'image actuellement affichée
     - Utilise la configuration définie via **"⭕ SET-UP PHOTOMÉTRIE"**
     - Calcule `mag_T1_G` à partir des comparateurs Gaia
     - Met à jour un CSV cumulatif dans `results`
   - **Utilisation** :
     - Cas comètes / objets diffus (réglages fins image par image)
     - Contrôle qualité avant génération finale ADES
   - **Prérequis** :
     - Image chargée et affichée
     - Astrométrie effectuée (pour un positionnement précis de T1)

2. **Photométrie batch (Astéroïdes)** :
   - **Bouton** : "📊 PHOTOMÉTRIE BATCH (Astéroïdes)"
   - **Fonctionnement** :
     - Traite automatiquement **toutes les images** du dossier chargé
     - Utilise la configuration photométrique validée
     - Applique les mêmes apertures à toutes les images
     - Génère une courbe de lumière complète
   - **Utilisation** :
     - Pour traiter une série complète d'observations
     - Après avoir validé la configuration sur une image de référence
     - Pour générer les données photométriques finales
   - **Prérequis** :
     - Configuration photométrique validée (via photométrie image par image)
     - Astrométrie effectuée sur toutes les images (ou en batch)
     - Éphémérides facultatives (fallbacks automatiques disponibles)
   - **Options d'optimisation batch (Zero-Aperture uniquement)** :
     - **"Toutes les images"** : Extrapolation zero-aperture complète sur toutes les images
       - Précision maximale sur toute la série
       - Plus lent mais meilleure qualité
       - Recommandé pour observations de haute précision
     - **"Optimisé (1 sur 5)"** : Extrapolation complète seulement sur 1 image sur 5
       - Les autres images utilisent l'astrométrie classique (sans extrapolation)
       - **5 fois plus rapide** que le mode "Toutes les images"
       - Précision légèrement réduite mais généralement suffisante
       - Recommandé pour grandes séries d'observations
     - Cette option n'affecte que l'astrométrie zero-aperture en mode batch
     - L'astrométrie classique traite toujours toutes les images de la même manière

#### Workflow recommandé

1. **Préparation** :
   - Chargez le dossier d'images calibrées
   - Récupérez les éphémérides de l'astéroïde
   - Effectuez l'astrométrie (zero-aperture ou classique) sur au moins une image de référence

2. **Configuration photométrique** (image par image) :
   - Naviguez vers une image de bonne qualité (bon seeing, bon SNR)
   - Cliquez sur l'image pour sélectionner T1 (l'astéroïde)
   - Cliquez sur "⭕ SET-UP PHOTOMÉTRIE"
   - Configurez les apertures (rayon, annulus)
   - Sélectionnez les étoiles de comparaison
   - Validez la configuration

3. **Traitement photométrique** :
   - **Comètes / image par image** : utilisez **« Photométrie image »** et avancez dans la série
   - **Astéroïdes** : utilisez "📊 PHOTOMÉTRIE BATCH (Astéroïdes)"
   - Le CSV de photométrie est compilé automatiquement dans `results`

4. **Validation** :
   - Vérifiez les magnitudes `mag_T1_G` et leurs incertitudes
   - Vérifiez la courbe de lumière générée
   - Examinez les résidus photométriques
   - Si nécessaire, ajustez la configuration et relancez

### Génération du rapport ADES

1. **Finalisation ADES** : ouvrez le **tableau zéro-ouverture** puis cliquez sur :
   - **"🧾 Création des rapports ADES"**
2. **Rapports générés** :
   - `ADES final` (standard)
   - `ADES final ZA0` (extrapolation ouverture nulle)
3. **Photométrie intégrée** :
   - `mag_T1_G` et `rmsMag_T1` sont injectés dans ADES/ADES_ZA0 si disponibles
   - `logSNR` est calculé depuis `rmsMag_T1` avec :
     - `SNR ≈ 1.0857 / sigma_mag`
     - `logSNR = log10(SNR)`
4. **Exports** :
   - Les rapports finaux sont écrits dans le dossier `results`.

---

### 5.1 Astrométrie Zero-Aperture

L'astrométrie **zero-aperture** est une méthode avancée qui permet d'obtenir une astrométrie très précise en minimisant les biais liés au rayon d'ouverture photométrique. La méthodologie est alignée sur le projet **Zero-Aperture-Astrometry** (https://github.com/bensharkey/Zero-Aperture-Astrometry).

#### Principe

La méthode zero-aperture fonctionne en :

1. Mesurant la position des étoiles de référence avec **plusieurs rayons d'ouverture** (6 ouvertures)
2. Pour chaque aperture : calcul d’un WCS, puis **résidus RA et Dec** (écart calcul–catalogue) en arcsec pour chaque étoile
3. **Extrapolation à rayon = 0** :
   - **RMS** : fit (pondéré) du RMS total en fonction de l’aperture, extrapolation à 0 → RMS zero-aperture
   - **Décalage systématique** : fit des **résidus moyens RA et Dec** en fonction de l’aperture, extrapolation à 0 → décalage à corriger ; un **WCS corrigé** est produit en appliquant ce décalage au WCS
4. Calcul d’un **WCS optimisé** (WCS corrigé zero-aperture lorsque disponible) et de **statistiques détaillées** (RMS global, RMS RA/Dec, corrélation, outliers, offsets résiduels à aperture 0, etc.)

Cette méthode est plus précise que la méthode classique (une seule ouverture basée sur le FWHM) mais est également plus lente.

#### Interface

Dans l'onglet **Photométrie Astéroïdes**, la section **"Paramètres Astrométrie"** contient :

- **FWHM [px]** : FWHM moyen des étoiles en pixels
- **Seuil [σ]** : Seuil de détection des étoiles (sigma du bruit)
- **Gaia Gmax** : Magnitude limite pour les étoiles Gaia (plus grand = plus d'étoiles)
- **Match [″]** : Rayon de matching en arcsec entre détection et Gaia
- **GPU** : Option pour utiliser le GPU si disponible (CuPy)
- **Méthode** : Combobox pour choisir entre :
  - **"Zero-Aperture (extrapolation)"** → Mode zero-aperture (plus précis)
  - **"Classique (FWHM)"** → Mode simple et rapide

#### Tutoriel : Astrométrie Zero-Aperture (image unique)

1. **Charger les images**
   - Cliquez sur "📁 Charger images" et sélectionnez un dossier contenant des FITS calibrés avec WCS initial

2. **(Recommandé) Récupérer les éphémérides**
   - Entrez l'ID de l'astéroïde et le code observatoire
   - Réglez le pas d'éphémérides (2m recommandé)
   - Cliquez sur "🔭 Récupérer Éphémérides"

3. **Choisir la méthode Zero-Aperture**
   - Dans "Paramètres Astrométrie", vérifiez/ajustez FWHM, Seuil [σ], Gaia Gmax, Match [″]
   - Dans le combobox "Méthode", sélectionnez **"Zero-Aperture (extrapolation)"**
   - (Optionnel) Cochez "Utiliser GPU" si disponible

4. **Lancer l'astrométrie**
   - Assurez-vous d'être sur l'image voulue (navigation via slider/boutons)
   - Cliquez sur "🧭 Lancer l'Astrométrie"
   - Le processus :
     - Détecte les étoiles dans l'image
     - Les match avec le catalogue Gaia DR3
     - Mesure les positions avec plusieurs apertures
     - Extrapole à zero-aperture
     - Calcule le WCS optimisé

5. **Résultat**
   - Le WCS de l'image est **mis à jour** (solution ajustée avec Gaia)
   - Des **mots-clés astrométriques** sont ajoutés au header FITS :
     - `ASTREF` : Catalogue de référence (Gaia DR3)
     - `ASTRRMS` : RMS global (arcsec)
     - `ASTRRMSR` : RMS RA (arcsec)
     - `ASTRRMSD` : RMS Dec (arcsec)
     - `ASTRMSC` : RMS méthode classique (arcsec)
     - `ASTRMSZ` : RMS zero-aperture (arcsec)
     - `ASTNREF` : Nombre d'étoiles de référence
     - `ASTNOUT` : Nombre d'outliers (>3σ)
     - `ASTMETHOD` : Méthode utilisée (`zero-aperture`)
     - `ASTR2` : R² de la régression zero-aperture (si disponible)
   - En interne, les statistiques peuvent inclure les **offsets résiduels à aperture 0** (RA/Dec en arcsec) et le **WCS corrigé** (après application du décalage extrapolé), utilisé comme WCS final lorsque la méthode zero-aperture est retenue
   - L'interface met à jour le WCS interne et rafraîchit l'affichage

#### Tutoriel : Batch Astrométrie Zero-Aperture

1. **Charger le dossier d'images**
   - Même étape que ci-dessus (bouton "📁 Charger images")

2. **Choisir Zero-Aperture comme méthode**
   - Dans "Méthode", laissez "Zero-Aperture (extrapolation)"

3. **Lancer le batch**
   - Cliquez sur "🧭 Batch Astrométrie (Thread)"
   - Le processus :
     - Parcourt toutes les images du dossier
     - Pour **1 image sur 5**, lance la zero-aperture complète
     - Pour les 4 autres, utilise une version optimisée (skip extrapolation) pour gagner du temps

4. **Résultat batch**
   - Chaque fichier FITS se voit ajouter une solution WCS et des mots-clés `ASTR*`
   - Le log indique pour chaque image : méthode utilisée, RMS estimé, nombre d'étoiles, etc.
   - Vous obtenez une **série d'images astrométrées**, prêtes pour :
     - La mesure de positions d'astéroïde
     - L'export MPC/ADES

#### Comparaison Zero-Aperture vs Classique

**Zero-Aperture (extrapolation)** :
- ✅ **Avantages** : Meilleure précision, réduit les biais systémiques, statistiques détaillées
- ❌ **Inconvénients** : Plus lent
- **À utiliser si** : Vous cherchez la meilleure précision possible (soumissions MPC), le champ contient suffisamment d'étoiles Gaia

**Classique (FWHM)** :
- ✅ **Avantages** : Rapide, simple
- ❌ **Inconvénients** : RMS généralement moins bon, moins robuste
- **À utiliser si** : Vous voulez un résultat rapide, contrôle de champ, test rapide

### 5.2 Courbe de lumière, modèle harmonique, détrendage et ALCDEF

Dans l’onglet **Photométrie Astéroïdes**, sous-onglet **2. Photométrie**, la colonne droite contient l’**analyse de courbe de lumière** (fichier temps / flux issu de la photométrie batch ou d’un export). Vous pouvez estimer une **période de rotation**, appliquer un **détrendage** lent, ajuster un **modèle** (Fourier ou cosinus), puis **valider** le jeu pour l’onglet **3. Publication** (voir [5.3](#53-onglet-publication-jeu-validé-s-alcdef)).

**Organisation (onglet 2)** — ordre vertical approximatif :

1. **Fichier courbe de lumière** : chemin, **Parcourir**, **Charger** — uniquement une série **temps (JD-UTC) + flux** (éventuellement erreur), telle que produite par la photométrie NPOAP (`light_curve.txt` ou `results.csv`).  
1b. **Modèle DAMIT (surimpression phase)** : chargement optionnel d’un second fichier `light_curve.txt` **téléchargé depuis le site DAMIT** (phase de rotation 0–1 et flux relatif) pour le **tracer par-dessus** la courbe de phase de vos observations (voir [5.2a](#52a-surimpression-modèle-damit)).  
2. **Détrendage** : méthode (voir ci-dessous), champs optionnels **SG pts** (Savitzky–Golay) et **Wōtan j** (fenêtre en jours), bouton **Appliquer**.  
3. **Périodogramme (Lomb-Scargle)** et **Modèle** (côte à côte).  
4. **Graphiques** (courbe temps + phase), redimensionnés pour occuper l’espace disponible.  
5. **Valider pour publication…** : enregistre le fichier modèle et alimente l’onglet 3 (détail en [5.3](#53-onglet-publication-jeu-validé-s-alcdef)).

#### 5.2a Surimpression modèle DAMIT

La base [**DAMIT**](https://damit.cuni.cz/) (*Database of Asteroid Models from Inversion Techniques*, Charles University) propose pour certains astéroïdes un fichier **`light_curve.txt`** à deux colonnes : **phase de rotation** (entre 0 et 1) et **flux relatif** (comme sur le graphique web du site), parfois avec une figure SVG équivalente.

**Ce fichier n’est pas une courbe d’observations** : il ne contient pas de dates julienne. Si vous l’ouvrez avec **« 1. Fichier courbe de lumière »** puis **Charger**, NPOAP affiche une erreur et vous invite à utiliser le cadre **« 1b. Modèle DAMIT (surimpression sur la courbe de phase) »**.

**Procédure recommandée**

1. Chargez et traitez d’abord **votre courbe** (JD + flux) comme d’habitude ; affichez la **courbe de phase** en fixant **\(P\)**, **Lomb-Scargle** si besoin, puis **Calculer le modèle** (Fourier ou cosinus manuel) pour que la phase affichée corresponde à votre analyse.  
2. Dans **1b**, cliquez sur **Charger light_curve.txt (DAMIT)…** et sélectionnez le fichier issu du site.  
3. Cochez **Afficher** : une courbe **orange en tirets** (*DAMIT (modèle)*) apparaît sur le **graphique de phase** (panneau du bas), sur le **même axe de phase** que les points et le modèle rouge.  
4. **Déphasage [0–1]** : translation en phase du modèle DAMIT pour faire coïncider les bosses avec vos observations. La valeur utilisée est un décalage additif sur la phase (modulo 1).  
5. **Auto (aligner aux observations)** : recherche automatique du déphasage qui minimise (au sens des moindres carrés) l’écart aux points, avec un **facteur d’échelle** optimal à chaque test ; après alignement, **Échelle auto (médiane)** est activée. Vous pouvez encore affiner le déphasage à la main.  
6. **Échelle auto (médiane)** : multiplie le modèle DAMIT pour rapprocher globalement les flux du modèle de ceux des observations (échelle dérivée des médianes). Si vous décochez cette option, le **Facteur manuel** s’applique.  
7. **Effacer modèle** retire la surimpression et les données DAMIT chargées.

**Remarques**

- La surimpression sert à **comparer visuellement** la forme de phase issue du modèle publié sur DAMIT avec vos données ; ce n’est pas un ajustement publié dans le fichier `*_asteroid_model.txt` ni dans l’export ALCDEF.  
- Après modification du **déphasage** ou du **facteur** uniquement, recliquez sur **Calculer le modèle** si le graphique ne se met pas à jour (les cases **Afficher** et **Échelle auto** déclenchent un redessin lorsque c’est possible).

#### Détrendage (séries longues)

Le cadre **« Détrendage (sur données brutes du fichier) »** recalcule la courbe **affichée** à partir des **données brutes** du fichier à chaque **Appliquer** ou à chaque **Charger** : les tendances lentes (extinction, conditions, etc.) peuvent être atténuées avant Lomb-Scargle et le modèle.

Méthodes proposées (liste déroulante) :

- **Aucun** : pas de modification (hors recopie depuis le brut).  
- **Polynôme degré 1, 2 ou 3** : retrait d’un polynôme en temps (centré), en conservant un niveau proche du flux médian ; les erreurs ne sont pas modifiées.  
- **Savitzky–Golay** : lissage puis soustraction additive ; la **fenêtre en points** peut rester vide (choix automatique) ou être saisie (nombre **impair** de points).  
- **Wōtan biweight** / **Wōtan median** : estimateurs robustes de type `wotan.flatten` sur la durée couverte ; la fenêtre en **jours** peut rester vide (valeur dérivée de l’étendue temporelle). Ces lignes n’apparaissent **que si** le paquet Python **`wotan`** est installé (voir `requirements.txt` et les remerciements).

#### Invalidation de la « validation publication »

Sans recharger le fichier, la validation pour l’onglet 3 est **annulée** dès que vous :

- **Appliquez** un autre détrendage ;  
- lancez **Lancer Lomb-Scargle** ;  
- cliquez sur **Calculer le modèle** ou **Ajuster (inversion)**.

Tout **nouveau chargement** de courbe (ou rechargement automatique après photométrie batch) réinitialise aussi la validation. Il faut alors **Valider pour publication** à nouveau.

#### Rôle de la section « 3. Modèle » (harmonique ou manuel)

Après avoir chargé une courbe et, en général, fixé une **période \(P\)** (jours) — souvent issue du **Lomb–Scargle** —, l’étape **3** ajuste une **variation périodique** du flux dans le temps, typique d’une **courbe de rotation** d’astéroïde.

**Deux modes** (libellé de l’encadré) :

1. **Série de Fourier + offset (P fixe)** — *modèle harmonique*  
   Moindres carrés **pondérés** (les erreurs de flux sont prises en compte) avec **\(P\) imposé** (champ **P (j)**). Le flux est modélisé par  
   \[
   \mathrm{flux} \approx b_0 + \sum_{k=1}^{N} \bigl[ a_k \cos(k\,\omega (t - t_{\mathrm{ref}})) + b_k \sin(k\,\omega (t - t_{\mathrm{ref}})) \bigr]
   \]  
   avec \(\omega = 2\pi / P\), et éventuellement une **dérive linéaire** (tendance lente : extinction, conditions, etc.) si la case **« Dérive linéaire (tendance / extinction) »** est cochée.

   - **Harmoniques \(N\)** (valeur 1 à 8, défaut souvent 3) : nombre de termes de Fourier à la fréquence fondamentale \(1/P\) et à ses harmoniques \(k/P\).  
   - **\(N = 1\)** : une seule sinusoïde de période \(P\).  
   - **\(N \geq 2\)** : formes **non sinusoïdales** (deux bosses, etc.). Si la fenêtre temporelle ne couvre qu’une **faible fraction** de \(P\), une seule harmonique est souvent **mal contrainte** ; augmenter \(N\) peut améliorer l’ajustement lorsque la forme n’est pas purement sinusoïdale. L’interface peut signaler un modèle **fragile** lorsque \(\Delta t / P\) est petit.

2. **Modèle manuel** \(F_0 \cdot (1 + A\cos(2\pi(t-t_0)/P))\)  
   Si la case **« Série de Fourier + offset (P fixe) »** est **décochée**, le modèle est un **cosinus simple** : niveau **\(F_0\)**, amplitude relative **\(A\)**, phase via **\(t_0\)** (JD). Les paramètres sont ceux des champs **t0**, **A**, **F0** (éventuellement complétés après d’autres actions comme **Inversion**).

#### Comment choisir \(N\)

| Situation | Conseil |
|-----------|--------|
| Peu de points ou moins d’environ **un quart de période** couverte | Rester à **\(N = 1\) ou 2\)** ; éviter des \(N\) élevés qui sur-ajustent. |
| Courbe à **une bosse** (à peu près sinusoïdale) | **\(N = 1\)** ou **2** suffit souvent. |
| **Deux bosses** ou forme nettement **non sinusoïdale** | Essayer **\(N = 2\)**, puis **3** ; comparer **χ²**, les **résidus** et la **courbe de phase**. |
| **\(N\)** grand (4–8) | Réservé à **beaucoup de points** et une **bonne couverture en phase** (souvent plusieurs rotations ou une large fraction de \(P\)). |

**Méthode pratique** : partir de **\(N = 2\) ou 3**, augmenter si le modèle sous-ajuste visiblement, diminuer si la courbe modèle devient trop « bruitée » entre les points sans gain de χ².

#### Interpréter **t0** et **A**

- **Mode manuel** (\(F_0(1 + A\cos(\ldots))\)) :  
  - **\(t_0\)** (JD) : époque où le **cosinus vaut +1**, donc **maximum de flux** (pour \(A > 0\)).  
  - **\(A\)** : **amplitude relative sans dimension sur le flux**, pas directement une amplitude en **magnitudes**. En petite amplitude, \(A \approx \Delta\mathrm{flux}/F_0\) ; en ordre de grandeur si \(A \ll 1\), \(\Delta m \approx 1{,}09 \times A\) (soit \(2{,}5/\ln 10\)).

- **Mode Fourier** (\(N \geq 1\)) : le modèle réel est une **somme** de cosinus et sinus ; il n’existe plus un unique « instant de maximum » comme pour une seule sinusoïde. NPOAP calcule tout de même un **\(t_0\) équivalent** et un **\(A\) approximatif** à partir de la **première harmonique** uniquement, pour remplir les champs et rapprocher l’affichage du modèle cosinus. Si **\(N \geq 2\)** et que la forme est irrégulière, **\(t_0\)** et **\(A\)** affichés sont une **simplification** : la **courbe rouge** du graphe suit le **modèle Fourier complet**.

#### Graphiques (temps et phase)

- **Courbe en fonction du temps (JD)** : ordre réel des observations, lacunes, dérives résiduelles ; sur de longues durées la forme de rotation peut paraître **compressée** sur l’axe temps.  
- **Courbe en phase** : données repliées modulo **\(P\)** ; en général plus lisible pour la **morphologie** de la rotation sur une longue série.  
- **Fourier** : la phase est construite pour **aligner la première harmonique** ; les points et le modèle sont repliables en phase de façon cohérente avec \(\alpha\) (déphasage du terme \(k=1\)).  
- **Manuel** : la phase est basée sur \((t - t_0) \bmod P\) : **\(t_0\)** fixe l’**origine de phase** à 0 sur l’axe « phase ».  
- **Modèle DAMIT** : si vous en avez chargé un ([5.2a](#52a-surimpression-modèle-damit)), la courbe orange apparaît **uniquement** sur ce graphique de phase, avec le même abscisse **phase 0–1** que les observations.

### 5.3 Onglet Publication (jeu validé, S-ALCDEF)

L’onglet **3. Publication** ne repose **pas** sur la courbe « en cours d’édition » dans l’onglet 2 : il utilise uniquement le **dernier jeu validé** explicitement.

#### Bouton « Valider pour publication » (onglet 2)

Après détrendage, Lomb-Scargle et modèle souhaités, cliquez sur **Valider pour publication…** (sous les graphiques). NPOAP :

- enregistre un fichier texte **`{ID}_asteroid_model.txt`** ( **`ID`** dérivé du champ **ID cible** de la session, sinon du **N° MPC** saisi en Publication, sinon `UNKNOWN` ; caractères spéciaux adaptés au nom de fichier) ;  
- emplacement par défaut : sous-dossier **`results`** du dossier d’images chargé, ou **`results`** à côté du fichier de courbe si aucun dossier de session, ou dossier choisi par dialogue ;  
- le fichier contient des **lignes d’en-tête** (date UTC, nombre de points, méthode de détrendage, paramètres de modèle **P**, **t0**, **A**, **F0**, χ² si disponibles) puis les colonnes **JD_UTC**, **flux**, **flux_err** (deux colonnes seulement si pas d’erreurs) ;  
- **fige** ce jeu en mémoire pour l’export ALCDEF et le rapport ;  
- bascule l’affichage vers l’onglet **3. Publication**.

Tant que vous n’avez pas validé, le rapport indique qu’il faut valider depuis l’onglet 2 ; le bouton **Enregistrer S-ALCDEF…** refuse l’export sans jeu validé.

#### Contenu de l’onglet 3

- Bandeau d’**état** (jeu actif ou rappel de validation).  
- Métadonnées (N° MPC, observateurs, MagBand, etc.) et **rapport** de soumission généré à partir du **jeu validé** uniquement.  
- **Enregistrer S-ALCDEF…** : fichier **Simple ALCDEF** (JD, magnitude, erreur) ; un commentaire peut mentionner le chemin du fichier `*_asteroid_model.txt`.  
- Les paramètres **\(P\), \(A\), \(t_0\)** dans le rapport reflètent le **modèle figé au moment de la validation** ; pour un Fourier à **\(N>1\)**, la ligne « cosinus » reste une **approximation** (voir ci-dessous).

#### ALCDEF et rapport de soumission (rappel)

- Le fichier **S-ALCDEF** contient surtout les **données** : **JD**, **magnitude**, **erreur** (format Simple ALCDEF — voir [alcdef.org](https://alcdef.org/)). Les paramètres **\(P\), \(A\), \(t_0\)** ne font pas partie du bloc **DATA** standard.
- Le **texte d’aide** peut mentionner un ajustement **« modèle cosinus »** avec **\(P\), \(A\), \(t_0\)** : si vous avez utilisé un **Fourier à \(N>1\)**, cette ligne résume une **approximation cosinus** ; la description scientifique complète repose sur le **tracé du modèle**, le **χ²** et les **résidus**.
- Pour les publications ou LCDB, une **amplitude en magnitude** est souvent donnée par **écart min–max** des magnitudes (ou une définition dédiée) ; NPOAP peut aussi indiquer une **amplitude indicative en mag** dans le rapport — **distincte** de **\(A\)** (flux relatif du modèle cosinus).

#### Vérifications avant envoi

- Contrôler les métadonnées (bande magnétique ALCDEF, observateurs, filtre) sur le formulaire / guide [ALCDEF](https://alcdef.org/).  
- Valider le fichier avec **ALCDEFVerify** lorsque c’est requis par la procédure de dépôt.

---

## 6. Photométrie Transitoires

L'onglet **Photométrie Transitoires** est conçu pour l'analyse rapide d'événements transitoires (novae, supernovae, variables, etc.).

### Workflow

Le processus est similaire à celui de la photométrie d'astéroïdes :

1. **Chargement d'image** : Ouvrez une image FITS
2. **Astrométrie** : Résolution du champ (en ligne ou locale)
3. **Sélection de la cible** : Cliquez sur l'objet transitoire
4. **Sélection des comps** : Choisissez des étoiles de comparaison
5. **Analyse** : Calcul automatique de la magnitude

### Différences principales

- **Optimisé pour la rapidité** : Interface simplifiée pour un traitement rapide
- **Pas de format ADES** : Les résultats sont exportés en CSV standard
- **Analyse de séries** : Traitement par lots possible

### Tutoriel

1. Chargez les images calibrées du transitoire
2. Marquez la source transitoire (clic ou coordonnées)
3. Choisissez des étoiles de comparaison stables
4. Lancez la photométrie et suivez la magnitude vs temps
5. Exportez la courbe de lumière pour archivage ou publication

---

## 7. Analyse de Données

L'onglet **Analyse de Données** propose des outils d'analyse avancés pour les courbes de lumière. Il contient **4 sous-onglets** :

- **A. Détermination Période** : Recherche de période orbitale
- **B. Recherche & Analyse TTV** : Analyse des Transit Timing Variations
- **C. Analyse Système Multiple** : Comparaison de plusieurs planètes
- **D. Simulation N-body** : Simulation gravitationnelle de systèmes planétaires

---

### 7.1 Sous-onglet A : Détermination Période

#### But

Déterminer la **période orbitale** d'une exoplanète et extraire les **mid-times de transit** à partir d'un fichier de courbe de lumière.

#### Fonctionnalités principales

- **Chargement de courbe de lumière** : Format CSV/TXT (temps, flux, erreur)
- **Périodogrammes** :
  - **Lomb-Scargle** : Détection de périodicité dans les données
  - **BLS (Box Least Squares)** : Détection de transits
  - **Plavchan** : Méthode alternative
- **Extraction automatique des mid-times** : Génération d'un CSV de mid-times pour l'analyse TTV
- **Affinage du temps de transit central (Tc)** : modèle **trapézoïdal** local (SciPy) depuis la fenêtre **« Voir LC avec transits »**, époque par époque
- **Détrendage optionnel [Wōtan](https://github.com/hippke/wotan)** : dans le même viewer de fit, avec **masque transit** dérivé du Tc courant et de la largeur **T14** (voir ci-dessous)

#### Affinage des mid-times (fit trapèze et option Wotan)

1. Après **Voir LC avec transits (marqueurs, export mid-time.csv par période)** : la liste à gauche associe une **couleur** à chaque période retenue sur le périodogramme.
2. **Valider tous les transits** calcule des mid-times par minimum local dans une fenêtre autour de l’épéméride ; un **clic** sur la courbe fixe ou corrige un point pour la période sélectionnée.
3. **Affiner Tc (fit trapèze)** ouvre une fenêtre : choisissez l’**époque**, la **demi-fenêtre** (j), la **largeur transit T14** (j) et l’ingress ; lancez le **fit** (trapèze + baseline affine, pente optionnelle). Le graphique compare données et modèle ; **Appliquer Tc** enregistre le résultat (colonne **Uncertainty** et source `fit` à l’export).
4. **Détrendage Wotan** (case à cocher) : si le module Python **`wotan`** est installé (`pip install wotan`, aussi listé dans `requirements.txt`), le segment est d’abord passé par **`wotan.flatten`** en excluant de l’estimation de tendance les points où **|t − Tc| ≤ 1,15 × T14/2** (T14 = largeur « contact à contact » saisie). Le flux est ensuite ramené vers **~1** hors transit (médiane). Paramètres : **fenêtre Wotan (j)** (vide = automatique selon l’étendue du segment), **méthode** (ex. biweight, median, mean, huber, lowess). Si `wotan` est absent, la case est désactivée et un message indique l’installation.
5. **Conseil** : un flux déjà normalisé (~1 hors transit) facilite le fit ; Wotan est surtout utile lorsque la LC présente une **dérive lente** ou des systématiques sur l’intervalle considéré.

#### Tutoriel

1. **Fichiers source (courbes de lumière)**
   - **Ajouter fichiers…** : sélectionnez un ou plusieurs `.txt` / `.csv` (formats LcTools, Kepler/TESS commentés, CSV Time/Flux, etc.).
   - **Ajouter depuis un dossier…** : ajoute en une fois tous les `.txt` et `.csv` du dossier (hors `concatenated_lightcurve.csv` et `mid-time.csv`).
   - **Concaténer Lightcurves** : fusionne la liste dans `concatenated_lightcurve.csv` (dossier du **premier** fichier de la liste), puis charge la LC pour les périodogrammes.

2. **Lancer la recherche de période**
   - Choisissez la méthode (Lomb-Scargle, BLS, ou Plavchan)
   - Ajustez Min / Max P (jours) puis lancez le calcul
   - Le graphique affiche les pics de période

3. **Identifier la période orbitale**
   - Examinez le périodogramme pour identifier le pic principal
   - La période proposée est affichée dans l'interface

4. **Extraire et affiner les mid-times**
   - Ouvrez **Voir LC avec transits**, puis **Valider tous les transits** et/ou des **clics** sur la courbe pour fixer les mid-times.
   - (Optionnel) Utilisez **Affiner Tc (fit trapèze)** et, si besoin, le **détrendage Wotan** pour raffiner chaque époque avant export.
   - **Exporter mid-time.csv par période** : le CSV (Epoch, Mid_time, O-C_minutes, Uncertainty, …) sert d’entrée au sous-onglet **B (TTV)** pour les O-C et l’analyse TTV.

---

### 7.2 Sous-onglet B : Recherche & Analyse TTV

#### But

Analyser les **Transit Timing Variations (TTV)** à partir des mid-times de transit. Les TTV peuvent révéler la présence de planètes supplémentaires dans le système.

#### Fonctionnalités principales

- **Chargement des mid-times** : CSV de mid-times ou tableau O-C
- **Calcul des O-C** : Calcul des résidus O-C (en jours/minutes)
- **Affichage de la courbe TTV** : Visualisation interactive avec TTVViewer
- **Fit MCMC sinusoïdal** : Ajustement d'un ou plusieurs sinusoïdes aux données TTV
- **Paramètres configurables** :
  - **Nombre de fréquences** : Nombre de sinusoïdes à ajuster
  - **nwalkers** : Nombre de marcheurs MCMC (défaut: 32)
  - **nsteps** : Nombre d'itérations MCMC (défaut: 5000)
  - **burn-in** : Fraction de burn-in (défaut: 0.25)
- **Visualisation avancée** :
  - Courbe TTV principale + panneau de résidus
  - Échelles dynamiques basées sur l'amplitude
  - RMS des résidus affiché
- **Génération de rapport TTV** : Rapport complet avec :
  - Amplitude TTV, Période TTV (P_ttv en époques & jours)
  - BIC nul vs BIC TTV, interprétation
  - Vérification physique (ratio P_ttv/P_orb)
  - Prédiction de résonances (1:2, 2:3, 3:2, 2:1, 3:1)

#### Tutoriel

1. **Charger les mid-times**
   - Cliquez sur "Charger CSV" et sélectionnez le fichier de mid-times généré en A
   - Vérifiez la **P_orb détectée** ou saisissez-la si nécessaire

2. **Afficher la courbe TTV**
   - Cliquez sur "Afficher O-C" pour voir la courbe TTV

3. **Configurer les paramètres de fit**
   - Cliquez sur le bouton **"⚙️ Param. Fit"**
   - Ajustez :
     - **Nombre de fréquences** : 1 pour un signal simple, 2-3 pour des signaux complexes
     - **nwalkers** : 32-64 (plus = meilleur échantillonnage mais plus lent)
     - **nsteps** : 5000-10000 (plus = meilleure convergence)
     - **burn-in** : 0.25 (fraction à ignorer au début)

4. **Lancer le fit TTV**
   - Cliquez sur "Lancer Fit TTV (MCMC)"
   - Le processus peut prendre quelques minutes selon les paramètres
   - Inspectez :
     - La sinusoïde ajustée sur la courbe TTV
     - Le panneau de résidus, RMS, etc.

5. **Générer le rapport TTV**
   - Cliquez sur "Générer Rapport TTV"
   - Le rapport contient :
     - Les paramètres du fit (amplitude, période, phase)
     - Les statistiques (BIC, RMS)
     - Les prédictions de résonances (périodes possibles des planètes perturbatrices)
   - Sauvegardez le fichier texte
   - Ce rapport servira pour l'onglet C (système multiple) et le transfert vers D

---

### 7.3 Sous-onglet C : Analyse Système Multiple

#### But

Comparer **plusieurs rapports TTV** (plusieurs planètes) pour analyser un système multi-planétaire.

#### Fonctionnalités principales

- **Chargement de plusieurs rapports TTV** : Fichiers `.txt` générés en B
- **Affichage dans une liste** : Chaque rapport affiche **P_orb** et **P_ttv** de la planète
- **Analyse des paires de planètes** :
  - **Ratio des périodes** (P₂/P₁)
  - Différences de phase (Δφ)
  - Indices de résonance potentielle
- **Génération d'un rapport système multiple** : Rapport comparatif (texte)
- **Bouton "Transférer vers Simulation N-body"** :
  - Tri des rapports par P_orb
  - Remplissage automatique de l'onglet D avec les planètes (P_orb, masses estimées)
  - Détection des périodes orbitales dupliquées avec avertissement

#### Tutoriel

1. **Ajouter des rapports**
   - Cliquez sur "➕ Ajouter Rapport" et sélectionnez plusieurs fichiers de rapport TTV
   - Chaque fichier correspond à une planète différente
   - Vérifiez dans la liste que **P_orb** et **P_ttv** sont cohérents pour chaque planète

2. **Analyser le système**
   - Cliquez sur "🚀 Analyser Système" pour générer un rapport comparatif
   - Le rapport contient :
     - Les ratios de périodes entre paires de planètes
     - Les différences de phase
     - Les indices de résonance

3. **Transférer vers Simulation N-body**
   - Cliquez sur "🔄 Transférer vers Simulation N-body"
   - Entrez la masse de l'étoile
   - Pour chaque planète :
     - Vérifiez la période orbitale (P_orb)
     - Confirmez/ajustez la masse estimée (Mjup)
   - **Attention** : Si deux planètes ont la même P_orb, une boîte de dialogue vous avertit
   - Passez à l'onglet D pour lancer la simulation

---

### 7.4 Sous-onglet D : Simulation N-body

#### But

Simuler numériquement les interactions gravitationnelles entre les planètes d'un système et calculer les TTV simulés.

#### Prérequis

- **rebound** : Bibliothèque Python pour simulations N-body (`pip install rebound`)
- **ultranest** (optionnel) : Pour le fitting N-body (`pip install ultranest`)

#### Fonctionnalités principales

- **Saisie de la masse de l'étoile** : En masses solaires (Msun)
- **Gestion d'une liste de planètes** :
  - **Ajout manuel** : Entrez m, P, e, i, ω, etc.
  - **Remplissage automatique** : Depuis les rapports TTV (onglet C)
  - **Édition/Suppression** : Modifiez ou supprimez des planètes
- **Lancement de simulations N-body** :
  - Intégration temporelle du système
  - Extraction de **temps de transit** simulés
  - Calcul des **TTV simulés**
- **(Optionnel) Fitting N-body** : Ajustement des paramètres pour reproduire les TTV observées
- **Visualisations** : Évolution orbitale, TTV simulés

#### Tutoriel

1. **Vérifier les dépendances**
   - Assurez-vous que `rebound` est installé
   - Si ce n'est pas le cas, l'onglet affiche un message avec la commande `pip install rebound`

2. **Remplir les planètes**
   - **Option 1** : Utilisez le bouton "Transférer vers Simulation N-body" en C
   - **Option 2** : Ajoutez manuellement des planètes (m, P, e, i, ω)

3. **Vérifier la liste des planètes**
   - Vérifiez que les **masses** sont plausibles
   - Vérifiez que les **périodes** sont différentes (NPOAP avertit en cas de P_orb dupliquée)

4. **Configurer la simulation**
   - Entrez la masse de l'étoile
   - (Si disponible) Configurez les paramètres de simulation (durée, pas de temps, etc.)

5. **Lancer la simulation**
   - Cliquez sur "Lancer Simulation N-body"
   - Le processus peut prendre du temps selon le nombre de planètes et la durée
   - Examinez :
     - Les trajectoires orbitales
     - Les TTV simulés

6. **Comparer avec les observations**
   - Comparez les TTV simulés avec les TTV observées (onglet B)
   - Ajustez les paramètres si nécessaire
   - Utilisez le fitting N-body (si disponible) pour optimiser automatiquement

---

## 8. Étoiles Binaires

L'onglet **Étoiles Binaires** permet de modéliser et d'analyser des systèmes d'étoiles binaires à éclipses en utilisant PHOEBE2.

### Prérequis

PHOEBE2 doit être installé pour utiliser cet onglet. Si ce n'est pas le cas, l'application affichera un message avec les instructions d'installation.

### Création d'un système

1. **Type de système** :
   - **Binaire à éclipses** : Système binaire classique
   - **Système de contact** : Étoiles en contact

2. **Créer un Bundle** : Cliquez sur "Créer un nouveau Bundle" pour initialiser un système

### Chargement de données observées

1. **Format CSV** : Préparez un fichier CSV avec les colonnes :
   - `time` : Temps (JD ou jours depuis une référence)
   - `flux` : Flux observé (normalisé ou absolu)
   - `error` : Erreur sur le flux (optionnel)

2. **Chargement** : Utilisez "Parcourir" pour sélectionner votre fichier, puis "Charger les données"

### Paramètres du modèle

Dans la section **Paramètres du modèle**, vous pouvez ajuster :

- **Période** : Période orbitale en jours
- **t0** : Époque de référence (JD)
- **Inclinaison** : Angle d'inclinaison du système en degrés

### Calcul du modèle

1. **Calculer le modèle** : Lance le calcul avec les paramètres actuels
2. **Visualisation** : La courbe de lumière calculée s'affiche avec les observations (si chargées)

### Ajustement des paramètres

Le bouton **Ajuster les paramètres** lance un processus d'optimisation pour trouver les paramètres qui correspondent le mieux aux observations.

**Note** : Cette fonctionnalité peut prendre du temps selon la complexité du système.

### Visualisation 3D

Le bouton **🎬 Visualisation 3D** ouvre une fenêtre interactive permettant de :

1. **Visualiser l'orbite** :
   - **Vue de face** : Vue dans le plan orbital (X-Y)
   - **Vue de côté** : Vue perpendiculaire montrant l'inclinaison (X-Z)
   - **Vue de dessus** : Vue depuis le pôle (Y-Z)

2. **Animation** :
   - **Play/Pause** : Contrôle de l'animation
   - **Reset** : Retour au début
   - **Vitesse** : Réglage de la vitesse d'animation (0.1x à 5.0x)

3. **Informations** : Affichage en temps réel des paramètres et de la phase orbitale

4. **Trajectoires** : Visualisation des orbites complètes des deux étoiles

### Sauvegarde et chargement

- **Sauvegarder le Bundle** : Enregistre votre modèle PHOEBE2 pour réutilisation ultérieure
- **Charger un Bundle** : Recharge un modèle précédemment sauvegardé

### Workflow recommandé

1. Créer un Bundle (système binaire ou contact)
2. Charger vos données observées (CSV)
3. Ajuster manuellement les paramètres de base (période, t0, inclinaison)
4. Calculer un premier modèle
5. Utiliser "Ajuster les paramètres" pour optimiser
6. Visualiser en 3D pour comprendre la géométrie du système
7. Sauvegarder le Bundle final

---

## 9. Easy Lucky Imaging

L'onglet **Easy Lucky Imaging** permet de traiter des images d'étoiles doubles/binaires avec des techniques de réduction d'images (méthodes REDUC) et de mesurer précisément la séparation et l'angle de position entre deux étoiles.

### 9.1 Mesure de séparation

La fonction **Mesure de séparation** permet de mesurer précisément la séparation angulaire (ρ) et l'angle de position (θ) entre deux étoiles sur une image astrométriée.

#### Prérequis

- **Image astrométriée** : L'image doit avoir été astrométriée (plate-solving) pour déterminer le Nord céleste et calculer les angles correctement
- **WCS valide** : Les informations WCS (World Coordinate System) doivent être présentes dans l'en-tête FITS

#### Fonctionnalités principales

- **Détermination automatique du Nord et de l'Est célestes** :
  - Utilise la fonction `get_image_orientation()` qui calcule précisément l'orientation réelle du Nord et de l'Est sur le capteur via le WCS
  - Utilise `proj_plane_pixel_scales()` pour obtenir les échelles pixel exactes
  - Projette des points célestes (Dec+ pour Nord, RA+ pour Est) sur le plan image
  - Calcule la rotation du capteur et la parité (Standard ou Inversée/Miroir/Sud)

- **Sélection interactive des étoiles** :
  - **Étoile 1 (Primaire)** : Cliquez sur la première étoile, le centroïde est affiné automatiquement
  - **Étoile 2 (Secondaire)** : Cliquez sur la seconde étoile, le centroïde est affiné en excluant la position de l'étoile 1

- **Affichage graphique** :
  - **Flèche Nord (N)** : Flèche blanche pointant vers le Nord céleste réel
  - **Flèche Est (E)** : Flèche blanche pointant vers l'Est céleste réel
  - **Ligne de séparation** : Ligne blanche reliant les deux étoiles
  - **Arc theta (θ)** : Arc pointillé blanc partant du Nord, passant par l'Est, et allant jusqu'à la séparation
    - Le sens de l'arc est **constant** : toujours du Nord vers l'Est, indépendamment de la valeur de θ
    - L'arc utilise directement les angles WCS calculés pour déterminer le sens

- **Mesures affichées** :
  - **Séparation (ρ)** : Séparation angulaire en secondes d'arc (si pixel scale disponible) ou en pixels
  - **Angle de position (θ)** : Angle de position depuis le Nord céleste vers la séparation, mesuré positivement vers l'Est (0° = Nord, 90° = Est, 180° = Sud, 270° = Ouest)

#### Tutoriel

1. **Préparer l'image** :
   - Assurez-vous que l'image a été astrométriée (onglet "Réduction de Données" → Astrométrie)
   - Chargez l'image dans l'onglet "Easy Lucky Imaging"

2. **Ouvrir la fenêtre de mesure** :
   - Cliquez sur le bouton **"Afficher Image pour Mesure"**
   - L'application vérifie automatiquement la présence du WCS
   - Si le WCS est absent, un message d'erreur s'affiche

3. **Sélectionner l'étoile 1 (Primaire)** :
   - Cliquez sur la première étoile (généralement la plus brillante)
   - Le centroïde est automatiquement affiné
   - Un marqueur blanc "+" et le label "Star 1" apparaissent

4. **Sélectionner l'étoile 2 (Secondaire)** :
   - Cliquez sur la seconde étoile
   - Le centroïde est affiné en excluant la position de l'étoile 1
   - Un marqueur blanc "+" et le label "Star 2" apparaissent

5. **Lire les mesures** :
   - **Séparation (ρ)** : Affichée au milieu de la ligne de séparation, en secondes d'arc (") ou en pixels
   - **Angle de position (θ)** : Affiché le long de l'arc pointillé, en degrés depuis le Nord vers l'Est
   - Les valeurs sont aussi affichées dans les champs textuels sous l'image

6. **Vérifier l'affichage** :
   - Vérifiez que la flèche Nord (N) pointe bien vers le Nord céleste
   - Vérifiez que la flèche Est (E) pointe bien vers l'Est céleste (à 90° du Nord dans le sens horaire céleste)
   - Vérifiez que l'arc θ part bien du Nord, passe par l'Est, et va jusqu'à la séparation

#### Ajustement du contraste

- Utilisez les sliders **Minimum** et **Maximum** pour ajuster le contraste de l'image
- Le bouton **"Réinitialiser (Médiane)"** remet les valeurs par défaut
- Les valeurs sont centrées sur la médiane de l'image pour un meilleur contraste

#### Réinitialisation

- Le bouton **"🔄 Réinitialiser"** efface toutes les sélections et permet de recommencer

#### Notes techniques

- **Méthode de calcul du Nord/Est** : La fonction `get_image_orientation()` :
  - Utilise `proj_plane_pixel_scales()` pour obtenir les échelles pixel en degrés
  - Projette un point vers le Nord (Dec+) et un point vers l'Est (RA+) depuis le centre de l'image
  - Convertit ces coordonnées célestes en pixels via `world_to_pixel()`
  - Calcule les angles dans le système matplotlib (0° = droite, sens anti-horaire)

- **Affichage de l'arc theta** :
  - L'arc utilise les angles `north_celestial_angle` et `east_celestial_angle` calculés via WCS
  - Le sens de l'arc est déterminé par la direction de l'Est depuis le Nord
  - L'arc part toujours du Nord, passe par l'Est, et continue jusqu'à la séparation
  - Ce sens est **constant** et indépendant de la valeur de θ

- **Précision** :
  - La précision dépend de la qualité de l'astrométrie (WCS)
  - La précision dépend aussi de l'affinage du centroïde (méthode 2D Gaussienne)
  - Pour les meilleures mesures, utilisez des images avec un bon SNR et un bon WCS

---

## 10. Analyse des outils et fichiers compatibles LcTools (TESS)

Cette section documente **uniquement** le travail d’analyse et d’interopérabilité avec les **exports de type LcTools** (courbes TESS normalisées en texte) : où ils sont produits dans NPOAP, quel format ils suivent, et comment le code les **reconnaît** ensuite pour concaténation ou analyse.

### 10.1 Où produire un export « type LcTools » dans l’interface

- Onglet principal **📚 Catalogues** → sous-onglet **🔭 Exoplanètes** → sous-onglet interne **📊 Courbe TESS : FITS → LcTools**.
- **Étape 1 — FITS édité** : à partir d’un fichier produit TESS `*_lc.fits`, masquage des cadences en transit (NaN sur le flux en transit) et enregistrement d’une copie FITS ; les priors de transit (T₀, P, durée) peuvent venir d’**exoplanet.eu** et/ou de la **NASA Exoplanet Archive**, ou d’une saisie manuelle ; fenêtres de transit éditables à la souris si l’option est activée (voir module `core/lc_transit_pick.py`).
- **Étape 2 — Normalisation LcTools → .txt** : à partir du FITS de l’étape 1, appel du script `scripts/tess_lc_fits_to_txt.py` en mode normalisé : temps en **BJD-TDB**, flux **divisé par la moyenne des points hors transit** (les cadences en transit sont exclues du calcul de baseline), avec options (qualité, chemins de sortie, etc.). Le journal de la zone de texte reprend les messages `[FITS édité]` et `[LcTools]`.
- **FITS « édités » et transits** : l’étape 1 met le flux à **NaN** pendant le transit ; sans valeurs d’origine, l’export ne gardait que les points finis. Le script retrouve le FITS brut si le nom se termine par **`_edited.fits`** (ex. `*_lc_edited.fits` SPOC, `*_fast-lc_edited.fits` HLSP TESS-fast) : **même nom sans `_edited`**, même dossier. Sinon : CLI `--flux-source-fits`.

En ligne de commande, le même comportement est décrit dans l’en-tête de `scripts/tess_lc_fits_to_txt.py` (arguments `--normalized-lctools`, `--norm-baseline`, `--pick-transit-windows`, `--planet-name`, `--flux-source-fits`, etc.).

### 10.2 Fusion des épémérides (catalogues → script TESS/LcTools)

Le module `core/transit_catalog_merge.py` sert au **flux Catalogues → TESS → LcTools** : il fusionne les informations de transit issues **d’exoplanet.eu** (TAP) et de la **NASA Exoplanet Archive** (via astroquery) pour fournir P, T₀ et la durée de transit (T₁₄) au script, avec des règles de priorité documentées dans le script (CLI > EU > NASA selon les paramètres).

### 10.3 Analyse automatique du format fichier (reconnaissance « style LcTools »)

Le module `core/lightcurve_tools.py` implémente la **lecture et l’homogénéisation** des courbes texte/CSV pour les usages aval (concaténation, périodogrammes, TTV, etc.).

1. **Fichiers au style export LcTools** (produits par `tess_lc_fits_to_txt.py --normalized-lctools`) :
   - Première ligne non vide : commentaire commençant par `#` contenant le mot **time** (ex. `#Time (BJD-TDB),...`) suivi des noms de colonnes séparés par des virgules.
   - Les données commencent à la ligne suivante ; le lecteur reconstruit les noms de colonnes normalisés (espaces → `_`, casse unifiée).

2. **Identification des colonnes temps et flux** : après normalisation des en-têtes, NPOAP cherche des **synonymes** :
   - **Temps** : sous-chaînes parmi `TIME`, `BJD`, `JD`, `BJD_TDB`, `BTJD`, `HJD`, `DATE`, `MID_TIME`, etc.
   - **Flux** : `DETRENDED_FLUX`, `PDCSAP_FLUX`, `CORR_FLUX`, `FLUX`, `SAP_FLUX`, `RAW_FLUX`, `NORMALIZED` (ce dernier correspond typiquement au flux normalisé issu de l’export LcTools).

3. **Si l’entête LcTools n’est pas détectée** : le code tente un CSV classique (ligne d’en-tête `Time`, `Flux`, …) puis, à défaut, l’extraction des lignes commentées **Kepler/TESS** du type `# Column N: …` pour mapper les indices de colonnes vers des noms reconnus. Sinon le fichier est ignoré avec un message d’avertissement explicite.

4. **Concaténation** (`concatenate_lightcurves`) : pour chaque fichier reconnu, les paires (temps, flux) sont extraites ; une **normalisation supplémentaire par la médiane du flux** est appliquée par fichier avant concaténation, afin d’harmoniser l’échelle entre segments.

### 10.4 Résumé pratique

| Élément | Rôle |
|--------|------|
| GUI Catalogues → Exoplanètes → TESS → LcTools | Enchaîne édition FITS + export `.txt` normalisé type LcTools |
| `scripts/tess_lc_fits_to_txt.py` | Conversion FITS TESS → texte ; mode LcTools : BJD-TDB + baseline hors transit |
| `core/transit_catalog_merge.py` | Fusion EU + NASA pour P, T₀, T₁₄ |
| `core/lc_transit_pick.py` | Édition interactive des fenêtres « en transit » |
| `core/lightcurve_tools.py` | **Analyse** des fichiers : détection entête LcTools, synonymes colonnes, fallback Kepler/TESS |

Pour les autres fonctions de l’onglet **Catalogues** (Vizier, Gaia MPC, binaires, MAST, etc.), se référer aux libellés et info-bulles dans l’interface ; elles ne sont pas détaillées dans cette section.

---

## 11. Spectroscopie

L'onglet **Spectroscopie** permet de charger, visualiser et analyser des spectres astronomiques. Il comprend deux types d'analyse :

- **Analyse de spectres d'étoiles** : Mesure de raies, normalisation de continuum, analyse spectrale standard
- **Analyse de galaxies avec Prospector** : Inférence de propriétés stellaires à partir de SED (Spectral Energy Distribution)

### Chargement de spectres

1. **Format supporté** : FITS ou texte (CSV)
2. **Chargement** : Utilisez le bouton "Charger spectre" pour sélectionner votre fichier
3. **Détection automatique** : L'application détecte automatiquement le format (FITS structuré, FITS standard, ASCII)

### Visualisation

- **Graphique flux vs longueur d'onde** : Affichage interactif du spectre
- **Zoom** : Utilisez les outils de zoom pour examiner des régions spécifiques
- **Navigation** : Parcourez le spectre avec les contrôles

### Outils d'analyse de spectres d'étoiles

#### Normalisation de continuum

1. **Bouton "Normaliser Continuum"** : Calcule et applique une normalisation du continuum
2. **Méthode** : Ajustement polynomial du continuum (gère automatiquement les valeurs non-finies)
3. **Résultat** : Le spectre normalisé est affiché avec le continuum ajusté

#### Analyse de raies spectrales

Une fois un spectre chargé et normalisé, vous pouvez :

- **Sélectionner une région** : Cliquez-glisser sur le graphique pour définir une région spectrale
- **Mesurer des quantités** :
  - **Équivalent de largeur** (EW) : Largeur équivalente d'une raie
  - **Flux de raie** : Flux intégré dans la raie
  - **Centroïde** : Position centrale de la raie
  - **FWHM** : Largeur à mi-hauteur de la raie

### Analyse de galaxies avec Prospector (Optionnel)

Prospector permet d'inférer les propriétés stellaires des galaxies à partir de leurs spectres et/ou de leurs données photométriques (SED).

#### Prérequis

- **Prospector installé** : Utilisez `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat` pour l'installation
- **FSPS** (optionnel) : Recommandé pour les fonctionnalités avancées

#### Statut de Prospector

L'onglet affiche automatiquement le statut de Prospector :
- **"Installed"** : Prospector est disponible et prêt à l'emploi
- **"Non installé"** : Cliquez sur "Installer Prospector" pour l'installer automatiquement

#### Fonctionnalités Prospector

1. **🌌 Inférer Propriétés Stellaires (SED)** :
   - Chargez un spectre de galaxie (FITS ou ASCII)
   - Prospector ajuste un modèle de population stellaire simple (SSP)
   - Infère : âge, métallicité, extinction de la poussière, redshift
   - Utilise l'inférence bayésienne (MCMC avec dynesty ou emcee)

2. **📊 Créer SED depuis Photométrie** :
   - Entrez des données photométriques (magnitudes ou flux dans différents filtres)
   - Créez une SED à partir de ces données
   - Combine avec des données spectroscopiques si disponibles

#### Tutoriel : Analyse avec Prospector

1. **Installer Prospector** (si non installé) :
   - Cliquez sur "Installer Prospector" dans la section "Analyse de Galaxies (Prospector)"
   - Attendez la fin de l'installation (peut prendre plusieurs minutes)
   - Redémarrez l'application

2. **Charger un spectre de galaxie** :
   - Cliquez sur "Charger spectre" et sélectionnez un fichier FITS ou ASCII
   - Le spectre s'affiche dans le graphique principal

3. **Normaliser le continuum** (optionnel mais recommandé) :
   - Cliquez sur "Normaliser Continuum"
   - Vérifiez que le continuum est correctement ajusté

4. **Inférer les propriétés stellaires** :
   - Cliquez sur "🌌 Inférer Propriétés Stellaires (SED)"
   - Configurez les paramètres du modèle (si des options sont disponibles)
   - Lancez l'inférence (peut prendre plusieurs minutes selon les paramètres)
   - Examinez les résultats : âge, métallicité, extinction, redshift

5. **Créer SED depuis photométrie** (alternative) :
   - Cliquez sur "📊 Créer SED depuis Photométrie"
   - Entrez les données photométriques (magnitudes ou flux par filtre)
   - La SED est créée et peut être utilisée pour l'inférence

### Export et sauvegarde

- Les spectres normalisés peuvent être exportés
- Les résultats d'analyse (EW, flux, etc.) sont affichés dans l'interface
- Les résultats Prospector sont affichés sous forme de résumé textuel

### Tutoriel rapide

**Pour les spectres d'étoiles** :
1. Chargez un spectre dans l'onglet **Spectroscopie**
2. Normalisez le continuum si nécessaire
3. Sélectionnez une région spectrale pour mesurer des raies
4. Examinez les quantités mesurées (EW, flux, FWHM, etc.)

**Pour les galaxies (Prospector)** :
1. Installez Prospector si nécessaire
2. Chargez un spectre de galaxie
3. Cliquez sur "🌌 Inférer Propriétés Stellaires (SED)"
4. Attendez les résultats de l'inférence bayésienne
5. Examinez les propriétés stellaires inférées

---

## 12. Observation de la nuit

L’onglet **🌙 Observation de la nuit** permet de préparer une soirée d’observation : téléchargement des catalogues locaux (MPC, AAVSO), calcul des éphémérides pour la date choisie, liste des cibles et graphique d’altitude sur la nuit astronomique, export JSON vers **NINA**.

### Catalogues et dossiers

- **MPC** : `NEA.txt` (astéroïdes géocroiseurs) et `AllCometEls.txt` (comètes), dans le dossier MPC indiqué (par défaut sous `~/.npoap/catalogues`).
- **AAVSO** : `index.csv` (binaires à éclipses uniquement) dans le dossier AAVSO/NASA indiqué.
- Le bouton **Téléchargement NASA Exoplanet Archive** (export CSV local `nasa_exoplanet_transits.csv`) a été **retiré de l’interface** ; les exoplanètes en transit sont interroguées **directement au moment du calcul** via le service TAP NASA.

### Exoplanètes en transit : plusieurs sources et fusion

Au calcul des éphémérides, les exoplanètes proviennent de sources combinées :

- **ExoClock** (module Python `exoclock`, si installé) : **priorité principale** pour la fiche retenue en cas de doublon (notamment les magnitudes).
- **NASA Exoplanet Archive** (TAP, table `pscomppars`, planètes en transit) : source complémentaire.
- **ETD / VarAstro** : intégration prévue ; l’API catalogue (`/api/Search/Exoplanets`) exige une **authentification**. Sans session valide, la source ETD est **désactivée automatiquement** (pas de requêtes inutiles).

Les doublons sont fusionnés avec une **clé canonique** sur le nom (insensible aux espaces, tirets, casse ; normalisation du type `WASP-12 b` → `WASP-12b`). Si la fiche prioritaire manque un champ (par ex. magnitude, période, milieu de transit, durée, profondeur), il est **complété** à partir de l’autre source lorsque possible.

À la fin du calcul, une boîte de dialogue peut indiquer le détail par source (nombre de cibles avant/après dédoublonnage).

### Affichage : liste et graphique

- **Liste « Objets observables »** : pour les exoplanètes, le nom est précédé de **`[ExoClock]`** ou **`[NASA]`** selon la source retenue pour la ligne ; la ligne est colorée (**rouge** / **bleu**) pour distinguer visuellement les sources.
- **Colonne « Mag vis »** : libellé utilisé pour les magnitudes affichées (exoplanètes : magnitude hôte/visuelle selon la source ; **ce n’est pas** la magnitude absolue **H** du fichier NEA pour les astéroïdes — voir ci-dessous).
- **Graphique d’altitude** : les courbes d’altitude et les schémas de transit des objets **cochés** utilisent les mêmes couleurs par source (**ExoClock en rouge**, **NASA en bleu**) pour rester cohérents avec la liste.

Sélectionnez les cibles avec la case **☑** dans la première colonne ; seules les cibles cochées sont surlignées sur le graphique (dont les schémas de transit pour les exoplanètes).

### Astéroïdes : filtre H et colonne « Mag vis »

Pour les astéroïdes issus de `NEA.txt`, la magnitude absolue **H** (champs fixes MPC) sert au **filtre** initial (objets avec H trop élevée exclus avant l’appel d’éphémérides). La valeur affichée dans **Mag vis** est extraite de la **réponse du service d’éphémérides MPC (MPES)** pour la date et le lieu configurés : il s’agit en principe d’une **magnitude apparente estimée** (souvent assimilée à une magnitude visuelle à l’éphéméride), et non de **H**.

### Export NINA

Les objets cochés peuvent être exportés en JSON pour **NINA** ; les coordonnées ICRS saisies manuellement peuvent également être exportées si l’option correspondante est activée.

---

## 13. Conseils généraux

### Organisation des données

- Gardez vos images organisées dans des dossiers clairs
- Utilisez des noms de fichiers cohérents
- Sauvegardez régulièrement vos configurations

### Qualité des données

- Vérifiez toujours la qualité des images avant analyse
- Utilisez un nombre suffisant d'étoiles de comparaison (au moins 3-5)
- Vérifiez que les étoiles de comparaison ne sont pas variables

### Export et sauvegarde

- Exportez régulièrement vos résultats
- Conservez les fichiers sources (images FITS)
- Documentez vos observations

### Performance

- Pour de grandes séries d'images, le traitement peut prendre du temps
- Utilisez l'astrométrie locale pour de meilleures performances si vous avez beaucoup d'images
- L'ajustement de modèles PHOEBE2 peut être long pour les systèmes complexes
- Pour l'astrométrie zero-aperture en batch, le processus optimise automatiquement (1 image sur 5 avec extrapolation complète)

### Astrométrie Zero-Aperture

- **Quand l'utiliser** : Pour des soumissions MPC sérieuses, lorsque vous avez besoin de la meilleure précision possible
- **Quand utiliser la méthode classique** : Pour des tests rapides, des contrôles de champ, ou lorsque la vitesse est prioritaire
- **Vérifiez les statistiques** : Les mots-clés `ASTR*` dans le header FITS vous donnent toutes les informations sur la qualité de l'astrométrie

### Analyse TTV

- **Nombre de fréquences** : Commencez avec 1, augmentez si le signal est complexe
- **Paramètres MCMC** : Augmentez `nwalkers` et `nsteps` pour une meilleure convergence si nécessaire
- **Interprétation des rapports** : Les prédictions de résonances sont des hypothèses, pas des confirmations

### Simulation N-body

- **Vérifiez les périodes** : Assurez-vous que les planètes ont des périodes différentes
- **Masses réalistes** : Utilisez des masses plausibles (typiquement 0.1-10 Mjup pour les exoplanètes)
- **Durée de simulation** : Simulez sur plusieurs périodes orbitales pour voir les effets à long terme

---

## Support

Pour toute question ou problème, consultez la documentation technique ou contactez l'équipe de développement.

---

**NPOAP - Manuel Utilisateur v1.0** (mise à jour : onglet Observation de la nuit, sources exoplanètes, libellé Mag vis)

*Ce manuel décrit l'utilisation de base de NPOAP. Pour les procédures d'installation, consultez le manuel d'installation séparé.*
