# Manuel d'Installation NPOAP

**NPOAP - Nouvelle Plateforme d'Observation et d'Analyse Photométrique**

Version 1.0 (révisions procédure : vérification `test_installation.py`, CuPy/`importlib.metadata`, ligne KBMOD dans `requirements.txt`)

**Responsable HOPS-modified** : J.P Vignes  
**Contact** : jeanpascal.vignes@gmail.com

---

## Table des matières

1. [Prérequis système](#1-prérequis-système)
2. [Installation de Python et Conda](#2-installation-de-python-et-conda)
3. [Création de l'environnement Conda](#3-création-de-lenvironnement-conda)
4. [Installation des dépendances Python](#4-installation-des-dépendances-python)
5. [Installation des dépendances optionnelles](#5-installation-des-dépendances-optionnelles)
6. [Installation d'Astrometry.net (optionnel)](#6-installation-dastrometrynet-optionnel)
7. [KBMOD – Détection Synthetic Tracking (optionnel, via WSL)](#7-kbmod--détection-synthetic-tracking-optionnel-via-wsl)
8. [Configuration initiale](#8-configuration-initiale)
9. [Vérification de l'installation](#9-vérification-de-linstallation)
10. [Lancement de l'application](#10-lancement-de-lapplication)
11. [Dépannage](#11-dépannage)
12. [Remerciements (Acknowledgments)](#12-remerciements-acknowledgments)

---

## 1. Prérequis système

### Système d'exploitation

NPOAP est compatible avec :
- **Windows 10/11** (testé et recommandé)
- **Linux** (Ubuntu 20.04+, Debian 11+, etc.)
- **macOS** (10.15+)

### Mémoire et espace disque

- **RAM** : Minimum 4 Go (8 Go recommandé)
- **Espace disque** : 2 Go minimum pour l'installation complète
- **Processeur** : Processeur multi-cœurs recommandé pour les traitements par lots

### Logiciels requis

- **Python 3.11.x** (3.11 recommandé ; validé avec SciPy 1.17.1 et PHOEBE 2.4.22)
- **Conda** (Miniconda ou Anaconda) - Recommandé pour la gestion des environnements
- **Git** (optionnel, pour cloner le dépôt)

---

## 2. Installation de Python et Conda

### Option A : Installation avec Conda (RECOMMANDÉ)

Conda facilite la gestion des dépendances et des environnements Python.

1. **Télécharger Miniconda** :
   - Windows : https://docs.conda.io/en/latest/miniconda.html
   - Sélectionnez la version Python 3.11
   - Exécutez l'installateur et suivez les instructions

2. **Vérifier l'installation** :
   ```cmd
   conda --version
   python --version
   ```

### Option B : Installation avec Python seul

Si vous préférez utiliser Python sans Conda :

1. **Télécharger Python** :
   - Site officiel : https://www.python.org/downloads/
   - Version 3.11+ recommandée
   - **Important** : Cochez "Add Python to PATH" lors de l'installation

2. **Installer pip** (généralement inclus) :
   ```cmd
   python --version
   pip --version
   ```

---

## 3. Création de l'environnement Conda

### Avec Conda (Recommandé)

1. **Ouvrir un terminal Anaconda Prompt** (Windows) ou terminal (Linux/macOS)

2. **Créer un nouvel environnement nommé `astroenv`** :
   ```cmd
   conda create -n astroenv python=3.11
   ```
   
   **Important** : N'oubliez pas l'option `-n` (ou `--name`) pour spécifier le nom de l'environnement. Sans cette option, conda affichera une erreur.

3. **Activer l'environnement** :
   ```cmd
   conda activate astroenv
   ```

### Avec Python seul (venv)

1. **Naviguer vers le répertoire NPOAP** :
   ```cmd
   cd C:\Users\VotreNom\Documents\NPOAP
   ```

2. **Créer un environnement virtuel** :
   ```cmd
   python -m venv venv
   ```

3. **Activer l'environnement** :
   - **Windows** :
     ```cmd
     venv\Scripts\activate
     ```
   - **Linux/macOS** :
     ```bash
     source venv/bin/activate
     ```

---

## 4. Installation des dépendances Python

### Méthode recommandée : requirements.txt

1. **S'assurer que l'environnement est activé** :
   ```cmd
   conda activate astroenv
   ```

2. **Mettre à jour pip** :
   ```cmd
   python -m pip install --upgrade pip
   ```

3. **Installer les dépendances** :
   ```cmd
   pip install -r requirements.txt
   ```

   Cette commande installe toutes les dépendances requises, incluant :
   - Les bibliothèques astronomiques (astropy, photutils, astroquery, **ccdproc**, **reproject**)
   - Les outils de calcul scientifique (numpy, scipy, pandas)
   - **reproject** : nécessaire pour l’alignement WCS des images en réduction (reprojection sur une grille commune)
   - **exotethys**, **h5py**, **click** : limb darkening / modules associés (ex. workflow ExoTETHyS)
   - **specutils**, **lightkurve** : spectroscopie et séries temporelles (selon les onglets utilisés)
   - **stdpipe**, **synphot** : pipelines photométriques / SED
   - **ezpadova** (dépôt Git, ligne `git+https://github.com/mfouesneau/ezpadova` dans le fichier) : isochrones PARSEC pour l’analyse d’amas
   - Les bibliothèques pour l'extraction de catalogues (requests pour MPC, TESS EBS, Exoplanet.eu)
   - Les outils d'analyse (emcee, statsmodels, pylightcurve)
   - Les outils optionnels ou avancés (phoebe, rebound, ultranest, etc.)

   **Note** : Les packages **ezpadova** (isochrones PARSEC pour l'onglet Analyse d'amas) et **Prospector** ne sont pas sur PyPI ; ils doivent être installés séparément depuis GitHub (voir section 5).

   **KBMOD** : le `requirements.txt` à la racine du dépôt peut contenir une ligne active **`kbmod @ git+https://github.com/dirac-institute/kbmod.git`**. Sous **Windows natif**, `pip install -r requirements.txt` échoue souvent à cette étape (compilation). **À faire** : commentez temporairement cette ligne dans `requirements.txt`, réexécutez `pip install -r requirements.txt`, puis installez KBMOD **uniquement dans WSL** selon **`docs/INSTALL_KBMOD_WSL.md`** et `requirements-kbmod.txt`. Les distributions « profil » (réduction, exoplanètes, etc.) peuvent fournir un `requirements.txt` sans KBMOD.

### Installation manuelle des dépendances

Si vous préférez installer manuellement :

```cmd
pip install "numpy>=2.4.0,<2.5"
pip install "scipy>=1.17.1,<1.18"
pip install pandas>=1.5.0
pip install "astropy>=7.2.0,<7.3"
pip install photutils>=1.6.0
pip install astroquery>=0.4.6
pip install reproject>=0.8
pip install "matplotlib>=3.10.0,<3.11"
pip install requests>=2.33.0
pip install Pillow>=12.1.0
pip install emcee>=3.1.0
pip install "phoebe>=2.4.22,<2.5"
pip install reportlab>=3.6.0
pip install pylightcurve>=4.0.0
pip install statsmodels>=0.13.0
```

### Installation sur Linux : dépendances système

Sur Linux, vous devrez peut-être installer certaines dépendances système :

**Ubuntu/Debian** :
```bash
sudo apt-get update
sudo apt-get install python3-tk python3-dev build-essential
```

**Fedora/CentOS** :
```bash
sudo dnf install python3-tkinter python3-devel gcc gcc-c++
```

---

## 5. Installation des dépendances optionnelles

### PHOEBE2 (Étoiles Binaires)

PHOEBE2 permet la modélisation d'étoiles binaires à éclipses.

#### Sur Windows

**ATTENTION** : PHOEBE2 nécessite Microsoft Visual C++ Build Tools pour compiler son extension C++ (`libphoebe`).

1. **Installer Microsoft Visual C++ Build Tools** :
   - Téléchargez depuis : https://visualstudio.microsoft.com/visual-cpp-build-tools/
   - Installez le workload **"Desktop development with C++"**
   - Composants essentiels :
     - MSVC v143 - VS 2022 C++ x64/x86 build tools
     - Windows 10/11 SDK
     - C++ CMake tools for Windows (recommandé)

2. **Installer PHOEBE2 dans l'environnement `astroenv` (Python 3.11)** :
   ```cmd
   conda activate astroenv
   pip install --force-reinstall phoebe
   ```

   - Sous Windows, NPOAP est testé avec **PHOEBE 2.4.22** et **SciPy 1.17.1** (voir `requirements.txt`).
   - Si l'installation échoue avec un message "Microsoft Visual C++ 14.0 or greater is required", vérifiez que les Build Tools sont bien installés, puis relancez la commande.

#### Sur Linux/macOS

```bash
conda activate astroenv
pip install phoebe
```

### Prospector (Analyse de Spectres de Galaxies - Optionnel)

Prospector est un outil Python pour inférer les propriétés stellaires à partir de données spectroscopiques et photométriques (SED - Spectral Energy Distribution). Il utilise l'inférence bayésienne avec des modèles de populations stellaires simples (SSP).

**IMPORTANT** : le paquet **astro-prospector** doit être installé depuis **GitHub** (`bd-j/prospector`), pas depuis PyPI. Le nom `prospector` sur PyPI désigne un **autre** logiciel (analyse statique de code Python).

#### Faut-il garder `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat` et `.ps1` ?

**Oui**, tant que vous utilisez NPOAP sous Windows avec l’environnement Conda **`astroenv`** : ces scripts installent **Prospector dans ce même Python** que `python main.py`. Le fichier **`.bat`** est un lanceur (double-clic, politique d’exécution PowerShell) ; le **`.ps1`** contient toute la logique et les paramètres (`-InstallFSPS`, `-CondaEnv autre_env`, `-SkipVerification`, `-ForceReinstall`, etc.). Vous pouvez n’utiliser que le `.ps1` si vous préférez, mais le `.bat` reste utile pour les utilisateurs qui n’ouvrent pas PowerShell à la main.

#### Voie A — Windows natif + Conda `astroenv` (recommandée pour NPOAP)

**Méthode recommandée** :

```cmd
INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat
```

Équivalent (PowerShell, depuis le dossier NPOAP) :

```powershell
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1
```

Avec tentative de compilation **FSPS** (nécessite **CMake** et **gfortran** sur le PATH Windows) :

```powershell
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 -InstallFSPS
```

Sans **gfortran**, l’installation de la roue **fsps** échoue en général : le script enchaîne avec `pip install … --no-deps`, crée **SPS_HOME** et des **fichiers stub** ; Prospector reste importable et une grande partie des fonctionnalités fonctionne.

**Vérification** :

```cmd
conda activate astroenv
python -c "import prospect; from prospect.models import SpecModel; print('Prospector OK!')"
```

#### Voie B — WSL2 (Linux) : FSPS / Prospector dans un venv dédié

Pour un environnement **Linux** avec **gfortran** et souvent une installation **fsps** plus simple, utilisez le dossier **`Installation_fsps/`** :

- **`Installation_fsps/prospector.bat`** (Windows) appelle **`Installation_fsps/install_prospector_wsl.sh`** dans WSL.
- Prérequis : **WSL2** et une distro type **Ubuntu** (`install_wsl.bat`, `install_ubuntu_wsl.bat` si besoin).

Cette voie installe un **venv séparé** sous Linux (`~/.local/share/npoap-prospector-wsl/venv` par défaut) : ce n’est **pas** le même interpréteur que **`astroenv`** sous Windows. NPOAP lancé depuis Windows n’utilise pas automatiquement ce venv ; la voie B sert surtout à du travail Prospector/FSPS en ligne de commande sous WSL. Certaines archives placent `prospector.bat` et `install_prospector_wsl.sh` à la **racine** du dossier : gardez les **deux** fichiers côte à côte et exécutez le `.bat`.

#### Voie C — Installation manuelle (Windows, dans `astroenv`)

```cmd
conda activate astroenv

pip install numpy>=1.20.0 scipy>=1.7.0 pandas>=1.3.0 astropy>=5.0.0
pip uninstall sedpy -y
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir
pip install dynesty>=2.0.0 dill>=0.3.0 h5py>=3.0.0 emcee>=3.1.0 corner>=2.2.3 matplotlib>=3.8.4
pip install git+https://github.com/bd-j/prospector.git --no-cache-dir
```

Si la dernière ligne échoue à cause de **fsps** : `pip install git+https://github.com/bd-j/prospector.git --no-deps` puis installer **fsps** séparément ou s’en tenir aux stubs (comme le script automatique).

**FSPS** : optionnel mais utile pour les modèles SED complets. Détails, ordre des paquets et dépannage : **`docs/PROTOCOLE_INSTALLATION_PROSPECTOR_WINDOWS.md`**.

### CuPy (Accélération GPU - Optionnel - NON REQUIS)

CuPy permet d'accélérer certains calculs astrométriques avec une carte graphique NVIDIA, notamment la détection d'étoiles dans l'onglet **Photométrie Astéroïdes** (astrométrie zero-aperture).

**Cartes compatibles** :
- **NVIDIA GeForce GTX 1660 Ti** ✅ (et autres cartes NVIDIA avec support CUDA)
- Architecture Turing (Compute Capability 7.5)
- Toutes les cartes NVIDIA récentes (GTX 1000+, RTX série)

**Prérequis** :
- Carte graphique NVIDIA avec support CUDA
- CUDA Toolkit installé (version 11.x ou 12.x recommandée)

**Étape 1 : Vérifier votre version CUDA**

1. **Ouvrir le gestionnaire de périphériques Windows** :
   - Appuyez sur `Win + X` et sélectionnez "Gestionnaire de périphériques"
   - Développez "Adaptateurs d'affichage"
   - Vérifiez que votre carte NVIDIA est reconnue (ex: "NVIDIA GeForce GTX 1660 Ti")

2. **Vérifier la version CUDA installée** :
   ```cmd
   nvidia-smi
   ```
   - La version CUDA apparaît en haut à droite (ex: "CUDA Version: 12.2")
   - Si `nvidia-smi` n'est pas reconnu, installez les **pilotes NVIDIA** depuis https://www.nvidia.com/drivers/

3. **Si CUDA n'est pas installé** :
   - Téléchargez **CUDA Toolkit** depuis https://developer.nvidia.com/cuda-downloads
   - Pour GTX 1660 Ti, **CUDA 11.8** ou **12.x** est recommandé
   - Installez le toolkit CUDA (suivez les instructions du site NVIDIA)

**Étape 2 : Installer CuPy**

Une fois CUDA installé, installez CuPy correspondant à votre version CUDA :

**Pour CUDA 12.x**  :
```cmd
pip install cupy-cuda12x
```

**Pour CUDA 13.x** (CUDA 13.0+) :
```cmd
pip install cupy-cuda13x
```

**Vérification de l'installation** :
```cmd
python -c "import cupy as cp; print(f'CuPy {cp.__version__} installé'); print(f'GPU disponible: {cp.cuda.is_available()}'); print(f'GPU: {cp.cuda.Device(0).compute_capability if cp.cuda.is_available() else \"N/A\"}')"
```

Si l'installation réussit, vous verrez :
```
CuPy 13.x.x installé
GPU disponible: True
GPU: (7, 5)
```

**Note** : CuPy est optionnel. L'application fonctionne sans GPU, mais sera plus lente pour certaines opérations d'astrométrie (détection d'étoiles dans les images). Pour la GTX 1660 Ti, l'accélération GPU peut améliorer les performances de 2-5x pour la détection d'étoiles sur de grandes images.

**Métadonnées pip / import CuPy** : certaines installations mélangées (conda/pip, wheel incomplet) font que `import cupy` lève une erreur du type `AttributeError: 'NoneType' object has no attribute 'get'` dans `cupy._detect_duplicate_installation`. NPOAP applique un **contournement au démarrage** dans `main.py` (filtrage des entrées `importlib.metadata` sans métadonnées), et l’onglet photométrie astéroïdes **désactive le GPU** proprement si CuPy ne se charge pas. Si le problème persiste, voir la section **Dépannage** (CuPy).

### ezpadova / PARSEC (Analyse d'amas – Optionnel)

L'onglet **Analyse d'amas** permet d'estimer l'âge et la distance d'un amas à partir du diagramme couleur-magnitude (G, G_BP, G_RP). Pour utiliser les **isochrones PARSEC** (modèle Padova), le package **ezpadova** doit être installé. ezpadova interroge le serveur CMD (stev.oapd.inaf.it) et fournit des isochrones en filtres Gaia EDR3.

**IMPORTANT** : ezpadova n'est **pas** disponible sur PyPI. Il faut l'installer depuis GitHub (Git doit être installé et accessible dans le PATH).

**Installation** :

```cmd
conda activate astroenv
pip install git+https://github.com/mfouesneau/ezpadova
```

**Sans ezpadova** : l'onglet Analyse d'amas fonctionne quand même ; l'estimation d'âge utilise alors une grille empirique de tour de courbe au lieu des isochrones PARSEC.

**Vérification** :
```cmd
python -c "import ezpadova; r = ezpadova.get_isochrones(logage=(9, 9, 1), MH=(0, 0, 1), photsys_file='gaiaEDR3'); print('ezpadova / PARSEC OK')"
```

---

## 6. Installation d'Astrometry.net (Optionnel) [ne pas confondre avec NOVA.Astrometry.net qui nécessite internet et une clé]

Astrometry.net permet l'astrométrie locale (sans connexion internet).

### Sur Windows (via WSL)

1. **Installer WSL (Windows Subsystem for Linux)** :
   ```powershell
   wsl --install
   ```
   Redémarrez votre ordinateur après l'installation.

2. **Installer Ubuntu dans WSL** :
   - Après le redémarrage, configurez Ubuntu
   - Créez un nom d'utilisateur et un mot de passe

3. **Installer Astrometry.net dans WSL** :
   ```bash
   sudo apt-get update
   sudo apt-get install astrometry.net
   ```

4. **Télécharger les index** (optionnel mais recommandé) :
   - Les index sont volumineux (plusieurs Go)
   - Consultez : http://astrometry.net/doc/readme.html#getting-index-files
   - Placez-les dans `/usr/share/astrometry/data/` dans WSL

### Sur Linux

```bash
sudo apt-get update
sudo apt-get install astrometry.net
```

### Sur macOS

```bash
brew install astrometry-net
```

**Note** : Astrometry.net est optionnel. L'application peut utiliser Astrometry.net en ligne (NOVA) si vous avez une clé API.

---

## 7. KBMOD – Détection Synthetic Tracking (optionnel, via WSL)

KBMOD permet la détection d'astéroïdes par *Synthetic Tracking* (empilement d'images déplacées selon une trajectoire). **KBMOD n'est pas installé sous Windows** (échec de compilation avec MSVC). NPOAP lance la détection **via WSL** : le script `scripts/kbmod_wsl_detect.py` s'exécute dans l'environnement Linux et écrit les candidats dans `kbmod_candidates.csv`, que NPOAP lit ensuite.

### Prérequis

- **WSL 2** avec une distribution Linux (Ubuntu recommandée).
- **GPU NVIDIA** (optionnel mais recommandé pour des temps de calcul raisonnables).
- **CUDA** installé sous WSL si vous utilisez le GPU.

### Installation sous WSL

Suivez le guide dédié : **`docs/INSTALL_KBMOD_WSL.md`**.

Résumé des étapes dans WSL :

1. Installer les outils de compilation et Python : `build-essential`, `cmake`, `python3-dev`, `python3-pip`.
2. (Optionnel) Installer le toolkit CUDA sous WSL.
3. Cloner et installer KBMOD : `git clone https://github.com/dirac-institute/kbmod.git && cd kbmod && pip3 install -e .`
4. Vérifier : `python3 -c "import kbmod.search; print('KBMOD OK')"`

### Utilisation dans NPOAP

Dans l'onglet **Astéroïdes**, après avoir chargé un dossier d'images FITS, cliquez sur **« Détection KBMOD (via WSL) »**. NPOAP appelle automatiquement `wsl python3 .../scripts/kbmod_wsl_detect.py` avec le dossier FITS (converti en chemin WSL) et affiche les candidats une fois le script terminé. Vous pouvez en sélectionner un comme cible T1 pour la photométrie.

**Note** : le `requirements.txt` racine peut lister **kbmod** ; sous Windows natif, commentez cette ligne pour finir l’installation pip, puis utilisez **`requirements-kbmod.txt`** et **`docs/INSTALL_KBMOD_WSL.md`** pour KBMOD **dans WSL** uniquement.

---

## 8. Configuration initiale

1. **télécharger NPOAP** :
   
   -  extrayez l'archive ZIP NPOAP_XXXXXX_XXXX.ZIP

2. **Vérifier la structure des fichiers** :
   ```
   NPOAP/
   ├── main.py
   ├── requirements.txt
   ├── config.py
   ├── gui/
   ├── core/
   └── utils/
   ```

3. **Configuration de la clé API Astrometry.net (optionnel)** :
   - Lancez l'application (voir section 9)
   - Allez dans l'onglet "Accueil"
   - Entrez votre clé API Astrometry.net (obtenue sur https://nova.astrometry.net/)
   - Cliquez sur "Sauvegarder"

4. **Installation de HOPS dans l'onglet Photométrie Exoplanètes** :
   - Ouvrez NPOAP, onglet **Photométrie Exoplanètes**
   - Cliquez sur **Installer / Réinstaller HOPS (ZIP)**
   - Sélectionnez l'archive `HOPS-modified.zip` (ou `hops-master.zip` si vous utilisez l'archive officielle)
   - L'installation est faite dans : `external_apps/hops/hops-master`
   - Cliquez ensuite sur **Lancer HOPS**

   **Note** : au premier lancement, les dépendances HOPS peuvent être installées automatiquement (ex. `exoclock`).

### Filtres photométriques (choix et installation)

Le choix du filtre doit rester cohérent entre acquisition, calibration et ajustement :

- **Réduction/calibration** : les flats doivent être pris avec un filtre compatible avec les lights (voir les contrôles dans l'onglet Réduction).
- **HOPS Data & Target** : le filtre sélectionné est utilisé ensuite par la photométrie et le fitting.
- **Limb-darkening (HOPS)** : il dépend du passband connu par `pylightcurve41`.

#### Cas des filtres Gaia (`Gaia G`, `Gaia BP`, `Gaia RP`)

- NPOAP expose ces filtres dans la liste HOPS.
- Dans la distribution NPOAP, les fichiers passband Gaia sont déjà fournis.
- Si les passbands Gaia sont installés dans l'environnement Python, HOPS calcule les LD avec `GAIA_G`, `GAIA_BP`, `GAIA_RP`.
- Sinon, NPOAP applique automatiquement un repli pour éviter l'échec du fitting :
  - `Gaia G` -> `JOHNSON_V`
  - `Gaia BP` -> `JOHNSON_B`
  - `Gaia RP` -> `COUSINS_R`

---

## 9. Vérification de l'installation

### Vérifier les dépendances Python

Le dépôt fournit le script **`test_installation.py`** à la **racine du projet** (à côté de `main.py`). Il vérifie les modules requis et optionnels (numpy, astropy, photutils, **reproject**, phoebe, **specutils**, **exotethys**, CuPy, Prospector, FSPS, modules locaux NPOAP, etc.) et affiche un résumé dans le terminal.

**Important** : au démarrage, ce script applique le même type de **correctif `importlib.metadata`** que `main.py`, afin que la vérification de **CuPy** ne plante pas sur des métadonnées pip incohérentes (`AttributeError` / `_detect_duplicate_installation`).

Depuis le dossier NPOAP, avec l’environnement activé :

```cmd
conda activate astroenv
cd C:\Users\VotreNom\Documents\NPOAP
python test_installation.py
```

En cas d’échec sur un module optionnel (CuPy, Prospector…), l’application peut tout de même fonctionner pour les autres onglets ; lisez les messages et indices affichés par le script.

### Vérifier Astrometry.net (si installé)

Dans WSL (Windows) ou terminal (Linux) :
```bash
solve-field --version
```

---

## 10. Lancement de l'application

**Important** : NPOAP doit être lancé dans l’environnement **astroenv** (ou un environnement contenant toutes les dépendances, dont **reproject** pour l’alignement). Les distributions préconstruites (réduction, exoplanètes, etc.) incluent un script qui active astroenv automatiquement.

### Méthode 1 : Distribution préconstruite (recommandé)

Si vous utilisez une distribution NPOAP (dossier « reduction », « exoplanets », « full », etc., généré par le build) :

1. **Prérequis** : [Miniconda](https://docs.conda.io/en/latest/miniconda.html) ou Anaconda installé (Python **3.11** recommandé pour l’environnement **astroenv**).
2. Ouvrez le dossier de la distribution (il contient `main.py` et un `requirements.txt` propre au **profil** : réduction, exoplanètes, complet, etc.).
3. **Première installation** : double-cliquez sur **`INSTALLER_NPOAP_ASTROENV_WINDOWS.bat`**. Le script crée l’environnement conda **`astroenv`** s’il n’existe pas (`python=3.11`), l’active, puis exécute **`pip install -r requirements.txt`** pour ce profil.  
   - *Rétro-compatibilité* : le profil « full » peut aussi fournir **`INSTALLER_DISTRIBUTION_FULL_WINDOWS.bat`**, qui appelle le même installateur.
4. **Lancements suivants** : double-cliquez sur **`lancement.bat`**. Ce script **exige** que **astroenv** soit utilisable (il ne retombe plus sur un Python « nu » du PATH). En cas d’erreur, réexécutez l’installateur ou ouvrez une **Anaconda Prompt** et vérifiez `conda activate astroenv` (voir Méthode 2).

Les anciennes distributions qui utilisaient un dossier **`venv`** local ne sont plus le mode par défaut : tout est centralisé dans **astroenv**.

### Méthode 2 : Projet complet (script fourni ou manuel)

Un script **`LANCER_NPOAP_ASTROENV.bat`** est fourni à la racine du projet :

```cmd
LANCER_NPOAP_ASTROENV.bat
```

**Ou en manuel** :

1. **Activer l'environnement** :
   ```cmd
   conda activate astroenv
   ```

2. **Lancer l'application** :
   ```cmd
   cd C:\Users\VotreNom\Documents\NPOAP
   python main.py
   ```

### Méthode 3 : Depuis Python (développement)

Pour un comportement identique au lancement utilisateur (journalisation, correctif métadonnées **avant** l’import de l’interface), préférez **`python main.py`**.

Si vous importez directement `MainWindow` dans un interpréteur ou un script personnalisé, vous ne bénéficiez pas du contournement CuPy/`importlib.metadata` défini en tête de **`main.py`** ; en cas d’erreur à l’import, lancez plutôt `main.py` ou appelez le même patch qu’au début de `main.py` avant `from gui.main_window import MainWindow`.

```python
from gui.main_window import MainWindow
import tkinter as tk
import utils.logging_handler as logging_handler

logging_handler.setup_logging()
root = tk.Tk()
app = MainWindow(root)
root.mainloop()
```

---

## 11. Dépannage

### Erreur : "ModuleNotFoundError: No module named 'tkinter'"

**Windows** : tkinter devrait être inclus. Si l'erreur persiste, réinstallez Python et cochez "tcl/tk and IDLE".

**Linux** :
```bash
sudo apt-get install python3-tk
```

**macOS** : tkinter devrait être inclus avec Python.

### Erreur lors de l'installation de PHOEBE2 sur Windows

**Erreur** : "Microsoft Visual C++ 14.0 or greater is required"

**Solution** :
1. Installez Microsoft Visual C++ Build Tools (voir section 5.1)
2. Redémarrez votre ordinateur
3. Réessayez : `pip install phoebe`

**Alternative** : Utilisez conda :
```cmd
conda install -c conda-forge phoebe
```

### Erreur : "CRASH WSL (Code 1)" lors de l'astrométrie locale

**Causes possibles** :
- Astrometry.net n'est pas installé dans WSL
- Les index ne sont pas téléchargés
- Problème de permissions

**Solution** :
1. Vérifiez l'installation dans WSL : `wsl solve-field --version`
2. Testez manuellement : `wsl solve-field image.fits`
3. Vérifiez les logs dans le dossier `logs/`

### Logs HOPS runtime

Pour diagnostiquer les erreurs HOPS pendant le déroulement des opérations, consultez :

`external_apps/hops/hops-master/runtime_logs/hops_runtime.log`

Le journal est configuré en mode **erreurs uniquement** (avec traces d'exception détaillées).

### Erreur KBMOD (via WSL) : "KBMOD not installed in this environment" ou "Script introuvable"

**Cause** : Le script `scripts/kbmod_wsl_detect.py` est exécuté sous WSL mais KBMOD n’est pas installé dans l’environnement Python utilisé par `wsl python3`, ou le script est absent.

**Solution** :
1. Vérifiez que le dossier `scripts/` contient `kbmod_wsl_detect.py` (à la racine du projet NPOAP).
2. Installez KBMOD **dans WSL** en suivant `docs/INSTALL_KBMOD_WSL.md`.
3. Dans WSL, testez : `python3 -c "import kbmod.search; print('KBMOD OK')"`.

### Erreur : "No module named 'reproject'" (alignement WCS)

**Cause** : Le module **reproject** est requis pour l’alignement des images (onglet Réduction → Aligner images WCS). Il doit être installé dans le même environnement que celui utilisé au lancement (généralement **astroenv**).

**Solution** :
```cmd
conda activate astroenv
pip install reproject
```
Puis relancez NPOAP avec **lancement.bat** (ou après `conda activate astroenv`). Si vous lancez par double-clic, assurez-vous que le script active bien astroenv (voir section 10).

### Erreur : "No module named 'readline'" (PHOEBE2 sur Windows)

Cette erreur est normale sur Windows. Elle est gérée automatiquement par l'application. Si vous la voyez, l'application devrait quand même fonctionner.

### Erreur CuPy : `'NoneType' object has no attribute 'get'` ou `_detect_duplicate_installation`

**Contexte** : au premier `import cupy`, la bibliothèque inspecte les paquets installés via `importlib.metadata`. Si certaines entrées ont **`metadata`** nul (installation pip/conda incohérente, paquet corrompu), CuPy peut lever `AttributeError` avant même que l’application s’ouvre.

**Que fait NPOAP** : `main.py` remplace temporairement l’itérateur `importlib.metadata.distributions` pour ignorer ces entrées ; l’onglet photométrie astéroïdes accepte aussi l’absence de CuPy (astrométrie en CPU).

**Si l’erreur apparaît encore** (script lancé sans passer par `main.py`, environnement très endommagé) :

1. Exécutez `python test_installation.py` (il applique le même correctif) pour diagnostiquer.
2. `pip check` ; corrigez les paquets en conflit signalés.
3. Désinstallez toutes les variantes CuPy puis n’en réinstallez qu’**une** :  
   `pip uninstall cupy cupy-cuda11x cupy-cuda12x cupy-cuda13x -y` puis par exemple `pip install cupy-cuda12x` selon votre CUDA (voir section CuPy).

### Erreur lors de l'installation de Prospector

**Erreur** : `No CMAKE_Fortran_COMPILER could be found` / échec de compilation de **fsps** sous Windows

**Cause** : il manque **gfortran** (Fortran). CMake et MSVC peuvent être présents sans Fortran.

**Conséquence habituelle** : l’installateur PowerShell enchaîne avec `pip install prospector --no-deps` et des stubs **SPS_HOME** ; Prospector peut quand même s’installer. Pour un **fsps** complet : installez gfortran (ex. MSYS2 / MinGW) et relancez avec `-InstallFSPS`, ou utilisez la voie **WSL** (`Installation_fsps/prospector.bat`).

**Erreur** : `ModuleNotFoundError: No module named 'sedpy.observate'`

**Cause** : `sedpy` a été installé depuis PyPI au lieu de GitHub.

**Solution** :
```cmd
pip uninstall sedpy -y
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir
```

**Erreur** : `ModuleNotFoundError: No module named 'astropy'`

**Solution** :
```cmd
pip install astropy>=5.0.0
```

**Erreur** : `ValueError: could not assign tuple of length 1 to structure with 10 fields`

**Cause** : Format incorrect du fichier stub FSPS.

**Solution** : Supprimez le fichier stub et relancez l'installation :
```cmd
del "%USERPROFILE%\.local\share\fsps\dust\Nenkova08_y010_torusg_n10_q2.0.dat"
# Puis relancez INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat
```

Pour plus de détails, consultez `docs/PROTOCOLE_INSTALLATION_PROSPECTOR_WINDOWS.md`.

### L'application est lente

**Solutions** :
- Utilisez CuPy si vous avez une carte graphique NVIDIA (voir section 5.2)
- Réduisez le nombre d'images traitées par lot
- Vérifiez que vous utilisez un processeur récent
- Fermez les autres applications lourdes

### Problèmes de connexion avec Astrometry.net

**Si l'astrométrie en ligne échoue** :
- Vérifiez votre connexion internet
- Vérifiez votre clé API dans l'onglet "Accueil"
- Consultez https://nova.astrometry.net/ pour le statut du service

### Erreur d'importation de modules locaux

**Erreur** : "ModuleNotFoundError: No module named 'core'"

**Solution** :
1. Assurez-vous d'être dans le répertoire NPOAP
2. Vérifiez que les dossiers `core/`, `gui/`, `utils/` existent
3. Vérifiez que `__init__.py` existe dans ces dossiers

---

## Récapitulatif des commandes essentielles

### Installation complète (Conda)

```cmd
# Créer l'environnement
conda create -n astroenv python=3.11
conda activate astroenv

# Installer les dépendances
pip install -r requirements.txt

# Optionnel : PHOEBE2
pip install phoebe

# Lancer l'application
python main.py
```

### Mise à jour

```cmd
conda activate astroenv
pip install --upgrade -r requirements.txt
```

### Installation de Prospector (optionnel)

```cmd
conda activate astroenv
REM Script Windows (installe dans astroenv) :
INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat

REM Alternative WSL : Installation_fsps\prospector.bat (venv Linux separe)

REM Ou installation manuelle dans astroenv :
pip install numpy>=1.20.0 scipy>=1.7.0 pandas>=1.3.0 astropy>=5.0.0
pip uninstall sedpy -y
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir
pip install dynesty>=2.0.0 dill>=0.3.0 h5py>=3.0.0 emcee>=3.1.0 corner>=2.2.3 matplotlib>=3.8.4
pip install git+https://github.com/bd-j/prospector.git --no-cache-dir
```

### Désinstallation

```cmd
conda deactivate
conda env remove -n astroenv
```

---

## Support et ressources

- **Documentation utilisateur** : `docs/MANUEL_UTILISATEUR.md`
- **Protocole d'installation Prospector** : `docs/PROTOCOLE_INSTALLATION_PROSPECTOR_WINDOWS.md`
- **Installation KBMOD sous WSL** : `docs/INSTALL_KBMOD_WSL.md`
- **Logs de l'application** : Dossier `logs/`
- **Configuration** : `config.py` et `config.json`

---

## 12. Remerciements (Acknowledgments)

NPOAP utilise de nombreuses bibliothèques et outils open-source de la communauté scientifique. Cette section résume les principaux remerciements ; la liste complète et détaillée figure dans **`docs/ACKNOWLEDGEMENTS.md`**.

### Bibliothèques Python principales

- **Astropy** : Bibliothèque Python pour l'astronomie (https://www.astropy.org/)
  - Utilisée pour la manipulation des images FITS, les coordonnées célestes, les transformations WCS, et les calculs astronomiques

- **photutils** : Outils de photométrie astronomique (https://photutils.readthedocs.io/)
  - Utilisée pour la détection d'étoiles, la photométrie par ouverture, et l'estimation du fond de ciel

- **astroquery** : Interface Python pour les services astronomiques en ligne (https://astroquery.readthedocs.io/)
  - Utilisée pour interroger les catalogues Gaia (Vizier) et les éphémérides JPL Horizons

- **NumPy, SciPy, Pandas** : Bibliothèques fondamentales pour le calcul scientifique en Python
  - Utilisées pour les calculs numériques, l'analyse de données, et la manipulation de tableaux

- **Matplotlib** : Bibliothèque de visualisation (https://matplotlib.org/)
  - Utilisée pour l'affichage graphique des images, courbes de lumière, et périodogrammes

- **PHOEBE2** : Bibliothèque Python pour la modélisation d'étoiles binaires à éclipses (https://phoebe-project.org/)
  - Utilisée pour la modélisation et l'analyse des systèmes binaires dans l'onglet "Étoiles Binaires"

- **emcee** : Bibliothèque MCMC pour l'ajustement de modèles (https://emcee.readthedocs.io/)
  - Utilisée pour l'analyse des variations de temps de transit (TTV)

- **statsmodels** : Bibliothèque d'analyse statistique (https://www.statsmodels.org/)
  - Utilisée pour les analyses statistiques avancées des courbes de lumière

- **pylightcurve** : Bibliothèque Python pour la modélisation et l'analyse de courbes de lumière d'exoplanètes (https://github.com/ucl-exoplanets/pylightcurve)
  - Développée par l'équipe UCL Exoplanets (University College London)
  - Utilisée pour la modélisation des transits d'exoplanètes, le calcul des coefficients d'assombrissement du limbe, et l'ajustement flexible de courbes de lumière multi-époques
  - Licence MIT

- **Pillow (PIL)** : Bibliothèque de traitement d'images (https://python-pillow.org/)
  - Utilisée pour le traitement et la manipulation d'images

- **reportlab** : Bibliothèque de génération de PDF (https://www.reportlab.com/)
  - Utilisée pour la génération de documentation PDF

- **requests** : Bibliothèque HTTP pour Python (https://requests.readthedocs.io/)
  - Utilisée pour les requêtes HTTP vers les services en ligne

- **specutils** : Bibliothèque Python pour l'analyse de données spectroscopiques (https://specutils.readthedocs.io/)
  - Utilisée pour la représentation, le chargement, la manipulation et l'analyse de spectres d'étoiles
  - Développée dans le cadre du projet Astropy

- **rebound** : Bibliothèque de simulation N-body pour systèmes planétaires (https://rebound.readthedocs.io/)
  - Utilisée pour les simulations gravitationnelles dans l'analyse de systèmes multiples (onglet Analyse de Données)
  - Permet de modéliser les interactions gravitationnelles entre planètes

- **ultranest** : Bibliothèque d'échantillonnage bayésien avec nested sampling (https://johannesbuchner.github.io/UltraNest/)
  - Utilisée pour le fitting bayésien des modèles N-body aux observations TTV
  - Alternative robuste aux méthodes MCMC pour l'estimation de paramètres

- **STDPipe** : Simple Transient Detection Pipeline (https://stdpipe.readthedocs.io/)
  - Auteur principal : Sergey Karpov
  - Utilisée pour la photométrie des transitoires : astrométrie automatique, soustraction d'images, détection de transitoires, photométrie calibrée
  - Fournit des méthodes avancées de détection (segmentation, DAOStarFinder, IRAFStarFinder)
  - Permet le téléchargement d'images de référence depuis Pan-STARRS, SDSS, DES
  - Licence MIT

### Services et catalogues externes

- **Gaia DR3** : Catalogue astrométrique de l'Agence Spatiale Européenne (ESA)
  - Utilisé comme référence astrométrique et photométrique (https://www.cosmos.esa.int/web/gaia)

- **Astrometry.net** : Service d'astrométrie automatique (http://astrometry.net/)
  - Utilisé pour la résolution de champs stellaires (plate solving) via NOVA (service en ligne) et solve-field (installation locale)

- **JPL Horizons** : Service d'éphémérides du Jet Propulsion Laboratory (NASA)
  - Utilisé pour obtenir les éphémérides des astéroïdes et planètes (https://ssd.jpl.nasa.gov/horizons/)

- **Vizier** : Service de catalogues astronomiques du Centre de données astronomiques de Strasbourg (CDS)
  - Utilisé pour interroger les catalogues Gaia et autres catalogues astronomiques (https://vizier.cds.unistra.fr/)

### Outils et frameworks

- **Python** : Langage de programmation (https://www.python.org/)
  - Langage principal de développement

- **Conda/Miniconda** : Système de gestion d'environnements et de paquets (https://docs.conda.io/)
  - Utilisé pour la gestion des environnements Python et des dépendances

- **Tkinter** : Bibliothèque d'interface graphique (incluse dans Python)
  - Utilisée pour l'interface utilisateur graphique

- **WSL (Windows Subsystem for Linux)** : Système pour exécuter Linux sur Windows (Microsoft)
  - Utilisé pour l'installation locale d'Astrometry.net sur Windows

### Méthodologie astrométrie zero-aperture

- **Zero-Aperture-Astrometry** (Ben Sharkey) : https://github.com/bensharkey/Zero-Aperture-Astrometry
  - Méthodologie d’extrapolation des positions astrométriques à aperture nulle (fit RA/Dec vs aperture, extrapolation à 0) ; NPOAP s’en inspire dans l’onglet Photométrie Astéroïdes

### Références bibliographiques

- **Farnocchia et al. (2022)** : "International Asteroid Warning Network Timing Campaign: 2019 XS", *Planetary Science Journal*, 3:156
  - Guide de référence pour les améliorations astrométriques implémentées dans NPOAP
  - DOI: https://doi.org/10.3847/PSJ/ac7224

- **Prša (2018)** : "Modeling and Analysis of Eclipsing Binary Stars: The Theory and Design Principles of PHOEBE"
  - Référence principale pour l'intégration de PHOEBE2 (DOI: 10.1088/978-0-7503-1287-5)

### Licences

NPOAP est développé en utilisant des bibliothèques sous diverses licences open-source (principalement BSD, MIT, et Apache 2.0). Nous respectons les licences de toutes les bibliothèques utilisées.

Pour plus d'informations sur les licences spécifiques, consultez :
- Les fichiers LICENSE ou COPYRIGHT de chaque bibliothèque
- Les sites web officiels des projets mentionnés ci-dessus

Pour les remerciements complets (bibliothèques, services, références, licences), consultez **`docs/ACKNOWLEDGEMENTS.md`**.

---

**NPOAP - Manuel d'Installation v1.0**

*Pour toute question ou problème, consultez la section Dépannage ou contactez l'équipe de développement.*

