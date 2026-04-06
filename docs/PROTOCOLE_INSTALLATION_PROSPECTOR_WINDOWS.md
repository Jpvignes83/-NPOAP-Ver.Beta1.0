# Protocole d'Installation Prospector pour Windows

## Vue d'ensemble

Ce document décrit le protocole complet d'installation de Prospector sur Windows, basé sur l'analyse approfondie des logs d'installation réussis et des erreurs rencontrées. Prospector est un outil Python pour inférer les propriétés stellaires à partir de données spectroscopiques et photométriques (SED - Spectral Energy Distribution).

**Prospector GitHub** : https://github.com/bd-j/prospector  
**Documentation Prospector** : https://prospect.readthedocs.io/

---

## Table des matières

1. [Enseignements Clés](#1-enseignements-clés)
2. [Ordre d'Installation Critique](#2-ordre-dinstallation-critique)
3. [Problèmes Identifiés et Solutions](#3-problèmes-identifiés-et-solutions)
4. [Scripts d'Installation](#4-scripts-dinstallation)
5. [Installation Manuelle](#5-installation-manuelle)
6. [Format du Fichier Stub FSPS](#6-format-du-fichier-stub-fsps)
7. [Vérification Post-Installation](#7-vérification-post-installation)
8. [Intégration dans l'Installation Globale NPOAP](#8-intégration-dans-linstallation-globale-npoap)
9. [Dépannage](#9-dépannage)
10. [Installation de FSPS (Optionnel)](#10-installation-de-fsps-optionnel)
11. [Références](#11-références)

---

## 1. Enseignements Clés

### Analyse des Logs d'Installation

L'analyse complète des logs d'installation WSL et des erreurs rencontrées a révélé plusieurs points critiques :

1. **Ordre d'installation est crucial** : Certaines dépendances doivent être installées avant d'autres
2. **Sources d'installation importantes** : `sedpy` et `prospector` doivent venir de GitHub, pas de PyPI
3. **Dépendances manquantes** : `astropy` doit être explicitement installé (requis par Prospector)
4. **Format des fichiers stub** : Le format exact des fichiers stub FSPS est critique (3 espaces, 129 lignes)
5. **Initialisation des sous-modules Git** : FSPS nécessite `git submodule update --init --recursive`

### Points Critiques

- ✅ **Ordre d'installation** : Dépendances de base → sedpy (GitHub) → Autres dépendances → FSPS/SPS_HOME → Prospector (GitHub)
- ✅ **sedpy depuis GitHub** : Le package PyPI n'a pas le module `observate`
- ✅ **astropy explicite** : Dépendance nécessaire mais parfois omise
- ✅ **FSPS optionnel** : Prospector fonctionne avec des fichiers stub si FSPS n'est pas installé
- ✅ **SPS_HOME défini tôt** : Doit être défini avant l'import de `prospect`

---

## 2. Ordre d'Installation Critique

L'ordre suivant est **essentiel** pour garantir une installation réussie :

### Étape 1 : Prérequis Système

- **Conda/Miniconda** : Installé et fonctionnel
- **Git** : Installé et dans le PATH (requis pour installer depuis GitHub)
- **Python 3.10 ou 3.11** : Dans l'environnement Conda

### Étape 2 : Dépendances Python de Base

```cmd
pip install numpy>=1.20.0 scipy>=1.7.0 pandas>=1.3.0 astropy>=5.0.0
```

**⚠️ IMPORTANT** : `astropy>=5.0.0` est **crucial** - Prospector en a besoin mais cette dépendance peut être omise si non explicitement listée.

### Étape 3 : Installation de sedpy depuis GitHub

```cmd
# Désinstaller sedpy de PyPI s'il existe
pip uninstall sedpy -y

# Installer depuis GitHub
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir
```

**⚠️ CRITIQUE** : Le package `sedpy` sur PyPI n'a **pas** le module `observate` requis par Prospector. Il **DOIT** être installé depuis GitHub.

**Vérification** :
```cmd
python -c "from sedpy import observate; print('OK')"
```

### Étape 4 : Autres Dépendances Prospector

```cmd
pip install dynesty>=2.0.0 dill>=0.3.0 h5py>=3.0.0 emcee>=3.1.0
```

### Étape 5 : Configuration FSPS (Optionnel ou Stubs)

#### Option A : Installation de FSPS (si CMake et gfortran sont disponibles)

```cmd
# Cloner python-fsps
git clone https://github.com/dfm/python-fsps.git
cd python-fsps

# CRUCIAL: Initialiser les sous-modules Git
git submodule update --init --recursive

# Installer
set FC=gfortran
pip install . --no-cache-dir
```

#### Option B : Création de fichiers stub (si FSPS n'est pas installé)

Les fichiers stub permettent à Prospector de s'importer même sans FSPS complet. Voir section [Format du Fichier Stub FSPS](#6-format-du-fichier-stub-fsps).

### Étape 6 : Installation de Prospector depuis GitHub

```cmd
pip install git+https://github.com/bd-j/prospector.git --no-cache-dir
```

**⚠️ IMPORTANT** : Le package `prospector` sur PyPI est un **autre outil** (analyse de code Python). La version astronomique **DOIT** être installée depuis GitHub.

---

## 3. Problèmes Identifiés et Solutions

### Problème 1 : `ModuleNotFoundError: No module named 'sedpy.observate'`

**Symptôme** :
```
ModuleNotFoundError: No module named 'sedpy.observate'
```

**Cause** : Installation de `sedpy` depuis PyPI au lieu de GitHub.

**Solution** :
```cmd
pip uninstall sedpy -y
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir
python -c "from sedpy import observate; print('OK')"
```

**Vérification** : `sedpy.observate` doit être importable.

### Problème 2 : `ModuleNotFoundError: No module named 'astropy'`

**Symptôme** :
```
ModuleNotFoundError: No module named 'astropy'
```

**Cause** : Dépendance non explicitement listée dans certaines installations.

**Solution** :
```cmd
pip install astropy>=5.0.0
```

**Prévention** : Toujours installer `astropy>=5.0.0` **avant** Prospector.

### Problème 3 : FSPS - Sous-modules Git non initialisés

**Symptôme** :
```
CMake Error at src/fsps/CMakeLists.txt:2 (message):
  The source code in the FSPS submodule was not found.
  Please run 'git submodule update --init' to initialize the submodule.
```

**Cause** : `git submodule update --init --recursive` non exécuté après le clonage de `python-fsps`.

**Solution** :
```cmd
cd python-fsps
git submodule update --init --recursive
```

**Prévention** : **TOUJOURS** exécuter cette commande après avoir cloné `python-fsps`.

### Problème 4 : Format du fichier stub FSPS incorrect

**Symptôme** :
```
ValueError: could not assign tuple of length 1 to structure with 10 fields.
```

**Cause** : Le fichier stub `Nenkova08_y010_torusg_n10_q2.0.dat` a un format incorrect :
- Séparateur incorrect (1 espace au lieu de 3 espaces)
- Nombre de colonnes incorrect
- Nombre de lignes insuffisant

**Solution** : Créer le fichier stub avec le format exact :
- **4 lignes d'en-tête** (commentaires)
- **125 lignes de données**
- **10 colonnes** : 1 numérique (`wave`) + 9 chaînes (`fnu_*`)
- **Séparateur** : exactement **3 espaces** (`'   '`)

Voir section [Format du Fichier Stub FSPS](#6-format-du-fichier-stub-fsps) pour le format exact.

### Problème 5 : `SpecModel` vs `SedModel`

**Symptôme** :
```
ImportError: cannot import name 'SedModel' from 'prospect.models'
```

**Cause** : Prospector utilise `SpecModel` dans les versions récentes (commit 49ef5a17e et suivants).

**Solution** :
```python
# Utiliser SpecModel (version récente)
from prospect.models import SpecModel

# SedModel est obsolète dans les versions récentes
```

**Détection** :
```cmd
python -c "from prospect.models import SpecModel; print('SpecModel disponible')"
```

### Problème 6 : `RuntimeError: You need to have the SPS_HOME environment variable`

**Symptôme** :
```
RuntimeError: You need to have the SPS_HOME environment variable
```

**Cause** : Variable d'environnement `SPS_HOME` non définie avant l'import de `prospect`.

**Solution** :
```cmd
# Définir SPS_HOME
setx SPS_HOME "%USERPROFILE%\.local\share\fsps"

# Ou dans PowerShell:
[System.Environment]::SetEnvironmentVariable("SPS_HOME", "$env:USERPROFILE\.local\share\fsps", "User")
```

**Prévention** : Définir `SPS_HOME` **AVANT** d'importer `prospect` ou d'installer Prospector.

### Problème 7 : `TypeError: expected str, bytes or os.PathLike object, not NoneType`

**Symptôme** :
```
TypeError: expected str, bytes or os.PathLike object, not NoneType
```

**Cause** : `SPS_HOME` est `None` lors de l'import de `prospect`, même si défini plus tard.

**Solution** : S'assurer que `SPS_HOME` est défini **dès le début du script d'installation**, avant toute tentative d'import de `prospect`.

---

## 4. Scripts d'Installation

### Script PowerShell (Recommandé)

**Fichier** : `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1`

**Utilisation** :
```powershell
# Installation de base (sans FSPS, avec stubs)
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1

# Installation avec FSPS (nécessite CMake et gfortran)
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 -InstallFSPS

# Spécifier un environnement Conda différent
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 -CondaEnv monenv

# Sans vérification finale (pour tests)
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 -SkipVerification

# Forcer la recréation des fichiers stub
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 -ForceReinstall
```

**Avantages** :
- Gestion d'erreurs robuste avec `try-catch`
- Vérifications détaillées à chaque étape
- Messages colorés pour une meilleure lisibilité
- Création automatique des fichiers stub avec format correct
- Vérification du format des fichiers stub existants

**Fonctionnalités** :
- Vérification automatique des prérequis (Conda, Git)
- Installation de toutes les dépendances dans le bon ordre
- Création automatique de `SPS_HOME` et des fichiers stub
- Installation optionnelle de FSPS (si CMake et gfortran disponibles)
- Vérification post-installation complète

### Script Batch (Alternative)

**Fichier** : `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat`

**Utilisation** :
```cmd
INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat
```

**Avantages** :
- Double-clic simple, pas besoin d’ajuster la politique d’exécution PowerShell à la main
- Passe les arguments au `.ps1` (ex. `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat -InstallFSPS`)

### Installation sous WSL2 (alternative pour FSPS)

**Fichiers** : `Installation_fsps/prospector.bat` et `Installation_fsps/install_prospector_wsl.sh` (à conserver **ensemble** dans le même dossier). Certaines livraisons placent la paire à la **racine** du ZIP : même principe.

**Usage** : après installation de **WSL2** + distro **Ubuntu** (ou équivalent), exécutez `prospector.bat` sous Windows. Le script Linux installe **gfortran**, **cmake**, un **venv** dédié sous `~/.local/share/npoap-prospector-wsl/`, **sedpy** depuis GitHub, puis **Prospector** (avec repli `--no-deps` + **fsps** si besoin).

**Rôle** : environnement **Linux** séparé ; utile pour compiler ou utiliser **fsps** plus facilement qu’en natif Windows. **NPOAP** lancé depuis Windows continue d’utiliser **`astroenv`** ; la voie WSL ne remplace pas `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.*` pour intégrer Prospector au même Python que l’interface graphique.

---

## 5. Installation Manuelle

Si les scripts automatisés échouent ou si vous préférez installer manuellement, suivez ces étapes dans l'ordre exact :

### Étape 1 : Activer l'environnement Conda

```cmd
conda activate astroenv
```

### Étape 2 : Installer les dépendances de base

```cmd
pip install numpy>=1.20.0 scipy>=1.7.0 pandas>=1.3.0 astropy>=5.0.0
```

**⚠️ Ne pas omettre `astropy`** : C'est une dépendance cruciale de Prospector.

### Étape 3 : Installer sedpy depuis GitHub

```cmd
# Vérifier si sedpy est installé
python -c "import sedpy" 2>nul
if %ERRORLEVEL% EQU 0 (
    echo Desinstallation de sedpy (PyPI)...
    pip uninstall sedpy -y
)

# Installer depuis GitHub
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir

# Vérifier que sedpy.observate est disponible
python -c "from sedpy import observate; print('OK')"
```

**Si la vérification échoue** : Réinstaller depuis GitHub.

### Étape 4 : Installer les autres dépendances

```cmd
pip install dynesty>=2.0.0 dill>=0.3.0 h5py>=3.0.0 emcee>=3.1.0
```

### Étape 5 : Configurer SPS_HOME et créer les fichiers stub

#### Définir SPS_HOME

```cmd
# Dans CMD
setx SPS_HOME "%USERPROFILE%\.local\share\fsps"

# Dans PowerShell
[System.Environment]::SetEnvironmentVariable("SPS_HOME", "$env:USERPROFILE\.local\share\fsps", "User")
```

**Important** : Redémarrer le terminal après avoir défini `SPS_HOME` avec `setx`.

#### Créer les répertoires

```cmd
mkdir "%USERPROFILE%\.local\share\fsps\dust"
mkdir "%USERPROFILE%\.local\share\fsps\sed"
```

#### Créer le fichier stub (Format exact)

Le fichier stub doit avoir le format suivant (voir section [Format du Fichier Stub FSPS](#6-format-du-fichier-stub-fsps) pour le contenu exact) :

**Fichier** : `%USERPROFILE%\.local\share\fsps\dust\Nenkova08_y010_torusg_n10_q2.0.dat`

**Format** :
- 4 lignes d'en-tête (commentaires, ignorées par `skip_header=4`)
- 125 lignes de données
- 10 colonnes : 1 colonne numérique (`wave`) + 9 colonnes de chaînes (`fnu_*`)
- Séparateur : exactement **3 espaces** (`'   '`)

**Création avec Python** (méthode recommandée) :

```python
import pathlib

sps_home = pathlib.Path.home() / '.local' / 'share' / 'fsps'
dust_dir = sps_home / 'dust'
dust_dir.mkdir(parents=True, exist_ok=True)

dust_file = dust_dir / 'Nenkova08_y010_torusg_n10_q2.0.dat'

with open(dust_file, 'w') as f:
    # 4 lignes d'en-tête
    f.write("# Nenkova08 AGN torus dust model - Stub file\n")
    f.write("# This is a stub file created automatically\n")
    f.write("# Replace with real FSPS data file for full functionality\n")
    f.write("# wave   fnu_5   fnu_10   fnu_20   fnu_30   fnu_40   fnu_60   fnu_80   fnu_100   fnu_150\n")
    
    # 125 lignes de données
    for i in range(125):
        wave = 1.0 + i * 0.1
        fnu_values = [f"{(j+1)*0.001 + i*0.0001:.6f}" for j in range(9)]
        line = f"{wave:.6f}   {'   '.join(fnu_values)}\n"
        f.write(line)

print(f"Fichier stub créé: {dust_file}")
```

### Étape 6 : Installer Prospector depuis GitHub

```cmd
# S'assurer que SPS_HOME est défini dans cette session
set SPS_HOME=%USERPROFILE%\.local\share\fsps

# Installer Prospector
pip install git+https://github.com/bd-j/prospector.git --no-cache-dir
```

### Étape 7 : Vérifier l'installation

```cmd
# Vérifier sedpy
python -c "from sedpy import observate; print('sedpy.observate OK')"

# Vérifier Prospector
python -c "import prospect; print('prospect OK')"

# Vérifier SpecModel
python -c "from prospect.models import SpecModel; print('SpecModel OK')"

# Vérifier FastStepBasis
python -c "from prospect.sources import FastStepBasis; print('FastStepBasis OK')"

# Vérifier fit_model
python -c "from prospect.fitting import fit_model; print('fit_model OK')"
```

---

## 6. Format du Fichier Stub FSPS

### Fichier : `Nenkova08_y010_torusg_n10_q2.0.dat`

### Format Attendu par Prospector

Prospector lit ce fichier avec `np.genfromtxt` avec les paramètres suivants :

```python
dtype=[('wave', 'f8'),          # Colonne 1: longueur d'onde (float64)
       ('fnu_5', '<U20'),       # Colonnes 2-10: valeurs fnu (strings, max 20 chars)
       ('fnu_10', '<U20'),
       ('fnu_20', '<U20'),
       ('fnu_30', '<U20'),
       ('fnu_40', '<U20'),
       ('fnu_60', '<U20'),
       ('fnu_80', '<U20'),
       ('fnu_100', '<U20'),
       ('fnu_150', '<U20')]
delimiter='   '  # EXACTEMENT 3 espaces (pas 1, pas 2, pas de tabulation)
skip_header=4    # Ignore les 4 premières lignes
```

### Structure du Fichier

```
# Ligne 1: Commentaire (ignorée)
# Ligne 2: Commentaire (ignorée)
# Ligne 3: Commentaire (ignorée)
# Ligne 4: Commentaire avec noms de colonnes (ignorée)
# Ligne 5: Première ligne de données (wave fnu_5 fnu_10 ...)
# Ligne 6: Deuxième ligne de données
# ...
# Ligne 129: Dernière ligne de données (125 lignes au total)
```

### Exemple de Ligne de Données

```
1.000000   0.001000   0.002000   0.003000   0.004000   0.005000   0.006000   0.007000   0.008000   0.009000
```

**Format** :
- Colonne 1 (`wave`) : Nombre décimal (ex: `1.000000`)
- Colonnes 2-10 (`fnu_*`) : Chaînes de caractères (ex: `0.001000`)
- **Séparateur** : Exactement **3 espaces** entre chaque colonne

### Vérification du Format

Pour vérifier que votre fichier stub a le bon format :

```python
import numpy as np

dust_file = r"C:\Users\VotreNom\.local\share\fsps\dust\Nenkova08_y010_torusg_n10_q2.0.dat"

# Essayer de lire avec les mêmes paramètres que Prospector
data = np.genfromtxt(
    dust_file,
    dtype=[('wave', 'f8'),
           ('fnu_5', '<U20'), ('fnu_10', '<U20'), ('fnu_20', '<U20'),
           ('fnu_30', '<U20'), ('fnu_40', '<U20'), ('fnu_60', '<U20'),
           ('fnu_80', '<U20'), ('fnu_100', '<U20'), ('fnu_150', '<U20')],
    delimiter='   ',  # 3 espaces
    skip_header=4
)

print(f"Format OK: {len(data)} lignes lues avec succès")
print(f"Exemple: wave={data[0]['wave']}, fnu_5={data[0]['fnu_5']}")
```

Si cette vérification réussit, le fichier stub a le bon format.

### Erreurs Courantes

❌ **Séparateur incorrect** : 1 espace au lieu de 3
```
1.000000 0.001000 0.002000 ...  # INCORRECT
```

✅ **Séparateur correct** : 3 espaces exactement
```
1.000000   0.001000   0.002000   ...  # CORRECT
```

❌ **Nombre de colonnes incorrect** : Moins ou plus de 10 colonnes
```
1.000000   0.001000   0.002000  # INCORRECT (seulement 3 colonnes)
```

✅ **Nombre de colonnes correct** : Exactement 10 colonnes
```
1.000000   0.001000   0.002000   0.003000   0.004000   0.005000   0.006000   0.007000   0.008000   0.009000  # CORRECT
```

❌ **Nombre de lignes insuffisant** : Moins de 129 lignes au total
```
# Fichier avec seulement 50 lignes de données  # INCORRECT
```

✅ **Nombre de lignes correct** : 129 lignes au total (4 en-têtes + 125 données)
```
# 4 lignes d'en-tête
# + 125 lignes de données = 129 lignes total  # CORRECT
```

---

## 7. Vérification Post-Installation

### Test Complet

Créez un script de test `test_prospector_installation.py` :

```python
#!/usr/bin/env python
"""Script de vérification complète de l'installation Prospector"""

print("=" * 60)
print("Vérification de l'installation Prospector")
print("=" * 60)
print()

errors = []
warnings = []

# Test 1: sedpy
print("[1/7] Vérification de sedpy...")
try:
    from sedpy import observate
    print("  ✓ sedpy.observate disponible")
except ImportError as e:
    print(f"  ✗ sedpy.observate non disponible: {e}")
    errors.append("sedpy.observate")

# Test 2: astropy
print("[2/7] Vérification de astropy...")
try:
    import astropy
    print(f"  ✓ astropy {astropy.__version__}")
except ImportError as e:
    print(f"  ✗ astropy non disponible: {e}")
    errors.append("astropy")

# Test 3: FSPS (optionnel)
print("[3/7] Vérification de FSPS...")
try:
    import fsps
    print(f"  ✓ FSPS {fsps.__version__}")
except ImportError:
    print("  ⚠ FSPS non installé (utilisation de fichiers stub)")
    warnings.append("FSPS non installé")

# Test 4: Prospector (prospect)
print("[4/7] Vérification de prospect...")
try:
    import prospect
    version = getattr(prospect, '__version__', 'version inconnue')
    print(f"  ✓ prospect {version}")
except ImportError as e:
    print(f"  ✗ prospect non disponible: {e}")
    errors.append("prospect")
    print("\n  STOP: Impossible de continuer sans prospect")
    exit(1)

# Test 5: SpecModel
print("[5/7] Vérification de SpecModel...")
try:
    from prospect.models import SpecModel
    print("  ✓ SpecModel disponible")
except ImportError as e:
    print(f"  ⚠ SpecModel non disponible: {e}")
    warnings.append("SpecModel non disponible")

# Test 6: FastStepBasis
print("[6/7] Vérification de FastStepBasis...")
try:
    from prospect.sources import FastStepBasis
    print("  ✓ FastStepBasis disponible")
except ImportError as e:
    print(f"  ⚠ FastStepBasis non disponible: {e}")
    warnings.append("FastStepBasis non disponible (FSPS peut être requis)")

# Test 7: fit_model
print("[7/7] Vérification de fit_model...")
try:
    from prospect.fitting import fit_model
    print("  ✓ fit_model disponible")
except ImportError as e:
    print(f"  ⚠ fit_model non disponible: {e}")
    warnings.append("fit_model non disponible")

# Test 8: SPS_HOME
print("[8/8] Vérification de SPS_HOME...")
import os
sps_home = os.environ.get('SPS_HOME')
if sps_home:
    print(f"  ✓ SPS_HOME défini: {sps_home}")
    # Vérifier que le fichier stub existe
    from pathlib import Path
    stub_file = Path(sps_home) / 'dust' / 'Nenkova08_y010_torusg_n10_q2.0.dat'
    if stub_file.exists():
        print(f"  ✓ Fichier stub existe: {stub_file}")
    else:
        print(f"  ⚠ Fichier stub manquant: {stub_file}")
        warnings.append("Fichier stub FSPS manquant")
else:
    print("  ⚠ SPS_HOME non défini (sera créé automatiquement si nécessaire)")
    warnings.append("SPS_HOME non défini")

# Résumé
print()
print("=" * 60)
print("Résumé")
print("=" * 60)

if errors:
    print(f"✗ Erreurs ({len(errors)}):")
    for err in errors:
        print(f"  - {err}")
    print()
    print("L'installation n'est pas complète. Veuillez corriger ces erreurs.")
    exit(1)
else:
    print("✓ Aucune erreur critique")

if warnings:
    print(f"⚠ Avertissements ({len(warnings)}):")
    for warn in warnings:
        print(f"  - {warn}")
    print()
    print("Ces avertissements n'empêchent pas l'utilisation de Prospector,")
    print("mais certaines fonctionnalités peuvent être limitées.")
else:
    print("✓ Aucun avertissement")

print()
print("=" * 60)
print("✓ Installation Prospector vérifiée avec succès!")
print("=" * 60)
```

**Exécution** :
```cmd
conda activate astroenv
python test_prospector_installation.py
```

### Tests Individuels

#### Test sedpy
```cmd
python -c "from sedpy import observate; print('sedpy.observate OK')"
```

#### Test Prospector
```cmd
python -c "import prospect; print('prospect OK')"
```

#### Test SpecModel
```cmd
python -c "from prospect.models import SpecModel; print('SpecModel OK')"
```

#### Test FastStepBasis
```cmd
python -c "from prospect.sources import FastStepBasis; print('FastStepBasis OK')"
```

#### Test fit_model
```cmd
python -c "from prospect.fitting import fit_model; print('fit_model OK')"
```

#### Test FSPS (si installé)
```cmd
python -c "import fsps; print(f'FSPS {fsps.__version__} OK')"
```

---

## 8. Intégration dans l'Installation Globale NPOAP

Le script `installation.bat` inclut maintenant une étape optionnelle pour installer Prospector. Cette étape :

1. Demande confirmation à l'utilisateur
2. Vérifie la présence de Git
3. Appelle `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat` si disponible
4. Sinon, effectue une installation manuelle

### Modification dans installation.bat

```batch
REM ===================================================================
REM ETAPE 7: Installation de Prospector (optionnel)
REM ===================================================================
echo %BLUE%=== ETAPE 7: Installation de Prospector (optionnel) ===%RESET%
...
```

L'utilisateur peut choisir d'installer Prospector lors de l'installation initiale de NPOAP, ou l'installer plus tard avec le script dédié.

---

## 9. Dépannage

### Erreur : `ModuleNotFoundError: No module named 'sedpy.observate'`

**Solution** :
```cmd
pip uninstall sedpy -y
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir
python -c "from sedpy import observate; print('OK')"
```

### Erreur : `ModuleNotFoundError: No module named 'astropy'`

**Solution** :
```cmd
pip install astropy>=5.0.0
```

### Erreur : `ValueError: could not assign tuple of length 1 to structure with 10 fields`

**Cause** : Format du fichier stub FSPS incorrect.

**Solution** :
1. Supprimer l'ancien fichier stub :
   ```cmd
   del "%USERPROFILE%\.local\share\fsps\dust\Nenkova08_y010_torusg_n10_q2.0.dat"
   ```
2. Relancer le script d'installation ou recréer le fichier avec le format correct (voir section 6).

### Erreur : `ImportError: cannot import name 'SedModel'`

**Solution** : Utiliser `SpecModel` au lieu de `SedModel` (versions récentes de Prospector) :
```python
from prospect.models import SpecModel  # Correct
# from prospect.models import SedModel  # Obsolète
```

### Erreur : `RuntimeError: You need to have the SPS_HOME environment variable`

**Solution** :
```cmd
# Définir SPS_HOME
setx SPS_HOME "%USERPROFILE%\.local\share\fsps"

# Redémarrer le terminal, puis vérifier
echo %SPS_HOME%
```

**Note** : Si vous utilisez `setx`, redémarrez le terminal. Pour la session actuelle, utilisez `set` :
```cmd
set SPS_HOME=%USERPROFILE%\.local\share\fsps
```

### Erreur : Prospector affiche "non-installé" après installation et redémarrage

**Causes possibles** :
1. Fichier stub FSPS avec format incorrect
2. `SPS_HOME` non défini dans le profil utilisateur
3. Erreur lors de l'import de `prospect` (vérifier les logs)

**Solution** :
1. Vérifier les logs dans `logs/npoap_*.log`
2. Vérifier que `SPS_HOME` est défini : `echo %SPS_HOME%`
3. Vérifier le format du fichier stub (voir section 6)
4. Relancer le script d'installation avec `-ForceReinstall` pour recréer les stubs

### Erreur : Installation de FSPS échoue

**Solutions** :
1. **Utiliser WSL** (recommandé) : FSPS s'installe plus facilement sous Linux
2. **Installer CMake et gfortran** : Voir `docs/INSTALLATION_CMAKE.md` et `docs/INSTALLATION_GFORTRAN.md`
3. **Utiliser les fichiers stub** : Prospector fonctionne sans FSPS complet

---

## 10. Installation de FSPS (Optionnel)

### Prérequis

- **CMake** : Version 3.10 ou supérieure
- **gfortran** : Compilateur Fortran (TDM-GCC recommandé sur Windows)
- **Git** : Pour cloner et initialiser les sous-modules

### Installation sur Windows (Complexe)

L'installation native de FSPS sur Windows est **complexe** et nécessite :
1. Installation de CMake
2. Installation de gfortran (TDM-GCC, MSYS2, ou via Conda)
3. Configuration des variables d'environnement
4. Compilation du code Fortran

**Pour une installation plus simple, utilisez WSL** (voir section ci-dessous).

### Installation via WSL (Recommandé)

L'installation de FSPS est **beaucoup plus simple** sous Linux/WSL :

```bash
# Dans WSL
conda activate astroenv

# Installer les outils de build
sudo apt-get update
sudo apt-get install -y build-essential gfortran cmake git

# Cloner python-fsps
git clone https://github.com/dfm/python-fsps.git
cd python-fsps

# CRUCIAL: Initialiser les sous-modules Git
git submodule update --init --recursive

# Installer FSPS
export FC=gfortran
pip install . --no-cache-dir
```

**Vérification** :
```bash
python -c "import fsps; print(f'FSPS {fsps.__version__} OK')"
```

### Installation Native Windows (Avancé)

Si vous devez absolument installer FSPS nativement sur Windows :

1. **Installer CMake** : Voir `docs/INSTALLATION_CMAKE.md`
2. **Installer gfortran** : Voir `docs/INSTALLATION_GFORTRAN.md` ou `docs/INSTALLATION_GFORTRAN_RAPIDE.md`
3. **Configurer les variables d'environnement** :
   ```cmd
   set CMAKE_Fortran_COMPILER=C:\TDM-GCC-64\bin\gfortran.exe
   set FC=C:\TDM-GCC-64\bin\gfortran.exe
   ```
4. **Installer FSPS** :
   ```cmd
   git clone https://github.com/dfm/python-fsps.git
   cd python-fsps
   git submodule update --init --recursive
   pip install . --no-cache-dir
   ```

**Note** : Même avec CMake et gfortran configurés, l'installation peut échouer sur Windows. WSL reste la méthode la plus fiable.

---

## 11. Références

### Repositories GitHub

- **Prospector** : https://github.com/bd-j/prospector
- **sedpy** : https://github.com/bd-j/sedpy
- **FSPS (python-fsps)** : https://github.com/dfm/python-fsps
- **FSPS (sources Fortran)** : https://github.com/cconroy20/fsps (sous-module de python-fsps)

### Documentation

- **Prospector Documentation** : https://prospect.readthedocs.io/
- **FSPS Documentation** : https://fsps.readthedocs.io/
- **sedpy** : https://github.com/bd-j/sedpy (pas de documentation en ligne spécifique)

### Articles et Références

- **Johnson et al. (2021)** : "Prospector: A modular code for inferring stellar population parameters" (article principal)
- **Conroy et al. (2009)** : "FSPS: Flexible Stellar Population Synthesis" (article principal sur FSPS)

### Scripts et Outils

- **Script PowerShell** : `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1`
- **Script Batch** : `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat`
- **Script de test** : `test_prospector_installation.py` (à créer)

### Support

Pour toute question ou problème :
1. Consultez les logs dans le dossier `logs/`
2. Vérifiez que toutes les étapes de ce protocole ont été suivies
3. Vérifiez les références GitHub pour les dernières mises à jour

---

## Résumé des Commandes Essentielles

### Installation Complète (Script Automatique)

```powershell
# Méthode recommandée
.\INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1
```

### Installation Manuelle (Ordre Critique)

```cmd
conda activate astroenv

# 1. Dépendances de base
pip install numpy>=1.20.0 scipy>=1.7.0 pandas>=1.3.0 astropy>=5.0.0

# 2. sedpy depuis GitHub
pip uninstall sedpy -y
pip install git+https://github.com/bd-j/sedpy.git --no-cache-dir

# 3. Autres dépendances
pip install dynesty>=2.0.0 dill>=0.3.0 h5py>=3.0.0 emcee>=3.1.0

# 4. Configurer SPS_HOME (optionnel si FSPS non installé)
setx SPS_HOME "%USERPROFILE%\.local\share\fsps"
# Créer les fichiers stub (voir section 6)

# 5. Prospector depuis GitHub
pip install git+https://github.com/bd-j/prospector.git --no-cache-dir

# 6. Vérification
python -c "import prospect; from prospect.models import SpecModel; print('OK')"
```

### Vérification Rapide

```cmd
python -c "from sedpy import observate; import prospect; from prospect.models import SpecModel; print('Prospector OK!')"
```

---

**NPOAP - Protocole d'Installation Prospector Windows**

*Dernière mise à jour : Basé sur l'analyse des logs d'installation WSL et des erreurs rencontrées*
