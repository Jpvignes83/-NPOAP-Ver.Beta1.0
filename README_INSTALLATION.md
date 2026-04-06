# Installation Automatique NPOAP - Windows

Ce guide décrit comment installer NPOAP automatiquement sur Windows 10/11.

## Prérequis

- Windows 10 ou Windows 11
- Connexion internet (pour télécharger les composants)
- Droits administrateur (obligatoires pour `installation.bat`)
- Environ 5 Go d'espace disque libre
- **Version Python cible** : Python **3.11.x** (l'installateur crée/attend un environnement Conda `astroenv` en 3.11)

## Installation rapide

1. **Téléchargez l'archive ZIP** de NPOAP
2. **Extrayez l'archive** dans un dossier de votre choix
3. **Clic droit sur `installation.bat`** → **"Exécuter en tant qu'administrateur"**
4. **Acceptez** l'accord d'utilisation, puis **indiquez un chemin d'installation absolu** (ex. `C:\Astro\NPOAP`). Aucun dossier par défaut n'est imposé.
5. **Suivez les étapes** affichées (Python, Miniconda, environnement `astroenv`, dépendances `requirements.txt`, copie des fichiers, etc.)
6. À la fin, la console liste les **scripts optionnels** ; le fichier **`LISTE_INSTALL_OPTIONNELS.txt`** (dans le même dossier que NPOAP) reprend la même liste.

## Ce que fait le script `installation.bat`

Le script automatise le **cœur** de l'installation sur Windows. Les outils lourds ou spécialisés (compilateur MSVC, WSL, Astrometry.net, KBMOD, CMake, Prospector) ne sont **plus** intégrés dans ce fichier : ils sont lancés à part, quand vous en avez besoin (voir section suivante).

### Étapes automatisées

1. **Accord d'utilisation** puis vérification des privilèges administrateur
2. **Sélection du dossier d'installation** (chemin absolu au choix de l'utilisateur)
3. **Vérifications système** : architecture (x64 recommandé) et **droits d'écriture** dans le dossier choisi
4. **Installation de Python 3.11** (facultatif si vous conservez une installation existante)
5. **Installation de Miniconda** (ou réutilisation de Conda déjà présent)
6. **Création ou réutilisation** de l'environnement Conda **`astroenv`**
7. **Installation des dépendances Python** depuis `requirements.txt` (via `pip`)
8. **Copie des fichiers NPOAP** vers le dossier d'installation
9. **Création de `LANCER_NPOAP.bat`** dans ce dossier
10. **Test** si `test_installation.py` est présent
11. **Protection en lecture seule** des fichiers sources (comportement inchangé si `protect_files.bat` est utilisé)

KBMOD n'est **pas** dans `requirements.txt` pour Windows ; pour le Synthetic Tracking sous Linux/WSL, voir `docs/INSTALL_KBMOD_WSL.md` et `install_kbmod_wsl.bat`.

## Composants optionnels (scripts séparés)

Ces actions se font **après** l'installation principale, depuis le dossier NPOAP (copie d'installation ou dossier du déploiement), en double-cliquant ou depuis `cmd` :

| Script | Usage |
|--------|--------|
| `install_msvc_build_tools.bat` | Ouvre la page **Visual C++ Build Tools** (MSVC). Utile si l'installation ou la compilation de paquets comme **PHOEBE2** échoue faute de toolchain C++. |
| `install_cmake.bat` | Aide à installer **CMake** (téléchargement officiel ou rappel `winget`). Requis pour certaines compilations (ex. KBMOD sous WSL/Linux). |
| `install_wsl.bat` | Lance `wsl --install` (peut demander un **redémarrage**). |
| `install_ubuntu_wsl.bat` | Ajoute la distribution **Ubuntu** dans WSL (`wsl --install -d Ubuntu`). |
| `install_astrometry_wsl.bat` | Installe **astrometry.net** dans la distro WSL par défaut (`apt-get`). |
| `install_kbmod_wsl.bat` | Ouvre la documentation **KBMOD** et les pistes d'installation sous WSL/Linux. |
| `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat` | Installe **Prospector** (astro-prospector) dans l’environnement Conda **`astroenv`** sous Windows : appelle le script PowerShell ci-dessous. |
| `INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1` | Même installation que le `.bat`, avec options (`-InstallFSPS`, `-CondaEnv`, etc.). Le `.bat` est un lanceur pratique ; le `.ps1` est le script réel. |
| `Installation_fsps/prospector.bat` | Lance sous **WSL2** le script Linux `Installation_fsps/install_prospector_wsl.sh` : **venv dédié** dans WSL, **gfortran**, tentative d’installation complète incluant **fsps** (souvent plus simple que sous Windows natif). **Indépendant** de Conda Windows : NPOAP continue d’utiliser `astroenv` pour l’interface ; cette voie sert surtout à un environnement Prospector/FSPS sous Linux. |

La liste courte est aussi dans **`LISTE_INSTALL_OPTIONNELS.txt`**.

### Prospector / FSPS : faut-il encore le `.bat` et le `.ps1` Windows ?

**Oui**, si vous voulez **Prospector dans le même Python que NPOAP** (`conda activate astroenv` → `python main.py`). C’est le cas le plus courant : exécutez **`INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat`** (ou le `.ps1` à la main). Sans **gfortran**, la compilation de **fsps** échoue souvent sous Windows : le script bascule alors sur **`pip install … --no-deps`** et des **fichiers stub** sous `SPS_HOME` ; Prospector reste utilisable pour une grande partie des fonctions.

Le dossier **`Installation_fsps/`** est une **alternative** pour un outillage **Prospector + FSPS** complet dans **WSL**, sans remplacer l’installateur Windows pour l’application NPOAP. Certaines livraisons ZIP placent `prospector.bat` et `install_prospector_wsl.sh` à la **racine** du dossier : même usage, en gardant les deux fichiers au même endroit.

Voir aussi **`docs/MANUEL_INSTALLATION.md`** (section Prospector) et **`docs/PROTOCOLE_INSTALLATION_PROSPECTOR_WINDOWS.md`**.

**Ordre usuel pour l'astrométrie locale** : `install_wsl.bat` → redémarrage si Windows le demande → `install_ubuntu_wsl.bat` → configuration du compte Ubuntu → `install_astrometry_wsl.bat`.

Sans WSL, vous pouvez utiliser **Astrometry.net en ligne** (NOVA, clé API) selon la configuration de NPOAP.

## Installation manuelle (si le script échoue)

Si le script d'installation automatique échoue, suivez **`docs/MANUEL_INSTALLATION.md`** pour une installation manuelle étape par étape.

## Après l'installation

### Lancer NPOAP

1. **Double-cliquez sur `LANCER_NPOAP.bat`** dans le dossier d'installation

OU

2. Ouvrez un terminal et exécutez :
   ```cmd
   call "chemin\vers\miniconda3\Scripts\activate.bat" astroenv
   cd C:\Chemin\Vers\Votre\NPOAP
   python main.py
   ```
   Remplacez les chemins par ceux indiqués à la fin de `installation.bat` (racine Conda et dossier NPOAP).

### Vérifier l'installation

Exécutez le script de test (si présent) :
```cmd
call "chemin\vers\miniconda3\Scripts\activate.bat" astroenv
cd C:\Chemin\Vers\Votre\NPOAP
python test_installation.py
```

## Installation d'Astrometry.net (WSL)

1. Assurez-vous que **WSL** et **Ubuntu** sont en place (`install_wsl.bat`, puis `install_ubuntu_wsl.bat` si besoin).
2. Exécutez **`install_astrometry_wsl.bat`** (une invite WSL peut demander le mot de passe `sudo`).

Ou manuellement dans WSL :

```bash
sudo apt-get update
sudo apt-get install -y astrometry.net
```

## Dépannage

### Le script demande toujours les droits administrateur

- Faites **clic droit** sur `installation.bat` → **"Exécuter en tant qu'administrateur"**.

### Python n'est pas détecté après installation

- Redémarrez le terminal ou PowerShell
- Vérifiez le PATH : `python --version`

### Conda n'est pas détecté après installation

- Fermez et rouvrez le terminal
- Ou : `conda init cmd.exe` puis nouveau terminal  
- Pour une session `cmd` sans init : `call "...\miniconda3\Scripts\activate.bat" astroenv` (comme dans `LANCER_NPOAP.bat`)

### Erreur lors de l'installation ou de l'utilisation de PHOEBE2

- Installez **Visual C++ Build Tools** : exécutez `install_msvc_build_tools.bat`, installez la charge utile **Desktop development with C++** (ou équivalent MSVC x64), puis réessayez dans l'environnement `astroenv` :  
  `python -m pip install --force-reinstall phoebe`  
  (ou réinstallez la dépendance concernée selon le message d'erreur).

### WSL demande un redémarrage

- Redémarrez l'ordinateur, puis poursuivez avec `install_ubuntu_wsl.bat` ou `install_astrometry_wsl.bat` selon votre objectif. Il n'est pas nécessaire de relancer `installation.bat` pour cela.

### Erreurs de dépendances Python (`pip`)

- Vérifiez la connexion internet
- Activez `astroenv`, puis : `python -m pip install -r requirements.txt`
- Consultez `docs/MANUEL_INSTALLATION.md` pour plus de détails

## Protection des fichiers

Les fichiers sources de NPOAP sont automatiquement protégés en lecture seule lors de l'installation pour éviter les modifications accidentelles.

### Scripts de protection

- **`protect_files.bat`** : Remet les fichiers en lecture seule
- **`unprotect_files.bat`** : Retire la protection (à utiliser avec précaution)

**Note** : Les fichiers de configuration (`config.json`) et les dossiers `logs/` et `output/` ne sont pas protégés et peuvent être modifiés librement.

**Important** : L'attribut « lecture seule » est une protection faible. Pour plus de détails, consultez `INFO_PROTECTION.md` si présent dans votre livraison.

## Support

Pour toute question ou problème, consultez :

- **`docs/MANUEL_INSTALLATION.md`** : guide d'installation détaillé
- **`docs/MANUEL_UTILISATEUR.md`** : guide d'utilisation
- Les logs dans le dossier `logs/`

---

**NPOAP - Installation automatique Windows**
