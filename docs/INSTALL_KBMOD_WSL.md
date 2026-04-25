# Installation de KBMOD sous WSL (Windows Subsystem for Linux)

NPOAP sur Windows n'installe pas KBMOD (compilation difficile sous MSVC).  
La detection KBMOD est lancee via WSL avec `scripts/kbmod_wsl_detect.py`, qui ecrit les candidats dans `kbmod_candidates.csv`.

## Prerequis

- WSL 2 avec une distribution Linux (Ubuntu recommandee)
- GPU NVIDIA (optionnel mais recommande)
- CUDA sous WSL si vous utilisez le GPU

## Etapes d'installation

### 1) Ouvrir WSL

Depuis PowerShell/CMD:

```powershell
wsl
```

### 2) Mettre a jour Linux

```bash
sudo apt update
sudo apt upgrade -y
```

### 3) Installer les outils de build

```bash
sudo apt install -y build-essential cmake python3-dev python3-pip python3-venv git
```

### 4) (Optionnel) Installer CUDA sous WSL

Suivez la doc NVIDIA: <https://docs.nvidia.com/cuda/wsl-user-guide/>

Exemple Ubuntu:

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
sudo apt install -y cuda-toolkit
```

Si besoin, redemarrez WSL:

1. Quittez WSL: `exit`
2. Dans **PowerShell/CMD Windows (hors WSL)**:

```powershell
wsl --shutdown
```

3. Rouvrez WSL (`wsl` ou `wsl -d Ubuntu`)

### 5) Cloner KBMOD avec sous-modules

```bash
cd ~
git clone --recursive https://github.com/dirac-institute/kbmod.git
cd kbmod
pip3 install --user -e .
```

Si vous avez deja clone sans `--recursive`:

```bash
cd ~/kbmod
git submodule update --init --recursive
pip3 install --user -e .
```

### 6) (Optionnel) Ajouter ~/.local/bin au PATH

```bash
export PATH="$HOME/.local/bin:$PATH"
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 7) Verification

```bash
python3 -c "import kbmod.search; from kbmod.core.psf import PSF; print('KBMOD OK')"
```

## Utilisation depuis NPOAP

Dans l'onglet Asteroides:

1. Chargez un dossier FITS
2. Cliquez "Detection KBMOD (via WSL)"
3. Lancez la detection
4. NPOAP execute `wsl python3 .../scripts/kbmod_wsl_detect.py ...`
5. Le script genere `kbmod_candidates.csv`

## Depannage rapide

- `No module named 'kbmod'`: KBMOD n'est pas installe dans le `python3` utilise par WSL
- `No FITS files in directory`: verifier le chemin WSL (`/mnt/c/...`)
- Erreurs pybind11/CMake: relancer `git submodule update --init --recursive`
- Python 3.13: preferer Python 3.10/3.11

