#!/usr/bin/env bash
# NPOAP — Installation Prospector + FSPS dans WSL (Ubuntu/Debian).
# Appelé par prospector.bat ; ne pas déplacer seul sans adapter le .bat.
set -euo pipefail

echo ""
echo "========================================"
echo " NPOAP — Prospector + FSPS (WSL/Linux)"
echo "========================================"
echo ""

if [[ -f /proc/version ]] && grep -qiE 'microsoft|wsl' /proc/version 2>/dev/null; then
  echo "Environnement WSL détecté."
else
  echo "Note: exécution hors WSL possible ; paquets système requis (gfortran, etc.)."
fi
echo ""

if ! command -v apt-get >/dev/null 2>&1; then
  echo "ERREUR: apt-get introuvable (distribution non Debian/Ubuntu)."
  echo "Installez manuellement : gfortran git build-essential python3 python3-venv python3-pip cmake pkg-config"
  exit 1
fi

echo "[1/4] Paquets système (sudo peut demander votre mot de passe)..."
sudo apt-get update -y
sudo apt-get install -y \
  gfortran \
  git \
  build-essential \
  python3 \
  python3-venv \
  python3-pip \
  cmake \
  pkg-config

echo ""
echo "[2/4] Environnement virtuel Python..."
VENV_ROOT="${NPOAP_PROSPECTOR_VENV:-$HOME/.local/share/npoap-prospector-wsl}"
mkdir -p "$VENV_ROOT"
python3 -m venv "$VENV_ROOT/venv"
# shellcheck source=/dev/null
source "$VENV_ROOT/venv/bin/activate"
pip install -U pip wheel setuptools

echo ""
echo "[3/4] Dépendances Python (sedpy depuis GitHub)..."
pip install "numpy>=1.20" "scipy>=1.7" "pandas>=1.3" "astropy>=5"
pip install "git+https://github.com/bd-j/sedpy.git"
pip install \
  "dynesty>=2.0.0" \
  "dill>=0.3.0" \
  "h5py>=3.0.0" \
  "emcee>=3.1.0" \
  "corner>=2.2.3" \
  "matplotlib>=3.8.4"

SPS_HOME="${SPS_HOME:-$HOME/.local/share/fsps}"
export SPS_HOME
mkdir -p "$SPS_HOME/dust" "$SPS_HOME/sed"

MARK="# NPOAP Prospector WSL (SPS_HOME)"
if [[ -f "$HOME/.bashrc" ]] && ! grep -qF "$MARK" "$HOME/.bashrc" 2>/dev/null; then
  {
    echo ""
    echo "$MARK"
    echo "export SPS_HOME=\"$SPS_HOME\""
  } >> "$HOME/.bashrc"
  echo "SPS_HOME ajouté à ~/.bashrc : $SPS_HOME"
else
  echo "SPS_HOME=$SPS_HOME (déjà configuré ou .bashrc absent — export manuel si besoin)"
fi

echo ""
echo "[4/4] Prospector (GitHub) + fsps..."
if pip install "git+https://github.com/bd-j/prospector.git" --no-cache-dir; then
  echo "Installation Prospector avec dépendances pip : OK"
else
  echo "Échec avec dépendances ; essai astro-prospector --no-deps puis fsps..."
  pip install "git+https://github.com/bd-j/prospector.git" --no-cache-dir --no-deps
  pip install "fsps" --no-cache-dir
fi

echo ""
python3 -c "from sedpy import observate; print('sedpy.observate : OK')"
python3 -c "import prospect; from prospect.models import SpecModel; print('prospect', getattr(prospect, '__version__', '?'), '- SpecModel : OK')"
if python3 -c "import fsps" 2>/dev/null; then
  python3 -c "import fsps; print('fsps', getattr(fsps, '__version__', '?'), ': OK')"
else
  echo "Attention : import fsps impossible (données SPS_HOME ou build). Prospector peut rester utilisable selon le modèle."
fi

echo ""
echo "========================================"
echo " Terminé."
echo "========================================"
echo ""
echo "Pour utiliser cet environnement dans WSL :"
echo "  source \"$VENV_ROOT/venv/bin/activate\""
echo ""
echo "Test rapide :"
echo "  python3 -c \"import prospect, fsps; print('OK')\""
echo ""
