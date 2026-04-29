#!/usr/bin/env bash
# Installation automatique KBMOD sous WSL (appele par install_kbmod_wsl.bat).
set -euo pipefail

KBMOD_REPO="${KBMOD_REPO:-https://github.com/dirac-institute/kbmod.git}"
KBMOD_DIR="${KBMOD_DIR:-$HOME/kbmod}"

die() {
  echo "[KBMOD WSL] ERREUR: $*" >&2
  exit 1
}

ensure_apt() {
  if [[ "${SKIP_APT:-0}" == "1" ]]; then
    echo "[KBMOD WSL] SKIP_APT=1 — paquets APT ignores."
    return 0
  fi
  echo "[KBMOD WSL] Installation des paquets APT (sudo peut demander votre mot de passe)..."
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential cmake python3-dev python3-pip python3-venv git
}

clone_or_update() {
  if [[ -d "$KBMOD_DIR/.git" ]]; then
    echo "[KBMOD WSL] Depot existant: $KBMOD_DIR — mise a jour..."
    git -C "$KBMOD_DIR" fetch --tags origin || true
    git -C "$KBMOD_DIR" pull --ff-only || echo "[KBMOD WSL] Attention: git pull a echoue (branche locale?). Continue avec les sous-modules."
    git -C "$KBMOD_DIR" submodule update --init --recursive
  else
    echo "[KBMOD WSL] Clone recursive dans $KBMOD_DIR ..."
    mkdir -p "$(dirname "$KBMOD_DIR")"
    git clone --recursive "$KBMOD_REPO" "$KBMOD_DIR"
  fi
}

pip_install_editable() {
  echo "[KBMOD WSL] Compilation / installation pip en mode editable..."
  python3 -m pip install --user --upgrade pip setuptools wheel || true
  if python3 -m pip install --user -e "$KBMOD_DIR"; then
    return 0
  fi
  echo "[KBMOD WSL] Nouvelle tentative avec --break-system-packages (PEP 668 sur certaines distributions)."
  python3 -m pip install --break-system-packages --user -e "$KBMOD_DIR"
}

ensure_local_bin_path() {
  local marker='# NPOAP KBMOD PATH'
  if [[ -f "$HOME/.bashrc" ]] && grep -qF "$marker" "$HOME/.bashrc" 2>/dev/null; then
    return 0
  fi
  echo "[KBMOD WSL] Ajout de ~/.local/bin au PATH dans ~/.bashrc"
  {
    echo ""
    echo "$marker"
    echo 'export PATH="$HOME/.local/bin:$PATH"'
  } >> "$HOME/.bashrc"
}

verify_import() {
  echo "[KBMOD WSL] Verification import..."
  python3 -c "import kbmod.search; from kbmod.core.psf import PSF; print('KBMOD OK')"
}

main() {
  echo "[KBMOD WSL] Repertoire KBMOD: $KBMOD_DIR"
  command -v git >/dev/null || die "'git' introuvable. Lancez sans SKIP_APT ou installez git."
  ensure_apt
  clone_or_update
  pip_install_editable
  ensure_local_bin_path
  verify_import
  echo ""
  echo "[KBMOD WSL] Termine avec succes."
  echo "Python utilise: $(command -v python3)"
  echo "Si NPOAP ne trouve pas kbmod, definissez KBMOD_WSL_PYTHON dans config.py sur ce chemin (voir docs/INSTALL_KBMOD_WSL.md)."
}

main "$@"
