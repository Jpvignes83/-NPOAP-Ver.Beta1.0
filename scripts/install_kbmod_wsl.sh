#!/usr/bin/env bash
# Installation automatique KBMOD sous WSL (appele par install_kbmod_wsl.bat).
# KBMOD_TARGET_PYTHON : defaut (vide ou "astroenv") = uniquement .../envs/astroenv/bin/python3
#                       (miniconda3, mambaforge, miniforge3, anaconda3).
#                       "system" = Python du PATH (ex. /usr/bin/python3), non recommande NPOAP.
#                       chemin absolu = interpreter explicite.
# KBMOD_PURGE_USER_KBMOD=1 : apres install conda, desinstalle kbmod (et coverage) du site-user /usr/bin/python3 si present.
set -euo pipefail

KBMOD_REPO="${KBMOD_REPO:-https://github.com/dirac-institute/kbmod.git}"
KBMOD_DIR="${KBMOD_DIR:-$HOME/kbmod}"
# Rempli par resolve_target_python()
PY=""

die() {
  echo "[KBMOD WSL] ERREUR: $*" >&2
  exit 1
}

use_pip_user() {
  [[ "$PY" == /usr/bin/python* ]] || [[ "$PY" == /usr/local/bin/python* ]]
}

resolve_target_python() {
  local raw="${KBMOD_TARGET_PYTHON:-}"

  if [[ "$raw" == "system" ]]; then
    if ! command -v python3 >/dev/null 2>&1; then
      die "python3 introuvable dans PATH (mot-cle system)."
    fi
    PY="$(command -v python3)"
    echo "[KBMOD WSL] Python systeme (mot-cle system) : $PY"
    return 0
  fi

  if [[ -n "$raw" && "$raw" != "astroenv" ]]; then
    if [[ ! -x "$raw" ]]; then
      die "KBMOD_TARGET_PYTHON n'est pas executable : $raw"
    fi
    PY="$raw"
    echo "[KBMOD WSL] Python cible (chemin explicite) : $PY"
    return 0
  fi

  for candidate in \
    "$HOME/miniconda3/envs/astroenv/bin/python3" \
    "$HOME/mambaforge/envs/astroenv/bin/python3" \
    "$HOME/miniforge3/envs/astroenv/bin/python3" \
    "$HOME/anaconda3/envs/astroenv/bin/python3"; do
    if [[ -x "$candidate" ]]; then
      PY="$candidate"
      echo "[KBMOD WSL] Environnement conda astroenv (defaut NPOAP) : $PY"
      return 0
    fi
  done
  die "conda env 'astroenv' introuvable. Creez-le (ex: conda create -n astroenv python=3.12) ou passez le chemin : KBMOD_TARGET_PYTHON=/home/VOUS/miniconda3/envs/astroenv/bin/python3. Pour forcer le Python systeme : KBMOD_TARGET_PYTHON=system"
}

ensure_apt() {
  if [[ "${SKIP_APT:-0}" == "1" ]]; then
    echo "[KBMOD WSL] SKIP_APT=1 - paquets APT ignores."
    return 0
  fi
  echo "[KBMOD WSL] Installation des paquets APT (sudo peut demander votre mot de passe)..."
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential cmake python3-dev python3-pip python3-venv git
}

clone_or_update() {
  if [[ -d "$KBMOD_DIR/.git" ]]; then
    echo "[KBMOD WSL] Depot existant: $KBMOD_DIR - mise a jour..."
    git -C "$KBMOD_DIR" fetch --tags origin || true
    git -C "$KBMOD_DIR" pull --ff-only || echo "[KBMOD WSL] Attention: git pull a echoue (branche locale?). Continue avec les sous-modules."
    git -C "$KBMOD_DIR" submodule update --init --recursive
  else
    echo "[KBMOD WSL] Clone recursive dans $KBMOD_DIR ..."
    mkdir -p "$(dirname "$KBMOD_DIR")"
    git clone --recursive "$KBMOD_REPO" "$KBMOD_DIR"
  fi
}

# KBMOD compile une extension C++/CMake ; sans ces variables, CMake peut choisir /usr/bin/python3
# alors que pip utilise conda -> wheel .cpython-312.so au lieu de .cpython-310.so (echec copie).
kbmod_with_cmake_python() (
  export PYTHON_EXECUTABLE="$PY"
  export Python3_EXECUTABLE="$PY"
  _pyroot="$(dirname "$(dirname "$PY")")"
  export CMAKE_ARGS="-DPython3_EXECUTABLE=$PY -DPYTHON_EXECUTABLE=$PY -DPython3_ROOT_DIR=$_pyroot"
  "$@"
)

pip_install_editable() {
  echo "[KBMOD WSL] Compilation / installation pip en mode editable avec : $PY"
  kbmod_with_cmake_python "$PY" -m pip install --upgrade pip setuptools wheel 2>/dev/null || true

  if use_pip_user; then
    local em
    em="$("$PY" -c "import sysconfig; from pathlib import Path; p = Path(sysconfig.get_path('stdlib')) / 'EXTERNALLY-MANAGED'; print(p if p.is_file() else '')")"
    local bsp=()
    if [[ -n "$em" ]]; then
      echo "[KBMOD WSL] PEP 668 detecte ($em) : pip avec --break-system-packages et --user."
      bsp=(--break-system-packages)
    fi
    if [[ "${#bsp[@]}" -eq 0 ]]; then
      if ! kbmod_with_cmake_python "$PY" -m pip install --user -e "$KBMOD_DIR"; then
        echo "[KBMOD WSL] Nouvelle tentative avec --break-system-packages."
        kbmod_with_cmake_python "$PY" -m pip install --break-system-packages --user -e "$KBMOD_DIR"
      fi
    else
      kbmod_with_cmake_python "$PY" -m pip install "${bsp[@]}" --user -e "$KBMOD_DIR"
    fi
  else
    echo "[KBMOD WSL] Installation dans l'environnement Python (conda/venv), sans --user."
    kbmod_with_cmake_python "$PY" -m pip install -e "$KBMOD_DIR"
  fi
}

# Numba 0.62+ charge coverage si present ; il faut coverage>=7.6.1 (coverage.types.Tracer).
# Voir https://github.com/numba/numba/issues/10239 - versions intermediaires cassent l'import.
pin_coverage_for_numba() {
  echo "[KBMOD WSL] Compatibilite numba : assurance coverage>=7.6.1 (API coverage.types.Tracer)."
  if use_pip_user; then
    local em
    em="$("$PY" -c "import sysconfig; from pathlib import Path; p = Path(sysconfig.get_path('stdlib')) / 'EXTERNALLY-MANAGED'; print(p if p.is_file() else '')")"
    local bsp=()
    [[ -n "$em" ]] && bsp=(--break-system-packages)
    "$PY" -m pip install "${bsp[@]}" --user "coverage>=7.6.1"
  else
    "$PY" -m pip install "coverage>=7.6.1"
  fi
}

ensure_local_bin_path() {
  export PATH="$HOME/.local/bin:$PATH"
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
  "$PY" -c "import kbmod.search; from kbmod.core.psf import PSF; print('KBMOD OK')"
}

maybe_purge_user_site_kbmod() {
  [[ "${KBMOD_PURGE_USER_KBMOD:-0}" == "1" ]] || return 0
  use_pip_user && return 0
  [[ -x /usr/bin/python3 ]] || return 0
  echo "[KBMOD WSL] KBMOD_PURGE_USER_KBMOD=1 : desinstallation eventuelle kbmod/coverage du site-user systeme..."
  /usr/bin/python3 -m pip uninstall -y kbmod 2>/dev/null || true
  /usr/bin/python3 -m pip uninstall -y coverage 2>/dev/null || true
}

main() {
  echo "[KBMOD WSL] Repertoire KBMOD: $KBMOD_DIR"
  command -v git >/dev/null || die "'git' introuvable. Lancez sans SKIP_APT ou installez git."
  resolve_target_python
  ensure_apt
  clone_or_update
  pip_install_editable
  pin_coverage_for_numba
  if use_pip_user; then
    ensure_local_bin_path
  fi
  maybe_purge_user_site_kbmod
  verify_import
  echo ""
  echo "[KBMOD WSL] Termine avec succes."
  echo "[KBMOD WSL] Dans config.py (NPOAP), mettez exactement :"
  echo "           KBMOD_WSL_PYTHON = \"$PY\""
  if ! use_pip_user; then
    echo "[KBMOD WSL] Pour retirer une ancienne KBMOD dans ~/.local (Python systeme) :"
    echo "           /usr/bin/python3 -m pip uninstall -y kbmod"
    echo "           Ou sous Windows avant le .bat : set KBMOD_PURGE_USER_KBMOD=1 puis relancer ce script."
  fi
}

main "$@"
