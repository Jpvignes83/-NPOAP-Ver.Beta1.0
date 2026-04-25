@echo off
REM NPOAP — Réinstalle ou met à jour uniquement SORA (sora-astro) dans « astroenv ».
REM L’installation standard inclut sora-astro dans requirements_install_core.txt ; ce script met a jour SORA seul.
REM Utile en dépannage ou pour forcer une version de sora-astro sans réinstaller tout NPOAP.

chcp 65001 >nul
echo ========================================
echo Installation SORA (sora-astro) dans astroenv
echo ========================================
echo.

echo Activation de l'environnement astroenv...
call conda activate astroenv

if errorlevel 1 (
    echo ERREUR: Impossible d'activer l'environnement astroenv.
    echo Créez-le par exemple avec : conda create -n astroenv python=3.11 -y
    pause
    exit /b 1
)

echo.
echo Installation via pip (peut prendre plusieurs minutes ; dépendances réseau)...
python -m pip install --upgrade pip
python -m pip install "sora-astro"

if errorlevel 1 (
    echo.
    echo ERREUR lors de l'installation de sora-astro.
    echo En cas d'échec sur Cartopy / GEOS, consultez la doc SORA et les prérequis système.
    pause
    exit /b 1
)

echo.
echo ========================================
echo Installation terminee avec succes.
echo ========================================
echo.
echo sora-astro est a jour dans astroenv. Relancez NPOAP depuis cet environnement.
echo.
pause
