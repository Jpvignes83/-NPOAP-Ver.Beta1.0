@echo off
setlocal EnableDelayedExpansion
REM ============================================================
REM Script de construction de distribution NPOAP
REM Double-cliquez sur ce fichier pour construire une distribution
REM ============================================================

title Construction Distribution NPOAP

REM Aller dans le repertoire du script (racine NPOAP)
cd /d "%~dp0"

REM Verifier que nous sommes au bon endroit
if not exist "main.py" (
    echo ============================================================
    echo ERREUR: main.py non trouve.
    echo Assurez-vous d'etre dans le repertoire NPOAP.
    echo ============================================================
    pause
    exit /b 1
)

if not exist "build\build_distribution.py" (
    echo ============================================================
    echo ERREUR: build\build_distribution.py non trouve.
    echo Le systeme de build n'est pas complet.
    echo ============================================================
    pause
    exit /b 1
)

REM Afficher le menu (ASCII uniquement pour compatibilite CMD)
echo.
echo ============================================================
echo          CONSTRUCTION DE DISTRIBUTION NPOAP
echo ============================================================
echo.
echo Profils disponibles:
echo.
echo   1. Reduction - Accueil + Reduction
echo   2. Exoplanetes - Accueil + Reduction + Exoplanetes
echo   3. Asteroides - Accueil + Reduction + Photometrie Asteroides
echo   4. Etoiles Doubles - Accueil + Etoiles Binaires + ELI
echo   5. Transitoires - Accueil + Reduction + Photometrie Transitoires
echo   6. Catalogues - Accueil + extraction catalogues + Analyse des Donnees
echo   7. Complet - Tous les modules
echo   0. Annuler
echo.
echo ============================================================
echo.

REM Demander le choix
set /p CHOIX="Votre choix (1-7 ou 0 pour annuler): "

REM Traiter le choix
if "%CHOIX%"=="1" (
    set PROFILE=reduction
    set PROFILE_NAME=Reduction
) else if "%CHOIX%"=="2" (
    set PROFILE=exoplanets
    set PROFILE_NAME=Exoplanetes
) else if "%CHOIX%"=="3" (
    set PROFILE=asteroids
    set PROFILE_NAME=Asteroides
) else if "%CHOIX%"=="4" (
    set PROFILE=binary_stars
    set PROFILE_NAME=Etoiles Doubles
) else if "%CHOIX%"=="5" (
    set PROFILE=transient
    set PROFILE_NAME=Transitoires
) else if "%CHOIX%"=="6" (
    set PROFILE=catalogues
    set PROFILE_NAME=Catalogues
) else if "%CHOIX%"=="7" (
    set PROFILE=full
    set PROFILE_NAME=Complet
) else if "%CHOIX%"=="0" (
    echo.
    echo Construction annulee.
    pause
    exit /b 0
) else (
    echo.
    echo ============================================================
    echo ERREUR: Choix invalide. Choisissez 1, 2, 3, 4, 5, 6, 7 ou 0.
    echo ============================================================
    pause
    exit /b 1
)

echo.
echo ============================================================
echo Construction de la distribution: %PROFILE_NAME%
echo ============================================================
echo.

REM Chercher Python (priorite aux executables concrets Conda/Miniconda)
set "PYTHON_CMD="
if exist "%ProgramData%\Miniconda3\envs\astroenv\python.exe" set "PYTHON_CMD="%ProgramData%\Miniconda3\envs\astroenv\python.exe""
if not defined PYTHON_CMD if exist "%ProgramData%\Miniconda3\python.exe" set "PYTHON_CMD="%ProgramData%\Miniconda3\python.exe""
if not defined PYTHON_CMD if exist "%USERPROFILE%\miniconda3\python.exe" set "PYTHON_CMD="%USERPROFILE%\miniconda3\python.exe""
if not defined PYTHON_CMD if exist "%LOCALAPPDATA%\miniconda3\python.exe" set "PYTHON_CMD="%LOCALAPPDATA%\miniconda3\python.exe""
if not defined PYTHON_CMD if exist "%ProgramData%\Anaconda3\python.exe" set "PYTHON_CMD="%ProgramData%\Anaconda3\python.exe""
if not defined PYTHON_CMD if exist "%USERPROFILE%\Anaconda3\python.exe" set "PYTHON_CMD="%USERPROFILE%\Anaconda3\python.exe""

if not defined PYTHON_CMD (
    py -3 -V >nul 2>&1
    if %ERRORLEVEL% EQU 0 (
        set "PYTHON_CMD=py -3"
    ) else (
        python -V >nul 2>&1
        if %ERRORLEVEL% EQU 0 (
            set "PYTHON_CMD=python"
        ) else (
            python3 -V >nul 2>&1
            if %ERRORLEVEL% EQU 0 (
                set "PYTHON_CMD=python3"
            ) else (
                echo ============================================================
                echo ERREUR: Python n'est pas trouve dans le PATH ni dans les chemins Conda standards.
                echo.
                echo Solutions:
                echo   1. Lancez: call "C:\ProgramData\Miniconda3\Scripts\activate.bat" astroenv
                echo   2. Ou installez Python: https://www.python.org/downloads/
                echo ============================================================
                pause
                exit /b 1
            )
        )
    )
)

echo Utilisation de: %PYTHON_CMD%
echo.

REM Aller dans le repertoire build
cd build

REM Executer le script Python
echo Demarrage de la construction...
echo.
%PYTHON_CMD% build_distribution.py %PROFILE%

REM Stocker le code de retour
set BUILD_STATUS=%ERRORLEVEL%

REM Revenir a la racine
cd ..

echo.
echo ============================================================

if %BUILD_STATUS% EQU 0 (
    echo.
    echo [OK] DISTRIBUTION CONSTRUITE AVEC SUCCES!
    echo.
    echo Repertoire: build\distributions\%PROFILE%
    echo Archive: build\distributions\%PROFILE%.zip
    echo.
    echo La distribution est prete a etre utilisee.
    echo.
) else (
    echo.
    echo [ERREUR] Echec de la construction de la distribution.
    echo.
    echo Code d'erreur: %BUILD_STATUS%
    echo.
    echo Verifiez les messages d'erreur ci-dessus.
    echo.
    echo Aide: build\TROUBLESHOOTING.md et build\GUIDE_UTILISATION.md
    echo.
)

echo ============================================================
echo.

REM Demander si on veut ouvrir le repertoire de build
if %BUILD_STATUS% EQU 0 (
    set /p OPEN="Ouvrir le repertoire de build? (O/N): "
    if /i "!OPEN!"=="O" (
        if exist "build\distributions\%PROFILE%" (
            explorer "build\distributions\%PROFILE%"
        )
    )
)

pause
