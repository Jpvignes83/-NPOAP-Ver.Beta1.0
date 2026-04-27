@echo off
chcp 65001 >nul
REM ============================================================
REM Lancement NPOAP - environnement conda astroenv
REM Ne pas utiliser conda activate seul: sans conda init cmd
REM On appelle activate.bat directement (comme LANCEMENT.bat).
REM ============================================================

title NPOAP - Lancement (astroenv)

cd /d "%~dp0"

set "CONDA_ROOT="
where conda >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    for /f "delims=" %%i in ('conda info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%ProgramData%\Miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%USERPROFILE%\miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%LOCALAPPDATA%\miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%LOCALAPPDATA%\miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%ProgramData%\Anaconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%ProgramData%\Anaconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%USERPROFILE%\anaconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%USERPROFILE%\anaconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT (
    echo ERREUR: conda introuvable. Installez Miniconda puis creez astroenv ^(installation.bat ou INSTALLER_NPOAP_ASTROENV_WINDOWS.bat^).
    pause
    exit /b 1
)

if not exist "%CONDA_ROOT%\Scripts\activate.bat" (
    echo ERREUR: activate.bat introuvable dans "%CONDA_ROOT%\Scripts"
    pause
    exit /b 1
)

call "%CONDA_ROOT%\Scripts\activate.bat" astroenv
if errorlevel 1 (
    echo ERREUR: impossible d'activer astroenv.
    echo Exemple : conda create -n astroenv python=3.11 -y
    pause
    exit /b 1
)

REM Mode par defaut: compatibilite (autorise user-site) pour eviter les blocages
REM si astroenv en ProgramData est non-ecrivable sans droits admin.
REM Pour forcer l'isolation stricte: set NPOAP_STRICT_ENV=1 avant lancement.
if /i "%NPOAP_STRICT_ENV%"=="1" (
    set "PYTHONNOUSERSITE=1"
    set "PIP_USER=0"
)

REM Backend Prospector: WSL par defaut (FSPS Linux), surchargeable via variable d'environnement.
if not defined NPOAP_PROSPECTOR_BACKEND set "NPOAP_PROSPECTOR_BACKEND=wsl"

echo Lancement de NPOAP...
echo.

if not exist "main.py" (
    echo Erreur: main.py non trouve dans %CD%
    pause
    exit /b 1
)

if not exist "logs" mkdir logs

python main.py

if errorlevel 1 (
    echo.
    echo ERREUR lors du lancement de NPOAP ^(code %ERRORLEVEL%^).
    pause
)
