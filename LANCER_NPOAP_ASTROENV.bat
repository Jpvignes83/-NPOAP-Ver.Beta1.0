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

REM Par defaut : isolation conda — pas de site-packages utilisateur Windows
REM ^(%AppData%\Roaming\Python\Python311\site-packages^).
REM Pour autoriser ce site ^(cas rare^) : definir NPOAP_ALLOW_USERSITE=1 avant ce script.
set "PYTHONNOUSERSITE=1"
set "PIP_USER="
if /i "%NPOAP_ALLOW_USERSITE%"=="1" (
    set "PYTHONNOUSERSITE="
)

REM Backend Prospector: WSL par defaut (FSPS Linux), surchargeable via variable d'environnement.
if not defined NPOAP_PROSPECTOR_BACKEND set "NPOAP_PROSPECTOR_BACKEND=wsl"

set "PY_ASTRO=%CONDA_ROOT%\envs\astroenv\python.exe"

echo Lancement de NPOAP...
echo.

if not exist "main.py" (
    echo Erreur: main.py non trouve dans %CD%
    pause
    exit /b 1
)

if not exist "%PY_ASTRO%" (
    echo ERREUR: interprete conda introuvable : %PY_ASTRO%
    pause
    exit /b 1
)

if not exist "logs" mkdir logs

"%PY_ASTRO%" main.py

if errorlevel 1 (
    echo.
    echo ERREUR lors du lancement de NPOAP ^(code %ERRORLEVEL%^).
    pause
)
