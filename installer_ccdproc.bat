@echo off
setlocal EnableDelayedExpansion
REM Script pour installer ccdproc dans l'environnement conda astroenv uniquement.

echo ========================================
echo Installation de ccdproc
echo ========================================
echo.

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

if not defined CONDA_ROOT (
    echo ERREUR: conda introuvable.
    pause
    exit /b 1
)

echo Activation de l'environnement astroenv...
call "%CONDA_ROOT%\Scripts\activate.bat" astroenv

if errorlevel 1 (
    echo ERREUR: Impossible d'activer l'environnement astroenv
    echo Assurez-vous que conda est installe et que l'environnement astroenv existe.
    pause
    exit /b 1
)

set "PYTHONNOUSERSITE=1"
set "PIP_USER="
set "PY_ASTRO=%CONDA_ROOT%\envs\astroenv\python.exe"
if not exist "!PY_ASTRO!" (
    echo ERREUR: !PY_ASTRO! introuvable
    pause
    exit /b 1
)

echo.
echo Installation de ccdproc^>=2.4.0...
"!PY_ASTRO!" -m pip install "ccdproc>=2.4.0"

if errorlevel 1 (
    echo.
    echo ERREUR lors de l'installation de ccdproc
    pause
    exit /b 1
)

echo.
echo ========================================
echo Installation terminee avec succes!
echo ========================================
echo.
echo ccdproc est maintenant installe dans l'environnement astroenv.
echo Vous pouvez relancer NPOAP pour utiliser le pre-processing complet.
echo.
pause
