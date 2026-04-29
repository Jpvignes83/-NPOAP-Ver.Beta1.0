@echo off
setlocal EnableDelayedExpansion
REM Script batch pour lancer les tests de l'application NPOAP
REM Utilise uniquement l'interprete conda astroenv ^(pas AppData\Roaming\Python311^).

echo Activation de l'environnement astroenv...

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

if not defined CONDA_ROOT (
    echo Erreur: conda introuvable.
    pause
    exit /b 1
)

call "%CONDA_ROOT%\Scripts\activate.bat" astroenv

if errorlevel 1 (
    echo Erreur: Impossible d'activer l'environnement astroenv
    echo Essayez : call "%CONDA_ROOT%\Scripts\activate.bat" astroenv
    pause
    exit /b 1
)

set "PYTHONNOUSERSITE=1"
set "PIP_USER="
set "PY_ASTRO=%CONDA_ROOT%\envs\astroenv\python.exe"
if not exist "!PY_ASTRO!" (
    echo Erreur: !PY_ASTRO! introuvable
    pause
    exit /b 1
)

echo.
echo Lancement des tests...
echo.

"!PY_ASTRO!" test_application_launch.py

if errorlevel 1 (
    echo.
    echo Les tests ont detecte des erreurs. Verifiez les details ci-dessus.
    pause
    exit /b 1
) else (
    echo.
    echo Tous les tests sont passes avec succes !
    pause
)
