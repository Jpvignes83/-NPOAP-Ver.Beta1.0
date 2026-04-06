@echo off
chcp 65001 >nul
REM ============================================================
REM Script de lancement NPOAP - full
REM Obligatoire : environnement conda « astroenv » (Python 3.11+).
REM Installation : INSTALLER_NPOAP_ASTROENV_WINDOWS.bat
REM ============================================================

title NPOAP - Lancement

cd /d "%~dp0"

set "CONDA_ROOT="
where conda >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    for /f "delims=" %%i in ('conda info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%USERPROFILE%\miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%USERPROFILE%\anaconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%USERPROFILE%\anaconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT (
    echo ERREUR: conda introuvable. Installez Miniconda puis creez astroenv :
    echo   INSTALLER_NPOAP_ASTROENV_WINDOWS.bat
    pause
    exit /b 1
)

if not exist "%CONDA_ROOT%\Scripts\activate.bat" (
    echo ERREUR: activate.bat introuvable.
    pause
    exit /b 1
)

call "%CONDA_ROOT%\Scripts\activate.bat" astroenv
if errorlevel 1 (
    echo ERREUR: impossible d'activer astroenv.
    echo Executez INSTALLER_NPOAP_ASTROENV_WINDOWS.bat dans ce dossier.
    pause
    exit /b 1
)

echo Lancement de NPOAP...
echo.

if not exist "main.py" (
    echo Erreur: main.py non trouve.
    pause
    exit /b 1
)

if not exist "logs" mkdir logs

python main.py

if errorlevel 1 (
    echo.
    echo Erreur lors du lancement ^(code %ERRORLEVEL%^).
    echo Dependances : pip install -r requirements.txt dans astroenv
    pause
)
