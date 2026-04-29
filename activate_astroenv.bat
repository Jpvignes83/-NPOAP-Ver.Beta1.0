@echo off
REM Activation astroenv pour scripts NPOAP (sans conda init dans cmd).
REM Definit PY_ASTRO = interprete Python de l'environnement (chemins pip explicites).
REM Pip et imports ignorent le site utilisateur Windows %AppData%\Roaming\Python\Python311 (PYTHONNOUSERSITE=1).

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
    echo [ERREUR] conda introuvable.
    exit /b 1
)

if not exist "%CONDA_ROOT%\Scripts\activate.bat" (
    echo [ERREUR] activate.bat introuvable dans "%CONDA_ROOT%\Scripts"
    exit /b 1
)

call "%CONDA_ROOT%\Scripts\activate.bat" astroenv
if errorlevel 1 exit /b 1

set "PYTHONNOUSERSITE=1"
set "PIP_USER="
set "PY_ASTRO=%CONDA_ROOT%\envs\astroenv\python.exe"

if not exist "%PY_ASTRO%" (
    echo [ERREUR] Environnement astroenv incomplet : "%PY_ASTRO%" introuvable.
    exit /b 1
)

exit /b 0
