@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul
cd /d "%~dp0"

set "LOG=%~dp0install_sora_analysis.log"
set "ENV_NAME=astroenv"

echo ============================================================ > "%LOG%"
echo NPOAP - Paquet sora-astro ^(analyse SORA^) - %date% %time% >> "%LOG%"
echo Repertoire: %CD% >> "%LOG%"
echo ============================================================ >> "%LOG%"

echo.
echo ============================================================
echo   NPOAP - Installation / mise a jour : sora-astro ^(analyse^)
echo ============================================================
echo   Environnement Conda : %ENV_NAME% ^(meme env que NPOAP^)
echo   Journal : %LOG%
echo ============================================================
echo.

set "CONDA_EXE="
if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "CONDA_EXE=%ProgramData%\Miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%ProgramData%\Anaconda3\Scripts\conda.exe" set "CONDA_EXE=%ProgramData%\Anaconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%LocalAppData%\miniconda3\Scripts\conda.exe" set "CONDA_EXE=%LocalAppData%\miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%UserProfile%\miniconda3\Scripts\conda.exe" set "CONDA_EXE=%UserProfile%\miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%UserProfile%\Miniconda3\Scripts\conda.exe" set "CONDA_EXE=%UserProfile%\Miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%UserProfile%\Anaconda3\Scripts\conda.exe" set "CONDA_EXE=%UserProfile%\Anaconda3\Scripts\conda.exe"

if not defined CONDA_EXE (
    echo ERREUR : conda.exe introuvable.
    echo ERREUR : conda.exe introuvable. >> "%LOG%"
    goto :end_fail
)

echo Conda : !CONDA_EXE!
echo Conda : !CONDA_EXE! >> "%LOG%"
echo.

set "ENV_LIST=%TEMP%\npoap_conda_envs_sora_%RANDOM%.txt"
"!CONDA_EXE!" info --envs >"%ENV_LIST%" 2>>"%LOG%"
findstr /i "%ENV_NAME%" "%ENV_LIST%" >nul 2>&1
set "ENV_FOUND=!errorlevel!"
del "%ENV_LIST%" >nul 2>&1
if not "!ENV_FOUND!"=="0" (
    echo ERREUR : l'environnement "%ENV_NAME%" est introuvable.
    echo Creez-le avec installation.bat ^(NPOAP^) ou : conda create -n %ENV_NAME% python=3.11 -y
    echo ERREUR : env %ENV_NAME% absent >> "%LOG%"
    goto :end_fail
)

echo Installation / mise a jour : pip, sora-astro ^>=0.3.2...
echo pip install sora-astro >> "%LOG%"
"!CONDA_EXE!" run -n %ENV_NAME% python -m pip install --upgrade pip 1>>"%LOG%" 2>>&1
if errorlevel 1 (
    echo AVERTISSEMENT : mise a jour pip a echoue ^(voir journal^).
)
"!CONDA_EXE!" run -n %ENV_NAME% python -m pip install "sora-astro>=0.3.2" 1>>"%LOG%" 2>>&1
if errorlevel 1 (
    echo ERREUR : pip install sora-astro a echoue. Voir %LOG%
    goto :end_fail
)

echo.
echo Test d'import...
"!CONDA_EXE!" run -n %ENV_NAME% python -c "import sora; from sora.prediction import prediction; print('SORA OK:', getattr(sora, '__version__', '?'))" 1>>"%LOG%" 2>>&1
if errorlevel 1 (
    echo AVERTISSEMENT : le test d'import a echoue. Voir %LOG% ^(Cartopy / dependances ?^)
) else (
    echo Import SORA reussi.
)

echo.
echo ============================================================
echo Termine. Redemarrez NPOAP si l'application etait ouverte.
echo ============================================================
echo.
goto :end_ok

:end_fail
echo.
pause
exit /b 1

:end_ok
pause
exit /b 0
