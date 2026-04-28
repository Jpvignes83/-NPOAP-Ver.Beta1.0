@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul
cd /d "%~dp0"

set "LOG=%~dp0install_pymovie_env.log"
set "ENV_NAME=pymovie"
set "PY_VER=3.10"

echo ============================================================ > "%LOG%"
echo NPOAP - Environnement Conda PyMovie / PyOTE - %date% %time% >> "%LOG%"
echo Repertoire: %CD% >> "%LOG%"
echo ============================================================ >> "%LOG%"

echo.
echo ============================================================
echo   NPOAP - Installation environnement Conda : PyMovie + PyOTE
echo ============================================================
echo   Environnement : %ENV_NAME%  ^|  Python : %PY_VER%
echo   Journal : %LOG%
echo ============================================================
echo.

REM --- Localiser conda (ordre courant Windows / NPOAP) ---
set "CONDA_EXE="
if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "CONDA_EXE=%ProgramData%\Miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%ProgramData%\Anaconda3\Scripts\conda.exe" set "CONDA_EXE=%ProgramData%\Anaconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%LocalAppData%\miniconda3\Scripts\conda.exe" set "CONDA_EXE=%LocalAppData%\miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%UserProfile%\miniconda3\Scripts\conda.exe" set "CONDA_EXE=%UserProfile%\miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%UserProfile%\Miniconda3\Scripts\conda.exe" set "CONDA_EXE=%UserProfile%\Miniconda3\Scripts\conda.exe"
if not defined CONDA_EXE if exist "%UserProfile%\Anaconda3\Scripts\conda.exe" set "CONDA_EXE=%UserProfile%\Anaconda3\Scripts\conda.exe"

if not defined CONDA_EXE (
    echo ERREUR : conda.exe introuvable.
    echo Installez Miniconda/Anaconda ou indiquez le chemin dans ce script.
    echo ERREUR : conda.exe introuvable. >> "%LOG%"
    goto :end_fail
)

echo Conda : !CONDA_EXE!
echo Conda : !CONDA_EXE! >> "%LOG%"
echo.

set "ENV_LIST=%TEMP%\npoap_conda_envs_%RANDOM%.txt"
"!CONDA_EXE!" info --envs >"%ENV_LIST%" 2>>"%LOG%"
findstr /i "%ENV_NAME%" "%ENV_LIST%" >nul 2>&1
del "%ENV_LIST%" >nul 2>&1
if errorlevel 1 (
    echo Creation de l'environnement %ENV_NAME% ^(Python %PY_VER%^)...
    echo Creation environnement %ENV_NAME% >> "%LOG%"
    "!CONDA_EXE!" create -n %ENV_NAME% python=%PY_VER% -y 1>>"%LOG%" 2>>&1
    if errorlevel 1 (
        echo ERREUR : conda create a echoue. Voir %LOG%
        goto :end_fail
    )
) else (
    echo Environnement %ENV_NAME% deja present : mise a jour des paquets pip.
    echo Env %ENV_NAME% existe deja >> "%LOG%"
)

echo.
echo Installation / mise a jour : pip, pymovie, pyote...
echo pip install >> "%LOG%"
"!CONDA_EXE!" run -n %ENV_NAME% python -m pip install --upgrade pip 1>>"%LOG%" 2>>&1
if errorlevel 1 (
    echo AVERTISSEMENT : mise a jour pip a echoue ^(voir journal^).
)
"!CONDA_EXE!" run -n %ENV_NAME% python -m pip install pymovie pyote 1>>"%LOG%" 2>>&1
if errorlevel 1 (
    echo ERREUR : pip install pymovie pyote a echoue. Voir %LOG%
    goto :end_fail
)

echo.
echo Test d'import rapide...
"!CONDA_EXE!" run -n %ENV_NAME% python -c "from pymovie import main; from pyoteapp import pyote; print('Imports OK')" 1>>"%LOG%" 2>>&1
if errorlevel 1 (
    echo AVERTISSEMENT : le test d'import a echoue ^(PyMovie peut quand meme fonctionner selon la version^). Voir %LOG%
) else (
    echo Imports OK.
)

echo.
echo ============================================================
echo Termine. Environnement : conda activate %ENV_NAME%
echo Racine type : %USERPROFILE%\.conda\envs\%ENV_NAME%
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
