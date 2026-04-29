@echo off
REM NPOAP - Installation automatique KBMOD sous WSL (clone, build, pip -e).
REM Sans argument : installe uniquement dans conda env astroenv (defaut NPOAP).
REM Arguments facultatifs :
REM   skip-apt           - ne pas lancer apt-get (build tools deja installes)
REM   astroenv           - idem defaut (conda .../envs/astroenv/bin/python3)
REM   system             - Python du PATH systeme (/usr/bin/python3), non recommande
REM   chemin/python3     - interpreter WSL exact
REM Exemples :
REM   install_kbmod_wsl.bat
REM   install_kbmod_wsl.bat skip-apt
setlocal EnableDelayedExpansion
chcp 65001 >nul
cd /d "%~dp0"

echo.
echo ============================================================
echo NPOAP - Installation KBMOD sous WSL
echo ============================================================
echo Documentation : docs\INSTALL_KBMOD_WSL.md
echo Script Linux   : scripts\install_kbmod_wsl.sh
echo.
echo Defaut : conda astroenv uniquement ^(pas le Python systeme WSL^).
echo Option   : skip-apt ^| system ^| ou chemin complet vers python3 sous WSL
echo.
echo Mot de passe sudo possible dans WSL pour apt install ^(sauf skip-apt^).
echo Purge ~/.local ^(anciennes installs systeme^) : avant ce .bat,  CMD : set KBMOD_PURGE_USER_KBMOD=1
echo.

where wsl >nul 2>&1
if errorlevel 1 (
  echo ERREUR : WSL introuvable. Installez WSL puis relancez ^(install_wsl.bat^).
  goto :finish_error
)

set "SH=%~dp0scripts\install_kbmod_wsl.sh"
if not exist "%SH%" (
  echo ERREUR : Script introuvable : %SH%
  goto :finish_error
)

for /f "delims=" %%i in ('wsl wslpath -a "%SH%" 2^>nul') do set "SHWSL=%%i"
if not defined SHWSL (
  echo ERREUR : wslpath a echoue. Verifiez que la distribution WSL par defaut demarre ^(wsl echo ok^).
  goto :finish_error
)

set "SKIP_APT=0"
set "KBMOD_TP="
:parse_args
if "%~1"=="" goto :parse_done
if /i "%~1"=="skip-apt" (
  set "SKIP_APT=1"
  shift
  goto :parse_args
)
set "KBMOD_TP=%~1"
shift
goto :parse_args
:parse_done

if not defined KBMOD_TP set "KBMOD_TP=astroenv"

set "KU=0"
if /i "!KBMOD_PURGE_USER_KBMOD!"=="1" set "KU=1"

if "!SKIP_APT!"=="1" echo Mode : SKIP_APT=1
echo KBMOD_TARGET_PYTHON=!KBMOD_TP!
if "!KU!"=="1" echo KBMOD_PURGE_USER_KBMOD=1 ^(desinstalle kbmod du Python systeme ~/.local si present^)
echo.

wsl.exe env SKIP_APT=!SKIP_APT! KBMOD_TARGET_PYTHON=!KBMOD_TP! KBMOD_PURGE_USER_KBMOD=!KU! bash "!SHWSL!"

if errorlevel 1 goto :finish_error
echo.
echo ============================================================
echo Installation KBMOD terminee avec succes.
echo ============================================================
goto :finish_ok

:finish_error
echo.
echo ============================================================
echo Installation KBMOD terminee avec erreur ^(code %ERRORLEVEL%^).
echo ============================================================
set ERR=1
goto :end

:finish_ok
set ERR=0

:end
if /i "%NPOAP_AUTO%"=="1" exit /b %ERR%
pause
exit /b %ERR%
