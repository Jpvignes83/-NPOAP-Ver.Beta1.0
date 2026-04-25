@echo off
REM NPOAP — Lance l'installation Prospector + FSPS dans WSL2 (Linux).
REM Necessite install_prospector_wsl.sh dans le meme dossier.
setlocal
cd /d "%~dp0"
set "SCRIPT_SH="

echo.
echo ========================================
echo  NPOAP - Prospector + FSPS via WSL2
echo ========================================
echo.

where wsl >nul 2>&1
if errorlevel 1 (
  echo ERREUR : la commande "wsl" est introuvable.
  echo Installez WSL2 depuis une invite PowerShell administrateur :
  echo   wsl --install
  echo Puis redemarrez si demande, installez Ubuntu depuis le Store si besoin.
  pause
  exit /b 1
)

if exist "%~dp0install_prospector_wsl.sh" (
  set "SCRIPT_SH=%~dp0install_prospector_wsl.sh"
)
if not defined SCRIPT_SH if exist "%~dp0Installation_fsps\install_prospector_wsl.sh" (
  set "SCRIPT_SH=%~dp0Installation_fsps\install_prospector_wsl.sh"
)

if not defined SCRIPT_SH (
  echo ERREUR : install_prospector_wsl.sh est introuvable dans :
  echo   %~dp0
  echo   %~dp0Installation_fsps\
  echo Placez prospector.bat et install_prospector_wsl.sh dans le meme dossier,
  echo ou gardez la structure standard du package NPOAP.
  pause
  exit /b 1
)

echo Demarrage du script Linux dans WSL...
echo Un mot de passe peut etre demande pour "sudo apt-get".
echo.

REM Convertir le chemin Windows du script en chemin Linux WSL (/mnt/c/...).
for /f "delims=" %%P in ('wsl wslpath -a "%SCRIPT_SH%"') do set "SCRIPT_WSL=%%P"
if not defined SCRIPT_WSL (
  echo ERREUR : impossible de convertir le chemin du script pour WSL.
  echo Script Windows: %SCRIPT_SH%
  pause
  exit /b 1
)

REM WSL 0.64+ : --cd accepte un chemin Windows (repertoire de ce .bat).
wsl --cd "%CD%" -- bash "%SCRIPT_WSL%"
set "ERR=%ERRORLEVEL%"

if not "%ERR%"=="0" (
  echo.
  echo Echec ^(code %ERR%^). Si "wsl --cd" n'est pas reconnu, mettez a jour WSL :
  echo   wsl --update
  echo Ou lancez depuis WSL :  bash "%CD%\install_prospector_wsl.sh"
  echo ^(apres conversion du chemin vers /mnt/c/...^)
)

echo.
pause
exit /b %ERR%
