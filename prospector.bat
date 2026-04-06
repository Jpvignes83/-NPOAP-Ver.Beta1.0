@echo off
REM NPOAP — Lance l'installation Prospector + FSPS dans WSL2 (Linux).
REM Necessite install_prospector_wsl.sh dans le meme dossier.
setlocal
cd /d "%~dp0"

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

if not exist "%~dp0install_prospector_wsl.sh" (
  echo ERREUR : install_prospector_wsl.sh est introuvable dans :
  echo   %~dp0
  echo Les deux fichiers doivent rester dans le meme dossier.
  pause
  exit /b 1
)

echo Demarrage du script Linux dans WSL...
echo Un mot de passe peut etre demande pour "sudo apt-get".
echo.

REM WSL 0.64+ : --cd accepte un chemin Windows (repertoire de ce .bat)
wsl --cd "%CD%" -- bash ./install_prospector_wsl.sh
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
