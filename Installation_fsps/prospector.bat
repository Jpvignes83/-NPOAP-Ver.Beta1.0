@echo off
REM NPOAP — Lance l'installation Prospector + FSPS dans WSL2 (Linux).
REM Necessite install_prospector_wsl.sh dans le meme dossier.
setlocal EnableDelayedExpansion
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

REM Recherche du .sh : %%~fI normalise les ".." ^(dossier Installation_fsps duplique, etc.^).
set "SCRIPT_SH="
for %%I in (
  "%~dp0install_prospector_wsl.sh"
  "%~dp0..\install_prospector_wsl.sh"
  "%~dp0..\..\install_prospector_wsl.sh"
  "%~dp0..\..\..\install_prospector_wsl.sh"
  "%~dp0..\..\..\..\install_prospector_wsl.sh"
  "%~dp0..\..\..\..\..\install_prospector_wsl.sh"
  "%~dp0..\..\..\..\..\..\install_prospector_wsl.sh"
  "%~dp0Installation_fsps\install_prospector_wsl.sh"
  "%~dp0..\Installation_fsps\install_prospector_wsl.sh"
  "%~dp0..\..\Installation_fsps\install_prospector_wsl.sh"
  "%~dp0..\..\..\Installation_fsps\install_prospector_wsl.sh"
) do (
  if exist "%%~fI" (
    set "SCRIPT_SH=%%~fI"
    goto :script_sh_found
  )
)

echo ERREUR : install_prospector_wsl.sh introuvable.
echo Dossier du prospector.bat : %~dp0
echo Verifiez que install_prospector_wsl.sh est dans ce dossier ou dans un parent
echo ^(jusqu'a 6 niveaux au-dessus^) ou sous ...\Installation_fsps\
echo.
echo Astuce : copiez les deux fichiers depuis le depot vers
echo   C:\NPOAP\Installation_fsps\
echo puis lancez C:\NPOAP\Installation_fsps\prospector.bat
pause
exit /b 1

:script_sh_found

REM Executer WSL depuis le repertoire du script ^(chemins relatifs eventuels^).
for %%J in ("%SCRIPT_SH%") do cd /d "%%~dpJ"

echo Demarrage du script Linux dans WSL...
echo Un mot de passe peut etre demande pour "sudo apt-get".
echo.

REM Convertir le chemin Windows du script en chemin Linux WSL (/mnt/c/...).
for /f "delims=" %%P in ('wsl wslpath -a "!SCRIPT_SH!"') do set "SCRIPT_WSL=%%P"
if not defined SCRIPT_WSL (
  echo ERREUR : impossible de convertir le chemin du script pour WSL.
  echo Script Windows: !SCRIPT_SH!
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
