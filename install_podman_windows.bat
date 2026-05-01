@echo off
chcp 65001 >nul
title NPOAP - Installation Podman (winget)
echo.
echo Installe Podman via winget (paquet RedHat.Podman).
echo Executez ce fichier en invite de commandes ouvre en Administrateur si winget le demande (dependances WSL, etc.).
echo.

set "WINGET=%LOCALAPPDATA%\Microsoft\WindowsApps\winget.exe"
if not exist "%WINGET%" (
  where winget >nul 2>&1
  if errorlevel 1 (
    echo ERREUR : winget introuvable.
    echo Installez "App Installer" depuis le Microsoft Store, ou ouvrez :
    echo https://podman.io/getting-started/installation
    pause
    exit /b 1
  )
  set "WINGET=winget"
)

"%WINGET%" install -e --id RedHat.Podman --accept-package-agreements --accept-source-agreements
set EC=%ERRORLEVEL%
echo.
if %EC%==0 (
  echo Podman installe. Si une NOUVELLE installation :
  echo   1^) Ouvrez un NOUVEAU terminal ^(PATH mis a jour^).
  echo   2^) podman machine init     ^(une seule fois ; peut prendre plusieurs minutes^)
  echo   3^) podman machine start    ^(a chaque session ou au demarrage si vous preferez^)
  echo Ensuite relancez scripts\exominer_copy_pipeline_data.bat pour les donnees ExoMiner.
) else (
  echo Code retour : %EC%
)
pause
exit /b %EC%
