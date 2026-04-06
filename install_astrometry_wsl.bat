@echo off
REM NPOAP - Installe astrometry.net dans la distro WSL par defaut (souvent Ubuntu).
REM Une fenetre WSL peut demander le mot de passe sudo.
setlocal
cd /d "%~dp0"
echo.
echo Mise a jour des paquets et installation astrometry.net via WSL...
wsl bash -c "sudo apt-get update && sudo apt-get install -y astrometry.net"
echo.
if %errorLevel% neq 0 (
    echo Echec : ouvrez une session WSL ^(wsl^), verifiez sudo, puis :
    echo   sudo apt-get update ^&^& sudo apt-get install -y astrometry.net
) else (
    echo OK : astrometry.net installe dans WSL.
)
pause
