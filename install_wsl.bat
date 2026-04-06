@echo off
REM NPOAP - Installe le sous-systeme Windows pour Linux (WSL).
REM Executez en administrateur si la commande echoue.
setlocal
cd /d "%~dp0"
echo.
echo Installation / mise a jour de WSL...
echo Un redemarrage peut etre demande par Windows.
echo.
wsl --install
echo.
echo Code de sortie : %errorLevel%
echo Ensuite : install_ubuntu_wsl.bat pour ajouter la distro Ubuntu si besoin.
pause
