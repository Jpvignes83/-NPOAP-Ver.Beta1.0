@echo off
REM NPOAP - Enregistre la distribution Ubuntu dans WSL (apres install_wsl.bat).
setlocal
cd /d "%~dp0"
echo.
echo Installation de la distribution Ubuntu dans WSL...
wsl --install -d Ubuntu
echo.
if %errorLevel% neq 0 (
    echo Si le message indique que la distribution existe deja, c'est normal.
    echo Lancez Ubuntu depuis le menu demarrer ou : wsl -d Ubuntu
)
pause
