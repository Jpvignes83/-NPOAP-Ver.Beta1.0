@echo off
REM NPOAP - KBMOD se compile sous Linux / WSL (CUDA, CMake, toolchain).
setlocal
cd /d "%~dp0"
echo.
echo Documentation : docs\INSTALL_KBMOD_WSL.md
echo Dependances pip dediees : requirements-kbmod.txt ^(a utiliser dans l'environnement Linux/WSL^).
if exist "%~dp0docs\INSTALL_KBMOD_WSL.md" (
    start "" "%~dp0docs\INSTALL_KBMOD_WSL.md"
) else (
    echo Fichier docs\INSTALL_KBMOD_WSL.md introuvable dans ce dossier.
)
echo.
pause
