@echo off
REM NPOAP - Lance l'installateur PowerShell Prospector (meme dossier que ce .bat).
setlocal
cd /d "%~dp0"
if not exist "%~dp0INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1" (
    echo ERREUR : INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 introuvable.
    pause
    exit /b 1
)
powershell -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1" %*
set "EX=%errorLevel%"
if not "%EX%"=="0" echo Code de sortie PowerShell : %EX%
pause
exit /b %EX%
