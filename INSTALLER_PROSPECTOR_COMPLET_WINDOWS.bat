@echo off
REM NPOAP - Lance l'installateur PowerShell Prospector (meme dossier que ce .bat).
setlocal EnableDelayedExpansion
cd /d "%~dp0"
if not exist "%~dp0INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1" (
    echo ERREUR : INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1 introuvable.
    pause
    exit /b 1
)
call :run_ps
set "EX=!ERRORLEVEL!"
if not "!EX!"=="0" echo Code de sortie PowerShell : !EX!
pause
exit /b !EX!

:run_ps
start "NPOAP Prospector FSPS" /wait powershell -NoLogo -NoProfile -ExecutionPolicy Bypass -NoExit -File "%~dp0INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1" -InstallFSPS %*
exit /b %ERRORLEVEL%
