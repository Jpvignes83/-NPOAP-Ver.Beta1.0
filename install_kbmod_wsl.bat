@echo off
REM NPOAP — Installation automatique KBMOD sous WSL (clone, build, pip -e).
setlocal
cd /d "%~dp0"

echo.
echo ============================================================
echo NPOAP — Installation KBMOD sous WSL
echo ============================================================
echo Documentation : docs\INSTALL_KBMOD_WSL.md
echo Script Linux   : scripts\install_kbmod_wsl.sh
echo.
echo Mot de passe sudo possible dans WSL pour apt install.
echo Option : install_kbmod_wsl.bat skip-apt   ^(ignore apt si deja fait^)
echo.

where wsl >nul 2>&1
if errorlevel 1 (
  echo ERREUR : WSL introuvable. Installez WSL puis relancez ^(install_wsl.bat^).
  goto :finish_error
)

set "SH=%~dp0scripts\install_kbmod_wsl.sh"
if not exist "%SH%" (
  echo ERREUR : Script introuvable : %SH%
  goto :finish_error
)

for /f "delims=" %%i in ('wsl wslpath -a "%SH%" 2^>nul') do set "SHWSL=%%i"
if not defined SHWSL (
  echo ERREUR : wslpath a echoue. Verifiez que la distribution WSL par defaut demarre ^(wsl echo ok^).
  goto :finish_error
)

if /i "%~1"=="skip-apt" (
  echo Mode : SKIP_APT=1
  wsl env SKIP_APT=1 bash "%SHWSL%"
) else (
  wsl bash "%SHWSL%"
)

if errorlevel 1 goto :finish_error
echo.
echo ============================================================
echo Installation KBMOD terminee avec succes.
echo ============================================================
goto :finish_ok

:finish_error
echo.
echo ============================================================
echo Installation KBMOD terminee avec erreur ^(code %ERRORLEVEL%^).
echo ============================================================
set ERR=1
goto :end

:finish_ok
set ERR=0

:end
if /i "%NPOAP_AUTO%"=="1" exit /b %ERR%
pause
exit /b %ERR%
