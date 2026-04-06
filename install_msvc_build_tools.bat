@echo off
REM NPOAP - Ouvre la page officielle Visual C++ Build Tools (MSVC).
REM Installez la charge utile "Desktop development with C++" ou equivalent x64 si besoin pour PHOEBE2 / compilation.
setlocal
cd /d "%~dp0"
echo.
echo Ouverture du site Visual C++ Build Tools...
echo Fermez cette fenetre quand vous avez termine.
echo.
start "" "https://visualstudio.microsoft.com/visual-cpp-build-tools/"
pause
