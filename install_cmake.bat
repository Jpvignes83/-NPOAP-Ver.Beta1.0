@echo off
REM NPOAP - Aide a installer CMake (compilation KBMOD et autres projets C++).
setlocal
cd /d "%~dp0"
echo.
echo Option 1 : page de telechargement officielle CMake.
start "" "https://cmake.org/download/"
echo.
echo Option 2 : depuis un terminal avec winget (sans ce script) :
echo   winget install Kitware.CMake
echo.
pause
