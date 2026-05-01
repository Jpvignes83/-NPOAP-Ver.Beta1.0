@echo off
chcp 65001 >nul
set "HERE=%~dp0"
set "PS1=%HERE%exominer_copy_pipeline_data.ps1"
if not exist "%PS1%" (
  echo Script introuvable : %PS1%
  exit /b 1
)
powershell -NoProfile -ExecutionPolicy Bypass -File "%PS1%"
exit /b %ERRORLEVEL%
