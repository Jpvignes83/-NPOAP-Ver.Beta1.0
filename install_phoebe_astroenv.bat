@echo off
setlocal

REM Installation one-click de PHOEBE2 dans le contexte astroenv NPOAP.
REM A lancer depuis NPOAP (double-clic ou terminal).

cd /d "%~dp0"

if not exist "activate_astroenv.bat" (
    echo [ERREUR] activate_astroenv.bat introuvable dans %CD%
    pause
    exit /b 1
)

call "activate_astroenv.bat"
if errorlevel 1 (
    echo [ERREUR] Activation astroenv echouee.
    pause
    exit /b 1
)

echo.
echo [INFO] Installation de PHOEBE ^(seul^)...
python -m pip install --prefer-binary "phoebe>=2.4.22,<2.5"
if errorlevel 1 (
    echo [ERREUR] Echec de l'installation PHOEBE.
    pause
    exit /b 1
)

echo.
echo [INFO] Verification de l'import PHOEBE...
python -c "import sys; m=type('M',(),{'add_history':lambda *a,**k:None,'read_history_file':lambda *a,**k:None,'write_history_file':lambda *a,**k:None,'set_history_length':lambda *a,**k:None,'get_history_length':lambda *a,**k:0,'set_completer':lambda *a,**k:None,'set_completer_delims':lambda *a,**k:None,'set_completion_display_matches_hook':lambda *a,**k:None,'parse_and_bind':lambda *a,**k:None}); sys.modules.setdefault('readline', m()); import phoebe; print('[OK] PHOEBE', phoebe.__version__)"
if errorlevel 1 (
    echo [ATTENTION] Installation terminee, mais verification import en echec.
    echo Essayez directement dans NPOAP ^(l'onglet PHOEBE applique deja le fallback readline sous Windows^).
    pause
    exit /b 1
)

echo.
echo [OK] PHOEBE est installe et verifie.
pause
exit /b 0
