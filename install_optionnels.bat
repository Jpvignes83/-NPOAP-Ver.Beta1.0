@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul

REM NPOAP - Installation automatique des composants optionnels
REM A lancer depuis le dossier NPOAP apres installation.bat

cd /d "%~dp0"

set "FAIL_COUNT=0"
set "MISS_COUNT=0"
set "DONE_WSL_PROSPECTOR=0"
set "AUTO_MODE=0"
if /i "%NPOAP_AUTO%"=="1" set "AUTO_MODE=1"
set "LOG_FILE=%~dp0install_optionnels.log"

echo ============================================================ > "%LOG_FILE%"
echo NPOAP - Installation optionnelle (demarrage: %date% %time%) >> "%LOG_FILE%"
echo Dossier: %CD% >> "%LOG_FILE%"
echo ============================================================ >> "%LOG_FILE%"

echo.
echo ============================================================
echo NPOAP - Installation des composants optionnels
echo ============================================================
echo Dossier: %CD%
echo Journal: %LOG_FILE%
if "!AUTO_MODE!"=="1" echo Mode auto: OUI ^(NPOAP_AUTO=1^)
echo.
echo Ce script lance automatiquement les installateurs optionnels.
echo Certains composants peuvent demander des droits administrateur ou un redemarrage.
echo.

call :run_bat_step "Visual C++ Build Tools / MSVC" "install_msvc_build_tools.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence
call :run_bat_step "CMake (build C++/KBMOD)" "install_cmake.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence
call :run_bat_step "Sous-systeme Windows pour Linux (WSL)" "install_wsl.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence
call :run_bat_step "Distribution Ubuntu dans WSL" "install_ubuntu_wsl.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence
call :run_bat_step "Astrometry.net dans WSL" "install_astrometry_wsl.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence
call :run_bat_step "KBMOD sous WSL/Linux (aide/doc)" "install_kbmod_wsl.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence

call :install_gfortran_windows
call :wait_next_step
if errorlevel 1 goto :end_sequence

REM Prospector Windows: on privilegie le .ps1 (plus complet), sinon fallback .bat.
call :run_prospector_windows
call :wait_next_step
if errorlevel 1 goto :end_sequence

call :run_bat_step "Prospector + FSPS via WSL" "Installation_fsps\prospector.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence
call :run_bat_step "SORA seul (reinstall astroenv)" "INSTALLER_SORA_ASTROENV.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence

call :install_sne_requirements
call :wait_next_step
if errorlevel 1 goto :end_sequence

:end_sequence

echo.>> "%LOG_FILE%"
echo ============================================================>> "%LOG_FILE%"
echo Fin: %date% %time% >> "%LOG_FILE%"
echo Echecs: !FAIL_COUNT! ; Fichiers manquants: !MISS_COUNT! >> "%LOG_FILE%"
echo ============================================================>> "%LOG_FILE%"

echo.
echo ============================================================
echo Installation optionnelle terminee
echo ============================================================
echo Echecs: !FAIL_COUNT!
echo Fichiers manquants: !MISS_COUNT!
echo Voir le journal: %LOG_FILE%
if "!AUTO_MODE!"=="1" echo Mode auto utilise: OUI
echo.
if not "!FAIL_COUNT!"=="0" (
    echo ATTENTION: certaines etapes ont echoue. Verifiez le journal.
)
if not "!MISS_COUNT!"=="0" (
    echo ATTENTION: certains scripts etaient absents. Verifiez le package.
)
echo.
pause
exit /b 0

:run_bat_step
set "STEP_LABEL=%~1"
set "STEP_FILE=%~2"
echo.
echo ------------------------------------------------------------
echo [OPTIONNEL] !STEP_LABEL!
echo Script: !STEP_FILE!
echo ------------------------------------------------------------

call :ask_install_or_skip "!STEP_LABEL!"
if errorlevel 1 (
    echo [SKIP] Ignore par l'utilisateur: !STEP_LABEL!
    echo [SKIP] Ignore par l'utilisateur: !STEP_LABEL!>> "%LOG_FILE%"
    exit /b 0
)

echo [START] !STEP_LABEL! ^(!STEP_FILE!^)>> "%LOG_FILE%"

if /i "!STEP_FILE!"=="install_wsl.bat" (
    call :has_ubuntu_wsl
    if "!errorlevel!"=="0" (
        echo [SKIP] Ubuntu WSL deja detectee - etape WSL ignoree
        echo [SKIP] Ubuntu WSL deja detectee - etape WSL ignoree>> "%LOG_FILE%"
        exit /b 0
    )
)

if /i "!STEP_FILE!"=="install_ubuntu_wsl.bat" (
    call :has_ubuntu_wsl
    if "!errorlevel!"=="0" (
        echo [SKIP] Ubuntu WSL deja installee - etape ignoree
        echo [SKIP] Ubuntu WSL deja installee - etape ignoree>> "%LOG_FILE%"
        exit /b 0
    )
)

if /i "!STEP_FILE!"=="Installation_fsps\prospector.bat" (
    if "!DONE_WSL_PROSPECTOR!"=="1" (
        echo [SKIP] Prospector WSL deja traite plus tot dans la sequence
        echo [SKIP] Prospector WSL deja traite plus tot dans la sequence>> "%LOG_FILE%"
        exit /b 0
    )
)

if not exist "!STEP_FILE!" (
    echo [SKIP] Fichier introuvable: !STEP_FILE!
    echo [SKIP] Fichier introuvable: !STEP_FILE!>> "%LOG_FILE%"
    set /a MISS_COUNT+=1
    exit /b 0
)

call "!STEP_FILE!"
set "RC=!errorlevel!"
if "!RC!"=="0" (
    echo [OK] !STEP_LABEL!
    echo [OK] !STEP_LABEL!>> "%LOG_FILE%"
) else (
    echo [ECHEC] !STEP_LABEL! ^(code !RC!^)
    echo [ECHEC] !STEP_LABEL! ^(code !RC!^)>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
)
exit /b 0

:has_ubuntu_wsl
REM Retourne 0 si une distribution Ubuntu est deja presente dans WSL, sinon 1.
where wsl >nul 2>&1
if errorlevel 1 exit /b 1

for /f "delims=" %%D in ('wsl -l -q 2^>nul') do (
    echo %%D | findstr /i "Ubuntu" >nul
    if not errorlevel 1 exit /b 0
)
exit /b 1

:install_gfortran_windows
echo.
echo ------------------------------------------------------------
echo [OPTIONNEL] gfortran Windows (pour FSPS natif)
echo ------------------------------------------------------------

call :ask_install_or_skip "gfortran Windows (pour FSPS natif)"
if errorlevel 1 (
    echo [SKIP] Ignore par l'utilisateur: gfortran Windows
    echo [SKIP] Ignore par l'utilisateur: gfortran Windows>> "%LOG_FILE%"
    exit /b 0
)

echo [START] gfortran Windows>> "%LOG_FILE%"

where gfortran >nul 2>&1
if not errorlevel 1 (
    for /f "delims=" %%G in ('where gfortran 2^>nul') do (
        echo [OK] gfortran detecte: %%G
        echo [OK] gfortran detecte: %%G>> "%LOG_FILE%"
        goto :gfortran_done
    )
)

echo [SKIP] gfortran non detecte - installation manuelle recommandee.
echo [SKIP] gfortran non detecte - installation manuelle recommandee.>> "%LOG_FILE%"
echo   1^)^ Ouvrez docs\MANUEL_INSTALLATION.md
echo   2^)^ Installez MinGW-w64/WinLibs puis rouvrez la console
echo   3^)^ Verifiez avec: where gfortran
start "" "https://winlibs.com/"
set /a MISS_COUNT+=1

:gfortran_done
exit /b 0

:run_prospector_windows
echo.
echo ------------------------------------------------------------
echo [OPTIONNEL] Prospector Windows (astroenv)
echo ------------------------------------------------------------

call :ask_install_or_skip "Prospector Windows (astroenv)"
if errorlevel 1 (
    echo [SKIP] Ignore par l'utilisateur: Prospector Windows
    echo [SKIP] Ignore par l'utilisateur: Prospector Windows>> "%LOG_FILE%"
    exit /b 0
)

echo [START] Prospector Windows>> "%LOG_FILE%"

if exist "INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1" (
    powershell -NoLogo -NoProfile -ExecutionPolicy Bypass -File "INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1"
    set "RC=!errorlevel!"
    if not defined RC set "RC=1"
    2>nul set /a RC_NUM=!RC!
    if not defined RC_NUM set "RC_NUM=1"
    if !RC_NUM! EQU 0 (
        echo [OK] Prospector Windows (.ps1)
        echo [OK] Prospector Windows (.ps1)>> "%LOG_FILE%"
    ) else (
        echo [ECHEC] Prospector Windows (.ps1) ^(code !RC_NUM!^)
        echo [ECHEC] Prospector Windows (.ps1) ^(code !RC_NUM!^)>> "%LOG_FILE%"
        set /a FAIL_COUNT+=1
        call :offer_wsl_prospector_fallback
    )
    exit /b 0
)

if exist "INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat" (
    call "INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat"
    set "RC=!errorlevel!"
    if "!RC!"=="0" (
        echo [OK] Prospector Windows (.bat)
        echo [OK] Prospector Windows (.bat)>> "%LOG_FILE%"
    ) else (
        echo [ECHEC] Prospector Windows (.bat) ^(code !RC!^)
        echo [ECHEC] Prospector Windows (.bat) ^(code !RC!^)>> "%LOG_FILE%"
        set /a FAIL_COUNT+=1
        call :offer_wsl_prospector_fallback
    )
    exit /b 0
)

echo [SKIP] INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1/.bat introuvable
echo [SKIP] Prospector Windows introuvable>> "%LOG_FILE%"
set /a MISS_COUNT+=1
exit /b 0

:install_sne_requirements
echo.
echo ------------------------------------------------------------
echo [OPTIONNEL] Profil SN Ia (requirements-cosmology-sne.txt)
echo ------------------------------------------------------------

call :ask_install_or_skip "Profil SN Ia (requirements-cosmology-sne.txt)"
if errorlevel 1 (
    echo [SKIP] Ignore par l'utilisateur: Profil SN Ia
    echo [SKIP] Ignore par l'utilisateur: Profil SN Ia>> "%LOG_FILE%"
    exit /b 0
)

echo [START] requirements-cosmology-sne.txt>> "%LOG_FILE%"

if not exist "requirements-cosmology-sne.txt" (
    echo [SKIP] requirements-cosmology-sne.txt introuvable
    echo [SKIP] requirements-cosmology-sne.txt introuvable>> "%LOG_FILE%"
    set /a MISS_COUNT+=1
    exit /b 0
)

set "CONDA_ROOT="
where conda >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    for /f "delims=" %%i in ('conda info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%ProgramData%\Miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)
if not defined CONDA_ROOT if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%USERPROFILE%\miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT=%%i"
)

if not defined CONDA_ROOT (
    echo [ECHEC] conda introuvable - pip SNE non lance
    echo [ECHEC] conda introuvable - pip SNE non lance>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
    exit /b 0
)

if not exist "%CONDA_ROOT%\Scripts\activate.bat" (
    echo [ECHEC] activate.bat introuvable sous %CONDA_ROOT%\Scripts
    echo [ECHEC] activate.bat introuvable sous %CONDA_ROOT%\Scripts>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
    exit /b 0
)

call "%CONDA_ROOT%\Scripts\activate.bat" astroenv
if errorlevel 1 (
    echo [ECHEC] activation astroenv impossible - pip SNE non lance
    echo [ECHEC] activation astroenv impossible - pip SNE non lance>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
    exit /b 0
)

python -m pip install --prefer-binary -r "requirements-cosmology-sne.txt"
set "RC=!errorlevel!"
if "!RC!"=="0" (
    echo [OK] requirements-cosmology-sne.txt installe
    echo [OK] requirements-cosmology-sne.txt installe>> "%LOG_FILE%"
) else (
    echo [ECHEC] pip install requirements-cosmology-sne.txt ^(code !RC!^)
    echo [ECHEC] pip install requirements-cosmology-sne.txt ^(code !RC!^)>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
)
exit /b 0

:wait_next_step
echo.
if "!AUTO_MODE!"=="1" (
    echo [AUTO] Continuation automatique de la sequence.
    exit /b 0
)
set "CONTINUE_OPTIONNELS="
set /p CONTINUE_OPTIONNELS="Continuer la sequence optionnelle ? (O/N): "
if /i "!CONTINUE_OPTIONNELS!"=="N" (
    echo Arret demande par l'utilisateur.
    echo [STOP] Sequence interrompue par l'utilisateur>> "%LOG_FILE%"
    exit /b 1
)
exit /b 0

:ask_install_or_skip
set "ASK_LABEL=%~1"
if "!AUTO_MODE!"=="1" (
    echo [AUTO] Installation de l'etape '!ASK_LABEL!'.
    exit /b 0
)
:ask_install_or_skip_loop
set "ASK_CHOICE="
set /p ASK_CHOICE="Installer l'etape '!ASK_LABEL!' ? (O/N): "
if "!ASK_CHOICE!"=="" exit /b 0
if /i "!ASK_CHOICE!"=="O" exit /b 0
if /i "!ASK_CHOICE!"=="N" exit /b 1
echo Reponse invalide. Tapez O pour installer ou N pour ignorer.
goto :ask_install_or_skip_loop

:offer_wsl_prospector_fallback
echo.
set "TRY_WSL_PROSPECTOR="
set /p TRY_WSL_PROSPECTOR="Echec Prospector Windows. Basculer maintenant vers Installation_fsps\prospector.bat ? (O/N): "
if /i not "!TRY_WSL_PROSPECTOR!"=="O" exit /b 0

if not exist "Installation_fsps\prospector.bat" (
    echo [SKIP] Installation_fsps\prospector.bat introuvable
    echo [SKIP] Installation_fsps\prospector.bat introuvable>> "%LOG_FILE%"
    set /a MISS_COUNT+=1
    exit /b 0
)

echo [START] Fallback Prospector WSL>> "%LOG_FILE%"
call "Installation_fsps\prospector.bat"
set "RC_FALLBACK=!errorlevel!"
if "!RC_FALLBACK!"=="0" (
    echo [OK] Fallback Prospector WSL termine
    echo [OK] Fallback Prospector WSL termine>> "%LOG_FILE%"
) else (
    echo [ECHEC] Fallback Prospector WSL ^(code !RC_FALLBACK!^)
    echo [ECHEC] Fallback Prospector WSL ^(code !RC_FALLBACK!^)>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
)
set "DONE_WSL_PROSPECTOR=1"
exit /b 0
