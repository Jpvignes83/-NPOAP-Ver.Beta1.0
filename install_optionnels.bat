@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul

REM NPOAP - Installation automatique des composants optionnels
REM A lancer depuis le dossier NPOAP apres installation.bat

cd /d "%~dp0"

set "FAIL_COUNT=0"
set "MISS_COUNT=0"
set "AUTO_MODE=0"
if /i "%NPOAP_AUTO%"=="1" set "AUTO_MODE=1"
set "LOG_FILE=%~dp0install_optionnels.log"

echo ============================================================ > "%LOG_FILE%"
echo NPOAP - Installation optionnelle (demarrage: %date% %time%) >> "%LOG_FILE%"
echo Dossier: %CD% >> "%LOG_FILE%"
echo ============================================================ >> "%LOG_FILE%"

if not defined INSTALL_COMPLET_LOG set "INSTALL_COMPLET_LOG=%~dp0install_complet.log"

echo.
echo ============================================================
echo NPOAP - Installation des composants optionnels
echo ============================================================
echo Dossier: %CD%
echo Journal optionnel: %LOG_FILE%
echo Journal complet ^(installation + optionnels^): !INSTALL_COMPLET_LOG!
if "!AUTO_MODE!"=="1" echo Mode auto: OUI ^(NPOAP_AUTO=1^)
echo.
echo Ce script lance automatiquement les installateurs optionnels.
echo Certains composants peuvent demander des droits administrateur ou un redemarrage.
echo.
echo Profil PyPI unifie : requirements_install_optionnels.txt ^(PHOEBE, rebound, ultranest, SORA, SNe Ia, STDPipe+GP^).
echo Si vous voyez encore "requirements-cosmology-sne.txt", mettez a jour ce .bat depuis le depot.
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
call :run_bat_step "KBMOD sous WSL/Linux (installation auto)" "install_kbmod_wsl.bat"
call :wait_next_step
if errorlevel 1 goto :end_sequence

call :install_gfortran_windows
call :wait_next_step
if errorlevel 1 goto :end_sequence

REM Prospector + FSPS : uniquement manuel ^(Installation_fsps\prospector.bat en WSL^), pas d'etape Windows dans cette sequence.

REM SORA ^(sora-astro^) : installe / met a jour avec les autres paquets PyPI ^(requirements_install_optionnels.txt^).
REM Mise a jour SORA seule : double-cliquez INSTALLER_SORA_ASTROENV.bat ou pip install -U sora-astro dans astroenv.

call :install_pip_optionnels
call :wait_next_step
if errorlevel 1 goto :end_sequence

:end_sequence

echo.>> "%LOG_FILE%"
echo ============================================================>> "%LOG_FILE%"
echo Fin: %date% %time% >> "%LOG_FILE%"
echo Echecs: !FAIL_COUNT! ; Fichiers manquants: !MISS_COUNT! >> "%LOG_FILE%"
echo ============================================================>> "%LOG_FILE%"

if not defined INSTALL_COMPLET_LOG set "INSTALL_COMPLET_LOG=%~dp0install_complet.log"
echo.>> "!INSTALL_COMPLET_LOG!"
echo ----- Contenu integral install_optionnels.log ^(%date% %time%^) ----- >> "!INSTALL_COMPLET_LOG!"
type "%LOG_FILE%" >> "!INSTALL_COMPLET_LOG!"
echo.>> "!INSTALL_COMPLET_LOG!"
echo Fin suite optionnelle dans install_complet.log: %date% %time% >> "!INSTALL_COMPLET_LOG!"

echo.
echo ============================================================
echo Installation optionnelle terminee
echo ============================================================
echo Echecs: !FAIL_COUNT!
echo Fichiers manquants: !MISS_COUNT!
echo Journal optionnel: %LOG_FILE%
echo Journal complet: !INSTALL_COMPLET_LOG!
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

:install_pip_optionnels
echo.
echo ------------------------------------------------------------
echo [OPTIONNEL] Paquets Python (requirements_install_optionnels.txt)
echo ------------------------------------------------------------

call :ask_install_or_skip "Paquets optionnels PyPI (requirements_install_optionnels.txt)"
if errorlevel 1 (
    echo [SKIP] Ignore par l'utilisateur: requirements_install_optionnels
    echo [SKIP] Ignore par l'utilisateur: requirements_install_optionnels>> "%LOG_FILE%"
    exit /b 0
)

echo [START] requirements_install_optionnels.txt>> "%LOG_FILE%"

if not exist "requirements_install_optionnels.txt" (
    echo [SKIP] requirements_install_optionnels.txt introuvable
    echo [SKIP] requirements_install_optionnels.txt introuvable>> "%LOG_FILE%"
    set /a MISS_COUNT+=1
    exit /b 0
)

set "CONDA_ROOT_OPT="
where conda >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    for /f "delims=" %%i in ('conda info --base 2^>nul') do set "CONDA_ROOT_OPT=%%i"
)
if not defined CONDA_ROOT_OPT if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%ProgramData%\Miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT_OPT=%%i"
)
if not defined CONDA_ROOT_OPT if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" (
    for /f "delims=" %%i in ('"%USERPROFILE%\miniconda3\Scripts\conda.exe" info --base 2^>nul') do set "CONDA_ROOT_OPT=%%i"
)

if not defined CONDA_ROOT_OPT (
    echo [ECHEC] conda introuvable - pip optionnels non lance
    echo [ECHEC] conda introuvable - pip optionnels non lance>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
    exit /b 0
)

if not exist "%CONDA_ROOT_OPT%\Scripts\activate.bat" (
    echo [ECHEC] activate.bat introuvable sous %CONDA_ROOT_OPT%\Scripts
    echo [ECHEC] activate.bat introuvable sous %CONDA_ROOT_OPT%\Scripts>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
    exit /b 0
)

call "%CONDA_ROOT_OPT%\Scripts\activate.bat" astroenv
if errorlevel 1 (
    echo [ECHEC] activation astroenv impossible - pip optionnels non lance
    echo [ECHEC] activation astroenv impossible - pip optionnels non lance>> "%LOG_FILE%"
    set /a FAIL_COUNT+=1
    exit /b 0
)

python -m pip install --upgrade "setuptools<81" 1> "%TEMP%\npoap_setuptools_fix.log" 2>&1
type "%TEMP%\npoap_setuptools_fix.log"
echo.>> "%LOG_FILE%"
echo ----- pip install --upgrade setuptools^<81 ----- >> "%LOG_FILE%"
type "%TEMP%\npoap_setuptools_fix.log" >> "%LOG_FILE%"

python -m pip install --prefer-binary -r "requirements_install_optionnels.txt" 1> "%TEMP%\npoap_pip_optionnels.log" 2>&1
set "RC_OPT=!errorlevel!"
type "%TEMP%\npoap_pip_optionnels.log"
echo.>> "%LOG_FILE%"
echo ----- pip install -r requirements_install_optionnels.txt ----- >> "%LOG_FILE%"
type "%TEMP%\npoap_pip_optionnels.log" >> "%LOG_FILE%"
if "!RC_OPT!"=="0" (
    echo [OK] requirements_install_optionnels.txt installe
    echo [OK] requirements_install_optionnels.txt installe>> "%LOG_FILE%"
    python -c "import stdpipe" 1>> "%TEMP%\npoap_stdpipe_check.log" 2>&1
    if errorlevel 1 (
        echo [WARN] stdpipe non importable apres pip - tentative: pip install --prefer-binary stdpipe
        echo [WARN] stdpipe non importable - relance pip install stdpipe>> "%LOG_FILE%"
        python -m pip install --prefer-binary stdpipe 1>> "%TEMP%\npoap_pip_stdpipe_retry.log" 2>&1
        type "%TEMP%\npoap_pip_stdpipe_retry.log"
        echo.>> "%LOG_FILE%"
        echo ----- pip install stdpipe ^(retry^) ----- >> "%LOG_FILE%"
        type "%TEMP%\npoap_pip_stdpipe_retry.log" >> "%LOG_FILE%"
        python -c "import stdpipe; print('stdpipe OK')" 1>> "%TEMP%\npoap_stdpipe_check2.log" 2>&1
        if errorlevel 1 (
            echo [ECHEC] stdpipe toujours non importable - voir logs ci-dessus et docs STDPipe
            echo [ECHEC] stdpipe import impossible apres retry>> "%LOG_FILE%"
            set /a FAIL_COUNT+=1
        ) else (
            echo [OK] stdpipe importable apres retry
            echo [OK] stdpipe importable apres retry>> "%LOG_FILE%"
        )
    )
) else (
    echo [ECHEC] pip install requirements_install_optionnels.txt ^(code !RC_OPT!^)
    echo [ECHEC] pip install requirements_install_optionnels.txt ^(code !RC_OPT!^)>> "%LOG_FILE%"
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
REM Supprimer espaces debut/fin (sinon " N " ou "non " etait refuse).
for /f "tokens=*" %%a in ("!ASK_CHOICE!") do set "ASK_CHOICE=%%a"
if not defined ASK_CHOICE (
    echo Reponse vide : tapez O ou oui pour installer, N ou non pour ignorer ^(Entree seule ne fait rien^).
    goto :ask_install_or_skip_loop
)
if /i "!ASK_CHOICE!"=="O" exit /b 0
if /i "!ASK_CHOICE!"=="OUI" exit /b 0
if /i "!ASK_CHOICE!"=="Y" exit /b 0
if /i "!ASK_CHOICE!"=="YES" exit /b 0
if /i "!ASK_CHOICE!"=="N" exit /b 1
if /i "!ASK_CHOICE!"=="NON" exit /b 1
if /i "!ASK_CHOICE!"=="NO" exit /b 1
echo Reponse invalide ^(!ASK_CHOICE!^). Tapez O / oui ou N / non.
goto :ask_install_or_skip_loop
