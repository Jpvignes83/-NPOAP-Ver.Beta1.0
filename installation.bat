@echo off
REM ===================================================================
REM NPOAP - Script d'installation automatique Windows
REM ===================================================================
REM Au double-clic, Windows lance souvent cmd en mode qui ferme la fenetre a la fin.
REM On rouvre dans une console /k pour garder la fenetre ouverte (logs, erreurs, pause).
if /i not "%~1"=="__NPOAP_INTERNAL__" (
    pushd "%~dp0" 2>nul
    start "Installation NPOAP" cmd.exe /k call "%~f0" __NPOAP_INTERNAL__
    popd 2>nul
    exit /b 0
)

setlocal enabledelayedexpansion

REM Toujours se placer dans le dossier du script ^(double-clic / start : evite xcopy depuis le mauvais repertoire^)
cd /d "%~dp0"

REM Pas de sequences ANSI : sortie texte brut uniquement (compatible toutes consoles)
set "GREEN="
set "RED="
set "YELLOW="
set "BLUE="
set "RESET="

REM Configuration
set VERSION=1.0
set MINICONDA_VERSION=latest
set ENV_NAME=astroenv
REM INSTALL_DIR : aucun chemin par defaut (chaque utilisateur choisit un dossier absolu)
set "INSTALL_DIR="

echo.
echo ============================================================
echo   NPOAP - Installation Automatique
echo   Version %VERSION%
echo ============================================================
echo.

REM ===================================================================
REM ACCORD D'UTILISATION
REM ===================================================================
echo.
echo ============================================================
echo   Setup V1.0
echo ============================================================
echo.
echo NPOAP Software Usage Agreement
echo ----------------------------------------
echo           (Rev.1.0)
echo.
echo As a user of the NPOAP software system (hereby referred to as "The System"^),
echo I agree that:
echo.
echo   1. Access to the system is restricted. As such, I shall not:
echo.
echo        A. Distribute the software to other users.
echo        B. Publish or disseminate the link to the software.
echo.
echo   2. The software is proprietary. As such I shall not change the source code
echo      except by written permission of the author, Jean-Pascal VIGNES.
echo.       jeanpascal.vignes@gmail.com
echo   3. The author shall not be held liable for damages that may occur while
echo      using the system.
echo.
echo ============================================================
echo.
:ASK_AGREEMENT
set /p ACCEPT_AGREEMENT="Do you accept the above agreement (y/n)?: "
if /i "%ACCEPT_AGREEMENT%"=="y" (
    echo.
    echo Agreement accepted. Proceeding with installation...
    echo.
    goto :AGREEMENT_ACCEPTED
) else if /i "%ACCEPT_AGREEMENT%"=="n" (
    echo.
    echo %RED%Installation cancelled. You must accept the agreement to proceed.%RESET%
    echo.
    pause
    exit /b 1
) else (
    echo.
    echo %YELLOW%Invalid input. Please enter 'y' for yes or 'n' for no.%RESET%
    echo.
    goto :ASK_AGREEMENT
)

:AGREEMENT_ACCEPTED

REM Vérifier les privilèges administrateur
net session >nul 2>&1
if %errorLevel% neq 0 (
    echo %RED%ERREUR: Ce script doit etre execute en tant qu'administrateur!%RESET%
    echo.
    echo Clic droit sur installation.bat ^> Executer en tant qu'administrateur
    pause
    exit /b 1
)

REM ===================================================================
REM ETAPE 1: Sélection du dossier d'installation
REM ===================================================================
echo %BLUE%=== ETAPE 1: Selection du dossier d'installation ===%RESET%
echo.
echo Indiquez un chemin absolu accessible en ecriture (ex. C:\Apps\NPOAP ou D:\Astro\NPOAP).
echo Aucun dossier personnel n'est impose par defaut.
echo Si les sources sont deja dans ce dossier ^(meme emplacement que installation.bat^), la copie sera
echo automatiquement ignoree ^(sinon xcopy signale une « copie cyclique »^).
echo.
:ASK_INSTALL_DIR
set "INSTALL_DIR_CUSTOM="
set /p INSTALL_DIR_CUSTOM="Chemin d'installation de NPOAP : "
if "!INSTALL_DIR_CUSTOM!"=="" (
    echo %RED%Le chemin ne peut pas etre vide. Saisissez un repertoire complet.%RESET%
    echo.
    goto :ASK_INSTALL_DIR
)
set "INSTALL_DIR=!INSTALL_DIR_CUSTOM!"
echo.
echo Dossier d'installation selectionne: !INSTALL_DIR!
echo.

if not exist "%INSTALL_DIR%" (
    mkdir "%INSTALL_DIR%"
    echo Dossier cree: %INSTALL_DIR%
) else (
    echo Dossier existe deja: %INSTALL_DIR%
)
echo.

REM ===================================================================
REM ETAPE 2: Vérification de la configuration système
REM ===================================================================
echo %BLUE%=== ETAPE 2: Verification de la configuration systeme ===%RESET%
echo.

REM MSVC, WSL, Ubuntu, Astrometry, KBMOD, CMake, Prospector : install_*.bat dans ce dossier
echo --- Controles generaux ---
echo Verification de l'architecture...
if /i "!PROCESSOR_ARCHITECTURE!"=="AMD64" (
    echo %GREEN%OK: Windows 64 bits (x64^).%RESET%
) else if /i "!PROCESSOR_ARCHITECTURE!"=="ARM64" (
    echo %YELLOW%ATTENTION: ARM64 - l'installateur Miniconda utilise la variante x86_64 ; verifiez la compatibilite.%RESET%
) else (
    echo %YELLOW%ATTENTION: Architecture !PROCESSOR_ARCHITECTURE! - une edition 64 bits est recommandee pour NPOAP.%RESET%
)
echo.

echo Test d'ecriture dans le dossier d'installation...
set "WRTEST=%INSTALL_DIR%\npoap_install_wr_test_%RANDOM%.tmp"
echo ok>"%WRTEST%" 2>nul
if exist "%WRTEST%" (
    del /f /q "%WRTEST%" >nul 2>&1
    echo %GREEN%OK: Le dossier est accessible en ecriture.%RESET%
) else (
    echo %RED%ERREUR: Impossible d'ecrire dans !INSTALL_DIR!%RESET%
    echo Verifiez les droits ou choisissez un autre chemin et relancez le script.
    pause
    exit /b 1
)
echo.

REM ===================================================================
REM ETAPE 3: Python (Miniconda + astroenv uniquement — pas d'installateur python.org)
REM ===================================================================
echo %BLUE%=== ETAPE 3: Python pour NPOAP ^(via Miniconda^) ===%RESET%
echo.
echo NPOAP s'execute dans l'environnement Conda "%ENV_NAME%" avec Python 3.11 ^(conda create^).
echo Un installateur separe depuis python.org n'est plus necessaire : Miniconda fournit deja
echo un interpreteur ^(base^) puis l'etape suivante cree %ENV_NAME% avec python=3.11.
echo.
echo Si vous souhaitez en plus un Python "systeme" hors Conda, installez-le manuellement depuis
echo https://www.python.org/ ^(optionnel, non requis pour ce script^).
echo.

REM ===================================================================
REM ETAPE 4: Installation de Miniconda
REM ===================================================================
echo %BLUE%=== ETAPE 4: Installation de Miniconda ===%RESET%
echo.

REM Detection conda : d'abord conda.exe dans Scripts ; puis where ^(priorite aux lignes se terminant par conda.exe^)
set "SKIP_CONDA="
set "CONDA_CMD="
set "CONDA_INST_ROOT="

if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "CONDA_CMD=%ProgramData%\Miniconda3\Scripts\conda.exe"
if not defined CONDA_CMD if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" set "CONDA_CMD=%USERPROFILE%\miniconda3\Scripts\conda.exe"
if not defined CONDA_CMD if exist "%LOCALAPPDATA%\miniconda3\Scripts\conda.exe" set "CONDA_CMD=%LOCALAPPDATA%\miniconda3\Scripts\conda.exe"
if not defined CONDA_CMD if exist "%ProgramData%\Anaconda3\Scripts\conda.exe" set "CONDA_CMD=%ProgramData%\Anaconda3\Scripts\conda.exe"
if not defined CONDA_CMD if exist "%USERPROFILE%\Anaconda3\Scripts\conda.exe" set "CONDA_CMD=%USERPROFILE%\Anaconda3\Scripts\conda.exe"
if not defined CONDA_CMD (
    where conda >nul 2>&1
    if !errorLevel! equ 0 (
        for /f "delims=" %%C in ('where conda 2^>nul') do (
            echo %%C | findstr /i /e "conda.exe" >nul
            if not errorlevel 1 (
                set "CONDA_CMD=%%C"
                goto :CONDA_CMD_DONE
            )
        )
        for /f "delims=" %%C in ('where conda 2^>nul') do (
            set "CONDA_CMD=%%C"
            goto :CONDA_CMD_DONE
        )
    )
)
:CONDA_CMD_DONE

if defined CONDA_CMD (
    echo !GREEN!Conda est deja present sur ce PC ^(detecte^).!RESET!
    echo Chemin : !CONDA_CMD!
    echo Version :
    set "CVER_OK=0"
    for /f "delims=" %%V in ('call "!CONDA_CMD!" --version 2^>nul') do (
        echo   %%V
        set "CVER_OK=1"
    )
    if "!CVER_OK!"=="0" echo !YELLOW!Impossible d'obtenir la version.!RESET!
    echo.
    echo O = telecharger et reinstaller Miniconda ^(installation silencieuse AllUsers sous ProgramData^).
    echo N ou Entree = conserver l'installation actuelle ^(pas de telechargement^).
    echo.
    set "REINSTALL_CONDA="
    set /p REINSTALL_CONDA="Votre choix ^(O/N^) : "
    if /i not "!REINSTALL_CONDA!"=="O" set "SKIP_CONDA=1"
)

if "!SKIP_CONDA!"=="1" if defined CONDA_CMD (
    for %%I in ("!CONDA_CMD!") do pushd "%%~dpI.." 2>nul
    if !errorLevel! equ 0 (
        set "CONDA_INST_ROOT=!CD!"
        popd
        set "PATH=!CONDA_INST_ROOT!;!CONDA_INST_ROOT!\Scripts;!CONDA_INST_ROOT!\Library\bin;%PATH%"
        echo Environnement : conda pour la suite du script depuis !CONDA_INST_ROOT!
        echo.
    ) else (
        for %%I in ("!CONDA_CMD!") do set "PATH=%%~dpI;%PATH%"
        echo %YELLOW%Impossible de trouver la racine Miniconda ; le dossier Scripts est ajoute au PATH.%RESET%
        echo.
    )
)

if not defined SKIP_CONDA (
    echo Telechargement de Miniconda...
    set CONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
    set "CONDA_INSTALLER=%TEMP%\Miniconda3-latest-Windows-x86_64.exe"
    if exist "!CONDA_INSTALLER!" del /f /q "!CONDA_INSTALLER!" >nul 2>&1
    
    REM Plusieurs methodes : BITS (sortie native Windows^), PowerShell + console, curl
    set "DL_OK=0"
    set "PSFULL=%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe"
    echo Vous devriez voir la progression ci-dessous (pourcentage BITS, barre PS, ou barre curl^).
    echo.
    where bitsadmin >nul 2>&1
    if !errorLevel! equ 0 (
        call set "BITSJOB=NPOAPmc%%RANDOM%%%%RANDOM%%"
        echo Methode A : BITS ^(bitsadmin^) - lignes de progression Windows...
        bitsadmin /transfer !BITSJOB! /download /priority NORMAL "!CONDA_URL!" "!CONDA_INSTALLER!"
        if !errorLevel! equ 0 if exist "!CONDA_INSTALLER!" set "DL_OK=1"
    )
    set "MCDLPS1=%~dp0installation_miniconda_download.ps1"
    if "!DL_OK!"=="0" if exist "!MCDLPS1!" (
        echo Methode B : PowerShell - barre ASCII + sortie rattachee a cette fenetre...
        if exist "!PSFULL!" (
            "!PSFULL!" -NoLogo -NoProfile -ExecutionPolicy Bypass -File "!MCDLPS1!" -Url "!CONDA_URL!" -Out "!CONDA_INSTALLER!"
        ) else (
            powershell -NoLogo -NoProfile -ExecutionPolicy Bypass -File "!MCDLPS1!" -Url "!CONDA_URL!" -Out "!CONDA_INSTALLER!"
        )
        if !errorLevel! equ 0 if exist "!CONDA_INSTALLER!" set "DL_OK=1"
    )
    if "!DL_OK!"=="0" if not exist "!MCDLPS1!" (
        echo %YELLOW%installation_miniconda_download.ps1 absent du dossier du .bat - methodes A/C seulement.%RESET%
    )
    if "!DL_OK!"=="0" (
        where curl >nul 2>&1
        if !errorLevel! equ 0 (
            echo Methode C : curl - barre sur stderr puis console...
            curl.exe --progress-bar -fSL --retry 3 --connect-timeout 30 -o "!CONDA_INSTALLER!" "!CONDA_URL!"
            if !errorLevel! equ 0 if exist "!CONDA_INSTALLER!" set "DL_OK=1"
        )
    )
    if not exist "!CONDA_INSTALLER!" set "DL_OK=0"
    REM Eviter les echecs silencieux (page HTML, fichier vide ; Miniconda fait typiquement ~80 Mo)
    if exist "!CONDA_INSTALLER!" for %%S in ("!CONDA_INSTALLER!") do (
        if %%~zS equ 0 set "DL_OK=0"
        if %%~zS LSS 5000000 (
            echo %RED%Le fichier telecharge est trop petit (%%~zS octets^) - echec probable du telechargement.%RESET%
            set "DL_OK=0"
        )
    )
    if "!DL_OK!"=="0" (
        if exist "!CONDA_INSTALLER!" del /f /q "!CONDA_INSTALLER!" >nul 2>&1
        echo %RED%ERREUR: Echec du telechargement de Miniconda.%RESET%
        echo Verifiez la connexion Internet, le proxy, ou telechargez l'installateur a la main depuis:
        echo   !CONDA_URL!
        pause
        exit /b 1
    )
    
    set "NSIS_RAN=0"
    set "NSIS_ERR=0"
    set "OFFER_CONDA_BASE_UPDATE=0"
    REM L'installateur NSIS AllUsers refuse un dossier deja rempli : eviter l'echec et le faux "succes".
    if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" (
        set "OFFER_CONDA_BASE_UPDATE=1"
        echo %YELLOW%Un Miniconda ^(tous utilisateurs^) est deja present :%RESET%
        echo   %ProgramData%\Miniconda3
        echo L'installateur silencieux ne peut pas reinstaller par-dessus ^(dossier non vide^).
        echo Desinstallation complete ^(Applications Windows^) : seulement si vous voulez tout refaire a zero.
        echo Sinon : apres l'etape des CGU des canaux, le script proposera une mise a jour de conda ^(environnement base^).
        echo.
        del /f /q "!CONDA_INSTALLER!" >nul 2>&1
    ) else (
        echo Installation de Miniconda en cours...
        echo %YELLOW%Ne fermez pas cette fenetre : cela peut prendre plusieurs minutes.%RESET%
        echo.
        echo %GREEN%Note:%RESET% Avec /InstallationType=AllUsers, /AddToPath=1 est ignore ^(comportement Anaconda^).
        echo.
        REM NSIS: l'option /D= doit etre en DERNIER. AllUsers -^> emplacement courant ProgramData.
        set "NSIS_RAN=1"
        "!CONDA_INSTALLER!" /InstallationType=AllUsers /RegisterPython=1 /S /D=%ProgramData%\Miniconda3
        set "NSIS_ERR=!errorLevel!"
        del /f /q "!CONDA_INSTALLER!" >nul 2>&1
        "%SystemRoot%\System32\timeout.exe" /t 3 /nobreak >nul
    )
    if "!NSIS_RAN!"=="1" if not "!NSIS_ERR!"=="0" if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "OFFER_CONDA_BASE_UPDATE=1"
    for /f "tokens=2*" %%A in ('reg query "HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment" /v PATH 2^>nul') do set "PATH=%%B"
    
    set "CONDA_EXE_OK=0"
    if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "CONDA_EXE_OK=1"
    if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" set "CONDA_EXE_OK=1"
    if "!CONDA_EXE_OK!"=="0" (
        echo %RED%ERREUR: Miniconda ne semble pas installe (conda.exe introuvable^).%RESET%
        if "!NSIS_RAN!"=="1" if not "!NSIS_ERR!"=="0" (
            echo L'installateur a peut-etre echoue ^(ex. ::error:: dossier non vide, droits^). Supprimez ou videz
            echo manuellement %ProgramData%\Miniconda3 en administrateur si besoin, ou desinstallez l'ancienne Miniconda.
        )
        echo Installez a la main en administrateur ou relancez le script.
        pause
        exit /b 1
    )
    if "!NSIS_RAN!"=="1" if not "!NSIS_ERR!"=="0" (
        echo %YELLOW%ATTENTION: L'installateur Miniconda a signale un echec ^(code !NSIS_ERR!^).%RESET%
        if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" (
            echo Une installation dans ProgramData est neanmoins detectee - le script continue avec celle-ci.
        ) else if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" (
            echo Conda est disponible depuis le profil utilisateur - le script continue.
        )
        echo.
    )
    
    if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "PATH=%ProgramData%\Miniconda3;%ProgramData%\Miniconda3\Scripts;%ProgramData%\Miniconda3\Library\bin;%PATH%"
    if exist "%ProgramData%\Anaconda3\Scripts\conda.exe" set "PATH=%ProgramData%\Anaconda3;%ProgramData%\Anaconda3\Scripts;%ProgramData%\Anaconda3\Library\bin;%PATH%"
    if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" set "PATH=%USERPROFILE%\miniconda3;%USERPROFILE%\miniconda3\Scripts;%USERPROFILE%\miniconda3\Library\bin;%PATH%"
    if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%ProgramData%\Miniconda3"
    if not defined CONDA_INST_ROOT if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%USERPROFILE%\miniconda3"
    
    conda --version >nul 2>&1
    if !errorLevel! equ 0 (
        if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" if "!NSIS_RAN!"=="0" (
            echo %GREEN%Conda : installation existante dans ProgramData conservee.%RESET%
        ) else if "!NSIS_RAN!"=="1" if "!NSIS_ERR!"=="0" (
            echo %GREEN%Miniconda installe avec succes!%RESET%
        ) else (
            echo %GREEN%Conda utilisable ^(version ci-dessous^).%RESET%
        )
        conda --version
    ) else (
        echo %YELLOW%ATTENTION: Conda n'est pas encore dans le PATH de cette session.%RESET%
        echo Reessayez: fermez cette fenetre puis relancez installation.bat.
        set CONDA_WARNING=1
    )
) else (
    echo %YELLOW%ETAPE 4 : telechargement Miniconda IGNORE - conda est deja detecte.%RESET%
    echo Pour afficher la progression, relancez et repondez O pour reinstaller Miniconda.
)
echo.

REM Racine installation conda : obligatoire pour activate.bat dans cmd ^(sans conda init^)
if not defined CONDA_INST_ROOT (
    if exist "%ProgramData%\Miniconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%ProgramData%\Miniconda3"
    if not defined CONDA_INST_ROOT if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%USERPROFILE%\miniconda3"
    if not defined CONDA_INST_ROOT if exist "%LOCALAPPDATA%\miniconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%LOCALAPPDATA%\miniconda3"
    if not defined CONDA_INST_ROOT if exist "%ProgramData%\Anaconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%ProgramData%\Anaconda3"
    if not defined CONDA_INST_ROOT if exist "%USERPROFILE%\Anaconda3\Scripts\conda.exe" set "CONDA_INST_ROOT=%USERPROFILE%\Anaconda3"
)
if not defined CONDA_INST_ROOT (
    echo %RED%ERREUR: impossible de determiner la racine conda ^(Scripts\conda.exe introuvable^).%RESET%
    pause
    exit /b 1
)

REM Conda recent : sans acceptation explicite des CGU des canaux Anaconda par defaut,
REM "conda create" en script echoue avec CondaToSNonInteractiveError.
echo Acceptation des CGU des canaux repo.anaconda.com ^(obligatoire pour conda create en mode script^)...
call conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>nul
call conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>nul
call conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2 2>nul
echo.

REM Mise a jour de conda dans "base" sans reinstaller l'installateur .exe (dossier deja rempli, etc.)
if "!OFFER_CONDA_BASE_UPDATE!"=="1" (
    echo %BLUE%--- Mise a jour optionnelle : environnement de base Conda ---%RESET%
    echo Met a jour l'outil %GREEN%conda%RESET% et les paquets du noyau ^(environnement %GREEN%base%RESET%^), sans reinstaller Miniconda.
    set "DUCMD="
    set /p DUCMD="Lancer : conda update -n base -c defaults conda -y ? (O/N) : "
    if /i "!DUCMD!"=="O" (
        echo.
        "%CONDA_INST_ROOT%\Scripts\conda.exe" update -n base -c defaults conda -y
        if !errorLevel! neq 0 (
            echo %YELLOW%conda update a echoue. Essayez en administrateur, ou en ligne de commande :%RESET%
            echo   call "%CONDA_INST_ROOT%\Scripts\activate.bat" base
            echo   conda update -n base -c defaults conda
        ) else (
            echo %GREEN%Mise a jour de conda ^(base^) terminee.%RESET%
        )
        echo.
    ) else (
        echo Mise a jour de base ignoree. Pour la faire plus tard : conda activate base ^& conda update -n base -c defaults conda
        echo.
    )
)

REM ===================================================================
REM ETAPE 5: Création de l'environnement Conda
REM ===================================================================
echo %BLUE%=== ETAPE 5: Creation de l'environnement Conda ===%RESET%
echo.

conda env list | findstr /C:"%ENV_NAME%" >nul 2>&1
if %errorLevel% equ 0 (
    echo L'environnement %ENV_NAME% existe deja.
    echo.
    set /p RECREATE_ENV="Voulez-vous le recreer? (O/N): "
    if /i "!RECREATE_ENV!"=="O" (
        echo Suppression de l'environnement existant...
        call conda env remove -n %ENV_NAME% -y
        echo Creation du nouvel environnement...
        call conda create -n %ENV_NAME% python=3.11 -y
    ) else (
        echo Utilisation de l'environnement existant.
    )
) else (
    echo Creation de l'environnement %ENV_NAME%...
    call conda create -n %ENV_NAME% python=3.11 -y
)

if %errorLevel% neq 0 (
    echo %RED%ERREUR: Echec de la creation de l'environnement%RESET%
    pause
    exit /b 1
)

echo %GREEN%Environnement %ENV_NAME% cree avec succes!%RESET%
echo.

REM ===================================================================
REM ETAPE 6: Copie du projet, normalisation UTF-8, pip ^(cœur PyPI uniquement^)
REM ===================================================================
echo %BLUE%=== ETAPE 6: Copie NPOAP, requirements UTF-8, installation pip ===%RESET%
echo.
echo Dependances installees depuis requirements_install_core.txt ^(PyPI seul, pas de git / ezpadova / PHOEBE / optionnels lourds^).
echo Pour tout requirements.txt ^(git+, etc.^): conda activate %ENV_NAME% ^& pip install -r requirements.txt
echo.

if not exist "requirements_install_core.txt" (
    echo %RED%ERREUR: requirements_install_core.txt introuvable dans le dossier du script!%RESET%
    pause
    exit /b 1
)
if not exist "requirements.txt" (
    echo %RED%ERREUR: requirements.txt introuvable dans le dossier du script!%RESET%
    pause
    exit /b 1
)

REM Liste d'exclusion avant xcopy ^(premier passage^)
if not exist exclude_files.txt (
    (
        echo __pycache__
        echo *.pyc
        echo .git
        echo logs\*.log
        echo *.log
        echo .conda
        echo .cache
    ) > exclude_files.txt
)

REM Reprise d'installation : retirer lecture seule ^(etape 9 ou protection precedente^), sinon ecriture requirements impossible
if exist "!INSTALL_DIR!" (
    attrib -R "!INSTALL_DIR!\*.*" /S >nul 2>&1
)

set "SKIP_PROJECT_COPY="
set "SRC_CMP=%~dp0"
if "!SRC_CMP:~-1!"=="\" set "SRC_CMP=!SRC_CMP:~0,-1!"
set "DST_CMP=!INSTALL_DIR!"
if "!DST_CMP:~-1!"=="\" set "DST_CMP=!DST_CMP:~0,-1!"
if /i "!SRC_CMP!"=="!DST_CMP!" set "SKIP_PROJECT_COPY=1"

if defined SKIP_PROJECT_COPY (
    echo Dossier cible identique au dossier des sources - copie ignoree ^(evite « copie cyclique »^).
) else (
    echo Copie des fichiers vers !INSTALL_DIR!...
    xcopy /E /I /Y /Q . "!INSTALL_DIR!" /EXCLUDE:exclude_files.txt
    set "XCOPY_EC=!errorLevel!"
    REM xcopy : 0 ok ; 1 aucun fichier a copier ; 2 arret utilisateur ; 4 erreur init ^(chemins, memoire^)
    if !XCOPY_EC! geq 4 (
        echo %YELLOW%ATTENTION: xcopy code !XCOPY_EC! - tentative robocopy...%RESET%
        robocopy "%~dp0." "!INSTALL_DIR!" /E /COPY:DAT /R:2 /W:3 /XD __pycache__ .git .conda .cache /XF *.pyc /NFL /NDL /NJH /NJS /NC /NS /NP
        if errorlevel 8 (
            echo %RED%ERREUR: copie vers !INSTALL_DIR! echouee ^(xcopy code !XCOPY_EC!, robocopy ^>=8^).%RESET%
        )
    ) else if !XCOPY_EC! equ 2 (
        echo %YELLOW%ATTENTION: xcopy interrompu ^(code 2^).%RESET%
    )
)

attrib -R "!INSTALL_DIR!\requirements.txt" >nul 2>&1
attrib -R "!INSTALL_DIR!\requirements_install_core.txt" >nul 2>&1

REM pip exige des fichiers requirements en UTF-8 ^(pas UTF-16 / Unicode Bloc-notes^)
set "NORM_PS=%~dp0normalize_requirements_utf8.ps1"
if exist "!NORM_PS!" (
    echo Normalisation UTF-8 des fichiers requirements pour pip...
    powershell -NoProfile -ExecutionPolicy Bypass -File "!NORM_PS!" "!INSTALL_DIR!\requirements_install_core.txt"
    powershell -NoProfile -ExecutionPolicy Bypass -File "!NORM_PS!" "!INSTALL_DIR!\requirements.txt"
) else (
    echo %YELLOW%ATTENTION: normalize_requirements_utf8.ps1 absent - encodage requirements non corrige.%RESET%
)

echo Activation de l'environnement %ENV_NAME% ^(activate.bat^)...
call "!CONDA_INST_ROOT!\Scripts\activate.bat" %ENV_NAME%

if !errorLevel! neq 0 (
    echo %RED%ERREUR: Impossible d'activer l'environnement%RESET%
    pause
    exit /b 1
)

echo Mise a jour de pip, setuptools, wheel...
python -m pip install --upgrade pip setuptools wheel

echo Installation des dependances ^(cœur PyPI^) depuis requirements_install_core.txt...
python -m pip install --prefer-binary -r "!INSTALL_DIR!\requirements_install_core.txt"

if !errorLevel! neq 0 (
    echo %YELLOW%ATTENTION: pip a signale une erreur. Verifiez les lignes ci-dessus ^(reseau, proxy, antivirus^).%RESET%
    echo Pour reessayer : conda activate %ENV_NAME% ^& cd /d "!INSTALL_DIR!" ^& python -m pip install --prefer-binary -r requirements_install_core.txt
) else (
    echo %GREEN%Dependances cœur installees avec succes!%RESET%
    echo %YELLOW%Optionnel : PHOEBE, ezpadova ^(git^), lightkurve, pip complet ^>^> pip install -r requirements.txt ^(Git requis^)%RESET%
)
echo.

REM ===================================================================
REM ETAPE 7: Creation des scripts de lancement
REM ===================================================================
echo %BLUE%=== ETAPE 7: Creation des scripts de lancement ===%RESET%
echo.

REM Script de lancement principal : LANCEMENT.bat (fourni a cote de installation.bat^)
if exist "%~dp0LANCEMENT.bat" (
    copy /Y "%~dp0LANCEMENT.bat" "%INSTALL_DIR%\LANCEMENT.bat" >nul
    echo %GREEN%Script de lancement copie: %INSTALL_DIR%\LANCEMENT.bat%RESET%
) else (
    echo %YELLOW%LANCEMENT.bat absent du package - generation d'un lanceur minimal vers %ENV_NAME%...%RESET%
    (
        echo @echo off
        echo call "!CONDA_INST_ROOT!\Scripts\activate.bat" %ENV_NAME%
        echo cd /d "%INSTALL_DIR%"
        echo python main.py
        echo pause
    ) > "%INSTALL_DIR%\LANCEMENT.bat"
    echo %GREEN%Script de lancement cree: %INSTALL_DIR%\LANCEMENT.bat%RESET%
)
echo.

REM ===================================================================
REM ETAPE 8: Test de l'installation
REM ===================================================================
echo %BLUE%=== ETAPE 8: Test de l'installation ===%RESET%
echo.

call "!CONDA_INST_ROOT!\Scripts\activate.bat" %ENV_NAME%
cd /d "%INSTALL_DIR%"

if exist "test_installation.py" (
    echo Execution des tests...
    python test_installation.py
    echo.
) else (
    echo %YELLOW%test_installation.py introuvable, test ignore%RESET%
)

REM ===================================================================
REM ETAPE 9: Protection des fichiers source
REM ===================================================================
echo %BLUE%=== ETAPE 9: Protection des fichiers source ===%RESET%
echo.

echo Protection des fichiers source contre les modifications...
if exist "protect_files.bat" (
    call protect_files.bat
) else (
    echo Mise en lecture seule des fichiers source...
    attrib +R "*.py" >nul 2>&1
    attrib +R "*.md" >nul 2>&1
    REM requirements*.txt non verrouilles : reinstall / pip / normalisation UTF-8
    attrib +R "core\*.*" /S >nul 2>&1
    attrib +R "gui\*.*" /S >nul 2>&1
    attrib +R "utils\*.*" /S >nul 2>&1
    attrib +R "docs\*.*" /S >nul 2>&1
    REM Ne pas protéger les fichiers nécessaires à l'exécution
    attrib -R "config.py" >nul 2>&1
    attrib -R "config.json" >nul 2>&1
    attrib -R "requirements.txt" >nul 2>&1
    attrib -R "requirements_install_core.txt" >nul 2>&1
    if exist "logs" attrib -R "logs\*.*" /S >nul 2>&1
    if exist "output" attrib -R "output\*.*" /S >nul 2>&1
    echo %GREEN%Fichiers proteges en lecture seule%RESET%
    echo %YELLOW%Note: Cette protection peut etre retiree par l'utilisateur%RESET%
)

echo.
echo ============================================================
echo   Installation terminee!
echo ============================================================
echo.
echo Pour lancer NPOAP, utilisez:
echo   %INSTALL_DIR%\LANCEMENT.bat
echo.
echo Ou depuis la ligne de commande:
echo   call "!CONDA_INST_ROOT!\Scripts\activate.bat" %ENV_NAME%
echo   cd %INSTALL_DIR%
echo   python main.py
echo.
echo ------------------------------------------------------------
echo   Composants optionnels (scripts dans le dossier NPOAP)
echo ------------------------------------------------------------
echo   Installation auto de tous les optionnels : install_optionnels.bat
echo   Visual C++ Build Tools / MSVC .... install_msvc_build_tools.bat
echo   CMake .......................... install_cmake.bat
echo   WSL ............................ install_wsl.bat
echo   Ubuntu (distro WSL) ............ install_ubuntu_wsl.bat
echo   Astrometry.net dans WSL ........ install_astrometry_wsl.bat
echo   KBMOD sous WSL/Linux ........... install_kbmod_wsl.bat
echo   Prospector (astroenv Windows) .. INSTALLER_PROSPECTOR_COMPLET_WINDOWS.bat
echo   idem (PowerShell + options) .... INSTALLER_PROSPECTOR_COMPLET_WINDOWS.ps1
echo   Prospector + FSPS via WSL ...... Installation_fsps\prospector.bat
echo   SORA seul (reinstall) .......... INSTALLER_SORA_ASTROENV.bat  (deja dans requirements_install_core.txt)
echo.
echo Ouvrez le fichier texte LISTE_INSTALL_OPTIONNELS.txt pour le meme tableau.
echo.
set "RUN_OPTIONNELS="
set /p RUN_OPTIONNELS="Lancer install_optionnels.bat maintenant pour executer toute la sequence optionnelle ? (O/N): "
if /i "!RUN_OPTIONNELS!"=="O" (
    if exist "install_optionnels.bat" (
        echo.
        echo Lancement de install_optionnels.bat...
        call "install_optionnels.bat"
    ) else (
        echo.
        echo %YELLOW%install_optionnels.bat introuvable dans %CD%.%RESET%
    )
) else (
    echo Composants optionnels ignores pour cette execution.
)
echo.
pause

