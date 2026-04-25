# ===================================================================
# NPOAP - Installation complète de Prospector pour Windows
# ===================================================================
# Ce script installe Prospector avec toutes ses dépendances
# Basé sur l'analyse complète des logs d'installation WSL
# ===================================================================

param(
    [string]$CondaEnv = "astroenv",
    [switch]$InstallFSPS = $false,
    [switch]$SkipVerification = $false,
    [switch]$ForceReinstall = $false
)

$ErrorActionPreference = "Stop"

# Fonctions utilitaires pour les messages colorés
function Write-Info { 
    param([string]$Message)
    Write-Host $Message -ForegroundColor Cyan 
}

function Write-Success { 
    param([string]$Message)
    Write-Host $Message -ForegroundColor Green 
}

function Write-Warning { 
    param([string]$Message)
    Write-Host $Message -ForegroundColor Yellow 
}

function Write-ErrorMsg { 
    param([string]$Message)
    Write-Host $Message -ForegroundColor Red 
}

function Write-Step {
    param([int]$Step, [int]$Total, [string]$Message)
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Blue
    Write-Host "[$Step/$Total] $Message" -ForegroundColor Blue
    Write-Host "========================================" -ForegroundColor Blue
    Write-Host ""
}

# En-tête
Write-Host ""
Write-Host "========================================" -ForegroundColor Blue
Write-Host "Installation complète: Prospector + FSPS" -ForegroundColor Blue
Write-Host "========================================" -ForegroundColor Blue
Write-Host ""
Write-Host "Ce script installe Prospector avec toutes ses dépendances"
Write-Host "Basé sur l'analyse des logs d'installation réussis"
Write-Host ""

# ===================================================================
# ÉTAPE 1: Vérification des prérequis
# ===================================================================
Write-Step -Step 1 -Total 7 -Message "Vérification des prérequis"

# Vérifier Conda (PATH ou emplacements standards Miniconda/Anaconda)
$CondaExe = $null
$condaPath = Get-Command conda -ErrorAction SilentlyContinue
if ($condaPath) {
    $CondaExe = $condaPath.Source
}

if (-not $CondaExe) {
    $candidates = @(
        "$env:ProgramData\Miniconda3\Scripts\conda.exe",
        "$env:USERPROFILE\miniconda3\Scripts\conda.exe",
        "$env:LOCALAPPDATA\miniconda3\Scripts\conda.exe",
        "$env:ProgramData\Anaconda3\Scripts\conda.exe",
        "$env:USERPROFILE\Anaconda3\Scripts\conda.exe"
    )
    foreach ($cand in $candidates) {
        if (Test-Path -LiteralPath $cand) {
            $CondaExe = $cand
            break
        }
    }
}

if (-not $CondaExe) {
    Write-ErrorMsg "ERREUR: Conda n'est pas installé ou n'est pas dans le PATH"
    Write-ErrorMsg "Installez Miniconda depuis: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
}
Write-Success "✓ Conda trouvé: $CondaExe"

# Vérifier Git (requis pour installer depuis GitHub)
$gitPath = Get-Command git -ErrorAction SilentlyContinue
if (-not $gitPath) {
    Write-ErrorMsg "ERREUR: Git n'est pas installé ou n'est pas dans le PATH"
    Write-ErrorMsg "Installez Git depuis: https://git-scm.com/download/win"
    Write-Host ""
    Write-Host "Après l'installation de Git, redémarrez ce script."
    exit 1
}
Write-Success "✓ Git trouvé: $($gitPath.Source)"

# Vérifier l'environnement Conda
Write-Info "Vérification de l'environnement Conda: $CondaEnv"
$envList = & $CondaExe env list 2>$null
if ($envList -match $CondaEnv) {
    Write-Success "✓ Environnement '$CondaEnv' existe"
} else {
    Write-Warning "L'environnement '$CondaEnv' n'existe pas. Création..."
    & $CondaExe create -n $CondaEnv python=3.11 -y
    if ($LASTEXITCODE -ne 0) {
        Write-ErrorMsg "ERREUR: Impossible de créer l'environnement Conda"
        exit 1
    }
    Write-Success "✓ Environnement '$CondaEnv' créé"
}

# Obtenir le chemin Python de l'environnement
# conda run peut renvoyer vide ou plusieurs lignes ; Test-Path refuse une chaine vide.
$rawPython = & $CondaExe run -n $CondaEnv python -c "import sys; print(sys.executable)" 2>$null
$pythonPath = $null
if ($null -ne $rawPython) {
    if ($rawPython -is [array]) {
        $pythonPath = $rawPython[-1].ToString().Trim()
    } else {
        $pythonPath = $rawPython.ToString().Trim()
    }
}

$pythonOk = $false
if (-not [string]::IsNullOrWhiteSpace($pythonPath)) {
    $pythonOk = Test-Path -LiteralPath $pythonPath
}
if (-not $pythonOk) {
    $condaBase = ((& $CondaExe info --base) | Out-String).Trim()
    if ([string]::IsNullOrWhiteSpace($condaBase)) {
        Write-ErrorMsg "ERREUR: conda info --base n'a pas renvoye de chemin"
        exit 1
    }
    $pythonPath = Join-Path $condaBase "envs\$CondaEnv\python.exe"
    if (-not (Test-Path -LiteralPath $pythonPath)) {
        Write-ErrorMsg "ERREUR: Impossible de trouver Python dans l'environnement '$CondaEnv'"
        Write-ErrorMsg "  Attendu: $pythonPath"
        exit 1
    }
}
Write-Info "Python: $pythonPath"

# Echecs d'import Python -> stderr. Sous Windows PowerShell 5.x, avec
# $ErrorActionPreference Stop, ce stderr devient NativeCommandError et
# interrompt le script meme avec 2>$null.
function Invoke-PythonCodeQuiet {
    param([Parameter(Mandatory)][string]$Code)
    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'SilentlyContinue'
    try {
        $null = & $pythonPath -c $Code 2>&1 | Out-Null
        return $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $prev
    }
}

function Invoke-PythonCodeOutput {
    param([Parameter(Mandatory)][string]$Code)
    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'SilentlyContinue'
    try {
        $chunks = @(& $pythonPath -c $Code 2>&1)
        $exit = $LASTEXITCODE
        $outLines = foreach ($c in $chunks) {
            if ($c -isnot [System.Management.Automation.ErrorRecord]) { [string]$c }
        }
        $errLines = foreach ($c in $chunks) {
            if ($c -is [System.Management.Automation.ErrorRecord]) { $c.ToString() }
        }
        return [pscustomobject]@{
            ExitCode   = $exit
            Stdout     = ($outLines -join "`n").Trim()
            StderrText = ($errLines -join "`n").Trim()
        }
    } finally {
        $ErrorActionPreference = $prev
    }
}

# ===================================================================
# ÉTAPE 2: Installation des dépendances Python de base
# ===================================================================
Write-Step -Step 2 -Total 7 -Message "Installation des dépendances Python de base"

$baseDeps = @(
    "numpy>=1.20.0",
    "scipy>=1.7.0",
    "pandas>=1.3.0",
    "astropy>=5.0.0"
)

foreach ($dep in $baseDeps) {
    Write-Info "  Installation de $dep..."
    & $pythonPath -m pip install $dep --quiet --upgrade
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "  ⚠ Erreur lors de l'installation de $dep"
    } else {
        Write-Success "  ✓ $dep installé"
    }
}

# ===================================================================
# ÉTAPE 3: Installation de sedpy depuis GitHub
# ===================================================================
Write-Step -Step 3 -Total 7 -Message "Installation de sedpy depuis GitHub"
Write-Warning "  IMPORTANT: sedpy doit être installé depuis GitHub (pas PyPI)"
Write-Warning "  La version PyPI n'a pas le module 'observate' requis par Prospector"

# Desinstaller sedpy PyPI si present (import Python valide sur une ligne)
Write-Info "  Verification de l'installation actuelle de sedpy..."
$sedpyCheck = "NOK"
if ((Invoke-PythonCodeQuiet -Code "import sedpy") -eq 0) { $sedpyCheck = "OK" }
if ($sedpyCheck -eq "OK") {
    # Eviter le faux echec PowerShell quand pip ecrit un warning sur stderr.
    $pipHasSedpy = Invoke-PythonCodeQuiet -Code "import importlib.metadata as m; m.version('sedpy'); print('OK')"
    if ($pipHasSedpy -eq 0) {
        Write-Info "  Desinstallation de sedpy (PyPI)..."
        $prevEA = $ErrorActionPreference
        $ErrorActionPreference = 'SilentlyContinue'
        try {
            & $pythonPath -m pip uninstall sedpy -y --quiet 2>$null
            $null = $LASTEXITCODE
        } finally {
            $ErrorActionPreference = $prevEA
        }
    } else {
        Write-Info "  sedpy importe mais n'est pas installe via pip dans cet environnement (pas de desinstallation)."
    }
}

# Installer depuis GitHub
Write-Info "  Installation depuis GitHub: git+https://github.com/bd-j/sedpy.git"
& $pythonPath -m pip install "git+https://github.com/bd-j/sedpy.git" --no-cache-dir
if ($LASTEXITCODE -ne 0) {
    Write-ErrorMsg "  ✗ Erreur lors de l'installation de sedpy depuis GitHub"
    Write-ErrorMsg "  Vérifiez que Git est correctement installé et dans le PATH"
    exit 1
}

# Vérifier que sedpy.observate est disponible
Write-Info "  Vérification de sedpy.observate..."
$observateResult = Invoke-PythonCodeOutput -Code "from sedpy import observate; print('OK')"
if (($observateResult.ExitCode -eq 0) -and ($observateResult.Stdout -eq "OK")) {
    Write-Success "  ✓ sedpy installé avec sedpy.observate disponible"
} else {
    Write-ErrorMsg "  ✗ sedpy.observate non disponible - installation échouée"
    Write-ErrorMsg "  Le module 'observate' est requis par Prospector"
    exit 1
}

# ===================================================================
# ÉTAPE 4: Installation des autres dépendances Prospector
# ===================================================================
Write-Step -Step 4 -Total 7 -Message "Installation des autres dépendances Prospector"

# Memes contraintes que pyproject.toml de bd-j/prospector (sauf fsps, compile souvent impossible sous Windows sans MSVC+gfortran)
$prospectorDeps = @(
    "dynesty>=2.0.0",
    "dill>=0.3.0",
    "h5py>=3.0.0",
    "emcee>=3.1.0",
    "corner>=2.2.3",
    "matplotlib>=3.8.4"
)

foreach ($dep in $prospectorDeps) {
    Write-Info "  Installation de $dep..."
    & $pythonPath -m pip install $dep --quiet --upgrade
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "  ⚠ Erreur lors de l'installation de $dep"
    } else {
        Write-Success "  ✓ $dep installé"
    }
}

# ===================================================================
# ÉTAPE 5: Configuration FSPS (optionnel ou création de stubs)
# ===================================================================
Write-Step -Step 5 -Total 7 -Message "Configuration FSPS"

# Définir SPS_HOME
$spsHome = Join-Path $env:USERPROFILE ".local\share\fsps"
$dustDir = Join-Path $spsHome "dust"
$sedDir = Join-Path $spsHome "sed"

# Créer les répertoires
Write-Info "  Création des répertoires FSPS..."
New-Item -ItemType Directory -Force -Path $dustDir | Out-Null
New-Item -ItemType Directory -Force -Path $sedDir | Out-Null
Write-Success "  ✓ Répertoires créés: $spsHome"

# Définir SPS_HOME dans l'environnement utilisateur
[System.Environment]::SetEnvironmentVariable("SPS_HOME", $spsHome, "User")
$env:SPS_HOME = $spsHome
Write-Success "  ✓ SPS_HOME défini: $spsHome"

if ($InstallFSPS) {
    Write-Info "  Installation de FSPS..."
    Write-Warning "  ATTENTION: FSPS nécessite CMake et gfortran sur Windows"
    Write-Warning "  L'installation peut échouer si ces outils ne sont pas configurés"
    
    # Vérifier CMake
    $cmakePath = Get-Command cmake -ErrorAction SilentlyContinue
    if (-not $cmakePath) {
        Write-Warning "  ⚠ CMake non trouvé - FSPS ne pourra pas être compilé"
        Write-Warning "  Installation des fichiers stub FSPS à la place..."
        $InstallFSPS = $false
    } else {
        Write-Success "  ✓ CMake trouvé: $($cmakePath.Source)"
        
        # Vérifier gfortran
        $gfortranPath = Get-Command gfortran -ErrorAction SilentlyContinue
        if (-not $gfortranPath) {
            Write-Warning "  ⚠ gfortran non trouvé - FSPS ne pourra pas être compilé"
            Write-Warning "  Installation des fichiers stub FSPS à la place..."
            Write-Warning "  Voir docs/INSTALLATION_GFORTRAN.md pour installer gfortran"
            $InstallFSPS = $false
        } else {
            Write-Success "  ✓ gfortran trouvé: $($gfortranPath.Source)"
            
            # Essayer d'installer FSPS
            Write-Info "  Installation de FSPS depuis GitHub..."
            $tempDir = New-TemporaryFile | ForEach-Object { 
                Remove-Item $_ -Force -ErrorAction SilentlyContinue
                New-Item -ItemType Directory -Path $_.FullName 
            }
            
            try {
                Push-Location $tempDir.FullName
                
                Write-Info "  Clonage de python-fsps depuis GitHub..."
                git clone https://github.com/dfm/python-fsps.git 2>&1 | Out-Null
                if ($LASTEXITCODE -ne 0) {
                    throw "Échec du clonage de python-fsps"
                }
                
                Set-Location python-fsps
                
                # CRUCIAL: Initialiser les sous-modules Git
                Write-Info "  Initialisation des sous-modules Git (CRUCIAL!)..."
                git submodule update --init --recursive
                if ($LASTEXITCODE -ne 0) {
                    throw "Échec de l'initialisation des sous-modules Git"
                }
                Write-Success "  ✓ Sous-modules Git initialisés"
                
                # Installer FSPS
                Write-Info "  Compilation et installation de FSPS (peut prendre 5-15 minutes)..."
                $env:FC = "gfortran"
                & $pythonPath -m pip install . --no-cache-dir
                if ($LASTEXITCODE -ne 0) {
                    throw "Échec de l'installation de FSPS"
                }
                
                Write-Success "  ✓ FSPS installé et compilé avec succès"
            } catch {
                Write-Warning "  ⚠ Erreur lors de l'installation de FSPS: $_"
                Write-Warning "  Installation des fichiers stub FSPS à la place..."
                $InstallFSPS = $false
            } finally {
                Pop-Location
                Remove-Item -Recurse -Force $tempDir.FullName -ErrorAction SilentlyContinue
            }
        }
    }
}

# FSPS absent : fichiers stub (logique en ASCII pour Windows PowerShell 5.x + encodage fichier)
if (-not $InstallFSPS) {
    Write-Info "  Creation des fichiers stub FSPS..."
    $dustFile = Join-Path $dustDir "Nenkova08_y010_torusg_n10_q2.0.dat"
    Write-Info "  Fichier stub: $dustFile"
    $threeSpace = '   '

    $needsRecreate = $true
    if ((Test-Path -LiteralPath $dustFile) -and -not $ForceReinstall) {
        try {
            $lines = @(Get-Content -LiteralPath $dustFile -ErrorAction Stop)
            if ($lines.Count -ge 129) {
                $dataLine = $lines[4].Trim()
                $dataLine2 = $lines[5].Trim()
                if ($dataLine -and -not $dataLine.StartsWith('#') -and $dataLine2 -and -not $dataLine2.StartsWith('#')) {
                    $cols = $dataLine -split $threeSpace
                    $cols2 = $dataLine2 -split $threeSpace
                    if ($cols.Count -eq 10 -and $cols2.Count -eq 10) {
                        $needsRecreate = $false
                        $nl = $lines.Count
                        Write-Success "  OK: stub deja valide (10 colonnes, $nl lignes)."
                    }
                }
            }
        } catch {
            Write-Warning "  Verification stub: $_"
        }
    }

    if ($needsRecreate) {
        if (Test-Path -LiteralPath $dustFile) {
            try {
                Remove-Item -LiteralPath $dustFile -Force -ErrorAction Stop
                Write-Info "  Ancien stub supprime."
            } catch {
                Write-Warning "  Suppression stub: $_"
            }
        }

        try {
            $stubContent = New-Object System.Text.StringBuilder
            [void]$stubContent.AppendLine("# Nenkova08 AGN torus dust model - Stub file")
            [void]$stubContent.AppendLine("# This is a stub file created automatically")
            [void]$stubContent.AppendLine("# Replace with real FSPS data file for full functionality")
            [void]$stubContent.AppendLine("# wave   fnu_5   fnu_10   fnu_20   fnu_30   fnu_40   fnu_60   fnu_80   fnu_100   fnu_150")

            for ($i = 0; $i -lt 125; $i++) {
                $wave = 1.0 + $i * 0.1
                $fnuValues = @()
                for ($j = 0; $j -lt 9; $j++) {
                    $fnuVal = "{0:F6}" -f (($j + 1) * 0.001 + $i * 0.0001)
                    $fnuValues += $fnuVal
                }
                $joined = $fnuValues -join $threeSpace
                $line = "{0:F6}   {1}" -f $wave, $joined
                [void]$stubContent.AppendLine($line)
            }

            $utf8NoBom = New-Object System.Text.UTF8Encoding $false
            [System.IO.File]::WriteAllText($dustFile, $stubContent.ToString(), $utf8NoBom)
            Write-Success "  OK: fichier stub ecrit (4 lignes en-tete + 125 lignes, separateur 3 espaces)."
        } catch {
            Write-ErrorMsg "  Erreur creation stub: $_"
            exit 1
        }
    }
}

# ===================================================================
# ÉTAPE 6: Installation de Prospector depuis GitHub
# ===================================================================
Write-Step -Step 6 -Total 7 -Message "Installation de Prospector depuis GitHub"
Write-Warning "  IMPORTANT: Prospector doit être installé depuis GitHub"
Write-Warning "  Le package 'prospector' sur PyPI est un autre outil (analyse de code Python)"

Write-Info "  Installation depuis: git+https://github.com/bd-j/prospector.git"
Write-Info "  Cela peut prendre quelques minutes..."

# S'assurer que SPS_HOME est défini avant l'installation
if (-not $env:SPS_HOME) {
    $env:SPS_HOME = $spsHome
}

# Tout installer (inclut fsps) si compilateurs OK ; sinon repli --no-deps (fsps absent, stubs SPS_HOME deja crees)
& $pythonPath -m pip install "git+https://github.com/bd-j/prospector.git" --no-cache-dir
$prospectorOk = ($LASTEXITCODE -eq 0)

if (-not $prospectorOk) {
    Write-Warning "  Echec pip avec dependances (souvent: fsps a compiler — CMake/nmake/Fortran absents sous Windows)."
    Write-Info "  Second essai: astro-prospector depuis Git sans dependances pip (fsps non installe)."
    Write-Warning "  Pour spectres FSPS complets: relancez avec -InstallFSPS apres MSVC + gfortran, ou compilez python-fsps."
    & $pythonPath -m pip install "git+https://github.com/bd-j/prospector.git" --no-cache-dir --no-deps
    if ($LASTEXITCODE -ne 0) {
        Write-ErrorMsg "  ✗ Erreur lors de l'installation de Prospector (meme avec --no-deps)"
        Write-ErrorMsg "  Verifiez Git, le reseau et les droits d'ecriture dans l'environnement Conda."
        exit 1
    }
}
Write-Success "  ✓ Prospector installe depuis GitHub"

# ===================================================================
# ÉTAPE 7: Vérification de l'installation
# ===================================================================
if (-not $SkipVerification) {
    Write-Step -Step 7 -Total 7 -Message "Vérification de l'installation"
    
    # Vérifier sedpy
    Write-Info "  Vérification de sedpy..."
    $sedpyR = Invoke-PythonCodeOutput -Code "from sedpy import observate; print('OK')"
    if (($sedpyR.ExitCode -eq 0) -and ($sedpyR.Stdout -eq "OK")) {
        Write-Success "  ✓ sedpy.observate disponible"
    } else {
        Write-ErrorMsg "  ✗ sedpy.observate non disponible"
    }
    
    # Vérifier FSPS (si installé)
    if ($InstallFSPS) {
        Write-Info "  Vérification de FSPS..."
        $fspsR = Invoke-PythonCodeOutput -Code "import fsps; print('OK')"
        if (($fspsR.ExitCode -eq 0) -and ($fspsR.Stdout -eq "OK")) {
            Write-Success "  ✓ FSPS disponible"
            $fspsVerR = Invoke-PythonCodeOutput -Code "import fsps; print(fsps.__version__)"
            if (($fspsVerR.ExitCode -eq 0) -and $fspsVerR.Stdout) {
                Write-Info "    Version: $($fspsVerR.Stdout)"
            }
        } else {
            Write-Warning "  ⚠ FSPS non disponible (mais les fichiers stub sont créés)"
        }
    } else {
        Write-Info "  FSPS: Utilisation des fichiers stub (FSPS non installé)"
        Write-Info "  Pour installer FSPS complètement, voir docs/INSTALLATION_FSPS.md"
    }
    
    # Vérifier Prospector
    Write-Info "  Vérification de Prospector..."
    $prospectR = Invoke-PythonCodeOutput -Code "import prospect; print('OK')"
    if (($prospectR.ExitCode -eq 0) -and ($prospectR.Stdout -eq "OK")) {
        Write-Success "  ✓ prospect importable"
        
        # Obtenir la version
        $pvR = Invoke-PythonCodeOutput -Code "import prospect; print(getattr(prospect, '__version__', 'version inconnue'))"
        if (($pvR.ExitCode -eq 0) -and $pvR.Stdout) {
            Write-Info "    Version: $($pvR.Stdout)"
        }
        
        # Vérifier SpecModel
        $specR = Invoke-PythonCodeOutput -Code "from prospect.models import SpecModel; print('OK')"
        if (($specR.ExitCode -eq 0) -and ($specR.Stdout -eq "OK")) {
            Write-Success "  ✓ SpecModel disponible"
        } else {
            Write-Warning "  ⚠ SpecModel non disponible (utilisez SedModel dans les anciennes versions)"
        }
        
        # Vérifier FastStepBasis
        $fastR = Invoke-PythonCodeOutput -Code "from prospect.sources import FastStepBasis; print('OK')"
        if (($fastR.ExitCode -eq 0) -and ($fastR.Stdout -eq "OK")) {
            Write-Success "  ✓ FastStepBasis disponible"
        } else {
            Write-Warning "  ⚠ FastStepBasis non disponible (FSPS peut être requis)"
        }
        
        # Vérifier fit_model
        $fitR = Invoke-PythonCodeOutput -Code "from prospect.fitting import fit_model; print('OK')"
        if (($fitR.ExitCode -eq 0) -and ($fitR.Stdout -eq "OK")) {
            Write-Success "  ✓ fit_model disponible"
        } else {
            Write-Warning "  ⚠ fit_model non disponible"
        }
    } else {
        Write-ErrorMsg "  ✗ prospect non importable"
        Write-ErrorMsg "  Vérifiez les messages d'erreur ci-dessus"
        
        $impFail = Invoke-PythonCodeOutput -Code "import prospect"
        $errorDetails = (@($impFail.StderrText, $impFail.Stdout) | Where-Object { $_ }) -join "`n"
        if ([string]::IsNullOrWhiteSpace($errorDetails)) { $errorDetails = "(aucun detail capture)" }
        Write-ErrorMsg "  Détails de l'erreur:"
        Write-ErrorMsg $errorDetails
        
        exit 1
    }
}

# ===================================================================
# Résumé final
# ===================================================================
Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "Installation terminée avec succès!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

Write-Info "Résumé de l'installation:"
Write-Info "  - Environnement Conda: $CondaEnv"
Write-Info "  - sedpy: Installé depuis GitHub (avec observate)"
Write-Info "  - Dépendances: dynesty, dill, h5py, emcee, astropy"
if ($InstallFSPS) {
    Write-Info "  - FSPS: Installé et compilé"
} else {
    Write-Info "  - FSPS: Fichiers stub créés (FSPS non installé)"
    Write-Info "    Pour installer FSPS complètement:"
    Write-Info "    - Voir docs/INSTALLATION_FSPS.md pour Windows"
    Write-Info "    - Ou utilisez WSL pour une installation plus simple"
}
Write-Info "  - Prospector: Installé depuis GitHub"
Write-Host ""

Write-Warning "IMPORTANT:"
Write-Warning "  1. Redémarrez l'application NPOAP pour utiliser Prospector"
Write-Warning "  2. SPS_HOME est défini dans votre profil utilisateur Windows"
Write-Warning "  3. Si vous avez des erreurs, consultez les logs dans le dossier 'logs/'"
Write-Host ""

Write-Info "Pour tester l'installation, exécutez:"
Write-Host "  conda activate $CondaEnv" -ForegroundColor Yellow
# Guillemets simples : evite que PowerShell interprete "from" comme mot-cle
Write-Host '  python -c "import prospect; from prospect.models import SpecModel; print(''Prospector OK!'')"' -ForegroundColor Yellow
Write-Host ""
