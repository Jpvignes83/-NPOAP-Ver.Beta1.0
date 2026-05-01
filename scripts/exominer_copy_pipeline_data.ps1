#Requires -Version 5.1
<#
.SYNOPSIS
  Copie exominer_pipeline/data (.keras + norm_stats) depuis ghcr.io/nasa/exominer vers le clone NPOAP.

.PARAMETER NoPause
  Ne pas attendre la touche Entree a la fin (CI). Env.: NPOAP_EXOMINER_COPY_SILENT=1

.NOTES
  Cible : <NPOAP>/external_apps/Exominer/exominer_pipeline/data
#>

param(
    [switch] $NoPause
)

$ErrorActionPreference = "Stop"

function Get-ContainerCli {
    $names = @("docker", "podman")
    foreach ($n in $names) {
        $cmd = Get-Command $n -ErrorAction SilentlyContinue
        if ($cmd) { return $cmd.Source }
    }
    $extra = @(
        (Join-Path $env:ProgramFiles "RedHat\Podman\podman.exe"),
        (Join-Path $env:ProgramFiles "Docker\Docker\resources\bin\docker.exe"),
        (Join-Path ${env:ProgramFiles(x86)} "Docker\Docker\resources\bin\docker.exe"),
        (Join-Path $env:ProgramFiles "Podman\podman.exe"),
        (Join-Path $env:LOCALAPPDATA "Programs\Podman\podman.exe")
    )
    foreach ($p in $extra) {
        if ($p -and (Test-Path -LiteralPath $p)) { return $p }
    }
    return $null
}

$exitCode = 0

try {
    $ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
    $NpoapRoot = (Resolve-Path (Join-Path $ScriptDir "..")).Path
    $DestPipeline = Join-Path $NpoapRoot "external_apps\Exominer\exominer_pipeline"
    $Image = "ghcr.io/nasa/exominer:latest"
    $SrcInContainer = "/app/exominer_pipeline/data"

    $cli = Get-ContainerCli
    if (-not $cli) {
        Write-Host ""
        Write-Host "ERREUR : docker.exe ou podman introuvable." -ForegroundColor Red
        Write-Host "Installez Docker Desktop (https://www.docker.com/products/docker-desktop/)"
        Write-Host "ou Podman pour Windows, demarrez-le, puis relancez ce script."
        Write-Host ""
        $exitCode = 1
        return
    }

    if (-not (Test-Path -LiteralPath $DestPipeline)) {
        Write-Host "ERREUR : dossier attendu absent : $DestPipeline" -ForegroundColor Red
        Write-Host 'Clonez ExoMiner sous external_apps/Exominer (bouton "Installer / MAJ" dans NPOAP).'
        $exitCode = 1
        return
    }

    $ctrName = "em-npoap-{0}" -f ([Guid]::NewGuid().ToString("N").Substring(0, 10))

    Write-Host ""
    Write-Host "ExoMiner - copie des poids / norm_stats depuis l'image NASA"
    Write-Host "  CLI        : $cli"
    Write-Host "  Image      : $Image"
    Write-Host "  Destination: $DestPipeline\data"
    Write-Host ""

    Write-Host "Telechargement de l'image (peut etre long la premiere fois)..."
    & $cli pull $Image
    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERREUR : docker pull a echoue (code $LASTEXITCODE)." -ForegroundColor Red
        $exitCode = $LASTEXITCODE
        if ($exitCode -eq 0) { $exitCode = 1 }
        return
    }

    Write-Host "Creation du conteneur temporaire $ctrName ..."
    & $cli create --name $ctrName $Image
    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERREUR : docker create a echoue (code $LASTEXITCODE)." -ForegroundColor Red
        $exitCode = $LASTEXITCODE
        if ($exitCode -eq 0) { $exitCode = 1 }
        return
    }

    try {
        Write-Host "Copie vers le clone..."
        & $cli cp "${ctrName}:${SrcInContainer}" $DestPipeline
        if ($LASTEXITCODE -ne 0) {
            throw "La commande cp a echoue (code $LASTEXITCODE)."
        }
        Write-Host ""
        Write-Host "Termine. Verifiez : $DestPipeline\data" -ForegroundColor Green
        Write-Host ""
    } finally {
        Write-Host "Suppression du conteneur temporaire..."
        & $cli rm $ctrName 2>$null
    }
}
catch {
    Write-Host ""
    Write-Host "ERREUR : $($_.Exception.Message)" -ForegroundColor Red
    if ($_.ScriptStackTrace) { Write-Host $_.ScriptStackTrace -ForegroundColor DarkGray }
    $exitCode = 1
}
finally {
    $silent = $NoPause -or ($env:NPOAP_EXOMINER_COPY_SILENT -eq "1")
    if (-not $silent) {
        Write-Host "-------------------------------------------------------------------"
        Read-Host 'Appuyez sur Entree pour fermer cette fenetre'
    }
}

exit $exitCode
