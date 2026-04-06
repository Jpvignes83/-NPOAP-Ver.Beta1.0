# Telechargement Miniconda — sortie console forcee (cmd lance souvent PowerShell sans hote interactif).
param(
    [Parameter(Mandatory = $true)][string]$Url,
    [Parameter(Mandatory = $true)][string]$OutPath
)
$ErrorActionPreference = 'Stop'

# Rattacher la sortie a la console du processus parent (fenetre cmd.exe).
try {
    Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class ParentConsole {
    [DllImport("kernel32.dll", SetLastError = true)]
    public static extern bool AttachConsole(uint dwProcessId);
    public const uint ATTACH_PARENT = 0xFFFFFFFF;
}
'@ -ErrorAction Stop
    [void][ParentConsole]::AttachConsole([ParentConsole]::ATTACH_PARENT)
    [Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
}
catch { }

function Write-Con([string]$Msg) {
    [Console]::Write($Msg)
    [Console]::Out.Flush()
}

try {
    Write-Con ("Telechargement Miniconda... (serveur Anaconda)`r`n")
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
    $dir = Split-Path -Parent $OutPath
    if ($dir -and -not (Test-Path -LiteralPath $dir)) {
        New-Item -ItemType Directory -Path $dir -Force | Out-Null
    }

    $req = [System.Net.HttpWebRequest]::Create($Url)
    $req.UserAgent = 'NPOAP-Windows-Installer'
    $resp = $req.GetResponse()
    $len = [int64]$resp.ContentLength
    $stream = $resp.GetResponseStream()
    $fs = [System.IO.File]::Create($OutPath)
    $buf = New-Object byte[] 65536
    [long]$total = 0
    $lastPct = -1
    $cr = [char]13

    try {
        while (($n = $stream.Read($buf, 0, $buf.Length)) -gt 0) {
            $fs.Write($buf, 0, $n)
            $total += $n
            if ($len -gt 0) {
                $p = [int](100L * $total / $len)
                if ($p -gt $lastPct) {
                    $lastPct = $p
                    $mbT = [math]::Round($total / 1MB, 1)
                    $mbL = [math]::Round($len / 1MB, 1)
                    $fill = [Math]::Min(20, [Math]::Floor($p * 20 / 100))
                    $bar = ('#' * $fill).PadRight(20, '-')
                    Write-Con "${cr}Miniconda  [$bar] ${p}%  (${mbT} / ${mbL} Mo)   "
                }
            }
            else {
                $mbT = [math]::Round($total / 1MB, 1)
                Write-Con "${cr}Miniconda  ${mbT} Mo telecharges (taille inconnue)...   "
            }
        }
    }
    finally {
        $fs.Flush()
        $fs.Close()
        $stream.Close()
        $resp.Close()
    }

    Write-Con "`r`nTermine.`r`n"
    exit 0
}
catch {
    try {
        [Console]::Error.WriteLine('')
        [Console]::Error.WriteLine($_.Exception.Message)
        [Console]::Error.Flush()
    }
    catch { }
    exit 1
}
