# NPOAP - Re-ecrit requirements.txt en UTF-8 sans BOM pour pip (Windows).
# Gere UTF-8 (avec ou sans BOM), UTF-16 LE/BE avec BOM, UTF-16 LE sans BOM.
param(
    [Parameter(Mandatory = $true)]
    [string] $Path
)
try {
    $Path = [System.IO.Path]::GetFullPath($Path)
} catch {
    Write-Warning "Chemin invalide: $Path"
    exit 1
}
if (-not (Test-Path -LiteralPath $Path)) {
    Write-Warning "Fichier introuvable: $Path"
    exit 0
}
$finfo = Get-Item -LiteralPath $Path -Force
if ($finfo.IsReadOnly) {
    $finfo.IsReadOnly = $false
}
$bytes = [IO.File]::ReadAllBytes($Path)
if ($bytes.Length -eq 0) { exit 0 }

function Write-Utf8NoBom([string] $Text) {
    $utf8 = New-Object System.Text.UTF8Encoding $false
    try {
        [IO.File]::WriteAllText($Path, $Text, $utf8)
    }
    catch [System.UnauthorizedAccessException] {
        $fi = Get-Item -LiteralPath $Path -Force
        $fi.IsReadOnly = $false
        [IO.File]::WriteAllText($Path, $Text, $utf8)
    }
}

[string] $text = $null

# UTF-8 avec BOM
if ($bytes.Length -ge 3 -and $bytes[0] -eq 0xEF -and $bytes[1] -eq 0xBB -and $bytes[2] -eq 0xBF) {
    $utf8 = New-Object System.Text.UTF8Encoding $false
    $text = $utf8.GetString($bytes, 3, $bytes.Length - 3)
    Write-Utf8NoBom $text
    exit 0
}
# UTF-16 LE BOM
if ($bytes.Length -ge 2 -and $bytes[0] -eq 0xFF -and $bytes[1] -eq 0xFE) {
    $text = [Text.Encoding]::Unicode.GetString($bytes, 2, $bytes.Length - 2)
    Write-Utf8NoBom $text
    exit 0
}
# UTF-16 BE BOM
if ($bytes.Length -ge 2 -and $bytes[0] -eq 0xFE -and $bytes[1] -eq 0xFF) {
    $text = [Text.Encoding]::BigEndianUnicode.GetString($bytes, 2, $bytes.Length - 2)
    Write-Utf8NoBom $text
    exit 0
}

# Heuristique UTF-16 LE sans BOM (ex. outil qui garde 23 00 23 00 pour "##")
$n = [Math]::Min(512, $bytes.Length)
if ($n -ge 4 -and ($n % 2) -eq 0) {
    $asciiLEPairs = 0
    for ($i = 0; $i -lt $n - 1; $i += 2) {
        $lo = $bytes[$i]
        $hi = $bytes[$i + 1]
        if ($hi -eq 0 -and (($lo -ge 0x20 -and $lo -le 0x7E) -or $lo -eq 0x09 -or $lo -eq 0x0A -or $lo -eq 0x0D)) {
            $asciiLEPairs++
        }
    }
    $half = [Math]::Floor($n / 2)
    if ($half -gt 0 -and ($asciiLEPairs / $half) -ge 0.72) {
        $text = [Text.Encoding]::Unicode.GetString($bytes)
        Write-Utf8NoBom $text
        exit 0
    }
}

# UTF-8 strict ; si invalide, tenter UTF-16 LE sur tout le fichier
try {
    $encStrict = New-Object System.Text.UTF8Encoding $false, $true
    $text = $encStrict.GetString($bytes)
    Write-Utf8NoBom $text
    exit 0
}
catch {
    $text = [Text.Encoding]::Unicode.GetString($bytes)
    Write-Utf8NoBom $text
    exit 0
}
