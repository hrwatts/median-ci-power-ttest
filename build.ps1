$ErrorActionPreference = "Stop"

if (Get-Command py -ErrorAction SilentlyContinue) {
    $pythonCommand = "py"
    $pythonArgs = @("-3.10")
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $pythonCommand = "python"
    $pythonArgs = @()
} else {
    throw "Python was not found. Install Python 3.10 or adjust build.ps1."
}

$pythonVersion = (& $pythonCommand @pythonArgs -c "import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}')").Trim()
if ($LASTEXITCODE -ne 0 -or $pythonVersion -ne "3.10") {
    throw "Python 3.10 is required, but '$pythonCommand' resolved to Python $pythonVersion. Install Python 3.10 or adjust build.ps1."
}

if (-not (Get-Command pdflatex -ErrorAction SilentlyContinue)) {
    throw "pdflatex was not found on PATH. Install a LaTeX distribution or add pdflatex to PATH."
}

if (-not (Get-Command bibtex -ErrorAction SilentlyContinue)) {
    throw "bibtex was not found on PATH. Install a LaTeX distribution or add bibtex to PATH."
}

& $pythonCommand @pythonArgs "scripts/median_power_ci.py"
Push-Location paper
try {
    pdflatex -interaction=nonstopmode -halt-on-error paper.tex
    bibtex paper
    pdflatex -interaction=nonstopmode -halt-on-error paper.tex
    pdflatex -interaction=nonstopmode -halt-on-error paper.tex
}
finally {
    Pop-Location
}
