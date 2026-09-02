#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# check_environment.sh — audit a local machine for the P. soli P1 pipeline
#
#   bash check_environment.sh                 # audit current shell / PATH
#   bash check_environment.sh -e psoli_qc_asm # audit inside a conda env
#   bash check_environment.sh -a              # audit ALL conda envs found
#
# Writes environment_report.md — send that file back and I will write the
# next step against what you actually have, rather than what you might have.
#
# Safe to run: reads only. Installs nothing, changes nothing.
# ---------------------------------------------------------------------------

REPORT="environment_report.md"
TARGET_ENV=""
ALL_ENVS=0

while getopts "e:ao:h" opt; do
  case $opt in
    e) TARGET_ENV="$OPTARG" ;;
    a) ALL_ENVS=1 ;;
    o) REPORT="$OPTARG" ;;
    h) echo "Usage: bash check_environment.sh [-e ENV] [-a] [-o report.md]"; exit 0 ;;
    *) echo "Usage: bash check_environment.sh [-e ENV] [-a] [-o report.md]"; exit 1 ;;
  esac
done

FOUND=0
MISSING=0
MISSING_LIST=""

# ---------------------------------------------------------------------------
# Generic version extractor. Tools disagree wildly on flags and on whether
# they print to stdout or stderr, so try the common forms and take the first
# line containing a digit.
# ---------------------------------------------------------------------------
get_version() {
  local tool="$1" out=""
  for flag in "--version" "-v" "-V" "version" "--help"; do
    out=$("$tool" $flag 2>&1 | grep -m1 -E '[0-9]+\.[0-9]+' | head -c 120)
    [ -n "$out" ] && break
  done
  [ -z "$out" ] && out="(installed, version not parsed)"
  # squeeze whitespace
  echo "$out" | tr -s ' \t' ' ' | sed 's/^ *//;s/ *$//'
}

check_tool() {
  local tool="$1" note="${2:-}"
  if command -v "$tool" >/dev/null 2>&1; then
    local ver path
    ver=$(get_version "$tool")
    path=$(command -v "$tool")
    printf "| %-22s | YES | %-58s |\n" "$tool" "$ver" >> "$REPORT"
    printf "  \033[32m[OK]\033[0m %-22s %s\n" "$tool" "$ver"
    FOUND=$((FOUND+1))
  else
    printf "| %-22s | **NO** | %-58s |\n" "$tool" "$note" >> "$REPORT"
    printf "  \033[31m[--]\033[0m %-22s missing\n" "$tool"
    MISSING=$((MISSING+1))
    MISSING_LIST="$MISSING_LIST $tool"
  fi
}

section() {
  echo "" >> "$REPORT"
  echo "### $1" >> "$REPORT"
  echo "" >> "$REPORT"
  echo "| Tool | Present | Version / note |" >> "$REPORT"
  echo "|---|---|---|" >> "$REPORT"
  echo ""
  echo "--- $1 ---"
}

# ===========================================================================
: > "$REPORT"
{
  echo "# Environment audit — P. soli P1 pipeline"
  echo ""
  echo "- Date: $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "- Host: $(hostname 2>/dev/null || echo unknown)"
  echo "- User: $(whoami 2>/dev/null || echo unknown)"
  echo "- Shell: $SHELL"
  echo "- Conda env active: ${CONDA_DEFAULT_ENV:-none}"
} >> "$REPORT"

echo "Auditing environment. This takes a minute; some tools are slow to report version."
echo

# ---------------------------------------------------------------------------
# System resources — these decide whether certain steps are feasible at all
# ---------------------------------------------------------------------------
echo "" >> "$REPORT"
echo "## 1. System resources" >> "$REPORT"
echo "" >> "$REPORT"

OS=$(uname -s)
CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo "?")
if [ "$OS" = "Darwin" ]; then
  RAM_GB=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))
else
  RAM_GB=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo "?")
fi
DISK=$(df -h . 2>/dev/null | awk 'NR==2{print $4" free of "$2}')

{
  echo "| Resource | Value |"
  echo "|---|---|"
  echo "| OS | $OS $(uname -r 2>/dev/null) |"
  echo "| CPU cores | $CORES |"
  echo "| RAM (GB) | $RAM_GB |"
  echo "| Disk here | $DISK |"
} >> "$REPORT"

echo "  OS: $OS | cores: $CORES | RAM: ${RAM_GB} GB | disk: $DISK"

# Feasibility flags
{
  echo ""
  if [ "$RAM_GB" != "?" ] && [ "$RAM_GB" -lt 16 ] 2>/dev/null; then
    echo "> **RAM below 16 GB.** SPAdes and GTDB-Tk will likely fail. GTDB-Tk alone needs ~70 GB."
  elif [ "$RAM_GB" != "?" ] && [ "$RAM_GB" -lt 80 ] 2>/dev/null; then
    echo "> **RAM below 80 GB.** Assembly and annotation are fine. GTDB-Tk needs ~70 GB and may not run here — use the web-based TYGS plus fastANI instead, which is sufficient for Phase 2."
  fi
} >> "$REPORT"

# ---------------------------------------------------------------------------
# Package managers
# ---------------------------------------------------------------------------
echo "" >> "$REPORT"
echo "## 2. Package managers and languages" >> "$REPORT"

section "Core"
check_tool conda "install Miniforge: https://github.com/conda-forge/miniforge"
check_tool mamba "strongly recommended; conda solves are 10-50x slower"
check_tool python3
check_tool perl "required by Prokka"
check_tool java "required by InterProScan, Bandage"
check_tool R "required by gggenomes, some plotting"
check_tool git
check_tool docker "optional; simplest route for PSORTb"
check_tool singularity "alternative to docker on HPC"
check_tool wget
check_tool curl

# Conda environments present
if command -v conda >/dev/null 2>&1; then
  {
    echo ""
    echo "**Conda environments found:**"
    echo ""
    echo '```'
    conda env list 2>/dev/null
    echo '```'
  } >> "$REPORT"
fi

# ---------------------------------------------------------------------------
# If asked, re-run inside a named env
# ---------------------------------------------------------------------------
if [ -n "$TARGET_ENV" ]; then
  echo ""
  echo "NOTE: to audit inside env '$TARGET_ENV', run:"
  echo "  conda activate $TARGET_ENV && bash $0 -o report_${TARGET_ENV}.md"
  echo ""
fi

# ---------------------------------------------------------------------------
# Phase-by-phase tool checks
# ---------------------------------------------------------------------------
echo "" >> "$REPORT"
echo "## 3. Pipeline tools by phase" >> "$REPORT"

section "Phase 1 — read QC (needed for Step 1, now)"
check_tool fastp
check_tool fastqc
check_tool multiqc
check_tool kraken2
check_tool seqkit
check_tool filtlong "long reads only"
check_tool NanoPlot "long reads only"

section "Phase 1 — assembly and assembly QC (Step 2)"
check_tool spades.py
check_tool unicycler
check_tool flye "long reads only"
check_tool medaka "long reads only"
check_tool trycycler "long reads only"
check_tool polypolish "hybrid only"
check_tool dnaapler
check_tool minimap2
check_tool bwa-mem2
check_tool samtools
check_tool quast.py
check_tool busco
check_tool checkm2
check_tool Bandage "GUI; optional but useful"

section "Phase 2 — taxonomy and phylogenomics"
check_tool fastANI
check_tool average_nucleotide_identity.py "pyani"
check_tool gtdbtk "needs ~70 GB RAM + 110 GB db"
check_tool orthofinder
check_tool mafft
check_tool trimal
check_tool iqtree2
check_tool datasets "ncbi-datasets-cli; for downloading comparators"

section "Phase 3 — annotation"
check_tool bakta
check_tool prokka
check_tool emapper.py "eggNOG-mapper"
check_tool interproscan.sh "standalone install, not conda"
check_tool exec_annotation "kofamscan"
check_tool run_dbcan
check_tool hmmsearch "HMMER"
check_tool diamond
check_tool blastp
check_tool antismash
check_tool barrnap
check_tool signalp6 "DTU licence required"

section "Phase 4-5 — PET degradome and catabolism"
check_tool jackhmmer "HMMER"
check_tool hmmbuild "HMMER"
check_tool foldseek "structure search; key for Phase 4.2"

section "Phase 7-9 — operons, synteny, pangenome, mobile elements"
check_tool panaroo
check_tool ppanggolin
check_tool clinker
check_tool progressiveMauve
check_tool isescan.py
check_tool mob_recon "MOB-suite"
check_tool integron_finder
check_tool genomad
check_tool cmscan "Infernal"
check_tool meme
check_tool transterm "TransTermHP"

section "Phase 10-11 — structure and metabolic model"
check_tool colabfold_batch "or use ColabFold notebooks / AlphaFold server"
check_tool pymol
check_tool vina "AutoDock Vina"
check_tool gmx "GROMACS"
check_tool gapseq
check_tool memote

# ---------------------------------------------------------------------------
# Databases
# ---------------------------------------------------------------------------
echo "" >> "$REPORT"
echo "## 4. Databases" >> "$REPORT"
echo "" >> "$REPORT"
echo "| Database | Env var | Value | Present |" >> "$REPORT"
echo "|---|---|---|---|" >> "$REPORT"
echo ""
echo "--- Databases ---"

check_db() {
  local label="$1" var="$2" val="${!2:-}"
  if [ -n "$val" ] && [ -e "$val" ]; then
    local sz
    sz=$(du -sh "$val" 2>/dev/null | cut -f1)
    echo "| $label | \`$var\` | $val | YES ($sz) |" >> "$REPORT"
    printf "  \033[32m[OK]\033[0m %-18s %s (%s)\n" "$label" "$val" "$sz"
  elif [ -n "$val" ]; then
    echo "| $label | \`$var\` | $val | **set but path missing** |" >> "$REPORT"
    printf "  \033[33m[!!]\033[0m %-18s set but path does not exist: %s\n" "$label" "$val"
  else
    echo "| $label | \`$var\` | (unset) | **NO** |" >> "$REPORT"
    printf "  \033[31m[--]\033[0m %-18s not configured\n" "$label"
  fi
}

check_db "Kraken2"   KRAKEN2_DB_PATH
check_db "CheckM2"   CHECKM2DB
check_db "Bakta"     BAKTA_DB
check_db "eggNOG"    EGGNOG_DATA_DIR
check_db "GTDB-Tk"   GTDBTK_DATA_PATH
check_db "BUSCO"     BUSCO_DOWNLOAD_PATH

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
{
  echo ""
  echo "## 5. Summary"
  echo ""
  echo "- Tools found: **$FOUND**"
  echo "- Tools missing: **$MISSING**"
  echo ""
  if [ "$MISSING" -gt 0 ]; then
    echo "### Missing"
    echo ""
    echo '```'
    for t in $MISSING_LIST; do echo "$t"; done
    echo '```'
  fi
  echo ""
  echo "---"
  echo ""
  echo "Send this file back to continue. Missing tools from later phases are not a"
  echo "problem now — only the Phase 1 read QC block is needed for Step 1."
} >> "$REPORT"

echo
echo "================================================================"
echo "  Found: $FOUND    Missing: $MISSING"
echo "  Report written to: $REPORT"
echo "================================================================"
echo
echo "Send $REPORT back and we continue from there."
