#!/usr/bin/env bash
# ==============================================================================
# 03_annotate.sh  --  Phase 3: structural + functional annotation
#
# Annotates all study genomes AND comparator type strains with IDENTICAL
# Prokka settings, so downstream presence/absence comparisons are valid.
#
# USAGE:
#   conda activate genomics_env
#   bash 03_annotate.sh setup     # fetch comparator genomes (needs `datasets`)
#   bash 03_annotate.sh run       # run Prokka on everything
#   bash 03_annotate.sh report    # summary table of gene counts
#
# WHY IDENTICAL SETTINGS MATTER:
#   Your central claim is that P1 has secretion/surface machinery that P2
#   lacks. An absence called by a different annotation pipeline is not an
#   absence, it is an artefact. Same tool, same version, same flags, always.
# ==============================================================================

set -uo pipefail

# ---------------------------------------------------------------- CONFIG -----
# Work on the LINUX filesystem, not /mnt/c -- Prokka does thousands of small
# file operations and the Windows mount makes that 5-10x slower.
PROJECT="${PROJECT:-$HOME/psoli_project}"
WINDIR="/mnt/c/Users/oluyo/Downloads/Secuenciación_Dra_Ninfa/Comparative genomics_Pseudomonas_Dra_ninfa"

GENOMES="$PROJECT/genomes"
ANNOT="$PROJECT/annotation"
LOGS="$PROJECT/logs"

THREADS="${THREADS:-6}"          # leave 2 cores for the OS on an 8-core box
GENUS_P1="Pseudomonas"
GENUS_P2="Pseudochrobactrum"

# Comparator type strains (accession:label:genus)
COMPARATORS=(
  "GCF_017848315.1:P_palmensis_BBB001T:Pseudomonas"
  "GCA_002806685.1:P_qingdaonensis_JJ3:Pseudomonas"
  "GCA_008386555.1:P_saccharolyticum_CCUG33852T:Pseudochrobactrum"
)

# ------------------------------------------------------------- FUNCTIONS -----
say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
ok()   { printf '    \033[0;32m[ok]\033[0m %s\n' "$*"; }
warn() { printf '    \033[0;33m[!]\033[0m  %s\n' "$*"; }
die()  { printf '\n\033[0;31m[FAIL]\033[0m %s\n\n' "$*"; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found. Activate the right conda env first."; }

# ----------------------------------------------------------------- SETUP -----
do_setup() {
  say "Creating project tree under $PROJECT"
  mkdir -p "$GENOMES" "$ANNOT" "$LOGS"
  ok "$PROJECT"

  say "Copying your two isolate assemblies off the Windows mount"
  if [[ -f "$WINDIR/Pseudomonas_sp_dr_ninfa.fasta" ]]; then
    cp "$WINDIR/Pseudomonas_sp_dr_ninfa.fasta" "$GENOMES/P1_Pseudomonas.fasta"
    ok "P1_Pseudomonas.fasta"
  else
    warn "P1 assembly not found at expected path -- copy it to $GENOMES/P1_Pseudomonas.fasta yourself"
  fi

  if [[ -f "$WINDIR/pseudochrobactrum_cleaned.fasta" ]]; then
    cp "$WINDIR/pseudochrobactrum_cleaned.fasta" "$GENOMES/P2_Pseudochrobactrum.fasta"
    ok "P2_Pseudochrobactrum.fasta"
  else
    warn "P2 assembly not found -- copy it to $GENOMES/P2_Pseudochrobactrum.fasta yourself"
  fi

  say "Downloading comparator type strains"
  need datasets
  local accs
  accs=$(printf '%s\n' "${COMPARATORS[@]}" | cut -d: -f1 | paste -sd,)

  ( cd "$GENOMES" && \
    datasets download genome accession "$accs" --include genome \
      --filename comparators.zip >/dev/null 2>&1 && \
    unzip -q -o comparators.zip -d comparators_raw ) || die "comparator download failed"

  for entry in "${COMPARATORS[@]}"; do
    local acc label
    acc="${entry%%:*}"
    label=$(echo "$entry" | cut -d: -f2)
    local src
    src=$(find "$GENOMES/comparators_raw" -path "*${acc}*" -name "*.fna" | head -1)
    if [[ -n "$src" ]]; then
      cp "$src" "$GENOMES/${label}.fasta"
      ok "${label}.fasta"
    else
      warn "could not locate downloaded file for $acc"
    fi
  done

  say "Genomes ready"
  ls -lh "$GENOMES"/*.fasta | awk '{printf "    %-45s %s\n", $9, $5}'
  echo
  echo "    Next:  bash 03_annotate.sh run"
}

# ------------------------------------------------------------------- RUN -----
annotate_one() {
  local fasta="$1" tag="$2" genus="$3"
  local out="$ANNOT/$tag"

  if [[ -s "$out/$tag.faa" ]]; then
    ok "$tag already annotated -- skipping (delete $out to redo)"
    return 0
  fi

  say "Annotating $tag  (genus $genus)"
  # --compliant  : INSDC-valid output, needed if you deposit at NCBI later
  # --rfam       : ncRNA detection; slower but you want it for regulatory context
  # --force      : overwrite partial output from an interrupted run
  prokka \
    --outdir "$out" \
    --prefix "$tag" \
    --locustag "$tag" \
    --genus "$genus" \
    --kingdom Bacteria \
    --gcode 11 \
    --cpus "$THREADS" \
    --mincontiglen 200 \
    --rfam \
    --compliant \
    --force \
    "$fasta" > "$LOGS/prokka_${tag}.log" 2>&1

  if [[ -s "$out/$tag.faa" ]]; then
    local n
    n=$(grep -c '^>' "$out/$tag.faa")
    ok "$tag done -- $n proteins"
  else
    warn "$tag FAILED -- see $LOGS/prokka_${tag}.log"
  fi
}

do_run() {
  need prokka
  mkdir -p "$ANNOT" "$LOGS"

  say "Prokka version"
  prokka --version 2>&1 | sed 's/^/    /'
  warn "Record this version number -- it goes in your Methods."

  [[ -f "$GENOMES/P1_Pseudomonas.fasta" ]] && \
    annotate_one "$GENOMES/P1_Pseudomonas.fasta" "P1_Pseudomonas" "$GENUS_P1"
  [[ -f "$GENOMES/P2_Pseudochrobactrum.fasta" ]] && \
    annotate_one "$GENOMES/P2_Pseudochrobactrum.fasta" "P2_Pseudochrobactrum" "$GENUS_P2"

  for entry in "${COMPARATORS[@]}"; do
    local label genus
    label=$(echo "$entry" | cut -d: -f2)
    genus=$(echo "$entry" | cut -d: -f3)
    [[ -f "$GENOMES/${label}.fasta" ]] && \
      annotate_one "$GENOMES/${label}.fasta" "$label" "$genus"
  done

  say "All annotation attempted"
  echo "    Next:  bash 03_annotate.sh report"
}

# ---------------------------------------------------------------- REPORT -----
do_report() {
  say "Annotation summary"
  local rpt="$ANNOT/annotation_summary.tsv"
  printf 'genome\tcontigs\tbases\tCDS\ttRNA\trRNA\n' > "$rpt"

  for d in "$ANNOT"/*/; do
    local tag
    tag=$(basename "$d")
    [[ -s "$d/$tag.txt" ]] || continue
    local contigs bases cds trna rrna
    contigs=$(awk -F': ' '/^contigs/{print $2}'  "$d/$tag.txt")
    bases=$(awk -F': ' '/^bases/{print $2}'      "$d/$tag.txt")
    cds=$(awk -F': ' '/^CDS/{print $2}'          "$d/$tag.txt")
    trna=$(awk -F': ' '/^tRNA/{print $2}'        "$d/$tag.txt")
    rrna=$(awk -F': ' '/^rRNA/{print $2}'        "$d/$tag.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$tag" "${contigs:-NA}" "${bases:-NA}" "${cds:-NA}" "${trna:-NA}" "${rrna:-NA}" >> "$rpt"
  done

  column -t -s $'\t' "$rpt" | sed 's/^/    /'
  echo
  ok "written to $rpt"

  echo
  say "Sanity checks before you move on"
  cat <<'EOF'
    1. CDS count should be roughly 1 gene per kb:
         ~5,200-5,500 CDS for the 5.76 Mb Pseudomonas genomes
         ~3,500-3,800 CDS for the 3.77 Mb Pseudochrobactrum
       Far outside that range means something went wrong.

    2. P1 and the two Pseudomonas type strains should have SIMILAR
       CDS counts. A big gap would suggest an assembly or annotation
       problem, not real biology.

    3. rRNA counts are usually undercounted in draft assemblies
       (rRNA operons collapse during assembly). Do not read anything
       into a low number here.

    Key output files for the next phase:
       annotation/<tag>/<tag>.faa   proteins   -> HMM / homology searches
       annotation/<tag>/<tag>.gff   coordinates-> operon & synteny work
       annotation/<tag>/<tag>.ffn   gene nt seqs
EOF
}

# ------------------------------------------------------------------ MAIN -----
case "${1:-}" in
  setup)  do_setup  ;;
  run)    do_run    ;;
  report) do_report ;;
  *)
    cat <<EOF

  03_annotate.sh -- Phase 3 annotation

    bash 03_annotate.sh setup     copy isolates + download comparators
    bash 03_annotate.sh run       run Prokka on all genomes (needs genomics_env)
    bash 03_annotate.sh report    summary table + sanity checks

  Run them in that order.

EOF
    ;;
esac
