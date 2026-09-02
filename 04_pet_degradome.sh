#!/usr/bin/env bash
# ==============================================================================
# 04_pet_degradome.sh  --  Phase 4: mining the PET degradome
#
# USAGE:
#   bash 04_pet_degradome.sh fetch    # download reference plastic-active enzymes
#   bash 04_pet_degradome.sh verify   # print what was actually retrieved
#
# NOTE ON ACCESSIONS:
#   Accession numbers are easy to get wrong (see: the P. brenneri incident).
#   The fetch step records the protein NAME returned by UniProt for every
#   accession, and `verify` prints them. Read that list before trusting the
#   database -- if a name does not match what the comment says it should be,
#   that entry is wrong and must be removed.
# ==============================================================================

set -uo pipefail

PROJECT="${PROJECT:-$HOME/psoli_project}"
REFDIR="$PROJECT/pet_refs"
DB="$REFDIR/pet_reference_enzymes.faa"
MANIFEST="$REFDIR/manifest.tsv"

say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
ok()   { printf '    \033[0;32m[ok]\033[0m %s\n' "$*"; }
warn() { printf '    \033[0;33m[!]\033[0m  %s\n' "$*"; }

# ------------------------------------------------------------------------------
# Reference enzymes: accession | short_label | what it should be
# Grouped by chemistry, because the groups mean different things biologically.
# ------------------------------------------------------------------------------
REFS=(
  # --- PET hydrolases, mesophilic (the ones that matter most for you: your
  #     strains came from ambient-temperature environments, not compost) ---
  "A0A0K8P6T7|IsPETase|PETase Ideonella sakaiensis 201-F6"
  "A0A0K8P8E7|MHETase|MHETase Ideonella sakaiensis 201-F6"

  # --- PET hydrolases, thermophilic (benchmarks; high activity but from
  #     actinomycetes/compost, so expect low identity to your Pseudomonas) ---
  "G9BY57|LCC|leaf-branch compost cutinase"
  "W0TJ64|Cut190|Saccharomonospora viridis AHK190 cutinase"
  "E9LVH8|TfCut2|Thermobifida fusca KW3 cutinase"

  # --- classical cutinase outgroup ---
  "P00590|FsCut|Fusarium solani pisi cutinase"

  # --- alpha/beta hydrolase reference points ---
  "P22862|EstA|Pseudomonas aeruginosa esterase A autotransporter"
  "Q9HXP6|LipA|Pseudomonas aeruginosa lipase A"
)

# ------------------------------------------------------------------------------
do_fetch() {
  command -v curl >/dev/null || { echo "curl not found"; exit 1; }
  mkdir -p "$REFDIR"
  : > "$DB"
  printf 'accession\tlabel\texpected\tretrieved_name\tlength\tstatus\n' > "$MANIFEST"

  say "Fetching reference enzymes from UniProt"
  for entry in "${REFS[@]}"; do
    local acc label expect
    acc="${entry%%|*}"
    label=$(echo "$entry" | cut -d'|' -f2)
    expect=$(echo "$entry" | cut -d'|' -f3)

    local tmp="$REFDIR/.tmp_${acc}.fa"
    if curl -sf "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "$tmp" && [[ -s "$tmp" ]]; then
      local hdr len
      hdr=$(head -1 "$tmp" | sed 's/^>//')
      len=$(grep -v '^>' "$tmp" | tr -d '\n' | wc -c)
      # rewrite header to a clean short label, keep original in manifest
      { printf '>%s_%s\n' "$label" "$acc"
        grep -v '^>' "$tmp"; } >> "$DB"
      printf '%s\t%s\t%s\t%s\t%s\tOK\n' "$acc" "$label" "$expect" "$hdr" "$len" >> "$MANIFEST"
      ok "$label ($acc) -- ${len} aa"
    else
      printf '%s\t%s\t%s\t-\t-\tFAILED\n' "$acc" "$label" "$expect" >> "$MANIFEST"
      warn "$label ($acc) -- FAILED to retrieve"
    fi
    rm -f "$tmp"
    sleep 0.3
  done

  local n
  n=$(grep -c '^>' "$DB" 2>/dev/null || echo 0)
  say "Retrieved $n / ${#REFS[@]} reference sequences"
  echo "    database : $DB"
  echo "    manifest : $MANIFEST"
  echo
  echo "    Next:  bash 04_pet_degradome.sh verify"
}

# ------------------------------------------------------------------------------
do_verify() {
  [[ -s "$MANIFEST" ]] || { echo "No manifest -- run 'fetch' first."; exit 1; }

  say "CHECK THESE BY EYE -- does each retrieved name match the expectation?"
  echo
  awk -F'\t' 'NR>1 {
    printf "  %-12s %-10s %s\n", $1, $2, $6
    printf "      expected : %s\n", $3
    printf "      retrieved: %s\n", substr($4, 1, 95)
    printf "      length   : %s aa\n\n", $5
  }' "$MANIFEST"

  local nfail
  nfail=$(awk -F'\t' 'NR>1 && $6=="FAILED"' "$MANIFEST" | wc -l)
  if [[ "$nfail" -gt 0 ]]; then
    warn "$nfail accession(s) failed to download -- these need replacing."
  fi

  cat <<'EOF'

    WHAT TO LOOK FOR:
      - Organism in the retrieved name should match the expectation.
      - IsPETase should be ~290 aa, MHETase ~600 aa, LCC ~290 aa,
        cutinases ~200-300 aa. A wildly different length means wrong entry.
      - Anything that does not match: tell me which, and I will replace it.

    Do not proceed to the search step until this list is clean.
EOF
}

case "${1:-}" in
  fetch)  do_fetch  ;;
  verify) do_verify ;;
  *) cat <<EOF

  04_pet_degradome.sh -- Phase 4

    bash 04_pet_degradome.sh fetch    download reference enzymes
    bash 04_pet_degradome.sh verify   check what came back

EOF
  ;;
esac
