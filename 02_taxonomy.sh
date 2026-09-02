#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# STEP 2 — Species assignment by average nucleotide identity
#
# The 2024 paper assigned P. soli from 98.59% 16S identity. That is below the
# 98.7% species threshold, so the assignment is provisional. ANI settles it.
#
# Thresholds:
#   ANI >= 95-96%  -> same species
#   ANI 90-95%     -> same genus, different species
#   ANI < 90%      -> more distant
#
# Run in two parts:
#   conda activate psoli_assess     # has ncbi-datasets-cli
#   bash 02_taxonomy.sh download
#
#   conda activate fastani          # your existing env
#   bash 02_taxonomy.sh ani
# ---------------------------------------------------------------------------
set -uo pipefail

P1="Pseudomonas_sp_dr_ninfa.fasta"
P2="pseudochrobactrum_cleaned.fasta"
REFDIR="02_taxonomy/references"
OUT="02_taxonomy"
THREADS=8
MODE="${1:-help}"

say()  { printf "\n\033[1;34m==> %s\033[0m\n" "$1"; }
ok()   { printf "  \033[32m[OK]\033[0m %s\n" "$1"; }
warn() { printf "  \033[33m[!!]\033[0m %s\n" "$1"; }

# Comparators. The Pseudomonas list covers the P. putida group plus the two
# strains that matter functionally: P. aestusnigri (source of PE-H, the only
# characterised Pseudomonas PET hydrolase) and P. umsongensis (terephthalate
# catabolism). The Pseudochrobactrum list covers the described species, with
# Ochrobactrum as an outgroup.
PSEUDOMONAS_SPP=(
  "Pseudomonas soli"
  "Pseudomonas putida"
  "Pseudomonas alloputida"
  "Pseudomonas mosselii"
  "Pseudomonas monteilii"
  "Pseudomonas plecoglossicida"
  "Pseudomonas entomophila"
  "Pseudomonas aestusnigri"
  "Pseudomonas umsongensis"
  "Pseudomonas fluorescens"
)

PSEUDOCHROBACTRUM_SPP=(
  "Pseudochrobactrum saccharolyticum"
  "Pseudochrobactrum asaccharolyticum"
  "Pseudochrobactrum kiredjianiae"
  "Pseudochrobactrum glaciei"
  "Pseudochrobactrum algeriensis"
  "Brucella anthropi"
)

# ===========================================================================
download_refs() {
  say "Downloading reference genomes from NCBI"

  if ! command -v datasets >/dev/null 2>&1; then
    echo "ERROR: 'datasets' not found. Run:  conda activate psoli_assess"
    exit 1
  fi

  mkdir -p "$REFDIR" "$OUT/tmp"

  for sp in "${PSEUDOMONAS_SPP[@]}" "${PSEUDOCHROBACTRUM_SPP[@]}"; do
    tag=$(echo "$sp" | tr ' ' '_')
    if ls "$REFDIR/${tag}"*.fna >/dev/null 2>&1; then
      ok "$sp already downloaded"
      continue
    fi
    echo "  fetching: $sp"
    rm -f "$OUT/tmp/${tag}.zip"
    datasets download genome taxon "$sp" \
      --reference --include genome \
      --filename "$OUT/tmp/${tag}.zip" > /dev/null 2>&1

    if [ -f "$OUT/tmp/${tag}.zip" ]; then
      rm -rf "$OUT/tmp/${tag}_x"
      unzip -qo "$OUT/tmp/${tag}.zip" -d "$OUT/tmp/${tag}_x" 2>/dev/null
      found=$(find "$OUT/tmp/${tag}_x" -name "*.fna" | head -1)
      if [ -n "$found" ]; then
        cp "$found" "$REFDIR/${tag}.fna"
        ok "$sp"
      else
        warn "$sp — no genome in archive"
      fi
      rm -rf "$OUT/tmp/${tag}_x" "$OUT/tmp/${tag}.zip"
    else
      warn "$sp — download failed (no reference genome available?)"
    fi
  done

  rm -rf "$OUT/tmp"
  echo
  ok "References in $REFDIR:"
  ls -1 "$REFDIR" | sed 's/^/      /'
  echo
  warn "IMPORTANT: --reference gives NCBI's representative genome, which is"
  echo "     usually but not always the type strain. Before publishing, verify"
  echo "     each accession is 'assembly from type material' on the NCBI page."
}

# ===========================================================================
run_ani() {
  say "Running fastANI"

  if ! command -v fastANI >/dev/null 2>&1; then
    echo "ERROR: fastANI not found. Run:  conda activate fastani"
    exit 1
  fi

  if [ ! -d "$REFDIR" ] || [ -z "$(ls -A "$REFDIR" 2>/dev/null)" ]; then
    echo "ERROR: no references. Run the download step first."
    exit 1
  fi

  mkdir -p "$OUT"
  ls "$REFDIR"/*.fna > "$OUT/ref_list.txt"
  echo "fastANI: $(fastANI --version 2>&1 | head -1)" > "$OUT/VERSIONS.txt"

  for q in "$P1" "$P2"; do
    tag=$(basename "$q" .fasta)
    echo "  querying: $tag"
    fastANI -q "$q" --rl "$OUT/ref_list.txt" \
      -o "$OUT/ani_${tag}.tsv" -t "$THREADS" \
      > "$OUT/fastani_${tag}.log" 2>&1

    echo
    echo "  --- $tag: top matches by ANI ---"
    if [ -s "$OUT/ani_${tag}.tsv" ]; then
      sort -k3 -nr "$OUT/ani_${tag}.tsv" | head -8 | \
        awk -F'\t' '{
          ref=$2; sub(/.*\//,"",ref); sub(/\.fna$/,"",ref); gsub(/_/," ",ref)
          verdict = ($3>=95) ? "SAME SPECIES" : ($3>=90 ? "same genus" : "distant")
          printf "      %-42s %6.2f%%  (%s/%s frags)  %s\n", ref, $3, $4, $5, verdict
        }'
    else
      warn "no hits above the fastANI reporting cutoff (~80% ANI)"
    fi
    echo
  done

  cat <<'EOF'

  How to read this:

    >= 95-96%   same species as that reference
    90-95%      same genus, but a different species -- and if nothing reaches
                95%, the strain may represent a new species
    < 90%       more distant

  fastANI reports nothing below roughly 80% ANI. An empty result is itself
  informative: it means no close relative in the comparator set.

  Confirm with the dDDH values from TYGS (tygs.dsmz.de). ANI and dDDH
  disagreeing at the boundary is common; both belong in the manuscript.
EOF
}

case "$MODE" in
  download) download_refs ;;
  ani)      run_ani ;;
  *)
    cat <<EOF
Usage:
  conda activate psoli_assess && bash 02_taxonomy.sh download
  conda activate fastani      && bash 02_taxonomy.sh ani
EOF
    ;;
esac
