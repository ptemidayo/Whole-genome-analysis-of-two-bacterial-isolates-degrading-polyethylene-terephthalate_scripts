#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# STEP 1 (revised) — Assess the two existing assemblies
#
# You have assemblies, not reads, so this replaces the original read-QC and
# assembly steps. It answers three questions:
#
#   1. How good are these assemblies?          (seqkit, QUAST)
#   2. Are they complete and uncontaminated?   (CheckM2, BUSCO)
#   3. Are the gene clusters we care about likely to be intact?  (contig sizes)
#
#   conda activate psoli_assess
#   bash 01_assess_assemblies.sh
#
# Run it from the folder containing the two .fasta files.
# ---------------------------------------------------------------------------
set -uo pipefail

P1="Pseudomonas_sp_dr_ninfa.fasta"
P2="pseudochrobactrum_cleaned.fasta"
OUT="01_assessment"
THREADS=8

say()  { printf "\n\033[1;34m==> %s\033[0m\n" "$1"; }
ok()   { printf "  \033[32m[OK]\033[0m %s\n" "$1"; }
warn() { printf "  \033[33m[!!]\033[0m %s\n" "$1"; }

# --- sanity checks --------------------------------------------------------
for f in "$P1" "$P2"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: cannot find $f in $(pwd)"
    echo "Run this script from the folder containing both .fasta files,"
    echo "or edit the P1= and P2= lines at the top."
    exit 1
  fi
done

mkdir -p "$OUT"
V="$OUT/VERSIONS.txt"
{
  echo "# Step 1 assembly assessment — $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  echo "# P1: $P1"
  echo "# P2: $P2"
  echo
} > "$V"

# ===========================================================================
say "1. Basic assembly statistics"
# ===========================================================================
# Pure-awk contig length extraction, so this works even before seqkit exists
contig_lengths() {
  awk '/^>/{if(n)print n"\t"L; n=substr($0,2); L=0; next}{L+=length($0)}END{if(n)print n"\t"L}' "$1"
}

if command -v seqkit >/dev/null 2>&1; then
  echo "seqkit: $(seqkit version 2>&1 | head -1)" >> "$V"
  seqkit stats -a -T "$P1" "$P2" > "$OUT/seqkit_stats.tsv" 2>/dev/null
  seqkit stats -a "$P1" "$P2" | tee "$OUT/seqkit_stats.txt"
else
  warn "seqkit not installed — using built-in fallback"
  echo "seqkit: NOT INSTALLED (awk fallback used)" >> "$V"
fi

echo
echo "  Assembly summary:"
{
  printf "%-40s %8s %12s %10s %10s %8s\n" "file" "contigs" "total_bp" "N50" "longest" "GC%"
  for f in "$P1" "$P2"; do
    contig_lengths "$f" | sort -k2 -nr | awk -v name="$(basename "$f")" -v file="$f" '
      {L[NR]=$2; total+=$2}
      END{
        half=total/2; run=0; n50=0
        for(i=1;i<=NR;i++){run+=L[i]; if(run>=half && n50==0) n50=L[i]}
        gc=0; at=0
        while((getline line < file)>0){
          if(substr(line,1,1)==">") continue
          gc+=gsub(/[GCgc]/,"",line)
          at+=gsub(/[ATat]/,"",line)
        }
        printf "%-40s %8d %12d %10d %10d %7.1f\n", name, NR, total, n50, L[1], 100*gc/(gc+at)
      }'
  done
} | tee "$OUT/assembly_summary.txt" | sed 's/^/    /'

echo
echo "  Contig length distribution (top 10 per assembly):"
for f in "$P1" "$P2"; do
  echo "  --- $(basename "$f")"
  contig_lengths "$f" | sort -k2 -nr | head -10 \
    | awk '{printf "      %-40s %10d bp\n", $1, $2}'
done

# How much of each assembly sits on contigs large enough to hold an operon?
# A tph or ped cluster spans roughly 8-15 kb; contigs below ~20 kb are
# unlikely to contain a full cluster with flanking context.
echo
echo "  Assembly fraction on contigs >= 20 kb (operon-capable):"
for f in "$P1" "$P2"; do
  contig_lengths "$f" | awk -v name="$(basename "$f")" '
    {total+=$2; if($2>=20000){big+=$2; n++}}
    END{printf "      %-40s %5.1f%%  (%d contigs)\n", name, (total?100*big/total:0), n+0}'
done

# ===========================================================================
say "2. QUAST — contiguity metrics"
# ===========================================================================
if command -v quast.py >/dev/null 2>&1; then
  echo "quast: $(quast.py --version 2>&1 | head -1)" >> "$V"
  quast.py "$P1" "$P2" -o "$OUT/quast" -t "$THREADS" \
    --labels "P1_Pseudomonas,P2_Pseudochrobactrum" \
    > "$OUT/quast.log" 2>&1
  ok "QUAST report: $OUT/quast/report.txt"
  echo
  sed -n '1,30p' "$OUT/quast/report.txt" 2>/dev/null | sed 's/^/    /'
else
  warn "quast.py not found — skipping (mamba install -n psoli_assess quast)"
fi

# ===========================================================================
say "3. CheckM2 — completeness and contamination"
# ===========================================================================
if command -v checkm2 >/dev/null 2>&1; then
  echo "checkm2: $(checkm2 --version 2>&1 | head -1)" >> "$V"

  if [ -z "${CHECKM2DB:-}" ]; then
    warn "CHECKM2DB not set. Download the database once (~3 GB):"
    echo "       checkm2 database --download --path \$HOME/checkm2_db"
    echo "       export CHECKM2DB=\$HOME/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
    echo "     Then re-run this script."
  else
    mkdir -p "$OUT/checkm2_input"
    cp "$P1" "$P2" "$OUT/checkm2_input/"
    checkm2 predict --threads "$THREADS" \
      --input "$OUT/checkm2_input" \
      --output-directory "$OUT/checkm2" \
      --force > "$OUT/checkm2.log" 2>&1
    if [ -f "$OUT/checkm2/quality_report.tsv" ]; then
      ok "CheckM2 complete"
      echo
      column -t -s $'\t' "$OUT/checkm2/quality_report.tsv" | cut -c1-110 | sed 's/^/    /'
    else
      warn "CheckM2 did not produce a report — see $OUT/checkm2.log"
    fi
    rm -rf "$OUT/checkm2_input"
  fi
else
  warn "checkm2 not found — skipping"
fi

# ===========================================================================
say "4. BUSCO — single-copy orthologue completeness"
# ===========================================================================
if command -v busco >/dev/null 2>&1; then
  echo "busco: $(busco --version 2>&1 | head -1)" >> "$V"
  for f in "$P1" "$P2"; do
    tag=$(basename "$f" .fasta)
    busco -i "$f" -m genome -l bacteria_odb10 \
      -o "busco_${tag}" --out_path "$OUT" -c "$THREADS" -f \
      > "$OUT/busco_${tag}.log" 2>&1
  done
  echo
  grep -h "C:" "$OUT"/busco_*/short_summary*.txt 2>/dev/null | sed 's/^/    /'
  ok "BUSCO summaries in $OUT/busco_*/"
else
  warn "busco not found in this environment (it is in your base env)"
  echo "     Either: conda activate base && busco ..."
  echo "     Or:     mamba install -n psoli_assess busco"
fi

# ===========================================================================
say "Summary"
# ===========================================================================
cat <<'EOF'

  What to look for:

    Completeness   >95%  good        <90%  the assembly is missing genes,
                                            and absence-of-gene claims become
                                            unsafe -- this matters a lot for us
    Contamination  <5%   good        >10%  mixed assembly, needs binning
    N50            >50 kb good       <20 kb gene clusters will be fragmented
    Total length   P. soli ~5.5 Mb expected
                   Pseudochrobactrum ~4.5-5.0 Mb expected

  The 3.7 Mb length of the Pseudochrobactrum assembly is worth watching. If
  CheckM2 completeness comes back low, that genome may be missing real content,
  and "gene X is absent" would then be unprovable rather than a finding.

  Files to send back:
    01_assessment/seqkit_stats.txt
    01_assessment/quast/report.txt
    01_assessment/checkm2/quality_report.tsv
    01_assessment/busco_*/short_summary*.txt

EOF
