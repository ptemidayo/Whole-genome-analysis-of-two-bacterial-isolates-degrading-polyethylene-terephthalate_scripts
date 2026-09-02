#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# inventory_data.sh — catalogue the sequencing data actually present
#
#   bash inventory_data.sh /path/to/Secuenciación_Dra_Ninfa
#
# Reads only. Reports file types, sizes, read counts, inferred platform,
# pairing, and estimated coverage. Writes data_inventory.md.
# ---------------------------------------------------------------------------
set -uo pipefail

DIR="${1:-.}"
OUT="data_inventory.md"
GENOME_MB=6.0     # assumed Pseudomonas genome size, Mb

[ ! -d "$DIR" ] && { echo "Not a directory: $DIR"; exit 1; }

: > "$OUT"
{
  echo "# Sequencing data inventory"
  echo ""
  echo "- Scanned: \`$DIR\`"
  echo "- Date: $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "- Assumed genome size for coverage: ${GENOME_MB} Mb"
} >> "$OUT"

echo "Scanning $DIR ..."
echo

# ---------------------------------------------------------------------------
# Directory tree, two levels
# ---------------------------------------------------------------------------
{
  echo ""
  echo "## Directory structure"
  echo ""
  echo '```'
  find "$DIR" -maxdepth 2 -type d 2>/dev/null | head -60
  echo '```'
} >> "$OUT"

# ---------------------------------------------------------------------------
# File type census
# ---------------------------------------------------------------------------
{
  echo ""
  echo "## File types"
  echo ""
  echo "| Extension | Count | Total size |"
  echo "|---|---|---|"
} >> "$OUT"

echo "--- File types ---"
for ext in fastq.gz fq.gz fastq fq fasta.gz fa.gz fasta fa fna gbk gbff gff bam sra pod5 fast5; do
  n=$(find "$DIR" -type f -name "*.${ext}" 2>/dev/null | wc -l)
  if [ "$n" -gt 0 ]; then
    sz=$(find "$DIR" -type f -name "*.${ext}" -print0 2>/dev/null | du -ch --files0-from=- 2>/dev/null | tail -1 | cut -f1)
    echo "| .$ext | $n | ${sz:-?} |" >> "$OUT"
    printf "  %-12s %3d files  %s\n" ".$ext" "$n" "${sz:-?}"
  fi
done

# ---------------------------------------------------------------------------
# FASTQ detail
# ---------------------------------------------------------------------------
{
  echo ""
  echo "## FASTQ files"
  echo ""
  echo "| File | Size | Reads | Mean len | Total bases | Est. cov |"
  echo "|---|---|---|---|---|---|"
} >> "$OUT"

echo
echo "--- FASTQ detail (counting reads; this takes a moment) ---"

HAVE_SEQKIT=0
command -v seqkit >/dev/null 2>&1 && HAVE_SEQKIT=1

find "$DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) 2>/dev/null \
| sort | while read -r f; do
  base=$(basename "$f")
  size=$(du -h "$f" 2>/dev/null | cut -f1)

  if [ "$HAVE_SEQKIT" -eq 1 ]; then
    stats=$(seqkit stats -T "$f" 2>/dev/null | tail -1)
    reads=$(echo "$stats" | cut -f4)
    avglen=$(echo "$stats" | cut -f7)
    bases=$(echo "$stats" | cut -f5)
  else
    # fallback: count lines / 4
    case "$f" in
      *.gz) reads=$(( $(zcat "$f" 2>/dev/null | wc -l) / 4 )) ;;
      *)    reads=$(( $(wc -l < "$f") / 4 )) ;;
    esac
    avglen="?"; bases="?"
  fi

  if [ "$bases" != "?" ] && [ -n "$bases" ]; then
    cov=$(awk -v b="$bases" -v g="$GENOME_MB" 'BEGIN{printf "%.0fx", b/(g*1000000)}')
  else
    cov="?"
  fi

  echo "| $base | $size | $reads | $avglen | $bases | $cov |" >> "$OUT"
  printf "  %-42s %6s  reads=%-12s cov=%s\n" "$base" "$size" "$reads" "$cov"
done

# ---------------------------------------------------------------------------
# Platform and pairing inference
# ---------------------------------------------------------------------------
{
  echo ""
  echo "## Platform inference"
  echo ""
  echo "First read header from each FASTQ (identifies the sequencer):"
  echo ""
  echo '```'
} >> "$OUT"

echo
echo "--- Read headers ---"
find "$DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) 2>/dev/null \
| sort | head -20 | while read -r f; do
  case "$f" in
    *.gz) hdr=$(zcat "$f" 2>/dev/null | head -1) ;;
    *)    hdr=$(head -1 "$f") ;;
  esac
  echo "$(basename "$f")" >> "$OUT"
  echo "  $hdr" >> "$OUT"
  printf "  %-42s %s\n" "$(basename "$f")" "$(echo "$hdr" | head -c 70)"
done
echo '```' >> "$OUT"

{
  echo ""
  echo "How to read the headers above:"
  echo ""
  echo "- \`@A00123:45:HXXXX:1:1101:...  1:N:0:INDEX\` -> Illumina, paired-end"
  echo "- \`@<uuid> runid=... sampleid=...\`            -> Oxford Nanopore"
  echo "- \`@m64012_...\`                               -> PacBio"
} >> "$OUT"

# ---------------------------------------------------------------------------
# Sample grouping — how many isolates are here?
# ---------------------------------------------------------------------------
{
  echo ""
  echo "## Inferred samples"
  echo ""
  echo '```'
} >> "$OUT"

echo
echo "--- Inferred sample names ---"
find "$DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null \
| xargs -n1 basename 2>/dev/null \
| sed -E 's/_(R?[12])(_001)?\.(fastq|fq)\.gz$//; s/\.(fastq|fq)\.gz$//' \
| sort -u | tee -a "$OUT" | sed 's/^/  /'
echo '```' >> "$OUT"

# ---------------------------------------------------------------------------
# Existing assemblies, if any
# ---------------------------------------------------------------------------
NFA=$(find "$DIR" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fasta.gz" -o -name "*.fna.gz" \) 2>/dev/null | wc -l)
if [ "$NFA" -gt 0 ]; then
  {
    echo ""
    echo "## Existing assemblies / FASTA"
    echo ""
    echo "| File | Size | Contigs |"
    echo "|---|---|---|"
  } >> "$OUT"
  echo
  echo "--- FASTA files found ---"
  find "$DIR" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fasta.gz" -o -name "*.fna.gz" \) 2>/dev/null \
  | sort | while read -r f; do
    size=$(du -h "$f" | cut -f1)
    case "$f" in
      *.gz) n=$(zcat "$f" 2>/dev/null | grep -c "^>") ;;
      *)    n=$(grep -c "^>" "$f") ;;
    esac
    echo "| $(basename "$f") | $size | $n |" >> "$OUT"
    printf "  %-42s %6s  %s contigs\n" "$(basename "$f")" "$size" "$n"
  done
fi

echo
echo "================================================================"
echo "  Inventory written to: $OUT"
echo "================================================================"
echo
echo "Send $OUT back and we start Step 1."
