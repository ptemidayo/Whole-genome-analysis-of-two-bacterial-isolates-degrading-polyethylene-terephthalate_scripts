#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# fix_setup.sh — resolve the three blockers found in the environment audit
#
#   bash fix_setup.sh            # report what needs doing, change nothing
#   bash fix_setup.sh --apply    # actually make the changes
#
# Blockers addressed:
#   1. WSL2 memory cap (reports how to fix; you must restart WSL from Windows)
#   2. Project living on /mnt/c (slow NTFS bridge) -> migrate to native ext4
#   3. Stale GTDBTK_DATA_PATH pointing at a nonexistent directory
#   4. Missing Phase 1 tools -> create conda env
# ---------------------------------------------------------------------------
set -uo pipefail

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

NATIVE_ROOT="$HOME/psoli_P1"
WIN_USER="oluyo"
WSLCONFIG="/mnt/c/Users/${WIN_USER}/.wslconfig"

say()  { printf "\n\033[1;34m==> %s\033[0m\n" "$1"; }
ok()   { printf "  \033[32m[OK]\033[0m %s\n" "$1"; }
warn() { printf "  \033[33m[!!]\033[0m %s\n" "$1"; }
act()  { printf "  \033[36m[->]\033[0m %s\n" "$1"; }

# ===========================================================================
say "1. WSL2 memory allocation"
# ===========================================================================
WSL_RAM=$(free -g | awk '/^Mem:/{print $2}')
echo "  WSL currently sees: ${WSL_RAM} GB"

HOST_RAM_BYTES=$(powershell.exe -NoProfile -Command \
  "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory" 2>/dev/null | tr -d '\r\n ')

if [ -n "$HOST_RAM_BYTES" ] && [ "$HOST_RAM_BYTES" -gt 0 ] 2>/dev/null; then
  HOST_RAM=$(( HOST_RAM_BYTES / 1073741824 ))
  echo "  Windows host has:   ${HOST_RAM} GB"
  # Leave ~4 GB for Windows itself
  SUGGEST=$(( HOST_RAM - 4 ))
  [ "$SUGGEST" -lt 8 ] && SUGGEST=8
else
  HOST_RAM="unknown"
  SUGGEST=12
  warn "Could not query Windows RAM. Assuming 16 GB host; adjust the number below."
fi

echo
echo "  Write this to  C:\\Users\\${WIN_USER}\\.wslconfig  :"
echo "  ---------------------------------------------"
cat <<EOF
  [wsl2]
  memory=${SUGGEST}GB
  processors=8
  swap=16GB
  localhostForwarding=true
EOF
echo "  ---------------------------------------------"

if [ "$APPLY" -eq 1 ]; then
  if [ -f "$WSLCONFIG" ]; then
    cp "$WSLCONFIG" "${WSLCONFIG}.bak.$(date +%s)"
    warn "Existing .wslconfig backed up"
  fi
  cat > "$WSLCONFIG" <<EOF
[wsl2]
memory=${SUGGEST}GB
processors=8
swap=16GB
localhostForwarding=true
EOF
  ok "Wrote $WSLCONFIG"
  echo
  warn "NOW DO THIS, from a Windows PowerShell or CMD prompt:"
  echo "      wsl --shutdown"
  echo "  Then reopen your WSL terminal and confirm with:  free -g"
  echo "  The change does NOT take effect until WSL is shut down and restarted."
else
  act "Run with --apply to write it, or create the file manually in Windows."
fi

# ===========================================================================
say "2. Filesystem location"
# ===========================================================================
CWD=$(pwd)
case "$CWD" in
  /mnt/c/*|/mnt/d/*)
    warn "You are on the Windows drive bridge: $CWD"
    echo "     I/O here is roughly 10-20x slower than native ext4."
    echo "     Assembly and annotation are I/O-bound; this matters."
    echo
    act "Migrate to: $NATIVE_ROOT"
    if [ "$APPLY" -eq 1 ]; then
      mkdir -p "$NATIVE_ROOT"
      ok "Created $NATIVE_ROOT"
      echo
      echo "  Copy your sequencing data across, e.g.:"
      echo "    mkdir -p $NATIVE_ROOT/00_raw"
      echo "    cp /mnt/c/Users/${WIN_USER}/Downloads/Secuenciación_Dra_Ninfa/**/*.fastq.gz \\"
      echo "       $NATIVE_ROOT/00_raw/"
      echo
      echo "  Keep the /mnt/c copy as your backup. Work from $NATIVE_ROOT."
    fi
    ;;
  *)
    ok "Working on native Linux filesystem: $CWD"
    ;;
esac

# ===========================================================================
say "3. Stale GTDB-Tk configuration"
# ===========================================================================
if [ -n "${GTDBTK_DATA_PATH:-}" ] && [ ! -e "${GTDBTK_DATA_PATH}" ]; then
  warn "GTDBTK_DATA_PATH set to a nonexistent path:"
  echo "     ${GTDBTK_DATA_PATH}"
  echo
  echo "  GTDB-Tk classify_wf needs 55-150 GB RAM depending on mode."
  echo "  It will not run on this machine even with the database present."
  echo "  Phase 2 will use fastANI + TYGS (web) instead, which is sufficient"
  echo "  for a species-level assignment and is standard in Pseudomonas taxonomy."
  echo
  act "Remove the stale export from your shell rc files:"
  grep -n "GTDBTK_DATA_PATH" "$HOME/.bashrc" "$HOME/.profile" 2>/dev/null \
    | sed 's/^/       /' || echo "       (not found in .bashrc/.profile)"
  if [ "$APPLY" -eq 1 ]; then
    for rc in "$HOME/.bashrc" "$HOME/.profile"; do
      if [ -f "$rc" ] && grep -q "GTDBTK_DATA_PATH" "$rc"; then
        cp "$rc" "${rc}.bak.$(date +%s)"
        sed -i '/GTDBTK_DATA_PATH/d' "$rc"
        ok "Removed stale export from $rc (backup kept)"
      fi
    done
  fi
else
  ok "No stale GTDB-Tk path"
fi

# ===========================================================================
say "4. Phase 1 environment"
# ===========================================================================
if conda env list 2>/dev/null | grep -q "psoli_qc_asm"; then
  ok "Environment psoli_qc_asm already exists"
else
  act "Create it with:"
  echo "       mamba env create -f env1_qc_assembly_lowmem.yml"
  echo
  echo "  Note: several tools (samtools, mafft, hmmer, diamond, blast, busco,"
  echo "  gtdbtk) are installed in your base environment. Leave them; the new"
  echo "  environment takes precedence when activated. Do not install more into"
  echo "  base -- it is how environments become unsolvable."
fi

# ===========================================================================
say "Summary"
# ===========================================================================
cat <<EOF

  Order of operations:

    1. Write .wslconfig, then from Windows:  wsl --shutdown
    2. Reopen WSL, confirm with:             free -g      (expect ~${SUGGEST} GB)
    3. Migrate project to:                   $NATIVE_ROOT
    4. Create the environment:               mamba env create -f env1_qc_assembly_lowmem.yml
    5. Inventory your data (see next message) and we proceed to Step 1.

EOF
