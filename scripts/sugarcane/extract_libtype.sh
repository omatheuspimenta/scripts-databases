#!/usr/bin/env bash
# =============================================================================
# extract_libtype.sh
# -----------------------------------------------------------------------------
# Scans all SRA run directories inside the RNA-seq results folder, reads the
# salmon/<sra_id>/cmd_info.json file for each run, extracts the `libType`
# field, and writes a summary CSV:
#
#   libtype_summary.csv
#   sra_id,libType
#   SRR123456,IU
#   SRR123457,ISR
#   ...
#
# Runs with a missing or unreadable cmd_info.json are skipped with a warning.
# Runs where jq cannot extract `libType` are also skipped with a warning.
#
# Usage:
#   ./extract_libtype.sh [--output /path/to/output.csv]
#
# Options:
#   --output FILE   Path for the CSV output file
#                   (default: <results_dir>/libtype_summary.csv)
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
readonly RESULTS_DIR="/mnt/hd02/sb100/sugarcane/rna/results"
DEFAULT_OUTPUT="${RESULTS_DIR}/libtype_summary.csv"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
OUTPUT_FILE="$DEFAULT_OUTPUT"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output)
            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --output requires a file path argument." >&2
                exit 1
            fi
            OUTPUT_FILE="$2"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            echo "Usage: $0 [--output /path/to/output.csv]" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2; }
err()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
if [[ ! -d "$RESULTS_DIR" ]]; then
    err "Results directory not found: $RESULTS_DIR"
    exit 1
fi

# Ensure jq is available (required for JSON parsing)
if ! command -v jq &>/dev/null; then
    err "'jq' is required but not found. Install it with: sudo apt install jq"
    exit 1
fi

# Ensure the output directory exists (create if needed)
output_dir="$(dirname "$OUTPUT_FILE")"
if [[ ! -d "$output_dir" ]]; then
    log "Creating output directory: $output_dir"
    mkdir -p "$output_dir"
fi

# ---------------------------------------------------------------------------
# Write CSV header
# Using a temp file so we don't overwrite the destination on partial failure
# ---------------------------------------------------------------------------
tmp_output="${OUTPUT_FILE}.tmp"
echo "sra_id,libType" > "$tmp_output"

log "Starting libType extraction."
log "  Source : $RESULTS_DIR"
log "  Output : $OUTPUT_FILE"
echo "---"

# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------
processed=0
skipped=0

# ---------------------------------------------------------------------------
# Main processing loop — iterate over run directories only
# ---------------------------------------------------------------------------
while IFS= read -r -d '' run_dir; do

    sra_id="$(basename "$run_dir")"

    # -----------------------------------------------------------------------
    # Locate cmd_info.json by dynamically finding the sample subdirectory.
    #
    # Inside salmon/ there are exactly two subdirectories:
    #   - deseq2_qc/   (always present, always excluded)
    #   - <sample_id>/ (the actual Salmon output; its name is an experiment/
    #                   sample ID that may differ from the SRA run ID)
    #
    # We iterate subdirectories sorted, skip "deseq2_qc", and take the first
    # remaining entry as the sample subdirectory.
    # -----------------------------------------------------------------------
    salmon_dir="${run_dir}/salmon"
    sample_subdir=""

    if [[ -d "$salmon_dir" ]]; then
        while IFS= read -r -d '' d; do
            dname="$(basename "$d")"
            if [[ "$dname" != "deseq2_qc" ]]; then
                sample_subdir="$dname"
                break
            fi
        done < <(find "$salmon_dir" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)
    fi

    # Resolve the full path to cmd_info.json
    cmd_info=""
    if [[ -n "$sample_subdir" ]]; then
        cmd_info="${salmon_dir}/${sample_subdir}/cmd_info.json"
        log "[$sra_id] Using salmon/${sample_subdir}/cmd_info.json"
    else
        warn "[$sra_id] No sample subdirectory found inside salmon/ (only deseq2_qc or directory missing) — skipping run."
        (( skipped++ )) || true
        continue
    fi

    # Validate that cmd_info.json exists and is non-empty
    if [[ ! -f "$cmd_info" ]]; then
        warn "[$sra_id] cmd_info.json not found: $cmd_info — skipping run."
        (( skipped++ )) || true
        continue
    fi

    if [[ ! -s "$cmd_info" ]]; then
        warn "[$sra_id] cmd_info.json is empty: $cmd_info — skipping run."
        (( skipped++ )) || true
        continue
    fi

    # -----------------------------------------------------------------------
    # Extract the libType field using jq.
    # -r   outputs raw string (no surrounding quotes)
    # -e   causes jq to exit with a non-zero status if the value is
    #      null or false, so we can detect missing keys cleanly.
    # -----------------------------------------------------------------------
    lib_type=""
    if ! lib_type="$(jq -re '.libType' "$cmd_info" 2>/dev/null)"; then
        warn "[$sra_id] Could not extract 'libType' from $cmd_info — skipping run."
        (( skipped++ )) || true
        continue
    fi

    # Sanity-check: libType should be a non-empty string
    if [[ -z "$lib_type" ]]; then
        warn "[$sra_id] 'libType' is empty in $cmd_info — skipping run."
        (( skipped++ )) || true
        continue
    fi

    # -----------------------------------------------------------------------
    # Append the result to the CSV temp file
    # Values are quoted to handle any unexpected whitespace or commas
    # -----------------------------------------------------------------------
    printf '%s,%s\n' "$sra_id" "$lib_type" >> "$tmp_output"

    log "[$sra_id] libType = $lib_type"
    (( processed++ )) || true

done < <(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

# ---------------------------------------------------------------------------
# Move the temp CSV to the final destination only if processing succeeded
# ---------------------------------------------------------------------------
mv "$tmp_output" "$OUTPUT_FILE"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "---"
log "Done."
log "  Processed : $processed"
log "  Skipped   : $skipped"
log "  Output    : $OUTPUT_FILE"

# Print a preview of the CSV if there are any results
if (( processed > 0 )); then
    echo ""
    echo "--- CSV preview (first 10 rows) ---"
    head -n 11 "$OUTPUT_FILE"
fi
