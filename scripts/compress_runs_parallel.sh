#!/usr/bin/env bash
# =============================================================================
# compress_runs.sh (Strict Parallel Edition)
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
export RESULTS_DIR="/mnt/hd02/sb100/sugarcane/rna/results"
export OUTPUT_DIR="/mnt/hd02/sb100/sugarcane/rna/matrices"
export XZ_LEVEL="-9"
export ARCHIVE_EXT=".tar.xz"
export TSV_BUNDLE_NAME="matrices.tar.xz"

# Parallel Settings
THREADS=44 # Leaves 4 cores free for OS tasks

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
export FORCE=false
export BUNDLE_TSV=true

for arg in "$@"; do
    case "$arg" in
        --force)         export FORCE=true        ;;
        --no-bundle-tsv) export BUNDLE_TSV=false  ;;
        *)
            echo "ERROR: Unknown argument: $arg" >&2
            echo "Usage: $0 [--force] [--no-bundle-tsv]" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Setup Counters & Directories
# ---------------------------------------------------------------------------
if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Results directory not found: $RESULTS_DIR" >&2
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Directory to hold empty files representing counts (sub-shells can't increment parent variables)
export COUNTER_DIR="$(mktemp -d)"
mkdir -p "$COUNTER_DIR/processed" "$COUNTER_DIR/skipped" "$COUNTER_DIR/errors"

cleanup_main() {
    [[ -d "$COUNTER_DIR" ]] && rm -rf "$COUNTER_DIR"
}
trap cleanup_main EXIT INT TERM

# ---------------------------------------------------------------------------
# Helpers (Exported for sub-shells)
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2; }
err()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }
export -f log warn err

# ---------------------------------------------------------------------------
# Main Worker Function (Exact match to original loop logic)
# ---------------------------------------------------------------------------
process_run() {
    # Ensure the subshell respects strict bash settings
    set -euo pipefail
    IFS=$'\n\t'

    local run_dir="$1"
    local sra_id="$(basename "$run_dir")"
    local archive="${OUTPUT_DIR}/${sra_id}${ARCHIVE_EXT}"

    log "Processing run: $sra_id"

    if [[ -f "$archive" ]] && [[ "$FORCE" == false ]]; then
        warn "[$sra_id] Archive already exists, skipping (use --force to overwrite): $archive"
        touch "$COUNTER_DIR/skipped/$sra_id"
        return 0
    fi

    # Subshell-local temp directory logic
    local current_tmpdir="$(mktemp -d)"
    trap '[[ -n "$current_tmpdir" && -d "$current_tmpdir" ]] && rm -rf "$current_tmpdir"' EXIT

    local include_args=()

    # 1. MultiQC
    local multiqc_path="${sra_id}/multiqc"
    if [[ -d "${RESULTS_DIR}/${multiqc_path}" ]]; then
        include_args+=( "$multiqc_path" )
    else
        warn "[$sra_id] multiqc/ directory not found — skipping."
    fi

    # 2. Pipeline info
    local params_json_rel=""
    while IFS= read -r -d '' f; do
        params_json_rel="${f#"${RESULTS_DIR}/"}"
        include_args+=( "$params_json_rel" )
    done < <(find "${RESULTS_DIR}/${sra_id}/pipeline_info" -maxdepth 1 -name 'params_*.json' -print0 2>/dev/null || true)

    if [[ -z "$params_json_rel" ]]; then
        warn "[$sra_id] No pipeline_info/params_*.json found — skipping."
    fi

    # 3. Salmon cmd_info.json
    local salmon_dir="${RESULTS_DIR}/${sra_id}/salmon"
    local sample_subdir=""
    if [[ -d "$salmon_dir" ]]; then
        while IFS= read -r -d '' d; do
            local dname="$(basename "$d")"
            if [[ "$dname" != "deseq2_qc" ]]; then
                sample_subdir="$dname"
                break
            fi
        done < <(find "$salmon_dir" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)
    fi

    if [[ -n "$sample_subdir" ]]; then
        local cmd_info_path="${sra_id}/salmon/${sample_subdir}/cmd_info.json"
        if [[ -f "${RESULTS_DIR}/${cmd_info_path}" ]]; then
            include_args+=( "$cmd_info_path" )
            log "[$sra_id] Found cmd_info.json in salmon/${sample_subdir}/"
        else
            warn "[$sra_id] salmon/${sample_subdir}/cmd_info.json not found — skipping."
        fi
    else
        warn "[$sra_id] No sample subdirectory found inside salmon/ (only deseq2_qc or directory missing) — skipping cmd_info.json."
    fi

    # 4. TSV files
    local tsv_files=()
    if [[ -d "$salmon_dir" ]]; then
        while IFS= read -r -d '' f; do
            tsv_files+=( "$f" )
        done < <(find "$salmon_dir" -maxdepth 1 -name '*.tsv' -print0 2>/dev/null || true)
    fi

    local tsv_bundle_ready=false

    if [[ ${#tsv_files[@]} -eq 0 ]]; then
        warn "[$sra_id] No .tsv files found in salmon/ — skipping TSV step."

    elif [[ "$BUNDLE_TSV" == true ]]; then
        local tsv_bundle_tmp="${current_tmpdir}/${TSV_BUNDLE_NAME}"
        log "[$sra_id] Bundling ${#tsv_files[@]} TSV file(s) into ${TSV_BUNDLE_NAME} ..."

        if tar --create --file="$tsv_bundle_tmp" --use-compress-program="xz -T1 ${XZ_LEVEL}" --directory="$salmon_dir" --ignore-failed-read "${tsv_files[@]##*/}" 2>/dev/null; then
            log "[$sra_id] TSV bundle created ($(du -sh "$tsv_bundle_tmp" | cut -f1))"
            tsv_bundle_ready=true
        else
            warn "[$sra_id] Failed to create TSV bundle — falling back to loose .tsv files."
            for f in "${tsv_files[@]}"; do
                include_args+=( "${f#"${RESULTS_DIR}/"}" )
            done
        fi
    else
        for f in "${tsv_files[@]}"; do
            include_args+=( "${f#"${RESULTS_DIR}/"}" )
        done
        log "[$sra_id] Including ${#tsv_files[@]} loose TSV file(s) (--no-bundle-tsv)."
    fi

    # Nothing to archive check
    if [[ ${#include_args[@]} -eq 0 ]] && [[ "$tsv_bundle_ready" == false ]]; then
        warn "[$sra_id] No files found to archive — skipping run."
        touch "$COUNTER_DIR/skipped/$sra_id"
        return 0
    fi

    # Archiving Phase
    local tmp_archive="${archive}.tmp"
    [[ -f "$tmp_archive" ]] && rm -f "$tmp_archive"
    local archive_ok=false

    if [[ "$tsv_bundle_ready" == true ]]; then
        log "[$sra_id] Creating archive with TSV bundle (${#include_args[@]} other path(s) + bundle)."
        
        local stage_dir="${current_tmpdir}/stage/${sra_id}/salmon"
        mkdir -p "$stage_dir"
        cp "${current_tmpdir}/${TSV_BUNDLE_NAME}" "${stage_dir}/${TSV_BUNDLE_NAME}"
        
        local staged_rel="${sra_id}/salmon/${TSV_BUNDLE_NAME}"
        local tar_cmd=( tar --create --file="$tmp_archive" --use-compress-program="xz -T1 ${XZ_LEVEL}" --ignore-failed-read )
        
        [[ ${#include_args[@]} -gt 0 ]] && tar_cmd+=( -C "$RESULTS_DIR" "${include_args[@]}" )
        tar_cmd+=( -C "${current_tmpdir}/stage" "$staged_rel" )

        if "${tar_cmd[@]}" 2>/dev/null; then archive_ok=true; fi
    else
        log "[$sra_id] Creating archive (${#include_args[@]} path(s))."
        if tar --create --file="$tmp_archive" --use-compress-program="xz -T1 ${XZ_LEVEL}" --directory="$RESULTS_DIR" --ignore-failed-read "${include_args[@]}" 2>/dev/null; then
            archive_ok=true
        fi
    fi

    # Result logic
    if [[ "$archive_ok" == true && -f "$tmp_archive" ]]; then
        mv "$tmp_archive" "$archive"
        log "[$sra_id] Archive created successfully: $archive ($(du -sh "$archive" | cut -f1))"
        touch "$COUNTER_DIR/processed/$sra_id"
    else
        err "[$sra_id] Compression failed. Removing incomplete archive."
        [[ -f "$tmp_archive" ]] && rm -f "$tmp_archive"
        touch "$COUNTER_DIR/errors/$sra_id"
    fi
}
export -f process_run

# ---------------------------------------------------------------------------
# Parallel Execution Loop
# ---------------------------------------------------------------------------
log "Starting compression of RNA-seq run directories. (PARALLEL MODE: $THREADS threads)"
log "  Source     : $RESULTS_DIR"
log "  Output     : $OUTPUT_DIR"
log "  Force      : $FORCE"
log "  Bundle TSV : $BUNDLE_TSV"
echo "---"

find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | \
    xargs -0 -n 1 -P "$THREADS" bash -c 'process_run "$@"' _

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
processed=$(ls -1q "$COUNTER_DIR/processed" | wc -l)
skipped=$(ls -1q "$COUNTER_DIR/skipped" | wc -l)
errors=$(ls -1q "$COUNTER_DIR/errors" | wc -l)

echo "---"
log "Done."
log "  Processed : $processed"
log "  Skipped   : $skipped"
log "  Errors    : $errors"

if (( errors > 0 )); then
    exit 1
fi
