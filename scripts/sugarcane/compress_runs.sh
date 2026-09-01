#!/usr/bin/env bash
# =============================================================================
# compress_runs.sh
# -----------------------------------------------------------------------------
# Iterates over SRA run directories inside the RNA-seq results folder and
# creates selective .tar.xz archives containing only the relevant output files:
#   - multiqc/                          (entire directory)
#   - pipeline_info/params_*.json
#   - salmon/<sample_id>/cmd_info.json  (sample_id auto-detected, != deseq2_qc)
#   - salmon/*.tsv                      (bundled into tsv_files.tar.xz by default)
#
# TSV bundling (default ON):
#   All .tsv files found directly inside salmon/ are first compressed together
#   into a single tsv_files.tar.xz in a temporary directory, then that bundle
#   is included in the final archive instead of the loose .tsv files.
#   This significantly reduces the final archive size.
#   Original .tsv files on disk are never modified.
#   Pass --no-bundle-tsv to include the loose .tsv files directly (old behaviour).
#
# Archives are written to the output (matrices) directory with maximum xz
# compression. Existing archives are skipped unless --force is passed.
#
# Usage:
#   ./compress_runs.sh [--force] [--no-bundle-tsv]
#
# Options:
#   --force          Overwrite existing .tar.xz archives (default: skip)
#   --no-bundle-tsv  Include loose .tsv files instead of a pre-bundled archive
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
readonly RESULTS_DIR="/mnt/hd02/sb100/sugarcane/rna/results"
readonly OUTPUT_DIR="/mnt/hd02/sb100/sugarcane/rna/matrices"
readonly XZ_LEVEL="-9"                       # Maximum xz compression
readonly ARCHIVE_EXT=".tar.xz"
readonly TSV_BUNDLE_NAME="matrices.tar.xz"  # Name of the TSV bundle inside the archive

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
FORCE=false
BUNDLE_TSV=true   # ON by default; disabled by --no-bundle-tsv

for arg in "$@"; do
    case "$arg" in
        --force)         FORCE=true       ;;
        --no-bundle-tsv) BUNDLE_TSV=false ;;
        *)
            echo "ERROR: Unknown argument: $arg" >&2
            echo "Usage: $0 [--force] [--no-bundle-tsv]" >&2
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
# Cleanup trap — ensures any per-run temp directory is always removed,
# even if the script exits early due to set -e or a signal.
# ---------------------------------------------------------------------------
CURRENT_TMPDIR=""
cleanup() {
    if [[ -n "$CURRENT_TMPDIR" && -d "$CURRENT_TMPDIR" ]]; then
        rm -rf "$CURRENT_TMPDIR"
    fi
}
trap cleanup EXIT INT TERM

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
if [[ ! -d "$RESULTS_DIR" ]]; then
    err "Results directory not found: $RESULTS_DIR"
    exit 1
fi

# Create output directory if it does not exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    log "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Verify required tools are available
for tool in tar xz find mktemp; do
    if ! command -v "$tool" &>/dev/null; then
        err "Required tool not found: $tool"
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
processed=0
skipped=0
errors=0

log "Starting compression of RNA-seq run directories."
log "  Source     : $RESULTS_DIR"
log "  Output     : $OUTPUT_DIR"
log "  Force      : $FORCE"
log "  Bundle TSV : $BUNDLE_TSV"
echo "---"

# Iterate only over immediate subdirectories (SRA run IDs).
# .tar.gz files and other non-directory items are intentionally ignored.
while IFS= read -r -d '' run_dir; do

    sra_id="$(basename "$run_dir")"
    archive="${OUTPUT_DIR}/${sra_id}${ARCHIVE_EXT}"

    log "Processing run: $sra_id"

    # -----------------------------------------------------------------------
    # Skip if archive already exists and --force was not requested
    # -----------------------------------------------------------------------
    if [[ -f "$archive" ]] && [[ "$FORCE" == false ]]; then
        warn "Archive already exists, skipping (use --force to overwrite): $archive"
        (( skipped++ )) || true
        continue
    fi

    # -----------------------------------------------------------------------
    # Create a per-run temp directory for intermediate files (TSV bundle).
    # The cleanup trap removes it automatically on exit or error.
    # -----------------------------------------------------------------------
    CURRENT_TMPDIR="$(mktemp -d)"

    # -----------------------------------------------------------------------
    # Build the list of files/directories to include.
    # Paths are collected relative to RESULTS_DIR so the archive preserves
    # the top-level run directory (e.g. SRR123456/multiqc/...).
    # -----------------------------------------------------------------------
    include_args=()   # Paths relative to RESULTS_DIR, passed to tar

    # 1. Entire multiqc/ directory
    multiqc_path="${sra_id}/multiqc"
    if [[ -d "${RESULTS_DIR}/${multiqc_path}" ]]; then
        include_args+=( "$multiqc_path" )
    else
        warn "[$sra_id] multiqc/ directory not found — skipping."
    fi

    # 2. pipeline_info/params_*.json  (only one file expected per run)
    params_json_rel=""
    while IFS= read -r -d '' f; do
        params_json_rel="${f#"${RESULTS_DIR}/"}"
        include_args+=( "$params_json_rel" )
    done < <(find "${RESULTS_DIR}/${sra_id}/pipeline_info" \
                  -maxdepth 1 -name 'params_*.json' \
                  -print0 2>/dev/null || true)

    if [[ -z "$params_json_rel" ]]; then
        warn "[$sra_id] No pipeline_info/params_*.json found — skipping."
    fi

    # 3. salmon/<sample_id>/cmd_info.json
    # -----------------------------------------------------------------------
    # Inside salmon/ there are exactly two subdirectories:
    #   - deseq2_qc/   (always present, always excluded)
    #   - <sample_id>/ (the actual Salmon output; its name is an experiment/
    #                   sample ID that may differ from the SRA run ID)
    #
    # We dynamically find the first subdirectory that is NOT "deseq2_qc" and
    # look for cmd_info.json inside it.
    # -----------------------------------------------------------------------
    salmon_dir="${RESULTS_DIR}/${sra_id}/salmon"
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

    if [[ -n "$sample_subdir" ]]; then
        cmd_info_path="${sra_id}/salmon/${sample_subdir}/cmd_info.json"
        if [[ -f "${RESULTS_DIR}/${cmd_info_path}" ]]; then
            include_args+=( "$cmd_info_path" )
            log "[$sra_id] Found cmd_info.json in salmon/${sample_subdir}/"
        else
            warn "[$sra_id] salmon/${sample_subdir}/cmd_info.json not found — skipping."
        fi
    else
        warn "[$sra_id] No sample subdirectory found inside salmon/ (only deseq2_qc or directory missing) — skipping cmd_info.json."
    fi

    # 4. TSV files from salmon/
    # -----------------------------------------------------------------------
    # Collect all .tsv files found directly inside salmon/.
    # Depending on BUNDLE_TSV:
    #
    #   true  (default) — compress them together into a temporary
    #                     tsv_files.tar.xz, then add that single bundle
    #                     to the final archive as:
    #                       <sra_id>/salmon/tsv_files.tar.xz
    #                     The originals on disk are never touched.
    #
    #   false           — add the loose .tsv files directly to the final
    #                     archive (original behaviour, --no-bundle-tsv flag).
    # -----------------------------------------------------------------------

    tsv_files=()
    if [[ -d "$salmon_dir" ]]; then
        while IFS= read -r -d '' f; do
            tsv_files+=( "$f" )
        done < <(find "$salmon_dir" -maxdepth 1 -name '*.tsv' -print0 2>/dev/null || true)
    fi

    tsv_bundle_ready=false   # Will be set true only if bundle succeeds

    if [[ ${#tsv_files[@]} -eq 0 ]]; then
        warn "[$sra_id] No .tsv files found in salmon/ — skipping TSV step."

    elif [[ "$BUNDLE_TSV" == true ]]; then
        # -----------------------------------------------------------------
        # Bundle mode: create tsv_files.tar.xz inside the temp directory,
        # storing only the basenames of the .tsv files (no leading path).
        # -----------------------------------------------------------------
        tsv_bundle_tmp="${CURRENT_TMPDIR}/${TSV_BUNDLE_NAME}"

        log "[$sra_id] Bundling ${#tsv_files[@]} TSV file(s) into ${TSV_BUNDLE_NAME} ..."

        # Change into salmon/ so only basenames are stored inside the bundle
        if tar \
            --create \
            --file="$tsv_bundle_tmp" \
            --use-compress-program="xz ${XZ_LEVEL}" \
            --directory="$salmon_dir" \
            --ignore-failed-read \
            "${tsv_files[@]##*/}"; then   # Strip directory → basename only

            log "[$sra_id] TSV bundle created ($(du -sh "$tsv_bundle_tmp" | cut -f1))"
            tsv_bundle_ready=true
        else
            warn "[$sra_id] Failed to create TSV bundle — falling back to loose .tsv files."
            # Fall back: include loose .tsv files so no data is lost
            for f in "${tsv_files[@]}"; do
                include_args+=( "${f#"${RESULTS_DIR}/"}" )
            done
        fi

    else
        # -----------------------------------------------------------------
        # No-bundle mode (--no-bundle-tsv): include loose .tsv files directly
        # -----------------------------------------------------------------
        for f in "${tsv_files[@]}"; do
            include_args+=( "${f#"${RESULTS_DIR}/"}" )
        done
        log "[$sra_id] Including ${#tsv_files[@]} loose TSV file(s) (--no-bundle-tsv)."
    fi

    # -----------------------------------------------------------------------
    # Nothing to archive — warn but do not treat as a hard error
    # -----------------------------------------------------------------------
    if [[ ${#include_args[@]} -eq 0 ]] && [[ "$tsv_bundle_ready" == false ]]; then
        warn "[$sra_id] No files found to archive — skipping run."
        (( skipped++ )) || true
        rm -rf "$CURRENT_TMPDIR"; CURRENT_TMPDIR=""
        continue
    fi

    # -----------------------------------------------------------------------
    # Create the final .tar.xz archive.
    #
    # When tsv_bundle_ready=true we need to merge two directory roots into one
    # compressed archive:
    #   - RESULTS_DIR  (multiqc, json files, cmd_info.json)
    #   - CURRENT_TMPDIR  (tsv_files.tar.xz, placed at <sra_id>/salmon/ inside
    #                      the archive via --transform)
    #
    # tar cannot append to a compressed file, so we pipe two uncompressed tar
    # streams into a single xz process:
    #
    #   { tar --create [pass A] ; tar --create [pass B] ; } | xz > archive
    #
    # When tsv_bundle_ready=false a single tar invocation is used as before.
    # -----------------------------------------------------------------------
    tmp_archive="${archive}.tmp"

    # Remove stale temp file from a previous interrupted run
    [[ -f "$tmp_archive" ]] && rm -f "$tmp_archive"

    archive_ok=false

    if [[ "$tsv_bundle_ready" == true ]]; then
        # -------------------------------------------------------------------
        # Single-pass approach using a staging tree.
        #
        # The two-pass pipe trick ({ tar; tar; } | xz) is unreliable: tar
        # stops reading at the first end-of-archive marker from Pass A, so
        # Pass B content (the TSV bundle) is silently dropped on extraction.
        #
        # Fix: copy the bundle into a staging directory that mirrors its
        # intended path inside the archive, then run ONE tar with two -C
        # anchors (GNU tar supports multiple -C flags):
        #
        #   tar -C "$RESULTS_DIR"           [run-level files]
        #       -C "$CURRENT_TMPDIR/stage"  <sra_id>/salmon/tsv_files.tar.xz
        #
        # Staging layout:
        #   $CURRENT_TMPDIR/stage/<sra_id>/salmon/tsv_files.tar.xz
        # -------------------------------------------------------------------
        log "[$sra_id] Creating archive with TSV bundle (${#include_args[@]} other path(s) + bundle)."

        tsv_bundle_tmp="${CURRENT_TMPDIR}/${TSV_BUNDLE_NAME}"

        # Build the staging subdirectory that mirrors the target archive path
        stage_dir="${CURRENT_TMPDIR}/stage/${sra_id}/salmon"
        mkdir -p "$stage_dir"
        cp "$tsv_bundle_tmp" "${stage_dir}/${TSV_BUNDLE_NAME}"

        # Relative path used as the tar argument (anchored to stage/)
        staged_rel="${sra_id}/salmon/${TSV_BUNDLE_NAME}"

        # Build tar arguments: first anchor RESULTS_DIR content, then stage
        tar_cmd=( tar --create --file="$tmp_archive"
                      --use-compress-program="xz ${XZ_LEVEL}"
                      --ignore-failed-read )

        if [[ ${#include_args[@]} -gt 0 ]]; then
            tar_cmd+=( -C "$RESULTS_DIR" "${include_args[@]}" )
        fi

        tar_cmd+=( -C "${CURRENT_TMPDIR}/stage" "$staged_rel" )

        if "${tar_cmd[@]}"; then
            archive_ok=true
        fi

    else
        # -------------------------------------------------------------------
        # Single-pass approach (no bundle or fallback)
        # -------------------------------------------------------------------
        log "[$sra_id] Creating archive (${#include_args[@]} path(s))."

        if tar \
            --create \
            --file="$tmp_archive" \
            --use-compress-program="xz ${XZ_LEVEL}" \
            --directory="$RESULTS_DIR" \
            --ignore-failed-read \
            "${include_args[@]}"; then
            archive_ok=true
        fi
    fi

    # -----------------------------------------------------------------------
    # Finalise or clean up
    # -----------------------------------------------------------------------
    if [[ "$archive_ok" == true && -f "$tmp_archive" ]]; then
        mv "$tmp_archive" "$archive"
        log "[$sra_id] Archive created successfully: $archive ($(du -sh "$archive" | cut -f1))"
        (( processed++ )) || true
    else
        err "[$sra_id] Compression failed. Removing incomplete archive."
        [[ -f "$tmp_archive" ]] && rm -f "$tmp_archive"
        (( errors++ )) || true
    fi

    # Clean up per-run temp directory
    rm -rf "$CURRENT_TMPDIR"
    CURRENT_TMPDIR=""

done < <(find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "---"
log "Done."
log "  Processed : $processed"
log "  Skipped   : $skipped"
log "  Errors    : $errors"

# Exit with a non-zero status if any run failed
if (( errors > 0 )); then
    exit 1
fi
