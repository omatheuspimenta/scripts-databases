#!/usr/bin/env bash
# Step 3 - assemble the length-scaled gene count matrix from the DISTRIBUTED export.
#
# Reads, for every run in the sample sheet:
#   <MATRIX_DIR>/<RUN>.tar.xz  ->  <RUN>/salmon/matrices.tar.xz
#                              ->  salmon.merged.gene_counts_length_scaled.tsv
# and column-binds the single data column of each (keyed by Experiment accession in the
# header) into out/counts_lengthscaled.tsv.gz, with columns renamed to Run accessions.
#
# The internal pipeline tree (.../rna/results/) is deliberately NOT touched.
#
# Streaming: each inner matrix (~14 MB, 190,974 lines) is extracted into a fresh empty
# directory, its 3rd column is written to a small per-run column file, and the extracted
# matrix is deleted immediately.
set -euo pipefail

ROOT="${ROOT:-/home/jorge/sweet_recycler_fixes/13_comment/example}"
MATRIX_DIR="${MATRIX_DIR:-/mnt/hd02/sb100/sugarcane/rna/matrices}"
SAMPLES="${SAMPLES:-$ROOT/out/samples_final.tsv}"
INNER="salmon/matrices.tar.xz"
TARGET="salmon.merged.gene_counts_length_scaled.tsv"
REF_MD5="${REF_MD5:-932ec246d946c970fa81c502711b4bc9}"

WORK="$ROOT/work/step3"
COLS="$WORK/cols"
EXTRACT="$WORK/extract"          # fresh empty dir per run; untrusted archive contents
CHECKSUMS="$ROOT/work/step3_checksums.tsv"
OUT_MATRIX="$ROOT/out/counts_lengthscaled.tsv.gz"
OUT_MAP="$ROOT/out/run2experiment.tsv"

rm -rf "$WORK"
mkdir -p "$COLS"
printf 'run\texperiment\tcol1_md5\torder_ok\tn_rows\n' > "$CHECKSUMS"
printf 'run\texperiment\n' > "$OUT_MAP"

n=0
tail -n +2 "$SAMPLES" | cut -f1,2 | while IFS=$'\t' read -r RUN EXP; do
    n=$((n+1))
    OUTER="$MATRIX_DIR/$RUN.tar.xz"
    [ -f "$OUTER" ] || { echo "FATAL: missing archive $OUTER" >&2; exit 1; }

    rm -rf "$EXTRACT"; mkdir -p "$EXTRACT"
    tar -xJf "$OUTER" -C "$EXTRACT" "$RUN/$INNER"
    tar -xJf "$EXTRACT/$RUN/$INNER" -C "$EXTRACT" "$TARGET"
    TSV="$EXTRACT/$TARGET"

    # 3. header check: third field must be the Experiment accession from the sample sheet
    HDR3=$(awk -F'\t' 'NR==1{print $3; exit}' "$TSV")
    NF1=$(awk -F'\t' 'NR==1{print NF; exit}' "$TSV")
    if [ "$HDR3" != "$EXP" ]; then
        echo "FATAL: $RUN header field 3 is '$HDR3' but sample sheet says '$EXP'" >&2; exit 1
    fi
    if [ "$NF1" != "3" ]; then
        echo "FATAL: $RUN matrix has $NF1 header fields, expected 3" >&2; exit 1
    fi

    # 4. gene-order checksum (identical command form as the coordinator used)
    MD5=$(tail -n +2 "$TSV" | cut -f1 | md5sum | cut -d' ' -f1)
    NROWS=$(( $(wc -l < "$TSV") - 1 ))
    if [ "$MD5" = "$REF_MD5" ]; then OK=yes; else OK=no; fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$RUN" "$EXP" "$MD5" "$OK" "$NROWS" >> "$CHECKSUMS"
    if [ "$OK" != "yes" ]; then
        echo "FATAL: $RUN column-1 md5 $MD5 != reference $REF_MD5; column bind invalid" >&2; exit 1
    fi

    # gene_id column, once, from the first run
    if [ ! -f "$WORK/gene_ids.col" ]; then
        { printf 'gene_id\n'; tail -n +2 "$TSV" | cut -f1; } > "$WORK/gene_ids.col"
    fi

    # 5. append the value column, renamed to the Run accession
    { printf '%s\n' "$RUN"; tail -n +2 "$TSV" | cut -f3; } > "$COLS/$RUN.col"
    printf '%s\t%s\n' "$RUN" "$EXP" >> "$OUT_MAP"

    # 6. delete the extracted files before the next run
    rm -rf "$EXTRACT"
    echo "[$n] $RUN $EXP md5=$OK rows=$NROWS"
done

# column bind in sample-sheet order
ORDER=()
while IFS=$'\t' read -r RUN _; do ORDER+=("$COLS/$RUN.col"); done < <(tail -n +2 "$SAMPLES" | cut -f1,2)
paste "$WORK/gene_ids.col" "${ORDER[@]}" | gzip -c > "$OUT_MATRIX"

echo "wrote $OUT_MATRIX"
