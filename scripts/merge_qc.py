#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pandas",
# ]
# ///
"""Merge RNA-seq sample metadata with run-level QC metrics.

Left-joins `rna_sugarcane.csv` (metadata, keyed on `Run`) with
`TableS2_run_level_qc.tsv` (QC metrics, keyed on `run`), keeping only a
fixed set of QC columns, and writes the combined table to CSV.

Usage:
    ./merge_qc.py
    ./merge_qc.py --rna path/to/rna_sugarcane.csv --qc path/to/TableS2_run_level_qc.tsv --out path/to/out.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

DEFAULT_RNA_PATH = Path(
    "rna_sugarcane.csv"
)
DEFAULT_QC_PATH = Path(
    "TableS2_run_level_qc.tsv"
)
DEFAULT_OUT_PATH = Path(
    "rna_sugarcane_qc.csv"
)

RNA_KEY_COL = "Run"
QC_KEY_COL = "run"

QC_COLS = [
    QC_KEY_COL,
    "qc_pass",
    "qc_flags",
    "cutadapt_reads_processed",
    "cutadapt_reads_with_adapters",
    "cutadapt_bp_processed",
    "cutadapt_bp_quality_trimmed",
    "cutadapt_bp_written",
    "cutadapt_pct_reads_with_adapters",
    "cutadapt_pct_bp_removed",
    "cutadapt_version",
    "detected_layout",
    "n_fastq_files",
    "reads_raw",
    "fastqc_raw_percent_gc",
    "fastqc_raw_avg_read_length",
    "fastqc_raw_percent_duplicates",
    "read_length_raw",
    "fastqc_trimmed_poor_quality_reads",
    "reads_trimmed",
    "fastqc_trimmed_percent_gc",
    "fastqc_trimmed_avg_read_length",
    "reads_retained_pct",
    "salmon_library_types",
    "frag_length_mean",
    "frag_length_sd",
    "num_processed",
    "num_mapped",
    "percent_mapped",
    "num_decoy_fragments",
    "num_dovetail_fragments",
    "num_valid_targets",
    "pct_decoy",
    "genes_total",
    "genes_detected_gt0",
    "genes_frac_gt0",
    "genes_detected_ge0p5",
    "genes_frac_ge0p5",
    "genes_detected_ge1",
    "genes_frac_ge1",
    "tx_total",
    "tx_detected_gt0",
    "tx_frac_gt0",
    "tx_detected_ge0p5",
    "tx_frac_ge0p5",
    "tx_detected_ge1",
    "tx_frac_ge1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rna", type=Path, default=DEFAULT_RNA_PATH,
        help="Path to rna_sugarcane.csv (sample metadata, keyed on 'Run')",
    )
    parser.add_argument(
        "--qc", type=Path, default=DEFAULT_QC_PATH,
        help="Path to TableS2_run_level_qc.tsv (QC metrics, keyed on 'run')",
    )
    parser.add_argument(
        "--out", type=Path, default=DEFAULT_OUT_PATH,
        help="Path to write the merged CSV",
    )
    parser.add_argument(
        "--allow-missing-qc-cols", action="store_true",
        help="Warn instead of failing if some expected QC columns are absent",
    )
    return parser.parse_args()


def load_table(path: Path, sep: str) -> pd.DataFrame:
    if not path.exists():
        log.error("File not found: %s", path)
        sys.exit(1)
    return pd.read_csv(path, sep=sep)


def select_qc_columns(qc_df: pd.DataFrame, allow_missing: bool) -> pd.DataFrame:
    missing = [c for c in QC_COLS if c not in qc_df.columns]
    if missing:
        msg = f"QC table is missing expected columns: {missing}"
        if allow_missing:
            log.warning(msg + " (continuing, use --allow-missing-qc-cols was set)")
        else:
            log.error(msg + " (use --allow-missing-qc-cols to continue anyway)")
            sys.exit(1)
    present = [c for c in QC_COLS if c in qc_df.columns]
    return qc_df[present]


def merge_and_validate(df: pd.DataFrame, qc_df: pd.DataFrame) -> pd.DataFrame:
    if RNA_KEY_COL not in df.columns:
        log.error("RNA table has no '%s' column", RNA_KEY_COL)
        sys.exit(1)

    dup_runs = qc_df[qc_df[QC_KEY_COL].duplicated()][QC_KEY_COL].unique()
    if len(dup_runs):
        log.warning(
            "QC table has %d duplicate run id(s), e.g. %s — merge will fan out rows",
            len(dup_runs), list(dup_runs[:5]),
        )

    n_before = len(df)
    merged = df.merge(
        qc_df,
        left_on=RNA_KEY_COL,
        right_on=QC_KEY_COL,
        how="left",
    )
    merged.drop(columns=[QC_KEY_COL], inplace=True)

    if len(merged) != n_before:
        log.warning(
            "Row count changed after merge (%d -> %d); check for duplicate run ids in QC table",
            n_before, len(merged),
        )

    unmatched = merged[merged[[c for c in QC_COLS if c != QC_KEY_COL and c in merged.columns]].isna().all(axis=1)]
    if len(unmatched):
        log.warning(
            "%d/%d rows had no matching QC data for '%s' values, e.g. %s",
            len(unmatched), len(merged), RNA_KEY_COL,
            list(unmatched[RNA_KEY_COL].head(5)),
        )

    return merged


def main() -> None:
    args = parse_args()

    log.info("Loading RNA metadata from %s", args.rna)
    df = load_table(args.rna, sep=",")

    log.info("Loading QC table from %s", args.qc)
    qc_df = load_table(args.qc, sep="\t")

    qc_df = select_qc_columns(qc_df, args.allow_missing_qc_cols)

    merged = merge_and_validate(df, qc_df)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, index=False)
    log.info("Wrote %d rows x %d cols to %s", *merged.shape, args.out)


if __name__ == "__main__":
    main()