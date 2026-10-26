#!/usr/bin/env python3
"""Step 1B: build gene -> GO map and GO term table for the Comment 13 worked example.

Deliverables (written to --outdir):
  gene2go.tsv   gene_id <TAB> go_id      gene_id in MATRIX convention (locus + ".v2.1")
  go_terms.tsv  go_id <TAB> name <TAB> namespace   (one row per non-obsolete term, plus alt_ids)

Plain stdlib only. Run with:  /usr/bin/python3 -I scripts/01_build_go_assets.py ...
"""

import argparse
import os
import re
import sys

# ---- defaults (stable locations; override on the command line) ----------------
ANNOT_DIR = "/home/jorge/sweet_recycler_fixes/sugarcane_reference/sugarcane/annotation"
DEF_ANNOT_INFO = os.path.join(
    ANNOT_DIR, "SofficinarumxspontaneumR570_771_v2.1.annotation_info.txt")
DEF_GO_TABLE = os.path.join(ANNOT_DIR, "go_table_locus_name.txt")
MATRIX_SUFFIX = ".v2.1"
GO_RE = re.compile(r"GO:\d+")


def parse_annotation_info(path):
    """locus -> set(GO:xxxxxxx) from the authoritative annotation_info file."""
    gene2go = {}
    n_rows = 0
    n_rows_with_go = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_locus = header.index("locusName")
        i_go = header.index("GO")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_locus:
                continue
            n_rows += 1
            locus = f[i_locus].strip()
            if not locus:
                continue
            go_field = f[i_go] if len(f) > i_go else ""
            terms = GO_RE.findall(go_field)
            if terms:
                n_rows_with_go += 1
                gene2go.setdefault(locus, set()).update(terms)
    return gene2go, n_rows, n_rows_with_go


def parse_go_table(path):
    """locus -> set(GO:xxxxxxx) from the pre-made quoted go_table_locus_name.txt."""
    gene2go = {}
    n_lines = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            n_lines += 1
            parts = [p.strip().strip('"') for p in line.split()]
            if len(parts) < 2:
                continue
            locus, raw = parts[0], parts[1]
            if not raw.isdigit():
                continue
            gene2go.setdefault(locus, set()).add("GO:" + raw.zfill(7))
    return gene2go, n_lines


def parse_obo(path):
    """Returns (rows, n_terms, n_obsolete, n_alt). rows = [(go_id, name, namespace)]."""
    rows = []
    n_terms = 0
    n_obsolete = 0
    n_alt = 0
    in_term = False
    cur = {}
    alts = []

    def flush():
        nonlocal n_terms, n_obsolete, n_alt
        if not cur:
            return
        n_terms += 1
        if cur.get("is_obsolete") == "true":
            n_obsolete += 1
            return
        gid = cur.get("id")
        name = cur.get("name", "")
        ns = cur.get("namespace", "")
        if not gid:
            return
        rows.append((gid, name, ns))
        for a in alts:
            rows.append((a, name, ns))
            n_alt += 1

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                if in_term:
                    flush()
                cur = {}
                alts = []
                in_term = (line.strip() == "[Term]")
                continue
            if not in_term or ":" not in line:
                continue
            key, val = line.split(":", 1)
            key = key.strip()
            val = val.split(" ! ")[0].strip()
            if key == "id" and "id" not in cur:
                cur["id"] = val
            elif key == "name" and "name" not in cur:
                cur["name"] = val
            elif key == "namespace" and "namespace" not in cur:
                cur["namespace"] = val
            elif key == "is_obsolete":
                cur["is_obsolete"] = val
            elif key == "alt_id":
                alts.append(val)
    if in_term:
        flush()
    return rows, n_terms, n_obsolete, n_alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot-info", default=DEF_ANNOT_INFO)
    ap.add_argument("--go-table", default=DEF_GO_TABLE)
    ap.add_argument("--obo", required=True, help="path to go-basic.obo")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--matrix-genes", default=None,
                    help="optional file: one matrix gene_id per line (first line may be header)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # --- gene2go from the authoritative annotation_info -----------------------
    g2g, n_rows, n_rows_with_go = parse_annotation_info(a.annot_info)
    pairs = sorted((g + MATRIX_SUFFIX, t) for g, ts in g2g.items() for t in ts)
    out_g2g = os.path.join(a.outdir, "gene2go.tsv")
    with open(out_g2g, "w", encoding="utf-8") as fh:
        fh.write("gene_id\tgo_id\n")
        for g, t in pairs:
            fh.write(g + "\t" + t + "\n")

    uniq_genes = set(g for g, _ in pairs)
    uniq_gos = set(t for _, t in pairs)
    print("[annotation_info] transcript rows: %d  rows with >=1 GO: %d" % (n_rows, n_rows_with_go))
    print("[gene2go.tsv] pairs: %d  unique gene_id: %d  unique go_id: %d"
          % (len(pairs), len(uniq_genes), len(uniq_gos)))

    # --- cross-check vs pre-made go_table ------------------------------------
    gt, gt_lines = parse_go_table(a.go_table)
    gt_pairs = set((g, t) for g, ts in gt.items() for t in ts)
    my_pairs = set((g, t) for g, ts in g2g.items() for t in ts)
    print("[go_table_locus_name.txt] lines: %d  loci: %d  distinct pairs: %d"
          % (gt_lines, len(gt), len(gt_pairs)))
    print("[xcheck] loci in both: %d  only annot_info: %d  only go_table: %d"
          % (len(set(g2g) & set(gt)), len(set(g2g) - set(gt)), len(set(gt) - set(g2g))))
    print("[xcheck] pairs in both: %d  only annot_info: %d  only go_table: %d"
          % (len(my_pairs & gt_pairs), len(my_pairs - gt_pairs), len(gt_pairs - my_pairs)))
    same = sum(1 for g in set(g2g) & set(gt) if g2g[g] == gt[g])
    print("[xcheck] loci with IDENTICAL GO sets: %d" % same)
    ex = [g for g in sorted(set(g2g) & set(gt)) if g2g[g] != gt[g]][:3]
    for g in ex:
        print("[xcheck] differing locus %s: annot_info=%d go_table=%d extra_in_annot=%s"
              % (g, len(g2g[g]), len(gt[g]), sorted(g2g[g] - gt[g])[:5]))

    # --- go_terms from the OBO ------------------------------------------------
    rows, n_terms, n_obsolete, n_alt = parse_obo(a.obo)
    out_terms = os.path.join(a.outdir, "go_terms.tsv")
    seen = set()
    written = 0
    with open(out_terms, "w", encoding="utf-8") as fh:
        fh.write("go_id\tname\tnamespace\n")
        for gid, name, ns in rows:
            if gid in seen:
                continue
            seen.add(gid)
            fh.write("%s\t%s\t%s\n" % (gid, name, ns))
            written += 1
    print("[obo] [Term] stanzas: %d  obsolete skipped: %d  alt_id rows added: %d"
          % (n_terms, n_obsolete, n_alt))
    print("[go_terms.tsv] rows written: %d" % written)

    # --- verification --------------------------------------------------------
    ns_of = {}
    for gid, name, ns in rows:
        ns_of.setdefault(gid, ns)
    unresolved = sorted(t for t in uniq_gos if t not in ns_of)
    print("[verify] go_ids in gene2go not resolvable in go_terms: %d  examples: %s"
          % (len(unresolved), unresolved[:10]))
    counts = {}
    for _, t in pairs:
        counts[ns_of.get(t, "UNRESOLVED")] = counts.get(ns_of.get(t, "UNRESOLVED"), 0) + 1
    for k in sorted(counts):
        print("[verify] pairs in namespace %s: %d" % (k, counts[k]))
    bp_loci = set(g for g, t in pairs if ns_of.get(t) == "biological_process")
    print("[verify] gene_ids with >=1 biological_process term: %d" % len(bp_loci))

    if a.matrix_genes:
        mg = set()
        with open(a.matrix_genes, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                v = line.strip()
                if v and v != "gene_id":
                    mg.add(v)
        print("[verify] matrix gene_ids: %d" % len(mg))
        print("[verify] gene2go gene_ids present in matrix: %d  absent: %d"
              % (len(uniq_genes & mg), len(uniq_genes - mg)))
        print("[verify] BP-annotated gene_ids present in matrix: %d" % len(bp_loci & mg))
    return 0


if __name__ == "__main__":
    sys.exit(main())
