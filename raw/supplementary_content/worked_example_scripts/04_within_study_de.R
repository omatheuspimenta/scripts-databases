#!/usr/bin/env Rscript
## STEP 4 - within-study differential expression (limma-voom), one model per BioProject.
## Usage: Rscript 04_within_study_de.R [ROOT]
## ROOT defaults to the example root; all paths are derived from it.

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT      <- if (length(args) >= 1) args[1] else "/home/jorge/sweet_recycler_fixes/13_comment/example"
COUNTS    <- file.path(ROOT, "out", "counts_lengthscaled.tsv.gz")
SAMPLES   <- file.path(ROOT, "out", "samples_final.tsv")
OUTDIR    <- file.path(ROOT, "out", "de")
WORKDIR   <- file.path(ROOT, "work", "step4")
FDR       <- 0.05
LFC       <- 1
dir.create(OUTDIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n", sep = "")

## ---- inputs -----------------------------------------------------------------
samples <- fread(SAMPLES, sep = "\t", colClasses = "character")
cnt     <- fread(COUNTS, sep = "\t")
say("counts: %d genes x %d columns (incl gene_id)", nrow(cnt), ncol(cnt))
genes   <- cnt[["gene_id"]]
mat     <- as.matrix(cnt[, -1, with = FALSE]); rownames(mat) <- genes
rm(cnt); invisible(gc())
## the matrix may hold a superset of the sheet (e.g. a BioProject excluded after
## it was built); take the sheet's runs, in the sheet's order.
stopifnot(all(samples$run %in% colnames(mat)))
mat <- mat[, samples$run, drop = FALSE]
stopifnot(identical(colnames(mat), samples$run))
say("matrix subset to the sample sheet: %d runs, order identical: TRUE", ncol(mat))
say("all-zero genes across all %d samples: %d", ncol(mat), sum(rowSums(mat) == 0))

## arm coding: control is the REFERENCE level -> coefficient 'armdrought'
samples$arm <- factor(samples$arm, levels = c("control", "drought"))
say("levels(samples$arm) = %s   [reference = %s]",
    paste(levels(samples$arm), collapse = ", "), levels(samples$arm)[1])

projects <- sort(unique(samples$bioproject))
keep_sets <- list(); fits <- list(); meta <- list(); normf <- list()

pdf(file.path(OUTDIR, "diagnostics.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1))

for (p in projects) {
  idx <- which(samples$bioproject == p)
  s   <- samples[idx, ]
  s$arm      <- factor(as.character(s$arm), levels = c("control", "drought"))
  s$genotype <- factor(s$genotype)
  est <- unique(s$genotype_estimable)
  stopifnot(length(est) == 1)
  if (est == "yes") {
    form   <- "~ arm + genotype"
    design <- model.matrix(~ arm + genotype, data = s)
  } else {
    form   <- "~ arm"
    design <- model.matrix(~ arm, data = s)
  }
  rk <- qr(design)$rank
  say("\n=== %s  n=%d (%dd/%dc)  genotypes=%d  design %s  rank %d / ncol %d",
      p, nrow(s), sum(s$arm == "drought"), sum(s$arm == "control"),
      nlevels(s$genotype), form, rk, ncol(design))
  say("    arm levels within project: %s", paste(levels(s$arm), collapse = ", "))
  say("    design columns: %s", paste(colnames(design), collapse = " | "))
  if (rk < ncol(design)) stop(sprintf("RANK DEFICIENT design for %s: rank %d < ncol %d", p, rk, ncol(design)))
  stopifnot("armdrought" %in% colnames(design))

  dge  <- DGEList(counts = mat[, idx, drop = FALSE], samples = s)
  keep <- filterByExpr(dge, design = design)
  say("    filterByExpr kept %d of %d genes", sum(keep), length(keep))
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)                       # TMM
  nf   <- data.table(run = s$run, arm = as.character(s$arm), genotype = as.character(s$genotype),
                     lib_size = dge$samples$lib.size, norm_factor = dge$samples$norm.factors)
  normf[[p]] <- nf
  say("    TMM norm factors: control median %.4f (n=%d), drought median %.4f (n=%d)",
      median(nf$norm_factor[nf$arm == "control"]), sum(nf$arm == "control"),
      median(nf$norm_factor[nf$arm == "drought"]), sum(nf$arm == "drought"))

  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "armdrought", number = Inf, sort.by = "P")
  tt  <- data.table(gene_id = rownames(tt), tt[, c("logFC","AveExpr","t","P.Value","adj.P.Val")])

  ## --- diagnostics: MDS coloured by arm, shaped by genotype; voom trend ------
  cols <- c(control = "#1b7fbf", drought = "#d1662a")[as.character(s$arm)]
  pchs <- c(15, 17, 19, 18, 8)[as.integer(s$genotype)]
  mds  <- plotMDS(v, col = cols, pch = pchs, cex = 1.5,
                  main = sprintf("%s  MDS (n=%d)", p, nrow(s)))
  legend("topright", bty = "n", cex = 0.75,
         legend = c(levels(s$arm), levels(s$genotype)),
         col = c("#1b7fbf", "#d1662a", rep("grey30", nlevels(s$genotype))),
         pch = c(19, 19, c(15, 17, 19, 18, 8)[seq_len(nlevels(s$genotype))]))
  d1 <- mds$x; d2 <- mds$y
  ldep <- log10(dge$samples$lib.size)
  ct <- function(d) {
    a <- suppressWarnings(cor(d, ldep, method = "spearman"))
    b <- suppressWarnings(cor(d, as.integer(s$arm), method = "spearman"))
    c(depth = a, arm = b)
  }
  c1 <- ct(d1); c2 <- ct(d2)
  say("    MDS dim1 Spearman: vs log10(libsize) %+.3f ; vs arm %+.3f", c1["depth"], c1["arm"])
  say("    MDS dim2 Spearman: vs log10(libsize) %+.3f ; vs arm %+.3f", c2["depth"], c2["arm"])
  ## separation of arms along dim1/dim2 (Wilcoxon, descriptive only)
  wt <- function(d) suppressWarnings(wilcox.test(d ~ s$arm)$p.value)
  say("    arm separation Wilcoxon p: dim1 %.4g ; dim2 %.4g", wt(d1), wt(d2))

  voom(dge, design, plot = TRUE)
  title(sub = sprintf("%s  voom mean-variance (%d genes)", p, sum(keep)), cex.sub = 0.9)

  keep_sets[[p]] <- rownames(dge)
  fits[[p]] <- tt
  meta[[p]] <- data.table(bioproject = p,
                          n_drought = sum(s$arm == "drought"), n_control = sum(s$arm == "control"),
                          design_formula = form, design_rank = rk, design_ncol = ncol(design),
                          n_genes_filtered = sum(keep))
}
dev.off()

## ---- common testable universe = intersection of all five filtered sets ------
universe <- Reduce(intersect, keep_sets)
n_proj   <- length(keep_sets)
say("\nfiltered sizes: %s", paste(sprintf("%s=%d", names(keep_sets), sapply(keep_sets, length)), collapse = "  "))
say("union of filtered sets: %d", length(Reduce(union, keep_sets)))
say("COMMON TESTABLE UNIVERSE (%d-way intersection): %d genes", n_proj, length(universe))
writeLines(universe, file.path(OUTDIR, "universe.txt"))

summ <- rbindlist(meta)
for (p in projects) {
  tt <- fits[[p]]
  tt[, in_universe := gene_id %in% universe]
  sig <- tt$adj.P.Val < FDR & abs(tt$logFC) > LFC
  summ[bioproject == p, `:=`(
    n_up            = sum(sig & tt$logFC > 0),
    n_down          = sum(sig & tt$logFC < 0),
    n_up_universe   = sum(sig & tt$logFC > 0 & tt$in_universe),
    n_down_universe = sum(sig & tt$logFC < 0 & tt$in_universe))]
  fwrite(tt, file.path(OUTDIR, sprintf("%s_toptable.tsv", p)), sep = "\t", quote = FALSE)
  say("\n--- %s  top 5 by adj.P.Val ---", p)
  print(head(tt[, .(gene_id, logFC, AveExpr, P.Value, adj.P.Val, in_universe)], 5))
  say("    sig (BH<%.2f & |logFC|>%g): up %d / down %d ; in universe: up %d / down %d",
      FDR, LFC, sum(sig & tt$logFC > 0), sum(sig & tt$logFC < 0),
      sum(sig & tt$logFC > 0 & tt$in_universe), sum(sig & tt$logFC < 0 & tt$in_universe))
}
summ[, design_ncol := NULL]
setcolorder(summ, c("bioproject","n_drought","n_control","design_formula","design_rank",
                    "n_genes_filtered","n_up","n_down","n_up_universe","n_down_universe"))
fwrite(summ, file.path(OUTDIR, "de_summary.tsv"), sep = "\t", quote = FALSE)
say("\n=== de_summary.tsv ===")
print(summ)

nf_all <- rbindlist(normf, idcol = "bioproject")
fwrite(nf_all, file.path(WORKDIR, "norm_factors.tsv"), sep = "\t", quote = FALSE)
say("\n=== PRJNA811175 per-sample TMM factors ===")
print(nf_all[bioproject == "PRJNA811175"][order(arm, -lib_size)])
say("\nsessionInfo versions: edgeR %s limma %s R %s",
    as.character(packageVersion("edgeR")), as.character(packageVersion("limma")),
    paste0(R.version$major, ".", R.version$minor))
