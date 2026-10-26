#!/usr/bin/env Rscript
## DIAGNOSTIC - PRJNA628529 is lane-split: 10 runs = 5 libraries x 2 lanes (L6, L7).
## The first-pass analysis treated the 10 runs as 10 independent replicates. This refits
## after summing technical replicates, as edgeR's sumTechReps is meant to be used.
##
## HISTORICAL: this is the diagnostic that justified excluding PRJNA628529. It needs the
## PRE-EXCLUSION sample sheet (90 runs / 5 BioProjects); against the released 80-run sheet
## it will find no PRJNA628529 rows and stop. See reports/design_taxonomy_and_628529.md.
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
ROOT <- "/home/jorge/sweet_recycler_fixes/13_comment/example"
P <- "PRJNA628529"; FDR <- 0.05; LFC <- 1
say <- function(...) cat(sprintf(...), "\n", sep = "")

samples <- fread(file.path(ROOT,"out","samples_final.tsv"), sep="\t", colClasses="character")
meta    <- fread(file.path(ROOT,"..","..","sweet_recycler_RNA.csv"), colClasses="character")
cnt     <- fread(file.path(ROOT,"out","counts_lengthscaled.tsv.gz"), sep="\t")
mat     <- as.matrix(cnt[,-1,with=FALSE]); rownames(mat) <- cnt[["gene_id"]]; rm(cnt); invisible(gc())
univ    <- readLines(file.path(ROOT,"out","de","universe.txt"))

s <- samples[bioproject == P]
m <- meta[match(s$run, meta$Run)]
s$lib      <- sub("_L[0-9]+$", "", m$`Library Name`)   # 201_L6 -> 201
s$lane     <- sub("^.*_(L[0-9]+)$", "\\1", m$`Library Name`)
s$biosample<- m$BioSample
s$arm      <- factor(s$arm, levels=c("control","drought"))
say("=== %s as published: %d runs ===", P, nrow(s))
print(s[, .(run, lib, lane, biosample, arm)])
say("\ndistinct libraries: %d   distinct BioSamples: %d", uniqueN(s$lib), uniqueN(s$biosample))
say("library -> arm:"); print(unique(s[, .(lib, arm)])[order(arm, lib)])
say("\nbiological replicates per arm: control %d, drought %d  <- the >=3-per-arm rule uses RUNS (4 and 6)",
    uniqueN(s[arm=="control", lib]), uniqueN(s[arm=="drought", lib]))

## how similar are the two lanes of the same library?
cpm_all <- cpm(mat[, s$run], log = TRUE, prior.count = 2)
expressed <- rowSums(mat[, s$run]) > 50
say("\n-- Pearson r on log-CPM (%d expressed genes) --", sum(expressed))
same <- c(); diff_same_arm <- c()
for (i in seq_len(nrow(s))) for (j in seq_len(nrow(s))) if (i < j) {
  r <- cor(cpm_all[expressed, i], cpm_all[expressed, j])
  if (s$lib[i] == s$lib[j]) same <- c(same, r)
  else if (s$arm[i] == s$arm[j]) diff_same_arm <- c(diff_same_arm, r)
}
say("   same library, different lane (n=%d): median %.5f  range %.5f-%.5f", length(same), median(same), min(same), max(same))
say("   different library, same arm (n=%d): median %.5f  range %.5f-%.5f",
    length(diff_same_arm), median(diff_same_arm), min(diff_same_arm), max(diff_same_arm))

fit_and_report <- function(label, counts, ss) {
  d <- model.matrix(~ arm, data = ss)
  dge <- DGEList(counts = counts, samples = ss)
  keep <- filterByExpr(dge, design = d)
  dge <- dge[keep, , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge)
  v <- voom(dge, d, plot = FALSE); fit <- eBayes(lmFit(v, d))
  tt <- topTable(fit, coef = "armdrought", number = Inf, sort.by = "none")
  inu <- rownames(tt) %in% univ
  sig <- tt$adj.P.Val < FDR & abs(tt$logFC) > LFC
  say("\n--- %s ---", label)
  say("  n = %d   residual df = %d   filtered genes = %d", nrow(ss), nrow(ss) - ncol(d), sum(keep))
  say("  median residual SD (sigma) = %.4f", median(fit$sigma))
  say("  DE: up %d  down %d   |   in universe: up %d  down %d  (%.1f%% of the 35,692)",
      sum(sig & tt$logFC>0), sum(sig & tt$logFC<0),
      sum(sig & tt$logFC>0 & inu), sum(sig & tt$logFC<0 & inu),
      100*sum(sig & inu)/length(univ))
  invisible(list(tt=tt, sig=sig, inu=inu))
}

a <- fit_and_report("AS PUBLISHED - 10 runs treated as 10 replicates", mat[, s$run], s)

## sum technical replicates (the two lanes of each library)
dge0 <- DGEList(counts = mat[, s$run])
dge_m <- sumTechReps(dge0, ID = s$lib)
sm <- unique(s[, .(lib, arm)])[match(colnames(dge_m), lib)]
stopifnot(identical(sm$lib, colnames(dge_m)))
b <- fit_and_report("LANES SUMMED - 5 libraries, the real biological n", dge_m$counts, sm)

say("\n--- overlap between the two gene lists (universe only) ---")
ga <- rownames(a$tt)[a$sig & a$inu]; gb <- rownames(b$tt)[b$sig & b$inu]
say("  published %d   lane-summed %d   shared %d   Jaccard %.3f",
    length(ga), length(gb), length(intersect(ga,gb)),
    length(intersect(ga,gb))/length(union(ga,gb)))

## what the core sets become if PRJNA628529 is dropped from the >=4-of-5 rule
core <- fread(file.path(ROOT,"out","core_sets.tsv"))
say("\ncore_sets.tsv columns: %s", paste(names(core), collapse=", "))
print(head(core, 3))
