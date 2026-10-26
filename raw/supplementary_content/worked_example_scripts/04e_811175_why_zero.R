#!/usr/bin/env Rscript
## Why does PRJNA811175 return zero DE genes?  Decompose t = logFC / SE.
## One common filter + normalisation; only the design changes, so every number below
## is comparable gene-for-gene.
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
ROOT <- "/home/jorge/sweet_recycler_fixes/13_comment/example"
FDR <- 0.05; LFC <- 1
say <- function(...) cat(sprintf(...), "\n", sep = "")

samples <- fread(file.path(ROOT,"out","samples_final.tsv"), sep="\t", colClasses="character")
meta    <- fread(file.path(ROOT,"..","..","sweet_recycler_RNA.csv"), colClasses="character")
cnt     <- fread(file.path(ROOT,"out","counts_lengthscaled.tsv.gz"), sep="\t")
mat     <- as.matrix(cnt[,-1,with=FALSE]); rownames(mat) <- cnt[["gene_id"]]; rm(cnt); invisible(gc())

s <- samples[bioproject == "PRJNA811175"]
m <- meta[match(s$run, meta$Run)]; stopifnot(identical(m$Run, s$run))
tr <- trimws(m$treatment)
s$time  <- factor(trimws(m$time))
s$level <- factor(ifelse(grepl("^Control$",tr),"control",
                  ifelse(grepl("moderate",tr),"moderate","severe")),
                  levels=c("control","moderate","severe"))
s$arm <- factor(s$arm, levels=c("control","drought")); s$genotype <- factor(s$genotype)

## one filter, one normalisation, shared by every model below
dfilt <- model.matrix(~ level + genotype + time, data = s)
dge <- DGEList(counts = mat[, s$run, drop=FALSE], samples = s)
keep <- filterByExpr(dge, design = dfilt)
dge <- dge[keep, , keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge)
say("common filtered gene set: %d genes, shared by every model below\n", nrow(dge))

fitof <- function(form) { d <- model.matrix(form, data=s)
  stopifnot(qr(d)$rank == ncol(d)); eBayes(lmFit(voom(dge, d, plot=FALSE), d)) }

## stdev.unscaled * sigma = SE of the coefficient
coefstats <- function(fit, cf, label) {
  tt <- topTable(fit, coef=cf, number=Inf, sort.by="none")
  j  <- which(colnames(fit$coefficients) == cf)
  se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  say("%-34s resid df %2d | median |logFC| %.3f | median SE %.3f | median |t| %.3f",
      label, as.integer(fit$df.residual[1]), median(abs(tt$logFC)), median(se), median(abs(tt$t)))
  say("%34s  FDR<0.05 alone: %6d | |logFC|>1 alone: %6d | BOTH: %6d | min BH %.4g",
      "", sum(tt$adj.P.Val < FDR), sum(abs(tt$logFC) > LFC),
      sum(tt$adj.P.Val < FDR & abs(tt$logFC) > LFC), min(tt$adj.P.Val))
  invisible(list(tt=tt, se=se))
}

say("=== 1. Is the problem the FDR threshold or the fold-change threshold? ===")
fP <- fitof(~ arm + genotype + time)
fL <- fitof(~ level + genotype + time)
P  <- coefstats(fP, "armdrought",      "POOLED  drought vs control")
Sv <- coefstats(fL, "levelsevere",     "LEVELS  severe vs control")
Md <- coefstats(fL, "levelmoderate",   "LEVELS  moderate vs control")
G  <- coefstats(fP, "genotypeSP81-3250","POSITIVE CONTROL genotype")

say("\n=== 2. On the genes severe-vs-control calls, what does the pooled test see? ===")
sig <- Sv$tt$adj.P.Val < FDR & abs(Sv$tt$logFC) > LFC
say("severe-significant genes: %d", sum(sig))
say("  severe : median |logFC| %.3f | median SE %.3f | median |t| %.3f | median BH %.3g",
    median(abs(Sv$tt$logFC[sig])), median(Sv$se[sig]), median(abs(Sv$tt$t[sig])), median(Sv$tt$adj.P.Val[sig]))
say("  pooled : median |logFC| %.3f | median SE %.3f | median |t| %.3f | median BH %.3g",
    median(abs(P$tt$logFC[sig])),  median(P$se[sig]),  median(abs(P$tt$t[sig])),  median(P$tt$adj.P.Val[sig]))
say("  ratios pooled/severe -> logFC %.3f | SE %.3f | t %.3f",
    median(abs(P$tt$logFC[sig]))/median(abs(Sv$tt$logFC[sig])),
    median(P$se[sig])/median(Sv$se[sig]),
    median(abs(P$tt$t[sig]))/median(abs(Sv$tt$t[sig])))
say("  of those %d genes, pooled |logFC| exceeds 1 for %d (%.1f%%)",
    sum(sig), sum(abs(P$tt$logFC[sig]) > LFC), 100*mean(abs(P$tt$logFC[sig]) > LFC))

say("\n=== 3. Is the residual variance inflated by the pooling? ===")
for (nm in c("~ arm + genotype", "~ arm + genotype + time", "~ level + genotype + time")) {
  f <- fitof(as.formula(nm))
  say("  %-28s median sigma %.4f | median s2.post %.4f | prior df %.2f",
      nm, median(f$sigma), median(f$s2.post), f$df.prior)
}

say("\n=== 4. Do the two genotypes respond in opposite directions? ===")
fI <- fitof(~ level * genotype + time)
ii <- topTable(fI, coef="levelsevere:genotypeSP81-3250", number=Inf, sort.by="none")
say("  severe x genotype interaction: %d genes at FDR<0.05 (min BH %.3g)",
    sum(ii$adj.P.Val < FDR), min(ii$adj.P.Val))
for (g in levels(s$genotype)) {
  sg <- s[genotype == g & level %in% c("control","severe")]
  sg$level <- factor(as.character(sg$level), levels=c("control","severe"))
  sg$time  <- factor(as.character(sg$time))
  d <- model.matrix(~ level + time, data=sg)
  dd <- dge[, match(sg$run, colnames(dge))]
  dd <- calcNormFactors(dd[filterByExpr(dd, design=d), , keep.lib.sizes=FALSE])
  f <- eBayes(lmFit(voom(dd, d, plot=FALSE), d))
  t2 <- topTable(f, coef="levelsevere", number=Inf, sort.by="none")
  say("  %-12s n=%d (%d severe / %d control) resid df %d -> severe vs control: %d up / %d down",
      g, nrow(sg), sum(sg$level=="severe"), sum(sg$level=="control"),
      as.integer(f$df.residual[1]),
      sum(t2$adj.P.Val<FDR & t2$logFC> LFC), sum(t2$adj.P.Val<FDR & t2$logFC< -LFC))
}

say("\n=== 5. Depth: compare against the shallower BioProject that works ===")
for (p in c("PRJNA811175","PRJNA808768")) {
  sp <- samples[bioproject == p]
  d <- if (uniqueN(sp$genotype) > 1) model.matrix(~ factor(arm, levels=c("control","drought")) + factor(genotype), data=sp)
       else model.matrix(~ factor(arm, levels=c("control","drought")), data=sp)
  g0 <- DGEList(counts = mat[, sp$run, drop=FALSE])
  k  <- filterByExpr(g0, design=d); g0 <- calcNormFactors(g0[k, , keep.lib.sizes=FALSE])
  say("  %-13s n=%2d | median lib %.1f M | genes passing filter %6d | median logCPM %.2f",
      p, nrow(sp), median(g0$samples$lib.size)/1e6, nrow(g0),
      median(rowMeans(cpm(g0, log=TRUE))))
}
