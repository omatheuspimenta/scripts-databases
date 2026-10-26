#!/usr/bin/env Rscript
## DIAGNOSTIC - does pooling heterogeneous drought levels into one arm dilute the signal?
## PRJNA811175 (moderate+severe x 3 timepoints) and PRJNA776107 (4/8/16/32 h) both have
## within-arm structure. Compare the pooled arm contrast against level-specific contrasts.
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/jorge/sweet_recycler_fixes/13_comment/example"
FDR <- 0.05; LFC <- 1
say <- function(...) cat(sprintf(...), "\n", sep = "")

samples <- fread(file.path(ROOT,"out","samples_final.tsv"), sep="\t", colClasses="character")
meta    <- fread(file.path(ROOT,"..","..","sweet_recycler_RNA.csv"), colClasses="character")
cnt     <- fread(file.path(ROOT,"out","counts_lengthscaled.tsv.gz"), sep="\t")
mat     <- as.matrix(cnt[,-1,with=FALSE]); rownames(mat) <- cnt[["gene_id"]]
rm(cnt); invisible(gc())
univ <- readLines(file.path(ROOT,"out","de","universe.txt"))

analyse <- function(P, make_factors, pooled_form, level_form, level_coefs) {
  s <- samples[bioproject == P]
  m <- meta[match(s$run, meta$Run)]
  stopifnot(identical(m$Run, s$run))
  s <- make_factors(s, m)
  s$arm <- factor(s$arm, levels=c("control","drought")); s$genotype <- factor(s$genotype)
  say("\n%s", strrep("=", 72)); say("%s   n=%d (%dd/%dc)", P, nrow(s), sum(s$arm=="drought"), sum(s$arm=="control"))
  print(table(level = s$level, arm = s$arm))

  base <- model.matrix(pooled_form, data = s)
  dge  <- DGEList(counts = mat[, s$run, drop=FALSE], samples = s)
  keep <- filterByExpr(dge, design = base)
  dge  <- dge[keep, , keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge)
  inu  <- rownames(dge) %in% univ
  say("filtered %d genes (%d in universe)", nrow(dge), sum(inu))

  fitP <- eBayes(lmFit(voom(dge, base, plot=FALSE), base))
  ttP  <- topTable(fitP, coef="armdrought", number=Inf, sort.by="none")

  dL   <- model.matrix(level_form, data = s)
  stopifnot(qr(dL)$rank == ncol(dL))
  fitL <- eBayes(lmFit(voom(dge, dL, plot=FALSE), dL))
  say("\n  pooled  %s : up %d down %d (universe up %d down %d)",
      paste(deparse(pooled_form),collapse=""),
      sum(ttP$adj.P.Val<FDR & ttP$logFC> LFC), sum(ttP$adj.P.Val<FDR & ttP$logFC< -LFC),
      sum(ttP$adj.P.Val<FDR & ttP$logFC> LFC & inu), sum(ttP$adj.P.Val<FDR & ttP$logFC< -LFC & inu))
  say("  levels  %s", paste(deparse(level_form),collapse=""))
  for (cf in level_coefs) {
    tt <- topTable(fitL, coef=cf, number=Inf, sort.by="none")
    sg <- tt$adj.P.Val<FDR & abs(tt$logFC)>LFC
    say("    %-22s up %5d down %5d (universe up %5d down %5d)   median |logFC| of its sig genes %.2f",
        cf, sum(sg & tt$logFC>0), sum(sg & tt$logFC<0),
        sum(sg & tt$logFC>0 & inu), sum(sg & tt$logFC<0 & inu),
        if (any(sg)) median(abs(tt$logFC[sg])) else NA_real_)
    ## dilution: how does the pooled contrast behave on THIS level's significant genes?
    if (any(sg)) {
      r <- ttP$logFC[sg] / tt$logFC[sg]
      say("      -> on those genes the POOLED logFC is a median %.2fx the level logFC (IQR %.2f-%.2f);",
          median(r), quantile(r,.25), quantile(r,.75))
      say("         pooled calls %d of them significant; pooled median BH on them %.3g",
          sum(ttP$adj.P.Val[sg]<FDR & abs(ttP$logFC[sg])>LFC), median(ttP$adj.P.Val[sg]))
    }
  }
  say("  Spearman cor(pooled t, %s t) = %+.3f", level_coefs[length(level_coefs)],
      cor(ttP$t, topTable(fitL, coef=level_coefs[length(level_coefs)], number=Inf, sort.by="none")$t,
          method="spearman"))

  ## clean, self-consistent variance partition (ordinary least squares on logCPM)
  E <- voom(dge, base, plot=FALSE)$E
  r2 <- function(form) { d <- model.matrix(form, data=s)
    f <- lmFit(E, d); 1 - rowSums((E - f$coefficients %*% t(d))^2) / rowSums((E - rowMeans(E))^2) }
  say("  median R2 (OLS): pooled %.3f ; levels %.3f", median(r2(pooled_form)), median(r2(level_form)))
}

## ---- PRJNA811175 : moderate / severe, 3 timepoints, n=1 per cell -------------
analyse("PRJNA811175",
  function(s, m) {
    tr <- trimws(m$treatment)
    s$time  <- factor(trimws(m$time))
    s$level <- factor(ifelse(grepl("^Control$",tr),"control",
                      ifelse(grepl("moderate",tr),"moderate","severe")),
                      levels=c("control","moderate","severe"))
    s },
  ~ arm + genotype + time, ~ level + genotype + time,
  c("levelmoderate","levelsevere"))

## ---- PRJNA776107 : 4 / 8 / 16 / 32 h after withholding water ----------------
analyse("PRJNA776107",
  function(s, m) {
    tr <- trimws(m$treatment)
    h  <- ifelse(grepl("^Control$",tr), "0", sub("^Drought For ([0-9]+) Hours$","\\1", tr))
    s$level <- factor(paste0("h", h), levels=paste0("h", c(0,4,8,16,32)))
    s },
  ~ arm + genotype, ~ level + genotype,
  c("levelh4","levelh8","levelh16","levelh32"))
