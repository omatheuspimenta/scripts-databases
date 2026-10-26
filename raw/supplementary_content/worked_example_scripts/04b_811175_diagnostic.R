#!/usr/bin/env Rscript
## DIAGNOSTIC - why does PRJNA811175 return zero DE genes?
## Re-fits the SAME counts with progressively richer designs to localise the cause.
## Usage: Rscript 04b_811175_diagnostic.R [ROOT]
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/home/jorge/sweet_recycler_fixes/13_comment/example"
P    <- "PRJNA811175"
FDR  <- 0.05; LFC <- 1
say  <- function(...) cat(sprintf(...), "\n", sep = "")

samples <- fread(file.path(ROOT,"out","samples_final.tsv"), sep="\t", colClasses="character")
meta    <- fread(file.path(ROOT,"..","..","sweet_recycler_RNA.csv"), colClasses="character")
cnt     <- fread(file.path(ROOT,"out","counts_lengthscaled.tsv.gz"), sep="\t")
mat     <- as.matrix(cnt[,-1,with=FALSE]); rownames(mat) <- cnt[["gene_id"]]
stopifnot(all(samples$run %in% colnames(mat)))   # matrix may hold a superset
univ <- readLines(file.path(ROOT,"out","de","universe.txt"))

s <- samples[bioproject == P]
m <- meta[match(s$run, meta$Run)]
stopifnot(identical(m$Run, s$run))

## recover the hidden within-arm factors from the released metadata
s$time      <- factor(trimws(m$time))                                    # 3 / 4 / 5 months after planting
tr          <- trimws(m$treatment)
s$severity  <- factor(ifelse(grepl("^Control$", tr), "control",
                      ifelse(grepl("moderate",  tr), "moderate", "severe")),
                      levels = c("control","moderate","severe"))
s$days      <- factor(ifelse(s$severity=="control", "0",
                             sub("^([0-9]+) days.*$","\\1", tr)))        # 0/30/60/90 d of treatment
s$arm       <- factor(s$arm, levels = c("control","drought"))
s$genotype  <- factor(s$genotype)

say("=== %s  n=%d  (%d drought / %d control) ===", P, nrow(s),
    sum(s$arm=="drought"), sum(s$arm=="control"))
say("\n-- what the pooled 'drought' arm actually contains --")
print(table(severity = s$severity, time = s$time))
print(table(arm = s$arm, time = s$time))
print(table(arm = s$arm, genotype = s$genotype))
print(table(severity = s$severity, genotype = s$genotype))
say("\ncells of severity x time x genotype (design cell sizes):")
print(addmargins(table(s$severity, s$time, s$genotype)))

## one filter / one normalisation, shared by every model, so only the design changes
base_design <- model.matrix(~ arm + genotype, data = s)
dge  <- DGEList(counts = mat[, s$run, drop = FALSE], samples = s)
keep <- filterByExpr(dge, design = base_design)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)
say("\nfilterByExpr kept %d genes; %d of them are in the common testable universe",
    sum(keep), sum(rownames(dge) %in% univ))

run_model <- function(label, form, coefs) {
  design <- model.matrix(form, data = s)
  rk <- qr(design)$rank
  say("\n---------------------------------------------------------------")
  say("MODEL %s :  %s", label, paste(deparse(form), collapse=""))
  say("  rank %d / ncol %d   residual df = %d", rk, ncol(design), nrow(s) - rk)
  if (rk < ncol(design)) { say("  RANK DEFICIENT - not estimable, skipped"); return(invisible(NULL)) }
  v   <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  say("  sigma (residual SD) median %.4f   prior df %.2f", median(fit$sigma), fit$df.prior)
  for (cf in coefs) {
    if (is.character(cf) && !(cf %in% colnames(design))) { say("  coef %s absent", cf); next }
    tt  <- topTable(fit, coef = cf, number = Inf, sort.by = "P")
    sig <- tt$adj.P.Val < FDR & abs(tt$logFC) > LFC
    inu <- rownames(tt) %in% univ
    say("  coef %-22s : up %5d  down %5d   (universe: up %5d down %5d)   min P %.3g  min BH %.3g",
        cf, sum(sig & tt$logFC>0), sum(sig & tt$logFC<0),
        sum(sig & tt$logFC>0 & inu), sum(sig & tt$logFC<0 & inu),
        min(tt$P.Value), min(tt$adj.P.Val))
  }
  invisible(fit)
}

## M0 - exactly the model used in the paper
run_model("0  (as published)", ~ arm + genotype, c("armdrought"))
## M1 - same contrast, but plant age / phenological stage blocked out
run_model("1  + time block",  ~ arm + genotype + time, c("armdrought"))
## M2 - stop pooling the two stress intensities
run_model("2  severity, time blocked", ~ severity + genotype + time,
          c("severitymoderate","severitysevere"))
## M3 - positive control: is there ANY signal this design can detect?
run_model("3  genotype (positive ctrl)", ~ arm + genotype, c("genotypeSP81-3250"))
run_model("4  time (positive ctrl)",     ~ arm + genotype + time, c("time4 months","time5 months"))

## how much variance does the within-arm structure carry?
say("\n---------------------------------------------------------------")
say("variance partition on the %d filtered genes (voom logCPM, base design)", nrow(dge))
v  <- voom(dge, base_design, plot = FALSE)
ss <- function(form) {
  d <- model.matrix(form, data = s)
  fit <- lmFit(v, d)
  res <- v$E - fit$coefficients %*% t(d)
  rowSums(res^2)
}
tot <- rowSums((v$E - rowMeans(v$E))^2)
r0  <- ss(~ arm + genotype)
r1  <- ss(~ arm + genotype + time)
r2  <- ss(~ severity + genotype + time)
say("  median R2  ~ arm+genotype            : %.3f", median(1 - r0/tot))
say("  median R2  ~ arm+genotype+time        : %.3f", median(1 - r1/tot))
say("  median R2  ~ severity+genotype+time   : %.3f", median(1 - r2/tot))
say("  median extra variance explained by time alone: %.3f", median((r0-r1)/tot))
