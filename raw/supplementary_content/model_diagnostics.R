#!/usr/bin/env Rscript
# Reviewer #4 comment 5 / Reviewer #2 comment 4: model specification, sample
# sizes, diagnostics, and whether repeated runs per BioProject and per genotype
# are adequately handled by the BioProject random intercept alone.
#
#   Rscript scripts/model_diagnostics.R        (run from 4_comment/)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(lme4)
})

OUT <- "tables"
GROUP_LEVELS <- c(
  "R570 (reference genotype)",
  "Saccharum hybrid cultivar (non-R570)",
  "Saccharum officinarum",
  "Saccharum spontaneum",
  "Tripidium arundinaceum (outgroup)"
)
RESPONSES <- c(percent_mapped = "mapping rate (%)",
               genes_detected_gt0 = "genes detected (TPM > 0)")

ann <- read_tsv("annotation/analysis_set.tsv", show_col_types = FALSE) %>%
  mutate(species_group = factor(species_group, levels = GROUP_LEVELS),
         log_depth = log10(num_mapped),
         bioproject = factor(bioproject),
         genotype_clean = factor(genotype_clean)) %>%
  filter(!is.na(percent_mapped), !is.na(num_mapped), num_mapped > 0)

model_df <- ann %>% filter(!is.na(fastqc_raw_avg_read_length)) %>% droplevels()
pass_df  <- model_df %>% filter(qc_pass == "PASS") %>% droplevels()

fit <- function(df, response, extra = "")
  lmer(as.formula(paste0(response,
        " ~ species_group + log_depth + fastqc_raw_avg_read_length + (1 | bioproject)",
        extra)),
       data = df, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4

# ── observations per group, per analysis set ──────────────────────────────────
census <- function(df, label)
  df %>% group_by(species_group) %>%
    summarise(n_runs = n(), n_bioprojects = n_distinct(bioproject),
              n_genotypes = n_distinct(genotype_clean), .groups = "drop") %>%
    mutate(analysis_set = label)

t27 <- bind_rows(
  census(ann,      "subset (Tables 16, 17, 20, 21)"),
  census(model_df, "mixed models (Tables 18, 19, 22)"),
  census(pass_df,  "sensitivity, qc_pass == PASS (Table 23)")
) %>% select(analysis_set, species_group, n_runs, n_bioprojects, n_genotypes)
write_tsv(t27, file.path(OUT, "SupplementaryTable27_model_sample_sizes.tsv"))
cat("== observations per species group ==\n")
print(as.data.frame(t27), row.names = FALSE)

cat(sprintf("\nsubset %d runs / %d projects; model set %d / %d; PASS set %d / %d\n",
            nrow(ann), n_distinct(ann$bioproject),
            nrow(model_df), n_distinct(model_df$bioproject),
            nrow(pass_df), n_distinct(pass_df$bioproject)))

# runs per BioProject and per genotype, to document the repeated-measures structure
rep_struct <- model_df %>%
  summarise(
    runs_per_project_median = median(table(bioproject)),
    runs_per_project_max = max(table(bioproject)),
    runs_per_genotype_median = median(table(genotype_clean)),
    runs_per_genotype_max = max(table(genotype_clean)),
    n_genotypes = n_distinct(genotype_clean))
gpp <- model_df %>% distinct(genotype_clean, bioproject) %>%
  count(genotype_clean, name = "n_projects")
cat(sprintf("\nruns/project median %.0f max %.0f; runs/genotype median %.0f max %.0f\n",
            rep_struct$runs_per_project_median, rep_struct$runs_per_project_max,
            rep_struct$runs_per_genotype_median, rep_struct$runs_per_genotype_max))
cat(sprintf("%d of %d genotypes occur in exactly one BioProject (%.1f%% of runs)\n",
            sum(gpp$n_projects == 1), nrow(gpp),
            100 * sum(model_df$genotype_clean %in%
                      gpp$genotype_clean[gpp$n_projects == 1]) / nrow(model_df)))

# ── diagnostics + genotype random effect ──────────────────────────────────────
diag_rows <- list(); comp_rows <- list()
for (resp in names(RESPONSES)) {
  m  <- fit(model_df, resp)
  mg <- fit(model_df, resp, " + (1 | genotype_clean)")

  r  <- resid(m); f <- fitted(m); rs <- r / sd(r)
  blup <- ranef(m)$bioproject[, 1]
  vc <- as.data.frame(VarCorr(m))
  v_proj <- vc$vcov[vc$grp == "bioproject"]; v_res <- vc$vcov[vc$grp == "Residual"]
  # Nakagawa marginal/conditional R2
  v_fix <- var(as.vector(model.matrix(m) %*% fixef(m)))
  # VIF from the fixed-effect design matrix
  X <- model.matrix(m)[, -1, drop = FALSE]
  vifs <- sapply(seq_len(ncol(X)), function(j)
    1 / (1 - summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared))
  names(vifs) <- colnames(X)
  # heteroscedasticity: regress squared residual on fitted value
  bp <- summary(lm(I(r^2) ~ f))
  vcg <- as.data.frame(VarCorr(mg))

  diag_rows[[resp]] <- data.frame(
    response = resp,
    n_obs = nobs(m), n_bioprojects = ngrps(m)[["bioproject"]],
    converged = is.null(m@optinfo$conv$lme4$messages),
    singular = isSingular(m),
    var_bioproject = round(v_proj, 4), var_residual = round(v_res, 4),
    icc_bioproject = round(v_proj / (v_proj + v_res), 4),
    r2_marginal = round(v_fix / (v_fix + v_proj + v_res), 4),
    r2_conditional = round((v_fix + v_proj) / (v_fix + v_proj + v_res), 4),
    resid_skewness = round(skew(r), 3), resid_kurtosis = round(kurt(r), 3),
    resid_abs_gt3sd_pct = round(100 * mean(abs(rs) > 3), 2),
    resid_vs_fitted_cor = round(cor(abs(r), f), 3),
    heterosced_p = signif(coef(bp)[2, 4], 3),
    blup_skewness = round(skew(blup), 3), blup_kurtosis = round(kurt(blup), 3),
    blup_normal_qq_cor = round(cor(sort(blup), qnorm(ppoints(length(blup)))), 4),
    max_vif = round(max(vifs), 2),
    var_genotype_added = round(vcg$vcov[vcg$grp == "genotype_clean"], 4),
    icc_genotype_added = round(vcg$vcov[vcg$grp == "genotype_clean"] /
                               sum(vcg$vcov), 4),
    aic_bioproject_only = round(AIC(m), 1), aic_plus_genotype = round(AIC(mg), 1))

  # species-group coefficients under both random-effect structures
  for (tm in grep("^species_group", names(fixef(m)), value = TRUE)) {
    ci  <- confint(m,  method = "Wald", parm = tm)
    cig <- confint(mg, method = "Wald", parm = tm)
    comp_rows[[paste(resp, tm)]] <- data.frame(
      response = resp, term = sub("^species_group", "", tm),
      est_bioproject_only = round(fixef(m)[[tm]], 1),
      ci_low = round(ci[1], 1), ci_high = round(ci[2], 1),
      est_plus_genotype = round(fixef(mg)[[tm]], 1),
      ci_low_g = round(cig[1], 1), ci_high_g = round(cig[2], 1))
  }
  cat(sprintf("\n== %s ==\n", RESPONSES[[resp]]))
  print(as.data.frame(t(diag_rows[[resp]])))
}

t28 <- bind_rows(diag_rows)
write_tsv(t28, file.path(OUT, "SupplementaryTable28_model_diagnostics.tsv"))
t29 <- bind_rows(comp_rows)
write_tsv(t29, file.path(OUT, "SupplementaryTable29_genotype_random_effect.tsv"))
cat("\n== species-group effects: BioProject intercept vs + genotype intercept ==\n")
print(as.data.frame(t29), row.names = FALSE)

# ── leave-one-BioProject-out influence on the species-group coefficients ──────
infl <- list()
for (resp in names(RESPONSES)) {
  base <- fixef(fit(model_df, resp))
  terms <- grep("^species_group", names(base), value = TRUE)
  projs <- levels(model_df$bioproject)
  mat <- matrix(NA_real_, length(projs), length(terms),
                dimnames = list(projs, terms))
  for (p in projs) {
    d <- droplevels(model_df[model_df$bioproject != p, ])
    if (n_distinct(d$species_group) < nlevels(model_df$species_group)) next
    co <- try(fixef(fit(d, resp)), silent = TRUE)
    if (inherits(co, "try-error")) next
    mat[p, ] <- co[terms] - base[terms]
  }
  infl[[resp]] <- data.frame(
    response = resp, term = sub("^species_group", "", terms),
    estimate = round(base[terms], 1),
    max_abs_change = round(apply(abs(mat), 2, max, na.rm = TRUE), 1),
    worst_bioproject = projs[apply(abs(mat), 2, which.max)],
    n_refits = colSums(!is.na(mat)))
}
t30 <- bind_rows(infl)
write_tsv(t30, file.path(OUT, "SupplementaryTable30_leave_one_project_out.tsv"))
cat("\n== leave-one-BioProject-out change in species-group coefficients ==\n")
print(as.data.frame(t30), row.names = FALSE)

cat("\ndone\n")
