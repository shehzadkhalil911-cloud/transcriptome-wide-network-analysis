## 06_crosstalk_analysis.R
## (A) Paired tumour vs normal pathway-score comparisons (Analysis 2 outputs).
## (B) Within-tissue ferroptosis-necroptosis correlations (Analysis 4A).
## (C) Differential coupling: patient-aware interaction models + robustness
##     (permutation, bootstrap, Fisher-z as secondary) (Analysis 4B).
## Outputs:
##   results/crosstalk/pathway_score_paired_comparisons.csv
##   results/crosstalk/within_tissue_correlations.csv
##   results/crosstalk/coupling_interaction_models.csv
##   results/crosstalk/coupling_robustness.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
ct_dir <- .dir_for("results_crosstalk")
data_proc <- .dir_for("data_processed")
set.seed(cfg$project$seed)
n_perm <- cfg$analysis$crosstalk$n_permutations
n_boot <- cfg$analysis$crosstalk$n_bootstrap

sc <- readRDS(file.path(data_proc, "pathway_scores_all.rds"))
s47 <- sc$scores_GSE70947
s47$tissue <- factor(s47$tissue, levels = c("normal", "tumor"))

pairs <- unique(s47$patient_id)
stopifnot(length(pairs) == 148)

methods <- c("singscore", "gsva")
pathways <- c("ferroptosis", "necroptosis")

# ---------------- (A) paired score comparisons ----------------
comp <- list()
for (m in methods) {
  for (pw in pathways) {
    coln <- paste0(pw, "_full_", m)
    tum <- s47[s47$tissue == "tumor", c("patient_id", coln)]
    nor <- s47[s47$tissue == "normal", c("patient_id", coln)]
    d <- merge(tum, nor, by = "patient_id")
    diff <- d[[paste0(coln, ".x")]] - d[[paste0(coln, ".y")]]
    tt <- t.test(diff)
    wt <- wilcox.test(diff, exact = FALSE)
    dz <- mean(diff) / sd(diff)
    ci <- tt$conf.int
    comp[[length(comp) + 1]] <- data.frame(
      pathway = pw, method = m, n_pairs = nrow(d),
      mean_tumour = mean(d[[paste0(coln, ".x")]]), mean_normal = mean(d[[paste0(coln, ".y")]]),
      mean_difference = mean(diff), sd_difference = sd(diff),
      effect_size_dz = dz,
      ci95_lower = ci[1], ci95_upper = ci[2],
      paired_t_pvalue = tt$p.value,
      paired_wilcoxon_pvalue = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}
comp <- do.call(rbind, comp)
# BH correction across the family of paired tests
comp$bh_fdr <- p.adjust(comp$paired_t_pvalue, method = "BH")
cat("Paired score comparisons:\n"); print(comp[, c("pathway", "method", "mean_difference", "paired_t_pvalue", "bh_fdr")])
.save_tsv(comp, "results_crosstalk", "pathway_score_paired_comparisons.csv")

# ---------------- (B) within-tissue correlations ----------------
cor_tab <- list()
for (m in methods) {
  fcol <- paste0("ferroptosis_full_", m)
  ncol_ <- paste0("necroptosis_full_", m)
  for (tis in c("normal", "tumor")) {
    sub <- s47[s47$tissue == tis, ]
    for (cm in c("pearson", "spearman")) {
      c1 <- cor_test_ci(sub[[fcol]], sub[[ncol_]], method = cm, n_boot = n_boot, seed = cfg$project$seed)
      cor_tab[[length(cor_tab) + 1]] <- data.frame(
        tissue = tis, method_score = m, correlation_type = cm,
        n = c1$n, r = c1$r, ci95_lower = c1$ci[1], ci95_upper = c1$ci[2],
        p_value = c1$p, stringsAsFactors = FALSE
      )
    }
  }
}
cor_tab <- do.call(rbind, cor_tab)
cor_tab$bh_fdr <- p.adjust(cor_tab$p_value, method = "BH")
cat("Within-tissue correlations:\n"); print(cor_tab[, c("tissue", "method_score", "correlation_type", "r", "p_value", "bh_fdr")])
.save_tsv(cor_tab, "results_crosstalk", "within_tissue_correlations.csv")

# ---------------- (C) differential coupling ----------------
suppressPackageStartupMessages(library(lme4))
inter <- list()
rob <- list()
for (m in methods) {
  fcol <- paste0("ferroptosis_full_", m)
  ncol_ <- paste0("necroptosis_full_", m)
  d <- s47[, c("patient_id", "tissue", fcol, ncol_)]
  d <- d[complete.cases(d), ]
  # grand-mean centre necroptosis score (not within-patient: that would make
  # it collinear with patient-by-tissue structure in fixed-effects fits)
  d[[ncol_]] <- d[[ncol_]] - mean(d[[ncol_]])
  dat <- d
  names(dat)[names(dat) == fcol] <- "Fer"
  names(dat)[names(dat) == ncol_] <- "Nec"

  # primary: mixed model with patient random intercept
  fit <- lmer(Fer ~ Nec * tissue + (1 | patient_id), data = dat)
  b <- fixef(fit)[["Nec:tissuetumor"]]
  se <- sqrt(diag(vcov(fit)))[["Nec:tissuetumor"]]
  pval <- 2 * pnorm(-abs(b / se))
  ci <- b + c(-1.96, 1.96) * se

  # fixed-effects lm equivalent (patient as factor) for robustness speed
  fit_lm <- lm(Fer ~ Nec * tissue + patient_id, data = dat)
  s <- summary(fit_lm)$coefficients
  b_lm <- s["Nec:tissuetumor", "Estimate"]
  se_lm <- s["Nec:tissuetumor", "Std. Error"]
  p_lm <- s["Nec:tissuetumor", "Pr(>|t|)"]
  ci_lm <- confint(fit_lm)["Nec:tissuetumor", ]

  inter[[length(inter) + 1]] <- data.frame(
    method = m, model = "lme4_lmer_rm_intercept", interaction_beta = b,
    se = se, p_value = pval, ci95_lower = ci[1], ci95_upper = ci[2],
    n_samples = nrow(dat), stringsAsFactors = FALSE
  )
  inter[[length(inter) + 1]] <- data.frame(
    method = m, model = "lm_fixed_patient", interaction_beta = b_lm,
    se = se_lm, p_value = p_lm, ci95_lower = ci_lm[1], ci95_upper = ci_lm[2],
    n_samples = nrow(dat), stringsAsFactors = FALSE
  )

  # robustness: within-patient tissue-swap permutation of interaction
  cat("permuting interaction (", n_perm, ") for ", m, " ...\n", sep = "")
  perm_b <- numeric(n_perm)
  set.seed(cfg$project$seed)
  for (k in seq_len(n_perm)) {
    dk <- dat
    pats <- unique(dk$patient_id)
    swapped <- sample(c(TRUE, FALSE), length(pats), replace = TRUE)
    names(swapped) <- pats
    swap_flag <- swapped[as.character(dk$patient_id)]
    dk$tissue <- ifelse(swap_flag,
                        ifelse(as.character(dk$tissue) == "tumor", "normal", "tumor"),
                        as.character(dk$tissue))
    dk$tissue <- factor(dk$tissue, levels = c("normal", "tumor"))
    fk <- lm(Fer ~ Nec * tissue + patient_id, data = dk)
    perm_b[k] <- coef(fk)[["Nec:tissuetumor"]]
  }
  perm_p <- (sum(abs(perm_b) >= abs(b_lm)) + 1) / (n_perm + 1)

  # bootstrap over patients
  set.seed(cfg$project$seed)
  boot_b <- numeric(n_boot)
  pats <- unique(dat$patient_id)
  for (k in seq_len(n_boot)) {
    ps <- sample(pats, length(pats), replace = TRUE)
    dk <- do.call(rbind, lapply(ps, function(p) {
      r <- dat[dat$patient_id == p, ]
      if (nrow(r) == 0) return(NULL)
      r
    }))
    fk <- try(lm(Fer ~ Nec * tissue + patient_id, data = dk), silent = TRUE)
    boot_b[k] <- if (inherits(fk, "try-error")) NA else coef(fk)[["Nec:tissuetumor"]]
  }
  boot_b <- boot_b[is.finite(boot_b)]
  boot_ci <- quantile(boot_b, c(0.025, 0.975))

  # Fisher-z secondary comparison (tumour vs normal correlation)
  rz_t <- atanh(cor(dat$Fer[dat$tissue == "tumor"], dat$Nec[dat$tissue == "tumor"]))
  rz_n <- atanh(cor(dat$Fer[dat$tissue == "normal"], dat$Nec[dat$tissue == "normal"]))
  n_t <- sum(dat$tissue == "tumor"); n_n <- sum(dat$tissue == "normal")
  fz <- (rz_t - rz_n) / sqrt(1 / (n_t - 3) + 1 / (n_n - 3))
  fz_p <- 2 * pnorm(-abs(fz))

  rob[[length(rob) + 1]] <- data.frame(
    method = m,
    permutation_interaction_beta = b_lm, permutation_p = perm_p,
    bootstrap_interaction_beta = mean(boot_b), bootstrap_ci_lower = boot_ci[1], bootstrap_ci_upper = boot_ci[2],
    fisher_z = fz, fisher_z_p = fz_p,
    stringsAsFactors = FALSE
  )
}
inter_tab <- do.call(rbind, inter)
inter_tab$bh_fdr <- p.adjust(inter_tab$p_value, method = "BH")
rob_tab <- do.call(rbind, rob)
cat("Interaction models:\n"); print(inter_tab)
cat("Robustness:\n"); print(rob_tab)
.save_tsv(inter_tab, "results_crosstalk", "coupling_interaction_models.csv")
.save_tsv(rob_tab, "results_crosstalk", "coupling_robustness.csv")
.log("Crosstalk analysis complete")
cat("DONE\n")