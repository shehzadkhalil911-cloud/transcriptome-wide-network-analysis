## 07_delta_analysis.R
## Patient-level tumour-minus-normal pathway score changes and their relationship.
## DeltaFerroptosis vs DeltaNecroptosis across patients (Analysis 5).
## Outputs:
##   results/crosstalk/patient_delta_scores.tsv
##   results/crosstalk/patient_delta_crosstalk.csv

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

tum <- s47[s47$tissue == "tumor", ]
nor <- s47[s47$tissue == "normal", ]
m <- merge(tum, nor, by = "patient_id", suffixes = c("_tumor", "_normal"))
cat("n paired patients:", nrow(m), "\n")

methods <- c("singscore", "gsva")
delta_tab <- data.frame(patient_id = m$patient_id, stringsAsFactors = FALSE)
for (meth in methods) {
  for (pw in c("ferroptosis", "necroptosis")) {
    cname <- paste0(pw, "_full_", meth)
    delta_tab[[paste0("Delta", pw, "_", meth)]] <-
      m[[paste0(cname, "_tumor")]] - m[[paste0(cname, "_normal")]]
  }
}
.save_tsv(delta_tab, "results_crosstalk", "patient_delta_scores.tsv")

# correlations
rows <- list()
for (meth in methods) {
  dF <- delta_tab[[paste0("Deltaferroptosis_", meth)]]
  dN <- delta_tab[[paste0("Deltanecroptosis_", meth)]]
  for (cm in c("pearson", "spearman")) {
    c1 <- cor_test_ci(dF, dN, method = cm, n_boot = n_boot, seed = cfg$project$seed)
    # permutation p
    set.seed(cfg$project$seed)
    perm_r <- numeric(n_perm)
    for (k in seq_len(n_perm)) {
      perm_r[k] <- suppressWarnings(cor(dF, sample(dN), method = cm))
    }
    perm_p <- (sum(abs(perm_r) >= abs(c1$r)) + 1) / (n_perm + 1)
    rows[[length(rows) + 1]] <- data.frame(
      method = meth, correlation_type = cm, n = c1$n, r = c1$r,
      ci95_lower = c1$ci[1], ci95_upper = c1$ci[2],
      p_value = c1$p, permutation_p = perm_p, stringsAsFactors = FALSE
    )
  }
}
delta_res <- do.call(rbind, rows)
delta_res$bh_fdr <- p.adjust(delta_res$p_value, method = "BH")
cat("Patient-level delta crosstalk:\n"); print(delta_res)
.save_tsv(delta_res, "results_crosstalk", "patient_delta_crosstalk.csv")

# leave-one-patient-out influence on Pearson delta correlation
influence <- list()
for (meth in methods) {
  dF <- delta_tab[[paste0("Deltaferroptosis_", meth)]]
  dN <- delta_tab[[paste0("Deltanecroptosis_", meth)]]
  r_all <- cor(dF, dN)
  out <- sapply(seq_along(dF), function(i) cor(dF[-i], dN[-i]))
  influence[[meth]] <- data.frame(
    method = meth, r_all = r_all, min_loo = min(out), max_loo = max(out),
    sd_loo = sd(out), stringsAsFactors = FALSE
  )
}
inf <- do.call(rbind, influence)
cat("Leave-one-out influence:\n"); print(inf)
.save_tsv(inf, "results_crosstalk", "delta_correlation_leave_one_out.csv")
.log("Delta analysis complete")
cat("DONE\n")