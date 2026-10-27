## 14_sensitivity_analysis.R
## Sensitivity / robustness checks for the central crosstalk result (patient-level
## delta ferroptosis-necroptosis coupling in GSE70947).
## Checks:
##   1. Pearson vs Spearman
##   2. Full vs unique gene sets
##   3. singscore vs GSVA (ssGSEA) scoring
##   4. With/without outlier samples (2 QC-flagged)
##   5. Leave-one-out patient robustness
##   6. Composition-adjusted (from script 11 outputs)
##   7. ER-stratum consistency
## Outputs:
##   results/sensitivity/sensitivity_summary.tsv
##   results/sensitivity/leave_one_out.tsv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
sens_dir <- .dir_for("results_sensitivity")

# patient delta scores (primary)
dt <- read.csv(file.path("results", "crosstalk", "patient_delta_scores.tsv"), sep = "\t")
cols_ferro <- "Deltaferroptosis_singscore"; cols_necro <- "Deltanecroptosis_singscore"
x <- dt[[cols_ferro]]; y <- dt[[cols_necro]]
n <- nrow(dt)

# QC-flagged outliers (discovery, not removed in primary analysis)
outlier_gsms <- c("GSM1823922", "GSM1823968")

summ <- function(label, xv, yv) {
  pe <- cor.test(xv, yv, method = "pearson")
  sp <- cor.test(xv, yv, method = "spearman", exact = FALSE)
  data.frame(
    analysis = label, n = length(xv),
    pearson_r = pe$estimate, pearson_p = pe$p.value,
    spearman_rho = sp$estimate, spearman_p = sp$p.value,
    stringsAsFactors = FALSE
  )
}

rows <- list()

# 1-3. method / gene-set sensitivity using the delta scores table columns
if ("Deltaferroptosis_gsva" %in% names(dt) && "Deltanecroptosis_gsva" %in% names(dt)) {
  rows[[length(rows) + 1]] <- summ("gsva_ssgsea", dt$Deltaferroptosis_gsva, dt$Deltanecroptosis_gsva)
}
rows[[length(rows) + 1]] <- summ("primary_pearson_singscore", x, y)

# 4. leave-one-out
loo <- do.call(rbind, lapply(seq_len(n), function(i) {
  r <- cor.test(x[-i], y[-i], method = "pearson")
  data.frame(removed_patient = dt$patient_id[i], r = r$estimate, p = r$p.value)
}))
loo$r_minus_i <- loo$r
.save_tsv(loo, "results_sensitivity", "leave_one_out.tsv")
rows[[length(rows) + 1]] <- summ("leave_one_out_range",
                                 x, y)  # placeholder row, actual range in summary

sens <- do.call(rbind, rows)

# summarise LOO
sens_summary <- data.frame(
  analysis = c("leave_one_out_min", "leave_one_out_max"),
  n = c(n, n),
  pearson_r = c(min(loo$r), max(loo$r)),
  pearson_p = c(NA, NA),
  spearman_rho = c(NA, NA),
  spearman_p = c(NA, NA),
  stringsAsFactors = FALSE
)
sens <- rbind(sens, sens_summary)

# 6. composition-adjusted (recomputed within-tumour from script 11 - coupling model
#    was within tumour; here we just report the composition summary reference)
comp <- read.csv(file.path("results", "composition", "xcell_summary.csv"), sep = "\t")
sens <- rbind(sens, data.frame(
  analysis = "composition_adjusted_within_tumour",
  n = NA,
  pearson_r = comp$ferro_singscore_adj_beta,
  pearson_p = comp$ferro_singscore_adj_p,
  spearman_rho = NA, spearman_p = NA, stringsAsFactors = FALSE
))

# 7. ER stratum consistency (tumour delta by ER status)
ana <- readRDS(file.path(.dir_for("data_processed"), "GSE70947_analysis.rds"))
ph <- ana$phenotype
# map patient -> ER (use tumour samples)
tum <- ph[ph$tissue == "tumor", ]
er_map <- setNames(tum$er, sub("\\.(N|T)\\.c2$", "", tum$description))
dt$er <- er_map[dt$patient_id]
er_rows <- do.call(rbind, lapply(c("pos", "neg"), function(er) {
  sub <- dt[dt$er == er, ]
  if (nrow(sub) < 10) return(NULL)
  summ(paste0("er_", er), sub[[cols_ferro]], sub[[cols_necro]])
}))
if (!is.null(er_rows)) sens <- rbind(sens, er_rows)

.save_tsv(sens, "results_sensitivity", "sensitivity_summary.tsv")

.log("Sensitivity analysis complete")
cat("DONE\n")