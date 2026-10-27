## 11_composition_analysis.R
## Cell-type composition analysis for the discovery cohort (GSE70947) using xCell.
## Steps:
##   1. Run xCell on the imputed gene-level matrix (gene symbols -> xCell signatures)
##   2. Report tumour vs adjacent-normal composition changes (unpaired t-test + BH)
##   3. Re-evaluate ferroptosis/necroptosis pathway coupling (within-tumour regression)
##      with adjustment for cell-type proportions to test robustness of crosstalk signal
## Outputs:
##   results/composition/xcell_scores_tumour.tsv, xcell_scores_normal.tsv
##   results/composition/cell_type_paired_de.tsv
##   results/composition/composition_adjusted_coupling.tsv
##   results/composition/xcell_summary.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
comp_dir <- .dir_for("results_composition")

# ---- prepare expression ----
ana <- readRDS(file.path(.dir_for("data_processed"), "GSE70947_analysis.rds"))
expr <- ana$expression_imputed
pheno <- ana$phenotype
stopifnot(all(colnames(expr) == rownames(pheno)))

tumour_cols <- rownames(pheno)[pheno$tissue == "tumor"]
normal_cols <- rownames(pheno)[pheno$tissue == "normal"]

# ---- xCell (cell types in rows, samples in cols) ----
library(xCell)
scr_t <- xCellAnalysis(expr[, tumour_cols])
scr_n <- xCellAnalysis(expr[, normal_cols])
common_ct <- intersect(rownames(scr_t), rownames(scr_n))

saveRDS(list(tumour = scr_t, normal = scr_n), file.path(comp_dir, "xcell_scores.rds"))

# write in sample x celltype orientation with explicit sample id column
write_scores <- function(scr, samples, fname) {
  d <- as.data.frame(t(scr))
  d <- data.frame(sample_id = samples, d, check.names = FALSE)
  .save_tsv(d, "results_composition", fname)
}
write_scores(scr_t, tumour_cols, "xcell_scores_tumour.tsv")
write_scores(scr_n, normal_cols, "xcell_scores_normal.tsv")

# ---- paired DE for cell types (tumour vs adjacent normal, same patients) ----
# map sample ids back to patient pairing key via phenotype$description (CM016.N.c2)
gsm_to_desc <- function(gsm_vec) ana$phenotype[gsm_vec, "description"]
pkey <- function(sid) sub("\\.(N|T)\\.c2$", "", sid)
tumour_key <- pkey(gsm_to_desc(tumour_cols))
normal_key <- pkey(gsm_to_desc(normal_cols))
pair_ok <- intersect(tumour_key, normal_key)
ct_de <- do.call(rbind, lapply(common_ct, function(ct) {
  vt <- scr_t[ct, match(pair_ok, tumour_key)]
  vn <- scr_n[ct, match(pair_ok, normal_key)]
  d <- vt - vn
  if (sd(d) == 0 || is.na(sd(d))) return(NULL)
  tt <- t.test(d)
  data.frame(celltype = ct, mean_tumour = mean(vt), mean_normal = mean(vn),
             mean_diff = mean(d), sd_diff = sd(d), t_stat = tt$statistic,
             df = tt$parameter, p_value = tt$p.value,
             stringsAsFactors = FALSE)
}))
ct_de$fdr <- p.adjust(ct_de$p_value, method = "BH")
ct_de <- ct_de[order(ct_de$p_value), ]
.save_tsv(ct_de, "results_composition", "cell_type_paired_de.tsv")

# ---- composition-adjusted coupling (within-tumour) ----
ps <- readRDS(file.path(.dir_for("data_processed"), "pathway_scores_all.rds"))
sc <- ps$scores_GSE70947
sc_t <- sc[sc$tissue == "tumor", ]
stopifnot(all(sc_t$sample_id %in% colnames(scr_t)))
tumour_cols <- sc_t$sample_id  # reorder to score order (GSM ids)

comps <- as.data.frame(t(scr_t))          # samples x celltypes
comps <- comps[tumour_cols, ]
comps <- comps[, setdiff(names(comps), c("ImmuneScore", "StromaScore", "MicroenvironmentScore"))]
ct_var <- apply(comps, 2, sd, na.rm = TRUE)
sel_ct <- names(ct_var)[ct_var > quantile(ct_var, 0.25, na.rm = TRUE)]
comp_adj <- scale(comps[, sel_ct, drop = FALSE])

coupling_lm <- function(dep, indep, comps2) {
  dat <- data.frame(y = dep, x = indep, comps2, check.names = FALSE)
  base <- lm(y ~ x, data = dat)
  adj <- lm(y ~ ., data = dat)
  base_b <- coef(base)[["x"]]; base_p <- summary(base)$coefficients["x", 4]
  adj_b <- coef(adj)[["x"]]; adj_p <- summary(adj)$coefficients["x", 4]
  c(base_beta = unname(base_b), base_p = unname(base_p),
    adj_beta = unname(adj_b), adj_p = unname(adj_p))
}

run_pair <- function(dep_name, indep_name) {
  r1 <- coupling_lm(sc_t[[dep_name]], sc_t[[indep_name]], comp_adj)
  r2 <- coupling_lm(sc_t[[indep_name]], sc_t[[dep_name]], comp_adj)
  data.frame(model = paste0(dep_name, " ~ ", indep_name), t(as.data.frame(r1)),
             model_rev = paste0(indep_name, " ~ ", dep_name), t(as.data.frame(r2)),
             check.names = FALSE)
}

out <- rbind(
  run_pair("ferroptosis_full_singscore", "necroptosis_full_singscore"),
  run_pair("ferroptosis_full_gsva", "necroptosis_full_gsva")
)
.save_tsv(out, "results_composition", "composition_adjusted_coupling.tsv")

summ <- data.frame(
  n_celltypes = nrow(ct_de),
  n_sig_celltypes_fdr05 = sum(ct_de$fdr < 0.05, na.rm = TRUE),
  top_celltype_diff = ct_de$celltype[which.max(abs(ct_de$mean_diff))],
  n_comp_celltypes_used = ncol(comp_adj),
  ferro_singscore_base_beta = out$base_beta[1],
  ferro_singscore_adj_beta = out$adj_beta[1],
  ferro_singscore_adj_p = out$adj_p[1],
  ferro_gsva_adj_beta = out$adj_beta[2],
  ferro_gsva_adj_p = out$adj_p[2]
)
.save_tsv(summ, "results_composition", "xcell_summary.csv")

.log("Composition analysis complete")
cat("DONE\n")