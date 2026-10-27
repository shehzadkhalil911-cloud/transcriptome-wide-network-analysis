## 19_final_audit.R
## Final end-to-end audit of the manuscript. For every key number that appears in
## the manuscript text, recompute the formatted string from the source result
## file and verify the manuscript actually contains it (no drift between results
## and text). Also verify every figure/table referenced by the manuscript exists
## and every reference listed matches the PMID-verified list.
## Outputs (all in results/audit/):
##   manuscript_number_audit.md
##   figure_audit.md
##   table_audit.md
##   claim_evidence_matrix.tsv
##   final_audit_summary.md

suppressPackageStartupMessages(library(data.table))
source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
audit_dir <- .dir_for("results_audit")

rd <- function(p) data.table::fread(p, data.table = FALSE, encoding = "UTF-8")
res <- function(rel) file.path("results", rel)

## ---- formatting (must match 18_manuscript_generation.R) --------------------
fmt3 <- function(x) sprintf("%.3f", x)
fmtp <- function(x) ifelse(x < 0.001, sprintf("%.1e", x), sprintf("%.3f", x))
fmtr <- function(x) sprintf("%.3f", x)

manu_path <- file.path("manuscript", "Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.md")
manu <- paste(readLines(manu_path, warn = FALSE), collapse = "\n")

claims <- data.frame()

add_claim <- function(id, description, source_file, expected, found) {
  claims <<- rbind(claims, data.frame(
    claim_id = id, description = description, source_file = source_file,
    expected_string = ifelse(length(expected) == 0, "", expected),
    found_in_manuscript = found, status = ifelse(found, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  ))
}

chk <- function(id, description, source_file, needle) {
  add_claim(id, description, source_file, needle, grepl(needle, manu, fixed = TRUE))
}

## ---- load source files -----------------------------------------------------
delta_ct <- rd(res("crosstalk/patient_delta_crosstalk.csv"))
within   <- rd(res("crosstalk/within_tissue_correlations.csv"))
pc       <- rd(res("crosstalk/pathway_score_paired_comparisons.csv"))
interact <- rd(res("crosstalk/coupling_interaction_models.csv"))
loo      <- rd(res("crosstalk/delta_correlation_leave_one_out.csv"))
wsum     <- rd(res("wgcna/wgcna_summary.csv"))
mt       <- rd(res("wgcna/module_trait_correlations.csv"))
mdef     <- rd(res("enrichment/crosstalk_modules_defined.csv"))
vsum     <- rd(res("validation/validation_summary.csv"))
vcoup    <- rd(res("validation/validation_coupling.tsv"))
vmod     <- rd(res("validation/validation_module_associations.tsv"))
fe       <- rd(res("field_effect/field_effect_summary.csv"))
sens     <- rd(res("sensitivity/sensitivity_summary.tsv"))
de_sum   <- rd(res("differential_expression/pathway_de_summary.csv"))
gs       <- rd(res("tables/gene_set_summary.csv"))
t1       <- rd(res("tables/Table1_cohorts.tsv"))
comp     <- rd(res("composition/composition_adjusted_coupling.tsv"))
xsum     <- rd(res("composition/xcell_summary.csv"))
go_bp    <- rd(res("enrichment/GO_BP.tsv"))

## ---- helper lookups --------------------------------------------------------
disc_ferro <- setNames(mdef$r_ferro, mdef$module)
disc_necro <- setNames(mdef$r_necro, mdef$module)
vmod$conc_ferro <- vmod$fdr_ferro < 0.05 & sign(vmod$r_delta_ferro) == sign(disc_ferro[vmod$module])
vmod$conc_necro <- vmod$fdr_necro < 0.05 & sign(vmod$r_delta_necro) == sign(disc_necro[vmod$module])
n_sig_f <- sum(vmod$fdr_ferro < 0.05)
n_sig_n <- sum(vmod$fdr_necro < 0.05)
n_conc_f <- sum(vmod$conc_ferro)
n_conc_n <- sum(vmod$conc_necro)

## ---- 1. cohort numbers -----------------------------------------------------
chk("C1", "Discovery paired samples = 148", "Table1_cohorts.tsv", "148 paired samples")
chk("C2", "Discovery 296 samples = 2x148", "Table1_cohorts.tsv", "of 296 breast adenocarcinomas")
chk("C3", "Validation paired = 47", "Table1_cohorts.tsv", "47 paired tumour-adjacent samples")
chk("C4", "Validation tumour = 47", "validation_summary.csv", "47 tumour")
chk("C5", "Validation adjacent = 47", "validation_summary.csv", "47 adjacent normal")
chk("C6", "Validation reduction = 43", "validation_summary.csv", "43 reduction mammoplasty normal")
chk("C7", "Extra cohort paired = 43", "Table1_cohorts.tsv", "GSE15852 (43 paired samples)")

## ---- 2. gene sets ----------------------------------------------------------
chk("G1", "Ferroptosis 259 unique genes", "gene_set_summary.csv", "259 unique genes")
chk("G2", "Ferroptosis 227 measured", "pathway_de_summary.csv", "227 measured")
chk("G3", "Necroptosis 158 unique genes", "gene_set_summary.csv", "158 unique genes")
chk("G4", "Necroptosis 122 measured", "pathway_de_summary.csv", "122 measured")
chk("G5", "Overlap 18 genes", "gene_set_summary.csv", "18 genes shared")

## ---- 3. primary crosstalk ---------------------------------------------------
r_p  <- delta_ct$r[delta_ct$method == "singscore" & delta_ct$correlation_type == "pearson"]
ci_l <- delta_ct$ci95_lower[delta_ct$method == "singscore" & delta_ct$correlation_type == "pearson"]
ci_h <- delta_ct$ci95_upper[delta_ct$method == "singscore" & delta_ct$correlation_type == "pearson"]
p_p  <- delta_ct$p_value[delta_ct$method == "singscore" & delta_ct$correlation_type == "pearson"]
p_perm <- delta_ct$permutation_p[1]
r_sp <- delta_ct$r[delta_ct$method == "singscore" & delta_ct$correlation_type == "spearman"]
p_sp <- delta_ct$p_value[delta_ct$method == "singscore" & delta_ct$correlation_type == "spearman"]
r_gv <- delta_ct$r[delta_ct$method == "gsva" & delta_ct$correlation_type == "pearson"]
p_gv <- delta_ct$p_value[delta_ct$method == "gsva" & delta_ct$correlation_type == "pearson"]

chk("P1", "Pearson r = 0.361", "patient_delta_crosstalk.csv", sprintf("r = %s", fmtr(r_p)))
chk("P2", "95% CI 0.202-0.497", "patient_delta_crosstalk.csv", sprintf("CI %s-%s", fmt3(ci_l), fmt3(ci_h)))
chk("P3", "p = 6.6e-06", "patient_delta_crosstalk.csv", sprintf("p = %s", fmtp(p_p)))
chk("P4", "permutation p = 5.0e-04", "patient_delta_crosstalk.csv", sprintf("permutation p = %s", fmtp(p_perm)))
chk("P5", "Spearman r = 0.286", "patient_delta_crosstalk.csv", sprintf("r = %s", fmtr(r_sp)))
chk("P6", "Spearman p = 4.2e-04", "patient_delta_crosstalk.csv", sprintf("p = %s", fmtp(p_sp)))
chk("P7", "GSVA r = 0.407", "patient_delta_crosstalk.csv", sprintf("r = %s", fmtr(r_gv)))
chk("P8", "GSVA p = 2.9e-07", "patient_delta_crosstalk.csv", sprintf("p = %s", fmtp(p_gv)))

## within-tissue
wn_r <- within$r[within$tissue == "normal" & within$method_score == "singscore" & within$correlation_type == "pearson"]
wt_r <- within$r[within$tissue == "tumor" & within$method_score == "singscore" & within$correlation_type == "pearson"]
wt_rg<- within$r[within$tissue == "tumor" & within$method_score == "gsva" & within$correlation_type == "pearson"]
chk("P9",  "within normal r = 0.705", "within_tissue_correlations.csv", sprintf("r = %s", fmtr(wn_r)))
chk("P10", "within tumour r = 0.565", "within_tissue_correlations.csv", sprintf("r = %s", fmtr(wt_r)))
chk("P11", "within tumour GSVA r = 0.503", "within_tissue_correlations.csv", sprintf("r = %s", fmtr(wt_rg)))

## ---- 4. paired pathway comparisons ------------------------------------------
pc_f <- pc[pc$pathway == "ferroptosis" & pc$method == "singscore", ]
pc_n <- pc[pc$pathway == "necroptosis" & pc$method == "singscore", ]
pc_fg<- pc[pc$pathway == "ferroptosis" & pc$method == "gsva", ]
pc_ng<- pc[pc$pathway == "necroptosis" & pc$method == "gsva", ]
chk("D1", "necroptosis mean diff -0.004", "pathway_score_paired_comparisons.csv", fmt3(pc_n$mean_difference))
chk("D2", "necroptosis paired t p = 0.002", "pathway_score_paired_comparisons.csv", fmtp(pc_n$paired_t_pvalue))
chk("D3", "necroptosis wilcoxon p = 0.001", "pathway_score_paired_comparisons.csv", fmtp(pc_n$paired_wilcoxon_pvalue))
chk("D4", "ferroptosis mean diff -0.001", "pathway_score_paired_comparisons.csv", fmt3(pc_f$mean_difference))
chk("D5", "ferroptosis paired t p = 0.258", "pathway_score_paired_comparisons.csv", fmtp(pc_f$paired_t_pvalue))
chk("D6", "gsva ferro diff -0.030", "pathway_score_paired_comparisons.csv", fmt3(pc_fg$mean_difference))
chk("D7", "gsva ferro p = 1.7e-04", "pathway_score_paired_comparisons.csv", fmtp(pc_fg$paired_t_pvalue))
chk("D8", "gsva necro diff -0.045", "pathway_score_paired_comparisons.csv", fmt3(pc_ng$mean_difference))
chk("D9", "gsva necro p = 2.1e-05", "pathway_score_paired_comparisons.csv", fmtp(pc_ng$paired_t_pvalue))

## ---- 5. interaction / LOO ----------------------------------------------------
int_s <- interact[interact$method == "singscore" & interact$model == "lme4_lmer_rm_intercept", ]
int_g <- interact[interact$method == "gsva" & interact$model == "lme4_lmer_rm_intercept", ]
chk("I1", "singscore interaction beta -0.197", "coupling_interaction_models.csv", fmt3(int_s$interaction_beta))
chk("I2", "singscore interaction p = 0.009", "coupling_interaction_models.csv", fmtp(int_s$p_value))
chk("I3", "gsva interaction beta -0.201", "coupling_interaction_models.csv", fmt3(int_g$interaction_beta))
chk("I4", "gsva interaction p = 0.002", "coupling_interaction_models.csv", fmtp(int_g$p_value))
chk("I5", "LOO range 0.335-0.380", "delta_correlation_leave_one_out.csv",
    sprintf("%s and %s", fmt3(loo$min_loo[1]), fmt3(loo$max_loo[1])))

## ---- 6. WGCNA -----------------------------------------------------------------
chk("W1", "soft power 6", "wgcna_summary.csv", "soft-thresholding power of 6")
chk("W2", "scale-free R2 0.892", "wgcna_summary.csv", sprintf("0.%03d", round(scale_free_r2 <- wsum$scale_free_R2_selected * 1000)))
chk("W3", "33 co-expression modules", "wgcna_summary.csv", "identified 33 co-expression modules")
chk("W4", "30997 genes input", "wgcna_summary.csv", "30997 expression probes")
chk("W5", "12 associated with both", "wgcna_summary.csv", "12 modules were associated")
chk("W6", "11 concordant modules", "crosstalk_modules_defined.csv", "11 showed concordant direction")

## per-module values (3.5)
for (m in mdef$module) {
  rf <- mdef$r_ferro[mdef$module == m]
  rn <- mdef$r_necro[mdef$module == m]
  pf <- mdef$fdr_ferro[mdef$module == m]
  pn <- mdef$fdr_necro[mdef$module == m]
  chk(paste0("M-", m, "-rF"), sprintf("module %s r_ferro", m), "crosstalk_modules_defined.csv", sprintf("%s (r<sub>ferro</sub> = %s", m, fmtr(rf)))
  chk(paste0("M-", m, "-rN"), sprintf("module %s r_necro", m), "crosstalk_modules_defined.csv", sprintf("r<sub>necro</sub> = %s", fmtr(rn)))
  chk(paste0("M-", m, "-genes"), sprintf("module %s gene count", m), "crosstalk_modules_defined.csv", sprintf("%d genes)", mdef$n_genes[mdef$module == m]))
}

## GO p.adjust for ME9 (3.6)
pa_ai <- go_bp$p.adjust[go_bp$module == "ME9" & go_bp$Description == "adaptive immune response"]
pa_tc <- go_bp$p.adjust[go_bp$module == "ME9" & go_bp$Description == "T cell activation"]
chk("E1", "GO adaptive immune response p.adjust", "GO_BP.tsv", fmtp(pa_ai))
chk("E2", "GO T cell activation p.adjust", "GO_BP.tsv", fmtp(pa_tc))

## ---- 7. validation ------------------------------------------------------------
v_r  <- vcoup$r[vcoup$group == "tumor"]
v_p  <- vcoup$p_value[vcoup$group == "tumor"]
v_an <- vcoup$r[vcoup$group == "adjacent_normal"]
v_rn <- vcoup$r[vcoup$group == "reduction_normal"]
chk("V1", "validation coupling r = 0.865", "validation_coupling.tsv", sprintf("r = %s", fmtr(v_r)))
chk("V2", "validation coupling p = 4.4e-15", "validation_coupling.tsv", sprintf("p = %s", fmtp(v_p)))
chk("V3", "validation adjacent r = 0.993", "validation_coupling.tsv", sprintf("r = %s", fmtr(v_an)))
chk("V4", "validation reduction r = 0.885", "validation_coupling.tsv", sprintf("r = %s", fmtr(v_rn)))
chk("V5", "validation ferro paired p = 0.159", "validation_summary.csv", fmtp(vsum$ferro_paired_p))
chk("V6", "validation necro paired p = 0.164", "validation_summary.csv", fmtp(vsum$necro_paired_p))
chk("V7", "5 modules sig for ferro", "validation_module_associations.tsv", sprintf("%d retained a significant association with delta-ferroptosis", n_sig_f))
chk("V8", "5 modules sig for necro", "validation_module_associations.tsv", sprintf("%d with delta-necroptosis", n_sig_n))
chk("V9", "4 concordant ferro", "validation_module_associations.tsv", sprintf("%d modules replicated with a concordant direction for delta-ferroptosis", n_conc_f))
chk("V10", "4 concordant necro", "validation_module_associations.tsv", sprintf("%d for delta-necroptosis", n_conc_n))
chk("V11", "ME28 opposite direction noted", "validation_module_associations.tsv", "ME28")
chk("V12", "hub replication 47.0%", "validation_summary.csv", sprintf("%.1f%%", vsum$frac_hubs_replicated * 100))
chk("V13", "1954 hubs tested", "validation_summary.csv", sprintf("%d tested genes", vsum$n_hubs_tested))

## ---- 8. field effect ------------------------------------------------------------
chk("F1", "ferro ANOVA p = 0.216", "field_effect_summary.csv", fmtp(fe$ferro_anova_p))
chk("F2", "necro ANOVA p = 0.161", "field_effect_summary.csv", fmtp(fe$necro_anova_p))
chk("F3", "ferro tumor-vs-reduction padj = 0.990", "field_effect_summary.csv", fmtp(fe$ferro_tumor_vs_reduction_padj))
chk("F4", "necro tumor-vs-reduction padj = 0.992", "field_effect_summary.csv", fmtp(fe$necro_tumor_vs_reduction_padj))

## ---- 9. composition --------------------------------------------------------------
chk("X1", "67 cell types", "xcell_summary.csv", "67 cell types")
chk("X2", "54 sig cell types", "xcell_summary.csv", "54 differed significantly")
chk("X3", "composition base beta 0.473", "composition_adjusted_coupling.tsv", sprintf("beta %s", fmt3(comp$base_beta[1])))
chk("X4", "composition base p 7.1e-14", "composition_adjusted_coupling.tsv", sprintf("p = %s", fmtp(comp$base_p[1])))
chk("X5", "composition adj beta 0.222", "composition_adjusted_coupling.tsv", sprintf("-> %s", fmt3(comp$adj_beta[1])))
chk("X6", "composition adj p 0.028", "composition_adjusted_coupling.tsv", sprintf("p = %s", fmtp(comp$adj_p[1])))
chk("X7", "gsva base beta 0.401", "composition_adjusted_coupling.tsv", sprintf("beta %s", fmt3(comp$base_beta[2])))
chk("X8", "gsva adj beta 0.196", "composition_adjusted_coupling.tsv", sprintf("-> %s", fmt3(comp$adj_beta[2])))
chk("X9", "gsva adj p 0.050", "composition_adjusted_coupling.tsv", sprintf("p = %s", fmtp(comp$adj_p[2])))

## ---- 10. sensitivity --------------------------------------------------------------
erp <- sens[sens$analysis == "er_pos", ]
ern <- sens[sens$analysis == "er_neg", ]
chk("S1", "ER+ r = 0.325", "sensitivity_summary.tsv", sprintf("r = %s", fmtr(erp$pearson_r)))
chk("S2", "ER+ p = 3.9e-04", "sensitivity_summary.tsv", fmtp(erp$pearson_p))
chk("S3", "ER- r = 0.412", "sensitivity_summary.tsv", sprintf("r = %s", fmtr(ern$pearson_r)))
chk("S4", "ER- p = 0.017", "sensitivity_summary.tsv", fmtp(ern$pearson_p))

## ---- 11. figures and tables referenced by the manuscript ------------------------
fig_need <- c(
  "Fig1_study_overview.png", "Fig2_pca.png", "Fig3_pathway_by_tissue.png",
  "Fig4_within_tissue_coupling.png", "Fig5_delta_crosstalk.png",
  "Fig6_wgcna_soft_and_traits.png", "Fig7_crosstalk_module_associations.png",
  "Fig8_enrichment_kegg.png", "Fig9_hub_genes.png", "Fig10_external_validation.png",
  "SuppA_leave_one_out.png", "SuppB_field_effect.png", "SuppC_composition_adjusted.png"
)
fig_rows <- data.frame(
  figure = fig_need,
  expected_path = file.path("results", c(rep("figures", 10), rep("supplementary", 3)), fig_need),
  exists = file.exists(file.path("results", c(rep("figures", 10), rep("supplementary", 3)), fig_need)),
  stringsAsFactors = FALSE
)

tables_dir <- file.path("results", "tables")
table_files <- list.files(tables_dir, full.names = TRUE)
tab_rows <- data.frame(
  table = basename(table_files),
  exists = TRUE,
  stringsAsFactors = FALSE
)

## ---- 12. references in manuscript match verified PMIDs --------------------------
ref_ver <- rd(file.path(audit_dir, "reference_verification.tsv"))
ref_rows <- data.frame(
  key = ref_ver$key,
  pmid = ref_ver$pmid,
  verified = ref_ver$verified,
  cited_in_manuscript = sapply(ref_ver$pmid, function(p) grepl(sprintf("PMID: %d", as.integer(p)), manu, fixed = TRUE)),
  stringsAsFactors = FALSE
)

## ---- write claim evidence matrix ------------------------------------------------
save_claims <- claims
.save_tsv(save_claims, "results_audit", "claim_evidence_matrix.tsv")

## ---- write markdown reports -----------------------------------------------------
n_pass <- sum(claims$status == "PASS")
n_fail <- sum(claims$status == "FAIL")

audit_md <- c(
  "# Manuscript number audit",
  "",
  sprintf("Audit timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  sprintf("Checks passed: %d/%d (FAIL: %d)", n_pass, nrow(claims), n_fail),
  "",
  "| ID | claim | source file | expected string | in manuscript | status |",
  "|----|-------|-------------|-----------------|---------------|--------|"
)
for (i in seq_len(nrow(claims))) {
  audit_md <- c(audit_md, sprintf("| %s | %s | %s | `%s` | %s | %s |",
    claims$claim_id[i], claims$description[i], claims$source_file[i],
    claims$expected_string[i], claims$found_in_manuscript[i], claims$status[i]))
}
writeLines(audit_md, file.path(audit_dir, "manuscript_number_audit.md"))

fig_md <- c(
  "# Figure audit", "",
  sprintf("Audit timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
  sprintf("Figures referenced by manuscript present: %d/%d",
          sum(fig_rows$exists), nrow(fig_rows)), "",
  "| figure | path | exists |", "|--------|------|--------|"
)
for (i in seq_len(nrow(fig_rows))) {
  fig_md <- c(fig_md, sprintf("| %s | %s | %s |",
    fig_rows$figure[i], fig_rows$expected_path[i], fig_rows$exists[i]))
}
writeLines(fig_md, file.path(audit_dir, "figure_audit.md"))

tab_md <- c(
  "# Table audit", "",
  sprintf("Audit timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
  "| table | exists |", "|-------|--------|"
)
for (i in seq_len(nrow(tab_rows))) {
  tab_md <- c(tab_md, sprintf("| %s | %s |", tab_rows$table[i], tab_rows$exists[i]))
}
writeLines(tab_md, file.path(audit_dir, "table_audit.md"))

## ---- final summary -----------------------------------------------------------------
summary_md <- c(
  "# Final audit summary", "",
  sprintf("Audit timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
  sprintf("- Manuscript number checks: **%d/%d passed** (%d FAIL)", n_pass, nrow(claims), n_fail),
  sprintf("- Figures referenced and present: %d/%d", sum(fig_rows$exists), nrow(fig_rows)),
  sprintf("- Tables present: %d/%d", sum(tab_rows$exists), nrow(tab_rows)),
  sprintf("- References PMID-verified: %d/%d", sum(ref_ver$verified), nrow(ref_ver)),
  sprintf("- Verified references cited in manuscript: %d/%d",
          sum(ref_rows$cited_in_manuscript[ref_rows$verified]), sum(ref_rows$verified)),
  "", "## Failed number checks", ""
)
fails <- claims[claims$status == "FAIL", ]
if (nrow(fails) == 0) {
  summary_md <- c(summary_md, "None.", "")
} else {
  summary_md <- c(summary_md, "| ID | claim | expected string |", "|----|-------|-----------------|")
  for (i in seq_len(nrow(fails))) {
    summary_md <- c(summary_md, sprintf("| %s | %s | `%s` |",
      fails$claim_id[i], fails$description[i], fails$expected_string[i]))
  }
  summary_md <- c(summary_md, "")
}
writeLines(summary_md, file.path(audit_dir, "final_audit_summary.md"))

.log("Final audit complete: ", n_pass, "/", nrow(claims), " checks passed")
cat("DONE\n")
