## 16_tables.R
## Assemble manuscript tables from computed result files.
## Tables:
##   Table 1  Cohort characteristics (discovery + validation)
##   Table 2  Differential expression summary (ferroptosis/necroptosis genes)
##   Table 3  Crosstalk results (paired comparisons, within-tissue coupling, interaction)
##   Table 4  Patient-level delta crosstalk (primary + sensitivity)
##   Table 5  Crosstalk modules (WGCNA) with gene counts and trait associations
##   Table 6  Hub genes (top 30) with evidence
##   Table 7  Enrichment top terms
##   Table 8  External validation results
##   Table 9  Field-effect analysis
## Outputs:
##   results/tables/Table1..Table9.tsv (long format, manuscript-ready)

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
tab_dir <- .dir_for("results_tables")

w <- readRDS(file.path(.dir_for("data_processed"), "wgcna_delta_result.rds"))

## ---- Table 1: cohort characteristics ----
mk_cohort_row <- function(gse, ana, n_pairs = NA) {
  ph <- ana$phenotype
  d <- data.frame(
    cohort = gse,
    platform = ana$platform,
    n_samples = nrow(ph),
    n_patients = length(unique(ph$patient_id)),
    n_paired = n_pairs,
    stringsAsFactors = FALSE
  )
  if ("er" %in% names(ph)) d$er_pos <- sum(ph$er == "pos", na.rm = TRUE)
  if ("her2" %in% names(ph)) d$her2_pos <- sum(ph$her2 == "pos", na.rm = TRUE)
  d
}
a1 <- readRDS(file.path(.dir_for("data_processed"), "GSE70947_analysis.rds"))
a2 <- readRDS(file.path(.dir_for("data_processed"), "GSE70905_analysis.rds"))
a3 <- readRDS(file.path(.dir_for("data_processed"), "GSE15852_analysis.rds"))
t1 <- rbind(
  mk_cohort_row("GSE70947", a1, 148),
  mk_cohort_row("GSE70905", a2, 47),
  mk_cohort_row("GSE15852", a3, 43)
)
.save_tsv(t1, "results_tables", "Table1_cohorts.tsv")

## ---- Table 2: DE summary ----
de_summ <- read.csv(file.path("results", "differential_expression", "pathway_de_summary.csv"), sep = "\t")
.save_tsv(de_summ, "results_tables", "Table2_de_summary.tsv")

## ---- Table 3: crosstalk ----
pc <- read.csv(file.path("results", "crosstalk", "pathway_score_paired_comparisons.csv"), sep = "\t")
wc <- read.csv(file.path("results", "crosstalk", "within_tissue_correlations.csv"), sep = "\t")
ic <- read.csv(file.path("results", "crosstalk", "coupling_interaction_models.csv"), sep = "\t")
.save_tsv(pc, "results_tables", "Table3a_paired_comparisons.tsv")
.save_tsv(wc, "results_tables", "Table3b_within_tissue_correlations.tsv")
.save_tsv(ic, "results_tables", "Table3c_interaction_models.tsv")

## ---- Table 4: delta crosstalk ----
dc <- read.csv(file.path("results", "crosstalk", "patient_delta_crosstalk.csv"), sep = "\t")
sens <- read.csv(file.path("results", "sensitivity", "sensitivity_summary.tsv"), sep = "\t")
.save_tsv(dc, "results_tables", "Table4a_delta_crosstalk.tsv")
.save_tsv(sens, "results_tables", "Table4b_sensitivity.tsv")

## ---- Table 5: crosstalk modules ----
mt <- read.csv(file.path("results", "wgcna", "module_trait_correlations.csv"), sep = "\t")
mod_summ <- read.csv(file.path("results", "hubs", "hub_summary.csv"), sep = "\t")
ct_mods <- c("ME3", "ME7", "ME9", "ME12", "ME17", "ME21", "ME22", "ME26", "ME28", "ME29", "ME30")
mt_ct <- mt[mt$module %in% ct_mods & mt$trait %in% c("DeltaFerroptosis_singscore", "DeltaNecroptosis_singscore"), ]
t5 <- merge(mod_summ, mt_ct[mt_ct$trait == "DeltaFerroptosis_singscore", c("module", "correlation", "p_value", "bh_fdr")],
            by = "module", all.x = TRUE)
names(t5)[names(t5) %in% c("correlation", "p_value", "bh_fdr")] <- c("r_ferro", "p_ferro", "fdr_ferro")
t5 <- merge(t5, mt_ct[mt_ct$trait == "DeltaNecroptosis_singscore", c("module", "correlation", "p_value", "bh_fdr")],
            by = "module", all.x = TRUE)
names(t5)[names(t5) %in% c("correlation", "p_value", "bh_fdr")] <- c("r_necro", "p_necro", "fdr_necro")
.save_tsv(t5, "results_tables", "Table5_crosstalk_modules.tsv")

## ---- Table 6: hub genes ----
hub <- read.csv(file.path("results", "hubs", "hub_gene_evidence.tsv"), sep = "\t")
hub_ct <- hub[hub$crosstalk_module, ]
hub_top <- head(hub_ct[order(hub_ct$hub_score, decreasing = TRUE), ], 30)
.save_tsv(hub_top[, c("gene_symbol", "module", "hub_score", "kWithin", "GS_ferro", "GS_necro",
                       "r_delta_ferro", "r_delta_necro", "neg_log10_fdr", "ferro_member", "necro_member")],
          "results_tables", "Table6_hub_genes.tsv")

## ---- Table 7: enrichment ----
for (src in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
  d <- read.csv(file.path("results", "enrichment", paste0(src, ".tsv")), sep = "\t")
  if (nrow(d) == 0) next
  d$p.adjust <- as.numeric(d$p.adjust)
  top <- head(d[order(d$p.adjust), ], 20)
  .save_tsv(top, "results_tables", paste0("Table7_", src, ".tsv"))
}

## ---- Table 8: validation ----
v1 <- read.csv(file.path("results", "validation", "validation_pathway_de.tsv"), sep = "\t")
v2 <- read.csv(file.path("results", "validation", "validation_coupling.tsv"), sep = "\t")
v3 <- read.csv(file.path("results", "validation", "validation_module_associations.tsv"), sep = "\t")
.save_tsv(v1, "results_tables", "Table8a_validation_pathway_de.tsv")
.save_tsv(v2, "results_tables", "Table8b_validation_coupling.tsv")
.save_tsv(v3, "results_tables", "Table8c_validation_modules.tsv")

## ---- Table 9: field effect ----
f1 <- read.csv(file.path("results", "field_effect", "field_effect_anova.tsv"), sep = "\t")
f2 <- read.csv(file.path("results", "field_effect", "group_pairwise.tsv"), sep = "\t")
f3 <- read.csv(file.path("results", "field_effect", "group_means.tsv"), sep = "\t")
.save_tsv(f1, "results_tables", "Table9a_field_effect_anova.tsv")
.save_tsv(f2, "results_tables", "Table9b_field_effect_pairwise.tsv")
.save_tsv(f3, "results_tables", "Table9c_field_effect_means.tsv")

.log("Tables assembled")
cat("DONE\n")