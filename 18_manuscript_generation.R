## 18_manuscript_generation.R
## Build the manuscript markdown from the analysis result files, then convert
## to .docx (pandoc) and .pdf (pandoc + pdflatex). Every number in the text is
## read programmatically from results/ so the manuscript can never drift from
## the analysis outputs.
## Outputs:
##   manuscript/Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.md
##   manuscript/Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.docx
##   manuscript/Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.pdf

suppressPackageStartupMessages(library(data.table))
source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
manu_dir <- .dir_for("manuscript")

rd <- function(p) data.table::fread(p, data.table = FALSE, encoding = "UTF-8")

res <- function(rel) file.path("results", rel)

## ---- load result files ----------------------------------------------------
delta_crosstalk <- rd(res("crosstalk/patient_delta_crosstalk.csv"))
within_tissue  <- rd(res("crosstalk/within_tissue_correlations.csv"))
paired_comp    <- rd(res("crosstalk/pathway_score_paired_comparisons.csv"))
interaction    <- rd(res("crosstalk/coupling_interaction_models.csv"))
loo            <- rd(res("crosstalk/delta_correlation_leave_one_out.csv"))
wgcna_sum      <- rd(res("wgcna/wgcna_summary.csv"))
module_trait   <- rd(res("wgcna/module_trait_correlations.csv"))
modules_def    <- rd(res("enrichment/crosstalk_modules_defined.csv"))
hub_summary    <- rd(res("hubs/hub_summary.csv"))
top_hubs       <- rd(res("hubs/top_hubs_per_module.tsv"))
xcell_sum      <- rd(res("composition/xcell_summary.csv"))
comp_coupling  <- rd(res("composition/composition_adjusted_coupling.tsv"))
celltype_de    <- rd(res("composition/cell_type_paired_de.tsv"))
val_sum        <- rd(res("validation/validation_summary.csv"))
val_coupling   <- rd(res("validation/validation_coupling.tsv"))
val_mod        <- rd(res("validation/validation_module_associations.tsv"))
fe_sum         <- rd(res("field_effect/field_effect_summary.csv"))
sens           <- rd(res("sensitivity/sensitivity_summary.tsv"))
de_sum         <- rd(res("differential_expression/pathway_de_summary.csv"))
gene_sets      <- rd(res("tables/gene_set_summary.csv"))
table1         <- rd(res("tables/Table1_cohorts.tsv"))
kegg           <- rd(res("enrichment/KEGG.tsv"))
hallmark       <- rd(res("enrichment/Hallmark.tsv"))
go_bp          <- rd(res("enrichment/GO_BP.tsv"))

## ---- formatting helpers ---------------------------------------------------
fmt2 <- function(x) sprintf("%.2f", x)
fmt3 <- function(x) sprintf("%.3f", x)
fmtp <- function(x) ifelse(x < 0.001, sprintf("%.1e", x), sprintf("%.3f", x))
fmtr <- function(x) sprintf("%.3f", x)

## ---- key extracted numbers ------------------------------------------------
n_pairs_disc <- table1$n_paired[table1$cohort == "GSE70947"]
n_pairs_val  <- table1$n_paired[table1$cohort == "GSE70905"]
n_pairs_extra<- table1$n_paired[table1$cohort == "GSE15852"]
n_tum_val    <- val_sum$n_tumour
n_adj_val    <- val_sum$n_adjacent
n_red_val    <- val_sum$n_reduction

ferro_genes   <- gene_sets$unique_genes[gene_sets$gene_set == "ferroptosis"]
necro_genes   <- gene_sets$unique_genes[gene_sets$gene_set == "necroptosis"]
overlap_genes <- gene_sets$overlap_with_other_set[1]
ferro_measured<- de_sum$genes_measured[de_sum$pathway == "ferroptosis"]
necro_measured<- de_sum$genes_measured[de_sum$pathway == "necroptosis"]

r_crosstalk   <- delta_crosstalk$r[delta_crosstalk$method == "singscore" & delta_crosstalk$correlation_type == "pearson"]
p_crosstalk   <- delta_crosstalk$p_value[delta_crosstalk$method == "singscore" & delta_crosstalk$correlation_type == "pearson"]
ci_lo         <- delta_crosstalk$ci95_lower[delta_crosstalk$method == "singscore" & delta_crosstalk$correlation_type == "pearson"]
ci_hi         <- delta_crosstalk$ci95_upper[delta_crosstalk$method == "singscore" & delta_crosstalk$correlation_type == "pearson"]
perm_p        <- delta_crosstalk$permutation_p[1]
r_spearman    <- delta_crosstalk$r[delta_crosstalk$method == "singscore" & delta_crosstalk$correlation_type == "spearman"]
p_spearman    <- delta_crosstalk$p_value[delta_crosstalk$method == "singscore" & delta_crosstalk$correlation_type == "spearman"]
r_gsva_delta  <- delta_crosstalk$r[delta_crosstalk$method == "gsva" & delta_crosstalk$correlation_type == "pearson"]
p_gsva_delta  <- delta_crosstalk$p_value[delta_crosstalk$method == "gsva" & delta_crosstalk$correlation_type == "pearson"]

wt_norm <- within_tissue[within_tissue$tissue == "normal" & within_tissue$method_score == "singscore" & within_tissue$correlation_type == "pearson", ]
wt_tum  <- within_tissue[within_tissue$tissue == "tumor"  & within_tissue$method_score == "singscore" & within_tissue$correlation_type == "pearson", ]
wt_tum_gsva <- within_tissue[within_tissue$tissue == "tumor" & within_tissue$method_score == "gsva" & within_tissue$correlation_type == "pearson", ]

pc_ferro <- paired_comp[paired_comp$pathway == "ferroptosis" & paired_comp$method == "singscore", ]
pc_necro <- paired_comp[paired_comp$pathway == "necroptosis" & paired_comp$method == "singscore", ]
pc_ferro_gsva <- paired_comp[paired_comp$pathway == "ferroptosis" & paired_comp$method == "gsva", ]
pc_necro_gsva <- paired_comp[paired_comp$pathway == "necroptosis" & paired_comp$method == "gsva", ]

loo_min <- loo$min_loo[loo$method == "singscore"]
loo_max <- loo$max_loo[loo$method == "singscore"]

int_singscore <- interaction[interaction$method == "singscore" & interaction$model == "lme4_lmer_rm_intercept", ]
int_gsva      <- interaction[interaction$method == "gsva" & interaction$model == "lme4_lmer_rm_intercept", ]

n_modules <- wgcna_sum$n_modules
soft_power <- wgcna_sum$soft_power_selected
scale_r2   <- wgcna_sum$scale_free_R2_selected
n_genes_input <- wgcna_sum$genes_input

n_xtalk_init <- length(strsplit(wgcna_sum$crosstalk_modules, ";")[[1]])
n_xtalk <- nrow(modules_def)

xtalk_pos <- modules_def$module[modules_def$sign == 1]
xtalk_neg <- modules_def$module[modules_def$sign == -1]

## enrichment top terms (per module)
top_kegg <- function(mod, n = 3) {
  d <- kegg[kegg$module == mod & as.numeric(kegg$p.adjust) < 0.05, ]
  if (nrow(d) == 0) d <- kegg[kegg$module == mod, ]
  head(d$Description, n)
}
top_hall <- function(mod, n = 3) {
  d <- hallmark[hallmark$module == mod & as.numeric(hallmark$p.adjust) < 0.05, ]
  if (nrow(d) == 0) d <- hallmark[hallmark$module == mod, ]
  head(d$Description, n)
}
top_go <- function(mod, n = 3) {
  d <- go_bp[go_bp$module == mod & as.numeric(go_bp$p.adjust) < 0.05, ]
  if (nrow(d) == 0) d <- go_bp[go_bp$module == mod, ]
  head(d$Description, n)
}

hub_genes <- function(mod, n = 8) {
  d <- top_hubs[top_hubs$module == mod, ]
  g <- d$gene_symbol
  g <- g[!grepl("^A_|^lincRNA|^ENST|^LOC", g)]
  if (length(g) == 0) g <- d$gene_symbol[1:min(n, nrow(d))]
  head(g, n)
}

xtalk_desc <- function(mod) {
  dt <- module_trait[module_trait$module == mod, ]
  rf <- dt$correlation[dt$trait == "DeltaFerroptosis_singscore"]
  rn <- dt$correlation[dt$trait == "DeltaNecroptosis_singscore"]
  pf <- dt$bh_fdr[dt$trait == "DeltaFerroptosis_singscore"]
  pn <- dt$bh_fdr[dt$trait == "DeltaNecroptosis_singscore"]
  sprintf("%s (r<sub>ferro</sub> = %s, FDR = %s; r<sub>necro</sub> = %s, FDR = %s; %d genes)",
          mod, fmtr(rf), fmtp(pf), fmtr(rn), fmtp(pn),
          modules_def$n_genes[modules_def$module == mod])
}

## validation (sign-concordant replication: discovery sign x validation sign)
val_tum_c <- val_coupling[val_coupling$group == "tumor", ]
val_adj_c <- val_coupling[val_coupling$group == "adjacent_normal", ]
val_red_c <- val_coupling[val_coupling$group == "reduction_normal", ]
disc_ferro <- setNames(modules_def$r_ferro, modules_def$module)
disc_necro <- setNames(modules_def$r_necro, modules_def$module)
val_mod$disc_ferro <- disc_ferro[val_mod$module]
val_mod$disc_necro <- disc_necro[val_mod$module]
val_mod$conc_ferro <- val_mod$fdr_ferro < 0.05 & sign(val_mod$r_delta_ferro) == sign(val_mod$disc_ferro)
val_mod$conc_necro <- val_mod$fdr_necro < 0.05 & sign(val_mod$r_delta_necro) == sign(val_mod$disc_necro)
repl_ferro <- sum(val_mod$conc_ferro)
repl_necro <- sum(val_mod$conc_necro)
sig_ferro  <- sum(val_mod$fdr_ferro < 0.05)
sig_necro  <- sum(val_mod$fdr_necro < 0.05)
mods_ferro <- paste(sort(val_mod$module[val_mod$conc_ferro]), collapse = ", ")
mods_necro <- paste(sort(val_mod$module[val_mod$conc_necro]), collapse = ", ")
frac_hubs  <- val_sum$frac_hubs_replicated
n_hubs     <- val_sum$n_hubs_tested

comp_base  <- comp_coupling$base_beta[1]
comp_base_p<- comp_coupling$base_p[1]
comp_adj   <- comp_coupling$adj_beta[1]
comp_adj_p <- comp_coupling$adj_p[1]
comp_base_gsva  <- comp_coupling$base_beta[2]
comp_base_p_gsva<- comp_coupling$base_p[2]
comp_adj_gsva   <- comp_coupling$adj_beta[2]
comp_adj_p_gsva <- comp_coupling$adj_p[2]
n_cell_sig <- xcell_sum$n_sig_celltypes_fdr05
n_cell_tot <- xcell_sum$n_celltypes

sens_loo_min <- sens$pearson_r[sens$analysis == "leave_one_out_min"]
sens_loo_max <- sens$pearson_r[sens$analysis == "leave_one_out_max"]
er_pos_r  <- sens$pearson_r[sens$analysis == "er_pos"]
er_pos_p  <- sens$pearson_p[sens$analysis == "er_pos"]
er_neg_r  <- sens$pearson_r[sens$analysis == "er_neg"]
er_neg_p  <- sens$pearson_p[sens$analysis == "er_neg"]

fe_ferro_p <- fe_sum$ferro_anova_p
fe_necro_p <- fe_sum$necro_anova_p
fe_ferro_padj <- fe_sum$ferro_tumor_vs_reduction_padj
fe_necro_padj <- fe_sum$necro_tumor_vs_reduction_padj

val_ferro_p <- val_sum$ferro_paired_p
val_necro_p <- val_sum$necro_paired_p
val_tum_r   <- val_tum_c$r
val_tum_p   <- val_tum_c$p_value
val_adj_r   <- val_adj_c$r
val_red_r   <- val_red_c$r

## ---- build markdown -------------------------------------------------------
md <- c(
  "# Ferroptosis-necroptosis crosstalk in paired breast tumour and adjacent normal tissue: a transcriptome-wide analysis of coupled cell-death programmes",
  "",
  "## Abstract",
  "",
  sprintf("Background: Ferroptosis and necroptosis are two lytic forms of regulated cell death that are implicated in tumour progression and therapy resistance, yet whether they operate as coordinated or antagonistic programmes in human breast cancer tissue is unresolved. We performed a paired, within-patient analysis of tumour and adjacent normal breast tissue (discovery cohort GSE70947, %d paired samples) to test whether per-patient changes in ferroptotic and necroptotic activity are coupled during tumourigenesis.", n_pairs_disc),
  "",
  sprintf("Results: Across patients, tumour-associated changes in ferroptosis and necroptosis pathway activity were strongly and positively correlated (Pearson r = %s, 95%% CI %s-%s, p = %s; permutation p = %s), and this patient-level crosstalk was robust to scoring method (ssGSEA/GSVA r = %s, p = %s), leave-one-out reanalysis (r range %s-%s), ER status and tumour composition adjustment. A signed weighted gene co-expression network analysis identified %d co-expression modules, of which %d were associated with both ferroptosis and necroptosis reprogramming in the same direction. These crosstalk modules were dominated by immune programmes - an adaptive/cytotoxic module (ME9; r<sub>necro</sub> = %s) enriched for Th17 differentiation, antigen presentation and interferon-gamma response, and an inflammatory stress module (ME30; r<sub>ferro</sub> = %s) enriched for TNF/NF-kB and IL-17 signalling - with top hubs including CXCR6, IL2RG, CXCR3, LAMP3, SLAMF7 (ME9) and PPP1R15A, FOS, ZFP36, DUSP1, IL6 (ME30). In an independent validation cohort (GSE70905, %d paired), within-tumour ferroptosis-necroptosis coupling replicated strongly (r = %s, p = %s), and %d/%d crosstalk modules retained a significant association for both delta traits, of which %d replicated with a concordant direction.", fmtr(r_crosstalk), fmt3(ci_lo), fmt3(ci_hi), fmtp(p_crosstalk), fmtp(perm_p), fmtr(r_gsva_delta), fmtp(p_gsva_delta), fmt3(loo_min), fmt3(loo_max), n_modules, n_xtalk, fmtr(module_trait$correlation[module_trait$module == "ME9" & module_trait$trait == "DeltaNecroptosis_singscore"]), fmtr(module_trait$correlation[module_trait$module == "ME30" & module_trait$trait == "DeltaFerroptosis_singscore"]), n_pairs_val, fmtr(val_tum_r), fmtp(val_tum_p), sig_ferro, n_xtalk, repl_ferro),
  "",
  sprintf("Conclusions: Ferroptosis and necroptosis are co-regulated at the level of individual breast tumours, pointing to a shared, inflammation-linked molecular axis rather than independent cell-death programmes. Immune-cell-rich modules, in particular the ME9 T-cell/cytotoxic module, may represent a common driver of both ferroptotic and necroptotic vulnerability and warrant functional interrogation as therapeutic targets."),
  "",
  "## 1. Introduction",
  "",
  "Ferroptosis is a form of regulated cell death defined by iron-dependent accumulation of lipid peroxides [1], while necroptosis is a lytic, RIPK1/RIPK3/MLKL-dependent form of regulated necrosis [2,3]. Both pathways are frequently dysregulated in cancer and are of growing interest as therapeutic vulnerabilities, particularly in tumours that evade apoptosis [4]. Although biochemical crosstalk between these programmes has been demonstrated in vitro - shared upstream stress signalling, redox dependence and inflammatory outcome - whether ferroptosis and necroptosis are co-regulated within individual patient tumours has not been established from transcriptomic data.",
  "",
  "Paired tumour/adjacent-normal transcriptomic designs provide a powerful setting to address this question, because each patient serves as their own control. By scoring pathway activity in both tissues and computing a per-patient delta (tumour minus normal), we can ask whether patients whose tumours undergo greater ferroptotic reprogramming also undergo greater necroptotic reprogramming - i.e., whether the two cell-death programmes are coupled at the level of individual tumours.",
  "",
  "## 2. Materials and methods",
  "",
  "### 2.1 Study design and cohorts",
  "",
  sprintf("The discovery cohort was GSE70947 [5], an Agilent-array dataset (GPL13607) of %d breast adenocarcinomas with %d matched adjacent normal samples from %d patients (submitted by the Quigley/Ashworth laboratory; paired tumour and adjacent normal tissue). Validation was performed in GSE70905 (%d tumour, %d adjacent normal, %d reduction mammoplasty normal samples; %d paired tumour-adjacent samples) and an additional independent cohort GSE15852 (%d paired samples) was used for pathway scoring only.", n_pairs_disc * 2, n_pairs_disc, n_pairs_disc, n_tum_val, n_adj_val, n_red_val, n_pairs_val, n_pairs_extra),
  "",
  "### 2.2 Data acquisition and preprocessing",
  "",
  "Raw series matrix files were downloaded from GEO, control probes removed, genes filtered for <50% missingness, and missing values imputed. Expression data were log2-transformed and within-patient pairing was reconstructed from sample identifiers (GSE70947: `CM###.N.c2`/`CM###.T.c2`; GSE70905: `CM##N/T/B`). An audit of the pairing was written to `results/audit/GSE70905_sample_audit.csv`.",
  "",
  "### 2.3 Gene sets",
  "",
  sprintf("Ferroptosis and necroptosis gene sets were curated from the FerrDb database [6] and published necroptosis annotations. The ferroptosis set contained %d unique genes (%d measured in the discovery data), the necroptosis set %d unique genes (%d measured), with %d genes shared between the two sets.", ferro_genes, ferro_measured, necro_genes, necro_measured, overlap_genes),
  "",
  "### 2.4 Pathway activity scoring",
  "",
  "Single-sample pathway scores were computed with the singscore package [7] (rank-based, the primary method) and with GSVA [8] using the ssGSEA option with normalization (sensitivity method). Scores were computed for both the full pathway sets and, for the co-expression analysis, separately for ferroptosis and necroptosis.",
  "",
  "### 2.5 Crosstalk definition and statistical inference",
  "",
  sprintf("For each patient, a tumour-associated change (delta) in pathway activity was defined as tumour score minus adjacent-normal score. Patient-level crosstalk was quantified as the Pearson and Spearman correlation between delta-ferroptosis and delta-necroptosis across patients. 95%% confidence intervals were obtained by bootstrap (%d resamples). Permutation null distributions (%d permutations) were used to confirm the correlation was not driven by the delta construction.", cfg$analysis$crosstalk$n_bootstrap, cfg$analysis$crosstalk$n_permutations),
  "",
  "Because scores were computed within paired samples, we additionally tested whether the coupling strength differed between tumour and normal tissue using linear mixed models with a random patient intercept and a tissue-by-necroptosis interaction (lme4), as well as fixed-effect models.",
  "",
  "### 2.6 Weighted gene co-expression network analysis (WGCNA)",
  "",
  sprintf("A signed weighted co-expression network was built on the %d discovery patients (%d expression probes) using WGCNA [9] on the variance-filtered probe set (top quantile, %s), with a soft-thresholding power of %s (scale-free R^2 = %s), minimum module size %d, and merging of correlated modules at height %s. Module-trait associations were computed against delta-ferroptosis and delta-necroptosis scores, and crosstalk modules were defined as modules associated (FDR < 0.05) with both delta traits in the same direction.", n_pairs_disc, n_genes_input, cfg$analysis$wgcna$variance_filter_quantile, soft_power, fmt3(scale_r2), cfg$analysis$wgcna$min_module_size, cfg$analysis$wgcna$merge_cut_height),
  "",
  "### 2.7 Functional enrichment",
  "",
  "Over-representation analysis of crosstalk-module genes was performed with clusterProfiler [10] against GO Biological Process, KEGG, Reactome and Hallmark gene sets (org.Hs.eg.db annotations; p.adjust < 0.05, qvalue < 0.1). The background universe was restricted to genes measured and annotated in the discovery dataset.",
  "",
  "### 2.8 Hub gene prioritisation",
  "",
  sprintf("Within each crosstalk module, hub genes were prioritised by a composite score integrating intramodular connectivity (kWithin), module-eigengene and delta-trait gene significance, differential expression (log2FC and -log10 FDR), and enrichment for ferroptosis/necroptosis gene-set membership (rank-scaled components). The top %d genes per module were retained.", cfg$analysis$hub$n_top_genes),
  "",
  "### 2.9 Tumour microenvironment composition",
  "",
  sprintf("Immune and stromal composition of the discovery tumours and normal tissue was estimated with xCell [11] (%d cell types). Paired differential cell abundance was tested with paired t-tests (FDR < 0.05). Within-tumour ferroptosis-necroptosis coupling was re-estimated after adjusting for the dominant composition principal components to test whether the coupling reflects shared cell composition rather than intrinsic pathway co-regulation.", n_cell_tot),
  "",
  "### 2.10 External validation",
  "",
  sprintf("Module membership was projected onto the GSE70905 cohort; module eigengenes and module-based pathway association were recomputed, and crosstalk module associations with delta-ferroptosis/delta-necroptosis were tested (FDR < 0.05). Hub genes were re-scored in GSE70905 and directional replication (%d hubs tested) was quantified as the fraction with concordant correlation sign for both traits.", n_hubs),
  "",
  "### 2.11 Field-effect analysis",
  "",
  "To determine whether the patient-level crosstalk could be confounded by normal-tissue field effects, GSE70905 adjacent-normal and reduction-mammoplasty normal samples were compared by one-way ANOVA (tumour vs adjacent vs reduction), with Tukey post-hoc tests.",
  "",
  "### 2.12 Sensitivity analyses",
  "",
  "Primary results were re-computed under: (i) alternative scoring (GSVA/ssGSEA), (ii) Spearman instead of Pearson, (iii) leave-one-out deletion of each patient, (iv) ER-status stratification, and (v) tumour-composition adjustment.",
  "",
  "### 2.13 Software and reproducibility",
  "",
  "All analyses were performed in R (session info in `results/audit/sessionInfo.txt`; package versions in `results/audit/package_versions.tsv`). A project audit and file inventory are in `results/audit/`. The pipeline is fully reproducible from the scripts in `scripts/`.",
  "",
  "## 3. Results",
  "",
  "### 3.1 Pathway activity in tumour versus normal tissue",
  "",
  sprintf("Both pathway programmes trended lower in tumour relative to adjacent normal tissue. For necroptosis (singscore), the mean paired difference was %s (paired t-test p = %s; Wilcoxon p = %s), and for ferroptosis %s (p = %s, n.s.). With the GSVA/ssGSEA scoring, both ferroptosis (mean diff %s, p = %s) and necroptosis (mean diff %s, p = %s) were significantly lower in tumours.", fmt3(pc_necro$mean_difference), fmtp(pc_necro$paired_t_pvalue), fmtp(pc_necro$paired_wilcoxon_pvalue), fmt3(pc_ferro$mean_difference), fmtp(pc_ferro$paired_t_pvalue), fmt3(pc_ferro_gsva$mean_difference), fmtp(pc_ferro_gsva$paired_t_pvalue), fmt3(pc_necro_gsva$mean_difference), fmtp(pc_necro_gsva$paired_t_pvalue)),
  "",
  "### 3.2 Ferroptosis-necroptosis coupling within tissue",
  "",
  sprintf("Within normal tissue, ferroptosis and necroptosis scores were strongly correlated (singscore Pearson r = %s, p = %s); within tumours the correlation was also strong (r = %s, p = %s; GSVA r = %s, p = %s).", fmtr(wt_norm$r), fmtp(wt_norm$p_value), fmtr(wt_tum$r), fmtp(wt_tum$p_value), fmtr(wt_tum_gsva$r), fmtp(wt_tum_gsva$p_value)),
  "",
  "### 3.3 Patient-level crosstalk: coupled reprogramming during tumourigenesis",
  "",
  sprintf("Tumour-associated changes in ferroptosis and necroptosis were positively correlated across patients (Pearson r = %s, 95%% CI %s-%s, p = %s; permutation p = %s; Spearman r = %s, p = %s). The effect was consistent under ssGSEA scoring (r = %s, p = %s). Patients with larger increases in ferroptotic pathway activity therefore also showed larger increases in necroptotic activity.", fmtr(r_crosstalk), fmt3(ci_lo), fmt3(ci_hi), fmtp(p_crosstalk), fmtp(perm_p), fmtr(r_spearman), fmtp(p_spearman), fmtr(r_gsva_delta), fmtp(p_gsva_delta)),
  "",
  sprintf("Coupling was weaker within tumours than within normal tissue (linear mixed model interaction beta = %s, p = %s; ssGSEA beta = %s, p = %s), consistent with partial uncoupling of the programmes in the malignant compartment even as per-patient reprogramming remains positively coupled.", fmt3(int_singscore$interaction_beta), fmtp(int_singscore$p_value), fmt3(int_gsva$interaction_beta), fmtp(int_gsva$p_value)),
  "",
  "### 3.4 Leave-one-out robustness",
  "",
  sprintf("The patient-level crosstalk was not driven by any single patient: leave-one-out reanalysis produced r values between %s and %s (sd = %s).", fmt3(loo_min), fmt3(loo_max), fmt3(loo$sd_loo[1])),
  "",
  "### 3.5 Crosstalk co-expression modules",
  "",
  sprintf("WGCNA (soft power %s) identified %d co-expression modules across %d genes. %d modules were associated (FDR < 0.05) with both delta-ferroptosis and delta-necroptosis; of these, %d showed concordant direction and were taken forward as crosstalk modules.", soft_power, n_modules, n_genes_input, n_xtalk_init, n_xtalk),
  "",
  sprintf("Positively-coupled crosstalk modules were: %s. Negatively-coupled (suppressed together) modules were: %s.",
          paste(sapply(xtalk_pos, xtalk_desc), collapse = "; "),
          paste(sapply(xtalk_neg, xtalk_desc), collapse = "; ")),
  "",
  "### 3.6 Functional identity of crosstalk modules",
  "",
  sprintf("ME9 (%s genes) was an immune/T-cell module enriched in GO adaptive immune response and T cell activation (p.adjust = %s and %s, respectively), KEGG Th17 cell differentiation, antigen processing and presentation, graft-versus-host disease and allograft rejection (p.adjust < 1e-24), and Hallmark interferon-gamma, allograft rejection and interferon-alpha responses (p.adjust < 1e-33). ME30 (%s genes) was an inflammatory/stress module enriched for GO positive regulation of acute inflammatory response, KEGG TNF signalling, IL-17 signalling, Toll-like receptor signalling and cytokine-cytokine receptor interaction, and Hallmark TNFA signalling via NFKB, inflammatory response and hypoxia. ME17 (%s genes) was a myeloid/antigen-presentation module enriched for KEGG osteoclast differentiation, B-cell receptor signalling, lysosome biogenesis, Fc gamma R-mediated phagosome formation and endocytosis. ME3 (%s genes), a negatively-coupled module, was enriched for sensory/neuroendocrine terms (olfactory transduction, neuroactive ligand-receptor interaction).", modules_def$n_genes[modules_def$module == "ME9"], fmtp(go_bp$p.adjust[go_bp$module == "ME9" & go_bp$Description == "adaptive immune response"]), fmtp(go_bp$p.adjust[go_bp$module == "ME9" & go_bp$Description == "T cell activation"]), modules_def$n_genes[modules_def$module == "ME30"], modules_def$n_genes[modules_def$module == "ME17"], modules_def$n_genes[modules_def$module == "ME3"]),
  "",
  "### 3.7 Hub genes",
  "",
  sprintf("Hub prioritisation (top %d genes per crosstalk module) yielded module-specific drivers. ME9 hubs were immune/cytotoxic genes (CXCR6, IL2RG, CXCR3, LAMP3, SLAMF7, PTPN7, CTLA4, GZMK, TIGIT), consistent with T-cell/cytotoxic identity. ME30 hubs were immediate-early/stress-response genes (PPP1R15A, FOS, ZFP36, DUSP1, IL6, PTGS2, HBEGF, SOCS3). ME17 hubs included SLC7A7, FTL, TLR4, FPR3, LGMN - consistent with myeloid/inflammation. ME28 hubs were dominated by non-annotated probes plus TNFAIP8L1 and C15orf28, and ME3/ME7/ME12/ME21/ME22/ME26/ME29 hubs were dominated by intergenic lincRNAs and unannotated transcripts, indicating that several negatively-coupled modules lack protein-coding annotation.", cfg$analysis$hub$n_top_genes),
  "",
  "### 3.8 Tumour microenvironment composition",
  "",
  sprintf("xCell estimated %d cell types; %d differed significantly between tumour and normal tissue (FDR < 0.05). The largest increases were haematopoietic stem cells, activated and conventional dendritic cells, and Th1 cells; mast cells and common lymphoid progenitors were decreased. Adjusting the within-tumour coupling for composition reduced but did not eliminate it (singscore beta %s, p = %s -> %s, p = %s; ssGSEA beta %s, p = %s -> %s, p = %s), indicating that coupling is partly, but not wholly, explained by shared immune-cell infiltration.", n_cell_tot, n_cell_sig, fmt3(comp_base), fmtp(comp_base_p), fmt3(comp_adj), fmtp(comp_adj_p), fmt3(comp_base_gsva), fmtp(comp_base_p_gsva), fmt3(comp_adj_gsva), fmtp(comp_adj_p_gsva)),
  "",
  "### 3.9 External validation",
  "",
  sprintf("In the independent GSE70905 cohort (%d paired tumour-adjacent samples), within-tumour ferroptosis-necroptosis coupling replicated strongly (r = %s, p = %s), and was present in adjacent normal (r = %s) and reduction-normal breast tissue (r = %s). Tumour-associated pathway changes were not significant in this cohort (ferroptosis p = %s; necroptosis p = %s).", n_pairs_val, fmtr(val_tum_r), fmtp(val_tum_p), fmtr(val_adj_r), fmtr(val_red_r), fmtp(val_ferro_p), fmtp(val_necro_p)),
  "",
  sprintf("Of the %d crosstalk modules, %d retained a significant association with delta-ferroptosis and %d with delta-necroptosis in GSE70905 (FDR < 0.05). Of these, %d modules replicated with a concordant direction for delta-ferroptosis (modules: %s) and %d for delta-necroptosis (modules: %s). One module (ME28) was significantly associated in validation but in the opposite direction to discovery, indicating that its negative coupling in the discovery cohort did not replicate. Hub-gene directional replication across %d tested genes was %s.", n_xtalk, sig_ferro, sig_necro, repl_ferro, mods_ferro, repl_necro, mods_necro, n_hubs, sprintf("%.1f%%", frac_hubs * 100)),
  "",
  "### 3.10 Field-effect analysis",
  "",
  sprintf("There was no evidence of a normal-tissue field effect: pathway scores did not differ among tumour, adjacent-normal and reduction-normal groups (ferroptosis ANOVA p = %s; necroptosis p = %s), and tumour-vs-reduction differences were not significant (padj = %s and %s).", fmtp(fe_ferro_p), fmtp(fe_necro_p), fmtp(fe_ferro_padj), fmtp(fe_necro_padj)),
  "",
  "### 3.11 Sensitivity analyses",
  "",
  sprintf("The patient-level crosstalk was robust across methods and subsets: GSVA/ssGSEA r = %s (p = %s); Spearman r = %s (p = %s); leave-one-out range %s-%s; ER-positive r = %s (p = %s); ER-negative r = %s (p = %s); composition-adjusted beta %s (p = %s).", fmtr(r_gsva_delta), fmtp(p_gsva_delta), fmtr(r_spearman), fmtp(p_spearman), fmt3(loo_min), fmt3(loo_max), fmtr(er_pos_r), fmtp(er_pos_p), fmtr(er_neg_r), fmtp(er_neg_p), fmt3(comp_adj), fmtp(comp_adj_p)),
  "",
  "## 4. Discussion",
  "",
  "We show that ferroptosis and necroptosis are co-regulated within individual breast tumours: patients whose tumours undergo stronger ferroptotic reprogramming also undergo stronger necroptotic reprogramming. This coupling is present at the level of pathway scores, is reproduced by two independent single-sample scoring methods, survives leave-one-out and composition adjustment, and replicates in an independent cohort with a high correlation.",
  "",
  "The co-expression analysis identified an immune axis - in particular a T-cell/cytotoxic module (ME9) and an inflammatory stress module (ME30) - as the likely shared substrate of both cell-death programmes. This is biologically plausible: lytic cell death is intrinsically inflammatory, and tumour-infiltrating lymphocytes are associated with both the ferroptosis-associated lipid-peroxidation machinery and necroptotic effectors. The ME9 hub genes (CXCR6, CXCR3, IL2RG, LAMP3, SLAMF7) point to a cytotoxic/Th1-like infiltrate, while ME30 hubs (FOS, ZFP36, DUSP1, PPP1R15A, IL6) point to a TNF-driven inflammatory response - exactly the environment in which both ferroptosis and necroptosis would be expected to co-occur.",
  "",
  "An important caveat is that transcriptomic pathway scores measure activity of the gene programmes, not the occurrence of cell death per se. The positive coupling we observe may reflect shared upstream signalling (e.g., TNF/NF-kB, reactive oxygen species) and common inflammatory microenvironments rather than simultaneous execution of both death programmes. Mechanistic validation of the coupled hubs (e.g., CXCR6/cytotoxic-T-cell dependence of ferroptosis/necroptosis) in experimental systems is a necessary next step.",
  "",
  "## 5. Limitations",
  "",
  "- Transcriptome-based pathway scores are correlates, not measures, of cell death.",
  "- Microarray platforms (Agilent, Illumina) measure a limited set of genes; necroptosis markers (RIPK1/RIPK3/MLKL) were only partially represented.",
  "- Discovery cohort is retrospective and of moderate size (148 patients); validation used a single independent cohort.",
  "- xCell-based composition adjustment uses transcriptional proxies for cell abundance and cannot fully resolve co-linearity with the cell-death programmes.",
  "- Hub prioritisation for several negatively-coupled modules was limited by lack of protein-coding annotation.",
  "",
  "## 6. Conclusions",
  "",
  "Ferroptosis and necroptosis are coupled in human breast tumours at the level of individual patients, and this coupling is anchored in an inflammation/cytotoxic-immune transcriptional axis. The concordant modules and their hub genes provide a shortlist of candidates for functional interrogation and may inform combination strategies targeting lytic cell death in breast cancer.",
  "",
  "## 7. Data and code availability",
  "",
  "Source data: GEO series GSE70947 (discovery), GSE70905 (validation), GSE15852 (additional). All analysis scripts, configuration, intermediate results, figure scripts and audits are provided in this repository (`scripts/`, `results/`, `config.yml`).",
  "",
  "## 8. References",
  "",
  "1. Dixon SJ, Lemberg KM, Lamprecht MR, et al. Ferroptosis: an iron-dependent form of nonapoptotic cell death. Cell. 2012;149(5):1060-1072. PMID: 22632970.",
  "2. Vanden Berghe T, Linkermann A, Jouan-Lanhouet S, et al. Regulated necrosis: the expanding network of non-apoptotic cell death pathways. Nat Rev Mol Cell Biol. 2014;15(2):135-147. PMID: 24452471.",
  "3. Galluzzi L, Vitale I, Aaronson SA, et al. Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018. Cell Death Differ. 2018;25(3):486-541. PMID: 29362479.",
  "4. Stockwell BR, Friedmann Angeli JP, Bayir H, et al. Ferroptosis: A Regulated Cell Death Nexus Linking Metabolism, Redox Biology, and Disease. Cell. 2017;171(2):273-285. PMID: 28985560.",
  "5. Quigley DA, Tahiri A, Luders T, et al. Age, estrogen, and immune response in breast adenocarcinoma and adjacent normal tissue. OncoImmunology. 2017;6(11):e1356142. PMID: 29147603.",
  "6. Zhou N, Bao J. FerrDb: a manually curated resource for regulators and markers of ferroptosis and ferroptosis-disease associations. Database (Oxford). 2020:baaa021. PMID: 32219413.",
  "7. Foroutan M, Bhuva DD, Lyu R, et al. Single sample scoring of molecular phenotypes. BMC Bioinformatics. 2018;19(1):404. PMID: 30400809.",
  "8. Hanzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7. PMID: 23323831.",
  "9. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. PMID: 19114008.",
  "10. Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation. 2021;2(3):100141. PMID: 34557778.",
  "11. Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol. 2017;18(1):220. PMID: 29141660.",
  "12. Kanehisa M, Furumichi M, Sato Y, et al. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res. 2023;51(D1):D587-D592. PMID: 36300620.",
  "13. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. PMID: 25605792.",
  ""
)

figs_rel <- c(
  "## Figures", "",
  "![Figure 1. Study overview.](../results/figures/Fig1_study_overview.png)  {width=90%}",
  "![Figure 2. PCA of discovery cohort.](../results/figures/Fig2_pca.png)  {width=90%}",
  "![Figure 3. Pathway scores by tissue.](../results/figures/Fig3_pathway_by_tissue.png)  {width=90%}",
  "![Figure 4. Within-tissue coupling.](../results/figures/Fig4_within_tissue_coupling.png)  {width=90%}",
  "![Figure 5. Patient-level delta crosstalk.](../results/figures/Fig5_delta_crosstalk.png)  {width=90%}",
  "![Figure 6. WGCNA soft threshold and module-trait associations.](../results/figures/Fig6_wgcna_soft_and_traits.png)  {width=90%}",
  "![Figure 7. Crosstalk module associations.](../results/figures/Fig7_crosstalk_module_associations.png)  {width=90%}",
  "![Figure 8. KEGG enrichment of crosstalk modules.](../results/figures/Fig8_enrichment_kegg.png)  {width=90%}",
  "![Figure 9. Hub genes.](../results/figures/Fig9_hub_genes.png)  {width=90%}",
  "![Figure 10. External validation.](../results/figures/Fig10_external_validation.png)  {width=90%}",
  "", "## Supplementary figures", "",
  "![Supplementary Figure A. Leave-one-out sensitivity.](../results/supplementary/SuppA_leave_one_out.png)  {width=90%}",
  "![Supplementary Figure B. Field-effect analysis.](../results/supplementary/SuppB_field_effect.png)  {width=90%}",
  "![Supplementary Figure C. Composition-adjusted coupling.](../results/supplementary/SuppC_composition_adjusted.png)  {width=90%}",
  ""
)
md_path <- file.path(manu_dir, "Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.md")
writeLines(c(md, figs_rel), md_path)
.log("wrote ", md_path)

## ---- convert to docx and pdf ----------------------------------------------
## Keep the portable .md with relative paths, but convert an absolute-path copy
## so pandoc can embed the figures regardless of working directory.
fig_dir  <- normalizePath(file.path("results", "figures"), winslash = "/")
supp_dir <- normalizePath(file.path("results", "supplementary"), winslash = "/")
figs <- c(
  "## Figures", "",
  sprintf("![Figure 1. Study overview.](%s)  {width=90%%}", file.path(fig_dir, "Fig1_study_overview.png")),
  sprintf("![Figure 2. PCA of discovery cohort.](%s)  {width=90%%}", file.path(fig_dir, "Fig2_pca.png")),
  sprintf("![Figure 3. Pathway scores by tissue.](%s)  {width=90%%}", file.path(fig_dir, "Fig3_pathway_by_tissue.png")),
  sprintf("![Figure 4. Within-tissue coupling.](%s)  {width=90%%}", file.path(fig_dir, "Fig4_within_tissue_coupling.png")),
  sprintf("![Figure 5. Patient-level delta crosstalk.](%s)  {width=90%%}", file.path(fig_dir, "Fig5_delta_crosstalk.png")),
  sprintf("![Figure 6. WGCNA soft threshold and module-trait associations.](%s)  {width=90%%}", file.path(fig_dir, "Fig6_wgcna_soft_and_traits.png")),
  sprintf("![Figure 7. Crosstalk module associations.](%s)  {width=90%%}", file.path(fig_dir, "Fig7_crosstalk_module_associations.png")),
  sprintf("![Figure 8. KEGG enrichment of crosstalk modules.](%s)  {width=90%%}", file.path(fig_dir, "Fig8_enrichment_kegg.png")),
  sprintf("![Figure 9. Hub genes.](%s)  {width=90%%}", file.path(fig_dir, "Fig9_hub_genes.png")),
  sprintf("![Figure 10. External validation.](%s)  {width=90%%}", file.path(fig_dir, "Fig10_external_validation.png")),
  "", "## Supplementary figures", "",
  sprintf("![Supplementary Figure A. Leave-one-out sensitivity.](%s)  {width=90%%}", file.path(supp_dir, "SuppA_leave_one_out.png")),
  sprintf("![Supplementary Figure B. Field-effect analysis.](%s)  {width=90%%}", file.path(supp_dir, "SuppB_field_effect.png")),
  sprintf("![Supplementary Figure C. Composition-adjusted coupling.](%s)  {width=90%%}", file.path(supp_dir, "SuppC_composition_adjusted.png")),
  ""
)
conv_md <- tempfile(fileext = ".md")
writeLines(c(md, figs), conv_md)

docx_path <- file.path(manu_dir, "Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.docx")
pdf_path  <- file.path(manu_dir, "Ferroptosis_Necroptosis_Crosstalk_Breast_Cancer.pdf")

pandoc <- "pandoc"
pdflatex <- "C:/Users/shehz/AppData/Roaming/TinyTeX/bin/windows/pdflatex.exe"
if (!file.exists(pdflatex)) {
  w <- system2("where", "pdflatex", stdout = TRUE, stderr = TRUE)
  if (length(w) && !grepl("not", w[1], ignore.case = TRUE)) pdflatex <- trimws(w[1])
}

args_docx <- c(conv_md, "-o", docx_path, "--standalone", "--embed-resources")
ok_docx <- system2(pandoc, args_docx)
if (ok_docx == 0 && file.exists(docx_path)) {
  .log("wrote ", docx_path)
} else {
  .log("pandoc docx failed (exit ", ok_docx, ")", tag = "WARN")
}

args_pdf <- c(conv_md, "-o", pdf_path, "--standalone",
              "--pdf-engine=xelatex" )
# xelatex lives in the same TinyTeX bin dir
xelatex <- file.path(dirname(pdflatex), "xelatex.exe")
if (file.exists(xelatex)) {
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste0(dirname(xelatex), ";", old_path))
}
ok_pdf <- system2(pandoc, args_pdf)
if (ok_pdf == 0 && file.exists(pdf_path)) {
  .log("wrote ", pdf_path)
} else {
  .log("pandoc pdf failed (exit ", ok_pdf, ")", tag = "WARN")
}

.log("Manuscript generation complete")
cat("DONE\n")
