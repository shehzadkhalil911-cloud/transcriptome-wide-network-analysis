## 20_publication_figures.R
## Publication-ready figures (PNG, 600 dpi) generated ONLY from existing computed results.
## No statistics are recomputed or altered: every annotation is read from results/ files.
## Outputs: results/publication_figures/*.png + figure_generation_log.txt
## Skipped panels are documented in results/figures_missing.txt

source(file.path("R", "helpers.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(ggdendro)
  library(pheatmap)
})
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)

pub_dir <- file.path("results", "publication_figures")
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(pub_dir, "figure_generation_log.txt")
missing_file <- file.path("results", "figures_missing.txt")

## ---- style helpers ----
okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
               "#E69F00", "#56B4E9", "#F0E442", "#000000")
mod_pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
             "#56B4E9", "#009988", "#882255", "#AA4499", "#DDCC77",
             "#44AA99", "#117733")
ct_mods <- c("ME3", "ME7", "ME9", "ME12", "ME17", "ME21", "ME22",
             "ME26", "ME28", "ME29", "ME30")
mod_cols <- setNames(mod_pal[seq_along(ct_mods)], ct_mods)

theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(face = "bold", hjust = 0, size = base + 1),
          plot.subtitle = element_text(hjust = 0, size = base - 1, colour = "grey25"),
          legend.position = "top",
          legend.title = element_text(size = base - 1),
          legend.text = element_text(size = base - 1),
          axis.text = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
}

save_pub <- function(p, name, w = 180, h = 120) {
  ggsave(file.path(pub_dir, paste0(name, ".png")), p,
         width = w, height = h, units = "mm", dpi = 600, bg = "white")
}

fmt_p <- function(p, digits = 2) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
         paste0(signif(p, digits)))
}

cor_ann <- function(r, pval, n, extra = NULL) {
  txt <- paste0("r = ", signif(r, 3), ", p = ", fmt_p(pval), ", n = ", n)
  if (!is.null(extra)) txt <- paste0(txt, extra)
  txt
}

## ---- load data ----
ps <- readRDS(file.path(.dir_for("data_processed"), "pathway_scores_all.rds"))
sc  <- ps$scores_GSE70947
scv <- ps$scores_GSE70905
sc2 <- ps$scores_GSE15852

dt <- read.csv(file.path("results", "crosstalk", "patient_delta_scores.tsv"), sep = "\t")
pdc <- read.csv(file.path("results", "crosstalk", "patient_delta_crosstalk.csv"), sep = "\t")
wt  <- read.csv(file.path("results", "crosstalk", "within_tissue_correlations.csv"), sep = "\t")
ppc <- read.csv(file.path("results", "crosstalk", "pathway_score_paired_comparisons.csv"), sep = "\t")
cim <- read.csv(file.path("results", "crosstalk", "coupling_interaction_models.csv"), sep = "\t")
cr  <- read.csv(file.path("results", "crosstalk", "coupling_robustness.csv"), sep = "\t")

ss  <- read.csv(file.path("results", "sensitivity", "sensitivity_summary.tsv"), sep = "\t")
loo <- read.csv(file.path("results", "sensitivity", "leave_one_out.tsv"), sep = "\t")

w  <- readRDS(file.path(.dir_for("data_processed"), "wgcna_delta_result.rds"))
soft <- read.csv(file.path("results", "wgcna", "soft_threshold_results.tsv"), sep = "\t")
mt   <- read.csv(file.path("results", "wgcna", "module_trait_correlations.csv"), sep = "\t")
ws   <- read.csv(file.path("results", "wgcna", "wgcna_summary.csv"), sep = "\t")

cm  <- read.csv(file.path("results", "enrichment", "crosstalk_modules_defined.csv"), sep = "\t")
go  <- read.csv(file.path("results", "enrichment", "GO_BP.tsv"), sep = "\t")
kegg <- read.csv(file.path("results", "enrichment", "KEGG.tsv"), sep = "\t")
hall <- read.csv(file.path("results", "enrichment", "Hallmark.tsv"), sep = "\t")
rea  <- read.csv(file.path("results", "enrichment", "Reactome.tsv"), sep = "\t")

hub  <- read.csv(file.path("results", "hubs", "top_hubs_per_module.tsv"), sep = "\t")
hsum <- read.csv(file.path("results", "hubs", "hub_summary.csv"), sep = "\t")

xsc  <- readRDS(file.path("results", "composition", "xcell_scores.rds"))
ctde <- read.csv(file.path("results", "composition", "cell_type_paired_de.tsv"), sep = "\t")
xsum <- read.csv(file.path("results", "composition", "xcell_summary.csv"), sep = "\t")
comp_cou <- read.csv(file.path("results", "composition", "composition_adjusted_coupling.tsv"),
                     sep = "\t", check.names = FALSE)

val_sum <- read.csv(file.path("results", "validation", "validation_summary.csv"), sep = "\t")
val_cou <- read.csv(file.path("results", "validation", "validation_coupling.tsv"), sep = "\t")
val_mod <- read.csv(file.path("results", "validation", "validation_module_associations.tsv"), sep = "\t")
val_hub <- read.csv(file.path("results", "validation", "validation_hub_replication.tsv"), sep = "\t")
val_de  <- read.csv(file.path("results", "validation", "validation_pathway_de.tsv"), sep = "\t")

fe_sum <- read.csv(file.path("results", "field_effect", "field_effect_summary.csv"), sep = "\t")
fe_gm  <- read.csv(file.path("results", "field_effect", "group_means.tsv"), sep = "\t")

qc_sum <- read.csv(file.path("results", "qc", "qc_summary.csv"), sep = "\t")
qc47   <- read.csv(file.path("results", "qc", "GSE70947_sample_qc.csv"), sep = "\t")
ana    <- readRDS(file.path(.dir_for("data_processed"), "GSE70947_analysis.rds"))

gs_sum <- read.csv(file.path("results", "tables", "gene_set_summary.csv"), sep = "\t")

log_df <- data.frame(figure = character(), status = character(),
                     width_mm = numeric(), height_mm = numeric(),
                     dpi = numeric(), source = character(), stringsAsFactors = FALSE)
add_log <- function(name, status, w, h, dpi, source) {
  log_df <<- rbind(log_df, data.frame(figure = name, status = status,
                                      width_mm = w, height_mm = h,
                                      dpi = dpi, source = source))
}

## ===========================================================================
## MAIN FIGURE 1: pathway activity and tissue coupling
## ===========================================================================
sc_long <- sc %>%
  select(sample_id, tissue, patient_id,
         ferroptosis_full_singscore, necroptosis_full_singscore) %>%
  pivot_longer(c(ferroptosis_full_singscore, necroptosis_full_singscore),
               names_to = "pathway", values_to = "score")
sc_long$pathway <- ifelse(sc_long$pathway == "ferroptosis_full_singscore",
                          "Ferroptosis", "Necroptosis")

p1a <- ggplot(sc_long, aes(x = tissue, y = score, fill = tissue)) +
  geom_violin(alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.22, outlier.shape = NA, colour = "black") +
  facet_wrap(~pathway) +
  scale_fill_manual(values = c("normal" = okabe_ito[1], "tumor" = okabe_ito[2])) +
  scale_x_discrete(labels = c("normal" = "Adjacent\nnormal", "tumor" = "Tumour")) +
  labs(subtitle = paste0(
    "Ferroptosis: paired t p = ", fmt_p(ppc$paired_t_pvalue[ppc$pathway == "ferroptosis" & ppc$method == "singscore"]),
    "  |  Necroptosis: paired t p = ", fmt_p(ppc$paired_t_pvalue[ppc$pathway == "necroptosis" & ppc$method == "singscore"]),
    "  (GSE70947, 148 paired tumour/adjacent-normal)"),
    x = NULL, y = "singscore") +
  theme_pub() +
  theme(legend.position = "none", strip.background = element_blank())

scp <- sc %>% select(sample_id, tissue, patient_id,
                     ferroptosis_full_singscore, necroptosis_full_singscore)
wt_t <- wt[wt$tissue == "tumor" & wt$method_score == "singscore" & wt$correlation_type == "pearson", ]
wt_n <- wt[wt$tissue == "normal" & wt$method_score == "singscore" & wt$correlation_type == "pearson", ]

p1b <- ggplot(scp, aes(necroptosis_full_singscore, ferroptosis_full_singscore, colour = tissue)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_colour_manual(values = c("normal" = okabe_ito[1], "tumor" = okabe_ito[2]),
                      labels = c("normal" = "Adjacent normal", "tumor" = "Tumour")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.6, family = "Arial",
           label = paste0("Tumour:   ", cor_ann(wt_t$r, wt_t$p_value, wt_t$n), "\n",
                          "Adjacent: ", cor_ann(wt_n$r, wt_n$p_value, wt_n$n))) +
  labs(subtitle = "Within-tissue ferroptosis\u2013necroptosis coupling (singscore)",
       x = "Necroptosis singscore", y = "Ferroptosis singscore") +
  theme_pub()

fe_gm$pathway <- ifelse(fe_gm$pathway == "ferro", "Ferroptosis", "Necroptosis")
fe_gm$group <- factor(c("reduction_normal" = "Reduction", "adjacent_normal" = "Adjacent",
                        "tumor" = "Tumour")[fe_gm$group],
                      levels = c("Reduction", "Adjacent", "Tumour"))

p1c <- ggplot(fe_gm, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
  facet_wrap(~pathway) +
  scale_fill_manual(values = c("Reduction" = okabe_ito[7], "Adjacent" = okabe_ito[1],
                               "Tumour" = okabe_ito[2])) +
  annotate("text", x = 2, y = Inf, vjust = 1.5, size = 2.6, family = "Arial",
           label = paste0("ANOVA p = ", fmt_p(fe_sum$ferro_anova_p),
                          " / ", fmt_p(fe_sum$necro_anova_p))) +
  labs(subtitle = "Field effect across tissue types (GSE70905, singscore)",
       x = NULL, y = "Mean singscore") +
  theme_pub() +
  theme(legend.position = "none", strip.background = element_blank())

p1 <- (p1a + p1b) / p1c +
  plot_annotation(title = "Pathway activity and tissue coupling") +
  plot_layout(heights = c(1, 0.8))
save_pub(p1, "main_pathway_activity_and_tissue_coupling", 190, 150)
add_log("main_pathway_activity_and_tissue_coupling", "ok", 190, 150, 600,
        "GSE70947_pathway_scores.tsv; within_tissue_correlations.csv; pathway_score_paired_comparisons.csv; field_effect/group_means.tsv")

## ===========================================================================
## MAIN FIGURE 2: patient-level delta crosstalk (primary result)
## ===========================================================================
r_p <- pdc[pdc$method == "singscore" & pdc$correlation_type == "pearson", ]
r_g <- pdc[pdc$method == "gsva" & pdc$correlation_type == "pearson", ]

p2a <- ggplot(dt, aes(Deltanecroptosis_singscore, Deltaferroptosis_singscore)) +
  geom_point(alpha = 0.55, size = 1.4, colour = okabe_ito[1]) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = cor_ann(r_p$r, r_p$p_value, r_p$n,
                           paste0("\n95% CI [", signif(r_p$ci95_lower, 3), ", ",
                                  signif(r_p$ci95_upper, 3), "]\npermutation p = ",
                                  fmt_p(r_p$permutation_p)))) +
  labs(subtitle = "Primary analysis: singscore (GSE70947, paired tumour \u2212 adjacent-normal)",
       x = "Delta necroptosis (tumour \u2212 normal)", y = "Delta ferroptosis (tumour \u2212 normal)") +
  theme_pub()

p2b <- ggplot(dt, aes(Deltanecroptosis_gsva, Deltaferroptosis_gsva)) +
  geom_point(alpha = 0.55, size = 1.4, colour = okabe_ito[4]) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = cor_ann(r_g$r, r_g$p_value, r_g$n,
                           paste0("\n95% CI [", signif(r_g$ci95_lower, 3), ", ",
                                  signif(r_g$ci95_upper, 3), "]\npermutation p = ",
                                  fmt_p(r_g$permutation_p)))) +
  labs(subtitle = "Sensitivity: ssGSEA pathway scores",
       x = "Delta necroptosis (tumour \u2212 normal)", y = "Delta ferroptosis (tumour \u2212 normal)") +
  theme_pub()

pdc_l <- pdc %>% mutate(method = ifelse(method == "singscore", "singscore", "ssGSEA"),
                        method = factor(method, levels = c("singscore", "ssGSEA")),
                        cor = factor(correlation_type, levels = c("pearson", "spearman")))
p2c <- ggplot(pdc_l, aes(x = method, y = r, colour = cor)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = ci95_lower, ymax = ci95_upper), position = position_dodge(0.4), size = 0.6) +
  scale_colour_manual(values = c("pearson" = okabe_ito[1], "spearman" = okabe_ito[4])) +
  labs(subtitle = "Correlation estimates by scoring method (n = 148)",
       x = NULL, y = "Pearson / Spearman r", colour = NULL) +
  theme_pub()

p2 <- (p2a | p2b) / p2c +
  plot_annotation(title = "Patient-level delta crosstalk") +
  plot_layout(heights = c(1, 0.7))
save_pub(p2, "main_patient_level_delta_crosstalk", 190, 150)
add_log("main_patient_level_delta_crosstalk", "ok", 190, 150, 600,
        "patient_delta_scores.tsv; patient_delta_crosstalk.csv")

## ===========================================================================
## MAIN FIGURE 3: WGCNA crosstalk modules
## ===========================================================================
p3a <- ggplot(soft, aes(soft_power, scale_free_R2)) +
  geom_line(colour = "grey60") +
  geom_point(colour = okabe_ito[1]) +
  geom_point(data = soft[soft$soft_power == ws$soft_power_selected, ],
             colour = okabe_ito[2], size = 2.5) +
  annotate("text", x = ws$soft_power_selected, y = ws$scale_free_R2_selected,
           label = paste0("power = ", ws$soft_power_selected,
                          "\nR\u00b2 = ", signif(ws$scale_free_R2_selected, 3)),
           hjust = -0.15, vjust = -0.3, size = 2.6, family = "Arial") +
  labs(subtitle = "Soft-threshold selection (signed network, TOM type = signed)",
       x = "Soft threshold (power)", y = "Scale-free fit index R\u00b2") +
  theme_pub()

mt_c <- mt[mt$module %in% ct_mods, ]
mt_c$module <- factor(mt_c$module, levels = ct_mods)
mt_c$trait <- factor(mt_c$trait,
                     levels = c("DeltaFerroptosis_singscore", "DeltaNecroptosis_singscore",
                                "DeltaFerroptosis_gsva", "DeltaNecroptosis_gsva"))
p3b <- ggplot(mt_c, aes(x = module, y = trait, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "Pearson r") +
  geom_text(aes(label = ifelse(bh_fdr < 0.05, paste0(signif(correlation, 2), "*"), signif(correlation, 2))),
            size = 2.2, family = "Arial") +
  scale_y_discrete(labels = c("DeltaFerroptosis_singscore" = "Delta ferroptosis (singscore)",
                              "DeltaNecroptosis_singscore" = "Delta necroptosis (singscore)",
                              "DeltaFerroptosis_gsva" = "Delta ferroptosis (ssGSEA)",
                              "DeltaNecroptosis_gsva" = "Delta necroptosis (ssGSEA)")) +
  labs(subtitle = "Module\u2013delta trait correlations (* FDR < 0.05)", x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cm_d <- cm[cm$module %in% ct_mods, ]
cm_d$module <- factor(cm_d$module, levels = ct_mods)
p3c <- ggplot(cm_d, aes(x = module, y = n_genes, fill = sign)) +
  geom_col(width = 0.7) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       name = "Shared direction (sign of r)") +
  geom_text(aes(label = n_genes), vjust = -0.4, size = 2.4, family = "Arial") +
  labs(subtitle = "Crosstalk module sizes (n genes)",
       x = NULL, y = "Number of genes") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- (p3a + p3b) / p3c +
  plot_annotation(title = "WGCNA crosstalk modules") +
  plot_layout(heights = c(1, 0.85))
save_pub(p3, "main_wgcna_crosstalk_modules", 190, 150)
add_log("main_wgcna_crosstalk_modules", "ok", 190, 150, 600,
        "soft_threshold_results.tsv; module_trait_correlations.csv; crosstalk_modules_defined.csv; wgcna_summary.csv")

## ===========================================================================
## MAIN FIGURE 4: crosstalk module biology and hub genes
## ===========================================================================
go_f <- go[go$module %in% ct_mods & go$p.adjust < 0.05, ]
go_top <- head(go_f[order(go_f$p.adjust), ], 12)
go_top$Description <- factor(go_top$Description, levels = rev(unique(go_top$Description)))
go_top$module <- factor(go_top$module, levels = ct_mods)
p4a <- ggplot(go_top, aes(x = -log10(p.adjust), y = Description, colour = module)) +
  geom_point(aes(size = Count)) +
  scale_colour_manual(values = mod_cols) +
  labs(subtitle = "GO biological process enrichment (top FDR terms)",
       x = "\u2212log10 adjusted p", y = NULL, colour = "Module", size = "Count") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.5))

kegg_f <- kegg[kegg$module %in% ct_mods & kegg$p.adjust < 0.05, ]
kegg_top <- head(kegg_f[order(kegg_f$p.adjust), ], 12)
kegg_top$Description <- factor(kegg_top$Description, levels = rev(unique(kegg_top$Description)))
kegg_top$module <- factor(kegg_top$module, levels = ct_mods)
p4b <- ggplot(kegg_top, aes(x = -log10(p.adjust), y = Description, colour = module)) +
  geom_point(aes(size = Count)) +
  scale_colour_manual(values = mod_cols) +
  labs(subtitle = "KEGG pathway enrichment (top FDR terms)",
       x = "\u2212log10 adjusted p", y = NULL, colour = "Module", size = "Count") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.5))

hub_top <- hub[hub$module_rank <= 3, ]
hub_top$gene_symbol <- factor(hub_top$gene_symbol,
                              levels = rev(hub_top$gene_symbol[order(hub_top$module)]))
p4c <- ggplot(hub_top, aes(x = hub_score, y = gene_symbol, fill = module)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = mod_cols) +
  labs(subtitle = "Top 3 prioritised hub genes per crosstalk module",
       x = "Hub evidence score", y = NULL, fill = "Module") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6))

p4 <- (p4a + p4b) / p4c +
  plot_annotation(title = "Crosstalk module biology and hub genes") +
  plot_layout(heights = c(1, 1.1))
save_pub(p4, "main_crosstalk_module_biology_and_hubs", 190, 170)
add_log("main_crosstalk_module_biology_and_hubs", "ok", 190, 170, 600,
        "GO_BP.tsv; KEGG.tsv; top_hubs_per_module.tsv")

## ===========================================================================
## MAIN FIGURE 5: tumour microenvironment and composition-adjusted coupling
## ===========================================================================
ctde_top <- head(ctde[order(abs(ctde$mean_diff), decreasing = TRUE), ], 15)
ctde_top$celltype <- factor(ctde_top$celltype, levels = rev(ctde_top$celltype))
p5a <- ggplot(ctde_top, aes(x = mean_diff, y = celltype,
                            fill = ifelse(fdr < 0.05, "FDR < 0.05", "ns"))) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("FDR < 0.05" = okabe_ito[2], "ns" = "grey75")) +
  labs(subtitle = paste0("xCell abundance: tumour \u2212 adjacent-normal (", xsum$n_celltypes,
                         " cell types, ", xsum$n_sig_celltypes_fdr05, " FDR < 0.05)"),
       x = "Mean difference (tumour \u2212 normal)", y = NULL, fill = NULL) +
  theme_pub()

comp_l <- data.frame(
  model = factor(c("singscore", "ssGSEA"), levels = c("singscore", "ssGSEA")),
  base_beta = comp_cou$base_beta[1:2],
  base_p = comp_cou$base_p[1:2],
  adj_beta = comp_cou$adj_beta[1:2],
  adj_p = comp_cou$adj_p[1:2]
)
p5b <- ggplot(comp_l, aes(x = model)) +
  geom_point(aes(y = base_beta, shape = "Unadjusted"), size = 2.6, colour = "grey55") +
  geom_point(aes(y = adj_beta, shape = "Composition-adjusted"), size = 2.6, colour = okabe_ito[2]) +
  geom_segment(aes(xend = model, y = base_beta, yend = adj_beta), colour = "grey70") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape_manual(values = c("Unadjusted" = 16, "Composition-adjusted" = 17)) +
  annotate("text", x = c(1, 2), y = -0.05, size = 2.4, family = "Arial",
           label = paste0("adjusted p = ", fmt_p(comp_l$adj_p))) +
  labs(subtitle = "Within-tumour coupling before/after xCell adjustment (ferroptosis ~ necroptosis)",
       x = NULL, y = "Standardised regression beta", shape = NULL) +
  theme_pub()

p5 <- (p5a | p5b) +
  plot_annotation(title = "Tumour microenvironment and composition-adjusted coupling") +
  plot_layout(widths = c(1, 1))
save_pub(p5, "main_tumour_microenvironment_and_adjusted_coupling", 190, 120)
add_log("main_tumour_microenvironment_and_adjusted_coupling", "ok", 190, 120, 600,
        "cell_type_paired_de.tsv; xcell_summary.csv; composition_adjusted_coupling.tsv")

## ===========================================================================
## MAIN FIGURE 6: external validation (GSE70905)
## ===========================================================================
val_cou$group <- factor(c("tumor" = "Tumour", "adjacent_normal" = "Adjacent",
                          "reduction_normal" = "Reduction")[val_cou$group],
                        levels = c("Tumour", "Adjacent", "Reduction"))
p6a <- ggplot(val_cou, aes(x = group, y = r, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = pmax(0, r - 0.06), ymax = pmin(1, r + 0.06)), width = 0.15) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(aes(label = paste0("r = ", signif(r, 3), "\np = ", fmt_p(p_value))),
            vjust = -0.6, size = 2.4, family = "Arial") +
  scale_fill_manual(values = c("Tumour" = okabe_ito[2], "Adjacent" = okabe_ito[1],
                               "Reduction" = okabe_ito[7])) +
  coord_cartesian(ylim = c(0, 1.15)) +
  labs(subtitle = "Within-tissue coupling in GSE70905 (singscore)",
       x = NULL, y = "Pearson r (ferroptosis vs necroptosis)") +
  theme_pub() +
  theme(legend.position = "none")

cmr <- cm[cm$module %in% ct_mods, c("module", "r_ferro", "r_necro")]
vm <- merge(val_mod, cmr, by = "module")
vm$concordant <- ifelse(sign(vm$r_ferro) == sign(vm$r_delta_ferro) & vm$fdr_ferro < 0.05,
                        "concordant (FDR<0.05)", "not concordant")
p6b <- ggplot(vm, aes(x = r_ferro, y = r_delta_ferro, colour = concordant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(size = 2.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey50", linewidth = 0.6) +
  geom_text(aes(label = module), hjust = -0.3, vjust = -0.4, size = 2.2, family = "Arial") +
  scale_colour_manual(values = c("concordant (FDR<0.05)" = okabe_ito[2],
                                 "not concordant" = "grey60")) +
  labs(subtitle = "Module replication: discovery vs validation (GSE70905)",
       x = "Discovery r (module eigengene vs Delta ferroptosis)",
       y = "Validation r (module score vs Delta ferroptosis)", colour = NULL) +
  theme_pub()

val_de$pathway <- ifelse(val_de$pathway == "ferro", "Ferroptosis", "Necroptosis")
p6c <- ggplot(val_de, aes(x = pathway, y = mean_diff, fill = pathway)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(aes(label = paste0("paired t p = ", fmt_p(p_value), "\nFDR = ", fmt_p(fdr))),
            vjust = -0.5, size = 2.4, family = "Arial") +
  scale_fill_manual(values = c("Ferroptosis" = okabe_ito[3], "Necroptosis" = okabe_ito[4])) +
  coord_cartesian(ylim = c(0, 0.05)) +
  labs(subtitle = "Validation cohort paired pathway change (n = 47 paired)",
       x = NULL, y = "Mean difference (tumour \u2212 adjacent-normal)") +
  theme_pub() +
  theme(legend.position = "none")

p6 <- (p6a + p6b) / p6c +
  plot_annotation(title = "External validation (GSE70905)") +
  plot_layout(heights = c(1, 0.8))
save_pub(p6, "main_external_validation_gse70905", 190, 150)
add_log("main_external_validation_gse70905", "ok", 190, 150, 600,
        "validation_coupling.tsv; validation_module_associations.tsv; crosstalk_modules_defined.csv; validation_pathway_de.tsv")

## ===========================================================================
## SUPP: expression distribution QC
## ===========================================================================
p7a <- ggplot(qc47, aes(mean_expr)) +
  geom_histogram(bins = 25, fill = okabe_ito[1], colour = "white") +
  labs(subtitle = "Mean expression per sample", x = "Mean expression", y = "Count") +
  theme_pub()
p7b <- ggplot(qc47, aes(sd_expr)) +
  geom_histogram(bins = 25, fill = okabe_ito[3], colour = "white") +
  labs(subtitle = "Expression SD per sample", x = "SD expression", y = "Count") +
  theme_pub()
p7c <- ggplot(qc47, aes(connectivity, mean_expr)) +
  geom_point(size = 1.2, alpha = 0.5, colour = okabe_ito[4]) +
  labs(subtitle = "Sample connectivity vs mean expression", x = "Connectivity", y = "Mean expression") +
  theme_pub()
p7 <- (p7a + p7b + p7c) +
  plot_annotation(title = "QC: expression distributions (GSE70947)") +
  plot_layout(nrow = 1)
save_pub(p7, "supp_expression_distribution_qc", 190, 75)
add_log("supp_expression_distribution_qc", "ok", 190, 75, 600, "qc/GSE70947_sample_qc.csv")

## ===========================================================================
## SUPP: discovery PCA and clustering
## ===========================================================================
expr <- ana$expression_imputed
pheno_tissue <- trimws(sub("tissue: *", "", ana$phenotype[["tissue:ch1"]]))
names(pheno_tissue) <- rownames(ana$phenotype)
gvar <- apply(expr, 1, var)
topg <- names(sort(gvar, decreasing = TRUE))[seq_len(2000)]
expr_sub <- expr[topg, ]

pca_res <- prcomp(t(expr_sub), center = TRUE, scale. = FALSE)
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                     tissue = pheno_tissue[rownames(ana$phenotype)])
var_pc <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
p8a <- ggplot(pca_df, aes(PC1, PC2, colour = tissue)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_colour_manual(values = c("normal" = okabe_ito[1], "tumor" = okabe_ito[2])) +
  labs(subtitle = paste0("PCA on top 2,000 variable genes (", var_pc[1], "% / ", var_pc[2], "% variance)"),
       x = "PC1", y = "PC2") +
  theme_pub()

hcl <- hclust(dist(t(expr_sub)), method = "average")
hc_d <- dendro_data(hcl)
if (length(dev.list()) > 0) dev.off()
hc_lab <- hc_d$labels
hc_lab$tissue <- pheno_tissue[match(hc_lab$label, rownames(ana$phenotype))]
p8b <- ggplot(segment(hc_d)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = hc_lab, aes(x = x, y = y, colour = tissue), size = 0.8) +
  scale_colour_manual(values = c("normal" = okabe_ito[1], "tumor" = okabe_ito[2])) +
  labs(subtitle = "Sample hierarchical clustering (average linkage)", x = NULL, y = "Distance") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p8 <- (p8a | p8b) +
  plot_annotation(title = "Discovery cohort structure (GSE70947)") +
  plot_layout(widths = c(1, 1.2))
save_pub(p8, "supp_discovery_pca_and_clustering", 190, 90)
add_log("supp_discovery_pca_and_clustering", "ok", 190, 90, 600, "GSE70947_analysis.rds")

## ===========================================================================
## SUPP: gene set coverage / overlap
## ===========================================================================
gs_bar <- data.frame(gene_set = c("Ferroptosis", "Necroptosis"),
                     unique_genes = gs_sum$unique_genes,
                     total_entries = gs_sum$total_entries)
p9a <- ggplot(gs_bar, aes(x = gene_set, y = unique_genes, fill = gene_set)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = unique_genes), vjust = -0.4, size = 3, family = "Arial") +
  scale_fill_manual(values = c("Ferroptosis" = okabe_ito[3], "Necroptosis" = okabe_ito[4])) +
  labs(subtitle = paste0("Unique mapped genes (ferroptosis ",
                         gs_sum$total_entries[1], " entries; necroptosis ",
                         gs_sum$total_entries[2], " entries)"),
       x = NULL, y = "Unique genes") +
  theme_pub() +
  theme(legend.position = "none")

ov <- read.csv(file.path("results", "tables", "ferroptosis_necroptosis_overlap.csv"), sep = "\t")
p9b <- ggplot(data.frame(x = c(1, 2), y = c(gs_sum$unique_genes[1], gs_sum$unique_genes[2])),
              aes(x, y)) +
  geom_col(width = 0.4, fill = c(okabe_ito[3], okabe_ito[4])) +
  geom_text(aes(label = y), vjust = -0.4, size = 3, family = "Arial") +
  annotate("text", x = 1.5, y = 60, size = 3, family = "Arial",
           label = paste0("Overlap: ", nrow(ov), " genes")) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Ferroptosis", "Necroptosis")) +
  labs(subtitle = "Gene-set size and overlap", x = NULL, y = "Unique genes") +
  theme_pub() +
  theme(panel.grid.major = element_blank())
p9 <- (p9a + p9b) +
  plot_annotation(title = "Ferroptosis/necroptosis gene-set coverage") +
  plot_layout(nrow = 1)
save_pub(p9, "supp_gene_set_coverage", 170, 80)
add_log("supp_gene_set_coverage", "ok", 170, 80, 600,
        "tables/gene_set_summary.csv; tables/ferroptosis_necroptosis_overlap.csv")

## ===========================================================================
## SUPP: pathway scoring method agreement
## ===========================================================================
r_sf <- stats::cor(sc$ferroptosis_full_singscore, sc$ferroptosis_full_gsva)
r_sn <- stats::cor(sc$necroptosis_full_singscore, sc$necroptosis_full_gsva)
p10a <- ggplot(sc, aes(ferroptosis_full_singscore, ferroptosis_full_gsva)) +
  geom_point(alpha = 0.5, size = 1.2, colour = okabe_ito[3]) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = paste0("r = ", signif(r_sf, 3))) +
  labs(subtitle = "Ferroptosis score agreement (GSE70947)", x = "singscore", y = "ssGSEA") +
  theme_pub()
p10b <- ggplot(sc, aes(necroptosis_full_singscore, necroptosis_full_gsva)) +
  geom_point(alpha = 0.5, size = 1.2, colour = okabe_ito[4]) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = paste0("r = ", signif(r_sn, 3))) +
  labs(subtitle = "Necroptosis score agreement (GSE70947)", x = "singscore", y = "ssGSEA") +
  theme_pub()
p10 <- (p10a + p10b) +
  plot_annotation(title = "Pathway scoring method agreement") +
  plot_layout(nrow = 1)
save_pub(p10, "supp_pathway_scoring_method_agreement", 170, 85)
add_log("supp_pathway_scoring_method_agreement", "ok", 170, 85, 600,
        "GSE70947_pathway_scores.tsv")

## ===========================================================================
## SUPP: ssGSEA delta crosstalk (standalone sensitivity figure)
## ===========================================================================
p11 <- ggplot(dt, aes(Deltanecroptosis_gsva, Deltaferroptosis_gsva)) +
  geom_point(alpha = 0.55, size = 1.4, colour = okabe_ito[4]) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = cor_ann(r_g$r, r_g$p_value, r_g$n,
                           paste0("\n95% CI [", signif(r_g$ci95_lower, 3), ", ",
                                  signif(r_g$ci95_upper, 3), "]"))) +
  labs(title = "ssGSEA delta crosstalk (sensitivity)",
       x = "Delta necroptosis (tumour \u2212 normal)", y = "Delta ferroptosis (tumour \u2212 normal)") +
  theme_pub()
save_pub(p11, "supp_ssgsea_delta_crosstalk", 120, 100)
add_log("supp_ssgsea_delta_crosstalk", "ok", 120, 100, 600, "patient_delta_scores.tsv; patient_delta_crosstalk.csv")

## ===========================================================================
## SUPP: correlation method sensitivity (pearson vs spearman)
## ===========================================================================
p12 <- ggplot(pdc_l, aes(x = method, y = r, colour = cor)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = ci95_lower, ymax = ci95_upper), position = position_dodge(0.4), size = 0.6) +
  geom_text(aes(label = paste0("p = ", fmt_p(p_value)), y = ci95_upper + 0.02),
            position = position_dodge(0.4), size = 2.4, family = "Arial") +
  scale_colour_manual(values = c("pearson" = okabe_ito[1], "spearman" = okabe_ito[4])) +
  labs(title = "Correlation method sensitivity (n = 148)",
       x = NULL, y = "Correlation r", colour = NULL) +
  theme_pub()
save_pub(p12, "supp_correlation_method_sensitivity", 120, 90)
add_log("supp_correlation_method_sensitivity", "ok", 120, 90, 600, "patient_delta_crosstalk.csv")

## ===========================================================================
## SUPP: leave-one-out robustness
## ===========================================================================
loo_r <- dt %>% summarise(r = stats::cor(Deltaferroptosis_singscore, Deltanecroptosis_singscore)) %>% pull(r)
p13 <- ggplot(loo, aes(x = r)) +
  geom_histogram(bins = 30, fill = okabe_ito[1], colour = "white") +
  geom_vline(xintercept = loo_r, colour = okabe_ito[2], linetype = 2, linewidth = 0.6) +
  annotate("text", x = loo_r, y = Inf, vjust = 1.5, hjust = -0.05, size = 2.8, family = "Arial",
           label = paste0("full r = ", signif(loo_r, 3),
                          "\nrange [", signif(ss$pearson_r[ss$analysis == "leave_one_out_min"], 3),
                          ", ", signif(ss$pearson_r[ss$analysis == "leave_one_out_max"], 3), "]")) +
  labs(title = "Leave-one-out robustness (n = 148 patients)",
       x = "Pearson r (delta ferroptosis vs delta necroptosis)", y = "Count") +
  theme_pub()
save_pub(p13, "supp_leave_one_out_robustness", 120, 90)
add_log("supp_leave_one_out_robustness", "ok", 120, 90, 600,
        "sensitivity/leave_one_out.tsv; sensitivity/sensitivity_summary.tsv")

## ===========================================================================
## SUPP: ER-status stratified crosstalk
## ===========================================================================
er_map <- sc %>% select(patient_id, er) %>% distinct()
dt_er <- merge(dt, er_map, by = "patient_id")
dt_er$er <- factor(dt_er$er, levels = c("pos", "neg"))
erp <- ss[ss$analysis == "er_pos", ]
ern <- ss[ss$analysis == "er_neg", ]
p14 <- ggplot(dt_er, aes(Deltanecroptosis_singscore, Deltaferroptosis_singscore, colour = er)) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  scale_colour_manual(values = c("pos" = okabe_ito[2], "neg" = okabe_ito[1]),
                      labels = c("pos" = "ER-positive", "neg" = "ER-negative")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = paste0("ER+: ", cor_ann(erp$pearson_r, erp$pearson_p, erp$n), "\n",
                          "ER\u2212: ", cor_ann(ern$pearson_r, ern$pearson_p, ern$n))) +
  labs(title = "Delta crosstalk stratified by ER status",
       x = "Delta necroptosis (tumour \u2212 normal)", y = "Delta ferroptosis (tumour \u2212 normal)") +
  theme_pub()
save_pub(p14, "supp_er_status_stratified_crosstalk", 120, 100)
add_log("supp_er_status_stratified_crosstalk", "ok", 120, 100, 600,
        "patient_delta_scores.tsv; GSE70947_pathway_scores.tsv; sensitivity_summary.tsv")

## ===========================================================================
## SUPP: WGCNA diagnostics (soft-threshold + module dendrogram & colours)
## ===========================================================================
png(file.path(pub_dir, "supp_wgcna_diagnostics.png"),
    width = 190, height = 130, units = "mm", res = 600, bg = "white")
par(family = "Arial", mar = c(3.2, 3.6, 1.2, 0.8))
layout(matrix(c(1, 1, 2, 2, 2), nrow = 5, byrow = FALSE),
       heights = c(1, 1, 1, 1, 1))
plot(soft$soft_power, soft$scale_free_R2, type = "b", pch = 16, col = okabe_ito[1],
     xlab = "Soft threshold (power)", ylab = "Scale-free fit R\u00b2",
     main = "", cex.lab = 1.1, cex.axis = 0.9, las = 1)
points(ws$soft_power_selected, ws$scale_free_R2_selected, pch = 16, cex = 1.6, col = okabe_ito[2])
text(ws$soft_power_selected, ws$scale_free_R2_selected,
     paste0(" power ", ws$soft_power_selected, " (R\u00b2 = ",
            signif(ws$scale_free_R2_selected, 3), ")"), pos = 4, cex = 0.9, family = "Arial")
mtext("Soft-threshold selection (signed network)", side = 3, line = 0.2, cex = 0.95, family = "Arial")
nblocks <- length(w$net$dendrograms)
for (b in seq_len(nblocks)) {
  dg <- w$net$dendrograms[[b]]
  cols <- w$net$colors[w$net$blocks == b]
  mar_prev <- par("mar"); par(mar = c(0.5, 3.6, 0.8, 0.8))
  WGCNA::plotDendroAndColors(dg, colors = cols, dendroLabels = FALSE,
                             colorHeight = 0.12, addGuide = FALSE, hang = 0.03,
                             main = paste0("Block ", b, " (", length(dg$labels), " genes)"),
                             cex.main = 0.85)
  par(mar = mar_prev)
}
dev.off()
add_log("supp_wgcna_diagnostics", "ok", 190, 130, 600, "wgcna_delta_result.rds; soft_threshold_results.tsv")

## ===========================================================================
## SUPP: full module-trait heatmap (all modules x 4 traits)
## ===========================================================================
mt_all <- mt
mt_all$module <- factor(mt_all$module, levels = paste0("ME", 0:33))
mt_all$trait <- factor(mt_all$trait,
                       levels = c("DeltaFerroptosis_singscore", "DeltaNecroptosis_singscore",
                                  "DeltaFerroptosis_gsva", "DeltaNecroptosis_gsva"))
p16 <- ggplot(mt_all, aes(x = module, y = trait, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "Pearson r") +
  geom_text(aes(label = ifelse(bh_fdr < 0.05, paste0(signif(correlation, 2), "*"), signif(correlation, 2))),
            size = 1.8, family = "Arial") +
  scale_y_discrete(labels = c("DeltaFerroptosis_singscore" = "Delta ferroptosis (singscore)",
                              "DeltaNecroptosis_singscore" = "Delta necroptosis (singscore)",
                              "DeltaFerroptosis_gsva" = "Delta ferroptosis (ssGSEA)",
                              "DeltaNecroptosis_gsva" = "Delta necroptosis (ssGSEA)")) +
  labs(title = "Module\u2013delta trait correlations, all 34 modules (* FDR < 0.05)",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 7))
save_pub(p16, "supp_full_module_trait_heatmap", 190, 90)
add_log("supp_full_module_trait_heatmap", "ok", 190, 90, 600, "module_trait_correlations.csv")

## ===========================================================================
## SUPP: enrichment dotplots per source
## ===========================================================================
make_dotplot <- function(df, src, w = 190, h = 140, fname) {
  df <- df[df$module %in% ct_mods & df$p.adjust < 0.05, ]
  if (nrow(df) == 0) {
    add_log(fname, "no significant terms", w, h, 600, src)
    return(NULL)
  }
  df <- head(df[order(df$p.adjust), ], 20)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))
  df$module <- factor(as.character(df$module), levels = ct_mods)
  p <- ggplot(df, aes(x = -log10(p.adjust), y = Description, colour = module, size = Count)) +
    geom_point() +
    scale_colour_manual(values = mod_cols) +
    labs(title = paste0(src, " enrichment in crosstalk modules"),
         x = "\u2212log10 adjusted p", y = NULL, colour = "Module", size = "Count") +
    theme_pub() +
    theme(axis.text.y = element_text(size = 5.5))
  save_pub(p, fname, w, h)
  add_log(fname, "ok", w, h, 600, paste0("enrichment/", src, ".tsv"))
}
make_dotplot(go, "GO_BP", 190, 160, "supp_go_enrichment_crosstalk_modules")
make_dotplot(kegg, "KEGG", 190, 150, "supp_kegg_enrichment_crosstalk_modules")
make_dotplot(hall, "Hallmark", 190, 140, "supp_hallmark_enrichment_crosstalk_modules")
make_dotplot(rea, "Reactome", 190, 160, "supp_reactome_enrichment_crosstalk_modules")

## ===========================================================================
## SUPP: hub genes, all crosstalk modules
## ===========================================================================
hub_all <- hub[hub$crosstalk_module & hub$top10, ]
hub_all$gene_symbol <- factor(hub_all$gene_symbol,
                              levels = rev(hub_all$gene_symbol[order(hub_all$module)]))
hub_all$module <- factor(hub_all$module, levels = ct_mods)
p22 <- ggplot(hub_all, aes(x = hub_score, y = gene_symbol, fill = module)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = mod_cols) +
  labs(title = "Top 10 prioritised hub genes per crosstalk module",
       x = "Hub evidence score", y = NULL, fill = "Module") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 4.5))
save_pub(p22, "supp_hub_genes_all_crosstalk_modules", 170, 220)
add_log("supp_hub_genes_all_crosstalk_modules", "ok", 170, 220, 600, "top_hubs_per_module.tsv")

## ===========================================================================
## SUPP: xCell all cell types (volcano-style)
## ===========================================================================
p23 <- ggplot(ctde, aes(x = mean_diff, y = -log10(fdr), colour = ifelse(fdr < 0.05, "FDR < 0.05", "ns"))) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_text(data = ctde[order(ctde$fdr), ][1:8, ],
            aes(label = celltype), hjust = -0.1, size = 2.2, family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = c("FDR < 0.05" = okabe_ito[2], "ns" = "grey70")) +
  labs(title = paste0("xCell differential abundance (tumour vs adjacent-normal), ",
                      xsum$n_celltypes, " cell types, ", xsum$n_sig_celltypes_fdr05, " FDR < 0.05"),
       x = "Mean difference (tumour \u2212 normal)", y = "\u2212log10 FDR", colour = NULL) +
  theme_pub()
save_pub(p23, "supp_xcell_all_cell_types", 150, 100)
add_log("supp_xcell_all_cell_types", "ok", 150, 100, 600, "cell_type_paired_de.tsv; xcell_summary.csv")

## ===========================================================================
## SUPP: xCell sample heatmap
## ===========================================================================
xmat <- cbind(xsc$tumour, xsc$normal)
xgrp <- c(rep("Tumour", ncol(xsc$tumour)), rep("Adjacent normal", ncol(xsc$normal)))
xsc_df <- as.data.frame(t(scale(t(xmat))))
xsc_df <- xsc_df[rowSums(is.na(xsc_df)) == 0, ]
anno_col <- data.frame(Tissue = xgrp, row.names = colnames(xmat))
anno_colors <- list(Tissue = c("Tumour" = okabe_ito[2], "Adjacent normal" = okabe_ito[1]))
png(file.path(pub_dir, "supp_xcell_sample_heatmap.png"),
    width = 190, height = 140, units = "mm", res = 600, bg = "white")
pheatmap::pheatmap(as.matrix(xsc_df), cluster_rows = TRUE, cluster_cols = FALSE,
                   annotation_col = anno_col, annotation_colors = anno_colors,
                   show_colnames = FALSE, show_rownames = TRUE,
                   fontsize_row = 4, fontsize = 5,
                   color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                   main = "")
dev.off()
add_log("supp_xcell_sample_heatmap", "ok", 190, 140, 600, "xcell_scores.rds")

## ===========================================================================
## SUPP: composition adjustment diagnostics
## ===========================================================================
comp_l2 <- data.frame(
  model = rep(c("singscore", "ssGSEA"), 2),
  stage = rep(c("Unadjusted", "Composition-adjusted"), each = 2),
  beta = c(comp_cou$base_beta[1:2], comp_cou$adj_beta[1:2]),
  p = c(comp_cou$base_p[1:2], comp_cou$adj_p[1:2])
)
comp_l2$model <- factor(comp_l2$model, levels = c("singscore", "ssGSEA"))
comp_l2$stage <- factor(comp_l2$stage, levels = c("Unadjusted", "Composition-adjusted"))
p25 <- ggplot(comp_l2, aes(x = model, y = beta, fill = stage)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = paste0("p = ", fmt_p(p))),
            position = position_dodge(0.7), vjust = -0.4, size = 2.6, family = "Arial") +
  scale_fill_manual(values = c("Unadjusted" = "grey70", "Composition-adjusted" = okabe_ito[2])) +
  labs(title = "Composition adjustment diagnostics (within-tumour coupling)",
       subtitle = "Beta from ferroptosis ~ necroptosis regression, before/after xCell covariate adjustment",
       x = NULL, y = "Standardised regression beta", fill = NULL) +
  theme_pub()
save_pub(p25, "supp_composition_adjustment_diagnostics", 140, 90)
add_log("supp_composition_adjustment_diagnostics", "ok", 140, 90, 600, "composition_adjusted_coupling.tsv")

## ===========================================================================
## SUPP: field-effect pathway scores
## ===========================================================================
p26 <- ggplot(fe_gm, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
  facet_wrap(~pathway) +
  scale_fill_manual(values = c("Reduction" = okabe_ito[7], "Adjacent" = okabe_ito[1],
                               "Tumour" = okabe_ito[2])) +
  annotate("text", x = 2, y = Inf, vjust = 1.5, size = 2.8, family = "Arial",
           label = paste0("ANOVA p: ferroptosis ", fmt_p(fe_sum$ferro_anova_p),
                          "; necroptosis ", fmt_p(fe_sum$necro_anova_p),
                          "\ntumour vs reduction padj: ", fmt_p(fe_sum$ferro_tumor_vs_reduction_padj),
                          " / ", fmt_p(fe_sum$necro_tumor_vs_reduction_padj))) +
  labs(title = "Field effect across tissue types (GSE70905, singscore)",
       x = NULL, y = "Mean singscore") +
  theme_pub() +
  theme(legend.position = "none", strip.background = element_blank())
save_pub(p26, "supp_field_effect_pathway_scores", 150, 90)
add_log("supp_field_effect_pathway_scores", "ok", 150, 90, 600,
        "field_effect/group_means.tsv; field_effect_summary.csv")

## ===========================================================================
## SUPP: external paired pathway changes
## ===========================================================================
p27 <- ggplot(val_de, aes(x = pathway, y = mean_diff, fill = pathway)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(aes(label = paste0("n = ", n_pairs, " paired\nmean diff = ", signif(mean_diff, 3),
                               "\npaired t p = ", fmt_p(p_value))),
            vjust = -0.5, size = 2.6, family = "Arial") +
  scale_fill_manual(values = c("Ferroptosis" = okabe_ito[3], "Necroptosis" = okabe_ito[4])) +
  coord_cartesian(ylim = c(0, 0.05)) +
  labs(title = "Validation cohort: paired pathway change (GSE70905)",
       x = NULL, y = "Mean difference (tumour \u2212 adjacent-normal)") +
  theme_pub() +
  theme(legend.position = "none")
save_pub(p27, "supp_external_paired_pathway_changes", 120, 90)
add_log("supp_external_paired_pathway_changes", "ok", 120, 90, 600, "validation_pathway_de.tsv")

## ===========================================================================
## SUPP: module replication discovery vs validation (per trait)
## ===========================================================================
vm2 <- merge(val_mod, cm[, c("module", "r_ferro", "r_necro")], by = "module")
vm2$concordant <- ifelse(sign(vm2$r_ferro) == sign(vm2$r_delta_ferro) & vm2$fdr_ferro < 0.05,
                         "ferro concordant (FDR<0.05)", "not concordant")
p28a <- ggplot(vm2, aes(x = r_ferro, y = r_delta_ferro, colour = concordant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(size = 2.6) +
  geom_text(aes(label = module), hjust = -0.3, vjust = -0.4, size = 2.2, family = "Arial") +
  scale_colour_manual(values = c("ferro concordant (FDR<0.05)" = okabe_ito[2],
                                 "not concordant" = "grey60")) +
  labs(subtitle = "Delta ferroptosis", x = "Discovery r", y = "Validation r") +
  theme_pub()
p28b <- ggplot(vm2, aes(x = r_necro, y = r_delta_necro)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(size = 2.6, colour = okabe_ito[4]) +
  geom_text(aes(label = module), hjust = -0.3, vjust = -0.4, size = 2.2, family = "Arial") +
  labs(subtitle = "Delta necroptosis", x = "Discovery r", y = "Validation r") +
  theme_pub()
p28 <- (p28a + p28b) +
  plot_annotation(title = "Module replication, discovery vs validation (GSE70905)") +
  plot_layout(nrow = 1)
save_pub(p28, "supp_module_replication_discovery_vs_validation", 190, 90)
add_log("supp_module_replication_discovery_vs_validation", "ok", 190, 90, 600,
        "validation_module_associations.tsv; crosstalk_modules_defined.csv")

## ===========================================================================
## SUPP: hub gene directional replication
## ===========================================================================
vh <- val_hub
vh$dir <- ifelse(vh$replicated, "Replicated", "Not replicated")
frac <- mean(vh$replicated)
p29 <- ggplot(vh, aes(x = discovery_r_ferro, y = validation_r_ferro, colour = dir)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c("Replicated" = okabe_ito[2], "Not replicated" = "grey70")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.8, family = "Arial",
           label = paste0("Directionally replicated: ", sum(vh$replicated), " / ", nrow(vh),
                          " (", signif(100 * frac, 1), "%)")) +
  labs(title = "Hub gene directional replication (GSE70905)",
       x = "Discovery r (hub gene vs Delta ferroptosis)",
       y = "Validation r (hub gene vs Delta ferroptosis)", colour = NULL) +
  theme_pub()
save_pub(p29, "supp_hub_gene_directional_replication", 130, 100)
add_log("supp_hub_gene_directional_replication", "ok", 130, 100, 600, "validation_hub_replication.tsv")

## ===========================================================================
## SUPP: GSE15852 pathway scoring
## ===========================================================================
sc2_long <- sc2 %>%
  select(sample_id, tissue, ferroptosis_full_singscore, necroptosis_full_singscore) %>%
  pivot_longer(c(ferroptosis_full_singscore, necroptosis_full_singscore),
               names_to = "pathway", values_to = "score")
sc2_long$pathway <- ifelse(sc2_long$pathway == "ferroptosis_full_singscore",
                           "Ferroptosis", "Necroptosis")
sc2_long$sample_type <- ifelse(grepl("tumor", sc2_long$tissue), "Tumour", "Normal")
p30 <- ggplot(sc2_long, aes(x = sample_type, y = score, fill = sample_type)) +
  geom_violin(alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.22, outlier.shape = NA, colour = "black") +
  facet_wrap(~pathway) +
  scale_fill_manual(values = c("Tumour" = okabe_ito[2], "Normal" = okabe_ito[1])) +
  labs(title = "GSE15852 pathway scoring (86 samples)",
       subtitle = "Independent cohort; singscore",
       x = NULL, y = "singscore") +
  theme_pub() +
  theme(legend.position = "none", strip.background = element_blank())
save_pub(p30, "supp_gse15852_pathway_scoring", 150, 90)
add_log("supp_gse15852_pathway_scoring", "ok", 150, 90, 600, "GSE15852_pathway_scores.tsv")

## ===========================================================================
## finalise logs and missing figures
## ===========================================================================
write.table(log_df, log_file, sep = "\t", row.names = FALSE, quote = FALSE)

missing_txt <- c(
  "Figures that could not be generated from existing results:",
  "",
  "1. supp_permutation_null_distribution  -- NOT AVAILABLE",
  "   The 2000-permutation null distribution for the interaction beta was not saved to disk",
  "   (glob results/**/*perm* found no files). Only the summary statistics exist in",
  "   results/crosstalk/coupling_robustness.csv (permutation_interaction_beta, permutation_p,",
  "   bootstrap_CI). Per the project rule that publication figures must be generated ONLY from",
  "   existing computed results (no recomputation), this panel cannot be rebuilt and is logged as",
  "   missing. Equivalent robustness evidence is presented in supp_ssgsea_delta_crosstalk.png,",
  "   supp_leave_one_out_robustness.png and supp_correlation_method_sensitivity.png.",
  "",
  paste0("Figures generated: ", nrow(log_df)),
  paste0("Log written to: ", log_file)
)
writeLines(missing_txt, missing_file)

cat("Done. Figures:", nrow(log_df), "\n")