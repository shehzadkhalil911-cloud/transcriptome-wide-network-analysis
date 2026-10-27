## 15_figures.R
## Generate all manuscript figures (Fig 1-10) and supplementary figures.
## Figures:
##   Fig 1  Study overview / sample flowchart
##   Fig 2  QC: PCA + sample distributions
##   Fig 3  Pathway scores by tissue (paired) + delta distributions
##   Fig 4  Within-tissue ferroptosis-necroptosis coupling (scatter)
##   Fig 5  Patient-level delta crosstalk (scatter, primary result)
##   Fig 6  WGCNA: soft-threshold + module-trait heatmap
##   Fig 7  Crosstalk module eigengene-delta correlations (bar)
##   Fig 8  Enrichment: top KEGG/Hallmark for crosstalk modules
##   Fig 9  Hub gene prioritisation (top hubs bar)
##   Fig 10 External validation (GSE70905 coupling + module replication)
##   Supp A  Leave-one-out sensitivity
##   Supp B  Field-effect (tumour vs adjacent vs reduction)
##   Supp C  Composition-adjusted coupling
## Outputs (PNG + PDF): results/figures/FigX_*.png / *.pdf

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
fig_dir <- .dir_for("results_figures")
supp_dir <- .dir_for("results_supp")

# helper to save both formats
save_fig <- function(base, w = 8, h = 6, dpi = 300) {
  ggsave(file.path(fig_dir, paste0(base, ".png")), width = w, height = h, dpi = dpi)
  ggsave(file.path(fig_dir, paste0(base, ".pdf")), width = w, height = h)
  ggsave(file.path(supp_dir, paste0(base, ".png")), width = w, height = h, dpi = dpi)
}

theme_ms <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(face = "bold", hjust = 0),
          legend.position = "top")
}

## ---- load data ----
ps <- readRDS(file.path(.dir_for("data_processed"), "pathway_scores_all.rds"))
sc <- ps$scores_GSE70947
ana <- readRDS(file.path(.dir_for("data_processed"), "GSE70947_analysis.rds"))

dt <- read.csv(file.path("results", "crosstalk", "patient_delta_scores.tsv"), sep = "\t")
w <- readRDS(file.path(.dir_for("data_processed"), "wgcna_delta_result.rds"))
mt <- read.csv(file.path("results", "wgcna", "module_trait_correlations.csv"), sep = "\t")
soft <- read.csv(file.path("results", "wgcna", "soft_threshold_results.tsv"), sep = "\t")
mods <- read.csv(file.path("results", "wgcna", "module_assignments.tsv"), sep = "\t")
hub <- read.csv(file.path("results", "hubs", "hub_gene_evidence.tsv"), sep = "\t")
val_sum <- read.csv(file.path("results", "validation", "validation_summary.csv"), sep = "\t")
val_cou <- read.csv(file.path("results", "validation", "validation_coupling.tsv"), sep = "\t")
val_mod <- read.csv(file.path("results", "validation", "validation_module_associations.tsv"), sep = "\t")
fe_sum <- read.csv(file.path("results", "field_effect", "field_effect_summary.csv"), sep = "\t")
sens <- read.csv(file.path("results", "sensitivity", "sensitivity_summary.tsv"), sep = "\t")
loo <- read.csv(file.path("results", "sensitivity", "leave_one_out.tsv"), sep = "\t")
comp_sum <- read.csv(file.path("results", "composition", "xcell_summary.csv"), sep = "\t")

## ---- Fig 1: study overview (schematic; data-driven counts) ----
tab1 <- data.frame(
  step = c("Discovery cohort", "Validation cohort", "Additional validation",
           "Gene sets", "Crosstalk modules"),
  n = c("GSE70947 (148 paired tumour/adjacent-normal)",
        "GSE70905 (47 tumour + 47 adjacent + 43 reduction)",
        "GSE15852 (43 paired)",
        "Ferroptosis 259; Necroptosis 158; overlap 18",
        "11 same-sign crosstalk modules (WGCNA)")
)
p1 <- ggplot(tab1, aes(x = step, y = 1, label = n)) +
  geom_tile(fill = "grey90", color = "white", height = 0.5) +
  geom_text(size = 4, hjust = 0) +
  coord_flip() + scale_x_discrete(limits = rev(tab1$step)) +
  labs(title = "Figure 1. Study overview") + theme_void()
save_fig("Fig1_study_overview", w = 10, h = 4)

## ---- Fig 2: QC (PCA + distribution) ----
ph <- ana$phenotype
pca_res <- prcomp(t(ana$expression_imputed), center = TRUE, scale. = FALSE)
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                     tissue = ph$tissue)
p2 <- ggplot(pca_df, aes(PC1, PC2, color = tissue)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("normal" = "#2e86ab", "tumor" = "#c0392b")) +
  labs(title = "Figure 2. Principal component analysis (GSE70947)") +
  theme_ms()
save_fig("Fig2_pca", w = 7, h = 6)

## ---- Fig 3: pathway scores by tissue (paired) ----
sc_long <- sc %>% select(sample_id, tissue, ferroptosis_full_singscore, necroptosis_full_singscore) %>%
  pivot_longer(c(ferroptosis_full_singscore, necroptosis_full_singscore),
               names_to = "pathway", values_to = "score")
sc_long$pathway <- ifelse(sc_long$pathway == "ferroptosis_full_singscore", "Ferroptosis", "Necroptosis")
p3 <- ggplot(sc_long, aes(x = tissue, y = score, fill = tissue)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_wrap(~pathway, scales = "free_y") +
  scale_fill_manual(values = c("normal" = "#2e86ab", "tumor" = "#c0392b")) +
  labs(title = "Figure 3. Pathway scores by tissue (GSE70947, paired)",
       x = NULL, y = "singscore") + theme_ms()
save_fig("Fig3_pathway_by_tissue", w = 8, h = 5)

## ---- Fig 4: within-tissue coupling scatter ----
sc_dot <- sc[sc$tissue == "tumor", ]
p4 <- ggplot(sc_dot, aes(necroptosis_full_singscore, ferroptosis_full_singscore)) +
  geom_point(alpha = 0.6, color = "#c0392b") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "Figure 4. Within-tumour ferroptosis-necroptosis coupling (GSE70947)",
       x = "Necroptosis singscore", y = "Ferroptosis singscore") + theme_ms()
save_fig("Fig4_within_tissue_coupling", w = 7, h = 6)

## ---- Fig 5: primary delta crosstalk ----
p5 <- ggplot(dt, aes(Deltanecroptosis_singscore, Deltaferroptosis_singscore)) +
  geom_point(alpha = 0.6, color = "#2e86ab") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "Figure 5. Patient-level delta crosstalk (GSE70947, n=148)",
       x = "Delta necroptosis (tumour - normal)", y = "Delta ferroptosis (tumour - normal)") + theme_ms()
save_fig("Fig5_delta_crosstalk", w = 7, h = 6)

## ---- Fig 6: WGCNA soft-threshold + module-trait heatmap ----
soft_plot <- soft[soft$power %in% c(1, 2, 4, 6, 8, 12, 16, 20, 24, 30), ]
p6a <- ggplot(soft_plot, aes(power, signed.R2)) +
  geom_line() + geom_point() +
  labs(x = "Soft threshold (power)", y = "Scale-free fit index R^2",
       title = "Figure 6A. Soft-threshold selection") + theme_ms()

mt_f <- mt[mt$trait %in% c("DeltaFerroptosis_singscore", "DeltaNecroptosis_singscore"), ]
mt_f$module <- factor(mt_f$module, levels = paste0("ME", 0:33))
p6b <- ggplot(mt_f, aes(x = module, y = trait, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  geom_text(aes(label = ifelse(bh_fdr < 0.05, "*", "")), size = 3) +
  scale_x_discrete(limits = paste0("ME", 0:33)) +
  labs(title = "Figure 6B. Module-trait correlations (singscore)",
       x = NULL, y = NULL, fill = "Pearson r") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_fig("Fig6_wgcna_soft_and_traits", w = 10, h = 6)

## ---- Fig 7: crosstalk module eigengene-delta associations ----
ct_mods <- c("ME3", "ME7", "ME9", "ME12", "ME17", "ME21", "ME22", "ME26", "ME28", "ME29", "ME30")
mod_f <- mt_f[mt_f$module %in% ct_mods, ]
mod_f$trait <- ifelse(mod_f$trait == "DeltaFerroptosis_singscore", "Delta ferroptosis", "Delta necroptosis")
p7 <- ggplot(mod_f, aes(x = module, y = correlation, fill = trait)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c("Delta ferroptosis" = "#7fbf7b", "Delta necroptosis" = "#af8dc3")) +
  labs(title = "Figure 7. Crosstalk module correlations with pathway delta scores",
       x = NULL, y = "Pearson r (module eigengene)") + theme_ms()
save_fig("Fig7_crosstalk_module_associations", w = 9, h = 5)

## ---- Fig 8: enrichment dotplot (KEGG, top terms across crosstalk modules) ----
kegg <- read.csv(file.path("results", "enrichment", "KEGG.tsv"), sep = "\t")
kegg$p.adjust <- as.numeric(kegg$p.adjust)
kegg_top <- kegg[kegg$module %in% ct_mods & kegg$p.adjust < 0.05, ]
kegg_top <- head(kegg_top[order(kegg_top$p.adjust), ], 15)
kegg_top$Description <- factor(kegg_top$Description, levels = rev(kegg_top$Description))
p8 <- ggplot(kegg_top, aes(x = p.adjust, y = Description, color = module)) +
  geom_point(aes(size = Count)) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Figure 8. Top KEGG pathways enriched in crosstalk modules",
       x = "Adjusted p-value (BH)", y = NULL, size = "Gene count") + theme_ms()
save_fig("Fig8_enrichment_kegg", w = 9, h = 6)

## ---- Fig 9: top hub genes ----
hub_ct <- hub[hub$crosstalk_module, ]
hub_top <- head(hub_ct[order(hub_ct$hub_score, decreasing = TRUE), ], 20)
hub_top$gene_symbol <- factor(hub_top$gene_symbol, levels = rev(hub_top$gene_symbol))
p9 <- ggplot(hub_top, aes(x = hub_score, y = gene_symbol, fill = module)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Figure 9. Top prioritized hub genes across crosstalk modules",
       x = "Hub evidence score", y = NULL) + theme_ms()
save_fig("Fig9_hub_genes", w = 8, h = 6)

## ---- Fig 10: external validation ----
val_cou$group <- c("Tumour" = "tumor", "Adjacent" = "adjacent_normal", "Reduction" = "reduction_normal")[val_cou$group]
p10a <- ggplot(val_cou, aes(x = group, y = r, fill = group)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Figure 10A. Within-tissue coupling in GSE70905",
       x = NULL, y = "Pearson r (ferroptosis vs necroptosis)") + theme_ms()

val_mod$sig <- ifelse(val_mod$fdr_ferro < 0.05 | val_mod$fdr_necro < 0.05, "sig", "ns")
p10b <- ggplot(val_mod, aes(x = module, y = r_delta_ferro, color = sig)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("sig" = "#c0392b", "ns" = "grey60")) +
  labs(title = "Figure 10B. Module delta-ferroptosis replication in GSE70905",
       x = NULL, y = "Validation r (module score vs Delta ferroptosis)") + theme_ms()
save_fig("Fig10_external_validation", w = 10, h = 6)

## ---- Supp A: leave-one-out ----
pA <- ggplot(loo, aes(x = r)) +
  geom_histogram(bins = 30, fill = "#2e86ab", color = "white") +
  geom_vline(xintercept = dt %>% {cor(.$Deltaferroptosis_singscore, .$Deltanecroptosis_singscore)}, color = "red", linetype = 2) +
  labs(title = "Supplementary A. Leave-one-out delta crosstalk correlations",
       x = "Pearson r", y = "Count") + theme_ms()
save_fig("SuppA_leave_one_out", w = 7, h = 5)

## ---- Supp B: field effect (group means) ----
fe <- read.csv(file.path("results", "field_effect", "group_means.tsv"), sep = "\t")
fe$pathway <- ifelse(fe$pathway == "ferro", "Ferroptosis", "Necroptosis")
fe$group <- c("reduction_normal" = "Reduction", "adjacent_normal" = "Adjacent", "tumor" = "Tumour")[fe$group]
fe$group <- factor(fe$group, levels = c("Reduction", "Adjacent", "Tumour"))
pB <- ggplot(fe, aes(x = group, y = mean, fill = group)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  facet_wrap(~pathway) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Supplementary B. Field-effect: pathway scores by tissue type (GSE70905)",
       x = NULL, y = "Mean singscore") + theme_ms()
save_fig("SuppB_field_effect", w = 8, h = 5)

## ---- Supp C: composition-adjusted coupling ----
comp_cou <- read.csv(file.path("results", "composition", "composition_adjusted_coupling.tsv"), sep = "\t")
comp_cou$comp_cou <- read.csv(file.path("results", "composition", "composition_adjusted_coupling.tsv"), sep = "\t")$adj_beta
pC <- ggplot(comp_cou, aes(x = model, y = base_beta, group = 1)) +
  geom_point(color = "grey60", size = 3) +
  geom_point(aes(y = adj_beta), color = "#c0392b", size = 3, shape = 17) +
  labs(title = "Supplementary C. Coupling before/after cell-type adjustment (tumour)",
       x = NULL, y = "Regression beta (ferroptosis ~ necroptosis)") + theme_ms()
save_fig("SuppC_composition_adjusted", w = 8, h = 5)

.log("Figures generated")
cat("DONE\n")