## 10_hub_prioritisation.R
## Multi-evidence hub prioritisation for WGCNA crosstalk modules.
## Evidence per gene:
##   1) kME / kWithin (intramodular hubness)
##   2) Gene significance vs DeltaFerroptosis and DeltaNecroptosis (singscore)
##   3) Differential expression (paired limma |logFC| and -log10 FDR)
##   4) Delta-correlation with each pathway trait (Pearson r to trait)
##   5) Membership in ferroptosis/necroptosis reference gene sets (bonus)
## Outputs:
##   results/hubs/hub_gene_evidence.tsv  - full evidence table (all analysed genes)
##   results/hubs/hub_gene_ranking.tsv   - composite-ranked genes within each crosstalk module
##   results/hubs/top_hubs_per_module.tsv - top-k hubs per module
##   results/hubs/hub_summary.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
hub_dir <- .dir_for("results_hubs")

gs <- read.csv(file.path("data", "gene_sets", "ferroptosis_genes_clean.txt"), header = FALSE)$V1
ns <- read.csv(file.path("data", "gene_sets", "necroptosis_genes_clean.txt"), header = FALSE)$V1
ferro_set <- toupper(gs); necro_set <- toupper(ns)

w <- readRDS(file.path(.dir_for("data_processed"), "wgcna_delta_result.rds"))
labels <- w$labels
Delta <- w$Delta
trait_df <- read.csv(file.path("results", "crosstalk", "patient_delta_scores.tsv"), sep = "\t")
trait_df <- trait_df[order(match(trait_df$patient, colnames(Delta))), ]

singscore_traits <- trait_df[, c("Deltaferroptosis_singscore", "Deltanecroptosis_singscore")]
stopifnot(nrow(singscore_traits) == ncol(Delta))

gsk <- read.csv(file.path("results", "wgcna", "module_gs_kme.tsv"), sep = "\t")
gsk <- gsk[!is.na(gsk$module_label) & gsk$module_label > 0, ]
gsk$kWithin <- sapply(seq_len(nrow(gsk)), function(i) {
  col <- paste0("kME", gsk$module_label[i])
  gsk[[col]][i]
})
gsk <- gsk[!is.na(gsk$kWithin), ]

de <- read.csv(file.path("results", "differential_expression", "GSE70947_paired_limma_all_genes.tsv"), sep = "\t")
de$abs_logFC <- abs(de$log2FC)
de$neg_log10_fdr <- -log10(pmax(de$bh_fdr, 1e-300))

# gene-level delta-trait correlations
genes <- rownames(Delta)
delta_corr <- lapply(c("Deltaferroptosis_singscore", "Deltanecroptosis_singscore"), function(tname) {
  tv <- singscore_traits[[tname]]
  rc <- vapply(seq_len(nrow(Delta)), function(i) suppressWarnings(cor(Delta[i, ], tv, method = "pearson")), numeric(1))
  names(rc) <- genes
  rc
})
names(delta_corr) <- c("r_delta_ferro", "r_delta_necro")

ev <- data.frame(
  gene_symbol = gsk$gene_symbol,
  module_label = gsk$module_label,
  module = paste0("ME", gsk$module_label),
  module_color = gsk$module_color,
  kWithin = gsk$kWithin,
  kME_own = gsk$kWithin,
  GS_ferro = gsk$GS_DeltaFerroptosis,
  GS_necro = gsk$GS_DeltaNecroptosis,
  stringsAsFactors = FALSE
)
ev$r_delta_ferro <- delta_corr$r_delta_ferro[ev$gene_symbol]
ev$r_delta_necro <- delta_corr$r_delta_necro[ev$gene_symbol]
ev <- merge(ev, de[, c("gene_symbol", "log2FC", "abs_logFC", "neg_log10_fdr")], by = "gene_symbol", all.x = TRUE)
ev$ferro_member <- ev$gene_symbol %in% ferro_set
ev$necro_member <- ev$gene_symbol %in% necro_set
ev$overlap_member <- ev$gene_symbol %in% intersect(ferro_set, necro_set)

rank_scale <- function(x) {
  r <- rank(x, na.last = "keep")
  r / max(r, na.rm = TRUE)
}

ev$s_kme <- rank_scale(ev$kWithin)
ev$s_gs_ferro <- rank_scale(abs(ev$GS_ferro))
ev$s_gs_necro <- rank_scale(abs(ev$GS_necro))
ev$s_de <- rank_scale(ev$neg_log10_fdr)
ev$s_dr_ferro <- rank_scale(abs(ev$r_delta_ferro))
ev$s_dr_necro <- rank_scale(abs(ev$r_delta_necro))

ev$hub_score <- ev$s_kme + ev$s_gs_ferro + ev$s_gs_necro + ev$s_de +
  ev$s_dr_ferro + ev$s_dr_necro
ev$hub_score_ferro <- ev$s_kme + ev$s_gs_ferro + ev$s_dr_ferro + ev$s_de
ev$hub_score_necro <- ev$s_kme + ev$s_gs_necro + ev$s_dr_necro + ev$s_de

ct_mods <- c("ME3", "ME7", "ME9", "ME12", "ME17", "ME21", "ME22", "ME26", "ME28", "ME29", "ME30")
ev$crosstalk_module <- ev$module %in% ct_mods

.save_tsv(ev, "results_hubs", "hub_gene_evidence.tsv")

ranked <- do.call(rbind, lapply(ct_mods, function(mod) {
  sub <- ev[ev$module == mod, ]
  sub <- sub[order(sub$hub_score, decreasing = TRUE), ]
  sub$module_rank <- seq_len(nrow(sub))
  sub$top10 <- sub$module_rank <= 10
  sub
}))
.save_tsv(ranked, "results_hubs", "hub_gene_ranking.tsv")

top <- do.call(rbind, lapply(ct_mods, function(mod) {
  sub <- ranked[ranked$module == mod & ranked$top10, ]
  head(sub[order(sub$hub_score, decreasing = TRUE), ], 10)
}))
.save_tsv(top, "results_hubs", "top_hubs_per_module.tsv")

# module-level summary
summ <- do.call(rbind, lapply(ct_mods, function(mod) {
  sub <- ranked[ranked$module == mod, ]
  data.frame(
    module = mod,
    n_genes = nrow(sub),
    n_ferro_members = sum(sub$ferro_member),
    n_necro_members = sum(sub$necro_member),
    n_overlap_members = sum(sub$overlap_member),
    mean_kWithin = round(mean(sub$kWithin, na.rm = TRUE), 3),
    mean_abs_gs_ferro = round(mean(abs(sub$GS_ferro), na.rm = TRUE), 3),
    mean_abs_gs_necro = round(mean(abs(sub$GS_necro), na.rm = TRUE), 3)
  )
}))
.save_tsv(summ, "results_hubs", "hub_summary.csv")

.log("Hub prioritisation complete")
cat("DONE\n")