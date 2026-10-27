## 08_wgcna_delta_expression.R
## WGCNA on the patient-level tumour-minus-normal delta expression matrix.
## Identifies modules of co-varying tumour-associated gene changes and tests
## their association with DeltaFerroptosis and DeltaNecroptosis.
## Outputs:
##   results/wgcna/delta_expression.rds
##   results/wgcna/soft_threshold_results.tsv
##   results/wgcna/module_assignments.tsv
##   results/wgcna/module_eigengenes.tsv
##   results/wgcna/module_trait_correlations.tsv
##   results/wgcna/module_gs_kme.tsv
##   results/wgcna/wgcna_summary.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
wgcna_dir <- .dir_for("results_wgcna")
data_proc <- .dir_for("data_processed")
set.seed(cfg$project$seed)
allowWGCNAThreads(nThreads = 2)

obj <- .load_analysis("GSE70947")
E_imp <- obj$expression_imputed
pn <- obj$phenotype
sid <- pn$`sample id:ch1`
patient <- factor(sub("\\.(N|T)\\.c2$", "", sid))
tissue <- pn$`tissue:ch1`

# patient-level delta expression (align patients by ID)
tum_idx <- which(tissue == "tumor"); nor_idx <- which(tissue == "normal")
ord_t <- tum_idx[order(patient[tum_idx])]; ord_n <- nor_idx[order(patient[nor_idx])]
stopifnot(all(patient[ord_t] == patient[ord_n]))
tum <- E_imp[, ord_t, drop = FALSE]
nor <- E_imp[, ord_n, drop = FALSE]
Delta <- tum - nor
colnames(Delta) <- as.character(patient[ord_t])
cat("Delta matrix dim:", paste(dim(Delta), collapse = "x"), "\n")

# remove genes with zero variance
v <- apply(Delta, 1, var)
Delta <- Delta[v > 0, ]
cat("genes after zero-variance removal:", nrow(Delta), "\n")

# variance filter: retain top (1 - variance_filter_quantile) fraction by variance
v <- apply(Delta, 1, var)
vfq <- cfg$analysis$wgcna$variance_filter_quantile
thr <- quantile(v, probs = vfq)
Delta <- Delta[v >= thr, ]
cat("genes retained after variance filter (quantile", vfq, "):", nrow(Delta), "\n")
v <- apply(Delta, 1, var)

# sample clustering for outlier report
hc <- hclust(dist(t(Delta)), method = "average")
png(file.path(wgcna_dir, "delta_sample_hclust.png"), width = 1000, height = 600, res = 130)
plot(hc, labels = FALSE, main = "Patient delta-expression clustering")
dev.off()

# soft threshold (on a deterministic subsample of 5000 genes for tractability)
powers <- c(seq(1, 10), seq(12, 30, by = 2))
set.seed(cfg$project$seed)
sub_genes <- sample(rownames(Delta), size = min(5000, nrow(Delta)))
cat("picking soft threshold on", length(sub_genes), "genes...\n")
sft <- pickSoftThreshold(t(Delta[sub_genes, ]), powerVector = powers, networkType = "signed", verbose = 0)
sft_tab <- data.frame(
  soft_power = sft$fitIndices$Power,
  scale_free_R2 = sft$fitIndices$SFT.R.sq,
  mean_connectivity = sft$fitIndices$mean.k.,
  median_connectivity = sft$fitIndices$median.k.,
  stringsAsFactors = FALSE
)
.save_tsv(sft_tab, "results_wgcna", "soft_threshold_results.tsv")

# choose power: first with R2 >= 0.85, else the one maximizing R2
ok <- sft_tab$scale_free_R2 >= 0.85
if (any(ok)) {
  power_sel <- sft_tab$soft_power[which(ok)[1]]
} else {
  power_sel <- sft_tab$soft_power[which.max(sft_tab$scale_free_R2)]
}
cat("selected soft power:", power_sel, " (R2 =", sft_tab$scale_free_R2[sft_tab$soft_power == power_sel], ")\n")

# blockwise network construction
mms <- cfg$analysis$wgcna$min_module_size
mch <- cfg$analysis$wgcna$merge_cut_height
net <- blockwiseModules(
  t(Delta), power = power_sel, TOMType = "signed", networkType = "signed",
  minModuleSize = mms, mergeCutHeight = mch, numericLabels = TRUE,
  maxBlockSize = 8000, verbose = 1, saveTOMs = FALSE,
  randomSeed = cfg$project$seed
)
labels <- net$colors
names(labels) <- rownames(Delta)
n_modules <- length(setdiff(unique(labels), 0))
cat("n modules (excl grey=0):", n_modules, "\n")
print(table(labels))

# module eigengenes
ME <- moduleEigengenes(t(Delta), colors = labels)$eigengenes
names(ME) <- gsub("ME", "ME", names(ME))

# module assignments
assign_df <- data.frame(
  gene_symbol = names(labels),
  module_label = as.integer(labels),
  module_color = labels2colors(labels),
  stringsAsFactors = FALSE
)
.save_tsv(assign_df, "results_wgcna", "module_assignments.tsv")

# traits
delta_tab <- read.csv(file.path("results", "crosstalk", "patient_delta_scores.tsv"),
                      sep = "\t", check.names = FALSE)
traits <- data.frame(
  patient_id = delta_tab$patient_id,
  DeltaFerroptosis_singscore = delta_tab$Deltaferroptosis_singscore,
  DeltaNecroptosis_singscore = delta_tab$Deltanecroptosis_singscore,
  DeltaFerroptosis_gsva = delta_tab$Deltaferroptosis_gsva,
  DeltaNecroptosis_gsva = delta_tab$Deltanecroptosis_gsva,
  stringsAsFactors = FALSE
)
rownames(traits) <- traits$patient_id
# align to Delta matrix columns
traits <- traits[colnames(Delta), ]

# module-trait correlations (Pearson) with FDR across tests
mod_names <- colnames(ME)
trait_cols <- setdiff(colnames(traits), "patient_id")
m_tab <- list()
for (mod in mod_names) {
  for (tr in trait_cols) {
    ct <- suppressWarnings(cor.test(ME[[mod]], traits[[tr]], method = "pearson"))
    m_tab[[length(m_tab) + 1]] <- data.frame(
      module = mod, trait = tr, correlation = ct$estimate, p_value = ct$p.value,
      n = length(traits[[tr]]), stringsAsFactors = FALSE
    )
  }
}
mt <- do.call(rbind, m_tab)
mt$bh_fdr <- p.adjust(mt$p_value, method = "BH")
# flag crosstalk modules: FDR<0.05 for both DeltaFerroptosis and DeltaNecroptosis (singscore primary)
prim <- mt[mt$trait %in% c("DeltaFerroptosis_singscore", "DeltaNecroptosis_singscore"), ]
prim_wide <- reshape(prim[, c("module", "trait", "bh_fdr")], idvar = "module", timevar = "trait", direction = "wide")
crosstalk_modules <- prim_wide$module[
  !is.na(prim_wide[[paste0("bh_fdr.", "DeltaFerroptosis_singscore")]]) &
    prim_wide[[paste0("bh_fdr.", "DeltaFerroptosis_singscore")]] < cfg$analysis$wgcna$module_trait_fdr &
    !is.na(prim_wide[[paste0("bh_fdr.", "DeltaNecroptosis_singscore")]]) &
    prim_wide[[paste0("bh_fdr.", "DeltaNecroptosis_singscore")]] < cfg$analysis$wgcna$module_trait_fdr
]
cat("Crosstalk modules (FDR<0.05 for both DeltaFerro & DeltaNecro, singscore):",
    paste(crosstalk_modules, collapse = ", "), "\n")

.save_tsv(mt, "results_wgcna", "module_trait_correlations.csv")
write.table(ME, file.path(wgcna_dir, "module_eigengenes.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)

# gene significance and module membership (kME) for DeltaFerroptosis/DeltaNecroptosis singscore
GS <- cor(t(Delta), traits$DeltaFerroptosis_singscore, use = "pairwise.complete.obs")
GSn <- cor(t(Delta), traits$DeltaNecroptosis_singscore, use = "pairwise.complete.obs")
kME <- signedKME(t(Delta), ME, outputColumnName = "kME")
gs_kme <- data.frame(
  gene_symbol = names(labels),
  module_label = as.integer(labels),
  module_color = labels2colors(labels),
  GS_DeltaFerroptosis = GS,
  GS_DeltaNecroptosis = GSn,
  stringsAsFactors = FALSE
)
gs_kme <- cbind(gs_kme, as.data.frame(kME))

# intramodular connectivity (within-module signed adjacency sum)
ik <- sapply(setdiff(unique(labels), 0), function(mod) {
  g <- names(labels)[labels == mod]
  subM <- t(Delta[g, ])
  adj <- adjacency(subM, power = power_sel, type = "signed")
  rowSums(adj)
})
kWithin <- rep(NA_real_, length(labels)); names(kWithin) <- names(labels)
for (mod in names(ik)) kWithin[names(ik[[mod]])] <- ik[[mod]]
gs_kme$kWithin <- kWithin[gs_kme$gene_symbol]

.save_tsv(gs_kme, "results_wgcna", "module_gs_kme.tsv")

summary <- data.frame(
  genes_input = nrow(Delta),
  n_patients = ncol(Delta),
  soft_power_selected = power_sel,
  scale_free_R2_selected = sft_tab$scale_free_R2[sft_tab$soft_power == power_sel],
  n_modules = n_modules,
  min_module_size = mms,
  merge_cut_height = mch,
  network_type = "signed",
  TOM_type = "signed",
  variance_filter_quantile = vfq,
  crosstalk_modules = paste(crosstalk_modules, collapse = ";"),
  stringsAsFactors = FALSE
)
.save_tsv(summary, "results_wgcna", "wgcna_summary.csv")
.save_tsv(assign_df, "results_wgcna", "module_assignments_export.tsv")

# save objects for downstream
saveRDS(list(Delta = Delta, labels = labels, ME = ME, traits = traits,
             power = power_sel, net = net, module_assignments = assign_df,
             crosstalk_modules = crosstalk_modules, gs_kme = gs_kme),
        file.path(data_proc, "wgcna_delta_result.rds"))
.log("WGCNA delta-expression analysis complete")
cat("DONE\n")