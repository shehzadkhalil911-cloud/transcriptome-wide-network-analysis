## 12_external_validation.R
## External validation of the discovery findings (GSE70947) in GSE70905
## (47 tumour + 47 adjacent-normal + 43 reduction-mammoplasty).
## Validates:
##   1. Cross-tissue pathway change (tumour vs adjacent normal, paired + unpaired)
##   2. Within-tissue ferroptosis/necroptosis coupling (normal vs tumour)
##   3. Crosstalk modules reconstructed from overlapping genes (gene-level module
##      scores), and association of module score with DeltaFerro/DeltaNecro
##   4. Hub gene replication: same-direction delta-correlation with pathway traits
## Outputs:
##   results/validation/validation_pathway_de.tsv
##   results/validation/validation_coupling.tsv
##   results/validation/validation_module_associations.tsv
##   results/validation/validation_hub_replication.tsv
##   results/validation/validation_summary.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
val_dir <- .dir_for("results_validation")

ana_val <- readRDS(file.path(.dir_for("data_processed"), "GSE70905_analysis.rds"))
expr_val <- ana_val$expression_imputed
ph_val <- ana_val$phenotype
stopifnot(all(colnames(expr_val) == rownames(ph_val)))

audit <- read.csv(file.path("results", "audit", "GSE70905_sample_audit.csv"), sep = "\t")
group <- setNames(audit$group, audit$sample_id)
ph_val$group <- group[rownames(ph_val)]

tumour_cols <- rownames(ph_val)[ph_val$group == "tumor"]
adjacent_cols <- rownames(ph_val)[ph_val$group == "adjacent_normal"]
reduction_cols <- rownames(ph_val)[ph_val$group == "reduction_normal"]
cat("tumour:", length(tumour_cols), " adjacent:", length(adjacent_cols), " reduction:", length(reduction_cols), "\n")

# patient pairing key from sample_id_raw (CM21N/CM21T -> CM21)
pkey <- function(sid) sub("(N|T|B)$", "", sid)
pair_t <- pkey(audit$sample_id_raw[match(tumour_cols, audit$sample_id)])
pair_a <- pkey(audit$sample_id_raw[match(adjacent_cols, audit$sample_id)])
pairs <- intersect(pair_t, pair_a)
cat("paired patients:", length(pairs), "\n")
# global maps: pair key -> tumour/adjacent sample id
names(pair_t) <- tumour_cols; names(pair_a) <- adjacent_cols
pt <- setNames(names(pair_t)[match(pairs, pair_t)], pairs)
pa <- setNames(names(pair_a)[match(pairs, pair_a)], pairs)
gs <- read.csv(file.path("data", "gene_sets", "ferroptosis_genes_clean.txt"), header = FALSE)$V1
ns <- read.csv(file.path("data", "gene_sets", "necroptosis_genes_clean.txt"), header = FALSE)$V1

score_matrix <- function(mat, gene_set, sig_sets = NULL) {
  common <- intersect(gene_set, rownames(mat))
  if (length(common) < 3) return(rep(NA, ncol(mat)))
  library(singscore)
  ranks <- rankGenes(mat)
  s <- simpleScore(ranks, upSet = common, centerScore = TRUE)
  as.numeric(s$TotalScore)
}
# match discovery score names
score_name_map <- list(
  ferro_full = "ferroptosis_full_singscore",
  necro_full = "necroptosis_full_singscore"
)

scores <- data.frame(
  sample_id = colnames(expr_val),
  stringsAsFactors = FALSE
)
scores$ferro <- score_matrix(expr_val, gs)
scores$necro <- score_matrix(expr_val, ns)
scores$group <- ph_val$group
scores$pair_key <- pkey(audit$sample_id_raw[match(scores$sample_id, audit$sample_id)])
rownames(scores) <- scores$sample_id

# ---- 1. cross-tissue pathway change (paired: tumour - adjacent) ----
ferro_by_id <- setNames(scores$ferro, scores$sample_id)
necro_by_id <- setNames(scores$necro, scores$sample_id)
de_out <- do.call(rbind, lapply(c("ferro", "necro"), function(pw) {
  vals <- if (pw == "ferro") ferro_by_id else necro_by_id
  diff <- vals[pt[pairs]] - vals[pa[pairs]]
  tt <- t.test(diff)
  data.frame(pathway = pw, n_pairs = length(diff), mean_diff = mean(diff),
             t_stat = tt$statistic, df = tt$parameter, p_value = tt$p.value)
}))
de_out$fdr <- p.adjust(de_out$p_value, method = "BH")
.save_tsv(de_out, "results_validation", "validation_pathway_de.tsv")

# ---- 2. within-tissue coupling ----
coupling <- do.call(rbind, lapply(c("tumor", "adjacent_normal", "reduction_normal"), function(grp) {
  sub <- scores[scores$group == grp, ]
  if (nrow(sub) < 5) return(NULL)
  ct <- cor.test(sub$ferro, sub$necro, method = "pearson")
  data.frame(group = grp, n = nrow(sub), r = ct$estimate, p_value = ct$p.value)
}))
coupling$fdr <- p.adjust(coupling$p_value, method = "BH")
.save_tsv(coupling, "results_validation", "validation_coupling.tsv")

# ---- 3. module reconstruction from overlapping genes ----
w <- readRDS(file.path(.dir_for("data_processed"), "wgcna_delta_result.rds"))
labels <- w$labels
module_of <- names(labels)[labels > 0]
labels2 <- labels[labels > 0]

ct_mods <- c("ME3", "ME7", "ME9", "ME12", "ME17", "ME21", "ME22", "ME26", "ME28", "ME29", "ME30")
mod_genes <- lapply(ct_mods, function(m) {
  lab <- as.integer(sub("ME", "", m))
  names(labels2)[labels2 == lab]
})
names(mod_genes) <- ct_mods

# Delta (tumour - adjacent) per gene for validation
common_genes <- intersect(rownames(expr_val), rownames(w$Delta))
expr_val_common <- expr_val[common_genes, ]
delta_val <- matrix(NA, nrow = length(common_genes), ncol = length(pairs))
rownames(delta_val) <- common_genes
colnames(delta_val) <- pairs
for (i in seq_along(pairs)) {
  tidx <- pt[[pairs[i]]]
  nidx <- pa[[pairs[i]]]
  delta_val[, i] <- expr_val_common[, tidx] - expr_val_common[, nidx]
}
# standardise per gene
delta_val_s <- t(scale(t(delta_val)))

# Delta traits from scores
delta_traits <- data.frame(
  pair = pairs,
  DeltaFerro = ferro_by_id[pt[pairs]] - ferro_by_id[pa[pairs]],
  DeltaNecro = necro_by_id[pt[pairs]] - necro_by_id[pa[pairs]]
)

mod_assoc <- do.call(rbind, lapply(ct_mods, function(m) {
  genes_in <- intersect(mod_genes[[m]], common_genes)
  genes_in <- genes_in[rowSums(is.na(delta_val_s[genes_in, , drop = FALSE])) == 0]
  if (length(genes_in) < 5) return(NULL)
  me <- colMeans(delta_val_s[genes_in, , drop = FALSE], na.rm = TRUE)
  ok <- is.finite(me) & is.finite(delta_traits$DeltaFerro) & is.finite(delta_traits$DeltaNecro)
  if (sum(ok) < 10) return(NULL)
  r1 <- cor(me[ok], delta_traits$DeltaFerro[ok])
  r2 <- cor(me[ok], delta_traits$DeltaNecro[ok])
  p1 <- cor.test(me[ok], delta_traits$DeltaFerro[ok])$p.value
  p2 <- cor.test(me[ok], delta_traits$DeltaNecro[ok])$p.value
  data.frame(module = m, n_overlapping_genes = length(genes_in),
             r_delta_ferro = r1, r_delta_necro = r2,
             p_ferro = p1, p_necro = p2)
}))
mod_assoc$fdr_ferro <- p.adjust(mod_assoc$p_ferro, method = "BH")
mod_assoc$fdr_necro <- p.adjust(mod_assoc$p_necro, method = "BH")
.save_tsv(mod_assoc, "results_validation", "validation_module_associations.tsv")

# ---- 4. hub replication ----
hub <- read.csv(file.path("results", "hubs", "hub_gene_evidence.tsv"), sep = "\t")
hub <- hub[hub$crosstalk_module, ]
hub$gene_symbol <- trimws(hub$gene_symbol)
hub_repl <- do.call(rbind, lapply(seq_len(nrow(hub)), function(i) {
  g <- hub$gene_symbol[i]
  if (!g %in% rownames(delta_val_s)) return(NULL)
  v <- delta_val_s[g, ]
  if (any(is.na(v))) return(NULL)
  r_ferro <- cor(v, delta_traits$DeltaFerro, use = "pairwise.complete.obs")
  r_necro <- cor(v, delta_traits$DeltaNecro, use = "pairwise.complete.obs")
  data.frame(gene = g, module = hub$module[i],
             discovery_r_ferro = hub$GS_ferro[i],
             discovery_r_necro = hub$GS_necro[i],
             validation_r_ferro = r_ferro,
             validation_r_necro = r_necro)
}))
hub_repl$same_dir_ferro <- sign(hub_repl$discovery_r_ferro) == sign(hub_repl$validation_r_ferro)
hub_repl$same_dir_necro <- sign(hub_repl$discovery_r_necro) == sign(hub_repl$validation_r_necro)
hub_repl$replicated <- hub_repl$same_dir_ferro & hub_repl$same_dir_necro
.save_tsv(hub_repl, "results_validation", "validation_hub_replication.tsv")

summ <- data.frame(
  n_tumour = length(tumour_cols), n_adjacent = length(adjacent_cols),
  n_reduction = length(reduction_cols), n_paired = length(pairs),
  ferro_paired_mean_diff = de_out$mean_diff[1], ferro_paired_p = de_out$p_value[1],
  necro_paired_mean_diff = de_out$mean_diff[2], necro_paired_p = de_out$p_value[2],
  coupling_tumor_r = coupling$r[coupling$group == "tumor"],
  coupling_tumor_p = coupling$p_value[coupling$group == "tumor"],
  n_reconstructed_modules = nrow(mod_assoc),
  n_sig_module_delta_ferro = sum(mod_assoc$fdr_ferro < 0.05, na.rm = TRUE),
  n_sig_module_delta_necro = sum(mod_assoc$fdr_necro < 0.05, na.rm = TRUE),
  n_hubs_tested = nrow(hub_repl),
  frac_hubs_replicated = mean(hub_repl$replicated, na.rm = TRUE)
)
.save_tsv(summ, "results_validation", "validation_summary.csv")

.log("External validation complete")
cat("DONE\n")