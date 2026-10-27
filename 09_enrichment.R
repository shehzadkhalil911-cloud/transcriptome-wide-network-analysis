## 09_enrichment.R
## Functional enrichment (GO BP, KEGG, Reactome, Hallmark) for WGCNA crosstalk
## modules and for the prioritized hub candidates.
## Universe = all measured genes in the discovery analysis.
## Outputs:
##   results/enrichment/GO_BP.tsv, KEGG.tsv, Reactome.tsv, Hallmark.tsv
##   results/enrichment/enrichment_summary.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
en_dir <- .dir_for("results_enrichment")
data_proc <- .dir_for("data_processed")
set.seed(cfg$project$seed)

w <- readRDS(file.path(data_proc, "wgcna_delta_result.rds"))
labels <- w$labels
Delta <- w$Delta
universe <- rownames(Delta)

# same-sign crosstalk modules: FDR<0.05 for both singscore traits, same correlation sign
mt <- read.csv(file.path("results", "wgcna", "module_trait_correlations.csv"), sep = "\t")
prim <- mt[mt$trait %in% c("DeltaFerroptosis_singscore", "DeltaNecroptosis_singscore"), ]
mods <- unique(prim$module)
ct <- list()
for (mod in mods) {
  r1 <- prim$correlation[prim$module == mod & prim$trait == "DeltaFerroptosis_singscore"]
  r2 <- prim$correlation[prim$module == mod & prim$trait == "DeltaNecroptosis_singscore"]
  f1 <- prim$bh_fdr[prim$module == mod & prim$trait == "DeltaFerroptosis_singscore"]
  f2 <- prim$bh_fdr[prim$module == mod & prim$trait == "DeltaNecroptosis_singscore"]
  if (f1 < 0.05 && f2 < 0.05 && sign(r1) == sign(r2)) {
    ct[[mod]] <- data.frame(module = mod, r_ferro = r1, r_necro = r2,
                            fdr_ferro = f1, fdr_necro = f2, sign = sign(r1))
  }
}
ct_tab <- do.call(rbind, ct)
ct_tab$module_label <- as.integer(sub("ME", "", ct_tab$module))
ct_tab$n_genes <- sapply(ct_tab$module_label, function(l) sum(labels == l))
cat("Same-sign crosstalk modules:\n"); print(ct_tab)
.save_tsv(ct_tab, "results_enrichment", "crosstalk_modules_defined.csv")

module_gene_lists <- lapply(ct_tab$module_label, function(l) names(labels)[labels == l])
names(module_gene_lists) <- ct_tab$module

# Restrict universe to genes annotated in org.Hs.eg.db (natural GO/Reactome universe)
keymap <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "SYMBOL", columns = "ENTREZID")
keymap <- keymap[!is.na(keymap$ENTREZID), ]
universe_annotated <- unique(keymap$SYMBOL)
cat("universe:", length(universe), " annotated:", length(universe_annotated), "\n")
flush.console()

# ---- GO BP ----
run_go <- function(gene_vector, universe, label) {
  ego <- try(enrichGO(gene = gene_vector, universe = universe, OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                      pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE), silent = TRUE)
  if (inherits(ego, "try-error")) return(NULL)
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  res <- as.data.frame(ego)
  res$module <- label
  res
}

all_go <- do.call(rbind, lapply(names(module_gene_lists), function(mn) {
  cat("GO BP module", mn, "...\n"); flush.console()
  run_go(module_gene_lists[[mn]], universe_annotated, mn)
}))
if (!is.null(all_go)) {
  all_go$gene_ratio <- all_go$GeneRatio
  all_go$bg_ratio <- all_go$BgRatio
  .save_tsv(all_go, "results_enrichment", "GO_BP.tsv")
} else {
  cat("GO BP: no enriched terms\n")
  .save_tsv(data.frame(), "results_enrichment", "GO_BP.tsv")
}

# ---- KEGG ----
run_kegg <- function(gene_vector, label) {
  kid <- try(bitr(gene_vector, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db), silent = TRUE)
  if (inherits(kid, "try-error") || nrow(kid) == 0) return(NULL)
  ek <- try(enrichKEGG(gene = unique(kid$ENTREZID), organism = "hsa",
                       pvalueCutoff = 1, qvalueCutoff = 1), silent = TRUE)
  if (inherits(ek, "try-error")) return(NULL)
  if (is.null(ek) || nrow(ek) == 0) return(NULL)
  res <- as.data.frame(ek)
  res$module <- label
  res
}
all_kegg <- do.call(rbind, lapply(names(module_gene_lists), function(mn) {
  cat("KEGG module", mn, "...\n"); flush.console()
  run_kegg(module_gene_lists[[mn]], mn)
}))
if (!is.null(all_kegg)) .save_tsv(all_kegg, "results_enrichment", "KEGG.tsv") else {
  cat("KEGG: none\n"); .save_tsv(data.frame(), "results_enrichment", "KEGG.tsv")
}

# ---- Reactome ----
suppressPackageStartupMessages(library(ReactomePA))
run_react <- function(gene_vector, label) {
  kid <- try(bitr(gene_vector, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db), silent = TRUE)
  if (inherits(kid, "try-error") || nrow(kid) == 0) return(NULL)
  er <- try(enrichPathway(gene = unique(kid$ENTREZID), organism = "human",
                          pvalueCutoff = 1, qvalueCutoff = 1), silent = TRUE)
  if (inherits(er, "try-error")) return(NULL)
  if (is.null(er) || nrow(er) == 0) return(NULL)
  res <- as.data.frame(er)
  res$module <- label
  res
}
all_react <- do.call(rbind, lapply(names(module_gene_lists), function(mn) {
  cat("Reactome module", mn, "...\n"); flush.console()
  run_react(module_gene_lists[[mn]], mn)
}))
if (!is.null(all_react)) .save_tsv(all_react, "results_enrichment", "Reactome.tsv") else {
  cat("Reactome: none\n"); .save_tsv(data.frame(), "results_enrichment", "Reactome.tsv")
}

# ---- Hallmark ----
hall <- msigdbr(species = "Homo sapiens", collection = "H")
hall_gs <- split(hall$gene_symbol, hall$gs_name)
ter <- data.frame(term = rep(names(hall_gs), lengths(hall_gs)), gene = unlist(hall_gs))
all_hall <- do.call(rbind, lapply(names(module_gene_lists), function(mn) {
  gene_vector <- module_gene_lists[[mn]]
  eh <- try(enricher(gene = gene_vector, TERM2GENE = ter, universe = universe,
                     pvalueCutoff = 1, qvalueCutoff = 1), silent = TRUE)
  if (inherits(eh, "try-error") || is.null(eh) || nrow(eh) == 0) return(NULL)
  res <- as.data.frame(eh); res$module <- mn; res
}))
if (!is.null(all_hall)) .save_tsv(all_hall, "results_enrichment", "Hallmark.tsv") else {
  cat("Hallmark: none\n"); .save_tsv(data.frame(), "results_enrichment", "Hallmark.tsv")
}

.log("Enrichment complete")
cat("DONE\n")