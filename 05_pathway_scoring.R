## 05_pathway_scoring.R
## Compute sample-level ferroptosis and necroptosis pathway scores.
## Two methods: singscore (primary, rank-based) and GSVA ssGSEA (sensitivity).
## Gene sets: full and overlap-removed (ferroptosis_unique, necroptosis_unique).
## Definitions are frozen here for both discovery and validation cohorts.
## Outputs:
##   data/processed/pathway_scores_all.rds
##   results/pathway_scores/GSE70947_pathway_scores.tsv
##   results/pathway_scores/GSE70905_pathway_scores.tsv
##   results/pathway_scores/GSE15852_pathway_scores.tsv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
score_dir <- .dir_for("results_scores")
gs_dir <- .dir_for("gene_sets")
data_proc <- .dir_for("data_processed")
set.seed(cfg$project$seed)

ferr <- readLines(file.path(gs_dir, "ferroptosis_genes_clean.txt"))
necr <- readLines(file.path(gs_dir, "necroptosis_genes_clean.txt"))
overlap <- readLines(file.path(gs_dir, "ferroptosis_necroptosis_overlap.txt"))
ferr_uniq <- setdiff(ferr, overlap)
necr_uniq <- setdiff(necr, overlap)

gene_sets <- list(
  ferroptosis_full = ferr,
  necroptosis_full = necr,
  ferroptosis_unique = ferr_uniq,
  necroptosis_unique = necr_uniq
)

score_one <- function(E_imp, gsets, gse, phenotype_extra = NULL) {
  # rank matrix (higher expression = higher rank) for singscore
  ranked <- rankGenes(E_imp)   # rows genes, cols samples; larger = higher
  res_list <- list()
  for (nm in names(gsets)) {
    gs <- gsets[[nm]]
    present <- intersect(gs, rownames(E_imp))
    cat(gse, nm, ": set size", length(gs), " measured", length(present), "\n")
    ss <- try(simpleScore(ranked, upSet = present, centerScore = TRUE), silent = TRUE)
    if (inherits(ss, "try-error")) {
      res_list[[paste0(nm, "_singscore")]] <- rep(NA_real_, ncol(E_imp))
    } else {
      res_list[[paste0(nm, "_singscore")]] <- as.numeric(ss$TotalScore)
    }
  }
  # GSVA ssgsea
  p <- ssgseaParam(exprData = E_imp, geneSets = gene_sets, normalize = TRUE)
  gsv <- gsva(p, verbose = FALSE)
  for (nm in names(gsets)) {
    res_list[[paste0(nm, "_gsva")]] <- as.numeric(gsv[nm, ])
  }
  df <- data.frame(sample_id = colnames(E_imp), do.call(cbind, res_list), stringsAsFactors = FALSE)
  df$gse <- gse
  if (!is.null(phenotype_extra)) {
    for (cn in names(phenotype_extra)) df[[cn]] <- phenotype_extra[[cn]]
  }
  df
}

g47 <- .load_analysis("GSE70947")
pn47 <- g47$phenotype
sid47 <- pn47$`sample id:ch1`
ph47 <- data.frame(
  tissue = pn47$`tissue:ch1`,
  patient_id = sub("\\.(N|T)\\.c2$", "", sid47),
  er = pn47$`er:ch1`,
  her2 = pn47$`her2:ch1`,
  age = suppressWarnings(as.numeric(pn47$`age_in_years:ch1`)),
  stringsAsFactors = FALSE
)
s47 <- score_one(g47$expression_imputed, gene_sets, "GSE70947", ph47)

g05 <- .load_analysis("GSE70905")
pn05 <- g05$phenotype
src05 <- pn05$source_name_ch1
grp05 <- ifelse(src05 == "breast adenocarcinoma", "tumor",
         ifelse(src05 == "breast tissue adjacent to tumor", "adjacent_normal",
         ifelse(src05 == "mammaplastic reduction", "reduction_normal", "other")))
ph05 <- data.frame(
  tissue = pn05$`tissue:ch1`,
  sample_type = pn05$sample_type,
  group = grp05,
  patient_id = sub("(N|T|B)$", "", pn05$description),
  er = pn05$`er:ch1`,
  her2 = pn05$`her2:ch1`,
  stringsAsFactors = FALSE
)
s05 <- score_one(g05$expression_imputed, gene_sets, "GSE70905", ph05)

g52 <- .load_analysis("GSE15852")
pn52 <- g52$phenotype
ph52 <- data.frame(
  tissue = pn52$`tissue:ch1`,
  sample_type = pn52$source_name_ch1,
  patient_id = pn52$description,
  stringsAsFactors = FALSE
)
s52 <- score_one(g52$expression_imputed, gene_sets, "GSE15852", ph52)

saveRDS(list(gene_sets = gene_sets, scores_GSE70947 = s47, scores_GSE70905 = s05, scores_GSE15852 = s52),
        file.path(data_proc, "pathway_scores_all.rds"))

.save_tsv(s47, "results_scores", "GSE70947_pathway_scores.tsv")
.save_tsv(s05, "results_scores", "GSE70905_pathway_scores.tsv")
.save_tsv(s52, "results_scores", "GSE15852_pathway_scores.tsv")
.log("Pathway scoring complete")
cat("DONE\n")