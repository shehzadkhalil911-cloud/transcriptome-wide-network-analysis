## 04_differential_expression.R
## Paired tumour vs matched normal differential expression using limma
## with an explicit patient design for GSE70947 (discovery).
## Outputs:
##   results/differential_expression/GSE70947_paired_limma_all_genes.tsv
##   results/differential_expression/ferroptosis_DE.tsv
##   results/differential_expression/necroptosis_DE.tsv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
de_dir <- .dir_for("results_de")
gs_dir <- .dir_for("gene_sets")
set.seed(cfg$project$seed)

obj <- .load_analysis("GSE70947")
E <- obj$expression          # gene-level, with NAs (limma handles)
pn <- obj$phenotype
sid <- pn$`sample id:ch1`
patient <- factor(sub("\\.(N|T)\\.c2$", "", sid))
tissue <- factor(pn$`tissue:ch1`, levels = c("normal", "tumor"))
stopifnot(all(colnames(E) == pn$geo_accession))

design <- model.matrix(~ patient + tissue)
colnames(design) <- make.names(colnames(design))
coef_interest <- "tissuetumor"
cat("design dim:", paste(dim(design), collapse = "x"), "\n")

fit <- lmFit(E, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = coef_interest, number = Inf, sort.by = "none")

# compute standard error of the coefficient of interest
se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_interest]
names(se) <- rownames(E)

res <- data.frame(
  gene_symbol = rownames(E),
  log2FC = tt$logFC,
  average_expression = tt$AveExpr,
  standard_error = unname(se[rownames(E)]),
  t_statistic = tt$t,
  raw_pvalue = tt$P.Value,
  bh_fdr = tt$adj.P.Val,
  stringsAsFactors = FALSE
)
# sanity: B statistic not needed
cat("n genes tested:", nrow(res), "\n")
cat("n DEG (FDR<0.05):", sum(res$bh_fdr < 0.05), "\n")
cat("n up (FDR<0.05):", sum(res$bh_fdr < 0.05 & res$log2FC > 0), "\n")
cat("n down (FDR<0.05):", sum(res$bh_fdr < 0.05 & res$log2FC < 0), "\n")

# save full table
.save_tsv(res, "results_de", "GSE70947_paired_limma_all_genes.tsv")

# pathway-specific tables
ferr <- readLines(file.path(gs_dir, "ferroptosis_genes_clean.txt"))
necr <- readLines(file.path(gs_dir, "necroptosis_genes_clean.txt"))
overlap <- readLines(file.path(gs_dir, "ferroptosis_necroptosis_overlap.txt"))

res$ferroptosis_member <- res$gene_symbol %in% ferr
res$necroptosis_member <- res$gene_symbol %in% necr
res$both_pathway_member <- res$gene_symbol %in% overlap

ferr_de <- res[res$ferroptosis_member, c("gene_symbol","log2FC","average_expression","standard_error","t_statistic","raw_pvalue","bh_fdr","necroptosis_member")]
necr_de <- res[res$necroptosis_member, c("gene_symbol","log2FC","average_expression","standard_error","t_statistic","raw_pvalue","bh_fdr","ferroptosis_member")]

cat("ferroptosis genes present in expression data:", nrow(ferr_de), "/", length(ferr), "\n")
cat("necroptosis genes present in expression data:", nrow(necr_de), "/", length(necr), "\n")
cat("ferroptosis genes DEG (FDR<0.05):", sum(ferr_de$bh_fdr < 0.05), "\n")
cat("necroptosis genes DEG (FDR<0.05):", sum(necr_de$bh_fdr < 0.05), "\n")

.save_tsv(ferr_de, "results_de", "ferroptosis_DE.tsv")
.save_tsv(necr_de, "results_de", "necroptosis_DE.tsv")

# Also save the full DE table with membership appended for reuse
.save_tsv(res, "results_de", "GSE70947_paired_limma_all_genes_with_membership.tsv")

# summary of pathway DE for downstream use
summary_de <- data.frame(
  pathway = c("ferroptosis", "necroptosis"),
  genes_in_set = c(length(ferr), length(necr)),
  genes_measured = c(nrow(ferr_de), nrow(necr_de)),
  deg_fdr005_up = c(sum(ferr_de$bh_fdr < 0.05 & ferr_de$log2FC > 0), sum(necr_de$bh_fdr < 0.05 & necr_de$log2FC > 0)),
  deg_fdr005_down = c(sum(ferr_de$bh_fdr < 0.05 & ferr_de$log2FC < 0), sum(necr_de$bh_fdr < 0.05 & necr_de$log2FC < 0)),
  stringsAsFactors = FALSE
)
.save_tsv(summary_de, "results_de", "pathway_de_summary.csv")
.log("Differential expression complete")
cat("DONE\n")