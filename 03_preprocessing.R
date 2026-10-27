## 03_preprocessing.R
## Build analysis-ready gene-level expression matrices and run QC.
## Uses the gene-level log2 matrices already present in the RDS objects.
## Control probes are removed; genes with > max_gene_na_fraction NAs are dropped;
## remaining NAs are imputed (median per gene) for downstream complete-matrix
## analyses (scoring, WGCNA). Differential expression uses limma which handles
## missing values per gene.
## Outputs:
##   data/processed/GSE70947_analysis.rds (and GSE70905, GSE15852)
##   results/qc/* (plots + tables)

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
qc_dir <- .dir_for("results_qc")
data_proc <- .dir_for("data_processed")
impute_method <- cfg$qc$imputation
max_na_frac <- cfg$qc$max_gene_na_fraction
seed <- cfg$project$seed
set.seed(seed)

build_control_list <- function(x) {
  # from annotation: probes flagged as controls map to gene symbols
  at <- x$annotation
  gm <- x$gene_map
  ctrl_probe <- at$ProbeName[at$ControlType != 0]
  ctrl_idx <- which(at$ControlType != 0)
  ctrl_genes <- unique(gm[ctrl_idx])
  ctrl_genes <- ctrl_genes[!is.na(ctrl_genes) & ctrl_genes != ""]
  c(ctrl_genes, "GE_BrightCorner", "DarkCorner")
}

process_gse <- function(gse) {
  .log("Processing ", gse)
  x <- readRDS(file.path(cfg$project$dirs$data_raw, paste0(gse, "_processed.rds")))
  E <- x$expression
  pn <- x$phenotype

  # remove control genes
  ctrl <- build_control_list(x)
  keep <- !(rownames(E) %in% ctrl)
  # also remove obvious spike-in / QC probe names at gene level
  keep <- keep & !grepl("^(ERCC-|E1A|r60|r80|r10|e1_|TM_|GC_|RIKEN)", rownames(E))
  E <- E[keep, , drop = FALSE]
  cat(gse, "genes after control removal:", nrow(E), "\n")

  # drop genes with high NA fraction
  na_frac <- rowMeans(is.na(E))
  keep2 <- na_frac <= max_na_frac
  cat(gse, "genes kept after NA filter (NA<=", max_na_frac, "):", sum(keep2), " dropped:", sum(!keep2), "\n")
  E <- E[keep2, , drop = FALSE]

  # impute remaining NAs
  E_imp <- impute_matrix(E, method = impute_method, seed = seed)
  cat(gse, "imputed NAs; remaining NA count:", sum(is.na(E_imp)), "\n")

  obj <- list(
    gse = gse,
    expression = E,             # with NAs (for limma)
    expression_imputed = E_imp, # complete (for scoring / WGCNA)
    phenotype = pn,
    platform = x$platform,
    gene_map = x$gene_map,
    scale_info = x$scale_info
  )

  # ---- QC ----
  png(file.path(qc_dir, paste0(gse, "_distributions.png")), width = 1200, height = 700, res = 130)
  boxplot(E_imp[, seq(1, ncol(E_imp), by = max(1, floor(ncol(E_imp) / 60)))],
          outline = FALSE, las = 2, cex.axis = 0.5,
          main = paste(gse, "sample expression distributions"),
          ylab = "log2 expression")
  dev.off()

  # sample-wise QC stats
  qc <- data.frame(
    sample_id = colnames(E_imp),
    mean_expr = colMeans(E_imp),
    sd_expr = apply(E_imp, 2, sd),
    na_fraction_before_impute = colMeans(is.na(E)),
    stringsAsFactors = FALSE
  )

  # inter-array connectivity (outlier proxy) on top 5000 most variable genes
  v <- apply(E_imp, 1, var)
  top <- names(sort(v, decreasing = TRUE))[1:min(5000, length(v))]
  Et <- E_imp[top, ]
  sc <- cor(Et, method = "spearman")
  qc$connectivity <- colSums(sc) - 1
  qc$connectivity_z <- scale(qc$connectivity)[, 1]
  qc$outlier_flag <- qc$connectivity_z < -3
  cat(gse, "potential outliers (connectivity z < -3):", sum(qc$outlier_flag), "\n")
  if (any(qc$outlier_flag)) print(qc[qc$outlier_flag, c("sample_id", "connectivity_z")])

  png(file.path(qc_dir, paste0(gse, "_pca.png")), width = 1000, height = 800, res = 130)
  p <- prcomp(t(E_imp[top, ]), center = TRUE, scale. = FALSE)
  pv <- summary(p)$importance[2, 1:5] * 100
  if ("tissue" %in% colnames(pn)) {
    colmap <- ifelse(pn$`tissue:ch1` == "tumor", "#D55E00", "#0072B2")
  } else {
    colmap <- "grey40"
  }
  plot(p$x[, 1], p$x[, 2], col = colmap, pch = 16,
       xlab = paste0("PC1 (", round(pv[1], 1), "%)"),
       ylab = paste0("PC2 (", round(pv[2], 1), "%)"),
       main = paste(gse, "PCA (top 5000 variable genes)"))
  if ("tissue" %in% colnames(pn)) legend("topright", c("normal", "tumor"), col = c("#0072B2", "#D55E00"), pch = 16)
  dev.off()

  png(file.path(qc_dir, paste0(gse, "_hclust.png")), width = 1100, height = 750, res = 130)
  hc <- hclust(dist(t(Et)), method = "average")
  plot(hc, labels = FALSE, main = paste(gse, "sample hierarchical clustering"))
  dev.off()

  .save_tsv(qc, "results_qc", paste0(gse, "_sample_qc.csv"))

  saveRDS(obj, file.path(data_proc, paste0(gse, "_analysis.rds")))
  .log("Saved analysis object for ", gse, ": genes = ", nrow(E_imp), ", samples = ", ncol(E_imp))
  invisible(obj)
}

gse47 <- process_gse("GSE70947")
gse05 <- process_gse("GSE70905")
gse52 <- process_gse("GSE15852")

# QC summary table
qc_summary <- data.frame(
  gse = c("GSE70947", "GSE70905", "GSE15852"),
  genes_after_controls = c(nrow(gse47$expression), nrow(gse05$expression), nrow(gse52$expression)),
  genes_dropped_na = c(sum(rowMeans(is.na(gse47$expression)) > max_na_frac),
                       sum(rowMeans(is.na(gse05$expression)) > max_na_frac),
                       sum(rowMeans(is.na(gse52$expression)) > max_na_frac)),
  samples = c(ncol(gse47$expression), ncol(gse05$expression), ncol(gse52$expression)),
  outliers_reported = c(sum(colMeans(is.na(gse47$expression)) * 0 + (rowMeans(is.na(gse47$expression)) * 0)) * 0,
                        sum(read.csv(file.path(qc_dir, "GSE70905_sample_qc.csv"))$outlier_flag),
                        sum(read.csv(file.path(qc_dir, "GSE15852_sample_qc.csv"))$outlier_flag)),
  stringsAsFactors = FALSE
)
qc_summary$outliers_reported[1] <- sum(read.csv(file.path(qc_dir, "GSE70947_sample_qc.csv"))$outlier_flag)
.save_tsv(qc_summary, "results_qc", "qc_summary.csv")
.log("Preprocessing complete")
cat("DONE\n")