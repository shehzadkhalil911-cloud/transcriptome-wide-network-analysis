## 02_gene_set_audit.R
## Audit the ferroptosis and necroptosis gene sets.
## Outputs:
##   data/gene_sets/ferroptosis_genes_clean.txt
##   data/gene_sets/necroptosis_genes_clean.txt
##   data/gene_sets/ferroptosis_necroptosis_overlap.txt
##   results/tables/gene_set_summary.csv
##   results/tables/ferroptosis_necroptosis_overlap.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
tab_dir <- .dir_for("results_tables")
gs_dir <- .dir_for("gene_sets")
if (!dir.exists(gs_dir)) dir.create(gs_dir, recursive = TRUE)

read_gs <- function(file) {
  d <- read_excel(file, sheet = 1)
  g <- as.character(d[[1]])
  g <- trimws(g)
  g <- g[!is.na(g) & g != ""]
  g
}

ferr <- read_gs(file.path(cfg$project$dirs$data_raw, "Ferroptosis Genes.xlsx"))
necr <- read_gs(file.path(cfg$project$dirs$data_raw, "Necroptosis Genes.xlsx"))

cat("Ferroptosis raw entries:", length(ferr), "\n")
cat("Necroptosis raw entries:", length(necr), "\n")

gs_summary <- data.frame(
  gene_set = c("ferroptosis", "necroptosis"),
  source_file = c("Ferroptosis Genes.xlsx", "Necroptosis Genes.xlsx"),
  source_version = c("not recorded", "not recorded"),
  total_entries = c(length(ferr), length(necr)),
  unique_genes = c(length(unique(ferr)), length(unique(necr))),
  duplicate_entries = c(length(ferr) - length(unique(ferr)), length(necr) - length(unique(necr))),
  identifier_type = c("HGNC gene symbol (as listed)", "HGNC gene symbol (as listed)"),
  overlap_with_other_set = NA_integer_,
  stringsAsFactors = FALSE
)

overlap <- intersect(unique(ferr), unique(necr))
gs_summary$overlap_with_other_set <- c(length(overlap), length(overlap))

cat("Overlap ferroptosis x necroptosis:", length(overlap), "genes\n")
cat("Overlap genes:", paste(overlap, collapse = ", "), "\n")

# cleaned unique sets
ferr_uniq <- unique(ferr)
necr_uniq <- unique(necr)
writeLines(sort(ferr_uniq), file.path(gs_dir, "ferroptosis_genes_clean.txt"))
writeLines(sort(necr_uniq), file.path(gs_dir, "necroptosis_genes_clean.txt"))
writeLines(sort(overlap), file.path(gs_dir, "ferroptosis_necroptosis_overlap.txt"))
cat("unique ferroptosis:", length(ferr_uniq), " unique necroptosis:", length(necr_uniq), "\n")

# overlap table
ov <- data.frame(
  gene = sort(overlap),
  ferroptosis_member = TRUE,
  necroptosis_member = TRUE
)
.save_tsv(ov, "results_tables", "ferroptosis_necroptosis_overlap.csv")
.save_tsv(gs_summary, "results_tables", "gene_set_summary.csv")

# compare against external reference sets (KEGG via msigdbr) for provenance note
## Download authoritative reference gene sets (recorded in download_manifest.csv)
audit_dir <- .dir_for("results_audit")
download_and_record <- function(url, filename, reason, dest_dir = file.path("data", "gene_sets", "external")) {
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest <- file.path(dest_dir, filename)
  if (!file.exists(dest)) {
    ok <- try(download.file(url, dest, mode = "wb", quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) {
      cat("DOWNLOAD FAILED:", url, "\n")
      return(NA_character_)
    }
  }
  md5 <- tools::md5sum(dest)
  dl <- data.frame(url = url, source = ifelse(grepl("kegg", url), "KEGG REST API",
                     ifelse(grepl("ferrdb|zhounan", url), "FerrDb", "other")),
                   date = as.character(Sys.Date()), filename = filename,
                   checksum = unname(md5), reason = reason, stringsAsFactors = FALSE)
  dl_old <- if (file.exists(file.path(audit_dir, "download_manifest.csv"))) {
    read.csv(file.path(audit_dir, "download_manifest.csv"), check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")
  } else data.frame(url=character(), source=character(), date=character(), filename=character(), checksum=character(), reason=character())
  dl_new <- rbind(dl_old, dl)
  .save_tsv(dl_new, "results_audit", "download_manifest.csv")
  dest
}

# KEGG ferroptosis hsa04216 and necroptosis hsa04217
fetch_kegg <- function(path_id) {
  lk <- readLines(url(paste0("https://rest.kegg.jp/link/hsa/path:hsa", path_id)))
  ids <- sub(paste0("^path:hsa", path_id, "\\thsa:"), "", lk)
  lst <- readLines(url("https://rest.kegg.jp/list/hsa"))
  id2sym <- setNames(sub(";.*$", "", sub("^[^\t]*\t[^\t]*\t[^\t]*\t", "", lst)),
                     sub("^hsa:", "", sub("\t.*$", "", lst)))
  id2sym <- id2sym[!grepl("^(CDS|ncRNA|tRNA|rRNA|pseudo)$", id2sym)]
  syms <- unname(id2sym[ids])
  unique(syms[!is.na(syms) & syms != ""])
}

kegg_ferr_file <- download_and_record("https://rest.kegg.jp/link/hsa/path:hsa04216",
                                      "kegg_hsa04216_link.txt", "KEGG ferroptosis pathway (hsa04216) gene list for sensitivity analysis")
kegg_necr_file <- download_and_record("https://rest.kegg.jp/link/hsa/path:hsa04217",
                                      "kegg_hsa04217_link.txt", "KEGG necroptosis pathway (hsa04217) gene list for sensitivity analysis")

kegg_ferr <- try(fetch_kegg("04216"), silent = TRUE)
kegg_necr <- try(fetch_kegg("04217"), silent = TRUE)
cat("KEGG ferroptosis (hsa04216) symbol count:", length(kegg_ferr), "\n")
cat("KEGG necroptosis (hsa04217) symbol count:", length(kegg_necr), "\n")
if (!inherits(kegg_ferr, "try-error")) writeLines(sort(kegg_ferr), file.path(gs_dir, "kegg_ferroptosis_ref.txt"))
if (!inherits(kegg_necr, "try-error")) writeLines(sort(kegg_necr), file.path(gs_dir, "kegg_necroptosis_ref.txt"))

# FerrDb ferroptosis gene list (extended download, July 2024 update)
ferrdb_url <- "https://www.zhounan.org/ferrdb/current/extdownload/ext_20240730.xlsx"
ferrdb_file <- download_and_record(ferrdb_url, "ferrdb_ext_20240730.xlsx",
                                   "FerrDb extended download (ferroptosis + cuproptosis + panoptosis regulators) for sensitivity analysis")
ferrdb_ferro <- try({
  fd <- readxl::read_excel(ferrdb_file, sheet = 1)
  # columns may include a 'RCD' or death-type column; ferroptosis genes
  cat("FerrDb xlsx dim:", dim(fd)[1], "x", dim(fd)[2], "\n")
  print(names(fd))
  syms <- character(0)
  for (cn in c("gene", "Gene", "symbol", "Symbol", "regulator", "Regulator")) {
    if (cn %in% names(fd)) syms <- c(syms, as.character(fd[[cn]]))
  }
  unique(trimws(syms[!is.na(syms)]))
}, silent = TRUE)
if (!inherits(ferrdb_ferro, "try-error") && length(ferrdb_ferro) > 0) {
  writeLines(sort(ferrdb_ferro), file.path(gs_dir, "ferrdb_ferroptosis_ref.txt"))
  cat("FerrDb extended extracted symbol count:", length(ferrdb_ferro), "\n")
} else {
  cat("FerrDb extended extraction failed:\n", as.character(ferrdb_ferro)[1], "\n")
}

# FerrDb V1 ferroptosis gene lists (driver / suppressor / marker)
ferrdb_base <- "https://www.zhounan.org/ferrdb/legacy/fordownload"
ferr_v1 <- character(0)
for (cat_f in c("driver", "suppressor", "marker")) {
  u <- file.path(ferrdb_base, paste0("1_info_", cat_f, ".csv"))
  f <- download_and_record(u, paste0("ferrdb_v1_", cat_f, ".csv"),
                           paste0("FerrDb V1 ferroptosis ", cat_f, " gene list for sensitivity analysis"))
  d <- try(read.csv(f, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(d, "try-error") && nrow(d) > 0) {
    # find gene symbol column
    cn <- intersect(c("Symbol", "symbol", "Gene", "gene", "Gene_symbol", "gene_symbol"), names(d))
    if (length(cn) > 0) {
      ferr_v1 <- c(ferr_v1, as.character(d[[cn[1]]]))
      cat("FerrDb V1", cat_f, "genes:", nrow(d), "\n")
    } else {
      cat("FerrDb V1", cat_f, "columns:", paste(names(d), collapse=", "), "\n")
    }
  }
}
ferr_v1 <- unique(trimws(ferr_v1[!is.na(ferr_v1) & ferr_v1 != ""]))
if (length(ferr_v1) > 0) {
  writeLines(sort(ferr_v1), file.path(gs_dir, "ferrdb_ferroptosis_ref.txt"))
  cat("FerrDb V1 ferroptosis gene count (driver+suppressor+marker):", length(ferr_v1), "\n")
  cat("overlap local ferroptosis vs FerrDb V1:", length(intersect(ferr_uniq, ferr_v1)), "\n")
}

.log("Gene-set audit complete")
cat("DONE\n")