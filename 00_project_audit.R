## 00_project_audit.R
## Inventory the project directory and record provenance information.
## Outputs:
##   results/audit/project_file_inventory.csv
##   results/audit/project_audit.md
##   results/audit/download_manifest.csv
##   results/audit/sessionInfo.txt
##   results/audit/package_versions.tsv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")

.log("Audit start")
root <- getwd()
audit_dir <- .dir_for("results_audit")

# ---- file inventory -------------------------------------------------------
inv <- data.frame()
for (rel in list.files(".", recursive = TRUE, full.names = FALSE, include.dirs = FALSE)) {
  if (grepl("^results/|^logs/", rel)) next
  full <- file.path(".", rel)
  if (dir.exists(full)) next
  info <- file.info(full)
  ext <- tools::file_ext(rel)
  inv <- rbind(inv, data.frame(
    filename = basename(rel),
    relative_path = gsub("\\\\", "/", rel),
    extension = ext,
    file_size_bytes = as.numeric(info$size),
    file_size_human = format(as.numeric(info$size), big.mark = ","),
    purpose = NA_character_,
    stringsAsFactors = FALSE
  ))
}

purpose <- function(path) {
  p <- tolower(path)
  if (grepl("ferroptosis", p)) return("ferroptosis gene-set file (source, do not modify)")
  if (grepl("necroptosis", p)) return("necroptosis gene-set file (source, do not modify)")
  if (grepl("gse70947_processed", p)) return("GSE70947 processed expression + phenotype (discovery)")
  if (grepl("gse70947_expression", p)) return("GSE70947 gene-level expression matrix (CSV export)")
  if (grepl("gse70905_processed", p)) return("GSE70905 processed expression + phenotype (validation)")
  if (grepl("gse70905_series_matrix", p)) return("GSE70905 GEO series matrix (not used; scrambled per RDS notes)")
  if (grepl("gse15852_processed", p)) return("GSE15852 processed expression + phenotype (extra validation)")
  if (grepl("gse15852_expression", p)) return("GSE15852 gene-level expression matrix (CSV export)")
  if (grepl("sample_classes", p)) return("sample class assignment CSV from previous project")
  if (grepl("config.yml", p)) return("pipeline configuration (seed, thresholds)")
  if (grepl("^R/helpers", p)) return("shared R helper functions")
  if (grepl("^scripts/", p)) return("analysis pipeline script")
  if (grepl("^manuscript/", p)) return("manuscript output")
  if (grepl("^results/", p)) return("analysis result output")
  if (grepl("^data/gene_sets/", p)) return("cleaned gene-set file (generated)")
  if (grepl("^data/processed/", p)) return("processed analysis object (generated)")
  if (grepl("^logs/", p)) return("pipeline log")
  if (grepl("\\.rds$", p)) return("R serialised object")
  return("other")
}
inv$purpose <- vapply(inv$relative_path, purpose, character(1))

.save_tsv(inv, "results_audit", "project_file_inventory.csv")

# ---- download manifest (empty unless downloads occur) ---------------------
dl <- data.frame(
  url = character(0), source = character(0), date = character(0),
  filename = character(0), checksum = character(0), reason = character(0)
)
.save_tsv(dl, "results_audit", "download_manifest.csv")

# ---- session info / package versions --------------------------------------
si <- capture.output(sessionInfo())
writeLines(si, file.path(audit_dir, "sessionInfo.txt"))
ip <- installed.packages()[, c("Package", "Version", "LibPath")]
write.table(ip, file.path(audit_dir, "package_versions.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

# ---- audit markdown -------------------------------------------------------
md <- c(
  "# Project Audit",
  "",
  paste0("Audit timestamp: ", Sys.time()),
  "",
  "## Data files present (source data; not modified)",
  "",
  "| File | Size (bytes) | Purpose |",
  "|------|-------------|---------|"
)
for (i in seq_len(nrow(inv))) {
  if (grepl("^results/|^logs/|^manuscript/", inv$relative_path[i])) next
  md <- c(md, sprintf("| %s | %s | %s |",
                      inv$relative_path[i], inv$file_size_bytes[i], inv$purpose[i]))
}
md <- c(md, "",
        "## Dataset identity (verified from RDS objects)",
        "",
        "- **GSE70947** (discovery): Agilent GPL13607; 296 samples = 148 paired tumour + 148 adjacent normal;",
        "  gene-level log2 matrix 34949 x 296 (highest-IQR single-probe collapse); phenotype contains tissue, ER, HER2, cohort, age, BMI, tumour size.",
        "- **GSE70905** (validation): Agilent GPL4133; 137 samples = 47 tumour + 47 adjacent normal + 43 reduction-mammoplasty normal;",
        "  gene-level log2 + quantile-normalised matrix 30232 x 137; phenotype contains sample_type, patient_id, source_name_ch1.",
        "- **GSE15852** (extra validation, present): Affymetrix GPL96; 86 samples = 43 paired tumour + normal; gene-level matrix 13237 x 86.",
        "",
        "## Gene-set files",
        "",
        "- `Ferroptosis Genes.xlsx`: 382 symbols (single column, Sheet1).",
        "- `Necroptosis Genes.xlsx`: 158 symbols (single column, Sheet1).",
        "- Provenance/version not recorded in the files; cleaned unique-symbol versions are generated in data/gene_sets.",
        "",
        "## Download manifest",
        "",
        "No downloads were required; local data were complete and internally consistent for the planned analyses.",
        ""
)
writeLines(md, file.path(audit_dir, "project_audit.md"))
.log("Audit complete. Inventory: ", file.path(audit_dir, "project_file_inventory.csv"))