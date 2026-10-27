## 01_metadata_and_pairing.R
## Verify cohort structure and pairing for GSE70947 (discovery) and GSE70905 (validation).
## Outputs:
##   results/audit/GSE70947_pairing_audit.csv
##   results/audit/GSE70905_sample_audit.csv
##   data/processed/pairing_GSE70947.rds
##   data/processed/pairing_GSE70905.rds

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
audit_dir <- .dir_for("results_audit")
data_proc <- .dir_for("data_processed")

# ---------------- GSE70947 (discovery) ----------------
.log("GSE70947 pairing audit")
x <- readRDS(file.path(cfg$project$dirs$data_raw, "GSE70947_processed.rds"))
pn <- x$phenotype
sid <- pn$`sample id:ch1`
gsm <- pn$geo_accession
tissue <- pn$`tissue:ch1`
# patient id from sample id pattern CM016.N.c2 / CM016.T.c2
patient <- sub("\\.(N|T)\\.c2$", "", sid)
stopifnot(all(!grepl("\\.(N|T)\\.c2$", patient)))
tissue_sid <- ifelse(grepl("\\.N\\.c2$", sid), "normal", "tumor")
cat("mismatch tissue vs sample-id pattern:", sum(tissue_sid != tissue), "\n")

paired <- data.frame(
  sample_id = gsm,
  geo_accession = gsm,
  sample_id_raw = sid,
  patient_id = patient,
  tissue = tissue,
  er = pn$`er:ch1`,
  her2 = pn$`her2:ch1`,
  cohort = pn$`cohort:ch1`,
  age_years = suppressWarnings(as.numeric(pn$`age_in_years:ch1`)),
  bmi = suppressWarnings(as.numeric(pn$`bmi:ch1`)),
  tumor_size_cm = suppressWarnings(as.numeric(pn$`tumor_size_cm:ch1`)),
  stringsAsFactors = FALSE
)

# build per-patient pairing table
tum <- paired[paired$tissue == "tumor", ]
nor <- paired[paired$tissue == "normal", ]
tab <- merge(
  data.frame(patient_id = tum$patient_id, tumor_sample_id = tum$geo_accession,
             tumor_phenotype = tum$tissue, tumor_er = tum$er, tumor_her2 = tum$her2,
             tumor_age = tum$age_years),
  data.frame(patient_id = nor$patient_id, normal_sample_id = nor$geo_accession,
             normal_phenotype = nor$tissue, normal_er = nor$er),
  by = "patient_id", all = TRUE
)
tab$pairing_status <- ifelse(!is.na(tab$tumor_sample_id) & !is.na(tab$normal_sample_id), "paired", "unpaired")

cat("n patients:", nrow(tab), "\n")
cat("n paired:", sum(tab$pairing_status == "paired"), "\n")
cat("n unpaired:", sum(tab$pairing_status == "unpaired"), "\n")
stopifnot(nrow(tab) == 148, all(tab$pairing_status == "paired"))

.save_tsv(tab, "results_audit", "GSE70947_pairing_audit.csv")
saveRDS(list(paired_samples = paired, pairing = tab), file.path(data_proc, "pairing_GSE70947.rds"))

# ---------------- GSE70905 (validation) ----------------
.log("GSE70905 sample audit")
y <- readRDS(file.path(cfg$project$dirs$data_raw, "GSE70905_processed.rds"))
pn5 <- y$phenotype
# group by source_name_ch1: tumour / adjacent normal / reduction
grp <- ifelse(pn5$source_name_ch1 == "breast adenocarcinoma", "tumor",
       ifelse(pn5$source_name_ch1 == "breast tissue adjacent to tumor", "adjacent_normal",
       ifelse(pn5$source_name_ch1 == "mammaplastic reduction", "reduction_normal", "other")))
# patient key: strip trailing N/T/B from raw sample id
raw <- pn5$description
pkey <- sub("(N|T|B)$", "", raw)
cat("sample_type vs source-derived group table:\n")
print(table(pn5$sample_type, grp))

audit5 <- data.frame(
  sample_id = pn5$geo_accession,
  sample_id_raw = raw,
  patient_id_raw = pkey,
  patient_id_field = pn5$patient_id,
  source_name = pn5$source_name_ch1,
  group = grp,
  tissue_field = pn5$`tissue:ch1`,
  er = pn5$`er:ch1`,
  her2 = pn5$`her2:ch1`,
  cohort = pn5$`cohort:ch1`,
  stringsAsFactors = FALSE
)
cat("group counts:\n"); print(table(audit5$group))

# paired tumour/adjacent-normal detection
tum5 <- audit5[audit5$group == "tumor", ]
adj5 <- audit5[audit5$group == "adjacent_normal", ]
pairs5 <- merge(
  data.frame(patient_key = tum5$patient_id_raw, tumor_sample_id = tum5$sample_id),
  data.frame(patient_key = adj5$patient_id_raw, adjacent_normal_sample_id = adj5$sample_id),
  by = "patient_key", all = TRUE
)
pairs5$pairing_status <- ifelse(!is.na(pairs5$tumor_sample_id) & !is.na(pairs5$adjacent_normal_sample_id), "paired", "unpaired")
cat("GSE70905 tumour-adjacent pairing: n paired =", sum(pairs5$pairing_status == "paired"),
    " n unpaired =", sum(pairs5$pairing_status == "unpaired"), "\n")

.save_tsv(audit5, "results_audit", "GSE70905_sample_audit.csv")
saveRDS(list(sample_audit = audit5, pairing = pairs5), file.path(data_proc, "pairing_GSE70905.rds"))

.log("Pairing audits complete")
cat("DONE\n")