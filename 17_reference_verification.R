## 17_reference_verification.R
## Verify every reference cited in the manuscript against PubMed (NCBI eutils).
## For each reference we record the PMID, whether it resolves, and the official
## title/journal/year as returned by NCBI so the manuscript bibliography is
## fact-checked (no fabricated references).
## Outputs:
##   results/audit/reference_verification.tsv
##   results/audit/reference_verification_report.md

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
audit_dir <- .dir_for("results_audit")

## candidate references (PMID or empty for software/database/URL references)
refs <- data.frame(
  key = c(
    "dixon2012", "stockwell2017", "galluzzi2018", "vandenberghe2014",
    "ferrdb2020", "quigley2017", "kanehisa2023", "langfelder2008",
    "hanzelmann2013", "foroutan2018", "aran2017", "wu2021",
    "ritchie2015"
  ),
  pmids = c(
    "22632970", "28985560", "29362479", "24452471",
    "32219413", "29147603", "36300620", "19114008",
    "23323831", "30400809", "29141660", "34557778",
    "25605792"
  ),
  expected_title = c(
    "Ferroptosis: an iron-dependent form of nonapoptotic cell death",
    "Ferroptosis: A Regulated Cell Death Nexus Linking Metabolism, Redox Biology, and Disease",
    "Molecular mechanisms of cell death: recommendations of the Nomenclature Committee on Cell Death 2018",
    "Regulated necrosis: the expanding network of non-apoptotic cell death pathways",
    "FerrDb: a manually curated resource for regulators and markers of ferroptosis and ferroptosis-disease associations",
    "Age, estrogen, and immune response in breast adenocarcinoma and adjacent normal tissue",
    "KEGG for taxonomy-based analysis of pathways and genomes",
    "WGCNA: an R package for weighted correlation network analysis",
    "GSVA: gene set variation analysis for microarray and RNA-seq data",
    "Single sample scoring of molecular phenotypes",
    "xCell: digitally portraying the tissue cellular heterogeneity landscape",
    "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data",
    "limma powers differential expression analyses for RNA-sequencing and microarray studies"
  ),
  kind = c(
    "primary_ferroptosis", "review_ferroptosis", "review_cell_death", "review_regulated_necrosis",
    "database", "cohort_source", "database", "software",
    "software", "software", "software", "software",
    "software"
  ),
  stringsAsFactors = FALSE
)

refs <- refs[refs$pmids != "", ]
refs$pmid <- as.integer(refs$pmids)

fetch_pubmed <- function(pmid) {
  url <- sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%d&retmode=json", pmid)
  ok <- tryCatch({
    js <- jsonlite::fromJSON(url)
    r <- js$result[[as.character(pmid)]]
    list(
      pmid = as.character(pmid),
      title = if (!is.null(r$title)) r$title else NA,
      journal = if (!is.null(r$fulljournalname)) r$fulljournalname else NA,
      pubdate = if (!is.null(r$pubdate)) r$pubdate else NA,
      authors = if (!is.null(r$authors)) paste(r$authors$name[1:min(3, length(r$authors$name))], collapse = "; ") else NA
    )
  }, error = function(e) list(pmid = as.character(pmid), title = NA, journal = NA, pubdate = NA, authors = NA))
  ok
}

rows <- lapply(seq_len(nrow(refs)), function(i) {
  info <- fetch_pubmed(refs$pmid[i])
  title_ok <- !is.na(info$title) &&
    (is.na(refs$expected_title[i]) || refs$expected_title[i] == "" ||
       grepl(substr(refs$expected_title[i], 1, 25), info$title, ignore.case = TRUE))
  data.frame(
    key = refs$key[i], pmid = refs$pmid[i], kind = refs$kind[i],
    expected_title = refs$expected_title[i],
    returned_title = info$title, journal = info$journal,
    pubdate = info$pubdate, authors = info$authors,
    title_matches_expected = title_ok,
    verified = !is.na(info$title) && title_ok,
    stringsAsFactors = FALSE
  )
})
ver <- do.call(rbind, rows)
.save_tsv(ver, "results_audit", "reference_verification.tsv")

report <- c(
  "# Reference verification report",
  "",
  sprintf("Verified %d/%d references against PubMed (generated %s).",
          sum(ver$verified), nrow(ver), format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "| key | PMID | verified | returned title | journal | pubdate |",
  "|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(ver))) {
  report <- c(report, sprintf("| %s | %s | %s | %s | %s | %s |",
    ver$key[i], ver$pmid[i], ifelse(ver$verified[i], "YES", "NO"),
    ver$returned_title[i], ver$journal[i], ver$pubdate[i]))
}
writeLines(report, file.path(audit_dir, "reference_verification_report.md"))

.log("Reference verification complete")
cat("DONE\n")