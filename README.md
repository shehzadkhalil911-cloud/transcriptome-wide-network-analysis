 Ferroptosis–necroptosis crosstalk in paired breast tumour and adjacent normal tissue

A fully reproducible, publication-grade transcriptomic study testing whether
ferroptosis and necroptosis are co-regulated (coupled) at the level of
individual breast tumours.

## Key result

Across 148 discovery patients (GSE70947), tumour-associated changes
(tumour − adjacent normal) in ferroptosis and necroptosis pathway activity were
positively correlated (Pearson r = 0.361, p = 6.6×10⁻⁶; permutation p = 5.0×10⁻⁴),
robust to scoring method (GSVA/ssGSEA r = 0.407), leave-one-out deletion
(r range 0.335–0.380), ER status and tumour-composition adjustment. A signed
WGCNA network identified 11 crosstalk modules anchored in an
inflammation/cytotoxic-immune axis (ME9 = T-cell/IFN-γ; ME30 = TNF/NF-κB/IL-17),
and the coupling replicated in an independent cohort (GSE70905, r = 0.865,
p = 4.4×10⁻¹⁵).

The full manuscript (markdown + Word + PDF, with all figures embedded) is in
`manuscript/` and is generated programmatically from the result files by
`scripts/18_manuscript_generation.R` — every number in the text is read from the
analysis outputs, so the manuscript cannot drift from the results.

## Cohorts

| cohort | platform | samples | design |
|--------|----------|---------|--------|
| GSE70947 | GPL13607 (Agilent) | 296 | 148 paired tumour + adjacent normal (discovery) |
| GSE70905 | GPL4133 (Agilent) | 137 | 47 tumour + 47 adjacent normal + 43 reduction normal (validation) |
| GSE15852 | GPL96 (Affymetrix) | 86 | 43 paired tumour + normal (additional) |

## Pipeline

Scripts run in order; each reads only prior outputs and `config.yml`.

| script | step | key outputs (`results/…`) |
|--------|------|---------------------------|
| `00_project_audit.R` | inventory source data | `audit/project_audit.md` |
| `01_metadata_and_pairing.R` | rebuild sample pairing | `audit/GSE*_*.csv` |
| `02_gene_set_audit.R` | ferroptosis/necroptosis gene sets | `tables/gene_set_summary.csv` |
| `03_preprocessing.R` | QC, imputation, log2, collapse | `data/processed/*.rds` |
| `04_differential_expression.R` | paired DE of pathway genes | `differential_expression/` |
| `05_pathway_scoring.R` | singscore (primary) + GSVA/ssGSEA | `pathway_scores/*.rds` |
| `06_crosstalk_analysis.R` | within-tissue coupling, interaction models | `crosstalk/within_tissue_correlations.csv`, `coupling_interaction_models.csv` |
| `07_delta_analysis.R` | patient-level delta crosstalk, permutation, bootstrap, LOO | `crosstalk/patient_delta_crosstalk.csv`, `delta_correlation_leave_one_out.csv` |
| `08_wgcna_delta_expression.R` | signed WGCNA on deltas | `wgcna/`, `enrichment/crosstalk_modules_defined.csv` |
| `09_enrichment.R` | GO/KEGG/Reactome/Hallmark ORA | `enrichment/*.tsv` |
| `10_hub_prioritisation.R` | composite hub scoring | `hubs/` |
| `11_composition_analysis.R` | xCell composition, adjusted coupling | `composition/` |
| `12_external_validation.R` | replicate coupling/modules/hubs in GSE70905 | `validation/` |
| `13_field_effect_analysis.R` | tumour vs adjacent vs reduction ANOVA | `field_effect/` |
| `14_sensitivity_analysis.R` | scoring/Spearman/LOO/ER/composition | `sensitivity/` |
| `15_figures.R` | Fig1–Fig10 + SuppA–C | `figures/`, `supplementary/` |
| `16_tables.R` | Table1–Table9 | `tables/` |
| `17_reference_verification.R` | fact-check all references vs PubMed | `audit/reference_verification_report.md` (13/13 verified) |
| `18_manuscript_generation.R` | build manuscript `.md`/`.docx`/`.pdf` | `manuscript/` |
| `19_final_audit.R` | cross-check every manuscript number/fig/table/ref | `audit/final_audit_summary.md`, `claim_evidence_matrix.tsv` |



## Requirements

R ≥ 4.2 with the packages in `R/helpers.R` (yaml, limma, WGCNA, GSVA, singscore,
clusterProfiler, org.Hs.eg.db, msigdbr, ggplot2, dplyr, tidyr, readxl, data.table,
jsonlite), plus pandoc and a LaTeX engine (xelatex) for the `.docx`/`.pdf`
manuscript conversion.
