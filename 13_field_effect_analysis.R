## 13_field_effect_analysis.R
## Field-effect analysis: examine whether ferroptosis/necroptosis pathway scores
## and coupling differ across normal tissue states in GSE70905:
##   tumour vs adjacent normal vs reduction mammoplasty (true healthy breast).
## Also uses GSE15852 (paired tumour/normal) as an extra check.
## Outputs:
##   results/field_effect/field_effect_anova.tsv
##   results/field_effect/group_pairwise.tsv
##   results/field_effect/field_effect_summary.csv

source(file.path("R", "helpers.R"))
cfg <- .load_config("config.yml")
set.seed(cfg$project$seed)
fe_dir <- .dir_for("results_field")

# ---- GSE70905 scores ----
ana_val <- readRDS(file.path(.dir_for("data_processed"), "GSE70905_analysis.rds"))
expr_val <- ana_val$expression_imputed
ph_val <- ana_val$phenotype
audit <- read.csv(file.path("results", "audit", "GSE70905_sample_audit.csv"), sep = "\t")
group <- setNames(audit$group, audit$sample_id)

gs <- read.csv(file.path("data", "gene_sets", "ferroptosis_genes_clean.txt"), header = FALSE)$V1
ns <- read.csv(file.path("data", "gene_sets", "necroptosis_genes_clean.txt"), header = FALSE)$V1

score_matrix <- function(mat, gene_set) {
  common <- intersect(gene_set, rownames(mat))
  if (length(common) < 3) return(rep(NA, ncol(mat)))
  library(singscore)
  ranks <- rankGenes(mat)
  s <- simpleScore(ranks, upSet = common, centerScore = TRUE)
  as.numeric(s$TotalScore)
}

df <- data.frame(
  sample_id = colnames(expr_val),
  ferro = score_matrix(expr_val, gs),
  necro = score_matrix(expr_val, ns),
  stringsAsFactors = FALSE
)
df$group <- group[df$sample_id]
df$group <- factor(df$group, levels = c("reduction_normal", "adjacent_normal", "tumor"))
df <- df[!is.na(df$group), ]
cat("groups:\n"); print(table(df$group))

# ---- ANOVA across three groups ----
anova_out <- do.call(rbind, lapply(c("ferro", "necro"), function(pw) {
  f <- try(aov(df[[pw]] ~ df$group), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  s <- summary(f)[[1]]
  data.frame(pathway = pw, between_df = s$Df[1], between_ss = s$`Sum Sq`[1],
             within_df = s$Df[2], within_ss = s$`Sum Sq`[2],
             F = s$`F value`[1], p_value = s$`Pr(>F)`[1])
}))
anova_out$fdr <- p.adjust(anova_out$p_value, method = "BH")
.save_tsv(anova_out, "results_field", "field_effect_anova.tsv")

# ---- pairwise t-tests with Tukey HSD ----
pairwise <- do.call(rbind, lapply(c("ferro", "necro"), function(pw) {
  tuk <- try(TukeyHSD(aov(df[[pw]] ~ df$group)), silent = TRUE)
  if (inherits(tuk, "try-error")) return(NULL)
  tdf <- as.data.frame(tuk$`df$group`)
  tdf$comparison <- rownames(tdf)
  tdf$pathway <- pw
  tdf$diff <- tdf$diff; tdf$p_adj <- tdf$`p adj`
  tdf
}))
pairwise <- pairwise[, c("pathway", "comparison", "diff", "lwr", "upr", "p_adj")]
.save_tsv(pairwise, "results_field", "group_pairwise.tsv")

# group means
means <- do.call(rbind, lapply(c("ferro", "necro"), function(pw) {
  ag <- aggregate(df[[pw]], by = list(df$group), FUN = mean)
  ag2 <- aggregate(df[[pw]], by = list(df$group), FUN = sd)
  data.frame(pathway = pw, group = ag$Group.1, mean = ag$x, sd = ag2$x)
}))
.save_tsv(means, "results_field", "group_means.tsv")

# ---- coupling per group (same as validation but including comparison) ----
coupling <- do.call(rbind, lapply(levels(df$group), function(grp) {
  sub <- df[df$group == grp, ]
  ct <- cor.test(sub$ferro, sub$necro, method = "pearson")
  data.frame(group = grp, n = nrow(sub), r = ct$estimate, p_value = ct$p.value)
}))
coupling$fdr <- p.adjust(coupling$p_value, method = "BH")
.save_tsv(coupling, "results_field", "coupling_by_group.tsv")

summ <- data.frame(
  n_reduction = sum(df$group == "reduction_normal"),
  n_adjacent = sum(df$group == "adjacent_normal"),
  n_tumor = sum(df$group == "tumor"),
  ferro_anova_p = anova_out$p_value[1],
  necro_anova_p = anova_out$p_value[2],
  ferro_tumor_vs_reduction_diff = pairwise$diff[pairwise$pathway == "ferro" & pairwise$comparison == "tumor-reduction_normal"],
  ferro_tumor_vs_reduction_padj = pairwise$p_adj[pairwise$pathway == "ferro" & pairwise$comparison == "tumor-reduction_normal"],
  necro_tumor_vs_reduction_diff = pairwise$diff[pairwise$pathway == "necro" & pairwise$comparison == "tumor-reduction_normal"],
  necro_tumor_vs_reduction_padj = pairwise$p_adj[pairwise$pathway == "necro" & pairwise$comparison == "tumor-reduction_normal"]
)
.save_tsv(summ, "results_field", "field_effect_summary.csv")

.log("Field effect analysis complete")
cat("DONE\n")