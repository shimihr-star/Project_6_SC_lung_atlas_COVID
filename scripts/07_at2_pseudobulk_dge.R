


# =============================================================================
# 07_at2_pseudobulk_dge.R
# AT2 Pseudobulk Differential Gene Expression Analysis
# COVID-19 vs Control
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load required packages
# -----------------------------------------------------------------------------

library(edgeR)


# -----------------------------------------------------------------------------
# 2. Define project paths
# -----------------------------------------------------------------------------

project_dir <- paste0(
  "G:/My Drive/Shimon/stuffs/Data_Science_New_way/",
  "Portfolio/Self_projects/Project_6_SC_lung_atlas_COVID"
)

counts_path <- file.path(
  project_dir,
  "data/results/at2_pseudobulk_counts.csv"
)

qc_path <- file.path(
  project_dir,
  "data/results/at2_pseudobulk_qc.csv"
)

clinical_path <- file.path(
  project_dir,
  "data/results/donor_metadata.csv"
)


# -----------------------------------------------------------------------------
# 3. Load AT2 pseudobulk counts
# -----------------------------------------------------------------------------

# Python exported donors as rows and genes as columns.
# edgeR requires genes as rows and donors as columns.

counts_df <- read.csv(
  counts_path,
  row.names = 1,
  check.names = FALSE
)

counts <- t(as.matrix(counts_df))

cat("\nAT2 count matrix dimensions:\n")
print(dim(counts))


# -----------------------------------------------------------------------------
# 4. Validate raw counts
# -----------------------------------------------------------------------------

stopifnot(
  is.numeric(counts),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts)),
  !anyDuplicated(colnames(counts))
)


# -----------------------------------------------------------------------------
# 5. Load AT2 sample QC
# -----------------------------------------------------------------------------

qc <- read.csv(
  qc_path,
  stringsAsFactors = FALSE
)

qc <- qc[
  match(colnames(counts), qc$donor_id),
  ,
  drop = FALSE
]

stopifnot(
  !anyNA(qc$donor_id),
  identical(qc$donor_id, colnames(counts))
)

cat("\nAT2 nuclei per donor:\n")
print(qc[, c("donor_id", "group", "n_nuclei")])


# -----------------------------------------------------------------------------
# 6. Load and align clinical metadata
# -----------------------------------------------------------------------------

clinical <- read.csv(
  clinical_path,
  stringsAsFactors = FALSE
)

clinical <- clinical[
  match(colnames(counts), clinical$donor_id),
  ,
  drop = FALSE
]

stopifnot(
  !anyNA(clinical$donor_id),
  identical(clinical$donor_id, colnames(counts)),
  identical(as.character(clinical$group), as.character(qc$group))
)

clinical$age <- as.numeric(clinical$age)
clinical$sex <- factor(clinical$sex)
clinical$group <- factor(
  clinical$group,
  levels = c("Control", "COVID-19")
)

cat("\nAT2 donors by group:\n")
print(table(clinical$group))

cat("\nMissing clinical information:\n")
print(colSums(is.na(clinical)))



# =============================================================================
# 7. Construct Adjusted Model and Normalize AT2 Counts
# =============================================================================

# -----------------------------------------------------------------------------
# 7.1. Select donors with complete clinical information
# -----------------------------------------------------------------------------

complete_samples <- complete.cases(
  clinical[, c("age", "sex", "group")]
)

clinical_adjusted <- droplevels(
  clinical[complete_samples, , drop = FALSE]
)

# Select matching count matrix columns
counts_adjusted <- counts[
  ,
  clinical_adjusted$donor_id,
  drop = FALSE
]

cat("\nDonors in adjusted AT2 model:\n")
print(table(clinical_adjusted$group))


# -----------------------------------------------------------------------------
# 7.2. Construct design matrix
# -----------------------------------------------------------------------------

design <- model.matrix(
  ~ age + sex + group,
  data = clinical_adjusted
)

stopifnot(
  identical(colnames(counts_adjusted),
            clinical_adjusted$donor_id),
  qr(design)$rank == ncol(design)
)

print(colnames(design))
print(dim(design))


# -----------------------------------------------------------------------------
# 7.3. Create edgeR object
# -----------------------------------------------------------------------------

dge <- DGEList(
  counts = counts_adjusted,
  group = clinical_adjusted$group
)


# -----------------------------------------------------------------------------
# 7.4. Filter low-expression genes
# -----------------------------------------------------------------------------

keep <- filterByExpr(
  dge,
  design = design
)

cat("\nGenes before filtering:", nrow(dge), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]


# -----------------------------------------------------------------------------
# 7.5. TMM normalization
# -----------------------------------------------------------------------------

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

cat("\nNormalized AT2 sample information:\n")
print(dge$samples)

cat("\nFinal AT2 DGEList dimensions:\n")
print(dim(dge))


cat("\nFinal AT2 DGEList dimensions:\n")
print(dim(dge))

# -----------------------------------------------------------------------------
# 7.6. Export TMM-normalized log2 CPM for Python heatmaps
# -----------------------------------------------------------------------------

# Calculate normalized log2 counts per million.
log_cpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)

# Transpose: donors as rows, genes as columns.
log_cpm_df <- as.data.frame(t(log_cpm))

# Save normalized expression.
output_path <- file.path(
  project_dir,
  "data/results/at2_log2cpm.csv"
)

write.csv(
  log_cpm_df,
  file = output_path,
  row.names = TRUE
)

cat("\nAT2 normalized expression dimensions:\n")
print(dim(log_cpm_df))

cat("\nSaved to:", output_path, "\n")

# =============================================================================
# 8. AT2 Differential Expression Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# 8.1. Estimate dispersion
# -----------------------------------------------------------------------------

dge <- estimateDisp(
  dge,
  design = design,
  robust = TRUE
)

cat("\nAT2 common dispersion:\n")
print(dge$common.dispersion)

cat("\nAT2 biological coefficient of variation:\n")
print(sqrt(dge$common.dispersion))


# -----------------------------------------------------------------------------
# 8.2. Fit quasi-likelihood negative binomial GLM
# -----------------------------------------------------------------------------

fit <- glmQLFit(
  dge,
  design = design,
  robust = TRUE
)


# -----------------------------------------------------------------------------
# 8.3. Test COVID-19 vs Control
# -----------------------------------------------------------------------------

qlf <- glmQLFTest(
  fit,
  coef = "groupCOVID-19"
)

dge_results <- topTags(
  qlf,
  n = Inf,
  adjust.method = "BH"
)$table


# -----------------------------------------------------------------------------
# 8.4. Summarize significant genes
# -----------------------------------------------------------------------------

fdr_threshold <- 0.05
logfc_threshold <- 1

significant <- (
  dge_results$FDR < fdr_threshold &
    abs(dge_results$logFC) > logfc_threshold
)

upregulated <- (
  significant & dge_results$logFC > 0
)

downregulated <- (
  significant & dge_results$logFC < 0
)

cat("\nAT2 DGE summary:\n")
cat("Genes tested:", nrow(dge_results), "\n")
cat("Genes with FDR < 0.05:",
    sum(dge_results$FDR < 0.05), "\n")
cat("Significant genes (FDR < 0.05, |logFC| > 1):",
    sum(significant), "\n")
cat("Upregulated in COVID-19:", sum(upregulated), "\n")
cat("Downregulated in COVID-19:", sum(downregulated), "\n")


# -----------------------------------------------------------------------------
# 8.5. Inspect top genes
# -----------------------------------------------------------------------------

cat("\nTop 10 AT2 differentially expressed genes:\n")
print(head(dge_results, 10))


# -----------------------------------------------------------------------------
# 8.6. Save results
# -----------------------------------------------------------------------------

write.csv(
  dge_results,
  file = file.path(
    project_dir,
    "data/results/at2_dge_edger.csv"
  ),
  row.names = TRUE
)


# =============================================================================
# 9. Verify AT2 and Macrophage Pseudobulk Input Files
# =============================================================================

# Load original macrophage pseudobulk counts.
macrophage_df <- read.csv(
  file.path(
    project_dir,
    "data/results/macrophage_pseudobulk_counts.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

# Load original AT2 pseudobulk counts.
at2_df <- read.csv(
  file.path(
    project_dir,
    "data/results/at2_pseudobulk_counts.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

# Check dimensions.
cat("Macrophage dimensions:", dim(macrophage_df), "\n")
cat("AT2 dimensions:", dim(at2_df), "\n")

# Check whether both matrices are identical.
cat(
  "Are the count matrices identical?",
  identical(macrophage_df, at2_df),
  "\n"
)

# Compare counts for one donor.
cat("\nTotal raw counts for C51ctr:\n")
cat("Macrophage:", sum(macrophage_df["C51ctr", ]), "\n")
cat("AT2:", sum(at2_df["C51ctr", ]), "\n")


# =============================================================================
# 9. AT2 Differential Expression Diagnostics
# =============================================================================

# -----------------------------------------------------------------------------
# 9.1. Define output directory
# -----------------------------------------------------------------------------

figures_dir <- file.path(
  project_dir,
  "figures",
  "at2_dge"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 9.2. Biological coefficient of variation (BCV) plot
# -----------------------------------------------------------------------------

png(
  filename = file.path(figures_dir, "at2_bcv.png"),
  width = 1000,
  height = 750,
  res = 120
)

plotBCV(dge)

dev.off()


# -----------------------------------------------------------------------------
# 9.3. MDS plot of AT2 pseudobulk samples
# -----------------------------------------------------------------------------

# Assign colors according to experimental group.
group_colors <- ifelse(
  clinical_adjusted$group == "COVID-19",
  "firebrick",
  "steelblue"
)

# Calculate MDS coordinates.
mds <- plotMDS(dge, plot = FALSE)

png(
  filename = file.path(figures_dir, "at2_mds.png"),
  width = 1100,
  height = 850,
  res = 120
)

plot(
  mds$x,
  mds$y,
  col = group_colors,
  pch = 19,
  cex = 1.2,
  xlab = "MDS Dimension 1",
  ylab = "MDS Dimension 2",
  main = "AT2 Pseudobulk MDS"
)

text(
  x = mds$x,
  y = mds$y,
  labels = clinical_adjusted$donor_id,
  pos = 3,
  cex = 0.65
)

legend(
  "topright",
  legend = c("Control", "COVID-19"),
  col = c("steelblue", "firebrick"),
  pch = 19,
  bty = "n"
)

dev.off()


# -----------------------------------------------------------------------------
# 9.4. Export MDS coordinates
# -----------------------------------------------------------------------------

mds_coordinates <- data.frame(
  donor_id = clinical_adjusted$donor_id,
  group = clinical_adjusted$group,
  n_nuclei = qc$n_nuclei[
    match(clinical_adjusted$donor_id, qc$donor_id)
  ],
  MDS1 = mds$x,
  MDS2 = mds$y
)

write.csv(
  mds_coordinates,
  file = file.path(
    project_dir,
    "data/results/at2_mds_coordinates.csv"
  ),
  row.names = FALSE
)

cat("\nAT2 diagnostic plots saved to:\n")
print(figures_dir)


# -----------------------------------------------------------------------------
# 9.5. Display diagnostic plots in RStudio
# -----------------------------------------------------------------------------

# Run these directly in the console to see the plots.
plotBCV(dge)

plot(
  mds$x,
  mds$y,
  col = group_colors,
  pch = 19,
  xlab = "MDS Dimension 1",
  ylab = "MDS Dimension 2",
  main = "AT2 Pseudobulk MDS"
)

text(
  mds$x,
  mds$y,
  labels = clinical_adjusted$donor_id,
  pos = 3,
  cex = 0.65
)



# =============================================================================
# 10. AT2 Sensitivity Analysis — Excluding Low-Nucleus Donors
# =============================================================================

# -----------------------------------------------------------------------------
# 10.1. Exclude L05cov and L11cov
# -----------------------------------------------------------------------------

exclude_donors <- c("L05cov", "L11cov")

clinical_sensitivity <- droplevels(
  clinical_adjusted[
    !clinical_adjusted$donor_id %in% exclude_donors,
    ,
    drop = FALSE
  ]
)

counts_sensitivity <- counts_adjusted[
  ,
  clinical_sensitivity$donor_id,
  drop = FALSE
]

design_sensitivity <- model.matrix(
  ~ age + sex + group,
  data = clinical_sensitivity
)

stopifnot(
  identical(
    colnames(counts_sensitivity),
    clinical_sensitivity$donor_id
  ),
  qr(design_sensitivity)$rank == ncol(design_sensitivity)
)

cat("\nSensitivity analysis donors:\n")
print(table(clinical_sensitivity$group))


# -----------------------------------------------------------------------------
# 10.2. Filter genes and normalize
# -----------------------------------------------------------------------------

dge_sensitivity <- DGEList(
  counts = counts_sensitivity,
  group = clinical_sensitivity$group
)

keep_sensitivity <- filterByExpr(
  dge_sensitivity,
  design = design_sensitivity
)

dge_sensitivity <- dge_sensitivity[
  keep_sensitivity,
  ,
  keep.lib.sizes = FALSE
]

dge_sensitivity <- calcNormFactors(
  dge_sensitivity,
  method = "TMM"
)


# -----------------------------------------------------------------------------
# 10.3. Fit edgeR model
# -----------------------------------------------------------------------------

dge_sensitivity <- estimateDisp(
  dge_sensitivity,
  design = design_sensitivity,
  robust = TRUE
)

fit_sensitivity <- glmQLFit(
  dge_sensitivity,
  design = design_sensitivity,
  robust = TRUE
)

qlf_sensitivity <- glmQLFTest(
  fit_sensitivity,
  coef = "groupCOVID-19"
)

results_sensitivity <- topTags(
  qlf_sensitivity,
  n = Inf,
  adjust.method = "BH"
)$table


# -----------------------------------------------------------------------------
# 10.4. Compare significant gene counts
# -----------------------------------------------------------------------------

count_significant <- function(results) {
  sum(
    results$FDR < 0.05 &
      abs(results$logFC) > 1
  )
}

cat("\nPrimary AT2 analysis:\n")
cat("Donors:", ncol(dge), "\n")
cat("Significant genes:",
    count_significant(dge_results), "\n")

cat("\nSensitivity analysis:\n")
cat("Donors:", ncol(dge_sensitivity), "\n")
cat("Significant genes:",
    count_significant(results_sensitivity), "\n")


# -----------------------------------------------------------------------------
# 10.5. Compare fold changes for shared genes
# -----------------------------------------------------------------------------

shared_genes <- intersect(
  rownames(dge_results),
  rownames(results_sensitivity)
)

comparison <- data.frame(
  gene = shared_genes,
  logFC_primary = dge_results[shared_genes, "logFC"],
  logFC_sensitivity = results_sensitivity[shared_genes, "logFC"],
  FDR_primary = dge_results[shared_genes, "FDR"],
  FDR_sensitivity = results_sensitivity[shared_genes, "FDR"]
)

fc_correlation <- cor(
  comparison$logFC_primary,
  comparison$logFC_sensitivity
)

cat("\nShared genes:", length(shared_genes), "\n")
cat("Log2FC correlation:", fc_correlation, "\n")


# -----------------------------------------------------------------------------
# 10.6. Compare top 10 genes
# -----------------------------------------------------------------------------

top_genes <- head(rownames(dge_results), 10)

top_gene_comparison <- comparison[
  match(top_genes, comparison$gene),
  ,
  drop = FALSE
]

cat("\nTop 10 AT2 genes — primary vs sensitivity:\n")
print(top_gene_comparison)


# -----------------------------------------------------------------------------
# 10.7. Save sensitivity results
# -----------------------------------------------------------------------------

write.csv(
  results_sensitivity,
  file = file.path(
    project_dir,
    "data/results/at2_dge_without_L05_L11.csv"
  ),
  row.names = TRUE
)

write.csv(
  comparison,
  file = file.path(
    project_dir,
    "data/results/at2_dge_sensitivity_comparison.csv"
  ),
  row.names = FALSE
)

