
# =============================================================================
# 06_pseudobulk_dge.R
# Pseudobulk Differential Gene Expression Analysis
# COVID-19 vs Control | Macrophages
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Install and load required packages
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}

library(edgeR)


# -----------------------------------------------------------------------------
# 2. Define project paths
# -----------------------------------------------------------------------------

project_dir <- "G:/My Drive/Shimon/stuffs/Data_Science_New_way/Portfolio/Self_projects/Project_6_SC_lung_atlas_COVID"

counts_path <- file.path(
  project_dir,
  "data/results/macrophage_pseudobulk_counts.csv"
)


# -----------------------------------------------------------------------------
# 3. Load pseudobulk count matrix
# -----------------------------------------------------------------------------

counts_df <- read.csv(
  counts_path,
  row.names = 1,
  check.names = FALSE
)

# Python exports donors as rows and genes as columns.
# edgeR requires genes as rows and donors as columns.

counts <- t(as.matrix(counts_df))


# -----------------------------------------------------------------------------
# 4. Verify count matrix
# -----------------------------------------------------------------------------

# Number of genes and donors
print(dim(counts))

# Donor identifiers
print(colnames(counts))

# Preview the first 5 genes and 5 donors
print(counts[1:5, 1:5])

# Verify that all values are valid raw counts
stopifnot(
  is.numeric(counts),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts))
)


# -----------------------------------------------------------------------------
# 5. Load and match donor metadata
# -----------------------------------------------------------------------------

metadata_path <- file.path(
  project_dir,
  "data/results/macrophage_pseudobulk_qc.csv"
)

metadata <- read.csv(
  metadata_path,
  stringsAsFactors = FALSE
)

# Inspect the available metadata
print(names(metadata))
print(head(metadata))

# Check that every donor appears exactly once
stopifnot(
  "donor_id" %in% names(metadata),
  !anyDuplicated(metadata$donor_id)
)

# Match metadata rows to count matrix columns
metadata <- metadata[
  match(colnames(counts), metadata$donor_id),
  ,
  drop = FALSE
]

# Verify that all donors matched
stopifnot(
  !anyNA(metadata$donor_id),
  identical(metadata$donor_id, colnames(counts))
)

# Inspect experimental groups
print(table(metadata$group))


print(names(metadata))


# -----------------------------------------------------------------------------
# 6. Load and validate clinical metadata
# -----------------------------------------------------------------------------

clinical_path <- file.path(
  project_dir,
  "data/results/donor_metadata.csv"
)

clinical <- read.csv(
  clinical_path,
  stringsAsFactors = FALSE
)

# Verify unique donor IDs
stopifnot(!anyDuplicated(clinical$donor_id))

# Match clinical metadata to the count matrix
clinical <- clinical[
  match(colnames(counts), clinical$donor_id),
  ,
  drop = FALSE
]

# Verify donor and group correspondence
stopifnot(
  !anyNA(clinical$donor_id),
  identical(clinical$donor_id, colnames(counts)),
  identical(
    as.character(clinical$group),
    as.character(metadata$group)
  )
)

# Convert variables to appropriate types
clinical$age <- as.numeric(clinical$age)

clinical$sex <- factor(clinical$sex)

clinical$group <- factor(
  clinical$group,
  levels = c("Control", "COVID-19")
)

# Inspect missing values
print(colSums(is.na(clinical)))

# Identify donors with complete clinical information
complete_samples <- complete.cases(
  clinical[, c("age", "sex", "group")]
)

cat("\nDonors available for adjusted model:",
    sum(complete_samples), "\n")

cat("\nDonors excluded from adjusted model:\n")
print(clinical$donor_id[!complete_samples])

# Construct design matrix for adjusted model
clinical_adjusted <- droplevels(
  clinical[complete_samples, , drop = FALSE]
)

design <- model.matrix(
  ~ age + sex + group,
  data = clinical_adjusted
)

# Verify that the design matrix has full rank
stopifnot(qr(design)$rank == ncol(design))

print(colnames(design))
print(dim(design))


# =============================================================================
# 7. Gene Filtering and TMM Normalization
# =============================================================================

# -----------------------------------------------------------------------------
# 7.1. Select donors for the adjusted model
# -----------------------------------------------------------------------------

# Select the 25 donors with complete age and sex metadata.
# L04cov is excluded from this model because its age is missing.
# The original count matrix remains unchanged.

counts_adjusted <- counts[
  ,
  clinical_adjusted$donor_id,
  drop = FALSE
]

# Verify that count columns match metadata rows
stopifnot(
  identical(colnames(counts_adjusted),
            clinical_adjusted$donor_id),
  ncol(counts_adjusted) == nrow(design)
)


# -----------------------------------------------------------------------------
# 7.2. Create edgeR DGEList object
# -----------------------------------------------------------------------------

# DGEList stores raw counts, sample information, and normalization factors.

dge <- DGEList(
  counts = counts_adjusted,
  group = clinical_adjusted$group
)


# -----------------------------------------------------------------------------
# 7.3. Filter low-expression genes
# -----------------------------------------------------------------------------

# filterByExpr identifies genes with sufficient expression
# for reliable statistical testing, accounting for the design matrix.

keep <- filterByExpr(
  dge,
  design = design
)

cat("\nGenes before filtering:", nrow(dge), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")

# Remove low-expression genes and recalculate library sizes
dge <- dge[keep, , keep.lib.sizes = FALSE]


# -----------------------------------------------------------------------------
# 7.4. Apply TMM normalization
# -----------------------------------------------------------------------------

# TMM (Trimmed Mean of M-values) adjusts effective library sizes
# to account for differences in RNA composition between samples.
# Raw gene counts are not transformed.

dge <- calcNormFactors(
  dge,
  method = "TMM"
)


# -----------------------------------------------------------------------------
# 7.5. Inspect normalization results
# -----------------------------------------------------------------------------

# Display sample information, library sizes, and normalization factors
print(dge$samples)

# Display final number of genes and donors
cat("\nFinal DGEList dimensions:\n")
print(dim(dge))

# Calculate effective library sizes
effective_lib_sizes <- dge$samples$lib.size *
  dge$samples$norm.factors

cat("\nEffective library sizes:\n")
print(effective_lib_sizes)

cat("\nEffective library sizes:\n")
print(effective_lib_sizes)

# -----------------------------------------------------------------------------
# 7.6. Export TMM-normalized log2 CPM for Python heatmaps
# -----------------------------------------------------------------------------

# 'dge' contains the 25 donors and 7,742 filtered genes,
# with TMM normalization factors calculated in Section 7.4.

log_cpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)

# Transpose: donors as rows, genes as columns.
log_cpm_df <- as.data.frame(t(log_cpm))

# Save to the same results directory used by Python.
output_path <- file.path(
  project_dir,
  "data/results/macrophage_log2cpm.csv"
)

write.csv(
  log_cpm_df,
  file = output_path,
  row.names = TRUE
)

cat("\nNormalized expression dimensions:\n")
print(dim(log_cpm_df))

cat("\nSaved to:", output_path, "\n")
# =============================================================================
# 8. Dispersion Estimation and Differential Expression Testing
# =============================================================================

# -----------------------------------------------------------------------------
# 8.1. Estimate gene-wise dispersion
# -----------------------------------------------------------------------------

# Estimate biological variability between donors.
# The design matrix adjusts for age, sex, and COVID-19 status.

dge <- estimateDisp(
  dge,
  design = design,
  robust = TRUE
)

# Inspect dispersion estimates
cat("\nCommon dispersion:\n")
print(dge$common.dispersion)

cat("\nBiological coefficient of variation (BCV):\n")
print(sqrt(dge$common.dispersion))


# -----------------------------------------------------------------------------
# 8.2. Fit negative-binomial quasi-likelihood model
# -----------------------------------------------------------------------------

# Fit the adjusted model for every retained gene.

fit <- glmQLFit(
  dge,
  design = design,
  robust = TRUE
)


# -----------------------------------------------------------------------------
# 8.3. Test COVID-19 versus Control
# -----------------------------------------------------------------------------

# Test the COVID-19 coefficient while adjusting for age and sex.

qlf <- glmQLFTest(
  fit,
  coef = "groupCOVID-19"
)


# -----------------------------------------------------------------------------
# 8.4. Extract differential expression results
# -----------------------------------------------------------------------------

# Extract all tested genes, ordered by statistical significance.
# FDR uses Benjamini-Hochberg multiple-testing correction.

dge_results <- topTags(
  qlf,
  n = Inf,
  adjust.method = "BH"
)$table

# Preview the most significant genes
cat("\nTop 10 differentially expressed genes:\n")
print(head(dge_results, 10))


# -----------------------------------------------------------------------------
# 8.5. Summarize significant genes
# -----------------------------------------------------------------------------

# Define significance thresholds:
# FDR < 0.05 and absolute log2 fold change > 1.

significant <- (
  dge_results$FDR < 0.05 &
    abs(dge_results$logFC) > 1
)

upregulated <- significant & dge_results$logFC > 1
downregulated <- significant & dge_results$logFC < -1

cat("\nDifferential expression summary:\n")
cat("Total genes tested:", nrow(dge_results), "\n")
cat("Significant genes:", sum(significant), "\n")
cat("Upregulated in COVID-19:", sum(upregulated), "\n")
cat("Downregulated in COVID-19:", sum(downregulated), "\n")


# -----------------------------------------------------------------------------
# 8.6. Save differential expression results
# -----------------------------------------------------------------------------

output_path <- file.path(
  project_dir,
  "data/results/macrophage_dge_edger.csv"
)

write.csv(
  dge_results,
  file = output_path,
  row.names = TRUE
)

cat("\nResults saved to:\n", output_path, "\n")





# Common dispersion
print(dge$common.dispersion)

# Top 10 differentially expressed genes
print(head(dge_results, 10))

# Number of significant genes (FDR < 0.05)
cat("FDR < 0.05:", sum(dge_results$FDR < 0.05), "\n")

# Number with FDR < 0.05 and |log2FC| > 1
cat(
  "Significant upregulated:",
  sum(dge_results$FDR < 0.05 & dge_results$logFC > 1),
  "\n"
)

cat(
  "Significant downregulated:",
  sum(dge_results$FDR < 0.05 & dge_results$logFC < -1),
  "\n"
)



# =============================================================================
# 9. Macrophage Differential Expression Diagnostics
# =============================================================================

# -----------------------------------------------------------------------------
# 9.1. Define output directory
# -----------------------------------------------------------------------------

figures_dir <- file.path(
  project_dir,
  "figures",
  "macrophage_dge"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 9.2. Biological coefficient of variation (BCV) plot
# -----------------------------------------------------------------------------

# The BCV plot shows how biological variability changes
# with average gene expression.
#
# Higher BCV indicates greater variability between donors.
# The common and trended dispersion estimates are displayed.

png(
  filename = file.path(figures_dir, "macrophage_bcv.png"),
  width = 1000,
  height = 750,
  res = 120
)

plotBCV(dge)

dev.off()


# -----------------------------------------------------------------------------
# 9.3. MDS plot of pseudobulk samples
# -----------------------------------------------------------------------------

# Multidimensional scaling (MDS) visualizes differences
# in gene expression profiles between donors.
#
# Samples close together have more similar expression profiles.
# This helps identify potential outliers and group separation.

group_colors <- ifelse(
  clinical_adjusted$group == "COVID-19",
  "firebrick",
  "steelblue"
)

png(
  filename = file.path(figures_dir, "macrophage_mds.png"),
  width = 1100,
  height = 850,
  res = 120
)

plotMDS(
  dge,
  col = group_colors,
  pch = 19,
  cex = 1.2,
  main = "Macrophage Pseudobulk MDS"
)

# Label donors
mds <- plotMDS(dge, plot = FALSE)

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
  MDS1 = mds$x,
  MDS2 = mds$y
)

write.csv(
  mds_coordinates,
  file = file.path(
    project_dir,
    "data/results/macrophage_mds_coordinates.csv"
  ),
  row.names = FALSE
)

cat("\nMacrophage diagnostic plots saved to:\n")
print(figures_dir)


# =============================================================================
# 10. Sensitivity Analysis — Excluding L01cov
# =============================================================================

# -----------------------------------------------------------------------------
# 10.1. Remove L01cov from the analysis
# -----------------------------------------------------------------------------

# L01cov was identified as an outlier in the MDS plot.
# We exclude it temporarily to evaluate its influence on DGE results.
# The original analysis remains unchanged.

exclude_donor <- "L01cov"

clinical_sensitivity <- droplevels(
  clinical_adjusted[
    clinical_adjusted$donor_id != exclude_donor,
    ,
    drop = FALSE
  ]
)

counts_sensitivity <- counts_adjusted[
  ,
  clinical_sensitivity$donor_id,
  drop = FALSE
]

# Rebuild the design matrix
design_sensitivity <- model.matrix(
  ~ age + sex + group,
  data = clinical_sensitivity
)

# Verify sample alignment and model rank
stopifnot(
  identical(
    colnames(counts_sensitivity),
    clinical_sensitivity$donor_id
  ),
  ncol(counts_sensitivity) == nrow(design_sensitivity),
  qr(design_sensitivity)$rank == ncol(design_sensitivity)
)

cat("\nDonors in sensitivity analysis:",
    ncol(counts_sensitivity), "\n")


# -----------------------------------------------------------------------------
# 10.2. Repeat gene filtering and normalization
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
# 10.3. Repeat dispersion estimation and DGE testing
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

cat("\nOriginal analysis:\n")
cat("Donors:", ncol(dge), "\n")
cat("Significant genes:", count_significant(dge_results), "\n")

cat("\nSensitivity analysis (without L01cov):\n")
cat("Donors:", ncol(dge_sensitivity), "\n")
cat("Significant genes:", count_significant(results_sensitivity), "\n")


# -----------------------------------------------------------------------------
# 10.5. Compare fold changes for shared genes
# -----------------------------------------------------------------------------

# Genes are matched by gene name, not by row position.
# Some genes may be filtered out in one analysis.

shared_genes <- intersect(
  rownames(dge_results),
  rownames(results_sensitivity)
)

comparison <- data.frame(
  gene = shared_genes,
  logFC_original = dge_results[shared_genes, "logFC"],
  logFC_sensitivity = results_sensitivity[shared_genes, "logFC"],
  FDR_original = dge_results[shared_genes, "FDR"],
  FDR_sensitivity = results_sensitivity[shared_genes, "FDR"]
)

# Calculate correlation between fold changes
fc_correlation <- cor(
  comparison$logFC_original,
  comparison$logFC_sensitivity,
  method = "pearson"
)

cat("\nShared genes:", length(shared_genes), "\n")
cat("Log2 fold-change correlation:", fc_correlation, "\n")


# -----------------------------------------------------------------------------
# 10.6. Examine the original top 10 genes
# -----------------------------------------------------------------------------

top_genes <- head(rownames(dge_results), 10)

top_gene_comparison <- comparison[
  match(top_genes, comparison$gene),
  ,
  drop = FALSE
]

cat("\nTop 10 genes — original vs sensitivity analysis:\n")
print(top_gene_comparison)


# -----------------------------------------------------------------------------
# 10.7. Save sensitivity analysis results
# -----------------------------------------------------------------------------

write.csv(
  results_sensitivity,
  file = file.path(
    project_dir,
    "data/results/macrophage_dge_without_L01cov.csv"
  ),
  row.names = TRUE
)

write.csv(
  comparison,
  file = file.path(
    project_dir,
    "data/results/macrophage_dge_sensitivity_comparison.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. Unadjusted Differential Expression Analysis — All 26 Donors
# =============================================================================

# -----------------------------------------------------------------------------
# 11.1. Construct the unadjusted design matrix
# -----------------------------------------------------------------------------

# Use all 26 donors, including L04cov.
# Only COVID-19 status is included as a predictor.
# Control is the reference group.

clinical_unadjusted <- clinical

design_unadjusted <- model.matrix(
  ~ group,
  data = clinical_unadjusted
)

# Verify donor alignment and design matrix rank
stopifnot(
  identical(colnames(counts), clinical_unadjusted$donor_id),
  nrow(design_unadjusted) == ncol(counts),
  qr(design_unadjusted)$rank == ncol(design_unadjusted)
)

cat("\nDonors in unadjusted model:", ncol(counts), "\n")
print(table(clinical_unadjusted$group))
print(colnames(design_unadjusted))


# -----------------------------------------------------------------------------
# 11.2. Gene filtering and TMM normalization
# -----------------------------------------------------------------------------

# Recalculate filtering and normalization using all 26 donors.

dge_unadjusted <- DGEList(
  counts = counts,
  group = clinical_unadjusted$group
)

keep_unadjusted <- filterByExpr(
  dge_unadjusted,
  design = design_unadjusted
)

dge_unadjusted <- dge_unadjusted[
  keep_unadjusted,
  ,
  keep.lib.sizes = FALSE
]

dge_unadjusted <- calcNormFactors(
  dge_unadjusted,
  method = "TMM"
)

cat("\nGenes retained in unadjusted model:",
    nrow(dge_unadjusted), "\n")


# -----------------------------------------------------------------------------
# 11.3. Estimate dispersion and fit the model
# -----------------------------------------------------------------------------

dge_unadjusted <- estimateDisp(
  dge_unadjusted,
  design = design_unadjusted,
  robust = TRUE
)

fit_unadjusted <- glmQLFit(
  dge_unadjusted,
  design = design_unadjusted,
  robust = TRUE
)


# -----------------------------------------------------------------------------
# 11.4. Test COVID-19 versus Control
# -----------------------------------------------------------------------------

qlf_unadjusted <- glmQLFTest(
  fit_unadjusted,
  coef = "groupCOVID-19"
)

results_unadjusted <- topTags(
  qlf_unadjusted,
  n = Inf,
  adjust.method = "BH"
)$table


# -----------------------------------------------------------------------------
# 11.5. Compare significant gene counts
# -----------------------------------------------------------------------------

# Use the same significance criteria in both analyses:
# FDR < 0.05 and |log2FC| > 1.

cat("\nAdjusted model (25 donors):\n")
cat("Significant genes:",
    count_significant(dge_results), "\n")

cat("\nUnadjusted model (26 donors):\n")
cat("Significant genes:",
    count_significant(results_unadjusted), "\n")


# -----------------------------------------------------------------------------
# 11.6. Compare fold changes for shared genes
# -----------------------------------------------------------------------------

# Gene filtering may differ between the models.
# Compare only genes tested in both analyses.

shared_genes_unadjusted <- intersect(
  rownames(dge_results),
  rownames(results_unadjusted)
)

comparison_unadjusted <- data.frame(
  gene = shared_genes_unadjusted,
  logFC_adjusted = dge_results[
    shared_genes_unadjusted, "logFC"
  ],
  logFC_unadjusted = results_unadjusted[
    shared_genes_unadjusted, "logFC"
  ],
  FDR_adjusted = dge_results[
    shared_genes_unadjusted, "FDR"
  ],
  FDR_unadjusted = results_unadjusted[
    shared_genes_unadjusted, "FDR"
  ]
)

fc_correlation_unadjusted <- cor(
  comparison_unadjusted$logFC_adjusted,
  comparison_unadjusted$logFC_unadjusted,
  method = "pearson"
)

cat("\nShared genes:",
    length(shared_genes_unadjusted), "\n")

cat("Fold-change correlation:",
    fc_correlation_unadjusted, "\n")


# -----------------------------------------------------------------------------
# 11.7. Compare the original top 10 genes
# -----------------------------------------------------------------------------

top_genes <- head(rownames(dge_results), 10)

top_gene_comparison_unadjusted <- comparison_unadjusted[
  match(top_genes, comparison_unadjusted$gene),
  ,
  drop = FALSE
]

cat("\nTop 10 genes — adjusted vs unadjusted:\n")
print(top_gene_comparison_unadjusted)


# -----------------------------------------------------------------------------
# 11.8. Save results
# -----------------------------------------------------------------------------

write.csv(
  results_unadjusted,
  file = file.path(
    project_dir,
    "data/results/macrophage_dge_unadjusted_26donors.csv"
  ),
  row.names = TRUE
)

write.csv(
  comparison_unadjusted,
  file = file.path(
    project_dir,
    "data/results/macrophage_dge_adjusted_vs_unadjusted.csv"
  ),
  row.names = FALSE
)

cat("\nUnadjusted sensitivity analysis completed.\n")

# =============================================================================
# 11.9. Unadjusted Model Using the Primary Analysis Gene Set
# =============================================================================

# Use the same genes tested in the adjusted model.
primary_genes <- rownames(dge_results)

# Select all 26 donors and the primary gene set.
counts_common <- counts[
  primary_genes,
  ,
  drop = FALSE
]

# Create a new edgeR object.
dge_common <- DGEList(
  counts = counts_common,
  group = clinical_unadjusted$group
)

# Recalculate TMM normalization.
dge_common <- calcNormFactors(
  dge_common,
  method = "TMM"
)

# Estimate dispersion.
dge_common <- estimateDisp(
  dge_common,
  design = design_unadjusted,
  robust = TRUE
)

# Fit the unadjusted model.
fit_common <- glmQLFit(
  dge_common,
  design = design_unadjusted,
  robust = TRUE
)

# Test COVID-19 versus Control.
qlf_common <- glmQLFTest(
  fit_common,
  coef = "groupCOVID-19"
)

# Extract results.
results_common <- topTags(
  qlf_common,
  n = Inf,
  adjust.method = "BH"
)$table

# Compare all original genes.
comparison_common <- data.frame(
  gene = primary_genes,
  logFC_adjusted = dge_results[primary_genes, "logFC"],
  logFC_unadjusted = results_common[primary_genes, "logFC"],
  FDR_adjusted = dge_results[primary_genes, "FDR"],
  FDR_unadjusted = results_common[primary_genes, "FDR"]
)

# Calculate fold-change correlation.
correlation_common <- cor(
  comparison_common$logFC_adjusted,
  comparison_common$logFC_unadjusted
)

cat("\nGenes compared:", nrow(comparison_common), "\n")
cat("Fold-change correlation:", correlation_common, "\n")

cat("Adjusted significant genes:",
    count_significant(dge_results), "\n")

cat("Unadjusted significant genes (common gene set):",
    count_significant(results_common), "\n")

# Examine the original top 10 genes.
top_genes <- head(rownames(dge_results), 10)

print(
  comparison_common[
    match(top_genes, comparison_common$gene),
    ,
    drop = FALSE
  ]
)

# Save results.
write.csv(
  results_common,
  file = file.path(
    project_dir,
    "data/results/macrophage_dge_unadjusted_common_genes.csv"
  ),
  row.names = TRUE
)




