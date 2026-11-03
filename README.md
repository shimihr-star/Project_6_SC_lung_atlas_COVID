# COVID-19 Lung Atlas — Single-Nucleus RNA-seq Analysis

Analysis of the publicly available [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524) lung atlas, examining donor-level cell-type composition and transcriptional changes associated with fatal COVID-19, with a focus on **macrophages** and **alveolar type 2 (AT2) epithelial cells**.

![Project pipeline overview](figures/pipeline_overview.png)

> **Figure note:** This infographic summarizes the completed workflow and selected numerical findings. Embedded plots and icons are schematic illustrations, not actual computed UMAPs, volcano plots, or enrichment figures. Numerical results are reported from the project analyses.

## Dataset

| Property | Description |
|---|---|
| Source | GEO [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524) |
| Donors | 26 total: 19 fatal COVID-19 cases and 7 controls |
| Integrated dataset | 96,915 nuclei × 31,969 genes |
| Modality | Single-nucleus RNA sequencing (snRNA-seq) |
| Biological replication | Donors, **not individual nuclei** |
| Adjusted DGE cohort | 25 donors (18 COVID-19, 7 controls); one COVID-19 donor excluded because age was missing |

## Analysis pipeline

| Step | Analysis | Status |
|---|---|---|
| 01 | Data acquisition, nucleus/gene QC, scVI/SOLO doublet detection | Complete |
| 02 | scVI integration, neighborhood graph, Leiden clustering, UMAP | Complete |
| 03 | Marker-gene analysis and manual cell-type annotation | Complete; annotations provisional |
| 04 | Donor-level cell-type composition, Mann–Whitney tests, BH-FDR correction, age/sex-adjusted analysis | Complete |
| 05 | Donor × cell-type pseudobulk aggregation, QC, edgeR differential expression, sensitivity analyses | Complete |
| 06 | Gene set enrichment analysis (MSigDB Hallmark gene sets) and cross-cell-type comparison | Complete |

### Methodological overview

**Preprocessing and integration.** Raw UMI counts were loaded into AnnData; nucleus QC and SOLO-based doublet filtering were applied. scVI was trained on highly variable genes with `biosample_id` as the batch covariate. The learned latent space was used for neighborhood graph construction, UMAP, and Leiden clustering. Marker-gene identification used a **separately normalized copy** of the AnnData object; the count matrix used downstream for pseudobulk was retained without normalization.

**Composition analysis.** Nuclei were counted by donor and annotated cell type. Donor-level percentages were compared using two-sided Mann–Whitney U tests with Benjamini–Hochberg correction. Additional age/sex-adjusted regression was performed for selected populations.

**Differential expression.** Raw UMI counts were summed within each donor × cell-type combination. Macrophage and AT2 pseudobulk counts were analyzed separately with edgeR using expression filtering, TMM normalization, robust dispersion estimation, and quasi-likelihood generalized linear models. The primary design was `~ age + sex + group`, with COVID-19 contrasted against control. The principal descriptive threshold was **FDR < 0.05 and |log2FC| > 1**. Sensitivity analyses assessed influential macrophage donors and low-nucleus-count AT2 donors.

**Functional enrichment.** Genes were ranked using the signed square root of the edgeR quasi-likelihood F statistic and evaluated with GSEA against 50 MSigDB Hallmark gene sets (1,000 permutations). A positive normalized enrichment score (NES) denotes enrichment toward genes associated with higher expression in COVID-19.

## Results

### 1. Donor-level cell-type composition

| Cell type | COVID-19 median | Control median | BH-adjusted p-value |
|---|---:|---:|---:|
| Macrophage | 25.283% | 6.356% | 0.0180 |
| AT2 | 7.432% | 27.201% | 0.0215 |
| Endothelial-like — unclassified | 2.017% | 1.123% | 0.0215 |

Macrophages accounted for a larger median share of recovered nuclei in COVID-19 lungs, while AT2 nuclei represented a smaller share. The endothelial-like label is provisional. **These are relative proportions of recovered nuclei, not absolute tissue cell abundances.**

### 2. Pseudobulk sample QC

Pseudobulk samples were defined by **donor × cell type**. All 26 donors contributed to the initial sample QC; the primary adjusted DGE model included the 25 donors with complete age data.

| Cell type | COVID-19 donors (QC) | Control donors (QC) | Median nuclei per COVID-19 donor | Median nuclei per control donor |
|---|---:|---:|---:|---:|
| Macrophage | 19 | 7 | 695 | 254 |
| AT2 | 19 | 7 | 208 | 1,033 |

Two COVID-19 donors had only 23 AT2 nuclei each; sensitivity analysis assessed their influence. Library size, detected genes, and exploratory PCA/MDS were examined as sample-level diagnostics.

### 3. Differential gene expression

| Cell type | Genes tested | Significant genes | Higher in COVID-19 | Lower in COVID-19 |
|---|---:|---:|---:|---:|
| Macrophage | 7,742 | **818** | 455 | 363 |
| AT2 | 9,188 | **1,580** | 891 | 689 |

Significance is reported at FDR < 0.05 and |log2FC| > 1. This is a post-testing effect-size filter, **not** a formal statistical test that the true absolute fold change exceeds twofold.

Selected macrophage-associated changes include `CEBPB`, `PKM`, and `JUND` (higher in COVID-19) and `CHIT1` and `NLRP3` (lower). In AT2, higher expression of `HLA-C`, `TAPBP`, `APOL6`, and `UNC93B1` is consistent with an immune-associated transcriptional signature. These gene-level observations are interpreted in conjunction with pathway results rather than as evidence of direct causal mechanisms.

The number of significant genes differs across cell types, but this **does not** establish that one population is more affected: expression levels, variance, nuclei recovered, and statistical power differ.

### 4. Hallmark pathway enrichment

Among 50 Hallmark gene sets, **22 were significant (FDR < 0.05) in at least one cell type**:

- **11 shared** between macrophages and AT2, including interferon-alpha and interferon-gamma responses, oxidative phosphorylation, glycolysis, mTORC1 signaling, reactive oxygen species, and p53 signaling.
- **7 significant only in macrophages**, including complement, adipogenesis/lipid-related signatures, protein secretion, and TGF-beta signaling.
- **4 significant only in AT2**, including allograft rejection, cholesterol homeostasis, and apoptosis-associated signatures.

| Selected shared Hallmark pathway | Macrophage NES | AT2 NES |
|---|---:|---:|
| Interferon alpha response | +1.77 | +2.37 |
| Interferon gamma response | +1.60 | +2.01 |
| Oxidative phosphorylation | +2.43 | +1.71 |
| Glycolysis | +1.73 | +1.71 |
| mTORC1 signaling | +2.33 | +2.09 |
| Reactive oxygen species | +1.97 | +1.65 |
| p53 pathway | +1.73 | +2.02 |

Notably, the AT2 TNF-alpha/NF-kappaB Hallmark signature showed **negative enrichment** (NES approximately −1.85), placing this gene set toward the control-associated end of the ranked list.

**Important interpretation:** A pathway significant in one population but not another is not, by itself, evidence of a statistically significant cell-type interaction. GSEA measures coordinated transcriptional associations, not direct protein activity or pathway flux.

## Biological interpretation

The analyses indicate **shared immune and metabolic transcriptional remodeling** in macrophage and AT2 populations in fatal COVID-19 lung tissue. Interferon-related and metabolic Hallmark signatures are prominent in both populations. Macrophages additionally show complement and lipid-metabolism-associated enrichment, whereas AT2 cells show antigen-presentation-related gene changes and distinct epithelial-associated pathway signatures.

These findings describe **associations in postmortem tissue**. They do not establish that viral infection directly causes each change or that the identified transcriptional programs represent functional pathway activation.

## Limitations

- **Observational, postmortem comparison:** disease severity, treatment, tissue injury, and sample collection may contribute to observed differences.
- **Small and unbalanced donor groups:** 19 fatal COVID-19 cases versus 7 controls; the primary adjusted DGE model uses 18 COVID-19 donors because age was unavailable for one donor.
- **Potential residual confounding:** age and sex were included in adjusted models, but other factors (including postmortem interval and tissue processing) may influence results.
- **Provisional cell-type labels:** marker-based annotations are not independent experimental validation. `Macrophage` and `Alveolar macrophage` are distinct annotated populations; DGE here refers specifically to `Macrophage`.
- **Relative cell composition:** percentages of recovered nuclei cannot be interpreted as absolute cell abundances.
- **Sensitivity to sample recovery:** two COVID-19 AT2 pseudobulk samples contain only 23 nuclei each; sensitivity analyses were performed.
- **Pathway interpretation:** enrichment differences between cell types are descriptive unless formally tested as interactions.

## Repository structure

```text
Project_6_SC_lung_atlas_COVID/
├── README.md
├── requirements.txt
├── .gitignore
├── scripts/
│   ├── 01_preprocess.py
│   ├── 02_integrate_cluster.py
│   ├── 03_annotate.py
│   ├── 04_cell_composition.py
│   ├── 05_pseudobulk.py
│   ├── 06_pseudobulk_dge.R
│   ├── 07_at2_pseudobulk_dge.R
│   ├── 08_dge_vizualization.py
│   └── 09_pathway_enrichment.py
├── src/
│   ├── config.py
│   ├── preprocessing.py
│   ├── integration.py
│   ├── annotation.py
│   ├── composition.py
│   └── pseudobulk.py
├── figures/
│   └── pipeline_overview.png
└── data/
    ├── raw_counts/       # not committed
    ├── processed/        # not committed
    └── results/          # generated analysis outputs
```

> File names above reflect the working analysis scripts; confirm the exact names in your local repository before committing. Other figures and generated result tables may also be present.

## Software and reproducibility

Python analyses use Scanpy, scVI-tools, PyTorch, pandas, NumPy, SciPy, scikit-learn, Matplotlib, and related libraries. Differential expression uses R and **edgeR**; GSEA uses **gseapy**. The working environment includes Scanpy 1.12.4, scVI-tools 1.5.1, and PyTorch 2.11 with CUDA 12.8. Export the exact Python and R package environments from the machines used to execute the analyses rather than treating this README as a lockfile.

Scripts are numbered by stage, with reusable functions in `src/`. Large raw datasets, processed `.h5ad` objects, and trained model checkpoints are excluded from Git. Original data are publicly available from [GEO GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524).

## Project status

**Core analysis complete.** Preprocessing, integration, cell-type annotation, donor-level composition, pseudobulk differential expression, sensitivity analyses, and Hallmark GSEA have been performed. Cell-type labels remain provisional, and biological conclusions are framed as observational associations.
