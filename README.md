<!-- COVER IMAGE -->
<p align="center">
  <img src="figures/cover_fig.png" alt="PBMC3k scVI Cover" width="80%">
</p>

# scVI + Scanpy Pipeline for PBMC3k  
### Deep Generative Modeling & Classical Single-Cell RNA-seq Analysis

This repository contains a complete, reproducible, and fully documented pipeline for **single-cell RNA sequencing (scRNA-seq)** using:

- **Scanpy** — classical statistical workflow  
- **scVI-tools** — deep generative modeling using Variational Autoencoders  
- **PyTorch** — GPU-accelerated backend  
- **UMAP / Leiden** — dimensionality reduction and clustering  
- **Differential Expression (Wilcoxon)**  

This project demonstrates a full **end-to-end analysis** of the canonical **PBMC3k** dataset, covering:

- preprocessing  
- QC  
- normalization  
- highly variable gene selection  
- scVI latent space learning  
- UMAP  
- clustering  
- DE analysis  
- manual cell-type annotation  
- publication-ready figures  

The repository is structured to be a professional example suitable for bioinformatics portfolios, graduate program applications, and research reproducibility.

---

# 1. Biological Background

## What are PBMCs?

PBMCs (**Peripheral Blood Mononuclear Cells**) are immune cells with a round nucleus.  
They include:

| Cell Type | Marker Genes | Biological Function |
|----------|--------------|---------------------|
| **T cells** | CD3D, IL7R, LTB | Adaptive immunity, antigen recognition |
| **B cells** | MS4A1, CD19, CD79A | Antibody production |
| **NK cells** | GNLY, NKG7, GZMB | Cytotoxic innate immune response |
| **Monocytes** | LST1, S100A8, CST3 | Phagocytosis, antigen presentation |
| **Platelets** | PPBP, SDPR | Clotting, wound healing |

The PBMC3k dataset (10x Genomics) is the most commonly used "hello world" dataset for single-cell pipelines.

---

# 2. What is Single-Cell RNA-seq?

**scRNA-seq** measures gene expression cell-by-cell.  

It enables:

- discovering **new cell types**  
- detecting **rare populations**  
- understanding **cell states**  
- mapping **immune activation/exhaustion**  
- studying **tumor microenvironment (TME)**  

The output is a sparse **cells × genes** matrix:
Example: 2700 cells × 20,000 genes


Most entries are zero due to technical dropout.

---

# 3. What is scVI?

**scVI (Single-Cell Variational Inference)** is a deep learning method that uses a **Variational Autoencoder (VAE)** to model scRNA-seq counts.

### Why scVI?

Traditional Scanpy/Seurat assume linearity and struggle with:

- dropout noise  
- sparsity  
- overdispersion  
- batch effects  

scVI provides:

- NB/zero-inflated likelihood modeling  
- Latent representation robust to noise  
- Batch correction  
- Improved cluster separation  
- Better differential expression  

The latent space is learned in **10–30 dimensions** and then visualized with UMAP.

---

# 4. Project Structure
melanoma_scvi/
│
├── data/
│ ├── pbmc3k_raw.h5ad
│ ├── pbmc3k_processed.h5ad
│ └── processed/
│ └── pbmc3k_scvi_processed.h5ad
│
├── notebooks/
│ └── 01_preprocess_scvi.ipynb
│
├── results/
│ ├── figures/
│ │ ├── pbmc3k_umap.png
│ │ ├── pbmc3k_umap_celltype.png
│ │ ├── pbmc3k_umap_leiden_scvi.png
│ │ ├── pbmc3k_heatmap_degs.png
│ │ ├── pbmc3k_dotplot_markers.png
│ │ └── cover_fig.png
│ │
│ ├── tables/
│ │ ├── pbmc3k_rank_genes_leiden_scvi.csv
│ │ ├── pbmc3k_rank_genes_leiden_scvi_wilcoxon.csv
│ │ └── pbmc3k_celltype_annotations.csv
│ │
│ ├── adata/
│ │ ├── pbmc3k_processed.h5ad
│ │ ├── pbmc3k_annotated.h5ad
│ │ └── pbmc3k_X_scVI.csv
│ │
│ └── models/
│ └── scvi_pbmc3k/
│ └── model.pt
│
├── figs/
│ └── cover_fig.png
│
├── envs/
│ └── sc-omics.yml
│
└── README.md


---

# 5. Pipeline Overview

## 5.1 Preprocessing (Scanpy)

Steps performed:

1. Load raw counts  
2. Basic QC  
3. Filter cells/genes  
4. Normalize total counts  
5. Log1p transform  
6. Select HVGs  
7. Scale the data  

Output: `pbmc3k_processed.h5ad`

## 5.2 scVI Training

1. Setup AnnData for scVI  
2. Train SCVI model with GPU  
3. Extract latent representation  
4. Save trained model  

Output: `results/models/scvi_pbmc3k/model.pt`

## 5.3 UMAP + Leiden Clustering

- UMAP on scVI latent space  
- Neighborhood graph  
- Leiden resolution 0.5  

Outputs:
- `pbmc3k_leiden.csv`
- `pbmc3k_umap.png`
- `pbmc3k_umap_leiden_scvi.png`

## 5.4 Differential Expression

Using Wilcoxon rank-sum:
- Compute DEGs for each cluster  
- Save tables as CSV  
- Plot heatmaps  

Outputs:
- CSV: `pbmc3k_rank_genes_leiden_scvi_wilcoxon.csv`
- Figures: heatmap, scatterplots

## 5.5 Cell-Type Annotation

Using canonical marker genes:

T cells: CD3D, IL32, LTB
B cells: MS4A1, CD79A, CD79B
Monocytes: LYZ, CST3, TYROBP
NK cells: GNLY, NKG7, GZMB
Platelets: PPBP, SDPR


Dotplot and UMAP show clear clusters.

---

# 6. Key Figures (Examples)

Figures generated in the pipeline:
results/figures/
├── pbmc3k_umap.png
├── pbmc3k_umap_celltype.png
├── pbmc3k_umap_leiden_scvi.png
├── pbmc3k_dotplot_markers.png
├── pbmc3k_heatmap_degs.png
└── cover_fig.png


---

# 7. Running the Pipeline

## 7.1 Create Environment

```bash
micromamba create -f envs/sc-omics.yml
micromamba activate sc-omics
7.2 Run Notebook

Open Jupyter:
jupyter lab
Then execute:
notebooks/01_preprocess_scvi.ipynb

8. Reproducibility

All intermediate files are saved:

Raw (H5AD)

Processed (Scanpy)

Latent space (CSV)

Annotated data

scVI model weights

Figures

DE tables

This ensures full reproducibility for reviewers or students.
9. References

Gayoso et al. (2022). scvi-tools: A library for deep probabilistic models in single-cell omics.

Wolf et al. (2018). Scanpy: scalable single-cell analysis.

10x Genomics PBMC3k dataset

McInnes et al. (2018). UMAP: Uniform Manifold Approximation and Projection.
10. Author

Leonardo Ferreira Morais
Veterinary Medicine → Bioinformatics / Computational Oncology
scRNA-seq • scVI • Deep Learning • Immune Microenvironment

GitHub: https://github.com/leomoraisf

