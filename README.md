# 🧬 Single-Cell RNA-seq Analysis with Seurat

A reproducible **single-cell RNA sequencing (scRNA-seq) analysis workflow in R using Seurat**, covering data preprocessing, quality control, normalization, dimensionality reduction, clustering, marker-gene identification, and visualization.

This repository is designed as a practical workflow for exploring **cellular heterogeneity and transcriptional states at single-cell resolution**.

---

## 🔬 Project Overview

Single-cell RNA sequencing allows researchers to measure gene expression at the level of individual cells, making it possible to identify distinct cell populations and characterize their molecular states.

This project demonstrates a typical scRNA-seq analysis workflow using **R and Seurat**.

### Workflow

```text
scRNA-seq Data
      ↓
Data Import
      ↓
Quality Control
      ↓
Filtering Low-Quality Cells
      ↓
Normalization
      ↓
Highly Variable Genes
      ↓
Scaling
      ↓
PCA
      ↓
Neighbor Graph
      ↓
Clustering
      ↓
UMAP / t-SNE
      ↓
Marker Gene Identification
      ↓
Cell-Type Annotation
      ↓
Biological Interpretation
```

---

## 🧪 Analysis Pipeline

### 1. Data Import

The workflow begins by loading single-cell gene-expression data into a Seurat object.

Seurat provides a framework for storing expression matrices together with cell-level metadata and analysis results.

---

### 2. Quality Control

Quality control is performed to identify and remove low-quality cells.

Common QC metrics include:

* Number of detected genes per cell
* Total RNA counts per cell
* Percentage of mitochondrial gene expression
* Identification of potential low-quality cells

Example visualization:

```r
VlnPlot(
  object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
```

---

### 3. Normalization

Gene expression counts are normalized to reduce differences caused by sequencing depth between cells.

```r
object <- NormalizeData(object)
```

---

### 4. Highly Variable Genes

Highly variable genes are identified to focus downstream analyses on genes that capture biological variation between cells.

```r
object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000
)
```

---

### 5. Scaling

Expression values are scaled before dimensionality reduction.

```r
object <- ScaleData(object)
```

---

### 6. Principal Component Analysis

PCA is used to reduce the dimensionality of the gene-expression matrix while retaining major sources of variation.

```r
object <- RunPCA(
  object,
  features = VariableFeatures(object)
)
```

PCA results can then be visualized to investigate the structure of the dataset.

---

### 7. Cell Clustering

