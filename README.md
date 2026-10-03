# 🧬 Single-Cell RNA-Sequencing Analysis with Seurat

### A practical R-based workflow for preprocessing, quality control, dimensionality reduction, clustering, visualization, and marker-gene analysis of single-cell RNA-seq data.

[![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r\&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-scRNA--seq-4B275F)](https://satijalab.org/seurat/)
[![Status](https://img.shields.io/badge/Status-Analysis%20Workflow-blue)]()
[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-green)]()

---

## 🔬 Overview

Single-cell RNA sequencing (scRNA-seq) enables researchers to investigate gene expression at the level of individual cells rather than measuring an average expression profile across a mixed population.

This repository contains an **R/Seurat-based single-cell RNA-seq analysis workflow** covering the major preprocessing and exploratory analysis steps:

```text
10X Genomics Data
       │
       ▼
Create Seurat Object
       │
       ▼
Quality Control
       │
       ▼
Cell Filtering
       │
       ▼
Normalization
       │
       ▼
Highly Variable Genes
       │
       ▼
Data Scaling
       │
       ▼
PCA
       │
       ▼
Elbow Plot
       │
       ▼
Nearest-Neighbor Graph
       │
       ▼
Cell Clustering
       │
       ▼
UMAP
       │
       ▼
Cluster Visualization
       │
       ▼
Marker Gene Identification
       │
       ▼
Basic Cell-Type Annotation
       │
       ▼
Save Seurat Object
```

---

# 🎯 Project Objectives

The purpose of this project is to demonstrate a complete introductory-to-intermediate scRNA-seq analysis workflow using **Seurat in R**.

The implemented workflow focuses on:

* Loading single-cell expression data
* Creating Seurat objects
* Performing quality control
* Filtering low-quality cells
* Normalizing gene expression
* Identifying highly variable genes
* Scaling expression data
* Performing PCA
* Evaluating principal components with an elbow plot
* Constructing a cell-neighborhood graph
* Clustering cells
* Visualizing cell populations using UMAP
* Identifying cluster-specific marker genes
* Performing basic cell-type annotation
* Saving the processed Seurat object

---

# 🧪 Implemented Analysis

## 1. Load Required Libraries

The workflow uses:

```r
library(Seurat)
library(dplyr)
library(ggplot2)
```

### Main packages

| Package     | Purpose                      |
| ----------- | ---------------------------- |
| **Seurat**  | Single-cell RNA-seq analysis |
| **dplyr**   | Data manipulation            |
| **ggplot2** | Data visualization           |

---

# 2. Load scRNA-seq Data

The repository contains workflows for loading both:

### 10X Genomics data

```r
data <- Read10X(
    data.dir = "your_data_file"
)
```

### Count matrix data

```r
counts <- read.csv(
    "path/to/counts_matrix.csv",
    row.names = 1
)
```

The resulting expression matrix is used to create a Seurat object.

---

# 3. Create Seurat Object

The expression data are converted into a Seurat object:

```r
seurat_object <- CreateSeuratObject(
    counts = data,
    project = "SingleCellRNASeq"
)
```

The Seurat object provides a structured container for:

* Gene-expression data
* Cell metadata
* Normalized data
* Variable features
* Dimensionality reductions
* Cluster assignments
* Downstream analysis results

---

# 4. Quality Control

The repository includes quality-control processing using common single-cell metrics.

### Mitochondrial gene percentage

```r
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
)
```

### QC filtering

The current workflow applies:

```r
seurat_obj <- subset(
    seurat_obj,
    subset =
        nFeature_RNA > 200 &
        nFeature_RNA < 2500 &
        percent.mt < 5
)
```

This filters cells based on:

* Minimum detected genes
* Maximum detected genes
* Mitochondrial gene percentage

> **Note:** These thresholds are parameters used in the current script and should be evaluated according to the specific dataset. They are not universal scRNA-seq thresholds.

---

# 5. Normalization

The workflow uses Seurat's standard normalization procedure:

```r
seurat_obj <- NormalizeData(seurat_obj)
```

Normalization helps account for differences in sequencing depth between cells.

---

# 6. Highly Variable Gene Identification

Highly variable genes are identified using:

```r
seurat_obj <- FindVariableFeatures(
    seurat_obj
)
```

These genes capture important variation across cells and are subsequently used in downstream analysis.

---

# 7. Data Scaling

The expression data are scaled using:

```r
seurat_obj <- ScaleData(
    seurat_obj
)
```

Scaling prepares the expression matrix for dimensionality-reduction analysis.

---

# 8. Principal Component Analysis

PCA is performed to reduce the dimensionality of the gene-expression dataset:

```r
seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj)
)
```

PCA transforms the high-dimensional expression space into a smaller number of principal components representing major sources of variation.

---

# 9. PCA Elbow Plot

The workflow includes an elbow plot to help evaluate the number of principal components:

```r
ElbowPlot(seurat_obj)
```

The elbow plot can help inform the selection of PCs for downstream neighborhood construction and clustering.

---

# 10. Cell-Cell Neighborhood Graph

A nearest-neighbor graph is constructed using the selected principal components:

```r
seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = 1:10
)
```

The current workflow uses the first **10 principal components**.

---

# 11. Cell Clustering

Cells are grouped into transcriptionally similar populations:

```r
seurat_obj <- FindClusters(
    seurat_obj,
    resolution = 0.5
)
```

The current implementation uses:

```text
Dimensions: 1–10 PCs
Resolution: 0.5
```

The resulting clusters represent groups of cells with similar gene-expression profiles.

---

# 12. UMAP Dimensionality Reduction

UMAP is performed using the selected principal components:

```r
seurat_obj <- RunUMAP(
    seurat_obj,
    dims = 1:10
)
```

UMAP provides a two-dimensional representation of the transcriptional structure of the dataset.

---

# 13. UMAP Cluster Visualization

Clusters are visualized using:

```r
DimPlot(
    seurat_obj,
    reduction = "umap",
    label = TRUE
)
```

The visualization displays:

* Individual cells
* Cluster structure
* Cluster labels
* Transcriptional relationships

The workflow also includes a titled UMAP visualization:

```r
DimPlot(
    seurat_object,
    reduction = "umap",
    label = TRUE
) +
ggtitle(
    "UMAP Plot of Single-Cell RNA-seq Data"
)
```

---

# 14. Cluster Marker Identification

The workflow identifies genes associated with individual clusters:

```r
cluster_markers <- FindAllMarkers(
    seurat_obj
)
```

This allows exploration of genes that distinguish different cell clusters.

Marker genes can subsequently be examined to help characterize cellular populations.

---

# 15. Marker Gene Visualization

The workflow includes `FeaturePlot()` for visualizing expression of selected genes:

```r
FeaturePlot(
    seurat_obj,
    features = c(
        "GeneA",
        "GeneB"
    )
)
```

`GeneA` and `GeneB` are currently placeholders and should be replaced with genes relevant to the dataset.

For example:

```r
FeaturePlot(
    seurat_obj,
    features = c(
        "CD3D",
        "MS4A1"
    )
)
```

would visualize the expression of example immune-cell markers if those genes are present in the dataset.

---

# 16. Basic Cell-Type Annotation

The repository contains a basic example of assigning cell types based on cluster IDs:

```r
seurat_obj$cell_type <- ifelse(
    seurat_obj$seurat_clusters == 0,
    "CellTypeA",
    "CellTypeB"
)
```

This is a **demonstration framework**, not a biological annotation of a specific dataset.

In a real analysis, cluster identities should be assigned using:

* Marker-gene expression
* Known biological markers
* Tissue context
* Experimental conditions
* Reference-based annotation where appropriate

---

# 17. Save the Seurat Object

The processed Seurat object can be saved for later analysis:

```r
saveRDS(
    seurat_obj,
    file = "seurat_obj.rds"
)
```

This allows the analysis to be continued without repeating the complete preprocessing workflow.

---

# 📊 Analysis Summary

| Analysis                   | Status |
| -------------------------- | :----: |
| 10X data import            |    ✅   |
| CSV count-matrix import    |    ✅   |
| Seurat object creation     |    ✅   |
| Mitochondrial QC           |    ✅   |
| Cell filtering             |    ✅   |
| Normalization              |    ✅   |
| Highly variable genes      |    ✅   |
| Scaling                    |    ✅   |
| PCA                        |    ✅   |
| Elbow plot                 |    ✅   |
| Nearest-neighbor graph     |    ✅   |
| Graph-based clustering     |    ✅   |
| UMAP                       |    ✅   |
| Cluster visualization      |    ✅   |
| Marker-gene identification |    ✅   |
| FeaturePlot                |    ✅   |
| Basic annotation framework |    ✅   |
| Save Seurat object         |    ✅   |

---

# 🧬 Current Workflow

```text
             scRNA-seq Data
                   │
          ┌────────┴────────┐
          │                 │
        10X             Count Matrix
          │                 │
          └────────┬────────┘
                   ▼
           Seurat Object
                   │
                   ▼
            Quality Control
                   │
                   ▼
             Cell Filtering
                   │
                   ▼
             Normalization
                   │
                   ▼
        Highly Variable Genes
                   │
                   ▼
                Scaling
                   │
                   ▼
                 PCA
                   │
                   ▼
             Elbow Plot
                   │
                   ▼
          Nearest Neighbors
                   │
                   ▼
              Clustering
                   │
                   ▼
                 UMAP
                   │
                   ▼
          Cluster Visualization
                   │
                   ▼
           Marker Identification
                   │
                   ▼
           Marker Visualization
                   │
                   ▼
        Basic Cell-Type Annotation
                   │
                   ▼
          Save Seurat Object
```

---

# 📁 Repository Structure

```text
single-Cell-RNA-Sequencing/
│
├── README.md
│
├── code.r
│
└── script.r
```

### `code.r`

Contains the workflow for:

* Loading 10X data
* Creating a Seurat object
* Normalization
* Variable feature detection
* Scaling
* PCA
* Clustering
* UMAP
* Cluster visualization

### `script.r`

Contains a more extended workflow including:

* Data import
* Quality control
* Mitochondrial percentage calculation
* Cell filtering
* Normalization
* Variable feature detection
* Scaling
* PCA
* Elbow plot
* Clustering
* UMAP
* Marker-gene identification
* Feature visualization
* Basic annotation
* Saving the Seurat object

---

# ⚙️ Requirements

The workflow requires:

* **R**
* **Seurat**
* **dplyr**
* **ggplot2**

Install packages with:

```r
install.packages("Seurat")
install.packages("dplyr")
install.packages("ggplot2")
```

Then load them:

```r
library(Seurat)
library(dplyr)
library(ggplot2)
```

---

# 🚀 Getting Started

Clone the repository:

```bash
git clone https://github.com/Bioinformatician-dev/single-Cell-RNA-Sequencing.git
```

Move into the repository:

```bash
cd single-Cell-RNA-Sequencing
```

Provide the path to your 10X dataset:

```r
data <- Read10X(
    data.dir = "path/to/your/10X/data"
)
```

Or provide a count matrix:

```r
counts <- read.csv(
    "path/to/counts_matrix.csv",
    row.names = 1
)
```

Then run the appropriate R script.

---

# ⚠️ Important Notes

### Dataset path

The scripts currently contain placeholder paths such as:

```r
"your_data_file"
```

and:

```r
"path/to/counts_matrix.csv"
```

These must be replaced with the location of the actual dataset.

### Marker genes

The example:

```r
c("GeneA", "GeneB")
```

is a placeholder and should be replaced with genes present in your dataset.

### Cell-type annotation

The current `CellTypeA` / `CellTypeB` assignment is an example of the annotation mechanism and should not be interpreted as validated biological annotation.

### QC thresholds

The current filtering parameters are:

```text
nFeature_RNA > 200
nFeature_RNA < 2500
percent.mt < 5
```

These values should be evaluated according to the characteristics of the dataset.

---

# 🔭 Future Extensions

The current repository establishes the core Seurat workflow. Possible future extensions include:

* [ ] SCTransform normalization
* [ ] Improved QC visualization
* [ ] Doublet detection
* [ ] Automated cell-type annotation
* [ ] Multi-sample integration
* [ ] Differential expression between experimental conditions
* [ ] GO enrichment
* [ ] KEGG/pathway analysis
* [ ] Cell-cell communication analysis
* [ ] Pseudotime / trajectory analysis
* [ ] Copy-number variation analysis
* [ ] scRNA-seq + scATAC-seq integration
* [ ] Spatial transcriptomics integration
* [ ] Interactive visualization

> These are **planned extensions**, not analyses currently implemented in this repository.

---

# 🧠 Biological Questions This Workflow Can Address

Once applied to an appropriate dataset, the workflow can help investigate:

* What cellular populations are present?
* Which cells have similar transcriptional profiles?
* Which genes distinguish clusters?
* Which marker genes characterize individual populations?
* What potential cell types correspond to identified clusters?
* How heterogeneous is the sampled tissue?

---

# 📚 Main Technologies

```text
R
│
├── Seurat
│   ├── scRNA-seq preprocessing
│   ├── PCA
│   ├── Clustering
│   ├── UMAP
│   └── Marker analysis
│
├── dplyr
│   └── Data manipulation
│
└── ggplot2
    └── Visualization
```

---

# 👩‍🔬 Author

## Salma Hafeez

**Bioinformatician | Computational Biology | Genomics | Transcriptomics**

GitHub:
**[Bioinformatician-dev](https://github.com/Bioinformatician-dev)**

Research interests include:

`Single-Cell RNA-seq` • `Transcriptomics` • `Cancer Genomics` • `Computational Biology` • `Machine Learning`

---

# ⭐ Project Status

**Status:** 🟢 Core scRNA-seq workflow implemented

The repository currently demonstrates the fundamental Seurat pipeline from expression-data loading through clustering, UMAP visualization, marker-gene identification, and basic annotation.

---

## 🧬 From Single Cells to Biological Insight

**Load → QC → Normalize → Reduce → Cluster → Visualize → Discover**
