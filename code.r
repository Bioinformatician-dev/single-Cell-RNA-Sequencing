# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the single-cell RNA-seq data
# Replace 'your_data_file' with the actual data file path
data <- Read10X(data.dir = "your_data_file")

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = data, project = "SingleCellRNASeq")

# Preprocess the data
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)

# Scale the data
seurat_object <- ScaleData(seurat_object)

# Perform PCA for dimensionality reduction
seurat_object <- RunPCA(seurat_object)

# Find clusters of cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Run UMAP for visualization
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Visualize the clusters
DimPlot(seurat_object, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Plot of Single-Cell RNA-seq Data")
