# Load required libraries
install.packages("Seurat")  # Uncomment if Seurat is not installed
library(Seurat)

# Step 1: Load your data
# Replace 'path/to/counts_matrix' with your actual data file
# Assuming 'counts' is a matrix with genes as rows and cells as columns
counts <- read.csv("path/to/counts_matrix.csv", row.names = 1)
seurat_obj <- CreateSeuratObject(counts = counts, project = "SingleCellRNAseq")

# Step 2: Quality Control
# Calculate percentage of mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Filter cells based on QC metrics
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Step 3: Normalization
seurat_obj <- NormalizeData(seurat_obj)

# Step 4: Identify Highly Variable Features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Step 5: Scaling the Data
seurat_obj <- ScaleData(seurat_obj)

# Step 6: PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Step 7: Elbow Plot
ElbowPlot(seurat_obj)

# Step 8: Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Step 9: UMAP Visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Step 10: Identify Cluster Markers
cluster_markers <- FindAllMarkers(seurat_obj)

# Step 11: Visualize Marker Genes
FeaturePlot(seurat_obj, features = c("GeneA", "GeneB"))  # Replace with actual gene names

# Step 12: Cell Type Annotation
# Here you can use known markers or reference datasets for annotation
# For example, assigning cell types based on clusters
seurat_obj$cell_type <- ifelse(seurat_obj$seurat_clusters == 0, "CellTypeA", "CellTypeB")  # Adjust as needed

# Step 13: Save the Seurat object
saveRDS(seurat_obj, file = "seurat_obj.rds")
