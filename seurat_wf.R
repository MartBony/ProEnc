library(dplyr)
library(Seurat)
library(patchwork)


#Directory of the 3 files with names barcodes, genes and matrix
dir = "/dir/to/data/" # contains the 3 source files
setwd(dir)

# Load dataset
data <- Read10X(data.dir = dir)

# Verify the data
#data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# HHLA2 <- data[c("HHLA2"),]
# TMIGD2 <- data[c("TMIGD2"),]
# plot(TMIGD2)
# plot(HHLA2)

# Sample down the matrix to make lighter on our computers for tests
dosample <- FALSE
if(dosample){
  data <- data[,sample(1:ncol(data), ncol(data)%/%10) ]
}


# Initialize Seurat object. Ignore "Avis"
SeurOBJ <- CreateSeuratObject(counts = data$`Gene Expression`, project = "TumeurTest", min.cells = 3, min.features = 200)

# SaveSeuratRds(SeurOBJ, "seurobj.rds")

SeurOBJ <- LoadSeuratRds("seurobj.rds")

# Check genes that were kept (Search for HHLA2)
genes <- data.frame(Features(SeurOBJ))







#Add mitochondrial RNA in samples
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SeurOBJ[["percent.mt"]] <- PercentageFeatureSet(SeurOBJ, pattern = "^MT-")

# # Visualize QC metrics as a violin plot
# VlnPlot(SeurOBJ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# Ncount_MT.plot <- FeatureScatter(SeurOBJ, feature1 = "nCount_RNA", feature2 = "percent.mt")
# Ncount_NFreature.plot <- FeatureScatter(SeurOBJ, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Ncount_MT.plot + Ncount_NFreature.plot




# % d'ARN mitochondrial maximal accepté
LimitMitoRNA = 10


SeurOBJ <- subset(SeurOBJ, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < LimitMitoRNA)
# VlnPlot(SeurOBJ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 

# normalize
SeurOBJ <- NormalizeData(SeurOBJ)







SeurOBJ <- FindVariableFeatures(SeurOBJ, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(SeurOBJ), 10)
# 
# # plot variable features with and without labels
# VarFeature.plot <- VariableFeaturePlot(SeurOBJ)
# LP.plot <- LabelPoints(plot = VarFeature.plot, points = top10, repel = TRUE)
# LP.plot




# Scale
all.genes <- rownames(SeurOBJ)
SeurOBJ <- ScaleData(SeurOBJ, features = all.genes)




# PCA on Seurat Object
SeurOBJ <- RunPCA(SeurOBJ, features = VariableFeatures(object = SeurOBJ))


# # Visualise most important coordinates of the two first principal components
# VizDimLoadings(SeurOBJ, dims = 1:2, reduction = "pca")
# 
# # Visualise first factorial axis
# DimPlot(SeurOBJ, reduction = "pca") + NoLegend()
# 
# # Visualise the variance dependding on the number of principle components
# ElbowPlot(SeurOBJ)


# Clusterize
SeurOBJ <- FindNeighbors(SeurOBJ, dims = 1:10)
SeurOBJ <- FindClusters(SeurOBJ, resolution = 1)

# UMAP
SeurOBJ <- RunUMAP(SeurOBJ, dims = 1:10)

# Plot clusters
DimPlot(SeurOBJ, reduction = "umap")


# Plot Mitococo RNA amount per cluster. Lower is better
# VlnPlot(SeurOBJ, features = "percent.mt") # Pas de cluster mauvais, tous autour de la même valeur


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# ⚠️ HUGE SPEEDUP WITH 👇
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')

SeurOBJ.markers <- FindAllMarkers(SeurOBJ, only.pos = TRUE)
# Table of most unique genes per cluster
diff.expressed.genes <- SeurOBJ.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# Plot expression of these per clusters
#VlnPlot(SeurOBJ, features = c("HHLA2"))
#VlnPlot(SeurOBJ, features = c("TMIGD2")) # CD28H


#FeaturePlot(SeurOBJ, features = c("HHLA2"))
#FeaturePlot(SeurOBJ, features = c("TMIGD2"))


source("../pre_annot/pre_annot.R", chdir=TRUE)

diff.expressed.genes <- mark_knowns(diff.expressed.genes)# optionnal 

type.annot.matrix <- get_annot_matrix(SeurOBJ, diff.expressed.genes)
type.avg.matrix <- get_avg_matrix(SeurOBJ, diff.expressed.genes)
display_heatmap(type.annot.matrix, "Nombre de gènes") + display_heatmap(type.avg.matrix, "Expression différentielle")


# Display annotations on UMAP
clusters.annot <- pre_labels(type.annot.matrix, seuil = 2)
clusters.annot

names(clusters.annot) <- levels(SeurOBJ)
LabeledSeurOBJ.alt <- RenameIdents(SeurOBJ, clusters.annot)
DimPlot(LabeledSeurOBJ.alt, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()

