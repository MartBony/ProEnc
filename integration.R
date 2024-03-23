# GSE154826 : CITEseq analysis of non-small-cell lung cancer lesions
# Load the datasets
library(harmony)
library(Seurat)
library(sctransform) # replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
library(dplyr)

do.sample <- FALSE
proj.name <- "lungAM"
dir = "./raw_data/" # Contains a list of files 
setwd(dir)

source("./utils.R")

merged.dir <- "./merged_seurat.rds"
images.dir <- "./out/"

print(paste("Integrating ", merged.dir))

if(file.exists(merged.dir)){
  merged_seurat <- LoadSeuratRds(merged.dir)
  
} else {
  batch_numbers <- c(1,2,3,4,5,6,7,8,9,10,11,14,15,16,17,18,19,22,23,24,25,26,27,28,
                     29,30,31,32,33,36,37,38,39,40,41,42,43,44,45,46,47,49,50,
                     94,95,230,231,244,245,292,293,307,308,309,342,343,344,345,346,347,348,350,351,668,
                     # not a matrix but list
                     48,51,52,53,54,55,56,88,89,90,91,92,93)
  
  # batch_numbers <- c(1,345, 53) #test
  samples.data <- read.csv("./samples_data.csv")
  
  seurobj_vec <- c()
  
  print("Importing")
  # Import by reading raw data or through file
  # Could use for(batch_num in batch_numbers){with append} but this is for progress bar
  seurobj_vec <- apply.to.all(batch_numbers, function(batch_num){
    seurobj_dir <- paste("./raw_data/seurobj_batch",batch_num,".rds", sep="")
    if(file.exists(seurobj_dir)){
      seurobj <- LoadSeuratRds(seurobj_dir)
      
    } else{
      print('')
      print(paste('Importing batch :', batch_num))
      data_dir <- paste("./raw_data/batch",batch_num, sep="")
      
      raw_data <- Read10X(data.dir = data_dir)
      if(class(raw_data) == "list"){ # 10X data sometimes conta other matrices, the first one is the count matrix
        raw_data <- raw_data[[1]]
        print("10X list data, extracted count matrix")
      }
      if(do.sample){
        raw_data <- raw_data[,sample(1:ncol(raw_data), ncol(raw_data)%/%6) ]
      }
      seurobj <- CreateSeuratObject(raw_data, project=paste(proj.name,".",batch_num, sep=""),
                                    min.cells = 3, 
                                    min.features = 200)
      
    
      SaveSeuratRds(seurobj, seurobj_dir)
      
    }

    #Add tissue metadata
    seurobj[[]]$tissue <- samples.data[samples.data$amp_batch_ID==batch_num,]$tissue[1] 
    seurobj[[]]$patient_ID <- samples.data[samples.data$amp_batch_ID==batch_num,]$patient_ID[1]

    return(seurobj)
    
  })

  print("Calculating HHLA2+ patients")
  patients.with.HHLA2 <- c()
  patients.with.HHLA2 <- apply.to.all(seurobj_vec, function(seurobj){
    if(any(row.names(seurobj[["RNA"]]$counts) == 'HHLA2') && sum(seurobj[["RNA"]]$counts["HHLA2",]) > 50){
      print(sum(seurobj[["RNA"]]$counts["HHLA2",]))
      return(seurobj[[]]$patient_ID[[1]])
    }
    return(NA)
  })

  print("Patients with HHLA2 :")
  print(patients.with.HHLA2)


  print("Filtering datasets")
  # Filter patients with no HHLA2 in any of their sample
  seurobj_vec <- apply.to.all(seurobj_vec, function(seurobj){ # Only take if one of this patient's sample expresses HHLA2
    if(seurobj$patient_ID[1] %in% patients.with.HHLA2){
      return(seurobj)
    }
    return(NA)
  })


  print("Calculating % of mitochondrial RNA")  
  ## Standard pre-processing workflow ## the goal is to remove unwanted cells from the dataset
  seurobj_vec <- apply.to.all(seurobj_vec, function(seurOBJ){
    seurOBJ[["percent.mt"]] <- PercentageFeatureSet(seurOBJ, pattern = "^MT-")
    return(seurOBJ)
  })
  
  # Do for a few
  VlnPlot(seurobj_vec[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  nFeature.min <- 200
  nFeature.max <- 3000
  percent.mt.max <- 12
  
  print("Subsetting")  
  ## Standard pre-processing workflow ## the goal is to remove unwanted cells from the dataset
  seurobj_vec <- apply.to.all(seurobj_vec, function(seurOBJ){
    seurOBJ <- subset(seurOBJ, subset = nFeature_RNA > nFeature.min & nFeature_RNA < nFeature.max & percent.mt < percent.mt.max)
    return(seurOBJ)
  })
  
  print("Running SCTransform")  
  # run sctransform
  seurobj_vec <- apply.to.all(seurobj_vec, function(seurOBJ){
    seurOBJ <- SCTransform(seurOBJ, vars.to.regress = "percent.mt", verbose=FALSE)
    return(seurOBJ)
  })
  
  
  # Find most variable features across samples to integrate
  integ_featuresallSCT <- SelectIntegrationFeatures(object.list = seurobj_vec, nfeatures = 20000) 
  print("HHLA2 %in% integ_featuresallSCT")
  print("HHLA2" %in% integ_featuresallSCT)
  print("KIR3DL3 %in% integ_featuresallSCT")
  print("KIR3DL3" %in% integ_featuresallSCT)
  print("TMIGD2 %in% integ_featuresallSCT")
  print("TMIGD2" %in% integ_featuresallSCT)
  
  # Merge normalized samples
  merged_seurat <- merge(x = seurobj_vec[[1]],
                               y = seurobj_vec[2:length(seurobj_vec)],
                               merge.data = TRUE)
  DefaultAssay(merged_seurat) <- "SCT"
  
  # Manually set variable features of merged Seurat object
  VariableFeatures(merged_seurat) <- integ_featuresallSCT

  
  # Calculate PCs using manually set variable features
  merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)
  #pb : utilise seulement le scale data T844
  
  #HARMONY
  merged_seurat <- RunHarmony(merged_seurat, 
    group.by.vars = c("orig.ident"), 
    reduction = "pca",  reduction.save = "harmony") #assay.use = "SCT"

    
  all.genes <- rownames(merged_seurat)
  print("HHLA2 %in% all.genes")
  print("HHLA2" %in% all.genes)
  print("KIR3DL3 %in% all.genes")
  print("KIR3DL3" %in% all.genes)
  print("TMIGD2 %in% all.genes")
  print("TMIGD2" %in% all.genes)
  
  
  # Usual Seurat workflow
  
  Elbow <- ElbowPlot(merged_seurat)
  save.plot.png(Elbow, paste(images.dir, "elbow.png"))

  numdims <- 16

  merged_seurat <- FindNeighbors(merged_seurat, dims = 1:numdims)
  merged_seurat <- FindClusters(merged_seurat, resolution = 0.2)
  
  
  merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", assay = all.genes, dims = 1:numdims)
    
  
  merged_seurat <- PrepSCTFindMarkers(merged_seurat, assay = "SCT")

  
  
  SaveSeuratRds(merged_seurat, merged.dir)

}  


# Visualisations
num_UMAP <- DimPlot(merged_seurat, reduction = "umap", label=TRUE, pt.size = 0.25, raster=FALSE)
save.plot.png(num_UMAP, paste(images.dir, "UMAP_numbered.png"))


manual_annot <- c("Macrophages","Lymphocytes T (bis)","Lymphocytes T","Lymphocytes T (bis2)",4
,"Cellules Tumorales Ciliées",6,"Lymphocytes B","Lymphocytes B (bis)",9
,"Natural Killer","Cellulles Endothéliales Environnantes","Granulocytes",13,"Cellules Tumorales")
get_cluster_name <- function(cluster_index){
  return(manual_annot[cluster_index+1])
}

names(manual_annot) <- levels(merged_seurat)
merged_seurat <- RenameIdents(merged_seurat, manual_annot)

UMAP_manual <- DimPlot(merged_seurat, reduction = "umap", label=TRUE, pt.size = 0.25, raster=FALSE)
save.plot.png(UMAP_manual, paste(images.dir, "UMAP_manualed.png"))


UMAP_manual_tissue <- DimPlot(merged_seurat, reduction = "umap", split.by="tissue", label=TRUE, pt.size = 0.25, raster=FALSE)
save.plot.png(UMAP_manual_tissue, paste(images.dir, "UMAP_manualed_tissue.png"), width=1200, height=600)


vln_HHLA <- VlnPlot(merged_seurat, features = c("HHLA2"))
vln_TMI <- VlnPlot(merged_seurat, features = c("TMIGD2"))
vln_KIR <- VlnPlot(merged_seurat, features = c("KIR3DL3"))
save.plot.png(vln_HHLA, paste(images.dir, "vln_HHLA.png"))
save.plot.png(vln_TMI, paste(images.dir, "vln_TMIGD2.png"))
save.plot.png(vln_KIR, paste(images.dir, "vln_KIR3DL3.png"))
  

feature_HHLA <- FeaturePlot(merged_seurat, features = c("HHLA2"), raster=FALSE) # split.by = "tissue",
feature_TMI <- FeaturePlot(merged_seurat, features = c("TMIGD2"), raster=FALSE) # split.by = "tissue",
feature_KIR <- FeaturePlot(merged_seurat, features = c("KIR3DL3"), raster=FALSE) # split.by = "tissue",
save.plot.png(feature_HHLA, paste(images.dir, "feature_HHLA.png"))
save.plot.png(feature_TMI, paste(images.dir, "feature_TMIGD2.png"))
save.plot.png(feature_KIR, paste(images.dir, "feature_KIR3DL3.png"))

dotplot <- DotPlot(merged_seurat, features = c("HHLA2", "TMIGD2", "KIR3DL3")) + RotatedAxis()
save.plot.png(dotplot, paste(images.dir, "dotplot.png"))

# Expressions différentielles pour HHLA, NK et autres
#nk.cells.id <- which(merged_seurat[[]]$seurat_clusters %in%c(11,13,15,31,35)) # Res 1.0
#nk_clusters <- c("5","12")
#nk_cells <- subset(merged_seurat, idents = nk_clusters)
HHLA_clusters <- get_cluster_name(5)
TMI_clusters <- c(get_cluster_name(1),get_cluster_name(2),get_cluster_name(3), get_cluster_name(8), get_cluster_name(10))
HHLA_cells <- subset(merged_seurat, idents = HHLA_clusters)
TMI_cells <- subset(merged_seurat, idents = TMI_clusters)

print("Expressions différentielles pour HHLA, NK et autres")

feature_tissue_HHLA <- FeaturePlot(merged_seurat, features = c("HHLA2"), split.by = "tissue", raster=FALSE)
save.plot.png(feature_tissue_HHLA, paste(images.dir, "feature_tissue_HHLA.png"), width=1200, height=600)
feature_tissue_TMI <- FeaturePlot(merged_seurat, features = c("TMIGD2"), split.by = "tissue", raster=FALSE)
save.plot.png(feature_tissue_TMI, paste(images.dir, "feature_tissue_TMI.png"), width=1200, height=600)




#nk_per_tissue <- DimPlot(nk_cells, reduction = "umap", group.by="tissue", pt.size = 0.25, raster=FALSE)
#save.plot.png(nk_per_tissue, paste(images.dir, "nk_per_tissue.png"))

#nk_tissue_markers <- FindMarkers(merged_seurat, assay = "SCT", ident.1="Tumor", group.by="tissue", subset.ident = nk_clusters,recorrect_umi = FALSE)
#write.csv2(nk_tissue_markers, paste(images.dir ,"markers_5_12_identTumor_GroupByTissue.csv"))
#saveRDS(nk_tissue_markers, paste(images.dir ,"markers_5_12_identTumor_GroupByTissue.rds"))

print("HHLA2, tumor vs normal in cluster 5 : ")
print(FindMarkers(merged_seurat, assay = "SCT", ident.1="Tumor", group.by="tissue", subset.ident=get_cluster_name(5), features="HHLA2", recorrect_umi = FALSE))
print("TMIGD2, tumor vs normal in cluster 2 : ")
print(FindMarkers(merged_seurat, assay = "SCT", ident.1="Tumor", group.by="tissue", subset.ident=get_cluster_name(2), features="TMIGD2", recorrect_umi = FALSE))
print("TMIGD2, tumor vs normal in cluster 10 : ")
print(FindMarkers(merged_seurat, assay = "SCT", ident.1="Tumor", group.by="tissue", subset.ident=get_cluster_name(10), features="TMIGD2", recorrect_umi = FALSE))


print("Calculating HHLA tissue markers")
HHLA_tissue_markers <- FindMarkers(merged_seurat, assay = "SCT", ident.1="Tumor", group.by="tissue", subset.ident = HHLA_clusters, recorrect_umi = FALSE)
write.csv2(HHLA_tissue_markers, paste(images.dir ,"markers_HHLA_identTumor_GroupByTissue.csv"))
saveRDS(HHLA_tissue_markers, paste(images.dir ,"markers_HHLA_identTumor_GroupByTissue.rds"))

print("Calculating TMIGD2 Cluster tissue markers")
nclusters <- length(levels(merged_seurat))
TMI_tissue_markers <- FindMarkers(merged_seurat, assay = "SCT", ident.1="Tumor", group.by = "tissue", subset.ident=TMI_clusters, recorrect_umi = FALSE)
write.csv2(TMI_tissue_markers, paste(images.dir ,"markers_TMI_identTumor_GroupByCluster.csv"))
saveRDS(TMI_tissue_markers, paste(images.dir ,"markers_TMI_identTumor_GroupByCluster.rds"))




# Marqueurs

print("Calculating cluster marker genes")
merged_seurat.markers <- FindAllMarkers(merged_seurat, only.pos = TRUE)
# Table of most unique genes per cluster
diff.expressed.genes <- merged_seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)






# Annotation

source("../pre_annot/pre_annot.R", chdir=TRUE)

diff.expressed.genes <- mark_knowns(diff.expressed.genes)# optionnal 
write.csv2(diff.expressed.genes,  paste(images.dir, "allmarkers_filtered_marked.csv"))
saveRDS(diff.expressed.genes, paste(images.dir ,"allmarkers_filtered_marked.rds"))


type.annot.matrix <- get_annot_matrix(merged_seurat, diff.expressed.genes)
type.avg.matrix <- get_avg_matrix(merged_seurat, diff.expressed.genes)
type.corresp.matrix <- get_corresp_matrix(merged_seurat, diff.expressed.genes)
type.modulated.matrix <- type.annot.matrix * type.avg.matrix
gene_heatmap <- display_heatmap(type.annot.matrix, "Nombre de gènes")
save.plot.png(gene_heatmap, paste(images.dir, "heatmap_gene.png"))
avg_heatmap <- display_heatmap(type.avg.matrix, "Expression différentielle")
save.plot.png(avg_heatmap, paste(images.dir, "heatmap_avg.png"))
display_heatmap(type.corresp.matrix, "% de correspondance")
modulated_heatmap <- display_heatmap(type.modulated.matrix, "Somme des log fold-change")
save.plot.png(modulated_heatmap, paste(images.dir, "heatmap_modulated.png"))



# Compute types from matrix of choice
clusters.annot <- pre_labels(type.annot.matrix, seuil = 2)
print(clusters.annot)

clusters.annot.alternatif <- pre_labels(type.modulated.matrix, seuil = 2)
clusters.annot.alternatif

names(clusters.annot) <- levels(merged_seurat)
merged_seurat.labeled <- RenameIdents(merged_seurat, clusters.annot)
gene_named_UMAP <- DimPlot(merged_seurat.labeled, reduction = "umap", label = TRUE, pt.size = 0.25) + 
  NoLegend() #+ labs(title = "Annotation par comparaison")
save.plot.png(gene_named_UMAP, paste(images.dir, "UMAP_gene_named.png"))


names(clusters.annot.alternatif) <- levels(merged_seurat)
merged_seurat.alt.labeled <- RenameIdents(merged_seurat, clusters.annot.alternatif)
modulated_named_UMAP <-DimPlot(merged_seurat.alt.labeled, reduction = "umap", label = TRUE, pt.size = 0.25) + 
  NoLegend() #+ labs(title = "Annotation par comparaison et expression")
save.plot.png(modulated_named_UMAP, paste(images.dir, "UMAP_modulated_named.png"))
