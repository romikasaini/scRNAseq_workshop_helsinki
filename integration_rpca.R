# This is a very rough/first draft version of integration script. Let's run/edit/explore together during hands-on session on Friday 7TH Oct 2022

library(Seurat)
library(ggplot2)
#setwd("/Users/admin_adhisadi/workspace/nemesyys_mobility/scRNAseq_workshop_helsinki") #your directory
set.seed(100101) #for UMAP


#I will be preprocessing the datasets with shortcut. Please do your preprocessing carefully. As you know this dataset has perfect quality I am just saving time here.
for (file in c("data1", "data2", "data3", "data4")) {
  seurat_data <-
    Read10X(data.dir = file)
  seurat_obj <-
    CreateSeuratObject(
      counts = seurat_data,
      min.features = 200,
      min.cells = 3,
      project = file,
      assay = "RNA"
    )
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  assign(file, seurat_obj)
}

all_data <- list(data1, data2, data3, data4)


for (i in 1:length(all_data)){
  all_data[[i]]<- subset(all_data[[i]], subset = nFeature_RNA > 200 & percent.mt < 20)
  all_data[[i]]$log10GenesPerUMI <- log10(all_data[[i]]$nFeature_RNA) / log10(all_data[[i]]$nCount_RNA)
  all_data[[i]]=subset(all_data[[i]], subset = log10GenesPerUMI > 0.8)
}

filtered_data <- all_data #we need this later for merge function

#normalization
for (i in 1:length(all_data)){
  all_data[[i]] <- NormalizeData(all_data[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

#identify highly variable features
for (i in 1:length(all_data)){
  all_data[[i]] <- FindVariableFeatures(all_data[[i]], selection.method = "vst", nfeatures = 2000)
}

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = all_data)
all_data <- lapply(X = all_data, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#INTEGRATION: This method is called Fast integration using reciprocal PCA (RPCA)
#I have used this method because it is fast (as the name suggests) and I had a very large dataset.
#But depending on your needs from the integration, you might need to use this or other methods. eg, if you have large batch effect.
#This publication, figure 5, summarizes major (all?) integration methods available so far and their plus and minus points. 
# https://www.nature.com/articles/s41592-021-01336-8

#We then identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
anchors_all_data <- FindIntegrationAnchors(object.list = all_data, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
integrate_data <- IntegrateData(anchorset = anchors_all_data)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrate_data) <- "integrated"

# Run the standard workflow for visualization and clustering
integrate_data <- ScaleData(integrate_data, verbose = FALSE)
integrate_data <- RunPCA(integrate_data, npcs = 30, verbose = FALSE)
integrate_data <- RunUMAP(integrate_data, reduction = "pca", dims = 1:30)
integrate_data <- FindNeighbors(integrate_data, reduction = "pca", dims = 1:30)
integrate_data <- FindClusters(integrate_data, resolution = 0.5)

# Visualization
p1 <- DimPlot(integrate_data, reduction = "umap", group.by = "orig.ident") 
p2 <- DimPlot(integrate_data, reduction = "umap", split.by = "orig.ident")
p2 
#Appears that there is integration, but each sample dominates some of the clusters. 
#This kind of cluster might mean biological difference or batch effect in your dataset.
#Note that rpca is not the best method for batch correction.

#But let's check how this would look if we did not integrated. We will have some idea on batch effect also if we check with merge function. 

# We will take filtered data here
for (i in 1:length(filtered_data)){
  filtered_data[[i]] <- FindVariableFeatures(filtered_data[[i]], selection.method = "vst", nfeatures = 2000)
}

for (i in 1:length(filtered_data)){
  merged_data <- Reduce(function(x,y) merge(x,y,add.cell.ids = c(x@project.name,y@project.name)) ,filtered_data)
  filtered_data[[i]] <- FindVariableFeatures(filtered_data[[i]], selection.method = "vst", nfeatures = 2000)
}
unique(Idents(merged_data)) #sample id
VlnPlot(merged_data, "nFeature_RNA", group.by = "orig.ident")

#normalize, variable genes
merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
var.genes <- SelectIntegrationFeatures(SplitObject(merged_data, split.by = "orig.ident"), nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, selection.method = "vst")
VariableFeatures(merged_data) <- var.genes
merged_data <- ScaleData(merged_data, features = VariableFeatures(merged_data))
# Do PCA on data including only the variable genes.
merged_data <- RunPCA(merged_data, features = VariableFeatures(merged_data), npcs = 20, ndims.print = 1:5, nfeatures.print = 5)
# Cluster the cells using the first ten principal components.
merged_data <- FindNeighbors(merged_data, reduction = "pca", dims = 1:10, k.param = 20)
merged_data <- FindClusters(merged_data, resolution = 0.5, algorithm = 1, random.seed = 100)
merged_data <- RunUMAP(merged_data, dims = 1:10, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

#let's check 
p6 <- DimPlot(merged_data, reduction = "umap", group.by = "orig.ident")
p6

#let's give a title so we can compare
p6 <- p6 + ggtitle("merged")
p1 <- p1 + ggtitle("rpca")
p1+p6

#let's discuss and write comments together
