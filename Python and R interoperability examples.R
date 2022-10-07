library(reticulate)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
ad <- import("anndata", convert = FALSE)
tmpDir="C:/Users/psergeev/Downloads/outs"
seurat_data <-Read10X(data.dir = file.path(tmpDir, "filtered_gene_bc_matrices", "GRCh38"))
seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 250)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet( seurat_obj, pattern = "^MT-")
seurat_obj$log10GenesPerUMI <- log10( seurat_obj$nFeature_RNA) / log10( seurat_obj$nCount_RNA)
seurat_obj=subset( seurat_obj, subset = log10GenesPerUMI > 0.8 & nFeature_RNA > 200)
seurat_obj=subset( seurat_obj, subset = percent.mt < ifelse(median(seurat_obj$percent.mt)*3 < 20,median(seurat_obj$percent.mt)*3, 20))
#remove 5% top outliers
seurat_obj=subset(seurat_obj, subset =nFeature_RNA <sort(seurat_obj$nFeature_RNA)[round(length(seurat_obj$nFeature_RNA)*0.95)])
counts <- GetAssayData(object = seurat_obj, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= as.numeric(summary(Matrix::rowSums(nonzero))[2])
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
seurat_obj <- CreateSeuratObject(filtered_counts, meta.data = seurat_obj@meta.data, project = "pmbc")
seurat_obj <- NormalizeData(seurat_obj, verbose = TRUE, scale.factor = 1000000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj<- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj<- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj<- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj<- RunUMAP(seurat_obj, reduction = "pca", dims = 1:10) #Define here the number of wanted PC's in dim= ! 

#Example of convertion of Seurat (R) object to AnnData (Python)
SaveH5Seurat(seurat_obj, overwrite = TRUE)
Convert("SeuratProject.h5Seurat", dest = "h5ad", overwrite = T)

sc=import("scanpy", convert = FALSE)
adata = sc$read_h5ad("SeuratProject.h5ad")
adata=sc$read_10x_mtx(file.path(tmpDir, "filtered_gene_bc_matrices", "GRCh38"))

sc$pl$violin(adata, c("log10GenesPerUMI", "percent_ribo", "percent_hb"),
             jitter=0.4, multi_panel=T)

sc$pl$umap(adata, color=c('CST3', 'NKG7', 'PPBP'))