library(Seurat)
path="C:/Users/psergeev/OneDrive/OneDrive - University of Helsinki/Project stuff/AP project (Juho, Romika)/scRNA myeloma/Samples/FH_6261_5/outs/filtered_feature_bc_matrix"
set.seed(100101) #for UMAP
seurat_data <-Read10X(data.dir = path)
seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 250) 
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scRNAseq)
  library(scater)
  library(flexmix)
  library(splines)
  library(biomaRt)
  library(miQC)
})
sce=as.SingleCellExperiment(seurat_obj)
mt_genes <- grepl("^mt-",  rownames(sce), ignore.case = T)
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sce <- addPerCellQC(sce, subsets = feature_ctrls)
model <- mixtureModel(sce)
sce <- filterCells(sce, model)
seurat_obj=seurat_obj[,colnames(sce)]
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj= PercentageFeatureSet(seurat_obj, "^RP[SL]", col.name = "percent_ribo")
seurat_obj= PercentageFeatureSet(seurat_obj, "^HB[^(P)]", col.name = "percent_hb")
seurat_obj=subset( seurat_obj, subset = log10GenesPerUMI > 0.8)
counts <- GetAssayData(object = seurat_obj, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >=10
filtered_counts <- counts[keep_genes, ]
df=NULL
df["before"]=nrow(seurat_obj)
# Reassign to filtered Seurat object
seurat_obj <- CreateSeuratObject(filtered_counts, meta.data = seurat_obj@meta.data)
seurat_obj <- NormalizeData(seurat_obj, verbose = TRUE, scale.factor = 1000000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

library(ensembldb)
library(AnnotationHub)
cc_file <- RCurl::getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
library(dplyr)
# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")
seurat_obj <- CellCycleScoring(seurat_obj,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)
DimPlot(seurat_obj,
        reduction = "pca",
        group.by= "Phase")
