#https://github.com/constantAmateur/SoupX

library(SoupX)
tmpDir="C:/Users/psergeev/Downloads/outs"

#load raw and filtered datasets
toc = Seurat::Read10X(file.path(tmpDir, "filtered_gene_bc_matrices", "GRCh38"),)
tod = Seurat::Read10X(file.path(tmpDir, "raw_gene_bc_matrices", "GRCh38"))

toc=toc[,names(clus1)]
tod=tod[,names(clus1)]

#create Soup object from those datasets
sc = SoupChannel(tod, toc)
sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc = estimateSoup(sc)
toc = sc$toc
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
# Calculate soup profile
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops, soupProf)

#add cluster info from 10x pipeline
clus=read.csv(file.path(tmpDir, "analysis", "clustering", "graphclust", "clusters.csv"), row.names = 1)

adjustCounts(sc, clusters = clus)
sc = setClusters(sc, setNames(clus$Cluster, names(clus)))

#add tSNE projection info for visualization purposes
tsn=read.csv(file.path(tmpDir, "analysis", "tsne", "2_components", "projection.csv"), row.names = 1)

sc = setDR(sc, tsn[colnames(sc$toc), c("TSNE.1", "TSNE.2")])
sc = autoEstCont(sc)
# Clean the data
out = adjustCounts(sc)

#plot examples
plotMarkerDistribution(sc)
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), 
                                      clusters = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)


plotChangeMap(sc, out, "S100A9")
