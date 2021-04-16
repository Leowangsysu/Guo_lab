library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
library(DoubletFinder)
library(parallel)
library(stringr)
library(ggsci)



#--------DoubletFinder---------------

sample_name.data <- Read10X(data.dir = "/filtered_feature_bc_matrix/")
sample_name <- CreateSeuratObject(counts = sample_name.data, project = "sample_name", min.cells = 3, min.features = 200)
sample_name[["percent.mt"]] <- PercentageFeatureSet(sample_name, pattern = "^MT-")
sample_name <- subset(sample_name, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & percent.mt < 25)
sample_name <- NormalizeData(sample_name)
sample_name <- FindVariableFeatures(sample_name, selection.method = "vst", nfeatures = 2000)
all.genes_sample_name <- rownames(sample_name)
sample_name <- ScaleData(sample_name, features = all.genes_sample_name)
sample_name <- RunPCA(sample_name, features = VariableFeatures(object = sample_name))
sample_name <- FindNeighbors(sample_name, dims = 1:15)
sample_name <- FindClusters(sample_name, resolution = 0.7)
sample_name <- RunUMAP(sample_name, dims = 1:15,umap.method = "uwot")


# pK Identification ---------------------------------------------------------------------------------------------------------
sweep.res.list_sample_name <- paramSweep_v3(sample_name, PCs = 1:15)
sweep.stats_sample_name <- summarizeSweep(sweep.res.list_sample_name, GT = FALSE)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------

homotypic.prop <- modelHomotypic(sample_name@active.ident)   
nExp_poi <- round(0.075*length(colnames(sample_name))) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sample_name <- doubletFinder_v3(sample_name, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
sample_name <- doubletFinder_v3(sample_name, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_443")

## Plot results --------------------------------------------------------------------------------------------------------------
sample_name@meta.data[,"DF_hi.lo"] <- sample_name@meta.data$DF.classifications_0.25_0.09_443
sample_name@meta.data$DF_hi.lo[which(sample_name@meta.data$DF_hi.lo == "Doublet" & sample_name@meta.data$DF.classifications_0.25_0.09_395 == "Singlet")] <- "Doublet_lo"
sample_name@meta.data$DF_hi.lo[which(sample_name@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

png(("/Doubletfinder_sample_name.png"), height = 15, width = 15, res = 400, units = "in")
DimPlot(object = sample_name, reduction = 'umap', pt.size = 3.5,
        label = FALSE, cols = NULL, group.by = "DF_hi.lo")
dev.off()
meta.data <- as.data.frame(sample_name@meta.data)
doublet_cellset_sample_name <- meta.data[which(!meta.data$DF_hi.lo == "Singlet"),] 
doublet_rm_cellset_sample_name <- rownames(doublet_cellset_sample_name)
write.table(doublet_rm_cellset_sample_name,
            file = "rm_Doubletlist.txt",
            sep = ",",
            row.names = T,
            quote = F)




