library(Seurat)  
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
library(DoubletFinder)
library(parallel)
library(stringr)
library(ggsci)


load("/liver_cluster.Rdata")

subset_cells_MP <- subset(subset_cells, idents = c(idents = c()))  ##select MP  
subset_cells_MP <- SplitObject(subset_cells_MP, split.by = "stim")
subset_cells_MP <- lapply(X = subset_cells_MP, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
subset_cells_MP <- FindIntegrationAnchors(object.list = subset_cells_MP, dims = 1:20)
subset_cells_MP.combined <- IntegrateData(anchorset = subset_cells_MP, dims = 1:20)
subset_cells_MP <- subset_cells_MP.combined 

DefaultAssay(subset_cells_MP) <- "integrated"

all.genes <- rownames(subset_cells_MP)
subset_cells_MP <- ScaleData(subset_cells_MP, features = all.genes)

subset_cells_MP$orig.ident <- factor(Map(f = function(x){x[2]}, strsplit(colnames(subset_cells_MP), split = "_")) %>% unlist(), levels = c("0H", "6H", "9H"))
subset_cells_MP <- RunPCA(subset_cells_MP, features = VariableFeatures(object = subset_cells_MP))
subset_cells_MP <- FindNeighbors(subset_cells_MP, dims = 1:15)
subset_cells_MP <- FindClusters(subset_cells_MP, resolution = 0.3)
subset_cells_MP <- RunUMAP(subset_cells_MP, dims = 1:15, umap.method = "uwot")

DefaultAssay(subset_cells_MP) <- "RNA"
all.genes <- rownames(subset_cells_MP)
subset_cells_MP <- ScaleData(subset_cells_MP, features = all.genes)



result <- mclapply(as.numeric(levels(subset_cells_MP@active.ident)),
                   FUN =  function(x) {FindMarkers(subset_cells_MP, ident.1 = x, ident.2 = NULL)},
                   mc.cores = 36)
RESULT <- result

roundN <- 1
while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
  if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
    print(recalculate_clusters)
    result1 <- mclapply(recalculate_clusters,
                        FUN =  function(x) {FindMarkers(subset_cells_MP, ident.1 = x, ident.2 = NULL)},
                        mc.cores = 35)
  }
  print(roundN + 1)
  for(i in 1:length(recalculate_clusters)){
    result[c(recalculate_clusters+1)[i]] <- result1[i]
  }
}

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result))
all_markers$cluster <- rep(levels(subset_cells_MP@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
subset_cells_MP.markers <- all_markers

subset_cells_MP.markers %>% TOP_N(50) -> top50
subset_cells_MP.markers <- subset_cells_MP.markers %>% TOP_N(5000)


 



#----------------GSVA---------------------use RNA slot not integrat

library(GSEABase)
library(Biobase)
library(limma)
library(RColorBrewer)
library(GSVA)
library(dplyr)
hallmark <- getGmt("h.all.v7.0.symbols.gmt")
c2.PID <- getGmt("c2.cp.pid.v7.1.symbols.gmt")
c2.BIO <- getGmt("c2.cp.biocarta.v7.1.symbols.gmt")

#---split---
colnames(x = subset_cells_MP[[]])
subset_cells_MP@meta.data$merge_ident <- paste(subset_cells_MP@meta.data$seurat_clusters,subset_cells_MP@meta.data$stim,sep="_")
Idents(object = subset_cells_MP) <- 'merge_ident'
levels(x = subset_cells_MP)

Average_exprs <- AverageExpression(subset_cells_MP)
es.max.c2.BIO <- gsva(Average_exprs$RNA %>% as.matrix(), c2.BIO, mx.diff = FALSE, verbose = FALSE, parallel.sz = 1)
write.table(es.max.c2.BIO,
            file = "MP_gsva_C2.BIO_split.csv",
            sep = ",",
            row.names = T,
            quote = F)


es.max.c2.PID <- gsva(Average_exprs$RNA %>% as.matrix(), c2.PID, mx.diff = FALSE, verbose = FALSE, parallel.sz = 1)
write.table(es.max.c2.PID,
            file = "MP_gsva_c2.PID_split.csv",
            sep = ",",
            row.names = T,
            quote = F)
es.max.c2 <- rbind(es.max.c2.BIO,es.max.c2.PID)
write.table(es.max.c2,
            file = "MP_gsva_c2_split.csv",
            sep = ",",
            row.names = T,
            quote = F)



##-------------------hcluster--------------------------------####

library(dendextend)
cluster.averages <- AverageExpression(object = subset_cells_MP, return.seurat = TRUE)


DEG <- read.table("MP_culster_all_DEGs.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F)
DEGs <- DEG[which(DEG$pct.1 > 0.2 & DEG$p_val_adj < 0.05 & DEG$avg_logFC > 0.25),]

spearman_cor <- cor(cluster.averages@assays$RNA@data %>% as.matrix() %>% `[`(DEGs$gene %>% unique(), ), method = "spearman")
row.names(spearman_cor) <- mapvalues(row.names(spearman_cor), from = subset_cells_MP$seurat_clusters %>% as.character, to = subset_cells_MP$integrated_snn_res.0.3 %>% as.character)
colnames(spearman_cor) <- mapvalues(colnames(spearman_cor), from = subset_cells_MP$seurat_clusters %>% as.character, to = subset_cells_MP$integrated_snn_res.0.3 %>% as.character)

distance_of_cluster <- as.dist(1-spearman_cor)

hc <- hclust(distance_of_cluster, method = "complete")
hc <- as.dendrogram(hc)
asggden <- as.ggdend(hc)



##################-----MHC ------########


MHC1_gene_sets <- read.table("MHC_CLASS_I_geneset.txt",
                             header = F,
                             row.names = NULL,
                             stringsAsFactors = F,
                             sep = "\t") %>% dplyr::select(1) %>% `[`(,1) 

MHC1_score <- AddModuleScore(object = subset_cells_MP, features = list(MHC1_gene_sets))
MHC1_score@meta.data$merge_ident <- paste(MHC1_score@meta.data$seurat_clusters,MHC1_score@meta.data$stim,sep="_")
MHC1_score_for_vln <- MHC1_score@meta.data %>% dplyr::select(merge_ident, Cluster1)


MHC11_gene_sets <- read.table("MHC_CLASS_II_geneset.txt",
                              header = T,
                              row.names = NULL,
                              stringsAsFactors = F,
                              sep = "\t") %>% dplyr::select(1) %>% `[`(,1) 

MHC11_score <- AddModuleScore(object = subset_cells_MP, features = list(MHC11_gene_sets))
MHC11_score@meta.data$merge_ident1 <- paste(MHC11_score@meta.data$seurat_clusters,MHC11_score@meta.data$stim,sep="_")
MHC11_score_for_vln <- MHC11_score@meta.data %>% dplyr::select(merge_ident1, Cluster1)

cols <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9), pal_futurama()(12), pal_aaas()(10))[-8]

merged_data <- cbind(x = MHC1_score_for_vln, y = MHC11_score_for_vln[, 2])
meltdata <- reshape::melt(merged_data)
meltdata$variable <- meltdata$variable %>% as.character()

meltdata$variable[meltdata$variable == "x.Cluster1"] <- "MHC I"
meltdata$variable[meltdata$variable == "y"] <- "MHC II"
write.table(meltdata,"/MHC-meltdata_MP.csv", col.names = T, row.names = F, sep = ",")


