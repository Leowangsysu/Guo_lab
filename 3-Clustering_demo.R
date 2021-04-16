library(Seurat)  
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
library(DoubletFinder)
library(parallel)
library(stringr)
library(ggsci)

paths <- "/2-Mapping/"
tmp <- c()

for (timepoint in c("sample_name1_cDNA", "sample_name2_cDNA", "sample_name3_cDNA")) {
  subset_cell <- Read10X(data.dir = paste0(paths, timepoint, "/outs/filtered_feature_bc_matrix/"))
  colnames(subset_cell) <- paste0(colnames(subset_cell), "_", strsplit(timepoint, split = "_") %>% unlist %>% `[`(2))
  tmp <- append(tmp, subset_cell)
}

tmp <- do.call(cbind, tmp)
subset_cell <- CreateSeuratObject(counts = tmp, project = "Liver", min.cells = 3, min.features = 200)

##rm doublet 
a <- read.table("rm_Doubletlist.txt")
a <- as.character(a$V1)
b <-  colnames(subset_cell)
b <- b[!(b%in%a)]
subset_cell <- subset(subset_cell,cells=b)

subset_cell[["percent.mt"]] <- PercentageFeatureSet(subset_cell, pattern = "^MT-")
subset_cell <- subset(subset_cell, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & percent.mt < 25) 


#integration#

subset_cell@meta.data$stim <- str_sub(rownames(subset_cell@meta.data),-2)
subset_cell <- SplitObject(subset_cell, split.by = "stim")
subset_cell <- lapply(X = subset_cell, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
subset_cell <- FindIntegrationAnchors(object.list = subset_cell, dims = 1:20)
subset_cell.combined <- IntegrateData(anchorset = subset_cell, dims = 1:20)
subset_cell <- subset_cell.combined 



### find the DEGs
DefaultAssay(subset_cells) <- "RNA"
all.genes <- rownames(subset_cells)
subset_cells <- ScaleData(subset_cells, features = all.genes)

result <- mclapply(as.numeric(levels(subset_cells@active.ident)),
                   FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                   mc.cores = 36)
RESULT <- result

roundN <- 1
while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
  if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
    print(recalculate_clusters)
    result1 <- mclapply(recalculate_clusters,
                        FUN =  function(x) {FindMarkers(subset_cells, ident.1 = x, ident.2 = NULL)},
                        mc.cores = 35)
  }
  print(roundN + 1)
  for(i in 1:length(recalculate_clusters)){
    result[c(recalculate_clusters+1)[i]] <- result1[i]
  }
}

all_markers <- do.call(rbind, result)
all_markers$gene <- unlist(mapply(rownames, result))
all_markers$cluster <- rep(levels(subset_cells@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
subset_cells.markers <- all_markers

  TOP_N <- function(x, n, pct.1 = 0.1, first = "avg_logFC", second = "p_val_adj", sig.padj = NULL, fc.threshold = c(-0.25, 0.25)){
    if(table(names(x) %in% c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))["TRUE"] == 7){
      reordered_daf <- data.frame()
      if (!is.null(pct.1) & (pct.1 > 0 & pct.1 < 1)) {
        x <- x[x$pct.1 >= pct.1, ]
      }
      if (!is.null(sig.padj) & is.numeric(sig.padj)) {
        x <- x[x$p_val_adj <= sig.padj, ]
      }
      if (length(fc.threshold) == 2) { 
        fc.threshold <- sort(fc.threshold)
        x <- x[(x$avg_logFC < fc.threshold[1] | x$avg_logFC > fc.threshold[2]), ]
      } else if (length(fc.threshold) == 1) {
        x <- x[x$avg_logFC >= fc.threshold, ]
      }
      
      for (i in unique(x$cluster)) {
        if(n > dim(x[x$cluster == i, ])[1]){
          message("FBI warning, n < ", n)
          tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first), get(second)))[1:(dim(x[x$cluster == i, ])[1]), ]
          reordered_daf <- rbind(reordered_daf, tmp)
        } else {
          tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first), get(second)))[1:n, ]
          reordered_daf <- rbind(reordered_daf, tmp)
        }
      }
    }
    else {
      message("be careful !")
    }
    return(reordered_daf)
  } 

subset_cells.markers %>% TOP_N(50) -> top50
subset_cells.markers <- subset_cells.markers %>% TOP_N(5000)


save(subset_cells, file = "/liver_cluster.Rdata")
save(subset_cells.markers, file = "/liver_cluster_DEG.Rdata")
