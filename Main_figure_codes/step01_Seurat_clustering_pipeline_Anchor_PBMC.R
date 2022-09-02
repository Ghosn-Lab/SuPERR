library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

setwd('E:\\Dropbox (GaTech)\\projects\\SuPERRSeq-part3\\Seurat_R_version')
# Load the BM dataset
pbmc.data.1 <- Read10X(data.dir = "..\\raw_data_version20200227\\GEX+ADT\\PBMC1")$`Gene Expression`
pbmc.data.2 <- Read10X(data.dir = "..\\raw_data_version20200227\\GEX+ADT\\PBMC2")$`Gene Expression`
pbmc.data.3 <- Read10X(data.dir = "..\\raw_data_version20200227\\GEX+ADT\\PBMC3")$`Gene Expression`

colnames(pbmc.data.1) = paste(colnames(pbmc.data.1), '_PBMC1', sep ='')
colnames(pbmc.data.2) = paste(colnames(pbmc.data.2), '_PBMC2', sep ='')
colnames(pbmc.data.3) = paste(colnames(pbmc.data.3), '_PBMC3', sep ='')

pbmc1 = CreateSeuratObject(counts = pbmc.data.1, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc2 = CreateSeuratObject(counts = pbmc.data.2, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc3 = CreateSeuratObject(counts = pbmc.data.3, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc1 <- PercentageFeatureSet(pbmc1, pattern = "^MT-", col.name = "percent.mt")
pbmc2 <- PercentageFeatureSet(pbmc2, pattern = "^MT-", col.name = "percent.mt")
pbmc3 <- PercentageFeatureSet(pbmc3, pattern = "^MT-", col.name = "percent.mt")

sample.list = list(pbmc1, pbmc2, pbmc3)
names(sample.list) = c("PBMC1", "PBMC2", "PBMC3")
for (i in 1:length(sample.list)) {
  sample.list[[i]] <- NormalizeData(sample.list[[i]], verbose = FALSE)
  sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

reference.list <- sample.list[c("PBMC1", "PBMC2", "PBMC3")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

DefaultAssay(pancreas.integrated) <- "integrated"

save(pancreas.integrated@assays, file = 'PBMC_integrated_seurat_object.RData')

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)

pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30, verbose = FALSE)
pancreas.integrated <- FindClusters(pancreas.integrated, verbose = FALSE)

write.csv(pancreas.integrated@reductions$umap@cell.embeddings, file = 'umap_cell_embeddings_Anchor_PBMC.csv', quote = FALSE)
write.csv(pancreas.integrated$seurat_clusters, file = 'Seurat_clustring_result_Anchor_PBMC_test.csv', quote = FALSE)
