library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
setwd('E:\\Dropbox (GaTech)\\projects\\SuPERRSeq-part3\\Seurat_R_version')

# Load the BM dataset
bm.data.1 <- Read10X(data.dir = "..\\raw_data_version20200227\\GEX+ADT\\BM2")$`Gene Expression`
bm.data.2 <- Read10X(data.dir = "..\\raw_data_version20200227\\GEX+ADT\\BM3")$`Gene Expression`
colnames(bm.data.1) = paste(colnames(bm.data.1), '_BM2', sep ='')
colnames(bm.data.2) = paste(colnames(bm.data.2), '_BM3', sep ='')

# sample integration
bm1 = CreateSeuratObject(counts = bm.data.1, min.cells = 3, min.features = 200)
bm2 = CreateSeuratObject(counts = bm.data.2, min.cells = 3, min.features = 200)
bm1 <- PercentageFeatureSet(bm1, pattern = "^MT-", col.name = "percent.mt")
bm2 <- PercentageFeatureSet(bm2, pattern = "^MT-", col.name = "percent.mt")
sample.list = list(bm1, bm2)
names(sample.list) = c("BM2", "BM3")
for (i in 1:length(sample.list)) {
  sample.list[[i]] <- NormalizeData(sample.list[[i]], verbose = FALSE)
  sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst",
                                             nfeatures = 2000, verbose = FALSE)
}
reference.list <- sample.list
sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
sample.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:30)
DefaultAssay(sample.integrated) <- "integrated"
save(sample.integrated, file = 'BM_integrated_seurat_object.RData')

# Run the standard workflow, genrating seurat clusters for visualization
sample.integrated <- ScaleData(sample.integrated, verbose = FALSE)
sample.integrated <- RunPCA(sample.integrated, npcs = 30, verbose = FALSE)
sample.integrated <- RunUMAP(sample.integrated, reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, dims = 1:30, verbose = FALSE)
sample.integrated <- FindClusters(sample.integrated, verbose = FALSE)
write.csv(sample.integrated@reductions$umap@cell.embeddings, file = 'umap_cell_embeddings_Anchor_BM.csv', quote = FALSE)
write.csv(sample.integrated$seurat_clusters, file = 'Seurat_clustring_result_Anchor_BM.csv', quote = FALSE)
