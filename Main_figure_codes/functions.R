seurat_pipeline <- function(pbmc, cells_of_each_main_cluster, mt.pct, Seurat_clustring_result_csv, cell_umap_loading_csv){
  ## ----setup, include=FALSE------------------------------------------------
  knitr::opts_chunk$set(
    tidy = TRUE,
    tidy.opts = list(width.cutoff = 120),
    message = FALSE,
    warning = FALSE
  )
  
  # Load the PBMC dataset from csv file
  specific_cells <- read.csv(cells_of_each_main_cluster, header = FALSE)
  specific_cells <- as.character(specific_cells$V1)
  pbmc <- subset(x = pbmc, cells = specific_cells)
  pbmc <- subset(pbmc, subset = percent.mt < mt.pct)

  ## ----var_features, fig.height=5, fig.width=11----------------------------
  pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'vst', nfeatures = 1000)
  
  ## ---- Z-score transformation -----------------
  all.genes <- rownames(x = pbmc)
  pbmc <- ScaleData(pbmc, verbose = FALSE)
  pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
  pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
  pbmc <- FindClusters(pbmc, verbose = FALSE)
  
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)  

  ## write the clustering result to files
  write.csv(pbmc$seurat_clusters, file = Seurat_clustring_result_csv, quote = FALSE, sep = ',')
  write.csv(pbmc@reductions$umap@cell.embeddings, file = cell_umap_loading_csv, quote = FALSE)
}


seurat_pipeline_noIG <- function(pbmc, cells_of_each_main_cluster, mt.pct, Seurat_clustring_result_csv, cell_umap_loading_csv){
  # This version excluds genes  IGHV, IGLV, IGKV special for clustering cells in Plasma lineage. 
  # and also regress out these genes using their mean values.
  ## ----setup, include=FALSE------------------------------------------------
  knitr::opts_chunk$set(
    tidy = TRUE,
    tidy.opts = list(width.cutoff = 120),
    message = FALSE,
    warning = FALSE
  )
  
  # Load the PBMC dataset from csv file
  # pbmc.data <- Read10X(data.dir = sc_Seq_bc_matrix)
  specific_cells <- read.csv(cells_of_each_main_cluster, header = FALSE)
  specific_cells <- as.character(specific_cells$V1)
  pbmc.data <- pbmc[["integrated"]]@data
  ind = match(specific_cells, colnames(pbmc.data))
  pbmc.data <- pbmc.data[, ind[!is.na(ind)]]
  ind = grepl('IGHV', row.names(pbmc.data))|grepl('IGLV', row.names(pbmc.data))|grepl('IGKV', row.names(pbmc.data))
  mat = as.matrix(pbmc.data[ind, ])
  variation = colMeans(mat)
  variation = data.frame(variation)
  names(variation) = c("HG")
  pbmc.data.noHG = data.frame(pbmc.data[!ind, ])
  pbmc <- CreateSeuratObject(counts = pbmc.data.noHG, project = "pbmc3k", min.cells = 0, min.features = 0)
  pbmc[["HG"]] <- variation
  pbmc <- FindVariableFeatures(object = pbmc,selection.method = 'vst', nfeatures = 2000)  
  #pbmc <- ScaleData(pbmc, vars.to.regress = "HG")
  pbmc <- ScaleData(pbmc) # no regress out
  pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
  pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)
  pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
  pbmc <- FindClusters(pbmc, verbose = FALSE)
  ## write the clustering result to files
  write.csv(pbmc$seurat_clusters, file = Seurat_clustring_result_csv, quote = FALSE, sep = ',')
  write.csv(pbmc@reductions$umap@cell.embeddings, file = cell_umap_loading_csv, quote = FALSE)
}




