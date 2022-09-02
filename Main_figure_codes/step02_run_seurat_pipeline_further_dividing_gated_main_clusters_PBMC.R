library(dplyr)
library(Seurat)
library(ggplot2)

setwd('E:/Dropbox (GaTech)/projects/SuPERRSeq-part3/Seurat_R_version')

source('functions.R')
input_dir = "../Main_gated_result/PBMC/singlets_each_main_gate"
result_dir = "GEX_subclusters_PBMC_AnchorVersion/singlets_version"
load('PBMC_integrated_seurat_object_v1.RData')
sc_data_main_clusters <- list.files(input_dir, recursive = FALSE, full.names = FALSE) 

mt.pct = rep(20, length(sc_data_main_clusters)) # mt<15 for CD4T 
for (cluster in sc_data_main_clusters){
  cells_of_each_main_cluster <- paste(input_dir, '/', cluster, sep ="")
  Seurat_clustring_result_csv <- paste(result_dir, "/Seurat_clustering_result_", cluster, sep = "")
  cell_umap_loading_csv <- paste(result_dir, "/Seurat_cell_umap_loading_", cluster, sep = "")
  Seurat_clustring_markers_csv <- paste(result_dir, "/Seurat_cell_clustering_markers_", cluster, sep = "")
  seurat_pipeline(pancreas.integrated, cells_of_each_main_cluster, mt.pct[cluster], Seurat_clustring_result_csv, cell_umap_loading_csv,Seurat_clustring_markers_csv)
  }
