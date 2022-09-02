library(dplyr)
library(Seurat)
library(ggplot2)

setwd('F:/projects/SuPERRSeq-part3/Seurat_R_version')

input_dir = "../Main_gated_result/BM/singlets_each_main_gate"
result_dir = "GEX_subclusters_BM_AnchorVersion/singlets_version"

load('BM_integrated_seurat_object.RData')
sc_data_main_clusters <- list.files(input_dir, recursive = FALSE, full.names = FALSE)
sc_data_main_clusters_1 = sc_data_main_clusters[] # leave out B cells
mitochondrial_percentage = rep(20, length(sc_data_main_clusters)) # CD4 14, HSPC 19, 

source('functions.R')
for (cluster in sc_data_main_clusters_1){
  cells_of_each_main_cluster <- paste(input_dir, '/', cluster, sep ="")
  Seurat_clustring_result_csv <- paste(result_dir, "/Seurat_clustering_result_", cluster, sep = "")
  cell_umap_loading_csv <- paste(result_dir, "/Seurat_cell_umap_loading_", cluster, sep = "")
  seurat_pipeline(pancreas.integrated, cells_of_each_main_cluster, mitochondrial_percentage[], Seurat_clustring_result_csv, cell_umap_loading_csv)

}

sc_data_main_clusters_2 = sc_data_main_clusters[] # B cells, no IG
for (cluster in sc_data_main_clusters_2){
  cells_of_each_main_cluster <- paste(input_dir, '/', cluster, sep ="")
  Seurat_clustring_result_csv <- paste(result_dir, "/Seurat_clustering_result_noIG_", cluster, sep = "")
  cell_umap_loading_csv <- paste(result_dir, "/Seurat_cell_umap_loading_noIG", cluster, sep = "")
  seurat_pipeline_noIG(pancreas.integrated, cells_of_each_main_cluster, mitochondrial_percentage[], Seurat_clustring_result_csv, cell_umap_loading_csv)
}
