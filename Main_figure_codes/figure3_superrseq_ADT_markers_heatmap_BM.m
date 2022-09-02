clc; clear all; close all
addpath(genpath('tools'));
%%
load raw_data_DSB_ADT_BM.mat
ab_data = raw_data_ADT; 
ab_cell_names = raw_cell_names_ADT;
anti_seq_gene_names = raw_gene_names_ADT;

ab_data = ab_data(~contains(anti_seq_gene_names, 'CD34'), :);
anti_seq_gene_names = anti_seq_gene_names(~contains(anti_seq_gene_names, 'CD34'));
%%
cluster_num = 0;
superrseq_clusters_adt_mean_matrix = [];

main_clusters = {'total_CD138low_plasma_cells', 'total_CD138high_plasma_cells', 'total_B_cells',  'HSPCs', 'Myeloid'};
for i = 1:length(main_clusters)
    table_content = read_xls_v2(['Seurat_R_version\GEX_subclusters_BM_AnchorVersion\singlets_version\Seurat_clustering_result_' main_clusters{i} '.csv'],',',1);
    tmp_cell_names = table_content(:, 1);
    tmp_cell_cluster = str2double(table_content(:, 2));
    tmp_cell_cluster_unique = unique(tmp_cell_cluster);
    for j = 1:length(tmp_cell_cluster_unique)
        cells_in_this_sub_cluster = tmp_cell_names(tmp_cell_cluster == j-1);
        column_indicator = ismember(ab_cell_names, cells_in_this_sub_cluster);
        data_tmp = ab_data(:, column_indicator);
        superrseq_clusters_adt_mean_matrix = [superrseq_clusters_adt_mean_matrix mean(data_tmp, 2)];
     end
    
    cluster_num = cluster_num + length(tmp_cell_cluster_unique);
end

%% only keep those markers not used for major lineage identification
ind_gene_left = ~ismember(anti_seq_gene_names, {'Total_Ig_trans','CD138','CD19','CD20','CD34','CD11c'});
superrseq_clusters_adt_mean_matrix = superrseq_clusters_adt_mean_matrix(ind_gene_left,:);
anti_seq_gene_names = anti_seq_gene_names(ind_gene_left);

%% expression level
load superrseq_cluster_order_BM.mat I

ab_subcluster_expression_matrix_normalization = normalize(superrseq_clusters_adt_mean_matrix, 2);

%%
% figure;
% imagesc(ab_subcluster_expression_matrix_normalization(:, I));
% set(gca, 'XTick', 1:cluster_num, 'XTickLabel', I, 'YTick', 1:cluster_num, 'YTickLabel', anti_seq_gene_names);
% dim_background = get(gca, 'position');
% % patching cell cluster labels
% start_point_x = dim_background(1);
% start_point_y = dim_background(2) + dim_background(4);
% total_length = dim_background(3);
% 
% for i = 1:length(main_clusters)
%     table_content = read_xls_v2(['Seurat_R_version\GEX_suclusters_PBMC\Seurat_clustering_result_sc_' main_clusters{i} '.csv'],',',1);
%     tmp_cell_cluster = str2double(table_content(:, 2));
%     tmp_cell_cluster_unique = unique(tmp_cell_cluster);
%     current_width = length(tmp_cell_cluster_unique)/cluster_num*total_length;
%     annotation('textbox',[start_point_x start_point_y current_width 0.04],'String', main_clusters{i}, ...
%         'BackgroundColor', [i/6*255 202 215]/255, 'FaceAlpha',1, 'Interpreter', 'none', 'VerticalAlignment', 'middle',...
%         'HorizontalAlignment', 'center');
%     start_point_x = start_point_x + current_width;
% end

%%
A = ab_subcluster_expression_matrix_normalization(:, I);
cg = clustergram(A);
set(cg, 'Cluster', 1, 'RowLabels', anti_seq_gene_names, 'ColumnLabels', I, 'Colormap',parula(256), 'ShowDendrogram', 'off');

%%
figure;
B = zeros(size(ab_subcluster_expression_matrix_normalization));
for i = 1:length(cg.RowLabels)
    B(i, :) = ab_subcluster_expression_matrix_normalization(ismember(anti_seq_gene_names, cg.RowLabels{i}), I);
end
%C = B; C(C>3) = 3;
imagesc(B); colorbar; 
caxis([-3, 4]);
set(gca, 'XTick', 1:cluster_num, 'XTickLabel', I, 'YTick', 1:length(cg.RowLabels), 'YTickLabel', strrep(cg.RowLabels, '_', ' '));
set(gcf, 'Position', 1.0e+03*[0.0026    0.2154    1.5320    0.6066])
% patching cell cluster labels
dim_background = get(gca, 'position');
start_point_x = dim_background(1);
start_point_y = dim_background(2) + dim_background(4);
total_length = dim_background(3);
panel_color = distinguishable_colors(10, 'k');

main_clusters_name = {'B cells', 'HSPCs', 'Myeloid'};

current_width = 4/cluster_num*total_length;
annotation('textbox',[start_point_x start_point_y current_width 0.04], 'String', 'Plasma cells', ...
    'BackgroundColor', panel_color(1, :), 'FaceAlpha',0.3, 'Interpreter', 'none', 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'w');
start_point_x = start_point_x + current_width;
for i = 3:length(main_clusters)
    table_content = read_xls_v2(['Seurat_R_version\GEX_subclusters_BM_AnchorVersion\singlets_version\Seurat_clustering_result_' main_clusters{i} '.csv'],',',1);
    tmp_cell_cluster = str2double(table_content(:, 2));
    tmp_cell_cluster_unique = unique(tmp_cell_cluster);
    current_width = length(tmp_cell_cluster_unique)/cluster_num*total_length;
    annotation('textbox',[start_point_x start_point_y current_width 0.04], 'String', main_clusters_name{i-2}, ...
        'BackgroundColor', panel_color(i-1, :), 'FaceAlpha',0.3, 'Interpreter', 'none', 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'w');
    start_point_x = start_point_x + current_width;
end

savepdf('superrseq_ADT_heatmap_BM.pdf')