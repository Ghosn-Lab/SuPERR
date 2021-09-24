clear all; close all; clc 
addpath(genpath('tools'))

%%
load raw_data_DSB_ADT_PBMC.mat
raw_data = full(raw_data_ADT); %
cell_names = raw_cell_names_ADT;
anti_seq_gene_names = raw_gene_names_ADT;
outlier_percentage = 10; % try 10
raw_data_umis = log10(raw_data(end, :)+1);
raw_data(end, :) = raw_data_umis;

%% Plasma vs Non_Plasma cells gating
figure(1);
before_gating_contour_plot(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans',  outlier_percentage);
set(gca, 'ylim', [-2, inf], 'xlim', [-9, 43])

h01 = drawpolygon; h02 = drawpolygon; 
figure0_gate1_position = h01.Position; figure0_gate2_position = h02.Position;

selected_indices_figure0_gate1 = manually_gating_v4(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans', figure0_gate1_position);
cells_plasma = cell_names(selected_indices_figure0_gate1);
data_plasma = raw_data(:, selected_indices_figure0_gate1);

selected_indices_figure0_gate2 = manually_gating_v4(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans', figure0_gate2_position);
cells_non_plasma = cell_names(selected_indices_figure0_gate2);
data_non_plasma = raw_data(:, selected_indices_figure0_gate2);

%% B cells gating 
figure(2)
before_gating_contour_plot(data_non_plasma, anti_seq_gene_names, 'CD20', 'CD19', outlier_percentage);
set(gca, 'xlim', [-20, 75], 'ylim', [-10, 30]);

h11 = drawpolygon; h12 = drawpolygon;
figure1_gate1_position = h11.Position; figure1_gate2_position = h12.Position;

selected_indices_figure1_gate1 = manually_gating_v4(data_non_plasma, anti_seq_gene_names, 'CD20', 'CD19', figure1_gate1_position);
cells_total_B = cells_non_plasma(selected_indices_figure1_gate1);
data_selected_figure1_B_cells = data_non_plasma(:, selected_indices_figure1_gate1);

selected_indices_figure1_gate2 = manually_gating_v4(data_non_plasma, anti_seq_gene_names, 'CD20', 'CD19', figure1_gate2_position);
cells_non_B = cells_non_plasma(selected_indices_figure1_gate2);
data_non_B = data_non_plasma(:, selected_indices_figure1_gate2);

%% NK cells
figure(3);
before_gating_contour_plot(data_non_B, anti_seq_gene_names, 'CD3', 'CD56', outlier_percentage);
set(gca, 'xlim', [-3, 20]);set(gca, 'ylim', [-3, 22]);

h21 = drawpolygon; h22 = drawpolygon;  h23 = drawpolygon;  h24 = drawpolygon;  
figure2_gate1_position = h21.Position; figure2_gate2_position = h22.Position; figure2_gate3_position = h23.Position; figure2_gate4_position = h24.Position; 

selected_indices_figure2_gate1 = manually_gating_v4(data_non_B, anti_seq_gene_names, 'CD3', 'CD56', figure2_gate1_position);
cells_toatl_NK = cells_non_B(selected_indices_figure2_gate1);
data_total_NK = data_non_B(:, selected_indices_figure2_gate1);

selected_indices_figure2_gate2 = manually_gating_v4(data_non_B, anti_seq_gene_names, 'CD3', 'CD56', figure2_gate2_position);
cells_non_NK = cells_non_B(selected_indices_figure2_gate2);
data_non_NK = data_non_B(:, selected_indices_figure2_gate2);

selected_indices_figure2_gate3 = manually_gating_v4(data_non_B, anti_seq_gene_names, 'CD3', 'CD56', figure2_gate3_position);
cells_NKT = cells_non_B(selected_indices_figure2_gate3);
data_NKT = data_non_B(:, selected_indices_figure2_gate3);

selected_indices_figure2_gate4 = manually_gating_v4(data_non_B, anti_seq_gene_names, 'CD3', 'CD56', figure2_gate4_position);
cells_T = cells_non_B(selected_indices_figure2_gate4);
data_T = data_non_B(:, selected_indices_figure2_gate4);

%% CD56hi
figure(4);
before_gating_contour_plot(data_total_NK, anti_seq_gene_names, 'CD16', 'CD56', outlier_percentage);
set(gca, 'xlim', [-3, 25]);set(gca, 'ylim', [-3, 22]);
h81 = drawpolygon; h82 = drawpolygon;  
figure8_gate1_position = h81.Position; figure8_gate2_position = h82.Position; 

selected_indices_figure8_gate1 = manually_gating_v4(data_total_NK, anti_seq_gene_names, 'CD16', 'CD56', figure8_gate1_position);
cells_CD56hi_NK = cells_toatl_NK(selected_indices_figure8_gate1);
data_CD56hi_NK = data_total_NK(:, selected_indices_figure8_gate1);

selected_indices_figure8_gate2 = manually_gating_v4(data_total_NK, anti_seq_gene_names, 'CD16', 'CD56', figure8_gate2_position);
cells_NK = cells_toatl_NK(selected_indices_figure8_gate2);
data_NK = data_total_NK(:, selected_indices_figure8_gate2);

%% Different monocytes
figure(5)
before_gating_contour_plot(data_non_NK, anti_seq_gene_names, 'CD14', 'CD16', outlier_percentage);
set(gca, 'xlim', [-3, 15], 'ylim', [-3, 25]);

h31 = drawpolygon; h32 = drawpolygon; h33 = drawpolygon; 
figure3_gate1_position = h31.Position;  figure3_gate2_position = h32.Position; 
figure3_gate3_position = h33.Position; 

selected_indices_figure3_gate1 = manually_gating_v4(data_non_NK, anti_seq_gene_names, 'CD14', 'CD16', figure3_gate1_position);
cells_non_classical_monocytes = cells_non_NK(selected_indices_figure3_gate1);

selected_indices_figure3_gate2 = manually_gating_v4(data_non_NK, anti_seq_gene_names, 'CD14', 'CD16', figure3_gate2_position);
cells_intermediate_monocytes = cells_non_NK(selected_indices_figure3_gate2);

selected_indices_figure3_gate3 = manually_gating_v4(data_non_NK, anti_seq_gene_names, 'CD14', 'CD16', figure3_gate3_position);
cells_classical_monocyctes = cells_non_NK(selected_indices_figure3_gate3);

%% CD4+ vs CD8+ T cells
figure(6)
before_gating_contour_plot(data_T, anti_seq_gene_names, 'CD8a', 'CD4', outlier_percentage);
set(gca, 'xlim', [-4, 35], 'ylim', [-3, 25]);

h41 = drawpolygon; h42 = drawpolygon; 
figure4_gate1_position = h41.Position; figure4_gate2_position = h42.Position; 

selected_indices_figure4_gate1 = manually_gating_v4(data_T, anti_seq_gene_names, 'CD8a', 'CD4', figure4_gate1_position);
cells_CD4hi_T = cells_T(selected_indices_figure4_gate1);
data_CD4hi_T = data_T(:, selected_indices_figure4_gate1);

selected_indices_figure4_gate2 = manually_gating_v4(data_T, anti_seq_gene_names, 'CD8a', 'CD4', figure4_gate2_position);
cells_CD8hi_T = cells_T(selected_indices_figure4_gate2);

%% Treg cells
figure(7)
before_gating_contour_plot(data_CD4hi_T, anti_seq_gene_names, 'CD127', 'CD25', outlier_percentage);
set(gca, 'xlim', [-5, 23], 'ylim', [-1, 10]);

h51 = drawpolygon; 
figure5_gate1_position = h51.Position; 

selected_indices_figure5_gate1 = manually_gating_v4(data_CD4hi_T, anti_seq_gene_names, 'CD127', 'CD25', figure5_gate1_position);
cells_treg = cells_CD4hi_T(selected_indices_figure5_gate1);

%%
save step02_anti_seq_manual_gating_PBMC_result.mat