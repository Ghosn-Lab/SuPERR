clear all; clear all; clc
addpath(genpath('tools'))

%%
load raw_data_DSB_ADT_BM.mat
raw_data = full(raw_data_ADT); %
cell_names = raw_cell_names_ADT;
anti_seq_gene_names = raw_gene_names_ADT;
outlier_percentage = 5;
raw_data_umis = log10(raw_data(end, :)+1);
raw_data(end, :) = raw_data_umis;
 
%% Plasma vs Non_Plasma cells gating
figure(1);
before_gating_contour_plot(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans',  outlier_percentage);
set(gca, 'xlim', [-20, 70]);

h01 = drawpolygon; h02 = drawpolygon; h03 = drawpolygon; 
figure0_gate1_position = h01.Position; figure0_gate2_position = h02.Position; figure0_gate3_position = h03.Position;

selected_indices_figure0_gate1 = manually_gating_v4(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans', figure0_gate1_position);
cells_CD138low_plasma = cell_names(selected_indices_figure0_gate1);

selected_indices_figure0_gate2 = manually_gating_v4(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans', figure0_gate2_position);
cells_CD138hi_plasma = cell_names(selected_indices_figure0_gate2);

selected_indices_figure0_gate3 = manually_gating_v4(raw_data, anti_seq_gene_names, 'CD138', 'Total_Ig_trans', figure0_gate3_position);
cells_non_plasma = cell_names(selected_indices_figure0_gate3);
data_non_plasma = raw_data(:, selected_indices_figure0_gate3);

%% B cells gating 
figure(2)
before_gating_contour_plot(data_non_plasma, anti_seq_gene_names, 'CD20', 'CD19', outlier_percentage);
set(gca, 'xlim', [-10, 60], 'ylim', [-3, 25]);

h11 = drawpolygon; h12 = drawpolygon;
figure1_gate1_position = h11.Position; figure1_gate2_position = h12.Position;

selected_indices_figure1_gate1 = manually_gating_v4(data_non_plasma, anti_seq_gene_names, 'CD20', 'CD19', figure1_gate1_position);
cells_B = cells_non_plasma(selected_indices_figure1_gate1);
data_B = data_non_plasma(:, selected_indices_figure1_gate1);

selected_indices_figure1_gate2 = manually_gating_v4(data_non_plasma, anti_seq_gene_names, 'CD20', 'CD19', figure1_gate2_position);
cells_non_B = cells_non_plasma(selected_indices_figure1_gate2);
data_non_B = data_non_plasma(:, selected_indices_figure1_gate2);

%% Projenitors
figure(3);
% before_gating_contour_plot(data_B, anti_seq_gene_names, 'CD10', 'CD34', outlier_percentage);
before_gating_contour_plot(data_B, anti_seq_gene_names, 'CD10', 'CD21', outlier_percentage);
set(gca, 'xlim', [-5, 45], 'ylim', [-5, 30]);

h21 = drawpolygon; h22 = drawpolygon; 
figure2_gate1_position = h21.Position; figure2_gate2_position = h22.Position; 

selected_indices_figure2_gate1 = manually_gating_v4(data_B, anti_seq_gene_names, 'CD10', 'CD21', figure2_gate1_position);
selected_indices_figure2_gate2 = manually_gating_v4(data_B, anti_seq_gene_names, 'CD10', 'CD21', figure2_gate2_position);

cells_naive_B = cells_B(selected_indices_figure2_gate1);
cells_CD10hi_projenitors = cells_B(selected_indices_figure2_gate2);

%% Non B subtypes
figure(4)
before_gating_contour_plot(data_non_B, anti_seq_gene_names, 'CD11c', 'CD34', outlier_percentage);
set(gca, 'xlim', [-10, 40], 'ylim', [-3,30]);

h31 = drawpolygon; h32 = drawpolygon;
figure3_gate1_position = h31.Position; figure3_gate2_position = h32.Position;
selected_indices_figure3_gate1 = manually_gating_v4(data_non_B, anti_seq_gene_names, 'CD11c', 'CD34', figure3_gate1_position);
selected_indices_figure3_gate2 = manually_gating_v4(data_non_B, anti_seq_gene_names, 'CD11c', 'CD34', figure3_gate2_position);

cells_HSPCs = cells_non_B(selected_indices_figure3_gate1);
cells_Myeloid = cells_non_B(selected_indices_figure3_gate2);

data_HSPCs = data_non_B(:, selected_indices_figure3_gate1);
data_Myeloid = data_non_B(:, selected_indices_figure3_gate2);

%%
save step02_anti_seq_manual_gating_BM_result.mat