function selected_indices = manually_gating_v3(data_last_figure, anti_gene_names, gene_x, gene_y, gate_coordinates)
% using data gated out from late step, do gating, all gates are polygons

data_this_step = [data_last_figure(ismember(anti_gene_names, gene_x), :); data_last_figure(ismember(anti_gene_names, gene_y), :)]';

%%
%this gate is a polygon
selected_indices = inpolygon(data_this_step(:, 1), data_this_step(:, 2), ...
    gate_coordinates(:, 1), gate_coordinates(:, 2));

end