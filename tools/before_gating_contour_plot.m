function selected_indices = before_gating_contour_plot(data_last_figure, anti_gene_names, gene_x, gene_y, outlier_percentage)
% using data gated out from last step, draw the contour plot

data_this_step = [data_last_figure(ismember(anti_gene_names, gene_x), :); data_last_figure(ismember(anti_gene_names, gene_y), :)]';

%figure;

FlowJo_contour2D(data_this_step(:, 1), data_this_step(:, 2), outlier_percentage);

xlabel(gene_x); ylabel(gene_y, 'Interpreter', 'none');

end