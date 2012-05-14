g = dlmread('../data/grid_0.dat');
four_cells = [1 1 -1 -1 ; 1 1 -1 -1 ; -1 -1 1 1 ; -1 -1 1 1];
n = size(g,1)/4;
f = repmat(four_cells, n, n);
fontsize = 20;

figure; imagesc(f);
axis off;
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/percolation/checkerboard_implementation.eps');



figure; imagesc(g);
axis off;
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/percolation/perturbed_implementation.eps');
