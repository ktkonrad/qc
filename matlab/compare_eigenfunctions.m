addpath('../vergini');

s = load_sta('../c/qust_low');
f_low = reshape(s(1,:,:), size(s,2), size(s,3))';

s = load_sta('../c/qust_high');
f_high = reshape(s(1,:,:), size(s,2), size(s,3))';

ratio = 20;
zoom = [160 210 5 45];
rect = [183 20 5 5];

figure; imagesc(f_low>0);
axis(zoom); axis off; axis square;
rectangle('Position', rect, 'EdgeColor', 'y', 'LineWidth', 3);
print('-deps2c', '../documents/thesis/figs/interpolation/eigenfunction_error_low.eps');

figure; imagesc(f_high>0);
axis((zoom-[1 1 1 1]).*ratio); axis off; axis square;
rectangle('Position', (rect-[1 1 0 0]).*ratio, 'EdgeColor', 'y', 'LineWidth', 3);
print('-deps2c', '../documents/thesis/figs/interpolation/eigenfunction_error_high.eps');
