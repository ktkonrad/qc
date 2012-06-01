addpath('../vergini');

s = load_sta('../c/qust_low');
f_low = reshape(s(1,:,:), size(s,2), size(s,3))';

s = load_sta('../c/qust_high');
f_high = reshape(s(1,:,:), size(s,2), size(s,3))';

ratio = 20;
zoom = [210 250 80 120];
rect = [229 98 5 5];

figure; imagesc(f_low>0);
axis off; axis equal; axis(zoom);
rectangle('Position', rect, 'EdgeColor', 'y', 'LineWidth', 3);
print('-deps2c', '../documents/thesis/figs/interpolation/eigenfunction_error_low.eps');

figure; imagesc(f_high>0);
axis off; axis equal; axis((zoom-[2 2 2 2]).*ratio);
rectangle('Position', (rect-[2.5 1.5 0 0]).*ratio, 'EdgeColor', 'y', 'LineWidth', 3);
print('-deps2c', '../documents/thesis/figs/interpolation/eigenfunction_error_high.eps');
