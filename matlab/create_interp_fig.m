addpath('../vergini/')

x = 52; y = 51;
rho = 30; % upsample ratio

s = load_sta('../c/qust');
f = reshape(s(1,:,:),size(s,2),size(s,3))';

figure; imagesc(f>0); axis off;
axis([x-4 x+6 y-4 y+6]);
print('-deps2c', '../documents/thesis/figs/interpolation/eigenfunction_zoom_domains.eps');
figure; imagesc(f); caxis([-2.5,2.5]); axis off;
axis([x-4 x+6 y-4 y+6]);
print('-deps2c', '../documents/thesis/figs/interpolation/eigenfunction_zoom.eps');

%%
in = read_stencil(sprintf('../c/interp_input_%d_%d.dat', x, y))';
figure; imagesc(in); axis off;
print('-deps2c', '../documents/thesis/figs/interpolation/interp_input.eps');
in_domains = 2*+(in > 0)-1;
in_domains([1 2 5 6], [1 6]) = -1;
in_domains([1 6], [1 2 5 6]) = -1;
figure; imagesc(in_domains); axis off;
print('-deps2c', '../documents/thesis/figs/interpolation/interp_input_domains.eps');

%%
out = read_dumped(sprintf('../c/interp_output_%d_%d.dat', x, y), rho+1, rho+1)';
figure; imagesc(out); axis off;
print('-deps2c', '../documents/thesis/figs/interpolation/interp_output.eps');
figure; imagesc(out>0); axis off;
print('-deps2c', '../documents/thesis/figs/interpolation/interp_output_domains.eps');

%%
N = 3;
large_out = zeros(N*rho+1);
for i = 1:N
    for j = 1:N
        A = read_dumped(sprintf('../c/interp_output_%d_%d.dat', ...
            x-ceil(N/2)+j, y-ceil(N/2)+i), rho+1, rho+1)';
        large_out((i-1)*rho+1:i*rho+1, (j-1)*rho+1:j*rho+1) = A;
    end
end
figure; imagesc(large_out); axis off;
line([30.5, 30.5], [1, N*rho+1], 'linewidth', 4);
line([60.5, 60.5], [1, N*rho+1], 'linewidth', 4);
line([1, N*rho+1], [30.5, 30.5], 'linewidth', 4);
line([1, N*rho+1], [60.5, 60.5], 'linewidth', 4);
print('-deps2c', '../documents/thesis/figs/interpolation/interp_output_neighborhood.eps');

figure; imagesc(large_out>0); axis off;
line([30.5, 30.5], [1, N*rho+1], 'linewidth', 4);
line([60.5, 60.5], [1, N*rho+1], 'linewidth', 4);
line([1, N*rho+1], [30.5, 30.5], 'linewidth', 4);
line([1, N*rho+1], [60.5, 60.5], 'linewidth', 4);
print('-deps2c', '../documents/thesis/figs/interpolation/interp_output_neighborhood_domains.eps');

