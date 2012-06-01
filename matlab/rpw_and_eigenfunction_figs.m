k = 200;
alpha = .5;
dx = alpha/k;
ppw = 2*pi/alpha;
n = 1/dx;

%%
rpw = rpw2dsample(n, ppw);

cmd = sprintf('../vergini/verg -l qust:2 -b 10 -s vepwoo:1.2:40:1.5 -k %f -V .01 -f %f -u', k, dx);
system(cmd);

%%
addpath('../vergini');
sta = load_sta('t');
efunc = reshape(sta(1,:,:), size(sta,2), size(sta,3))';

figure; imagesc(rpw); axis off; axis equal;
print('-deps2c', '../documents/thesis/figs/intro/rpw.eps')

figure; imagesc(efunc(1:n,1:n)); axis off; axis equal;
print('-deps2c', '../documents/thesis/figs/intro/eigenfunction.eps')
