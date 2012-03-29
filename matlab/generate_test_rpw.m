k = 200;
ppw = 5;

n = 20;
upsample = 20;


%%
f = rpw2dsample(n*upsample, ppw*upsample);
g = f((1:n)*upsample, (1:n)*upsample);
g_averaged = filter2(fspecial('average', upsample), f);
figure; imagesc(f>0);
figure; imagesc(g>0);



% k = 2*pi / lambda
% ppw = lambda / dx = 2*pi / k*dx

%%

dlmwrite('../data/rpw.dat', g, 'precision', '%.16g');
dlmwrite('../data/rpw_hi.dat', f, 'precision', '%.16g');

%%
sta=zeros(1,n,n);
sta(1,:,:) = g;
addpath('../vergini/');
save_sta(sta, k, '../data/rpw');

sta(1,:,:) = g_averaged;
save_sta(sta, k, '../data/rpw_averaged');