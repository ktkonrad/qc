n = 20; % coarse grid is n x n
upsample = 20; % fine grid is n*upsample x n*upsample
dx = 1 / n; % for coarse grid
ppw = 5; % for coarse grid
k = 2*pi / (ppw * dx);

% k = 2*pi / lambda
% ppw = lambda / dx = 2*pi / k*dx
% k = 2*pi / (ppw *dx)

%% generate a rpw
f = rpw2dsample(n*upsample, ppw*upsample);
g = f((1:n)*upsample, (1:n)*upsample);
%g_averaged = filter2(fspecial('average', upsample), f);
figure; imagesc(f>0);
figure; imagesc(g>0);


%% output rpw to text files
dlmwrite('../data/rpw.dat', g, 'precision', '%.16g');
dlmwrite('../data/rpw_hi.dat', f, 'precision', '%.16g');

%% output rpw to sta
sta=zeros(1,n,n);
sta(1,:,:) = g;
addpath('../vergini/');
save_sta(sta, k, '../data/rpw');

%sta(1,:,:) = g_averaged;
%save_sta(sta, k, '../data/rpw_averaged');

%% read rpw from sta and text files
addpath('../vergini/');
sta = load_sta('../data/rpw');
g = reshape(a, 20, 20);
figure; imagesc(b);
figure; imagesc(b>0);
lo_res = dlmread('../data/rpw.dat');
hi_res = dlmread('../data/rpw_hi.dat');
assert(all(all(abs(lo_res - b) < 1e-6))); % make sure the files are the same rpw

%% run count
cmd = sprintf('../c/count -f ../data/rpw.sta_bin -d %f -M 9 -u 20 -k %f', dx, k);
system(cmd);
c = load('../data/counted.dat')';
c(c>100000) = 0;
figure; imagesc(c);
figure; imagesc(hi_res>0);


%% compare interpolated with fine grid
sub = hi_res(260:280, 160:180);

interp = fliplr(read_dumped('../data/interpolated_12_7.dat', upsample+1, upsample+1));
figure; imagesc(interp); title('interpolated');
figure; imagesc(sub); title('high resolution');
figure; imagesc(interp>0); title('interpolated');
figure; imagesc(sub>0); title('high resolution');
figure; imagesc(interp - sub); title('error');