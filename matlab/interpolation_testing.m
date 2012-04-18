n = 20; % coarse grid is n x n
upsample = 20; % fine grid is n*upsample x n*upsample
dx = 1 / n; % for coarse grid
ppw = 5; % for coarse grid
k = 2*pi / (ppw * dx);

M = 8; % bessel function order

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
g = reshape(sta, n, n);
figure; imagesc(g);
figure; imagesc(g>0);
lo_res = dlmread('../data/rpw.dat');
hi_res = dlmread('../data/rpw_hi.dat');
assert(all(all(abs(lo_res - g) < 1e-6))); % make sure the files are the same rpw

%% generate a nice figure
trouble = [8, 13];
trouble_rect = [trouble - 0.5, 2, 2];
figure;

h = subplot(2,2,4);
imagesc(lo_res>0);
axis off;
set(h, 'pos', get(h, 'pos') + [-0.1,-0.1,0.1,0.1]);
rectangle('Position', trouble_rect, 'EdgeColor', 'red', 'LineWidth', 3);

h = subplot(2,2,3);
imagesc(lo_res);
axis off;
set(h, 'pos', get(h, 'pos') + [-0.1,-0.1,0.1,0.1]);
rectangle('Position', trouble_rect, 'EdgeColor', 'red', 'LineWidth', 3);

h = subplot(2,2,2);
imagesc(hi_res>0);
axis off;
set(h, 'pos', get(h, 'pos') + [-0.1,-0.1,0.1,0.1]);
rectangle('Position', trouble_rect*upsample, 'EdgeColor', 'red', 'LineWidth', 3);

h = subplot(2,2,1);
imagesc(hi_res);
axis off;
set(h, 'pos', get(h, 'pos') + [-0.1,-0.1,0.1,0.1]);
rectangle('Position', trouble_rect*upsample, 'EdgeColor', 'red', 'LineWidth', 3);

%% generate 4 separate figures
trouble = [8, 13];
trouble_rect = [trouble - 0.5, 2, 2];

figure;
imagesc(lo_res>0);
axis off;
rectangle('Position', trouble_rect, 'EdgeColor', 'red', 'LineWidth', 3);
daspect([1,1,1]);

figure;
imagesc(lo_res);
axis off;
rectangle('Position', trouble_rect, 'EdgeColor', 'red', 'LineWidth', 3);
daspect([1,1,1]);

figure;
imagesc(hi_res>0);
axis off;
rectangle('Position', trouble_rect*upsample, 'EdgeColor', 'red', 'LineWidth', 3);
daspect([1,1,1]);

figure;
imagesc(hi_res);
axis off;
rectangle('Position', trouble_rect*upsample, 'EdgeColor', 'red', 'LineWidth', 3);
daspect([1,1,1]);

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

%% generate an interpolation of the whole thing
interpolated_all = zeros(n*upsample);
addpath('../interpolation');
sten = stencil();
[xout, yout] = meshgrid(-0.5:1/upsample:0.5);
outpoints = horzcat(xout(:), yout(:));
interp_mat = interp_matrix(k*dx, sten, outpoints, M);    
for r=3:n-3
    for c=3:n-3
        interpolated_all((1:upsample+1)+(r*upsample), (1:upsample+1)+(c*upsample)) = flipud(interpolate(g, r, c, interp_mat));
    end
end
figure; imagesc(interpolated_all);