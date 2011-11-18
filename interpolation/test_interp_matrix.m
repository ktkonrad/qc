n = 4; % course grid is n x n
M = 8; % use bessel functions up to order M
k = 400;
dx = .001;
upsample = 30;

xmin = -(n-1)/2*dx;
xmax = (n-1)/2*dx;

x = xmin:dx:xmax;
xs = meshgrid(x);
ys = flipud(meshgrid(x)');
points = [xs(:) ys(:)]; % coarse grid

x2 = -.5*dx:dx/upsample:.5*dx;
xs2 = meshgrid(x2);
ys2 = flipud(meshgrid(x2)');
points2 = [xs2(:) ys2(:)]; % fine grid

cs = normrnd(0,1,[2, 10]);
theta_ks = 2*pi*rand(1,10);
k_vecs = [k * cos(theta_ks) ; k * sin(theta_ks)];

f_vals = rpw(cs, k_vecs, points); % random plane waves
f_vals2 = rpw(cs, k_vecs, points2);

[BAplus, preconditioner] = interp_matrix(k, points, points2, M);
interpolated = BAplus * f_vals;
error_norm = max(abs(interpolated - f_vals2))
interpolated2 = bessel_interp2(k, points, f_vals, M, points2);
error_norm2 = max(abs(interpolated2 - f_vals2))

figure;imagesc(reshape(f_vals, sqrt(length(f_vals)), sqrt(length(f_vals))))
figure;imagesc(reshape(f_vals2, sqrt(length(f_vals2)), sqrt(length(f_vals2))))
figure;imagesc(reshape(interpolated, sqrt(length(interpolated)), sqrt(length(interpolated))))
figure;imagesc(reshape(interpolated2, sqrt(length(interpolated2)), sqrt(length(interpolated2))))
%figure;imagesc(reshape(abs(f_vals2 - interpolated), sqrt(length(interpolated)), sqrt(length(interpolated))))