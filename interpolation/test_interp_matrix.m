n = 3; % course grid is n x n
M = 7; % use bessel functions up to order M
k = 400;
dx = .01;
upsample = 20;

xmin = -(n-1)/2*dx;
xmax = (n-1)/2*dx;

x = xmin:dx:xmax;
xs = meshgrid(x);
ys = flipud(meshgrid(x)');
points = [xs(:) ys(:)]; % coarse grid

x2 = xmin:dx/upsample:xmax;
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
cond(BAplus)