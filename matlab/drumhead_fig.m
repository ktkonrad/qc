%% constants
m = 2; n = 3;
lambda = 11.6198; % 2,3 mode
t = 1;
A = 0; B = 1; C= 0; D = 1;

r1 = 0.442; % radius of first circular nodal line
r2 = 0.724; % radius of second circular nodal line

%% polar
dr = 0.001;
dtheta = pi/240;
r = 0:dr:1;
theta = 0:dtheta:2*pi;
[rs, thetas] = meshgrid(r, theta);
[xs, ys] = pol2cart(thetas, rs);

%% rectangular
dx = 0.01;
x = -1:dx:1;
y = -1:dx:1;
[xs, ys] = meshgrid(x, y);
[rs, thetas] = cart2pol(xs, ys);


%% plotting
zs = real((A*cos(lambda*t) + B*sin(lambda*t)).*besselj(m,lambda*rs).*(C*cos(m*thetas) + D*sin(m*thetas)));

figure;
surf(xs,ys,zs, 'LineStyle', 'none');
axis off;
hold on;

% x axis
plot([-1, 1], [0, 0], '-k', 'LineWidth', 3);

% y-axis
plot([0, 0], [-1, 1], '-k', 'LineWidth', 3);

% inner cicular nodal line
[r1s, th1s] = meshgrid(r1, theta);
[x1s, y1s] = pol2cart(th1s, r1s);
plot(x1s, y1s, '-k', 'LineWidth', 3);

% outer circulr nodal line
[r2s, th2s] = meshgrid(r2, theta);
[x2s, y2s] = pol2cart(th2s, r2s);
plot(x2s, y2s, '-k', 'LineWidth', 3);

% boundary
[rbs, thbs] = meshgrid(1, theta);
[xbs, ybs] = pol2cart(thbs, rbs);
plot(xbs, ybs, '-k', 'LineWidth', 3);


%n = find(abs(zs) < 1e-3);
%plot(xs(n),ys(n),'.k');