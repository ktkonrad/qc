%        ns   y0    x1  cx
geom = '1 4 1 -1  2 1 1 1 0  0 -1 1  2 -1 -1 -1 0  0 1 -1'; % stadium
%      ccw  x0   arc  y1  cy

% ccw: counterclockwise [1 or -1]
% ns: number of segments
% (x0, y0) where to start
% type: 0 -> line segment, 1 -> circle, 2 -> arc
% (cx, cy) center of circle

geom = '1 4 0 0  0 2 0  2 1 1 1 0  0 0 1  0 0 0'; % quarter statium

cmd = sprintf('echo %s | ./spray -m 1000', geom);    % send 
[s,out] = system([cmd ' 2>/dev/null']);             % kill stderr output
L = strread(out, '%f', 1);
[bx,by,bnx,bny] = strread(out, '%f%f%f%f', 'headerlines', 1);
cmd = sprintf('echo %s | ./spray -T 0:%.15g:50', geom, 4/7);  % single traj
[s,out] = system([cmd ' 2>/dev/null']);
[x,y,vx,vy,obj,l] = strread(out, '%f%f%f%f%d%f');
figure; plot(bx, by, 'b-'); axis equal; axis off; hold on; plot(x, y, 'k-');
%figure; semilogy(abs(vx.^2+vy.^2-1));   % grows by factor 10 each arc bounce!
%print -depsc2 stadium_orbit.eps

%% plot a second nearby trajectory
cmd = sprintf('echo %s | ./spray -T 0:%.15g:50', geom, 4/7+.00001);  % single traj
[s,out] = system([cmd ' 2>/dev/null']);
[x2,y2,vx,vy,obj,l] = strread(out, '%f%f%f%f%d%f');
hold on; plot(x2, y2, 'r-');
%print('-deps2c', '../documents/thesis/figs/classical/stadium_divergent_orbits_early.eps');

%% plot eigenfuction
addpath('../vergini');
[s, ks] = load_sta('../c/t');
f = reshape(s(1,:,:), size(s,2), size(s,3))';
k_1 = ks(1);
k = 70;
scale = k_1/k;
dx = 0.7/k;
x = 0:dx:2; y = 0:dx:1;
figure; imagesc(x, y, f); caxis([-2.5, 2.5]);
set(gca, 'ydir', 'normal'); axis equal; axis off; hold on;
plot(bx, by, 'b-');
plot(bx.*scale, by.*scale, 'k-');
%print -depsc2 ../documents/thesis/figs/classical/stadium_eigenfunction.eps

%% plot nodal domains
f = dlmread('counted.dat');
f(abs(f)>10000) = 0;
k = 700.106;
dx = 0.7/k;
x = 0:dx:2; y = 0:dx:1;
figure; imagesc(x, y, (f==-1 | f==0)*2 - (f==0)); set(gca, 'ydir', 'normal'); axis equal; axis off; hold on; plot(bx, by, 'b-');
colormap(flipud(hot));
print -depsc2 ../documents/thesis/figs/classical/stadium_eigenfunction_largest_nodal_domain.eps
