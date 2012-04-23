%        ns   y0    x1  cx
geom = '1 4 1 -1  2 1 1 1 0  0 -1 1  2 -1 -1 -1 0  0 1 -1'; % stadium
%      ccw  x0   arc  y1  cy

% ccw: counterclockwise [1 or -1]
% ns: number of segments
% (x0, y0) where to start
% type: 0 -> line segment, 1 -> circle, 2 -> arc
% (cx, cy) center of circle

geom = '1 4 0 0  0 2 0  2 1 1 1 0  0 0 1  0 0 0'; % quarter statium

N = 50;
epsilon = 1e-5;

cmd = sprintf('echo %s | ./spray -m 1000', geom);    % send 
[s,out] = system([cmd ' 2>/dev/null']);             % kill stderr output
L = strread(out, '%f', 1);
[bx,by,bnx,bny] = strread(out, '%f%f%f%f', 'headerlines', 1);
cmd = sprintf('echo %s | ./spray -T 0:%.15g:%d', geom, 4/7, N);  % single traj
[s,out] = system([cmd ' 2>/dev/null']);
[x,y,vx,vy,obj,l] = strread(out, '%f%f%f%f%d%f');
figure; plot(bx, by, 'b-'); axis equal; axis off; hold on; plot(x, y, 'k-');

cmd = sprintf('echo %s | ./spray -T 0:%.15g:%d', geom, 4/7+epsilon, N);  % single traj
[s2,out2] = system([cmd ' 2>/dev/null']);
[x2,y2,vx2,vy2,obj2,l2] = strread(out2, '%f%f%f%f%d%f');
plot(x2, y2, 'r-');

%% calculate Lyuapunov exponent

distances = sqrt((x-x2).^2 + (y-y2).^2 + (vx-vx2).^2 + (vy-vy2).^2);
figure; plot(log(distances));

last_i = find(distances > 1, 1);
t = (1:last_i)';
d = log(distances(t));

p = polyfit(t, d, 1);
hold on; semilogy(1:last_i, p(1) * t + p(2), 'r-');
text(40, -6, ['\lambda = ' p(1)]);