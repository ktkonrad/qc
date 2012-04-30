%% storage space
%area = 0.61 % qugrs
area = 1 + pi/4 % qust
k = 800;
dk = 100;
alpha = 0.5;

efuncs = area / pi * k * dk;
gbs = efuncs * (alpha / k)^-2 * 4 / 2^30