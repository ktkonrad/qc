%% storage space
area = 0.61; % qugrs
%area = 1 + pi/4; % qust
k = 2010;
dk = 10;
alpha = .7;

efuncs = area / pi * k * dk;
gbs = efuncs * (alpha / k)^-2 * 4 / 2^30