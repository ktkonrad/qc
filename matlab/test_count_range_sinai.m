addpath('~/qc/matlab');
addpath('~/qc/vergini');
opts.kdx = exp(-.5);
opts.v = 2;
[counts,ks,wtms] = count_range_sinai(300, 302, .1, .1, ...
                   '-l qugrs:1:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1', opts);
save('counts.mat');
