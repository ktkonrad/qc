addpath('/altair/konrad/qc/matlab');
addpath('/altair/konrad/qc/vergini');
opts.kdx = exp(-.5);
opts.v = 2;
[counts,ks,wtms] = count_range(200, 400, .1, .1, '-l qust:2 10 -s vepwoo:1.2:40:1.5', opts);
save('counts.mat');