addpath('/altair/konrad/qc/matlab');
addpath('/altair/konrad/qc/vergini');
opts.dx = .001;
opts.v = 2;
[counts,ks,trouble_counts,wtms] = count_range2(100, 500, .1, .1, '-l qust:2 10 -s vepwoo:1.2:40:1.5', opts);
save('counts.mat');