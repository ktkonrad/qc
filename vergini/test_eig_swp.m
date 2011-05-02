% test the eig_swp function

%sys = '-l qugrs:1:.4:.7 -s oyooo:1.3:7:1';
sys = '-l qurf:0:.2:.2 -s oyooo:1.5:6:1';
opts.b = 15;
opts.v = 0;

%ks = 49:0.01:51;

ks = 39.1:0.01:39.6;

[mu] = eig_swp(ks, sys, opts);

figure; plot(ks, 1./sqrt(mu), '.');

