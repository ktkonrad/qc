% test effect of epsilon on GEP
%
% 8/24/04 barnett

sys = '-l qurf:0:.2:.2 -s oyooo:2:5:1';
k = 50.1219;
opts.v = 0;
opts.wei = 0;
opts.ne = 100;

zs = logspace(-17,-8,20);
nz = numel(zs);
mu = zeros(opts.ne, nz);
for i=1:nz
  opts.z = zs(i);
  mu(:,i) = eig_swp(k, sys, opts);
  r(i) = numel(find(mu(:,i)~=0));   % rank
end

figure;
loglog(zs, mu, '+');
