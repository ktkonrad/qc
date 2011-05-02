% function [mean_m, counts, bin_ctrs] = band_prof(k, m, dk, kmax {, k_off})
%
% compute band profile of matrix m (which has already squared elements)
%
% barnett 11/23/03
% 2/19/04 added optional offset for left edge of bins.

function [mean_m, counts, bin_ctrs] = band_prof(k, m, dk, kmax, k_off)

if nargin<5
  k_off = 0;       % default offset
end

% define bins
bin_edg = k_off + (0:dk:(kmax+dk)); % bin edges
nb = length(bin_edg) - 1;
bin_ctrs = bin_edg(1:nb) + dk/2;
if k_off==0
  bin_edg(1) = 1e-6; % don't count diagonal, even if k_off=0
end

% make k-differences array
ne = length(k);
kap = kron(ones(size(k)), k') - kron(k, ones(size(k))');
%figure;imagesc(k, k, kap');

% make index containing which bin each entry in m belongs in
[c, index] = histc(kap, bin_edg);
counts = histc(reshape(kap, [ne*ne 1]), bin_edg);
counts = counts(1:nb)';  % drop the last count, which is # not in bins

% avg m in each bin
mean_m = zeros(1, nb);
for i=1:nb
  mean_m(i) = sum(m(find(index==i)))/counts(i);
end
