function [ks, kos, ten, nrm, ne] = load_sum(head)

% head can include .sum now too, otherwise it adds it
if exist(head)~=2
  head = [head '.sum'];
end

[i, ks, kos, ten, nrm] = textread(head, '%d%f%f%f%f', 'headerlines', 6);
ne = length(i);
