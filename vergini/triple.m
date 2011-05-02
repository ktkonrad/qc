% function [m] = triple(sta, mask, V, dx, ne {, diag})
%
% Compute matrix of integrals M_ij = int_omega d^2r psi_i(r) V(r) psi_j(r)
% using grid representation and domain mask, for first ne wavefuncs
% in sta array. Seems to only give 10 Mflops on fricka for larger nx*ny,
% but 40 Mflops for smaller.
%
% If diag=1, only return diagonal entries, size(m) = [1, ne].
%
% barnett 11/23/03

function [m] = triple(sta, mask, V, dx, ne, diag)

if nargin==5
  diag = 0;
end
if diag
  disp('triple overlaps: diag only');
  m = zeros(ne, 1);
else
  disp('triple overlaps:');
  m = zeros(ne, ne);
end

ne_full = size(sta, 1);
nx = size(sta, 2);
ny = size(sta, 3);
V = reshape(V(:).*mask(:), [1 nx*ny]);

frac_done = 0; % for full m case
tic;
for i=1:ne
  f = (i-1)/ne;
  if ~diag & floor(20*f*(2-f))~=floor(20*frac_done) 
    disp(sprintf('%.1f%%', 100*f*(2-f)));
  end
  frac_done = f*(2-f);
  if diag
    m(i) = dx*dx*sum(sta(i,:).*V.*sta(i,:));
  else
    for j=i:ne
      m(i,j) = dx*dx*sum(sta(i,:).*V.*sta(j,:));
      m(j,i) = m(i,j);
    end
  end
end
if diag
  n = nx*ny*ne;
else
  n = nx*ny*ne*(ne+1)/2;
end
disp(sprintf('%d mults and %d adds done in %f secs = %.1f Mflops',...
	     2*n, n, toc, 1e-6*3*n/toc));
