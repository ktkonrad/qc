% function [k, t, per, ngr, sta] = efunc_window(sys, k {, delta {, delta_hi}})
%
% sys:   string giving billiard and basis options.
% k:        center wavenumber
% delta:    wavenumber window half-width

function [k, t, per, ngr, sta] = efunc_window(sys, k, delta_lo, delta_hi)

if nargin==3
  delta_hi = delta_lo;
end

if nargout==2 % k and t only
  
  [s,w] = system(sprintf('verg %s -k %f -V %f:%f', ...
			 sys, k, delta_lo, delta_hi));
  
  % load t.sum
  
  
end

