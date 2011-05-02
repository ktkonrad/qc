function [list, vals] = src(ns, t, f, e)
% function [list, vals] = src(ns, t, f, e)
%
% ns = # samples (to prevent overrun)
% t = time of impulsive source
% f = smaple freq
% e = gaussian std dev factor = impulse width in sample units
% list = which values in signal to add to
% vals = values to add to this list

s = ceil(e*e);
to = t*f;
list = floor(to - s):ceil(to + s);
list = list(find(list<=ns & list>=1));   % prevent overrun
vals = exp(-(list-to).^2/(2*e*e));
