% function [h] = show_sum(head, linetype)
%
% show summary info on tension and norm, returns graphics handle.

% Vergini package.
%
% Barnett 9/4/03

function [h] = show_sum(head, linetype)

if nargin<2
  linetype = '+-';
end

% load
[i, k, g, t, n] = textread(strcat(head,'.sum'), '%d %f %f %f %f', ...
			       'headerlines', 6);

% plot
h = semilogy(k,t,linetype);


