% function [I_sum, area_sum, max_rr] = mominert(head)
%
% Estimate moment of inertia from bdry colloc point geom data
% Needs files head.m_pts, head.rn and head.cf to get perim_colloc
% Also returns estimate of area, and max value of r^2 on bdry.
%
% 8/2/04 barnett

function [I_sum, area_sum, max_rr] = mominert(head)

[x, y, nx, ny] = textread(strcat(head,'.m_pts'), '%f %f %f %f', ...
			  'headerlines', 4);
M = size(x,1); % # boundary points
[rn dummy xoff dx] = load_1d(head, 'rn');  % get rn, xoff etc
rn = rn';
if (size(rn)~=size(x))
  disp(sprintf('inconsistent bdry pts numbers: M for x = %d, M for rn = %d',...
               M, numel(rn)));
end
%[perim, area, perim_colloc] = load_props(head);
%dx = perim_colloc/M;

% approximate integrals...
area_sum = dx*sum(rn)/2;        % can use to check accuracy
rr = x.^2 + y.^2;
I_sum = dx*sum(rn.*rr)/4;
max_rr = max(rr);
