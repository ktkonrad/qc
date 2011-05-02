% function [perim, area, perim_colloc] = load_props(head)
%
% returns full perimeter, area, and perimeter of pure collocation
% bdry. The first (full) perim value is to be used for Weyl law.
%
% 2/1/04 barnett

function [perim, area, perim_colloc] = load_props(head)

name = sprintf('%s.cf', head);
fid = fopen(name, 'r');
[a, count] = fscanf(fid, '%f', 8);
fclose(fid);
if count~=8
  disp(sprintf('problem reading %s!', name));
  return;
end
perim_colloc = a(7);
area = a(8);

[x, y, nx, ny] = textread(strcat(head,'.m_pts'), '%f %f %f %f', ...
			  'headerlines', 4);
M = size(x,1); % # boundary points
perim = perim_colloc + sqrt(x(1)^2+y(1)^2) + sqrt(x(M)^2+y(M)^2);
% add lengths along symmetry
% lines (C4 symm only)

% end
