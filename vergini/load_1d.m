% function [vecs, vals, xoff, dx, param1, param2] = load_1d(head, suf)
%
% Reads in a multiple-1d-array viewer file. vecs = ne rows.
% vals = col vec.
%
% 1/4/04 barnett.

function [vecs, vals, xoff, dx, param1, param2] = load_1d(head, suf)

fid = fopen(sprintf('%s.%s', head, suf), 'r');
[a, count] = fscanf(fid, '%f', 8);
if a(1)~=1
  disp('load_1d: Not a 1d viewer file!');
  return;
end
M = a(2);
ne = a(3);
xoff = a(4);
dx = a(5);
if a(6)~=2
  disp('load_1d: Not 2 params!');
  return;
end
param1 = a(7);
param2 = a(8);

[a, count] = fscanf(fid, '%f', ne*(2+M));
if count~=ne*(2+M)
  disp('load_1d: insufficient data!');
  return
end
a = reshape(a, [2+M, ne])';
if find(a(:,1)'~=1:ne)
  disp('load_1d: state number mismatch!');
  return
end
vals = a(:,2);
vecs = a(:,3:(2+M));

fclose(fid);
