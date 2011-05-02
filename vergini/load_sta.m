% function [sta, k, ne, nx, ny] = load_sta(head)
%
% Read file head.sta_bin into sta array (ne by nx by ny), and k (ne).
% Note convention: nx,ny are numbers of samples (1...nx, 1...ny).
% barnett 11/17/03
% Endian sensing added barnett 10/17/04

function [sta, k, ne, nx, ny] = load_sta(head)

% binary read in...
fid = fopen(sprintf('%s.sta_bin', head), 'r');
txt = fread(fid, 4, 'uchar');
ne = fread(fid, 1, 'int32');
if (ne<0 | ne>1e5)
  disp(sprintf('load_sta: ne = %d, suspect wrong endian, retrying...', ne));
  fclose(fid);
  fid = fopen(sprintf('%s.sta_bin', head), 'r', 'ieee-be');
  txt = fread(fid, 4, 'uchar');
  ne = fread(fid, 1, 'int32');
end
nx = fread(fid, 1, 'int32');
ny = fread(fid, 1, 'int32');
disp(sprintf('load_sta: ne = %d (nx=%d, ny=%d)', ne, nx, ny));
k = fread(fid, ne, 'double');
n = ne*nx*ny;
[sta, count] = fread(fid, n, 'float32');
if count~=n
  disp('load_sta: not enough binary data!');
end
fclose(fid);

sta = reshape(sta, [ne nx ny]);
