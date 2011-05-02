% function [status] = save_sta(sta, ks, head)
%
% saves in sta_bin format. sta(1..ne, 1..nx, 1..ny)
%
% eg to write out a .dir file from low_efuncs, use:
% clear; load ~/bdry/sarnak/qm_diag_VB/psi.dir1.10.100.mat
% save_sta(permute(psi, [3 1 2]), sqrt(E), 'lowest');

function [status] = save_sta(sta, ks, head)

ne = size(sta, 1);
nx = size(sta, 2);
ny = size(sta, 3);
status = 1;

% binary read out...
fid = fopen(sprintf('%s.sta_bin', head), 'w');

disp(sprintf('save_sta: ne = %d (nx=%d, ny=%d)...', ne, nx, ny));
c = fwrite(fid, sprintf('bin\n'), 'uchar');
if c~=4
  status = 0;
  return;
end
c = fwrite(fid, [ne nx ny], 'int32');
if c~=3
  status = 0;
  return;
end
c = fwrite(fid, ks, 'double');
if c~=ne
  status = 0;
  return;
end
c = fwrite(fid, sta, 'single');
if c~=ne*nx*ny
  status = 0;
  return;
end

fclose(fid);
disp('done.');
