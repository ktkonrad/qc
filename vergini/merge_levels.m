% function [ks, err, bdry, g] = merge_levels(ks1, ks2, err1, err2 ...
%    {, bdry1, bdry2 {, g1, g2}})
%
% merge 2 sets of level information
%
% 6/27/04 barnett

function [ks, err, bdry, g] = merge_levels(ks1, ks2, err1, err2, bdry1, ...
                                           bdry2, g1, g2)

ks = [ks1; ks2];
err.ten = [err1.ten; err2.ten];
err.nrm = [err1.nrm; err2.nrm];
err.kos = [err1.kos; err2.kos];

if nargout>2   % want bdry info
  bdry = bdry1;  % basic info
  if numel(bdry1.rn)~=numel(bdry2.rn)
    disp('merge_levels: mismatch in rn sizes in bdry.rn!');
  end
  bdry.ngr = [bdry1.ngr; bdry2.ngr];
  bdry.per = [bdry1.per; bdry2.per];
end

if nargout>3   % want grid info
  g = g1;
  if g1.dx~=g2.dx
    disp('merge_levels: mismatch in dx in grid.dx!');
  end
  g.sta = [g1.sta; g2.sta];
end
