% function [mu] = eig_swp(ks, sys, opts)
%
% Generate sweep of eigenvalues of a quadratic form
%
% ks = set of wavenumbers to use
% sys = system string for verg
%
% opts.ne sets number of eigenvalues to record
% opts.b sets collocation pt density
% opts.Mo fixes bdry output and norming density
% opts.v = -2: silent, -1: quiet, 0: default, 1: verbose, 2: verg verbose.
% opts.wei = 0: w=1,   wei = 1: w=1/rn
% opts.z sets epsilon for singular GEP truncation
% opts.m is present writes geom info to files for last k of sweep
%
% 5/19/04 barnett

function [mu] = eig_swp(ks, sys, opts)

verb = '';
if isfield(opts, 'v')
  if opts.v > 1           % verbose
    verb = '-v ';
  elseif opts.v <= 0       % quiet
    verb = '-q ';
  end
else
  opts.v = 0;
end
sys = [verb sys];  % include string

if ~isfield(opts, 'head')
  opts.head = './t';  % default
end
head = ['./' opts.head];  % ensure accessing only cwd

if ~isfield(opts, 'b')
  opts.b = 15;  % default
end
b = opts.b;

if ~isfield(opts, 'wei')
  opts.wei = 0;  % default
end
wei = opts.wei;

if ~isfield(opts, 'ne')
  opts.ne = 10;  % default
end
ne = opts.ne;

% this trick lets verg decide defaults...
zstr = ' ';
if isfield(opts, 'z')
  zstr = sprintf(' -z %g', opts.z);
end

mu = zeros(ne, numel(ks));
for i = 1:numel(ks) % ................. main loop over k
  k = ks(i);
  if (i==numel(ks) & isfield(opts, 'm'))
    zstr = [zstr ' -m'];                  % write out geom
  end
  cmd = sprintf('verg %s -o %s -b %g -k %.10g -d -E %d:%d:%d%s', ...
                sys, head, b, k, wei, 1, ne, zstr);
  if opts.v>=0
    disp(cmd);
  end
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed!');
    return;
  end

  % load in eigvals...
  mu(:,i) = load_sum(head);

end % ...................................................
