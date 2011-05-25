% function [counts, ks, wtms, prop, err, bdry, grid] = ...
%    count_range(k_lo, k_hi, delta_lo, delta_hi, sys, opts);
%
% Generate set of Dirchlet levels and eigenfunctions with fixed basis params
%
% opts.b changes collocation pt density from default 15
% opts.kdx changes k*dx from default exp(0.5)
% opts.head changes default file head 't'
% opts.v = -2: silent, -1: quiet, 0: default, 1: verbose, 2: verg verbose.
%
% only outputs bdry values if asked for in output args
% only calculates and outputs grid values if asked for in output args
%
% prop contains billiard properties, CPU time etc
%
% Based on gen_levels.m by Alex Barnett
%
% Kyle Konrad 5/12/2011

function [counts, ks, wtms] = ...
    count_range(k_lo, k_hi, delta_lo, delta_hi, sys, opts)

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

if ~isfield(opts, 'kdx')
  opts.kdx = exp(-.5);  % default
end
kdx = opts.kdx;


ks = []; counts = []; wtms = [];
ne = 0;

tic; % start timer
  
%.................................main loop over k...................
for k=(k_lo+delta_lo):(delta_lo+delta_hi):(k_hi+delta_lo)

    dx = kdx / k;
    
    outopts = sprintf('-f %g:0 %s', dx);
    
    cmd = sprintf('verg -q %s -o %s -b %g -k %g -V %g:%g %s', ...
                  sys, head, b, k, delta_lo, delta_hi, outopts);
    if opts.v>=0
      disp(cmd);
    end
    
    [s] = system(cmd);
    if (s~=0)
      disp('verg crashed!');
      return;
    end    

    ks_t = load_sum(head);
    j = find(ks_t<=k_hi);
    if opts.v>=-1
      disp(sprintf('\tfound %d', numel(j)));
    end
    ne = ne + numel(j);

    % only works for quarter stadium currently
    cmd = sprintf('../c/count -f %s.sta_bin -l qust:2 -k %g -d %g', head, k, dx);

    if opts.v>=0
      disp(cmd);
    end
    
    [s o] = system(cmd);

    out = regexp(o, '(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', 'tokens');

    if (s ~= 0)
        disp('count crashed!');
        disp(o);
        return;
    end

    counts = [counts cellfun(@(x) str2double(x(1)), out)];
    ks = [ks cellfun(@(x) str2double(x(2)), out)];
    wtms = [wtms cellfun(@(x) str2double(x(3)), out)];


end % ..................................................
  
prop.toc = toc;
if opts.v>=-1
  disp(sprintf('total time = %f minutes', toc/60));
end

% end
