% function [ks, prop, err, bdry, grid] = ...
%    gen_levels(k_lo, k_hi, delta_lo, delta_hi, sys, opts);
%
% Generate set of Dirchlet levels and eigenfunctions with fixed basis params
%
% opts.b changes collocation pt density from default 15
% opts.Mo changes bdry output and norming density from default 6
% opts.dx changes dx from default 0.01
% opts.head changes default file head 't'
% opts.sweep_dE uses sweep (of step size dE) rather than vergini
% opts.v = -2: silent, -1: quiet, 0: default, 1: verbose, 2: verg verbose.
%
% only outputs bdry values if asked for in output args
% only calculates and outputs grid values if asked for in output args
%
% prop contains billiard properties, CPU time etc
%
% 6/26/04 barnett
% 6/28/04 changed so if interrupt, data already present in output structs
%         (don't leave packaging until the end)
% 6/19/05 option to use sweep search instead of vergini (much slower!)

function [ks, prop, err, bdry, grid] = ...
    gen_levels(k_lo, k_hi, delta_lo, delta_hi, sys, opts)

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

want_bdry = nargout>3;
want_grid = nargout>4;
if opts.v>=-1
  disp(sprintf('want bdry = %d', want_bdry));
  disp(sprintf('want grid = %d', want_grid));
end

if ~isfield(opts, 'head')
  opts.head = './t';  % default
end
head = ['./' opts.head];  % ensure accessing only cwd

if ~isfield(opts, 'b')
  opts.b = 15;  % default
end
b = opts.b;

if ~isfield(opts, 'dx')
  opts.dx = 0.01;  % default
end
dx = opts.dx;

if ~isfield(opts, 'sweep_dE')
  opts.sweep_dE = 0;    % means use vergini, no sweep
end
sweep_dE = opts.sweep_dE;
  
outopts = '';
maskoutopts = '';
if want_grid
  if (sweep_dE==0)
    outopts = [outopts sprintf('-g %g:0', dx)];
  else
    outopts = [outopts sprintf('-f %g:0', dx)];
  end
  maskoutopts = sprintf('-g %g:10', dx);   % for mask grid
end

% output geom, and mask if wanted, via dummy task...
cmd = sprintf('verg %s -o %s -k %g -T 1:1 -m %s', ...
              sys, head, k_hi, maskoutopts);
if opts.v>=-1
  disp(cmd);
end
[s] = system(cmd);
if (s~=0)
  disp('verg crashed doing geom/mask output!');
  return;
end
% load properties
[prop.perim, prop.area, prop.perim_colloc] = load_props(head);

% need to define a default Mo so it stays fixed while k changes...
if ~isfield(opts, 'Mo')
  bo = 6;   % default
  opts.Mo = floor(bo*k_hi*prop.perim_colloc/(2*pi) + 1);
  if opts.v>=-1
    disp(sprintf('choosing Mo = %d, based on default b = %f', opts.Mo, bo));
  end
end
Mo = opts.Mo;

if want_bdry
  outopts = [outopts sprintf(' -p %g', Mo)];
  maskoutopts = [maskoutopts sprintf(' -p %g', Mo)];
  cmd = sprintf('verg %s -o %s -k %g -T 1:1 -m %s', ...
                sys, head, k_hi, maskoutopts);
  if opts.v>=-1
    disp(cmd);
  end
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed doing geom for rn output!');
    return;
  end
  [bdry.rn dummy bdry.xoff bdry.dx] = load_1d(head, 'rn');  % get rn, xoff etc
end

if want_grid
  [grid.mask, dummy_k, dummy_ne, nx, ny] = load_sta(sprintf('%s.mask', ...
                                                    head));       % load mask
  grid.xs = (0:nx-1)*dx;  % grid of cell-center x,y values
  grid.ys = (0:ny-1)*dx;
  grid.nx = nx;
  grid.ny = ny;
  grid.dx = dx;
end

ks = []; err.kos = []; err.ten = []; err.nrm = [];
if want_bdry
  bdry.per = []; bdry.ngr = [];
end
if want_grid
  grid.sta = [];
end
ne = 0;

tic; % start timer
  
if sweep_dE==0
    
  %.................................main loop over k...................
  for k=(k_lo+delta_lo):(delta_lo+delta_hi):(k_hi+delta_lo)
    
    cmd = sprintf('verg %s -o %s -b %g -k %g -V %g:%g %s', ...
                  sys, head, b, k, delta_lo, delta_hi, outopts);
    if opts.v>=0
      disp(cmd);
    end
    [s] = system(cmd);
    if (s~=0)
      disp('verg crashed!');
      return;
    end    
    
    [ks_t, kos_t, ten_t, nrm_t, ne_t] = load_sum(head);
    j = find(ks_t<=k_hi);
    if opts.v>=-1
      disp(sprintf('\tfound %d', numel(j)));
    end
    ne = ne + numel(j);
    
    ks = [ks; ks_t(j)];
    err.kos = [err.kos; kos_t(j)];
    err.ten = [err.ten; ten_t(j)];
    err.nrm = [err.nrm; nrm_t(j)];
    if want_bdry        % load ngr, per
      per_t = load_1d(head, 'per');
      ngr_t = load_1d(head, 'ngr');
      bdry.per = [bdry.per; per_t(j,:)];
      bdry.ngr = [bdry.ngr; ngr_t(j,:)];
    end
    if want_grid        % load sta, sum, ngr, per
      [sta_t, ks_t, ne_t, nx, ny] = load_sta(head);
      grid.sta = [grid.sta; sta_t(j,:,:)];
    end
  end % ..................................................

else
  
  %............................... main sweep over E using search........
  for E = k_lo^2:sweep_dE:k_hi^2
    
    klo = sqrt(E);
    khi = sqrt(E+sweep_dE);
    topts.sys = sys; topts.wei = 0; topts.v = '-v';
    oopts = optimset('fminbnd'); oopts = optimset(oopts, 'TolFun', 1e-6); % -15
    oopts = optimset(oopts, 'TolX', 1e-4);
    oopts = optimset(oopts, 'display', 'off');
    [kmin, tmin, flag, outp] = fminbnd('ten', klo, khi, oopts, topts);
    % hack to ignore case of active constraints (hits window edge):
    %if (abs(kmin-klo)>1e-2 & abs(kmin-khi)>1e-2)
    if (flag==1 & (tmin/kmin^2)<3e-3 & abs(kmin^2-klo^2)>1e-2 & ...
        abs(kmin^2-khi^2)>1e-2)
      % hack to say if found a state: cutoff causes missing, duplicates, etc
      if opts.v>=-1
        disp(sprintf('found 1 at: Emin=%g, t(E)_min=%g', kmin^2, sqrt(tmin)));
      end
      cmd = sprintf('verg %s -o %s -b %g -k %g -E %d:1:1 %s', ...
                    sys, head, b, kmin, topts.wei, outopts);
      if opts.v>=0
        disp(cmd);
      end
      [s] = system(cmd);
      if (s~=0)
        disp('verg crashed!');
        return;
      end    
      ks = [ks; kmin];
      err.kos = [err.kos; kmin];
      err.ten = [err.ten; tmin];
      
      if want_bdry
        per_t = load_1d(head, 'per');
        ngr_t = load_1d(head, 'ngr');
        bdry.per = [bdry.per; per_t(1,:)];
        bdry.ngr = [bdry.ngr; ngr_t(1,:)];
      end
      if want_grid
        [sta_t, ks_t, ne_t, nx, ny] = load_sta(head);
        grid.sta = [grid.sta; sta_t(1,:,:)];
      end
      ne = ne+1;
    end
    
  end
end
  
prop.toc = toc;
if opts.v>=-1
  disp(sprintf('total time = %f minutes', toc/60));
end

% end
