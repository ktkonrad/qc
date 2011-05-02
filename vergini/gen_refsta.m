% Generate reference set of dirchlet eigenfunctions.
%
% currently bodry output uses Mo=M, set in sys argument
%
% g is struct giving sta output on grid, with xs, ys, mask.
% Only calcs sta output on grid if g needed.
%
% 3/19/04 barnett

function [ks, bdry, prop, err, grid] = ...
    gen_refsta(k_lo, k_hi, delta_lo, delta_hi, delta_ref, sys, head, b_g, Mo)

want_sta = nargout>4;
if want_sta
  if nargin<8
    b_g = 6;
  end
  dx = 2*pi/(k_hi*b_g); % compute dx using b
end
if nargin<7
  head = 't';     % temp head file
end
if nargin<9
  Mo = 500;     % choose default Mo=500 (need perim to know this)
end

sta = []; ks = []; kos = []; ten = []; nrm = []; per = []; ngr = [];
ne = 0;
first_time = 1;

tic;
for k=(k_lo+delta_lo):(delta_lo+delta_hi):(k_hi-delta_hi) %.... loop over k

  if first_time % output geom too, including mask
    mask_flag = '-m';
    first_time = 0;
  else
    mask_flag = '';
  end
  if want_sta
    cmd = sprintf('verg %s -k %g -o %s -p %d -g %g:10 -R %g:%g:%g %s', sys, ...
                  k, head, Mo, dx, delta_lo, delta_hi, delta_ref, mask_flag);
  else
    cmd = sprintf('verg %s -k %g -o %s -p %d -R %g:%g:%g', sys, k, ...
		head, Mo, delta_lo, delta_hi, delta_ref);
  end
  disp(cmd);
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed!');
    return;
  end
  
  % load sta, sum, ngr, per
  if want_sta
    [sta_t, ks_t, ne_t, nx, ny] = load_sta(head);
    sta = [sta; sta_t];
  end
  [ks_t, kos_t, ten_t, nrm_t, ne_t] = load_sum(head);
  per_t = load_1d(head, 'per');
  ngr_t = load_1d(head, 'ngr');

  ks = [ks; ks_t];
  kos = [kos; kos_t];
  ten = [ten; ten_t];
  nrm = [nrm; nrm_t];
  per = [per; per_t];
  ngr = [ngr; ngr_t];
  ne = ne + ne_t;
end % ..................................................
disp(sprintf('total time = %f minutes', toc/60));

% load properties
[prop.perim, prop.area] = load_props(head);

% ---------------- package stuff --------------------
if want_sta
  [grid.mask] = load_sta(sprintf('%s.mask', head));       % load mask
  grid.xs = (0:nx-1)*dx;  % grid of cell-center x,y values
  grid.ys = (0:ny-1)*dx;
  grid.nx = nx;
  grid.ny = ny;
  grid.dx = dx;
  grid.sta = sta;
  save_sta(sta, ks, head);   % write all states back onto head file
end

bdry.ngr = ngr;
bdry.per = per;
[bdry.rn dummy bdry.xoff bdry.dx] = load_1d(head, 'rn');  % get rn, xoff etc

% also write per, ngr back to head file... ?

err.ten = ten;
err.nrm = nrm; 

% end
