% test vergini modules, peanut

theta = 0.4;
k = 50;
delta_lo = 0.1;
delta_hi = 0.2;
delta_ref = 0.005;
eta = 2.5;
kD = 3;
head = 't';
dx = 1/k;
perim = 3.16; % approx
b = 6;
oM = b * k*perim/(2*pi);
%disp(sprintf('dx corresponds to b = %g', 2*pi*k*dx));    % bug in b here!
sys = sprintf('-l qugrs:1:%g:%g -b 10 -s oyooo:%g:%g:1', ...
	      theta, -(pi/2 + theta), eta, kD);

sta = [];
ks = [];
kos = [];
ten = [];
nrm = [];
ne = 0;
first_time = 1;

tic;
for k=50:(delta_lo+delta_hi):52 %................... loop over k

  if first_time % output geom too, including mask
    mask_flag = '-m';
    first_time = 0;
  else
    mask_flag = '';
  end
  cmd = sprintf('verg %s -k %g -o %s -g %g:10 -R %g:%g:%g %s', sys, k, ...
		head, dx, delta_lo, delta_hi, delta_ref, mask_flag);
  disp(cmd);
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed!');
    return;
  end
  
  % load sta, sum
  [sta_t, ks_t, ne_t, nx, ny] = load_sta(head);
  [ks_t, kos_t, ten_t, nrm_t, ne_t] = load_sum(head);
  
  sta = [sta; sta_t];
  ks = [ks; ks_t];
  kos = [kos; kos_t];
  ten = [ten; ten_t];
  nrm = [nrm; nrm_t];
  ne = ne + ne_t;
end % ..................................................
disp(sprintf('total time = %f minutes', toc/60));

% load mask
[mask] = load_sta(sprintf('%s.mask', head));
xs = (0:nx-1)*dx;
ys = (0:ny-1)*dx;
if 0
  figure;
  imagesc(xs, ys, squeeze(mask(1,:,:))'); caxis([0 1]); colorbar;
  show_geom(head, gcf);
end

if 0
  figure;
  imagesc(xs, ys, squeeze(sta(1,:,:))'); caxis([-2 2]);
  show_geom(head, gcf);
end
  
if 0 % do ten vs delta plot
  figure;
  loglog(abs(ks-kos), ten, '+');
  title('tension vs delta');
end

%save qugrs.mat
%save_sta(sta, ks, 't');

% test orthonormality
%ne = 50; % choose subset
k = ks(1:ne);
o = triple(sta, mask, ones([1 nx*ny]), dx, ne);
figure;
imagesc(k, k, log10(abs(o - eye(ne)))'); caxis([-8 -1]); colorbar;
title('log10 of orthonormality error');

% compute overlaps with some function
yi = kron(ys, ones(1,nx));
xi = kron(ones(1,ny), xs);
V = xi.*xi + yi.*yi;
if 0
  figure;
  imagesc(xs, ys, reshape(V(1,:).*mask(1,:), [nx ny])');
  show_geom(head, gcf); colorbar;
end
m = triple(sta, mask, V, dx, ne);
figure;
imagesc(k, k, (m.*m)'); colorbar;

% compute Cqm_ii(kappa)
dk = 0.1;
[cqmvv, counts, kap] = band_prof(k, m, dk, 2);
figure; errorbar(kap, cqmvv, cqmvv./sqrt(counts), '-');
title('C^{(qm)}_{VV}(\kappa)'); xlabel('\kappa'); ylabel('C(\kappa)');


% end
