


% mass chop from part of billiard below a line segment
%
% barnett 1/17/04

%function [k, mass] = 


k_lo = 1500;
k_hi = 1501;
sys = '-l qugrs:1:0.4:0.7 -s oyooo:1.5:7:1';
head = 'c';
b = 10;
dx = 0.01;
d_lo = 0.1;  % vergini window % 0.05
d_hi = 0.1;
verb = '';   % -q

ks = []; mass = []; ten = []; nrm = [];
ne = 0;
first_time = 1;

tic;
for k=k_lo+d_lo:d_lo+d_hi:k_hi-d_hi % .................. collect states
  
  if first_time % output geom too, including mask
    mask_flag = '-m';
    first_time = 0;
  else
    mask_flag = '';
  end
  % output flags:    -f %g:0 -e 1 -p 0
  cmd = sprintf(['verg %s %s -b %g -k %g -o %s -4 0.85 -C %g:%g:0.5:0.75'...
		 ' %s'], verb, sys, b, k, head, d_lo, d_hi, mask_flag);
  disp(cmd);
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed!');
    return;
  end

  [ks_t, mass_t, ten_t, nrm_t, ne_t] = load_sum(head);
  ks = [ks; ks_t];
  mass = [mass; mass_t];
  ten = [ten; ten_t];
  nrm = [nrm; nrm_t];
  ne = ne + ne_t;
  ten_t
  
end % ..................................................
disp(sprintf('total time = %f minutes', toc/60));
totaltime = toc;

%[mask] = load_sta([head '.mask']);
%xs = (0:nx-1)*dx;
%ys = (0:ny-1)*dx;

if 0  % geom plots
  show_geom(head);
  show_geom([head '.chop'], gcf);
end


figure;
subplot(2,1,1);
semilogy(ks, ten, '+-');
subplot(2,1,2);
plot(ks, mass, '+');
v = axis;
axis([v(1) v(2) 0 1.0]);

if 1 % Weyl plot
  [perim, area] = load_props(head);
  weyl(k_lo, ks, perim, area);
end

mean(mass)
std(mass)

%save mass_chop200.mat

if 0  % plot all k variances
  nk = 4;
  mc(4) = load('mass_chop2000.mat', 'ne', 'ks', 'mass');
  mc(3) = load('mass_chop1500.mat', 'ne', 'ks', 'mass');
  mct = load('mass_chop1501.mat', 'ne', 'ks', 'mass');
  mc(3).ne = mc(3).ne + mct.ne;
  mc(3).ks = [mc(3).ks; mct.ks];
  mc(3).mass = [mc(3).mass; mct.mass];
  mc(2) = load('mass_chop1000.mat', 'ne', 'ks', 'mass');
  mc(1) = load('mass_chop500.mat', 'ne', 'ks', 'mass');
  clear v
  for i=1:nk
    v(i) = var(mc(i).mass);
    sig_v(i) = v(i)/sqrt(mc(i).ne);
    k_mean(i) = mean(mc(i).ks);
  end
  figure;
  errorbar(k_mean, v, sig_v);
  figure;
  loglog(k_mean, v, 'o');
  line([k_mean; k_mean], [v+sig_v; v-sig_v]);
  hold on; loglog(k_mean, 0.28*k_mean.^-1, ':'); hold off;
  axis([min(k_mean) max(k_mean) 1e-4 1e-3]);
end


