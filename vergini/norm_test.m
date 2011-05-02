% test area norm formulation in verg.cc

% no symmetry, test convergence with bx ( = b for area integration)
%sys = '-l rf:0:0.2:0:0 -s oyo:4:3:1';  % smooth, no symm.
%sys = '-l grs:1.4:0.3:0.7:0:0 -s oyo:4:3:1';  % corners, no symm.
sys = '-l qurf:0:0.2:0.2 -s oyoee:4:3:1';  % smooth, C_4 symm.
%sys = '-l qugrs:1.4:0.3:0.7 -s oyooo:4:3:1';  % corners, C_4 symm.
k = 20;
head = 'n';
seed = 1;                 % chooses random coefficient wavefunction
bxs = [5 7 10 15 22 30 40];
nbxs = length(bxs);
n = zeros(1,nbxs);
nrm = n;

for i = 1:nbxs
  b = bxs(i);                  % use same b as bx, converge together
  dx = 2*pi/(k*bxs(i));
  cmd = sprintf('verg -q %s -b %g -p 0 -k %g -o %s -g %g:6 -N %d -m', ...
		sys, b, k, head, dx, seed);
  disp(cmd);
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed!');
    return;
  end
  % load sta, sum, mask
  [sta, ks, ne, nx, ny] = load_sta(head);
  [ks, kos, ten, nrm(i), ne] = load_sum(head);
  [mask] = load_sta(sprintf('%s.mask', head));
  % take area norm over grid using mask...
  n(i) = triple(sta, mask, ones([1 nx*ny]), dx, ne);
end

% compare
figure;
subplot(2,1,1);plot(bxs, [nrm; n], '+-');xlabel('b and bx');
title(sprintf('norm test, k=%g, sys: %s', k, sys));
ylabel('area norms');legend('bdry integral', 'area integral');
subplot(2,1,2);loglog(bxs, abs(1 - nrm./n), '+-'); xlabel('b and bx');
ylabel('relative error in area norm');
hold on; loglog(bxs, bxs.^-2, '--'); hold off;
legend('rel err', '2nd order');

% how each version is converging...
nrm/nrm(nbxs) - 1
n/n(nbxs) - 1

% Done!

% Could do an ensemble average over random wavefuncs, but mask is
% slow to generate, so better to modify NORM_TEST task to output
% multiple funcs.

