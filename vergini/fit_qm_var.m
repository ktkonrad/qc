% fit mass-chop quantum variance data.

if 1
clear
load mc_all.mat
end

cl_mea = 0.5500;

if 0
% plot raw data...
ns = find(ks<500);
figure; set(gca, 'fontsize', 14);
h = plot(ks(ns).^2, mass(ns), '.');
set(h,'markersize', 1);
h1 = hline(cl_mea, 'r-');
%xlabel('k_j \equiv E_j^{1/2}'); ylabel('(\phi_j, A \phi_j)');
xlabel('E_n'); ylabel('<\phi_n, A \phi_n>');
%title('raw matrix element data and classical mean (A=masschop)');
legend([h, h1], 'quantum expectations', 'classical mean <A>');
text(1.2e5, 0.35, '11615 quantum states');
text(1.0e5, 0.28, 'No exceptional values -> QUE upheld');

print -dpsc2 talk/fig/raw_mc.eps
end

i = 1;  % will label chunks
nt = length(ks); % tot number states

% low set...
nch = 1000;   % chunk size
nl = length(find(ks<500));
for n=1:nch:nl  % aim for equal weight to each chunk
  ns = n:min(n+nch-1, nl);
  nm(i) = length(ns); % number in this chunk
  km(i) = mean(ks(ns));
  summ(i) = sum(mass(ns));
  sumsqm(i) = sum(mass(ns).^2);
  varm(i) = var(mass(ns));
  sigvarm(i) = varm(i)/sqrt(nm(i));
  i = i+1;
end

% isolated higher sets...
kcs = [700 1000 1500 2000 2500 3000 3500 4000];  % list of center k values
kwin = 100;        % enough to separate chunks
nkcs = length(kcs);

for j=1:nkcs
  ns = find(abs(ks-kcs(j))<kwin);
  nm(i) = length(ns);
  km(i) = mean(ks(ns));
  summ(i) = sum(mass(ns));
  sumsqm(i) = sum(mass(ns).^2);
  varm(i) = var(mass(ns));
  sigvarm(i) = varm(i)*sqrt(2/nm(i));
  i = i+1;
end

nc = i;
disp(sprintf('%d chunks collected', nc));

if 0
pfig = figure;
errorbar(km, varm, sigvarm); xlabel('k'); ylabel('qm var');
end

% fit...
mf = find(km>400);   % -------------- select which low chunks to omit
xi = [0.5 1 1];
opt = optimset('display', 'none');
%[x, f] = fminsearch('objfunc', xi, opt, ks, mass)
[x, f] = fminsearch('objfunc_chunk', ...
		    xi, opt, km(mf), summ(mf), sumsqm(mf), nm(mf))
mea = x(1);
pow = x(2);
c = x(3);
%pow = 1;
%c = 0.3;
% loglog errorbar hack...
figure; set(gca, 'fontsize', 14);
loglog(km.^2, varm, 'o');
line([km.^2; km.^2], [varm+sigvarm; varm-sigvarm]);
hold on; plot(km.^2, (c/1000)*(km./1000).^-pow, '--'); hold off;
%hold on; plot(km, c*(km./1000).^-1, ':'); hold off;
xlabel('E'); ylabel('V_A(E)');
%title('best power fit to qm variance, dotted:\gamma=1');
axis([1e4 2.02e7 7e-5 2.2e-3]);
text(1e6, 1.1e-3, 'Power law: V_A(E) = a E^{-\gamma}', 'fontsize', 14);
text(1e6, 7e-4, 'fitted \gamma = 0.48 \pm 0.01', 'fontsize', 14, ...
     'color', [0 0 1]);
text(4e4, 2e-4, '\sim 25000 quantum states');

print -depsc2 talk/fig/power_law.eps

% overlay semiclassical prediction............. no fitted params here!
[perim, area] = load_props('c');
Co = 0.109;                        % from cw_cl_masschop.m at omega=0
prefac = 2 * Co/area;  % 2 for TRS
if 1
hold on; h = plot(km.^2, prefac*(km).^-1, 'r-'); hold off;
legend(h, 'semiclassical (no fitted parameters)');
%print -depsc2 talk/fig/power_law_sc.eps
end

%return;

if 1
% likelihood plot: (hold mean const)
pows = 0.92:0.001:1.01;
cs = 0.325:0.001:0.37;
npows = length(pows)
ncs = length(cs);
clear o
for i=1:npows
  for j=1:ncs
    pow = pows(i);
    c = cs(j);
    o(i,j) = objfunc_chunk([mea pow c], km(mf), summ(mf), sumsqm(mf), nm(mf));
    %objfunc([mea pow c], ks, mass);
  end
end
mmo = min(min(o));
%figure;
%surf(pows, cs, exp(mmo-o)');
%figure(pfig);
%clabel(ch,hh);
axes('position', [.19 0.18 .35 .35]); set(gca, 'fontsize', 10);
[ch, hh] = contour(pows/2, cs, exp(mmo-o)');
%clabel(ch,hh);
%colorbar;
%xlabel('\gamma'); ylabel('a for k_0=1000'); title('likelihood at fixed \mu');
% add semiclassical prediction:
hold on; h = errorbar(1/2, prefac, prefac/500, 'r.'); hold off;
set(h,'markersize',30);
%text(0.97, 0.358, 'semiclassical estimate');
text(0.475, 0.36, 'semiclassical estimate', 'color', [1 0 0]);
xlabel('\gamma'); ylabel('a');
title('likelihood in parameter space (\gamma,a)          ');

print -depsc2 talk/fig/power_law_sc.eps

%print -depsc2 talk/fig/maxlik.eps
end

return;  % ----------------------------------------------------------------

% mean accuracy
ms = 0.549:0.00002:0.551;
nms = length(ms);
for i=1:nms
  om(i) = objfunc([ms(i) x(2) x(3)], ks, mass);
end
mo = min(om);
figure; set(gca, 'fontsize', 14);
plot(ms, exp(mo-om), '-');
h = vline(cl_mea, 'r-'); legend(h, 'classical mean <A>');

% check the model looks good... should give gaussian, look at scatter
mea = x(1);
pow = x(2);
c = x(3);
if 0 % use semiclass est
  mea = cl_mea;
  pow = 1.0;
  c = prefac;
end
ns = find(ks>400);             % select range to histogram
nns = length(ns);
vs = (c/1000)./(ks(ns)./1000).^pow;
zs = (mass(ns)-mea)./sqrt(vs);
figure;
h=plot(zs, '.');
set(h, 'markersize',1);
% histogram it:
dz = 0.2;
zcs = -5:dz:5;    % bin centers
n = hist(zs, zcs);
figure; set(gca, 'fontsize', 14);
errorbar(zcs, n, sqrt(n), '.');
hold on; plot(zcs, nns*dz*exp(-zcs.*zcs/2)/sqrt(2*pi), '-'); hold off;
%title('distribution of qm variance (after fitted 3 params)');
xlabel('z'); ylabel('counts per bin');
text(3, 1600, 'using best fit a,\gamma');

print -depsc2 talk/fig/hist.eps
