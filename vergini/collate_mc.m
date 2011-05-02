% Collate mass chop data into 1 big .mat file

% start cold
ks = []; mass = []; ten = []; nrm = [];
ne = 0;

% add shell-script calculated variety...
mcs = {'mc1s_all100', 'mc2s_all200', 'mc3s_all300', 'mc4s_all400',...
      'mco_all650', 'mcs_all950', 'mc2h_all2500', 'mc1h_all3000', ...
       'mc3h_all3500', 'mc4h_all4000'};
%mcs = {'t_all100'};
%mcs = {'mc_all100', 'mc_all200', 'mc_all300', 'mc_all400'};
%mcs = {'mco_all650'};
for mc = mcs
  mc = char(mc);
  disp(['collating ' mc '...']);
  [ks_t, mass_t, ten_t, nrm_t, ne_t] = load_mc(mc);
  ks = [ks; ks_t];
  mass = [mass; mass_t];
  ten = [ten; ten_t];
  nrm = [nrm; nrm_t];
  ne = ne + ne_t;
end

% add .mat file variety...
%mats = {};
mats = {'mass_chop1500', 'mass_chop1501', 'mass_chop2000'};

for mc = mats
  mc = char(mc);
  disp(['collating ' mc '...']);
  tmp = load(mc);
  ks = [ks; tmp.ks];
  mass = [mass; tmp.mass];
  ten = [ten; tmp.ten];
  nrm = [nrm; tmp.nrm];
  ne = ne + tmp.ne;
end
% use last one to get props...
%perim = tmp.perim;
%area = tmp.area;
[perim, area] = load_props('c');

% check via plot...
figure;
h = plot(ks, mass, '.');
set(h,'markersize', 1);

% monotonicity of ks...
figure;
plot(ks);

% weyl tests on chunks (don't know why so wobbly for low k - missing+duplicates)
weyl(100, ks(find(ks<500)), perim, area);
print -depsc2 talk/fig/mc100_weyl.eps
weyl(650, ks(find(ks>650 & ks<750)), perim, area);
print -depsc2 talk/fig/mc650_weyl.eps
weyl(650, ks(find(ks>650 & ks<750)), perim, area, 0.15);
print -depsc2 talk/fig/mc650_Nweyl.eps


% demo the duplicates sepaterated by uniform amounts (due to cubic
% k correction): edge at 4.6e-4
figure; plot(diff(ks(find(ks<500)))); v=axis; axis([v(1) v(2) 0 1e-3]);

if 0
weyl(950, ks(find(ks>950 & ks<1050)), perim, area);

figure; plot(diff(ks(find(ks>950 & ks<1050))));
v=axis; axis([v(1) v(2) 0 1e-3]);
ns = find(ks>950 & ks<1050);
nl = ns(1);
nh = ns(length(ns));  % nl, nh gives index range for k~1000 states
nc = nl - 1 + find(diff(ks(nl:nh))<1.1e-4); % set of duplicated states
figure; plot(ks(nc), [mass(nc) mass(nc+1)], '+');
end

weyl(2500, ks(find(ks>2500 & ks<2600)), perim, area);
weyl(3000, ks(find(ks>3000 & ks<3100)), perim, area);
weyl(3500, ks(find(ks>3500 & ks<3600)), perim, area);
weyl(4000, ks(find(ks>4000 & ks<4100)), perim, area);

figure; semilogy(ten);

save mc_all.mat ne ks mass ten nrm perim area

