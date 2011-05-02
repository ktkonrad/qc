% make complete list of levels of bulging-out qugrs

opts.Mo = 500;
opts.v = 0;

bil = '-l qurf:0:.05:.02 ';
sys = [bil '-s oyooo:4:2:1 -u -4 1'];
[ks, prop, err, bdry] = gen_levels(3, 30, .1, .1, sys, opts);
weyl(0, ks, prop.perim, prop.area);
%drawnow;
%figure; semilogy(ks, err.ten, '+');

sys = [bil '-s oyooo:2.5:4:1 -u -4 1'];
%sys = [bil '-s oyooo:2.5:3:1 -u -4 1'];
[ks_t, prop_t, err_t, bdry_t] = gen_levels(30, 100, .1, .1, sys, opts);  
[ks, err, bdry] = merge_levels(ks, ks_t, err, err_t, bdry, bdry_t);
%return;

sys = [bil '-s oyooo:1.5:5:1 -u -4 1'];
%sys = [bil '-s oyooo:2:4:1 -u -4 1'];
[ks_t, prop_t, err_t, bdry_t] = gen_levels(100,200,.1,.1,sys,opts);  
[ks, err, bdry] = merge_levels(ks, ks_t, err, err_t, bdry, bdry_t);

%sys = [bil '-s oyooo:1.5:6:1 -u -4 1'];
%[ks_t, prop_t, err_t, bdry_t] = gen_levels(200,300,.1,.1,sys,opts);  
%[ks, err, bdry] = merge_levels(ks, ks_t, err, err_t, bdry, bdry_t);

ne = numel(ks);
%save qurf0-300.mat ne ks err bdry prop opts sys
%save qurfh.0-200.mat ne ks err bdry prop opts sys

weyl(0, ks, prop.perim, prop.area);
Mo = opts.Mo;
Es = ks.^2;  
% compute Q
[Q] = triple(bdry.ngr, ones([1 Mo]), bdry.rn, bdry.dx, ne);
%save Q.qurfh.mat Q ne


