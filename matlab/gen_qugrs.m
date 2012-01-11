  sys = '-l qugrs:1:0.4:0.7 -s oyooo:2:5:1 -u -4 1';
  opts.v = 0;
  [ks, prop, err, bdry] = gen_levels(100,200,.1,.1,sys,opts);  

  sys = '-l qugrs:1:0.4:0.7 -s oyooo:2.5:4:1 -u -4 1';
  opts.v = 0;
  opts.Mo = 404;
  [ks_t, prop_t, err_t, bdry_t] = gen_levels(30,100,.1,.1,sys,opts);  

  [ks, err, bdry] = merge_levels(ks_t, ks, err_t, err, bdry_t, bdry);
  
  sys = '-l qugrs:1:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1';
  opts.v = 0;
  opts.Mo = 404;
  [ks_t, prop_t, err_t, bdry_t] = gen_levels(200,300,.1,.1,sys,opts);  

  [ks, err, bdry] = merge_levels(ks, ks_t, err, err_t, bdry, bdry_t);
  Mo = size(bdry.rn);
  
  sys = '-l qugrs:1:0.4:0.7 -s oyooo:4:2:1 -u -4 1';
  [ks_t, prop_t, err_t, bdry_t] = gen_levels(5,30,.1,.1,sys,opts);  
  [ks, err, bdry] = merge_levels(ks_t, ks, err_t, err, bdry_t, bdry);
  
  ne = numel(ks);
  save qugrs0-300.mat ne ks err bdry prop opts sys
