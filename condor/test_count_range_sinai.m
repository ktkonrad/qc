SYSLIBS='/usr/lib:/usr/lib/x86_64-linux-gnu';
LLP=strcat(SYSLIBS,':',getenv('LD_LIBRARY_PATH'));
setenv('LD_LIBRARY_PATH',LLP);

opts.kdx = exp(-.5);
opts.v = 2;
[counts,ks,wtms] = count_range_sinai(300, 302, .1, .1, ...
                   '-l qugrs:1:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1', opts);
save('counts.mat');
