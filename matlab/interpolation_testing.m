addpath('../vergini/');
a = load_sta('../data/rpw');
b = reshape(a, 20, 20);
figure; imagesc(b);
figure; imagesc(b>0);
lo_res = dlmread('../data/rpw_averaged.dat');
hi_res = dlmread('../data/rpw_hi.dat');
assert(all(all(abs(lo_res - b) < 1e-6))); % make sure the files are the same rpw

%%
upsample = 20;

cmd = '../c/count -f ../data/rpw_averaged.sta_bin -d 0.01 -M 9 -u 20 -k 30';
system(cmd);
c = load('../data/counted.dat')';
c(c>100000) = 0;
figure; imagesc(c);
figure; imagesc(hi_res>0);



%%
sub = hi_res(260:280, 160:180);

interp = fliplr(read_dumped('../data/interpolated_12_7.dat', upsample+1, upsample+1));
figure; imagesc(interp);
figure; imagesc(sub);