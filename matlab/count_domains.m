k0s = 200:50:400;
counts = zeros(30 * length(ks), 1);
ks = zeros(30 * length(ks), 1);
times = zeros(length(ks), 1);

for i = 1:length(ks);
    dx = .1 / ks(i);
    cmd = sprintf('../vergini/verg -l qust:2 -b 10 -s vepwoo:1.2:40:1.5 -k %g -V 0.05 -f %g:1 -m > /dev/null 2> /dev/null', ks(i), dx);

    tic;
    [status, output] = system(cmd);
    vergtime = toc;

    if (status ~= 0)
        disp(['error: ' output]);
        continue;
    end

    cmd = sprintf('../c/count -f t.sta_bin -l qust:2 -d %g -k %g', dx, ks(i));
    tic;
    [status, output] = system(cmd);%'../c/count -f t.sta_bin -m t.mask.sta_bin');
    counttime = toc;

    times(i) = vergtime/counttime;

    if (status ~= 0)
        disp(['error: ' output]);
        continue;
    end
    count = regexp(output, '\d+\s+\d+', 'match');
    counts(i) = str2double(count{1});
end

figure;
plot(log10(dxs'*ks), counts');
xlabel('log_{10}(k*dx)');
ylabel('nodal domain count');
legend('k = 200','k = 250','k = 300','k = 350','k = 400');

figure;
plot(log10(dxs), times');
xlabel('log_{10}(dx)');
ylabel('runtime ratio: verg / count');
legend('k = 200','k = 250','k = 300','k = 350','k = 400');
