ks = 400;%200:50:400;
dxs = [.0001];%[.1,.08,.05,.02,.01,.008,.005,.002,.001,.0005,.0002,.0001];
counts = zeros(length(ks), length(dxs));

for i = 1:length(ks);
    for j = 1:length(dxs);
        cmd = sprintf('../vergini/verg -l qust:2 -b 10 -s vepwoo:1.2:40:1.5 -k %g -V 0.1 -f %g:1 -m > /dev/null 2> /dev/null', ks(i), dxs(j));
        
        tic;
        [status, output] = system(cmd);
        vergtime = toc;
        
        if (status ~= 0)
            disp(['error: ' output]);
        end
        
        cmd = sprintf('../c/count -f t.sta_bin -l qust:2 -d %g -k %g', dxs(j), ks(i));
        tic;
        [status, output] = system(cmd);%'../c/count -f t.sta_bin -m t.mask.sta_bin');
        counttime = toc;
        
        disp(vergtime/counttime);
        
        if (status ~= 0)
            disp(['error: ' output]);
        end
        count = regexp(output, '\d+', 'match');
        counts(i,j) = str2double(count{1});
    end
end

plot(log10(dxs'*ks), counts');
xlabel('log_{10}(k*dx)')
ylabel('nodal domain count')
legend('k = 200','k = 250','k = 300','k = 350','k = 400')