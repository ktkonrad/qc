ws = 200; %window size
vars = [];
for i=1:ws:length(counts)
    vars= [vars std(counts(i:min(i+ws,length(counts))))^2/ks(min(i+ws/2,length(ks)))^2];
end

figure;
plot(vars);