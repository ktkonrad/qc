ws = 500; %window size
vars = [];
for i=1:ws:length(counts)
    vars= [vars std(counts(i:min(i+ws,length(counts))))^2/(vol/(4*pi)*ks(min(round(i+ws/2),length(ks)))^2)];
end

figure;
plot(1:ws:length(counts),vars);