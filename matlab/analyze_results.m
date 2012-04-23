%% read stats
filename = '../condor/qugrs_600_to_1200.out';
stats = dlmread(filename);
ks = stats(:,1);
counts = stats(:,2);

%% compute mean
mean_predicted = (3*sqrt(3) - 5)/pi;

% qugrs shape
a = 1;
t1=.4;
t2=.7;
R1 = a/sin(t1);
R2=1/sin(t2);
area = a - R1^2*(2*t1 - sin(2*t1))/4 - R2^2*(2*t2 - sin(2*t2))/4;

scaled_counts = 4*pi*counts./(area*ks.^2);
mean_measured = mean(scaled_counts);

mean_error = abs(mean_measured - mean_predicted) / mean_predicted;

figure;
plot(ks, scaled_counts, '.');
hold on;
plot(ks, mean_predicted, 'r-');
plot(ks, mean_measured, 'k-');

%% compute variance
variance_predicted = 18/pi^2 + 4*sqrt(3)/pi - 25/(2*pi);

ws = 200; %window size
vars = [];
for i=1:ws:length(counts)
    vars= [vars std(counts(i:min(i+ws,length(counts))))^2/ks(min(i+ws/2,length(ks)))^2];
end

figure;
plot(vars);
hold on;
plot(1:length(vars), variance_predicted, 'r-');