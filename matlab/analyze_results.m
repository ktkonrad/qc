%% read stats
filename = '../condor/counts.txt';
stats = dlmread(filename);
ks = stats(:,1);
counts = stats(:,2);

% constants
alpha = 0.5;
fontsize = 20;
area = 0.614; % qust 0.4:0.7

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
ws = 200; %window size
means = [];
k_windows = [];
for i=1:ws:length(counts)
    idx = i:min(i+ws,length(counts));
    k = ks(min(i+ws/2,length(ks)));
    k_windows = [k_windows k];
    means = [means mean(scaled_counts(idx))];
end

mean_error = abs(mean(means) - mean_predicted) / mean_predicted;

figure;
plot(ks, scaled_counts, '.');
hold on;
plot([min(ks), max(ks)], [mean_predicted, mean_predicted], 'r-', 'LineWidth', 3);
plot(k_windows, means, 'k-', 'LineWidth', 3);
xlabel('k', 'FontSize', fontsize);
ylabel('\nu(k)/N(k)', 'FontSize', fontsize);
legend('data', 'predicted mean', 'measured mean');
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/results/mean.eps');

%% compute variance
variance_predicted = 18/pi^2 + 4*sqrt(3)/pi - 25/(2*pi);

ws = 200; %window size
vars = [];
k_windows = [];
for i=1:ws:length(counts)
    idx = i:min(i+ws,length(counts));
    k = ks(min(i+ws/2,length(ks)));
    k_windows = [k_windows k];
    vars= [vars 4*pi*var(counts(idx))/(area*k^2)];
end

figure;
plot([min(ks), max(ks)], [variance_predicted, variance_predicted], 'r-', 'LineWidth', 3);
hold on;
plot(k_windows, vars, 'k-', 'LineWidth', 3);
xlabel('k', 'FontSize', fontsize);
ylabel('\sigma^{2}(k)/N(k)', 'FontSize', fontsize);
legend('predicted mean', 'measured mean');
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/results/variance.eps');

%% check interp counts
interp_counts = stats(:,4);

% interpolations per domain
%figure; plot(interp_counts./counts);

% interpolations per pixel
dxs = alpha ./ ks;
pixels = area ./ dxs.^2;
figure; plot(ks, interp_counts./pixels, '.');
hold on;
m = mean(interp_counts./pixels);
plot([min(ks), max(ks)], [m, m], 'r-', 'LineWidth', 3);
xlabel('k', 'FontSize', fontsize);
ylabel('interpolations per pixel', 'FontSize', fontsize);
legend('data', 'mean'); 
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/results/interp_counts.eps');
