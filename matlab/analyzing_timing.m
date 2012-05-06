%% read data
t = dlmread('../data/timing_data_parsed.dat');
ks = t(:,1);
solve_t = t(:,2);
verg_t = t(:,3);
count_t = t(:,4);
total_t = solve_t + verg_t + count_t;

%% plots
figure;
semilogy(ks, solve_t);
hold on;
semilogy(ks, verg_t, 'r');
semilogy(ks, count_t, 'k');
l=legend('solve', 'evaluate', 'count');
set(l, 'FontSize', 20);
set(gca, 'FontSize', 20);
xlabel('k');
ylabel('time (s)');

%% timing predictions
figure;
plot(ks, total_t, '.');
coeffs = polyfit(ks, total_t, 2);
fitted = coeffs(1)*ks.^2 + coeffs(2)*ks + coeffs(3);
hold on;
plot(ks, fitted, 'r');

kmin = 1000;
kmax = 1200;
krange = kmin:kmax;
hours = sum(coeffs(1)*krange.^2 + coeffs(2)*krange + coeffs(3)) / 3600