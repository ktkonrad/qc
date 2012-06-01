%% read data
t = dlmread('../results/old/qust_700_to_900_timing.txt');
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
print('-deps2c', '../documents/thesis/figs/timing/qust_700_to_900_partial_timing.txt');

%% timing predictions
figure;
plot(ks, total_t, '.');
coeffs = polyfit(ks, total_t, 2);
fitted = coeffs(1)*ks.^2 + coeffs(2)*ks + coeffs(3);
hold on;
plot(ks, fitted, 'r');

kmin = 700;
kmax = 900;
krange = kmin:kmax;
hours = sum(coeffs(1)*krange.^2 + coeffs(2)*krange + coeffs(3)) / 3600