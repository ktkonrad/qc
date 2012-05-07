%% read data
name = 'qugrs_1000_to_1200_sizes';
sizes = dlmread(['../results/' name '.txt']);
%% scaling for percolation
k = 115;
n = 2*sqrt(2)*k/(pi);
dx = 1/(n/4);
area = 1;
j1 = 2.4048; % bessel function zero
s_min = pi*(j1/k)^2; % smallest possible nodal domain
sizes = sizes*dx^2/s_min; % scale them

%% general
tau = 187/91;
N=10000; % number of buckets
ds = (max(sizes)-min(sizes))/N;
s = linspace(min(sizes), max(sizes), N);
freqs = hist(sizes,s);
s = s(freqs > 0);
freqs = freqs(freqs > 0)/(numel(sizes)*ds);
n = numel(s);

figure;

% histogram
s_plot = reshape([s-ds/2;s+ds/2],1,n*2);
freqs_plot = reshape([freqs;freqs],1,n*2);
loglog(s_plot,freqs_plot, 'k', 'LineWidth', 2)
hold on;

% theoretical
f = s.^-tau;
loglog(s,f, 'LineWidth', 3);

% best fit
fit_range = [10^1 10^3];
fit_idx = s > fit_range(1) & s < fit_range(2);
%p = polyfit(log(s),log(freqs),1);
p = lsqlin(log(s(fit_idx))', log(freqs(fit_idx)), [], [], 0, 0); % same as above but force zero intercept
fit = s.^p(1);
loglog(s, fit, 'r-.', 'LineWidth', 3);

fontsize = 20;
set(gca, 'FontSize', fontsize);
xlabel('s/s_{min}');
ylabel('n');

print('-deps2c', ['../documents/thesis/figs/results/' name '.eps']);