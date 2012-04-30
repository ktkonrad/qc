%% percolation
sizes = dlmread('../c/perc_sizes.dat');
k = 115;
n = 2*sqrt(2)*k/(pi);
dx = 1/(n/4);
area = 1;
j1 = 2.4048; % bessel function zero
s_min = pi*(j1/k)^2; % smallest possible nodal domain
sizes = sizes*dx^2/s_min; % scale them


%% eigenfunctions
sizes = dlmread('../c/all_sizes.dat'); % they're already scaled
sizes = sizes(:);
%%
sizes = sizes(sizes > 1);

%% general
tau = 187/91;
N=100000;
ds = (max(sizes)-min(sizes))/N;
s = linspace(min(sizes), max(sizes), N);
freqs = hist(sizes,s);
figure;

s_plot = reshape([s-ds/2;s+ds/2],1,N*2);
freqs_plot = reshape([freqs;freqs],1,N*2);
loglog(s_plot,freqs_plot/length(sizes)/ds, 'k', 'LineWidth', 2)

hold on;
f = s.^-tau;
loglog(s,f, 'LineWidth', 3);

fontsize = 20;
set(gca, 'FontSize', fontsize);
xlabel('s/s_{min}');
ylabel('n');
print('-deps2c', '../documents/thesis/figs/results/sizes.eps');