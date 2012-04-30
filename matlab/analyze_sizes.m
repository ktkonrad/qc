%% percolation
sizes = dlmread('../c/perc_sizes.dat');
k = 115;
n = 2*sqrt(2)*k/(pi);
dx = 1/(n/4);
area = 1;
sizes = sizes*dx^2;


%% eigenfunctions
sizes = dlmread('../c/sizes_all.dat');
dx = .00099;
k = 1001;
area = .614; % qust 0.4:0.7
sizes = sizes*dx^2;

%% general
tau = 187/91;
j1 = 2.4048; % bessel function zero
s_min = pi*(j1/k)^2; % smallest possible nodal domain
N=200;
ds = (max(sizes)-min(sizes))/N;
s = min(sizes):ds:max(sizes); % = linspace(min(sizes), max(sizes), N)
f = s.^-tau;
N_bar = area*k^2/(4*pi);
n_bar = .0624*N_bar;

% scale f
%int_f = sum(f) * ds;
%f = f * (length(sizes) / int_f);

figure; loglog(s/s_min,f/n_bar);
hold on;
[freqs, bin_centers] = hist(sizes/s_min,100);
loglog(bin_centers,freqs/n_bar, 'k')