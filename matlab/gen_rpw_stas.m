addpath('../vergini');

%% specify alpha
alpha = logspace(-1.5, 0.4, 6);
k = ones(1, numel(alpha));
dx = alpha;
n0 = 10; % number of 1-d sample points for largest dx
n = n0 ./ dx .* max(dx);
L = n .* dx; % should be constant

ppw = round(2*pi./L.*n); % must be integral

% recompute after rounding
n = n0 .* ppw ./ min(ppw);
L = 2*pi./ppw .* n; % width of window
dx = L./n;
k = 2*pi*n./(L.*ppw);  % should all be 1.0

%% specify ppw
ppw = [2:10 12:2:18 20:6:44 50:10:80];
n0 = 20;
n = n0 .* ppw ./ ppw(1);
L = 2*pi./ppw .* n; % width of window
dx = L./n;
k = 2*pi*n./(L.*ppw);

N = 100;
interp_errors = zeros(numel(ppw), N);
no_interp_errors = zeros(numel(ppw), N);
for j=1:N
%% generate rpws and save as sta_bin
[fs, xs] = rpw2dsample_multi(n0, ppw);

for i=1:numel(fs)
    %figure; imagesc(xs{i}, xs{i}, fs{i}>0);
    sta = zeros(1,n(i),n(i));
    sta(1,:,:) = fs{i};
    save_sta(sta, k(i), sprintf('../data/rpw_%.8f', dx(i)));
end

%% count nodal domains
system('../scripts/count_rpws.py ../data');

%% plot results
interp_stats = dlmread('rpw_counts_interp.txt');
no_interp_stats = dlmread('rpw_counts_no_interp.txt');
alpha = interp_stats(:,1) .* interp_stats(:,2);
interp_counts = interp_stats(:,3);
no_interp_counts = no_interp_stats(:,3);

true_count = interp_counts(1);
interp_errors(:,j) = abs(interp_counts - true_count) ./ true_count;
no_interp_errors(:,j) = abs(no_interp_counts - true_count) ./ true_count;

%figure;
%plot(alpha, interp_errors);
%hold on;
%plot(alpha, no_interp_errors, 'r');
end