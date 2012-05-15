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
ppw = [2:10 12:2:18 20:6:44 50:10:90 100:20:160];
%ppw = [2 200 260];
n0 = 100;
n = n0 .* ppw ./ ppw(1);
L = 2*pi./ppw .* n; % width of window
dx = L./n;
k = 2*pi*n./(L.*ppw); % should be all ones
alpha = k.*dx;

%% run 100 times
N = 100;
N = 68;
interp_counts= zeros(numel(ppw), N);
no_interp_counts= zeros(numel(ppw), N);
interp_errors = zeros(numel(ppw), N);
no_interp_errors = zeros(numel(ppw), N);
for j=1:N
[fs, xs] = rpw2dsample_multi(n0, ppw);

for i=1:numel(fs)
    %figure; imagesc(xs{i}, xs{i}, fs{i}>0);
    sta = zeros(1,n(i),n(i));
    sta(1,:,:) = fs{i};
    save_sta(sta, k(i), sprintf('../data/rpw_%.8f', dx(i)));
end

system('../scripts/count_rpws.py ../data');

interp_stats = dlmread('rpw_counts_interp.txt');
no_interp_stats = dlmread('rpw_counts_no_interp.txt');
alpha = interp_stats(:,1) .* interp_stats(:,2);
interp_counts(:,j) = interp_stats(:,3);
no_interp_counts(:,j) = no_interp_stats(:,3);

true_count = interp_counts(numel(ppw));
interp_errors(:,j) = abs(interp_counts(:,j) - true_count) ./ true_count;
no_interp_errors(:,j) = abs(no_interp_counts(:,j) - true_count) ./ true_count;

save('rpw_interp2.mat', 'interp_counts', 'no_interp_counts', 'interp_errors', 'no_interp_errors');
end

%% plot errors
fontsize = 20;
figure;

loglog(alpha, mean(no_interp_errors,2), 'r');
hold on;
loglog(alpha, mean(interp_errors,2));
xlabel('alpha', 'FontSize', fontsize);
ylabel('error', 'FontSize', fontsize);
legend('no interpolation', 'interpolation');
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/interpolation/rpw_errors.eps');

%% plot interp frequency
figure;loglog(alpha, mean(interp_counts,2)./(n.^2)')
xlabel('alpha', 'FontSize', fontsize);
ylabel('interpolation frequency', 'FontSize', fontsize);
set(gca, 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/interpolation/rpw_interp_frequencies.eps');
