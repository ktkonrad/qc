addpath('../vergini');

%% specify ppw
ppw = 2.^(2:8);% 100 150 200 300];
ppw = [10 200];
n0 = 50;
n = n0 .* ppw ./ ppw(1);
L = 2*pi./ppw .* n; % width of window
dx = L./n;
k = 2*pi*n./(L.*ppw); % should be all ones
alpha = k.*dx;

%% run 100 times
N = 1;
no_interp_counts= zeros(numel(ppw), N);
no_interp_errors = zeros(numel(ppw), N);
for j=1:N
[fs, xs] = rpw2dsample_multi(n0, ppw, 7);

for i=1:numel(fs)
    figure; imagesc(xs{i}, xs{i}, fs{i}>0);
    %sta = zeros(1,n(i),n(i));
    %sta(1,:,:) = fs{i};
    %save_sta(sta, k(i), sprintf('../data/rpw_%.8f', dx(i)));
end
 figure; imagesc(xs{2}(1:20:end,1:20:end), xs{2}(1:20:end,1:20:end), fs{2}(1:20:end,1:20:end)>0);
%%
 system('../scripts/count_rpws.py ../data');

no_interp_stats = dlmread('rpw_counts_no_interp.txt');
alpha = no_interp_stats(:,1) .* no_interp_stats(:,2);
no_interp_counts(:,j) = no_interp_stats(:,3);

true_count = no_interp_counts(numel(ppw));
no_interp_errors(:,j) = abs(no_interp_counts(:,j) - true_count) ./ true_count;

save('rpw_interp2.mat', 'alpha', 'no_interp_counts', 'no_interp_errors');
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
