load('rpw_interp2.mat');

true_counts=repmat(no_interp_counts(numel(alpha),:), numel(alpha),1);
no_interp_errors = no_interp_counts-true_counts;
no_interp_rel_errors = abs(no_interp_errors)./true_counts;

figure; loglog(alpha,mean(no_interp_rel_errors(:,1:35),2));
fontsize = 20;
set(gca, 'FontSize', fontsize);
xlabel('alpha', 'FontSize', fontsize);
ylabel('mean relative error', 'FontSize', fontsize);
print('-deps2c', '../documents/thesis/figs/interpolation/rpw_mean_relative_count_error.eps');