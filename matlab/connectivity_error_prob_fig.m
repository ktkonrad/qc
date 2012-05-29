alpha = logspace(-3,.3,50);
p = arrayfun(@(x) ambiguity_prob(x), alpha);

loglog(alpha, p);
fontsize = 20;
set(gca, 'FontSize', fontsize);
xlabel('\alpha');
ylabel('p(sign change)');
print('-deps2c', '../documents/thesis/figs/interpolation/analytic_connectivity_error_prob.eps');