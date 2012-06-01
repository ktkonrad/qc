alpha = logspace(log10(2.3e-2),log10(2.4),100);
p = arrayfun(@(x) ambiguity_prob(x), alpha);

figure;
loglog(alpha, p);
fontsize = 20;
set(gca, 'FontSize', fontsize);
xlabel('\alpha');
ylabel('p(ambiguity)');
print('-deps2c', '../documents/thesis/figs/interpolation/theoretical_ambiguity_prob.eps');
