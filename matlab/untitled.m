alpha = logspace(-3,.3,50);
p = arrayfun(@(x) ambiguity_prob(x), alpha);

loglog(alpha, p);
fontsize = 20;
set(gca, 'FontSize', fontsize);
xlabel('\alpha');
ylabel('p(error)');