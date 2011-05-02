% statistics on eigenproblem accuracy

eta = 3;
kD = 3;
sys = sprintf('-l qugrs:1:0.2:0.4 -b 10 -s oyooo:%g:%g:1', eta, kD);

[k, t] = efunc_window(sys, 20, 0.2);
