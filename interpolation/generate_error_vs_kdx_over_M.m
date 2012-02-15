dxs = 10.^(-1:-1:-5); %[.1, .05, .01, .005, .001, .0005, .0001, .00005, .00001];
k = 200;
n = 4.5;
upsample = 10;
Ms = 6:11;

error_norms = zeros(length(Ms), length(dxs));
conds = zeros(length(Ms), length(dxs));

i = 1;
for dx = dxs
    [error_norms(:,i), conds(:,i)] = plane_wave_interp_c(k, dx, Ms, 0);
    %[error_norms(:,i), conds(:,i)] = plane_wave_interp(k, dx, n, upsample, Ms);
    i = i+1;
end

figure;

%%% plot error_norms vs k*dx at each M %%%
loglog(k*dxs, error_norms);
legend(arrayfun(@(x) strcat('M = ', num2str(x)), Ms, 'UniformOutput', false));
xlabel('k*dx');
ylabel('||errors||_{\infty}');

%%% 3d plot error_norms vs k*dx vs M %%%
figure;surf(log10(k*dxs), Ms, log10(error_norms));xlabel('log(k*dx)');ylabel('M'),zlabel('log(error)');
figure;surf(log10(k*dxs), Ms, log10(conds));xlabel('log(k*dx)');ylabel('M'),zlabel('log(cond)');
