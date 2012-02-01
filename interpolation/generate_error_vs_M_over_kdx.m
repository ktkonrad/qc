dxs = 10.^(-1:-1:-8); %[.1, .05, .01, .005, .001, .0005, .0001, .00005, .00001];
k = 200;
n = 4;
upsample = 10;
Ms = 1:20;

error_norms = zeros(length(Ms), length(dxs));
conds = zeros(length(Ms), length(dxs));

i = 1;
for dx = dxs
    [error_norms(:,i), conds(:,i)] = plane_wave_interp(k, dx, n, upsample, Ms);
    i = i+1;
end

figure;

%%% plot error_norms vs M at each k*dx %%%
semilogy(error_norms);
legend(arrayfun(@(x) strcat('k*dx = ', num2str(x)), k*dxs, 'UniformOutput', false));
xlabel('M');
ylabel('||errors||_{\infty}');

figure;
semilogy(conds);
legend(arrayfun(@(x) strcat('k*dx = ', num2str(x)), k*dxs, 'UniformOutput', false));
xlabel('M');
ylabel('conds');