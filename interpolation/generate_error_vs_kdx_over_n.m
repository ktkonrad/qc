dxs = 10.^(-1:-1:-8);
k = 200;
ns = 2:7;
upsample = 10;
M = 7;

error_norms = zeros(length(dxs), length(ns));

i = 0;

for n=ns
     i = i+1;
     j = 0;
    for dx = dxs
        j = j+1;
        error_norms(j,i) = plane_wave_interp(k, dx, n, upsample, M);
    end
end

figure;

%%% plot error_norms vs k*dx %%%
loglog(k*dxs, error_norms);
xlabel('k*dx');
ylabel('||errors||_{\infty}');
legend(arrayfun(@(x) strcat('n = ', num2str(x)), ns, 'UniformOutput', false));
