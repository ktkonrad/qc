n = 100;
x = linspace(0,1,n);
y=x;

k_mag = 3; % magnitude of k
k_th = pi/6; % direction of k

k_vec = k_mag*[cos(k_th) sin(k_th)];

[xs,ys] = meshgrid(x,y);
x_stacked = xs(:);
y_stacked = ys(:);
points = [x_stacked y_stacked];

%% random fourier coeffs
f = zeros(N);
N=100;
for k=1:N
    for l=1:N
        c = rand(4,1);
        f = f + c(1) * cos(2*pi*k*xs) .* cos(2*pi*l*ys) + ...
            c(2) * cos(2*pi*k*xs) .* sin(2*pi*l*ys) + ...
            c(3) * sin(2*pi*k*xs) .* cos(2*pi*l*ys) + ...
            c(4) * sin(2*pi*k*xs) .* sin(2*pi*l*ys);
    end
end


figure;
surf(x,y,f);

figure;
imagesc(f>0);

%%

%pw = real(exp(1i*k_mag*2*pi*(k_vec(1)*xs+k_vec(2)*ys)));
pw = rpw2dsample(n, 12);
figure; surf(pw);

N=100;
c = zeros(N,N,4);
for k=1:N
    for l=1:N
        c(k,l,1) = sum(sum((cos(2*pi*k*xs) .* cos(2*pi*l*ys)) .* pw))*dx^2;
        c(k,l,2) = sum(sum((cos(2*pi*k*xs) .* sin(2*pi*l*ys)) .* pw))*dx^2;
        c(k,l,3) = sum(sum((sin(2*pi*k*xs) .* cos(2*pi*l*ys)) .* pw))*dx^2;
        c(k,l,4) = sum(sum((sin(2*pi*k*xs) .* sin(2*pi*l*ys)) .* pw))*dx^2;
    end
end

figure;
for i=1:4
subplot(2,2,i);
plot(c(:,:,i));
end
