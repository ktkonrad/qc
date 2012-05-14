k = 100:100:1100; % don't go above 1100 or you have to swap to disk (w/ 6GB ram)
alpha = 0.5;

N = 100; % number of repetitions

n = numel(L);
counts = zeros(n,N);
scaled_counts = zeros(n,N);

for i=1:n
    for j=1:N
        [counts(i,j), scaled_counts(i,j)] = count_rpw(alpha, k(i));
    end
end

figure;
plot(L, mean(scaled_counts));