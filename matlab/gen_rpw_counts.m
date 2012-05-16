k = 100:100:1100; % don't go above 1100 or you have to swap to disk (w/ 6GB ram)
alpha = 0.5;

N = 100; % number of repetitions

n = numel(k);
%%
counts = zeros(n,N);
scaled_counts = zeros(n,N);
%%
for i=11
    for j=1:N
        [counts(i,j), scaled_counts(i,j)] = count_rpw(alpha, k(i));
    end
end

save('rpw_counts.mat', 'counts', 'scaled_counts')
%% reshape for plotting with analyze_results script

t=repmat(k,N,1);
t(:);
ks = t(:);
c=counts';
counts = c(:);
s=scaled_counts';
scaled_counts = s(:);
