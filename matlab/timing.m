n = 20; % number of trials
s = 4; % increment to grid size
times = zeros(n,2);

for i=1:n
    a = random_percolation(s*i,s*i);
    
    tic;
    %count_nodal_domains(a);
    times(i,1) = toc;
    
    tic;
    count_nodal_domains_recursive(a);
    times(i,2) = toc;
end

%figure;
%plot(4*(1:n),times(:,1));
figure;
plot(4*(1:n),times(:,2));