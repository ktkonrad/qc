%% read stats
stats = dlmread('../c/counts.txt');
ks = stats(:,1);
counts = stats(:,2);

%% compute mean
% qugrs shape
a = 1;
t1=.4;
t2=.7;
R1 = a/sin(t1);
R2=1/sin(t2);
area = a - R1^2*(2*t1 - sin(2*t1))/4 - R2^2*(2*t2 - sin(2*t2))/4;

scaled_counts = 4*pi*counts./(area*ks.^2);
mean(scaled_counts)

