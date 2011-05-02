% see if matlab lin alg fast enough.

% time bd q pwoo 500 qust 1000 2.0 v 300 1e-14 20 4 0.1 0 1000
% Takes 25 sec.

% approx equiv timing in matlab:
%       40 sec (if everything vectorized, including fill of A)
%       290 sec (if N decisions are made per fill of A, using matlabspeed_fillA)
% Almost all time spend in writing out bdry functions.

% Conclusion: Matlab is about as fast as C/LAPACK if vectorize everything, but
% if need decisions, eg arbitrary different basis functions in same
% basis set, it will be too slow. (This could be got around be
% grouping the basis functions by type - messy?).

% barnett 8/22/03

tic;

N=500;
M=1000;

flag = rand(N,1);
a = kron((1:M)',1:N);
%A = sin(a) + cos(a);
A = matlabspeed_fillA(a, flag);
disp('A filled');
F = A'*A;
B = sin(a/10) + cos(a/5) + sin(a) + cos(a/4);
G = B'*B;
disp('F,G built');
size(F)

[V,d] = eig(F,G);
disp('eig done');

% note slowdown if vectors needed:
% tic;[V,d] = eig(F);toc          3.1572 sec
% tic;eig(F);toc                  0.7340 sec


p = zeros(1,M);
for i=1:20
  a = kron((1:M)',1:N);
  A = matlabspeed_fillA(a, flag);
  p = A*rand(N,1);
  B = sin(a/10) + cos(a/5) + sin(a) + cos(a/4);
  p = B*rand(N,1);
  disp(i);
end
disp('bdry values done');
toc
