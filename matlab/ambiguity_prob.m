function [p] = ambiguity_prob(alpha)
% numerically integrate a 3d gaussian to compute probability of an
% ambiguity. Uses a linear transformation from A = {a_0, a_1, a_2} to
% X = {x_0, x_1, x_2}

mu_A = zeros(3,1); % mean in A basis
S_A = diag([1;2;2]); % covariance matrix in A basis

beta = alpha/2;
beta_0 = besselj(0, beta);
beta_1 = besselj(1, beta);
beta_2 = besselj(2, beta);

% transform [a_0;a_1;a_2] to [x_0;x_1;x_2] using
%
% x_0 = a_0
% x_1 = a_1 + beta_0/beta_1 * a_0
% x_2 = a_2 + (a_0*beta_0 + a_1*beta_1)/beta_2
%
% this transforms integration region to [-\infty,0]^3

% J(i,j) = partial x_i / partial a_j
J = [-1 0 0 ; beta_0/beta_1 1 0 ; beta_0/beta_2 beta_1/beta_2 1]; % jacobian
S_X = inv(J*S_A*J'); % covariance in X basis. worked this out by hand. probably standard for coordinate transformations
mu_X = J*mu_A; % mean in X basis. still just zero

p = 2*mvncdf([0;0;0], mu_X, S_X);

end