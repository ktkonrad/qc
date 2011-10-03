% return a matrix A where
% A[i,j] = J_{j-1}(k*r_i) * exp(1i*j*theta_i)
% 1 <= i <= m (where m == length(r))
% 0 <= j-1 <= N

function [A] = bessel_matrix(k, r, theta, M)
    % here k is the norm of the vector k
    % r and theta are column vectors
    A = bessel([0 1:M 1:M], k*r) .* [ones(length(theta),1) sin(theta*(1:M)) cos(theta*(1:M))];
end