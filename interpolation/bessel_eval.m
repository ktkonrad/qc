function [vals] = bessel_eval(coeffs, k, r, theta)
    M = (length(coeffs) - 1) / 2;
    vals = bessel_matrix(k, r, theta, M) * coeffs;
end