% bessel interpolate f sampled at points with value f_vals
% use bessel functions up to order M
% return introplated function evalutated at eval_points

function [interpolated, condition] = bessel_interp2(k, points, f_vals, M, eval_points)
    % k scalar
    % points px2
    % f_vals column vectors
    % M order of highest bessel func
    % eval_points 2xN (x,y) to evaluate interpolation function at
    
    [theta, r] = cart2pol(points(:,1), points(:,2));
    bessel_vals = bessel_matrix(k, r, theta, M);
    [coeffs, condition] = bessel_coeffs(f_vals, bessel_vals);
    [eval_theta, eval_r] = cart2pol(eval_points(:,1), eval_points(:,2));
    interpolated = bessel_eval(coeffs, k, eval_r, eval_theta);    
end