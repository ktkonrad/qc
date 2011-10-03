% compute coefficients of Bessel expansion

function [coeffs, condition] = bessel_coeffs(f_vals, bessel_vals)
    preconditioner = max(bessel_vals);
    preconditioned = bessel_vals ./ repmat(preconditioner, size(bessel_vals, 1), 1);
    coeffs = lscov(preconditioned, f_vals) ./ preconditioner';
    condition = cond(preconditioned);
end