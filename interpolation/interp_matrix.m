% return a matrix BA^+
% where A^+ is the pseudoinverse of A
% where A is matrix of J bessel coeffs on course grid
% and B is matrix of J bessel coeffs on fine grid

function [C] = interp_matrix(k, points, points2, M)
    % here k is the norm of the vector k
    % r and theta are column vectors
    
    A = bessel_matrix2(k, points, M);
    %preconditioner = diag(1./max(A));
    %preconditioned = A * preconditioner;
    
    %A_cond = cond(A)
    %preconditioned_cond = cond(preconditioned)
    
    B = bessel_matrix2(k, points2, M);
    C = B*(pinv(A));
end
