function M = interp_matrix(k, points, points2, M)
rtyp = 3*sqrt(max(sum(points2.^2,2)));  % furthest dist of fine grid
    A = bessel_matrix2(k, points, M, rtyp);
    B = bessel_matrix2(k, points2, M, rtyp);
    M = B*pinv(A);
end