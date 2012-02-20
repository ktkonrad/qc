k = 200;
dx = 0.001;
M = 9;
rtyp = 3*sqrt(0.5)*dx;


cmd = sprintf('./c/interp %0.8f %d', k*dx, M);
disp(cmd);
[status, result] = system(cmd);
if status ~= 0
    fprintf('error: interp exited with status %d\n\t%s\n', status, result);
    return;
end

A = read_dumped('A.dat', 24, 19);
V = read_dumped('V.dat', 19, 19);
U = read_dumped('U.dat', 24, 19);
S = diag(read_dumped('S.dat', 19, 1));
interp = read_dumped('interp.dat', 121, 24);
B = read_dumped('B.dat', 121, 19);
A_plus = read_dumped('A_plus.dat', 19, 24);

points = stencil();
x2 = (-0.5:1/upsample:0.5)*dx;
[xs2, ys2] = meshgrid(x2, x2);
points2 = [xs2(:) ys2(:)];

A_mat = bessel_matrix2(k, points*dx, M, rtyp);
[U_mat, S_mat, V_mat] = svd(A_mat, 0);
B_mat = bessel_matrix2(k, points2, M, rtyp);
A_plus_mat = pinv(A_mat);
interp_mat = B_mat*A_plus_mat;

fprintf('A: %g\n', norm(A - A_mat));
fprintf('B: %g\n', norm(B - B_mat));
fprintf('U: %g\n', norm(U - U_mat));
fprintf('S: %g\n', norm(S - S_mat));
fprintf('V: %g\n', norm(V - V_mat));
fprintf('A+: %g\n', norm(A_plus - A_plus_mat));
fprintf('interp: %g\n', norm(interp - interp_mat));

