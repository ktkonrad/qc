function [f_vals] = plot_bessel(coeffs, k, dx)
xs = meshgrid(-.5:.1:.5);
ys = flipud(meshgrid(-.5:.1:.5)');

[r, theta] = cart2pol(xs(:), ys(:));

f_vals = reshape(real(bessel_eval(coeffs, k, r*dx, theta)), 11, 11);

imagesc(xs(:), ys(:), f_vals)
colorbar;

end

