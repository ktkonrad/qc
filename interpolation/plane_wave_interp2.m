function [error_norms] = plane_wave_interp2(k, ingrids, outgrid, M, cs, theta_ks, plot)
% k is wavenumber
% ingrids is cell array of grids, each N_j x 2  1 <= j <= l
% use Bessel fns up to degree M
% cs: 2xp matrix of coeffs for plane waves
% theta_ks: 1xp matrix of directions of plane waves
% error_norms: lx1

    if nargin < 4
        M = 8;
    end
    if nargin < 5
        cs = normrnd(0,1,[2, 10]);
    end
    if nargin < 6
        theta_ks = 2*pi*rand(1,10);
    end
    if nargin < 7
        plot = 0;
    end

    error_norms = zeros(length(ingrids), 1);
    k_vecs = [k * cos(theta_ks) ; k * sin(theta_ks)];
    f_vals2 = rpw(cs, k_vecs, outgrid);

    for j=1:length(ingrids)
        f_vals = rpw(cs, k_vecs, ingrids{j});
        [theta, r] = cart2pol(ingrids{j}(:,1), ingrids{j}(:,2));
        interpolated = bessel_interp(k, r, theta, f_vals, M, outgrid);
        %figure;imagesc(reshape(interpolated, sqrt(length(interpolated)), sqrt(length(interpolated))))
        errors = reshape(interpolated - f_vals2, sqrt(length(f_vals2)), sqrt(length(f_vals2)));
        error_norms(j) = max(max(abs(errors)));
  
        if plot
            xmin = min(ingrids{j}(:,1));
            xmax = max(ingrids{j}(:,2));
            figure;
            imagesc([xmin,xmax], [xmin,xmax], log10(abs(errors)));
            hold on;
            caxis([-16, 0]);
            scatter(ingrids{j}(:,1), ingrids{j}(:,2), 500, '+', 'k')
        end
    end
end