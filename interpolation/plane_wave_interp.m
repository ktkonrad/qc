function [error_norms] = plane_wave_interp(k, dx, n, upsample, Ms, cs, theta_ks, plot)
% k is wavenumber
% dx is grid spacing
% sample a grid of n^2 points
% upsample interpolated fn by upsample
% use Bessel fns up to degree M for M in Ms
% cs: 2xp matrix of coeffs for plane waves
% theta_ks: 1xp matrix of directions of plane waves

    if nargin < 5
        Ms = 1:20;
    end
    if nargin < 6
        cs = normrnd(0,1,[2, 10]);
    end
    if nargin < 7
        theta_ks = 2*pi*rand(1,10);
    end
    if nargin < 8
        plot = 0;
    end

    xmin = -(n-1)/2*dx;
    xmax = (n-1)/2*dx;

    x = xmin:dx:xmax;
    xs = meshgrid(x);
    ys = flipud(meshgrid(x)');
    points = [xs(:) ys(:)];
    [theta, r] = cart2pol(points(:,1), points(:,2));

    x2 = xmin:dx/upsample:xmax;
    xs2 = meshgrid(x2);
    ys2 = flipud(meshgrid(x2)');
    points2 = [xs2(:) ys2(:)];

    k_vecs = [k * cos(theta_ks) ; k * sin(theta_ks)];

    f_vals = rpw(cs, k_vecs, points);
    f_vals2 = rpw(cs, k_vecs, points2);

    error_norms = zeros(length(Ms), 1);
    conds = zeros(length(Ms), 1);

    if plot
        figure;
    end
    i = 0;
    for M=Ms
        i = i + 1;
        [interpolated, conds(i)] = bessel_interp(norm(k), r, theta, f_vals, M, points2);
        errors = reshape(interpolated - f_vals2, sqrt(length(f_vals2)), sqrt(length(f_vals2)));
        error_norms(i) = max(max(abs(errors)));

        if plot
            subplot(5,10,M);
            imagesc([xmin,xmax], [xmin,xmax], log10(abs(errors)));
            caxis([-16, 0]);
            line(repmat(x,2,1), repmat([xmin;xmax],1,length(x)), 'color', 'k');
            line(repmat([xmin;xmax],1,length(x)), repmat(x,2,1), 'color', 'k');
        end
    end

    
    if plot
        figure;
        subplot(2,1,1);
        semilogy(1:length(Ms), error_norms, '-');
        xlabel('M');
        ylabel('||errors||_{\infty}');

        subplot(2,1,2);
        semilogy(1:length(Ms), conds);
        xlabel('M');
        ylabel('cond(A)');
    end
end