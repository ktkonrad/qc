function [error_norms, conds] = plane_wave_interp_c(k, dx, Ms, plot)

    n = 4.5;
    upsample = 10;
    
    if nargin < 4
        plot = 0;
    end

    cs = normrnd(0,1,[2, 10]);
    theta_ks = 2*pi*rand(1,10);
    k_vecs = [k * cos(theta_ks) ; k * sin(theta_ks)];

    xmin = -(floor(n)-1)/2*dx;
    xmax = (floor(n)-1)/2*dx;

    points = stencil() *dx;

    x2 = (-0.5:1/upsample:0.5)*dx;
    [xs2, ys2] = meshgrid(x2, x2);
    points2 = [xs2(:) ys2(:)];

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
        
        cmd = sprintf('./c/interp %0.8f %d', k*dx, M);
        disp(cmd);
        [status, result] = system(cmd);
        if status ~= 0
            fprintf('error: interp exited with status %d\n\t%s\n', status, result);
            continue;
        end  
        interp = read_dumped('interp.dat', 121, 24);
        interp_mat = interp_matrix(k, points, points2, M);

        norm(interp-interp_mat)
        
        %dlmwrite(['interp_M=' int2str(M) '_dx=' num2str(dx) '.dat'], interp);
        interpolated = interp * f_vals;
        conds(i) = cond(interp);
        
        errors = reshape(interpolated - f_vals2, sqrt(length(f_vals2)), sqrt(length(f_vals2)));
        error_norms(i) = max(max(abs(errors)));

        if plot
            subplot(2,3,M-min(Ms)+1);
            imagesc([-0.5,0.5], [-0.5,0.5], log10(abs(errors)));
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