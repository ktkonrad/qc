function [count, scaled_count] = count_rpw(alpha, k)
    addpath('../vergini');
    % rescale to k = 1
    L = k; % length of one side of domain
    k = 1;
    dx = alpha/k;
    ppw = 2*pi/alpha; % sample points per wavelength
    area = L^2;
    n = L/alpha;
    
    count = -1;

    mem = (4*n)^2*8 / 2^30; % rpw2sample uses this much memory (in GB)
    
    if (mem > 5)
        printf('this will use more than 5GB');
        return
    end
    
    sta = zeros(1,n,n);
    sta(1,:,:) = rpw2dsample(n, ppw);
    name = sprintf('../data/rpw_%.8f', dx);
    save_sta(sta, k, name);
    cmd = sprintf('../c/count -f %s.sta_bin -k %f -a %f -M 9 -u 30 -q', name, k, alpha);
    [status, out] = system(cmd);
    if status
        printf('count failed');
        return
    end
    [k,dx,count,small,interp,bdry,edge] = strread(out,'%f%f%d%d%d%d%d','delimiter',',');
    scaled_count = count / (area*k^2/(4*pi));
end