% function [nd] = count_nodal_domains(grid)
%
% count number of nodal domains over a sampled grid
%
% Kyle Konrad 3/28/2011

function [nd] = count_nodal_domains(grid)
   
    [ny, nx] = size(grid);

    counted = zeros(ny, nx); % keep track of whether a point has been counted in a nodal domain

    nd = 0; % number of nodal domains
    i = 1;
    while i <= nx * ny
        i = find(~counted, 1);
    
        if isempty(i)
            break
        end
    
        nd = nd + 1;
        findDomain(i);
    end

    %figure;
    %imagesc(counted);
    
    % search outward from point k in grid to boundaries of nodal domain
    % mark all points in counted
    function [] = findDomain(j)
        if j < 0 || j > nx * ny
            return
        end
        
        current_sign = sign(grid(j));
        to_check = zeros(nx*ny,1);
        to_check(1) = j;
        last = 1;
        p = 1;
        
%         if current_sign == 0
%             nd = nd - 1;
%             return
%         end
        
        while p <= last;
            k = to_check(p);
            counted(k) = nd;
            
            %left
            if k > ny && sign(grid(k-ny)) == current_sign && ~counted(k - ny)
                last = last + 1;
                to_check(last) = k-ny;
            end

            %above
            if mod(k, ny) ~= 1 && sign(grid(k-1)) == current_sign && ~counted(k - 1)
                last = last + 1;
                to_check(last) = k-1;
            end

            %right
            if k <= (nx-1)*ny && sign(grid(k+ny)) == current_sign && ~counted(k + ny)
                last = last + 1;
                to_check(last) = k+ny;
            end

            %below
            if mod(k, ny) > 0 && sign(grid(k+1)) == current_sign && ~counted(k + 1)
                last = last + 1;
                to_check(last) = k+1;
            end
            
            p = p + 1;
        end
        
    end


end