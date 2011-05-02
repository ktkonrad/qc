% function [nd] = count_nodal_domains(grid)
%
% count number of nodal domains over a sampled grid
%
% Kyle Konrad 3/28/2011

function [nd] = count_nodal_domains_recursive(grid)
   
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

    figure;
    imagesc(counted);
    
    % search outward from point j in grid to boundaries of nodal domain
    % mark all points in counted
    function [] = findDomain(j)
        if j < 0 || j > nx * ny
            return
        end
        
        counted(j) = nd;
        current_sign = sign(grid(j));
        
%         if current_sign == 0
%             nd = nd - 1;
%             return
%         end
        
        
        %left
        if j > ny && sign(grid(j-ny)) == current_sign && ~counted(j - ny)
            findDomain(j - ny);
        end
            
        %above
        if mod(j, ny) ~= 1 && sign(grid(j-1)) == current_sign && ~counted(j - 1)
            findDomain(j - 1)
        end
        
        %right
        if j <= (nx-1)*ny && sign(grid(j+ny)) == current_sign && ~counted(j + ny)
            findDomain(j + ny);
        end
        
        %below
        if mod(j, ny) > 0 && sign(grid(j+1)) == current_sign && ~counted(j + 1)
            findDomain(j + 1)
        end
        
    end


end