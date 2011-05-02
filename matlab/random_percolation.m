% function [grid] = random_percolation(m,n)
%
% generate m x n grid using random percolation model
% (m and n will be rounded down to a multiple of 4)
%
% Kyle Konrad 3/29/2011

function [grid] = random_percolation(m,n)

    cell = [1,1,-1,-1;1,1,-1,-1;-1,-1,1,1;-1,-1,1,1];
    grid = repmat(cell, m/4, n/4);

    for i=1:m/4*2
        for j=1:n/4*2-1
            if rand(1)<.5;
                grid(2*i,2*j) = -2 * (mod(i, 2) ==  mod(j, 2)) + 1;
            else
                grid(2*i,2*j+1) = 2 * (mod(i, 2) ==  mod(j, 2)) - 1;
            end
        end
    end

end