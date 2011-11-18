% evaluate the solution to the helmholtz equation given by
% c1*cos(kr) + c2*sin(kr)
% at r

function [v] = rpw(cs, ks, r)
    % cs 2xp
    % ks 2xp
    % r Nx2
    % v Nx1
    v = zeros(size(r,1), 1);
    
    for i=1:size(cs,2)
        v = v + cs(1,i) * cos(r*ks(:,i)) + cs(2,i) * sin(r*ks(:,i));
    end
    v = v ./ (sqrt(sum(sum(abs(cs)))/4));
end