function [A] = read_stencil(filename)
    addpath('../interpolation');
    vals = dlmread(filename);
    A = zeros(6);
    N = 24;
    s = stencil() + repmat([3.5, 3.5], N, 1);
    s = [s(:,1) 7-s(:,2)];
    for i=1:N
        A(s(i,1),s(i,2)) = vals(i);
    end
end