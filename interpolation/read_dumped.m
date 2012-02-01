function [A] = read_dumped(filename, m, n)
    A = dlmread(filename);
    A = reshape(A, n, m)';
end