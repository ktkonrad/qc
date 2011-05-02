function [F] = circle_blobs(t,a,R,ppw,n,dk,e)
%
% this func spits out real F array, given blob-ampl data (t, a(t)), t in S^1.
% Ampl a gives the sum over the F array of each blob (there's no dt weights).
% Necessary since filling complex F is several times slower than real!
% ppw = spatial points per wavelength, e = gaussian err param, dk = wave# grid
% n = grid size, R = wave# magnitude (eg 1)
% barnett 12/1/06

F = zeros(n,n);
if numel(R)==1, R = R*ones(size(t)); end     % vectorize the radii
N = numel(t); tic;
for i=1:N       % loop over angular directions, add a gaussian blob at each
  ki = R(i)*[cos(t(i)) sin(t(i))];                            % wavevector
  [l1,v1] = src(n, dk + ppw/2 + ki(1), 1/dk, e);   % find coords to add blob
  [l2,v2] = src(n, dk + ppw/2 + ki(2), 1/dk, e);
  F(l1,l2) = F(l1,l2) + (a(i)/(2*pi*e*e))*kron(v1',v2); % add blob
end
disp(sprintf('circle_blobs: %d blobs size %dx%d in %g sec',...
             N,numel(l1),numel(l2),toc));
%end
