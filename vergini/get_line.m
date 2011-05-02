% function [l, chop_area] = get_line(head, a, t1, t2)
%
% return line start and end from .chop.m_pts file

function [l, chop_area] = get_line(head, a, t1, t2)

[x, y, nx, ny] = textread([head '.chop.m_pts'], '%f %f %f %f', ...
			  'headerlines', 4);
M = size(x,1); % # boundary points

% extrapolate to exact Gamma intersection points...
l(1) = x(1) - 0.5*(x(M)-x(1))/(M-1);
l(2) = y(1) - 0.5*(y(M)-y(1))/(M-1);
l(3) = x(M) + 0.5*(x(M)-x(1))/(M-1);
l(4) = y(M) + 0.5*(y(M)-y(1))/(M-1);

if nargin>=4
  % compute chopped area of qugrs...
  R2 = 1/sin(t2);
  trapez = (l(3)-l(1))*(l(2)+l(4))/2;   % trapezium part
  t = asin(l(4)/R2); % theta
  half_segment = R2*R2*(t - sin(2*t)/2)/2;
  chop_area = trapez - half_segment;
else
  chop_area = 0;
end
