% function [N] = show_geom(head {, opts})
%
% show billiard and basis geometry.
%
% opts.fig gives figure window to add to
% opts.basis (default 1) flag to say if add basis.
% opts.sense (default 1) flag to say if label start, end
% opts.box (default 1) show bound box
% opts.i (default 0) >0 to say label the bdry pts & basis pts, skip amount
% opts.bdry (default 1) flag to say if show bdry pts & normals
% opts.leg (default 1) flag to say if show legend
%
% Vergini package.
%
% Barnett 8/21/03
% 11/17/03 added fig
% 1/17/04 allowed b_pts file to be missing without raising error
% 7/23/04 opts added
% 2/17/05 more opts. Bug: no elegant way to make legend reflect all opts.

function [N] = show_geom(head, opts)

if (nargin<2) % defaults
  opts.basis = 1;
  opts.sense = 1;
  opts.box = 1;
end
if ~isfield(opts, 'basis')
  opts.basis = 1;
end
if ~isfield(opts, 'sense')
  opts.sense = 1;
end
if ~isfield(opts, 'box')
  opts.box = 1;
end
if ~isfield(opts, 'bdry')
  opts.bdry = 1;
end
if ~isfield(opts, 'leg')
  opts.leg = 1;
end
if ~isfield(opts, 'i')
  opts.i = 0;
end

% load everything:

if opts.box
  % box
  [boxx, boxy] = textread([head '.bound_box'], '%f %f', 'headerlines', 4);
end  
if opts.bdry
  % matching points and normal vectors...
  [x, y, nx, ny] = textread(strcat(head,'.m_pts'), '%f %f %f %f', ...
			  'headerlines', 4);
  M = size(x,1); % # boundary points
end
  
basis = opts.basis;
if (basis)
  % basis...
  basis = exist([head '.b_pts'], 'file');
  if basis
    [t, bx, by, bnx, bny] = textread([head '.b_pts'], '%d %f %f %f %f', ...
                                     'headerlines', 4);
    N = size(t,1); % basis size
  end
end

% plot:
if (~isfield(opts, 'fig'))
  figure;
  if opts.box
    % box...
    h_box = line(boxx', boxy', 'color', [0.7 0.7 0.7]);
  end
else
  figure(opts.fig);
  set(gca, 'ydir', 'normal');
  if opts.box
    hold on; h_box = line(boxx', boxy', 'color', [0.7 0.7 0.7]);
    hold off;
  end
end

if opts.bdry
  % boundary...
  hold on; h_m = plot(x, y, 'k.');
  if opts.i
    for i=1:opts.i:M 
      text(x(i), y(i), sprintf('%d', i));
    end
  end
  hold off;
  axis equal
  n_len = 0.5;   % 0.5, length to show normals
  h_n = line([x x+n_len*nx]', [y y+n_len*ny]','color',[0 0 1]);
end

if opts.sense
  text(x(1), y(1), 'start');
  text(x(M), y(M), 'end');
end

if basis
  % basis...
  b_len = 0.5;   % length to show basis directions
  oyon_len = 0.02;  % length to show oyo basis normals
  for j=1:N
    switch t(j)
      
     case {0,1,2,3},   % RPW, (RE, IM, O-O, E-E)
      h_b = line([0;b_len*bnx(j)], [0;b_len*bny(j)], 'color', [0 1 0]);
     case 4,         % EPW RE
      h_b = line([bx(j);bx(j)+b_len*bnx(j)], [by(j);by(j)+b_len*bny(j)], ...
		 'color', [0 1 0]);
     case 5,         % EPW IM
      h_b = line([bx(j);bx(j)+b_len*bnx(j)], [by(j);by(j)+b_len*bny(j)], ...
		 'color', [0 1 1]);
     case {6,7,8}         % Y0 and its odd-odd, even-even versions
      hold on; plot(bx(j), by(j), 'k.', 'markersize', 15); hold off;
      h_b = line([bx(j);bx(j)+oyon_len*bnx(j)], ...
		 [by(j);by(j)+oyon_len*bny(j)], ...
		 'color', [1 1 0]);
      if opts.i
        for i=1:opts.i:N 
          text(bx(i), by(i), sprintf('%d', i));
        end
      end
      
     otherwise,
      disp(sprintf('show_geom: basis element %d has invalid sfunc type %d',...
		   j, t(j)));
      
    end
  end
  if opts.leg
    legend([h_box h_m h_n(1) h_b], 'bound box', 'boundary \Gamma', ...
           'outward normal', 'basis');
    title(sprintf('Vergini: box, billard, basis geom of %s', head));
  end
else
  if opts.leg
    % adjust labels since no basis shown...
    legend([h_box h_m h_n(1)], 'bound box', 'boundary \Gamma', ...
           'outward normal');
    title(sprintf('Vergini: box, billard geom of %s', head));
  end
end

% end
