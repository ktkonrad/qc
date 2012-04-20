% test spray calling, barnett July 2006

% TEST GEOM ----------------------------------------------------------
% geom = '+1 2 0 1     0  0 -1     2  0 1 -1 0'; % half lemon (ccw sense)
% half mushroom... (cw sense)
%geom = '-1 5 0 0   0 0 -1  0 -.5 -1  0 -.5 1.5  -2 1 0 -.5 0  0 0 0';
% qusta:2
geom = '-1 4 0 1   0 1 1  -2 2 0 1 0  0 0 0  0 0 1';

cmd = sprintf('echo %s | ./spray -m 200', geom);    % send 
[s,out] = system([cmd ' 2>/dev/null']);             % kill stderr output

L = strread(out, '%f', 1);
[bx,by,bnx,bny] = strread(out, '%f%f%f%f');   % like textread for char array
figure; plot(bx, by, '.'); hold on; plot([bx bx+bnx]', [by by+bny]', '-');
axis equal

% TEST GEOM & TRAJ ----------------------------------------------------------
geom = '-1 5 0 0   0 0 -1  0 -.5 -1  0 -.5 1.5  -2 1 0 -.5 0  0 0 0'; %hamush
%geom = '+1 2 0 1     0  0 -1     2  0 1 1 0'; % trunc ball - bad |v| growth
%geom = '+1 4 1 1   -2 -1 1 0 4  -2 -1 -1 -5 0  -2 1 -1 0 -6  -2 1 1 3 0'; %grs

cmd = sprintf('echo %s | ./spray -m 200', geom);
[s,out] = system([cmd ' 2>/dev/null']);             % kill stderr output

L = strread(out, '%f', 1);
[bx,by,bnx,bny] = strread(out, '%f%f%f%f');   % like textread for char array
cmd = sprintf('echo %s | ./spray -T 1.79:0.705:10', geom);    % send 
%cmd = sprintf('echo %s | ./spray -T .117:.223:1000', geom);
%cmd = sprintf('echo %s | ./spray -T .1.975411:0.3882601:10', geom);
[s,out] = system([cmd ' 2>/dev/null']);

[x,y,vx,vy,obj,l] = strread(out, '%f%f%f%f%d%f');
figure; plot(x, y, '.-'); axis equal; hold on;
plot(bx, by, 'g.', 'markersize', 1); plot([x x+.1*vx]', [y y+.1*vy]', 'r-');
title('single trajectory'); xlabel('x'); ylabel('y');

% problem coords: q = 1.97541138165694 p = 0.388260139190508
figure; semilogy(abs(vx.^2+vy.^2-1));   % grows by factor 10 each arc bounce!
 

% TEST SPRAY --------------------------------------------------------------

geom = '-1 5 0 0   0 0 -1  0 -.5 -1  0 -.5 1.5  -2 1 0 -.5 0  0 0 0'; %hamush
cmd = sprintf('echo %s | ./spray -S 2.8:100:.01:10', geom);    % send 
[s,out] = system([cmd ' 2>/dev/null']);

[p,j,l] = strread(out, '%f%d%f');

% TEST POs --------------------------------------------------------------

geom = '-1 5 0 0   0 0 -1  0 -.5 -1  0 -.5 1.5  -2 1 0 -.5 0  0 0 0'; %hamush
cmd = sprintf('echo %s | ./spray -P 2.8:100:.01:.5:10', geom);    % send 
[s,out] = system([cmd ' 2>/dev/null']);

[p,j,l,dr2,dv2] = strread(out, '%f%d%f%f%f')



% GENERATE CLASS RETURN TIME PIC -----------------------------------7/27/06
%geom = '-1 5 0 0   0 0 -1  0 -.5 -1  0 -.5 1.5  -2 1 0 -.5 0  0 0 0'; %hamush
geom = '-1 4 0 1   0 1 1  -2 2 0 1 0  0 0 0  0 0 1';    % qusta:2
cmd = sprintf('echo %s | ./spray -q -m 1', geom);
[s,out] = system([cmd ' 2>/dev/null']);

L = strread(out, '%f', 1);
[bx,by,bnx,bny] = strread(out, '%f%f%f%f');            % clear stdin

%qs = 0:0.01:L-2;    % 0:0.003:L-2                       for hamush
qs = 0:0.02:L;
nq = numel(qs);
maxnrets = 10000;    % 10000;
rets = NaN*ones(maxnrets, nq); prets = rets; drrets = rets; dvrets = rets;
nrets = ones(1, nq);
for m=1:nq
  q = qs(m) %q = qs(m)+1        % convert from full q to hamush coord q
  cmd = sprintf('echo %s | ./spray -q -S %.15g:10000:.02:15', geom, q); 
  %cmd = sprintf('echo %s | ./spray -q -P %.15g:10000:.03:.03:40', geom, q); 
  [s,out] = system([cmd ' 2>/dev/null']);
  [p,j,l,dr2,dv2] = strread(out, '%f%d%f%f%f');
  nrets(1, m) = numel(l);
  nl = min(numel(l),maxnrets);             % number of l's to keep
  rets(1:nl,m) = l(1:nl); prets(1:nl,m) = p(1:nl);
  drrets(1:nl,m) = sqrt(dr2(1:nl)); dvrets(1:nl,m) = sqrt(dv2(1:nl));
end                            % don't paste beyond here
sum(nrets)       % how many collected total
figure; plot(qs, rets(1:min(max(nrets),maxnrets),:), 'b.', 'markersize', 1);
xlabel('q'); ylabel('l'); axis tight;
%save mu_spray1e4_.003_30.mat
%save mu_spray1e5_.01_20.mat
save qusta2_spray1e4_.02_18.mat


% convert to greyscale density image giving density of returns in [l,l+dl]
% for each q...
dl = 0.02;
ls = 0:dl:15;
nl = numel(ls);
%C = histc(rets, ls);
C = histc(rets + 0./(dvrets<0.02) + 0./(drrets<0.02), ls);
sum(sum(C))
figure;
imagesc(qs, ls, -sqrt(C)); set(gca,'ydir','normal');xlabel('q');ylabel('l');
colormap(gray(256)); c=caxis; caxis(c/3);
title('qusta:2 POs: dr=.02 dv=.02');
%print -deps2 qusta2_spray1e4_.02_15.eps


% density image in PSOS, within some return length range...
dp = 0.005;
ps = -1:dp:1; np = numel(ps);
C = histc(prets + 0./(dvrets<0.02) + 0./(rets<6), ps); sum(sum(C))
figure;
imagesc(qs, ps, -sqrt(C)); set(gca,'ydir','normal');xlabel('q');ylabel('p');
colormap(gray(256)); c=caxis; caxis(c/3);
title('qusta:2 POs: l<6 dv=.02');


% 3D plot in (q,p,l) space
figure;
plot3(repmat(qs, [maxnrets 1]), prets + 0./(dvrets<0.02), rets, '.', 'markersize', 1);
xlabel('q');ylabel('p');zlabel('l'); axis vis3d

% pare down and reformat as list of (q,l) pairs...
[r,c] = find(rets<18);
qrets = qs(c);
clear r c
rets = squeeze(rets(find(rets<18)))';
figure; plot(qrets, rets, 'm.', 'markersize', 1);
xlabel('q'); ylabel('l'); axis tight;
save mu_spray1e4_.003_18r.mat


%end

