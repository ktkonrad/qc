% investigate optimal C4 coeff for billiard shape
%
% 1/25/04 barnett

k = 300;
sys = '-l qugrs:1:0.4:0.7 -s oyooo:1.5:7:1';
head = 'c';
b = 10;
d_lo = 0.0;
d_hi = 0.2;
verb = '';   % -q

c4s = 0.5:0.05:1.0;
nc4s = length(c4s);

mne = 10;  % max ne, to fill array
ks = zeros(mne, nc4s); ten =ks; nmr = ks;

for i = 1:nc4s   % loop over C4 coeff .................................
  c4 = c4s(i);
  cmd = sprintf(['verg %s %s -b %g -k %g -o %s -V %g:%g -4 %g -p %d'], ...
		verb, sys, b, k, head, d_lo, d_hi, c4, 2*k);
  disp(cmd);
  [s] = system(cmd);
  if (s~=0)
    disp('verg crashed!');
    return;
  end
  [ks_t, kos, ten_t, nrm_t, ne] = load_sum(head);
  ks(1:ne,i) = ks_t;
  ten(1:ne,i) = ten_t;
  nrm(1:ne,i) = nrm_t;
  
end % .................................................................

figure;  % general tension vs delta
loglog(abs(ks - kos(1)), ten, '-+');
xlabel('\delta'); ylabel('tension');

figure;  % show tension minima vs C4
ns = 5; % number of stable states
plot(c4s, ten(1:ns,:)./kron(min(ten(1:ns,:),[],2), ones(1,nc4s)), '+-');
xlabel('C4'); ylabel('tension normalised by min value');
legnum(ks(1:ns,1) - kos(1));

figure;  % show computed delta vs C4
ns = 6;
plot(c4s, (ks(1:ns,:) - kron(ks(1:ns,1), ones(1,nc4s))), '+-');
axis([min(c4s) max(c4s) -1e-4 1e-4]);
xlabel('C4'); ylabel('computed \delta');
legnum(ks(1:ns,1) - kos(1));
