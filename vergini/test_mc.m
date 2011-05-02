% test different mass-chop methods against each other.
%
% 3/22/04 barnett


k_lo = 100;
k_hi = 101;
sys = '-l qurf:0:0.2:0.2 -s oyooo:2.5:5:1';
head = 'mct';
b = 10;             % try different b values ................
dx = 0.01;
d_lo = 0.1;  % vergini window % 0.05
d_hi = 0.2;
verb = '-q';   % -q

for method=1:2
  
  ks = []; mass = []; ten = []; nrm = [];
  ne = 0;
  first_time = 1;
  
  tic;    
  for k=k_lo+d_lo:d_lo+d_hi:k_hi-d_hi % .................. collect states
    
    if first_time % output geom too, including mask
      mask_flag = '-m';
      first_time = 0;
    else
      mask_flag = '';
    end
    % output flags:    -f %g:0 -e 1 -p 0
    m = '';
    if method==2
      m = '-2';
    end
    cmd = sprintf(['verg %s %s -b %g -k %g -o %s -4 2 -C'...
                   '%g:%g:0.5:0.75 %s %s -p 0'], verb, sys, b, k, head, ...
                  d_lo, d_hi, mask_flag, m);
    disp(cmd);
    [s] = system(cmd);
    if (s~=0)
      disp('verg crashed!');
      return;
    end
    
    [ks_t, mass_t, ten_t, nrm_t, ne_t] = load_sum(head);
    ks = [ks; ks_t];
    mass = [mass; mass_t];
    ten = [ten; ten_t];
    nrm = [nrm; nrm_t];
    ne = ne + ne_t;
    ten_t
    
  end % ..................................................
  disp(sprintf('total time = %f minutes', toc/60));
  totaltime = toc;
  
  massm(:,method) = mass;
  nrmm(:,method) = nrm;
end

figure; plot(ks, massm(:,1)-massm(:,2), '+');
% errors between 2 different gen-BC formulae die quadratically, 1/b^2.
% about 1/1000 rel err at b=10.
% Only 1% faster, for small systems (k=100) tested!


