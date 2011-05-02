% function [ks, mass, ten, nrm, ne, perim, area, per, ngr, xoff,
% dx] = load_mc(tarhead {, direc})
%
% load in possibly-tar'ed mass_chop datafiles
% 2/16/04 included chop-part only per, ngr if present
% 2/19/04: if direc set, use existing untared data files, don't erase them


function [ks, mass, ten, nrm, ne, perim, area, per, ngr, xoff, dx] ...
    = load_mc(tarhead, direc)

ks = []; mass = []; ten = []; nrm = []; per = []; ngr = [];
perim = 0; area = 0;
ne = 0;

untarflag = (nargin<2);

if untarflag
  % unpack tarred data
  if ~exist([tarhead '.tar'], 'file')
    disp(['load_mc: ' tarhead '.tar not found!']);
    return;
  end
  
  % use the tar head to name a tmp directory
  tmp = ['tmp_' tarhead];
  
  system(['rm -Rf ' tmp]);
  system(['mkdir ' tmp]);
  system(['cp ' tarhead '.tar ' tmp]);
  [status, tarlist] = system(['cd ' tmp '; tar xvf ' tarhead '.tar']);
else
  % use already-unpacked data (tarhead argument not used)
  tmp = direc;
end
a = dir(tmp);
  
for i=1:length(a)
  if ~isempty(strfind(a(i).name, '.cf'))  % if cf file, get props
    head = [tmp '/' a(i).name(1:length(a(i).name)-3)];  % remove suffix
    [perim, area] = load_props(head);
  end
  if ~isempty(strfind(a(i).name, '.sum'))  % if .sum file (not use findstr!)
    %[tmp '/' a(i).name]
    head = [tmp '/' a(i).name(1:length(a(i).name)-4)];  % remove suffix
    [ks_t, mass_t, ten_t, nrm_t, ne_t] = load_sum(head);
    ks = [ks; ks_t];
    mass = [mass; mass_t];
    ten = [ten; ten_t];
    nrm = [nrm; nrm_t];
    if exist([head '.chop.per'], 'file')
      [per_t dummy xoff dx] = load_1d(head, 'chop.per');  % get xoff etc
      per = [per; per_t];
    end
    if exist([head '.chop.ngr'], 'file')
      ngr_t = load_1d(head, 'chop.ngr');
      ngr = [ngr; ngr_t];
    end
    ne = ne + ne_t;
  end
end

if untarflag
  % clear up tmp direc (comment for debugging)
  system(['rm -Rf ' tmp]);
end

% now sort into inc k order
[ks, ind] = sort(ks);
mass = mass(ind);
ten = ten(ind);
nrm = nrm(ind);

% end



