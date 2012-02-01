dxs = 10.^(-1:-1:-8);
dxs = 10.^(-1:-.25:-5);
k = 200;
M = 10;

grids = cell(5,1);
names = cell(5,1);

grids{1} = [ -.5   .5
              .5   .5
             -.5  -.5
              .5  -.5 ]; % 2x2
names{1} = '2x2';
         
grids{2} = vertcat(grids{1}, ...
           [ -.5  1.5
              .5  1.5
            -1.5   .5
            -1.5  -.5
             -.5 -1.5
              .5 -1.5
             1.5   .5
             1.5  -.5 ]); % 4x4 without corners
names{2} = '4x4 no corners';
         
grids{3} = vertcat(grids{2}, ...
           [ -1.5  1.5
              1.5  1.5
             -1.5 -1.5
              1.5 -1.5 ]); % 4x4
names{3} = '4x4';
          
grids{4} = vertcat(grids{3}, ...
           [ -.5  2.5
              .5  2.5
            -2.5   .5
            -2.5  -.5
             -.5 -2.5
              .5 -2.5
             2.5   .5
             2.5  -.5 ]); % 4x4 with 2 extra on each side
names{4} = '4x4 +2/side';
         
grids{5} = vertcat(grids{4}, ...
           [-1.5  2.5
             1.5  2.5
            -2.5  1.5
            -2.5 -1.5
            -1.5 -2.5
             1.5 -2.5
             2.5  1.5
             2.5 -1.5 ]); % 6x6 without corners
names{5} = '6x6 no corners';            
         
outx = -.5:.1:.5;
outxs = meshgrid(-.5:.1:.5);
outys = flipud(meshgrid(-.5:.1:.5)');
outgrid = [outxs(:), outys(:)];    
error_norms = zeros(length(dxs), length(grids));

i = 0;
for dx=dxs
    i = i + 1;
    error_norms(i,:) = plane_wave_interp2(k, cellfun(@(x) dx*x, grids, 'UniformOutput', false), dx*outgrid, M);
end

figure;
loglog(k*dxs, error_norms);
xlabel('k*dx');
ylabel('||errors||_{\infty}');
legend(names);