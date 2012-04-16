function [grids, names] = all_stencils()
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
end