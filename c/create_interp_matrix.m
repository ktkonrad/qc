function [] = create_interp_matrix(k, dx, M, upsample, outfile)
    addpath('../interpolation');
    points = [ -.5   .5
                .5   .5
               -.5  -.5
                .5  -.5
               -.5  1.5
                .5  1.5
              -1.5   .5
              -1.5  -.5
               -.5 -1.5
                .5 -1.5
               1.5   .5
               1.5  -.5
              -1.5  1.5
               1.5  1.5
              -1.5 -1.5
               1.5 -1.5
               -.5  2.5
                .5  2.5
              -2.5   .5
              -2.5  -.5
               -.5 -2.5
                .5 -2.5
               2.5   .5
               2.5  -.5] * dx; % 4x4 plus 2 extra per side
	       
    x2 = -.5*dx:dx/upsample:.5*dx;
    xs2 = meshgrid(x2);
    ys2 = flipud(meshgrid(x2)');
    points2 = [xs2(:) ys2(:)]; % fine grid
               
    A = interp_matrix(k, points, points2, M);
    dlmwrite(outfile, A);
   
    exit(0);
end
