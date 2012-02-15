function points = stencil()

    n = 4;

    xmin = -(floor(n)-1)/2;
    xmax = (floor(n)-1)/2;

    x = xmin:1:xmax;
    [xs, ys] = meshgrid(x);
    points = [xs(:) ys(:)];
    points = vertcat(points, [ -.5  2.5
              .5  2.5
            -2.5   .5
            -2.5  -.5
             -.5 -2.5
              .5 -2.5
             2.5   .5
             2.5  -.5 ]);