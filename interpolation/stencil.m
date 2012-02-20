function stencil_points = stencil()
% 4x4 + 2 stencil in row major order

    npoints = 24;
    stencil_points = zeros(npoints, 2);
    row_widths = [2, 4, 6, 6, 4, 2];

    i = 1;
    for row=1:6
        row_edge = (row_widths(row) - 1) / 2;
        y = 3.5 - row;
        for x=-row_edge:row_edge
            stencil_points(i,:) = [x, y];
            i = i + 1;
        end
    end