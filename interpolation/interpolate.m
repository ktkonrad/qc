function [interpolated] = interpolate(grid, row, col, interp_mat)
    
    
    inpoints = [ grid(row-2,col)
                 grid(row-2,col+1)
                 grid(row-1,col-1)
                 grid(row-1,col)
                 grid(row-1,col+1)
                 grid(row-1,col+2)
                 grid(row,col-2)
                 grid(row,col-1)
                 grid(row,col)
                 grid(row,col+1)
                 grid(row,col+2)
                 grid(row,col+3)
                 grid(row+1,col-2)
                 grid(row+1,col-1)
                 grid(row+1,col)
                 grid(row+1,col+1)
                 grid(row+1,col+2)
                 grid(row+1,col+3)
                 grid(row+2,col-1)
                 grid(row+2,col)
                 grid(row+2,col+1)
                 grid(row+2,col+2)
                 grid(row+3,col)
                 grid(row+3,col+1) ];

     outpoints = interp_mat * inpoints;
     n = sqrt(size(interp_mat,1));
     interpolated = reshape(outpoints, n, n);
             
end