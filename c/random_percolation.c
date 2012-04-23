#include <stdio.h>
#include <stdlib.h>

/*
generates m x n grid using random percolation model
(m and n will be rounded down to a multiple of 4)

inputs: grid - 2d array to be filled
        m    - number of rows in grid
        n    - number of columns in grid

outputs: grid - 2d array filled according to random percolation model
*/
void randomPercolation(double **grid, int m, int n) {
  int i, j, k, r;
  for (i = 0 ; i < m-3 ; i+=4) {
    for (j = 0 ; j < n-3 ; j+=4) {
      grid[i][j] = 1;
      grid[i][j+1] = 1;
      grid[i][j+2] = -1;
      grid[i][j+3] = -1;
      grid[i+1][j] = 1;
      grid[i+1][j+1] = 1;
      grid[i+1][j+2] = -1;
      grid[i+1][j+3] = -1;
      grid[i+2][j] = -1;
      grid[i+2][j+1] = -1;
      grid[i+2][j+2] = 1;
      grid[i+2][j+3] = 1;
      grid[i+3][j] = -1;
      grid[i+3][j+1] = -1;
      grid[i+3][j+2] = 1;
      grid[i+3][j+3] = 1;
    }
  }

  k = 0;
  for (i = 0 ; i < m/4*2 ; i++) {
    for (j = 0 ; j < n/4*2 ; j++) {
      if (k == 0) { // generate new random number r
        k = 32;
        r = rand();
      }
      if ((r >> --k) & 1) { // check kth bit of r
        grid[2*i+1][2*j+1] = 2 * ((i ^ j) & 1) - 1; // connect top left to bottom right
      }
      else {
	grid[2*i+1][2*j+2] = -2 * ((i ^ j) & 1) + 1; // connect top right to bottom left
      }
    }
  }  
}
/*
i\j  0 1 2 3 4 5 6 7
     0   1   2   3
0 0  + + - - + + - - 
1    + + - - + + - - 
2 1  - - + + - - + + 
3    - - + + - - + + 
4 2  + + - - + + - - 
5    + + - - + + - - 
6 3  - - + + - - + + 
7    - - + + - - + + 

 */
