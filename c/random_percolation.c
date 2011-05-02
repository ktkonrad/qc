#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
generates m x n grid using random percolation model
(m and n will be rounded down to a multiple of 4)

inputs: grid - 2d array to be filled
        m    - number of rows in grid
        n    - number of columns in grid

outputs: grid - 2d array filled according to random percolation model
*/
void randomPercolation(double **grid, int m, int n) {
  srand(time(NULL));
  int i;
  int j;
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
  
  for (i = 0 ; i < m/4*2 ; i++) {
    for (j = 0 ; j < n/4*2 ; j++) {
      if (rand() % 2 == 0)
	grid[2*i+1][2*j+1] = -2 * (i % 2 == j % 2) + 1;
      else
	grid[2*i+1][2*j+2] = 2 * (i % 2 == j % 2) - 1;
    }
  }
  
}
