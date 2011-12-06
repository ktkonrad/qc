/*
functions to count nodal domains

Kyle Konrad
3/29/2011

*/

#include "count_nodal_domains.h"
#include "util/stack2.h"
#include "util.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define SIGN(x) ((x) == 0 ? 0 : ((x) > 0 ? 1 : -1))

/*
count the number of nodal domains in grid
precondition: grid must be ny x nx

inputs:
        grid - function values sampled at grid points
	mask - array that defines boundaries of grid (no boundaries if NULL)
        nx   - number of samples in x-direction for grid
        ny   - number of samples in y-direction for grid

output: return value - count of nodal domains
*/
int countNodalDomains(double **grid, char **mask, int ny, int nx) {
  int i, j;

  int **counted = (int **)malloc(ny * sizeof(int *)); // keep track of whether a point has been counted in a nodal domain
  for (i = 0 ; i < ny ; i++) {
    counted[i] = (int *)malloc(nx * sizeof(int));
    for (j = 0 ; j < nx ; j++) {
      counted[i][j] = 0;
    }
  }

  if (mask != NULL)
    applyMask(grid, counted, mask, ny, nx);

  array2file(grid, ny, nx, "../data/masked.dat");


  int nd = 0; // count of nodal domains
 
  i = 0;
  j = -1;
  while (findNextUnseen(counted, &i, &j, ny, nx)) {
    nd++;
    findDomain(grid, counted, i, j, nd, ny, nx);
  }

  intArray2file(counted, ny, nx, "../data/counted.dat");

  for (i = 0 ; i < ny ; i++)
    free(counted[i]);
  free(counted);
  
  return nd;
}

/*
find next unseen value of grid
precondition: *i and *j nonnull

inputs:
        **counted - 2d array indiciating which points have been counted
        *i        - row to start searching from 
        *j        - column to start searching from
	ny        - rows in counted
        nx        - columns in counted

outputs: *i - row of next zero found
         *j - column of next zero found
	 return value - boolean indicating whether something was found
*/
int findNextUnseen(int **counted, int *i, int *j, int ny, int nx) {
  int r, c;

  for (r = *i ; r < ny ; r++) {
    for (c = 0 ; c < nx ; c++) {
      if (r == *i && c <= *j)
	continue;
      if (counted[r][c] == 0) {
	*i = r;
	*j = c;
	return 1;
      }
    }
  }
  
  return 0; // all points have been seen
}


/*
DEPRECATED

mark nodal domain containing grid[i][j] as counted
recursive version

inputs:
	i - row of initial point in grid
        j - column of initial point in grid

outputs:
         none

void findDomainRecursive(int i, int j) {
  counted[i][j] = 1;
  int currentSign = SIGN(grid[i][j]);
  
  // left
  if (j >= 1 && SIGN(grid[i][j-1]) == currentSign && !counted[i][j-1])
    findDomainRecursive(i, j - 1);

  // above
  if (i >= 1 && SIGN(grid[i-1][j]) == currentSign && !counted[i-1][j])
   findDomainRecursive(i - 1, j);

  // right
  if (j < nx-1 && SIGN(grid[i][j+1]) == currentSign && !counted[i][j+1])
    findDomainRecursive(i, j + 1);

  // below
  if (i < ny-1 && SIGN(grid[i+1][j]) == currentSign && !counted[i+1][j])
    findDomainRecursive(i + 1, j);
}
*/

/*
mark nodal domain containing grid[i][j] as counted
non-recursive version

inputs:
        grid    - 2d array of function values
	counted - 2d array indicating whether a point has been counted
	i       - row of initial point in grid
        j       - column of initial point in grid
	nd      - number of current nodal domain
	ny      - rows in grid
	nx      - columns in grid

precondition: grid and counted are ny x nx

outputs:
         none
*/
void findDomain(double **grid, int **counted, int i, int j, int nd, int ny, int nx) {

  stack *s = newStack();
  push(s, j, i);
  
  int x, y;
  int currentSign;

  int trouble; // flag to track whether we have found a trouble spot

  while (pop(s, &x, &y)) {
    
    currentSign = SIGN(grid[i][j]);
    counted[y][x] = currentSign == 1 ? nd : -nd;
    
    trouble = 0;

    // left
    if (x >= 1 && !isinf(grid[y][x-1])) {
      if (SIGN(grid[y][x-1]) == currentSign) {
	if(!counted[y][x-1]) {
	  push(s, x - 1, y);
	}
      }
    }

    // above
    if (y >= 1 && !isinf(grid[y-1][x])) {
      if (SIGN(grid[y-1][x]) == currentSign) {
	if(!counted[y-1][x]) {
	  push(s, x, y - 1);
	}
      }
    }

    // right
    if (x < nx-1 && !isinf(grid[y][x+1])) {
      if (SIGN(grid[y][x+1]) == currentSign) {
	if(!counted[y][x+1]) {
	  push(s, x + 1, y);
	}
      } else {
	trouble = 1;
      }
    }

    // below
    if (y < ny-1 && !isinf(grid[y+1][x])) {
      if (SIGN(grid[y+1][x]) == currentSign) {
	trouble = 0;
	if(!counted[y+1][x]) {
	  push(s, x, y + 1);
	}
      } else {
	trouble &= 1;
      }
    } else {
      trouble = 0;
    }

    // check below & right if we may have a trouble spot
    if (trouble && x < nx-1 && y < ny-1 && !isinf(grid[y+1][x+1]) && SIGN(grid[y+1][x+1]) == currentSign) { // we have at trouble spot
      //printf("trouble at (%d, %d)\n", x, y);
      // TODO: interpolate and put results in trouble array      
    }

  }
  destroyStack(s);
}
