/*
functions to count nodal domains

Kyle Konrad
3/29/2011

*/

#include "count_nodal_domains.h"
#include "util/stack2.h"
#include "util/interp_matrix.h"
#include "util/util.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <gsl/gsl_blas.h>

#define SIGN(x) ((x) > 0 ? 1 : -1)

// orthogonal connection directions
#define LEFT 1
#define RIGHT 2
#define ABOVE 4
#define BELOW 8

// orthonognal connection getting macros
#define GET_LEFT(c) (c & LEFT)
#define GET_RIGHT(c) (c & RIGHT)
#define GET_ABOVE(c) (c & ABOVE)
#define GET_BELOW(c) (c & BELOW)

// orthonognal connection setting macros
#define SET_LEFT(c) c |= LEFT
#define SET_RIGHT(c) c |= RIGHT
#define SET_ABOVE(c) c |= ABOVE
#define SET_BELOW(c) c |= BELOW

// orthonognal connection resetting macros
#define RESET_LEFT(c) c &= ~LEFT
#define RESET_RIGHT(c) c &= ~RIGHT
#define RESET_ABOVE(c) c &= ~ABOVE
#define RESET_BELOW(c) c &= ~BELOW

// diagonal connection directions
#define BR_CONNECTED (INT_MAX - 1) // connected to below right
#define AR_CONNECTED (INT_MAX - 2) // connected to above right
#define AL_CONNECTED (INT_MAX - 3) // connected to above left
#define BL_CONNECTED (INT_MAX - 4) // connected to below left
#define BR_DISCONNECTED (INT_MAX - 5) // disconnected from below right
#define AR_DISCONNECTED (INT_MAX - 6) // disconnected from above right
#define AL_DISCONNECTED (INT_MAX - 7) // disconnected from above left
#define BL_DISCONNECTED (INT_MAX - 8) // disconnected from below left

// diagonal connection getting macros
#define GET_BR(c) (c == BR_CONNECTED)
#define GET_AR(c) (c == AR_CONNECTED)
#define GET_AL(c) (c == AL_CONNECTED)
#define GET_BL(c) (c == BL_CONNECTED)

// diagonal connection setting macros
#define SET_BR(c) (c = BR_CONNECTED)
#define SET_AR(c) (c = AR_CONNECTED)
#define SET_AL(c) (c = AL_CONNECTED)
#define SET_BL(c) (c = BL_CONNECTED)

// diagonal connection resetting macros
#define RESET_BR(c) (c = BR_DISCONNECTED)
#define RESET_AR(c) (c = AR_DISCONNECTED)
#define RESET_AL(c) (c = AL_DISCONNECTED)
#define RESET_BL(c) (c = BL_DISCONNECTED)

// special values
#define UNCOUNTED INT_MAX
#define MASKED INT_MIN

// boolean functions to check for counted or masked
#define IS_COUNTED(c) (c < BL_DISCONNECTED)
#define IS_MASKED(c) (c == MASKED)
#define IS_INTERPOLATED(c) (c >= BL_DISCONNECTED && c != UNCOUNTED)

/*
count the number of nodal domains in grid
precondition: grid must be ny x nx

inputs:
        grid   - function values sampled at grid points
	mask   - array that defines boundaries of grid (no boundaries if NULL)
        nx     - number of samples in x-direction for grid
        ny     - number of samples in y-direction for grid
        k        - wavenumber of eigenfunction being interpolated
        dx       - sampled resoultion of eigenfunction
        M        - highest order bessel function to do
	upsample - upsampling ratio for interpolation
output: return value - count of nodal domains
*/
int countNodalDomainsInterp(double **grid, char **mask, int ny, int nx, double k, double dx, int M, int upsample, interp_stats *stats) {
  int i, j;
  int rc;

  gsl_matrix *interp = create_interp_matrix(k*dx, M, upsample);

  int **counted = imatrix(ny, nx);
  for (i = 0 ; i < ny ; i++) {
    for (j = 0 ; j < nx ; j++) {
      counted[i][j] = UNCOUNTED;
    }
  }

  if (mask != NULL) {
    applyMask(grid, counted, mask, ny, nx);
  }

  // array2file(grid, ny, nx, "../data/masked.dat");

  int nd = 0; // count of nodal domains
 
  i = 0;
  j = 0;
  while (findNextUnseen(counted, &i, &j, ny, nx)) {
    nd++;
    findDomainInterp(grid, counted, i, j, nd, ny, nx, upsample, interp, stats);
  }

  // TODO: verbosity global to control this output
  //intArray2file(counted, ny, nx, "../data/counted.dat");
			       
  free_imatrix(counted);

  gsl_matrix_free(interp);
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
      if (!IS_COUNTED(counted[r][c])) {
	*i = r;
	*j = c;
	return 1;
      }
    }
  }
  
  return 0; // all points have been seen
}

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
         updates stats
*/
void findDomainInterp(double **grid, int **counted, int i, int j, int nd, int ny, int nx, int upsample, gsl_matrix *interp, interp_stats *stats) {
  stack *s = newStack();
  push(s, j, i);
  
  int x, y;
  int currentSign;
  char connections; // keep track of what (x,y) is connected to
  int size = 0;

  while (pop(s, &x, &y)) {
    size++;
    currentSign = SIGN(grid[i][j]);
    connections = 0;

    // orthongal directions
    // left
    if (x >= 1 && !IS_MASKED(grid[y][x-1])) {
      if (SIGN(grid[y][x-1]) == currentSign) {
	SET_LEFT(connections);
	if(!IS_COUNTED(counted[y][x-1])) {
	  push(s, x - 1, y);
	}
      }
    }

    // above
    if (y >= 1 && !IS_MASKED(grid[y-1][x])) {
      if (SIGN(grid[y-1][x]) == currentSign) {
	SET_ABOVE(connections);
	if(!IS_COUNTED(counted[y-1][x])) {
	  push(s, x, y-1);
	}
      }
    }

    // right
    if (x < nx-1 && !IS_MASKED(grid[y][x+1])) {
      if (SIGN(grid[y][x+1]) == currentSign) {
	SET_RIGHT(connections);
	if(!IS_COUNTED(counted[y][x+1])) {
	  push(s, x+1, y);
	}
      }
    }

    // below
    if (y < ny-1 && !IS_MASKED(grid[y+1][x])) {
      if (SIGN(grid[y+1][x]) == currentSign) {
	SET_BELOW(connections);
	if(!IS_COUNTED(counted[y+1][x])) {
	  push(s, x, y+1);
	}
      }
    }
    // end orthongal directions

    // diagonal directions
    // above left
    if (x > 0 && y > 0 && !GET_ABOVE(connections) && !GET_LEFT(connections)) {
      if (!IS_MASKED(grid[y-1][x-1])) {
	if (SIGN(grid[y-1][x-1]) == currentSign) {
	  if (!IS_INTERPOLATED(counted[y][x])) { // check if interpolation results are memoized
	    interpolate(grid, counted, y-1, x-1, ny, nx, upsample, interp, stats);
	  }
	  if (counted[y][x] == AL_CONNECTED) {
	    push(s, x-1, y-1);
	  }
	}
      }      
    }

    // below left
    if (x > 0 && y < ny-1 && !GET_BELOW(connections) && !GET_LEFT(connections)) {
      if (!IS_MASKED(grid[y+1][x-1])) {
	if (SIGN(grid[y+1][x-1]) == currentSign) {
	  if (!IS_INTERPOLATED(counted[y][x])) {
	    interpolate(grid, counted, y, x-1, ny, nx, upsample, interp, stats);
	  }
	  if (counted[y][x] == BL_CONNECTED) {
	    push(s, x-1, y+1);
	  }
	}
      }      
    }

    // below right
    if (x < nx-1 && y < ny-1 && !GET_BELOW(connections) && !GET_RIGHT(connections)) {
      if (!IS_MASKED(grid[y+1][x+1])) {
	if (SIGN(grid[y+1][x+1]) == currentSign) {
	  if (!IS_INTERPOLATED(counted[y][x])) {
	    interpolate(grid, counted, y, x, ny, nx, upsample, interp, stats);
	  }
	  if (counted[y][x] == AL_CONNECTED) {
	    push(s, x+1, y+1);
	  }
	}
      }      
    }

    // above right
    if (x < nx-1 && y > 0 && !GET_ABOVE(connections) && !GET_RIGHT(connections)) {
      if (!IS_MASKED(grid[y-1][x+1])) {
	if (SIGN(grid[y-1][x+1]) == currentSign) {
	  if (!IS_INTERPOLATED(counted[y][x])) {
	    interpolate(grid, counted, y-1, x, ny, nx, upsample, interp, stats);
	  }
	  if (counted[y][x] == AL_CONNECTED) {
	    push(s, x+1, y-1);
	  }
	}
      }      
    }
    counted[y][x] = currentSign == 1 ? nd : -nd;
  }
  if (size < SMALL_DOMAIN_SIZE) {
    stats->small_domain_count++;
  }
  destroyStack(s);
}

/*
  interpolate a region of the grid and put the results in counted
  
  precondition: grid and counted are ny x nx

  inputs:
        grid     - 2d array of function values
	counted  - 2d array indicating whether a point has been counted
	i        - row of initial point in grid
        j        - column of initial point in grid
	ny       - rows in grid
	nx       - columns in grid
	upsample - upsampling rate
	interp   - precalculated 24 x (upsapmle+1)^2 matrix to do interpolation
	
  outputs:
        stores connection values in counted[j:j+1][i:i+1]
*/
#define IDX(x, y) ((x)*n+(y))
void interpolate(double **grid, int **counted, int i, int j, int ny, int nx, int upsample, gsl_matrix *interp, interp_stats *stats) {
  int k;
  int x, y;
  int n = upsample + 1;
  int rc;
  int currentSign;
  gsl_vector *interp_input, *interp_output;
  int tl_br_connected = 0; // flag indicated whether (x,y) is connected to (x+1,y+1)
  int **interp_counted = imatrix(ny, nx);
  stack *s = newStack();


  // validate size of interp matrix
  if (interp->size1 != n*n || interp->size2 != NUM_STENCIL_POINTS) {
    ERROR("invalid size for interpolation matrix");
    exit(INTERP_ERR);
  }

  // check we're not too close to the edge
  if (i < 2 || i >= ny - 2 || j < 2 || j >= nx - 2) {
    ERROR("trouble spot near edge: (x,y) = (%d,%d)", j, i);
    stats->edge_trouble_count++;
    tl_br_connected = rand() % 2; // we can't interpolate so guess randomly which way it's connected
  } else {
    

    stats->interp_count++;
    interp_input = gsl_vector_alloc(NUM_STENCIL_POINTS);
    interp_output = gsl_vector_alloc(interp->size1);


    /*
      populate the vector of stencil points
      stencil points are laid out as follows:


      y\x  -2  -1   0   1   2   3
      +--------------------------
      -2 |          00  01
      -1 |      02  03  04  05
      0 |  06  07  08  09  10  11
      1 |  12  13  14  15  16  17
      2 |      18  19  20  21
      3 |          22  23
  
      point 8, the top left point of the inner square, is at (x,y)=(0,0)
    */
    gsl_vector_set(interp_input, 0 , grid[i-2][j]);
    gsl_vector_set(interp_input, 1 , grid[i-2][j+1]);
    gsl_vector_set(interp_input, 2 , grid[i-1][j-1]);
    gsl_vector_set(interp_input, 3 , grid[i-1][j]);
    gsl_vector_set(interp_input, 4 , grid[i-1][j+1]);
    gsl_vector_set(interp_input, 5 , grid[i-1][j+2]);
    gsl_vector_set(interp_input, 6 , grid[i][j-2]);
    gsl_vector_set(interp_input, 7 , grid[i][j-1]);
    gsl_vector_set(interp_input, 8 , grid[i][j]);
    gsl_vector_set(interp_input, 9 , grid[i][j+1]);
    gsl_vector_set(interp_input, 10, grid[i][j+2]);
    gsl_vector_set(interp_input, 11, grid[i][j+3]);
    gsl_vector_set(interp_input, 12, grid[i+1][j-2]);
    gsl_vector_set(interp_input, 13, grid[i+1][j-1]);
    gsl_vector_set(interp_input, 14, grid[i+1][j]);
    gsl_vector_set(interp_input, 15, grid[i+1][j+1]);
    gsl_vector_set(interp_input, 16, grid[i+1][j+2]);
    gsl_vector_set(interp_input, 17, grid[i+1][j+3]);
    gsl_vector_set(interp_input, 18, grid[i+2][j-1]);
    gsl_vector_set(interp_input, 19, grid[i+2][j]);
    gsl_vector_set(interp_input, 20, grid[i+2][j+1]);
    gsl_vector_set(interp_input, 21, grid[i+2][j+2]);
    gsl_vector_set(interp_input, 22, grid[i+3][j]);
    gsl_vector_set(interp_input, 23, grid[i+3][j+1]);
  
    // check that we are not near the boundary
    // TODO: figure out what to do if we are near the boundary
    for (k = 0 ; k < interp_input->size ; k++) {
      if (IS_MASKED(gsl_vector_get(interp_input, k))) {
        ERROR("trouble spot near boundary: (x,y) = (%d,%d)", j, i);
        stats->boundary_trouble_count++;
        break;
      }
    }

    // interp_output = interp * interp_input
    rc = gsl_blas_dgemv(CblasNoTrans, 1, interp, interp_input, 0, interp_output);
    if (rc) {
      ERROR("interpolation failed. gsl_blas_dgemv returned %d", rc);
      exit(INTERP_ERR);
    }
  
    push(s, 0, 0);
    currentSign = SIGN(gsl_vector_get(interp_output, 0));

    while (pop(s, &x, &y)) {
      if (x == n && y == n) { // we got to the bottom right
        tl_br_connected = 1;
        break;
      }
    
      interp_counted[y][x] = 1;
    
      // left
      if (x >= 1) {
        if (SIGN(gsl_vector_get(interp_output, IDX(x-1,y))) == currentSign && !interp_counted[y][x-1]) {
          push(s, x-1, y);
        }
      }

      // above
      if (y >= 1) {
        if (SIGN(gsl_vector_get(interp_output, IDX(x,y-1))) == currentSign && !interp_counted[y-1][x]) {
          push(s, x, y-1);
        }
      }

      // right
      if (x < n-1) {
        if (SIGN(gsl_vector_get(interp_output, IDX(x+1,y))) == currentSign && !interp_counted[y][x+1]) {
          push(s, x+1, y);
        }
      }

      // below
      if (y < n-1) {
        if (SIGN(gsl_vector_get(interp_output, IDX(x,y+1))) == currentSign && !interp_counted[y+1][x]) {
          push(s, x, y+1);
        }
      }
    }

    /*
    // debug output
    char filename[50];
    sprintf(filename, "../data/interpolated_%d_%d.dat", j, i);
    FILE *outfile = fopen(filename, "w");
    gsl_vector_fprintf(outfile, interp_output, "%.16g");
    fclose(outfile);
    */

    // cleanup
    destroyStack(s);
    gsl_vector_free(interp_input);
    gsl_vector_free(interp_output);
    free_imatrix(interp_counted);

  }

  // set appropriate values in counted
  if (tl_br_connected) {
    if (!IS_COUNTED(counted[i][j])) SET_BR(counted[i][j]);
    if (!IS_COUNTED(counted[i+1][j+1])) SET_AL(counted[i+1][j+1]);
    if (!IS_COUNTED(counted[i+1][j])) RESET_AR(counted[i+1][j]);
    if (!IS_COUNTED(counted[i][j+1])) RESET_BL(counted[i][j+1]);
  } else {
    if (!IS_COUNTED(counted[i][j])) RESET_BR(counted[i][j]);
    if (!IS_COUNTED(counted[i+1][j+1])) RESET_AL(counted[i+1][j+1]);
    if (!IS_COUNTED(counted[i+1][j])) SET_AR(counted[i+1][j]);
    if (!IS_COUNTED(counted[i][j+1])) SET_BL(counted[i][j+1]);
  }
}
