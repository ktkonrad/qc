/*
functions to count nodal domains

Kyle Konrad
3/29/2011

*/

#include "count_nodal_domains.h"
#include "util/stack2.h"
#include "util/interp_matrix.h"
#include "util/util.h"
#include "util/count_util.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <gsl/gsl_blas.h>

#define SIGN(x) fabs(x)

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

extern int verb;


/*
count the number of nodal domains in grid
precondition: grid must be ny x nx

inputs:
        grid     - function values sampled at grid points
	counted  - array to store nodal domain numbers in. mask should be applied to this already
        nx       - number of samples in x-direction for grid
        ny       - number of samples in y-direction for grid
        k        - wavenumber of eigenfunction being interpolated
        dx       - sampled resoultion of eigenfunction
        M        - highest order bessel function to do
        upsample - upsampling ratio for interpolation
        sizefile - file to write sizes of domains to
output: return value - count of nodal domains
*/
int countNodalDomainsInterp(double **grid, int **counted, int ny, int nx, double alpha, int M, int upsample_ratio, interp_stats *stats, FILE *sizefile) {
  int nodal_domain_count;
  bit_array_t *upsampled_grid;
  
  upsampled_grid = upsample(grid, counted, ny, nx, alpha, M, upsample_ratio, stats);
  nodal_domain_count = countNodalDomains(upsampled_grid, counted, sizefile);
  free_bit_array(upsampled_grid);
  return nodal_domain_count;
}

/*
count the number of nodal domains in grid
precondition: grid must be ny x nx

inputs:
        grid   - function values sampled at grid points
	counted  - array to store nodal domain numbers in. mask should be applied to this already
        sizefile - file to write size of nodal domains to
output: return value - count of nodal domains
*/
int countNodalDomains(bit_array_t *signs, int **counted, FILE *sizefile) {
  int i, j;
  int rc;
  int nx, ny;
  nx = signs->nx;
  ny = signs->ny;

  #ifdef DEBUG
    bit_array2file(signs, "../data/signs.dat");
  #endif
  
  int nd = 0; // count of nodal domains
  int size; // area of last nodal domain (in pixels)

  i = 0;
  j = 0;
  while (findNextUnseen(counted, &i, &j, ny, nx)) {
    nd++;
    size = findDomain(signs, counted, i, j, nd, ny, nx);
    if (sizefile) {
      fprintf(sizefile, "%d ", size);
    }
  }

  #ifdef DEBUG
    intArray2file(counted, ny, nx, "../data/counted.dat");
  #endif
  
  return nd;
}

/*
count the number of nodal domains in grid
precondition: grid must be ny x nx

inputs:
        grid    - function values sampled at grid points
        counted - array to store nodal domain counts in. mask should already by applied to this
        nx      - number of samples in x-direction for grid
        ny      - number of samples in y-direction for grid
        sizefile - file to write size of nodal domains to
output: return value - count of nodal domains
*/
int countNodalDomainsNoInterp(double **grid, int **counted, int ny, int nx, FILE *sizefile) {
  int i, j;
  int rc;
  
  int nd = 0; // count of nodal domains
  int size; // area of last nodal domain (in pixels)

  i = 0;
  j = 0;
  while (findNextUnseen(counted, &i, &j, ny, nx)) {
    nd++;
    size = findDomainNoInterp(grid, counted, i, j, nd, ny, nx);
    if (sizefile) {
      fprintf(sizefile, "%d ", size);
    }
  }

  #ifdef DEBUG
  intArray2file(counted, ny, nx, "../data/counted.dat");
  #endif
  
  free_imatrix(counted);

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
         return value: area of domain (in pixels)
         updates stats
*/
int findDomain(bit_array_t *signs, int **counted, int i, int j, int nd, int ny, int nx) {
  stack *s = newStack();
  push(s, j, i);
  
  int x, y;
  int currentSign;
  char connections; // keep track of what (x,y) is connected to
  int size = 0;

  while (pop(s, &x, &y)) {
    size++;
    currentSign = bit_array_get(signs, j, i);

    // orthongal directions
    // left
    if (x >= 1 && !IS_MASKED(counted[y][x-1])) {
      if (bit_array_get(signs, x-1, y) == currentSign) {
	if(!IS_COUNTED(counted[y][x-1])) {
	  push(s, x - 1, y);
	}
      }
    }

    // above
    if (y >= 1 && !IS_MASKED(counted[y-1][x])) {
      if (bit_array_get(signs, x, y-1) == currentSign) {
	if(!IS_COUNTED(counted[y-1][x])) {
	  push(s, x, y-1);
	}
      }
    }

    // right
    if (x < nx-1 && !IS_MASKED(counted[y][x+1])) {
      if (bit_array_get(signs, x+1, y) == currentSign) {
	if(!IS_COUNTED(counted[y][x+1])) {
	  push(s, x+1, y);
	}
      }
    }

    // below
    if (y < ny-1 && !IS_MASKED(counted[y+1][x])) {
      if (bit_array_get(signs, x, y+1) == currentSign) {
	if(!IS_COUNTED(counted[y+1][x])) {
	  push(s, x, y+1);
	}
      }
    }

    counted[y][x] = currentSign ? nd : -nd;
  }
  destroyStack(s);
  return size;
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
         return value: area of domain (in pixels)
         updates stats
*/
int findDomainNoInterp(double **grid, int **counted, int i, int j, int nd, int ny, int nx) {
  stack *s = newStack();
  push(s, j, i);
  
  int x, y;
  int currentSign;
  char connections; // keep track of what (x,y) is connected to
  int size = 0;

  while (pop(s, &x, &y)) {
    size++;
    currentSign = SIGN(grid[i][j]);

    // orthongal directions
    // left
    if (x >= 1 && !IS_MASKED(grid[y][x-1])) {
      if (SIGN(grid[y][x-1]) == currentSign) {
        if(!IS_COUNTED(counted[y][x-1])) {
          push(s, x - 1, y);
        }
      }
    }

    // above
    if (y >= 1 && !IS_MASKED(grid[y-1][x])) {
      if (SIGN(grid[y-1][x]) == currentSign) {
        if(!IS_COUNTED(counted[y-1][x])) {
          push(s, x, y-1);
        }
      }
    }

    // right
    if (x < nx-1 && !IS_MASKED(grid[y][x+1])) {
      if (SIGN(grid[y][x+1]) == currentSign) {
        if(!IS_COUNTED(counted[y][x+1])) {
          push(s, x+1, y);
        }
      }
    }

    // below
    if (y < ny-1 && !IS_MASKED(grid[y+1][x])) {
      if (SIGN(grid[y+1][x]) == currentSign) {
        if(!IS_COUNTED(counted[y+1][x])) {
          push(s, x, y+1);
        }
      }
    }

    counted[y][x] = currentSign == 1 ? nd : -nd;
  }
  destroyStack(s);
  return size;
}

/*
  interpolate a region of the grid and put the results in counted
  
  precondition: grid, counted, and mask are ny x nx

  inputs:
        grid     - 2d array of function values
        counted  - 2d array of to apply mask to
	ny       - rows in grid
	nx       - columns in grid
        alpha    - k*dx
	upsample - upsampling rate
        M        - highest order bessel function to upsample with
	
  outputs:
        returns an ((ny-1)*upsample)+1 x ((nx-1)*upsample)+1 bit array of signs on upsampled grid
*/
 
bit_array_t *upsample(double **grid, int **counted, int ny, int nx, double alpha, int M, int upsample_ratio, interp_stats *stats) {
  gsl_matrix *interp = create_interp_matrix(alpha, M, upsample_ratio);
  bit_array_t *upsampled = new_bit_array((ny-1)*upsample_ratio+1, (nx-1)*upsample_ratio+1);
  MALLOC_CHECK(upsampled);
  void (*bit_array_update_fn)(bit_array_t *, int, int); // function pointer
  int r,c,x,y;
  interp_workspace *w = new_interp_workspace(upsample_ratio);

  for (r = 0 ; r < ny - 1 ; r++) {
    for (c = 0 ; c < nx - 1 ; c++) {
      if (SIGN(grid[r][c]) == SIGN(grid[r][c+1]) == SIGN(grid[r+1][c]) == SIGN(grid[r+1][c+1])) {
        bit_array_update_fn = grid[r][c] > 0 ? &bit_array_set : &bit_array_reset;
        for (y = r*upsample_ratio ; y <= (r+1)*upsample_ratio ; y++) {
          for (x = c*upsample_ratio ; x <= (c+1)*upsample_ratio ; x++) {
            (*bit_array_update_fn)(upsampled,x,y);
          }
        }
      } else {
        interpolate(grid, upsampled, r, c, ny, nx, upsample_ratio, interp, stats, w);
      }
    }
  }
  free_interp_workspace(w);
  return upsampled;
}



/*
  interpolate a region of the grid and put the results in counted.
  
  precondition: grid and counted are ny x nx

  inputs:
        grid      - 2d array of function values
	upsampled - 2d array of upsampled function values (to be filled with results)
	i         - row of initial point in grid
        j         - column of initial point in grid
	ny        - rows in grid
	nx        - columns in grid
	upsample  - upsampling rate
	interp    - precalculated 24 x (upsapmle+1)^2 matrix to do interpolation
	
  outputs:
        stores upsampled function value signs in upsampled
*/
#define IDX(x, y) ((x)*n+(y))
#define Y_REFLECT_IDX
#define X_REFLECT_IDX
void interpolate(double **grid, bit_array_t *upsampled, int i, int j, int ny, int nx, int upsample, gsl_matrix *interp, interp_stats *stats, interp_workspace *w) {
  int k;
  int x, y;
  int n = upsample + 1;
  int rc;
  int currentSign;
  
  // validate size of interp matrix
  if (interp->size1 != n*n || interp->size2 != NUM_STENCIL_POINTS) {
    ERROR("invalid size for interpolation matrix");
    exit(INTERP_ERR);
  }

  // check we're not too close to the edge
  // TODO: interpolate on edges by reflecting
  if (i < 2 || i >= ny - 3 || j < 2 || j >= nx - 3) {
    ERROR("trouble spot near edge: (x,y) = (%d,%d)", j, i);
    stats->edge_trouble_count++;
    return;
  }
  stats->interp_count++;

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
  gsl_vector_set(w->input, 0 , grid[Y_REFLECT_IDX(i-2)][X_REFLECT_IDX(j)]);
  gsl_vector_set(w->input, 1 , grid[Y_REFLECT_IDX(i-2)][X_REFLECT_IDX(j+1)]);
  gsl_vector_set(w->input, 2 , grid[Y_REFLECT_IDX(i-1)][X_REFLECT_IDX(j-1)]);
  gsl_vector_set(w->input, 3 , grid[Y_REFLECT_IDX(i-1)][X_REFLECT_IDX(j)]);
  gsl_vector_set(w->input, 4 , grid[Y_REFLECT_IDX(i-1)][X_REFLECT_IDX(j+1)]);
  gsl_vector_set(w->input, 5 , grid[Y_REFLECT_IDX(i-1)][X_REFLECT_IDX(j+2)]);
  gsl_vector_set(w->input, 6 , grid[Y_REFLECT_IDX(i)][X_REFLECT_IDX(j-2)]);
  gsl_vector_set(w->input, 7 , grid[Y_REFLECT_IDX(i)][X_REFLECT_IDX(j-1)]);
  gsl_vector_set(w->input, 8 , grid[Y_REFLECT_IDX(i)][X_REFLECT_IDX(j)]);
  gsl_vector_set(w->input, 9 , grid[Y_REFLECT_IDX(i)][X_REFLECT_IDX(j+1)]);
  gsl_vector_set(w->input, 10, grid[Y_REFLECT_IDX(i)][X_REFLECT_IDX(j+2)]);
  gsl_vector_set(w->input, 11, grid[Y_REFLECT_IDX(i)][X_REFLECT_IDX(j+3)]);
  gsl_vector_set(w->input, 12, grid[Y_REFLECT_IDX(i+1)][X_REFLECT_IDX(j-2)]);
  gsl_vector_set(w->input, 13, grid[Y_REFLECT_IDX(i+1)][X_REFLECT_IDX(j-1)]);
  gsl_vector_set(w->input, 14, grid[Y_REFLECT_IDX(i+1)][X_REFLECT_IDX(j)]);
  gsl_vector_set(w->input, 15, grid[Y_REFLECT_IDX(i+1)][X_REFLECT_IDX(j+1)]);
  gsl_vector_set(w->input, 16, grid[Y_REFLECT_IDX(i+1)][X_REFLECT_IDX(j+2)]);
  gsl_vector_set(w->input, 17, grid[Y_REFLECT_IDX(i+1)][X_REFLECT_IDX(j+3)]);
  gsl_vector_set(w->input, 18, grid[Y_REFLECT_IDX(i+2)][X_REFLECT_IDX(j-1)]);
  gsl_vector_set(w->input, 19, grid[Y_REFLECT_IDX(i+2)][X_REFLECT_IDX(j)]);
  gsl_vector_set(w->input, 20, grid[Y_REFLECT_IDX(i+2)][X_REFLECT_IDX(j+1)]);
  gsl_vector_set(w->input, 21, grid[Y_REFLECT_IDX(i+2)][X_REFLECT_IDX(j+2)]);
  gsl_vector_set(w->input, 22, grid[Y_REFLECT_IDX(i+3)][X_REFLECT_IDX(j)]);
  gsl_vector_set(w->input, 23, grid[Y_REFLECT_IDX(i+3)][X_REFLECT_IDX(j+1)]);
  
  // check if we are near the boundary
  // if we are near the boundary we can still interpolate
  //   be there is a continuation of the eigenfunction for ~1 wavelength outside the boundary
  for (k = 0 ; k < w->input->size ; k++) {
    if (IS_MASKED(gsl_vector_get(w->input, k))) {
      ERROR("trouble spot near boundary: (x,y) = (%d,%d)", j, i);
      stats->boundary_trouble_count++;
      break;
    }
  }

  // interp_output = interp * w->input
  rc = gsl_blas_dgemv(CblasNoTrans, 1, interp, w->input, 0, w->output);
  if (rc) {
    ERROR("interpolation failed. gsl_blas_dgemv returned %d", rc);
    exit(INTERP_ERR);
  }

  // write signs into upsampled
  for (y = 0 ; y <= upsample ; y++) {
    for (x = 0 ; x <= upsample ; x++) {
      if (gsl_vector_get(w->output, IDX(x,y)) > 0) {
        bit_array_set(upsampled, j*upsample + x, i*upsample + y);
      } else {
        bit_array_reset(upsampled, j*upsample + x, i*upsample + y);
      }
    }
  }
}

interp_workspace *new_interp_workspace(int upsample_ratio) {
  interp_workspace *w = (interp_workspace *)malloc(sizeof(interp_workspace));
  MALLOC_CHECK(w);
  w->input = gsl_vector_alloc(NUM_STENCIL_POINTS);
  MALLOC_CHECK(w->input);
  w->output = gsl_vector_alloc((upsample_ratio+1)*(upsample_ratio+1));
  MALLOC_CHECK(w->output);
  return w;
}

void free_interp_workspace(interp_workspace *w) {
  free(w->input);
  free(w->output);
  free(w);
}
