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

#define SIGN(x) (((x)>0)-((x)<0))

#ifdef DEBUG
int **domain_numbers;
#endif

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
int countNodalDomainsInterp(double **grid, bit_array_t *counted, int ny, int nx, double alpha, int M, int upsample_ratio, interp_stats *stats, FILE *sizefile) {
  int nodal_domain_count;
  bit_array_t *upsampled_grid;
  
  upsampled_grid = upsample(grid, ny, nx, alpha, M, upsample_ratio, stats);
  // TODO(KK): It seems like we need to "upsample" the count array also...
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
int countNodalDomains(bit_array_t *signs, bit_array_t *counted, FILE *sizefile) {
  int i, j;
  int rc;
  int nx, ny;
  nx = signs->nx;
  ny = signs->ny;

#ifdef DEBUG
  bit_array2file(signs, "../data/signs.dat");
  domain_numbers = imatrix(ny, nx);
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
  intArray2file(domain_numbers, ny, nx, "../data/counted.dat");
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
int countNodalDomainsNoInterp(double **grid, bit_array_t *counted, int ny, int nx, FILE *sizefile) {
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
int findNextUnseen(bit_array_t *counted, int *i, int *j, int ny, int nx) {
  int r, c;

  for (r = *i ; r < ny ; r++) {
    for (c = 0 ; c < nx ; c++) {
      if (r == *i && c <= *j)
        continue;
      if (!bit_array_get(counted, c, r)) {
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
int findDomain(bit_array_t *signs, bit_array_t *counted, int i, int j, int nd, int ny, int nx) {
  stack *s = newStack();
  push(s, j, i);
  bit_array_set(counted, j, i);

  int x, y;
  int currentSign;
  int size = 0;

  while (pop(s, &x, &y)) {
    size++;
    currentSign = bit_array_get(signs, j, i);
    
#ifdef DEBUG
    domain_numbers[y][x] = nd;
#endif

    // orthongal directions
    // left
    if (x >= 1 && !bit_array_get(counted, x-1, y)) {
      if (bit_array_get(signs, x-1, y) == currentSign) {
        bit_array_set(counted, x-1, y);
        push(s, x - 1, y);
      }
    }

    // above
    if (y >= 1 && !bit_array_get(counted, x, y-1)) {
      if (bit_array_get(signs, x, y-1) == currentSign) {
        bit_array_set(counted, x, y-1);
        push(s, x, y-1);
      }
    }

    // right
    if (x < nx-1 && !bit_array_get(counted, x+1, y)) {
      if (bit_array_get(signs, x+1, y) == currentSign) {
        bit_array_set(counted, x+1, y);
        push(s, x+1, y);
      }
    }

    // below
    if (y < ny-1 && !bit_array_get(counted, x, y+1)) {
      if (bit_array_get(signs, x, y+1) == currentSign) {
        bit_array_set(counted, x, y+1);
        push(s, x, y+1);
      }
    }
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
int findDomainNoInterp(double **grid, bit_array_t *counted, int i, int j, int nd, int ny, int nx) {
  stack *s = newStack();
  push(s, j, i);
  bit_array_set(counted, j, i);

  int x, y;
  int currentSign;
  int size = 0;

  while (pop(s, &x, &y)) {
    size++;
    currentSign = SIGN(grid[i][j]);

    // orthongal directions
    // left
    if (x >= 1 && !bit_array_get(counted, x-1, y)) {
      if (SIGN(grid[y][x-1]) == currentSign) {
        bit_array_set(counted, x-1, y);
        push(s, x - 1, y);
      }
    }

    // above
    if (y >= 1 && !bit_array_get(counted, x, y-1)) {
      if (SIGN(grid[y-1][x]) == currentSign) {
        bit_array_set(counted, x, y-1);
        push(s, x, y-1);
      }
    }

    // right
    if (x < nx-1 && !bit_array_get(counted, x+1, y)) {
      if (SIGN(grid[y][x+1]) == currentSign) {
        bit_array_set(counted, x+1, y);
        push(s, x+1, y);
      }
    }

    // below
    if (y < ny-1 && !bit_array_get(counted, x, y+1)) {
      if (SIGN(grid[y+1][x]) == currentSign) {
        bit_array_set(counted, x, y+1);
        push(s, x, y+1);
      }
    }
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
 
bit_array_t *upsample(double **grid, int ny, int nx, double alpha, int M, int upsample_ratio, interp_stats *stats) {
  gsl_matrix *interp = create_interp_matrix(alpha, M, upsample_ratio);
  bit_array_t *upsampled = new_bit_array((ny-1)*upsample_ratio+1, (nx-1)*upsample_ratio+1);
  MALLOC_CHECK(upsampled);
  int r,c,x,y;
  interp_workspace *w = new_interp_workspace(upsample_ratio, stats);
  for (r = 0 ; r < ny - 3 ; r++) {
    for (c = 0 ; c < nx - 3 ; c++) {
      interpolate(grid, upsampled, r, c, ny, nx, upsample_ratio, interp, w);
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
void interpolate(double **grid, bit_array_t *upsampled, int i, int j, int ny, int nx, int upsample_ratio, gsl_matrix *interp, interp_workspace *w) {
  int k;
  int x, y, l;
  int n = upsample_ratio + 1;
  int rc;
  int currentSign;
  int refl_i[STENCIL_WIDTH], refl_j[STENCIL_WIDTH], refl_i_sign[STENCIL_WIDTH], refl_j_sign[STENCIL_WIDTH];

  // validate size of interp matrix
  if (interp->size1 != n*n || interp->size2 != NUM_STENCIL_POINTS) {
    ERROR("invalid size for interpolation matrix");
    exit(INTERP_ERR);
  }

  w->stats->interp_count++;

  fill_interp_input(grid, i, j, w, ny, nx, refl_i, refl_j, refl_i_sign, refl_j_sign);
  
  // interp_output = interp * w->input
  rc = gsl_blas_dgemv(CblasNoTrans, 1, interp, w->input, 0, w->output);
  if (rc) {
    ERROR("interpolation failed. gsl_blas_dgemv returned %d", rc);
    exit(INTERP_ERR);
  }
  
  // write signs into upsampled
  for (y = 0 ; y < n ; y++) {
    for (x = 0 ; x < n ; x++) {
      if (gsl_vector_get(w->output, IDX(x,y)) > 0) {
        bit_array_set(upsampled, j*upsample_ratio + x, i*upsample_ratio + y);
      } else {
        bit_array_reset(upsampled, j*upsample_ratio + x, i*upsample_ratio + y);
      }
    }
  }

  // debug output
  /*
    char interp_input_file[100];
    char interp_output_file[100];
    sprintf(interp_input_file, "interp_input_%d_%d.dat", j, i);
    sprintf(interp_output_file, "interp_output_%d_%d.dat", j, i);
    dump_vector(w->input, interp_input_file);
    dump_vector(w->output, interp_output_file);
  */
}

/*
  fill w->interp_input with grid values around (x,y)
  reflect eigenfunction across edge of data
  interpolate across domain boundary
  because there is a continuation of the eigenfunction for ~1 wavelength outside the boundary

  precondition: each int * argument is an array of length STENCIL_WIDTH
*/
// TODO: move all these params into interp_workspace struct
//         also move interp_stats to workspace
void fill_interp_input(double **grid, int i, int j, interp_workspace *w, int ny, int nx, int *refl_i, int *refl_j, int *refl_i_sign, int *refl_j_sign) {
  int l;

  if (i < 2 || j < 2) { // reflection case
    w->stats->edge_trouble_count++;
    for (l = 0 ; l < STENCIL_WIDTH ; l++) {
      refl_i_sign[l] = 1;
      refl_i[l] = i - 2 + l;
      // across x-axis
      if (refl_i[l] < 0) {
        refl_i_sign[l] *= -1;
        refl_i[l] *= -1;
      }
      refl_j_sign[l] = 1;
      refl_j[l] = j - 2 + l;
      // across y-axis
      if (refl_j[l] < 0) {
        refl_j_sign[l] = -1;
        refl_j[l] *= -1;
      }
    }
    gsl_vector_set(w->input, 0 , refl_i_sign[0]*refl_j_sign[2]*grid[refl_i[0]][refl_j[2]]);
    gsl_vector_set(w->input, 1 , refl_i_sign[0]*refl_j_sign[3]*grid[refl_i[0]][refl_j[3]]);
    gsl_vector_set(w->input, 2 , refl_i_sign[1]*refl_j_sign[1]*grid[refl_i[1]][refl_j[1]]);
    gsl_vector_set(w->input, 3 , refl_i_sign[1]*refl_j_sign[2]*grid[refl_i[1]][refl_j[2]]);
    gsl_vector_set(w->input, 4 , refl_i_sign[1]*refl_j_sign[3]*grid[refl_i[1]][refl_j[3]]);
    gsl_vector_set(w->input, 5 , refl_i_sign[1]*refl_j_sign[4]*grid[refl_i[1]][refl_j[4]]);
    gsl_vector_set(w->input, 6 , refl_i_sign[2]*refl_j_sign[0]*grid[refl_i[2]][refl_j[0]]);
    gsl_vector_set(w->input, 7 , refl_i_sign[2]*refl_j_sign[1]*grid[refl_i[2]][refl_j[1]]);
    gsl_vector_set(w->input, 8 , refl_i_sign[2]*refl_j_sign[2]*grid[refl_i[2]][refl_j[2]]);
    gsl_vector_set(w->input, 9 , refl_i_sign[2]*refl_j_sign[3]*grid[refl_i[2]][refl_j[3]]);
    gsl_vector_set(w->input, 10, refl_i_sign[2]*refl_j_sign[4]*grid[refl_i[2]][refl_j[4]]);
    gsl_vector_set(w->input, 11, refl_i_sign[2]*refl_j_sign[5]*grid[refl_i[2]][refl_j[5]]);
    gsl_vector_set(w->input, 12, refl_i_sign[3]*refl_j_sign[0]*grid[refl_i[3]][refl_j[0]]);
    gsl_vector_set(w->input, 13, refl_i_sign[3]*refl_j_sign[1]*grid[refl_i[3]][refl_j[1]]);
    gsl_vector_set(w->input, 14, refl_i_sign[3]*refl_j_sign[2]*grid[refl_i[3]][refl_j[2]]);
    gsl_vector_set(w->input, 15, refl_i_sign[3]*refl_j_sign[3]*grid[refl_i[3]][refl_j[3]]);
    gsl_vector_set(w->input, 16, refl_i_sign[3]*refl_j_sign[4]*grid[refl_i[3]][refl_j[4]]);
    gsl_vector_set(w->input, 17, refl_i_sign[3]*refl_j_sign[5]*grid[refl_i[3]][refl_j[5]]);
    gsl_vector_set(w->input, 18, refl_i_sign[4]*refl_j_sign[1]*grid[refl_i[4]][refl_j[1]]);
    gsl_vector_set(w->input, 19, refl_i_sign[4]*refl_j_sign[2]*grid[refl_i[4]][refl_j[2]]);
    gsl_vector_set(w->input, 20, refl_i_sign[4]*refl_j_sign[3]*grid[refl_i[4]][refl_j[3]]);
    gsl_vector_set(w->input, 21, refl_i_sign[4]*refl_j_sign[4]*grid[refl_i[4]][refl_j[4]]);
    gsl_vector_set(w->input, 22, refl_i_sign[5]*refl_j_sign[2]*grid[refl_i[5]][refl_j[2]]);
    gsl_vector_set(w->input, 23, refl_i_sign[5]*refl_j_sign[3]*grid[refl_i[5]][refl_j[3]]);
  } else {

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


    gsl_vector_set(w->input, 0 , grid[i-2][j]);
    gsl_vector_set(w->input, 1 , grid[i-2][j+1]);
    gsl_vector_set(w->input, 2 , grid[i-1][j-1]);
    gsl_vector_set(w->input, 3 , grid[i-1][j]);
    gsl_vector_set(w->input, 4 , grid[i-1][j+1]);
    gsl_vector_set(w->input, 5 , grid[i-1][j+2]);
    gsl_vector_set(w->input, 6 , grid[i][j-2]);
    gsl_vector_set(w->input, 7 , grid[i][j-1]);
    gsl_vector_set(w->input, 8 , grid[i][j]);
    gsl_vector_set(w->input, 9 , grid[i][j+1]);
    gsl_vector_set(w->input, 10, grid[i][j+2]);
    gsl_vector_set(w->input, 11, grid[i][j+3]);
    gsl_vector_set(w->input, 12, grid[i+1][j-2]);
    gsl_vector_set(w->input, 13, grid[i+1][j-1]);
    gsl_vector_set(w->input, 14, grid[i+1][j]);
    gsl_vector_set(w->input, 15, grid[i+1][j+1]);
    gsl_vector_set(w->input, 16, grid[i+1][j+2]);
    gsl_vector_set(w->input, 17, grid[i+1][j+3]);
    gsl_vector_set(w->input, 18, grid[i+2][j-1]);
    gsl_vector_set(w->input, 19, grid[i+2][j]);
    gsl_vector_set(w->input, 20, grid[i+2][j+1]);
    gsl_vector_set(w->input, 21, grid[i+2][j+2]);
    gsl_vector_set(w->input, 22, grid[i+3][j]);
    gsl_vector_set(w->input, 23, grid[i+3][j+1]);
  }
}
    
// TODO: move to util
interp_workspace *new_interp_workspace(int upsample_ratio, interp_stats *stats) {
  interp_workspace *w = (interp_workspace *)malloc(sizeof(interp_workspace));
  MALLOC_CHECK(w);
  w->input = gsl_vector_alloc(NUM_STENCIL_POINTS);
  MALLOC_CHECK(w->input);
  w->output = gsl_vector_alloc((upsample_ratio+1)*(upsample_ratio+1));
  MALLOC_CHECK(w->output);
  w->stats = stats;
  return w;
}

void free_interp_workspace(interp_workspace *w) {
  free(w->input);
  free(w->output);
  free(w);
}
