/*
  Utility functions for nodal_domain_count.c and nodal_domain_driver.c

  Kyle Konrad
  4/4/2011
*/

#include "util.h"
#include <gsl/gsl_matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <wait.h>
#include <limits.h>

/*
  create a matrix to do interpolation based on bessel functions
  
  input:
         k        - wavenumber of eigenfunction being interpolated
         dx       - sampled resoultion of eigenfunction
         M        - highest order bessel function to do
	 upsample - factor to upsample by
         m        - matrix to put values in

  output:
          return value - status (0 for success)
	  m            - interpolation matrix

  PRECONDITON: m is (upsample+1)^2 x 24
*/
// DEPRECATED
int fillInterpMatrix(double k, double dx, int M, int upsample, gsl_matrix *m) {
  char interpfile[100];
  char *executable = "../matlab/create_interp_matrix.sh"; // hacky way to allow calling this from other directories
  char *args[7];
  int pid;
  int rc;
  int i;

  sprintf(interpfile, "../data/interp_matrix_k=%0.4f_dx=%0.4f_M=%d.dat", k, dx, M);
  
  if (access(interpfile, R_OK) == -1) { // only create the matrix if we don't already have it
    pid = vfork();
    if (pid == 0) { // child
      for (i = 1 ; i < 5 ; i++) {
	args[i] = (char *)malloc(100*sizeof(char));
      }
      args[0] = executable;
      sprintf(args[1], "%f", k);
      sprintf(args[2], "%f", dx);
      sprintf(args[3], "%d", M); 
      sprintf(args[4], "%d", upsample);
      args[5] = interpfile;
      args[6] = NULL;
      execv(executable, args);
      ERROR("execv failed");
      for (i = 1 ; i < 5 ; i++) {
	free(args[i]);
      }    
      return INTERP_INIT_ERR;
    }
    // parent
    wait(&rc);
    if (rc != 0) {
      ERROR("failed to create interpolation matrix. child process returned %d", rc);
      return INTERP_INIT_ERR;
    }
  }

  FILE *interp_matrix_file = fopen(interpfile, "r");
  rc = gsl_matrix_fscanf(interp_matrix_file, m);
  fclose(interp_matrix_file);
  if (rc) {
    ERROR("gsl_matrix_fscanf failed. return code was %d", rc);
    return IO_ERR;
  }
  
  return 0;
}

/*
  read a sta_bin file into grid
  only reads middle eigenfunction in file
  Adapted from Alex Barnett's viewer.c

  input:
         grid - array to read data into
         file - name of file to read
	 m - where to put number of rows in grid
	 n - where to put number of columns in grid

  output:
          return value - grid (NULL if something failed)
	  m            - number of rows in grid
	  n            - number of columns in grid
*/
double **readOneSta(char *file, int *m, int *n) {
  FILE *fp;

  if (strcmp(file, "stdin") == 0)
    fp = stdin;
  else {
    fp = fopen(file, "r");
    if (fp == NULL) {
      ERROR("failed to open %s", file);
      return NULL;
    }
  }

  char c;
  int n_e, nx, ny;
  double temp_double;
  float temp_float;

  fscanf(fp,"%c",&c);
  if (c == 'b') {
    fseek(fp, 3, SEEK_CUR);
  }
  else {
    if (fp != stdin)
      fclose(fp);
    ERROR("incorrect file format in %s", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("failed to read n_e in %s", file);
    return NULL;
  }

  if (fread(&nx, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("failed to read nx in %s", file);
    return NULL;
  }

  if (fread(&ny, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("failed to read ny in %s", file);
    return NULL;
  } 

  if (n_e < 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("no eigenfunctions in %s", file);
    return NULL;
  }

  double **grid = createGrid(ny, nx);
  *m = ny;
  *n = nx;

  int i;

  for (i = 0 ; i < n_e ; i++) {
    if (fread(&temp_double,sizeof(double),1,fp) != 1) { // this is the energy - we don't care about it here
      if (fp != stdin)
	fclose(fp);
      ERROR("failed to read E_1 in %s", file);
      return NULL;
    }
  }

  int j;

  for (i = 0 ; i < nx * ny ; i++) {
    for (j = 0 ; j < n_e ; j++) {
      if (fread(&temp_float, 4, 1, fp) != 1) {
	if (fp != stdin) {
	  fclose(fp);
	}
	ERROR("failed to read data in %s", file);
	return NULL;
      }
      
      if (j == n_e / 2) { //only read 1 eigenfunction
	grid[i/nx][i%nx] = (double)temp_float;
      }
    }
  }

  if (fp != stdin)
    fclose(fp);
  return grid;
}

/*
  read a sta_bin file into grid
  Adapted from Alex Barnett's viewer.c

  input:
         grid - array to read data into
         file - name of file to read
	 ne   - where to put number of eigenfunctions in file
	 m    - where to put number of rows in grid
	 n    - where to put number of columns in grid
	 k    - where to put k of data read
	 l    - read the lth eigenfunction (first is l=0)

  output:
          return value - grid (NULL if something failed)
	  m            - number of rows in grid
	  n            - number of columns in grid
	  k            - wavenumber of eigenfunction
*/
double **readSta(char *file, int *ne, int *m, int *n, double *k, int l) {
  FILE *fp;

  if (strcmp(file, "stdin") == 0)
    fp = stdin;
  else {
    fp = fopen(file, "r");
    if (fp == NULL) {
      ERROR("failed to open %s", file);
      return NULL;
    }
  }

  char c;
  int n_e, nx, ny;
  double temp_double;
  float temp_float;

  fscanf(fp,"%c",&c);
  if (c == 'b') {
    fseek(fp, 3, SEEK_CUR);
  }
  else {
    if (fp != stdin)
      fclose(fp);
    ERROR("incorrect file format in %s", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("failed to read n_e in %s", file);
    return NULL;
  }

  if (l < 0 || l >= n_e) {
    if (fp != stdin)
      fclose(fp);
    ERROR("invalid eigenfunction number: %d", l);
    return NULL;
  }

  if (fread(&nx, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("failed to read nx in %s", file);
    return NULL;
  }

  if (fread(&ny, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("failed to read ny in %s", file);
    return NULL;
  } 

  if (n_e < 1) {
    if (fp != stdin)
      fclose(fp);
    ERROR("no eigenfunctions in %s", file);
    return NULL;
  }

  double **grid = createGrid(ny, nx);
  *ne = n_e;
  *m = ny;
  *n = nx;

  int i;
  // TODO: minor efficiency improvement: use fseek instead of reading everything
  for (i = 0 ; i < n_e ; i++) {
    if (fread(&temp_double,sizeof(double),1,fp) != 1) {
      if (fp != stdin)
	fclose(fp); 
      ERROR("failed to read k_%d in %s", i, file);
      return NULL;
    }
    if (i == l)
      *k = temp_double;
  }

  int j;
  // TODO: minor efficiency improvement: use fseek instead of reading everything
  for (i = 0 ; i < nx * ny ; i++) {
    for (j = 0 ; j < n_e ; j++) {
      if (fread(&temp_float, 4, 1, fp) != 1) {
	if (fp != stdin)
	  fclose(fp);
	ERROR("failed to read data in %s", file);
	return NULL;
      }
      
      if (j == l) { //only read lth eigenfunction
	grid[i/nx][i%nx] = (double)temp_float;
      }
    }
  }

  if (fp != stdin)
    fclose(fp);
  return grid;
}


/*
  read a mask from sta_bin file
  currently only reads single eigenfunction files

  input:
         grid - array to read data into
         file - name of file to read
	 m - where to put number of rows in grid
	 n - where to put number of columns in grid

  output:
          return value - grid (NULL if something failed)
	  m            - number of rows in grid
	  n            - number of columns in grid
*/
char **readMask(char *file, int *m, int *n) {
  FILE *fp = fopen(file, "r");
  if (fp == NULL) {
    ERROR("failed to open %s", file);
    return NULL;
  }

  char c;
  int n_e, nx, ny;
  double temp_double;
  float temp_float;

  fscanf(fp,"%c",&c);
  if (c == 'b') {
    fscanf(fp,"%c",&c);
    fscanf(fp,"%c",&c);
    fscanf(fp,"%c",&c);
  }
  else {
    fclose(fp);
    ERROR("incorrect file format in %s", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    ERROR("failed to read n_e in %s", file);
    return NULL;
  }

  if (fread(&nx, (size_t)4,1,fp) != 1) {
    ERROR("failed to read nx in %s", file);
    return NULL;
  }

  if (fread(&ny, (size_t)4,1,fp) != 1) {
    ERROR("failed to read ny in %s", file);
    return NULL;
  } 

  if (n_e != 1) {
    ERROR("more than one eigenfunction in %s", file);
    return NULL;
  }

  char **mask = createMask(ny, nx);
  *m = ny;
  *n = nx;

  if (fread(&temp_double,sizeof(double),1,fp) != 1) { // this is the energy - we don't care about it here
    ERROR("failed to read E_1 in %s", file);
    return NULL;
  }

  int i;
  for (i = 0 ; i < nx * ny ; i++) {
    if (fread(&temp_float,4,n_e,fp) != (unsigned int)n_e) {
      ERROR("failed to read data in %s", file);
      return NULL;
    }
    mask[i/nx][i%nx] = (char)temp_float;
  }

  fclose(fp);
  return mask;
}

/*
malloc an ny x nx array of doubles
input:
       ny   - y dimension of array
       nx   - x dimension of array
*/
double **createGrid(int ny, int nx) {
  double **grid;
  int i;
  grid = (double **)calloc(ny, sizeof(double *));
  if (!grid) {
    ERROR("failed to allocate matrix");
    return NULL;
  }
  grid[0] = (double *)calloc(ny*nx, sizeof(double));
  if (!grid[0]) {
    ERROR("failed to allocate matrix");
    return NULL;
  }
  for(i = 1 ; i < ny ; i++) {
    grid[i] = grid[i-1] + nx;
  }
  return grid;
}

/*
malloc an ny x nx array of chars
input:
       ny   - y dimension of array
       nx   - x dimension of array
output:
       returns - uninitialized 2d array
*/
char **createMask(int ny, int nx) {
  char **mask;
  int i;
  mask = (char **)calloc(ny, sizeof(char *));
  if (!mask) {
    ERROR("failed to allocate matrix");
    return NULL;
  }
  mask[0] = (char *)calloc(ny*nx, sizeof(char));
  if (!mask[0]) {
    ERROR("failed to allocate matrix");
    return NULL;
  }
  for(i = 1 ; i < ny ; i++) {
    mask[i] = mask[i-1] + nx;
  }
  return mask;
}

/*
  free memory used by a grid

  input:
         grid - array to be freed
         ny - y dimension of grid
*/
void destroyGrid(double **grid) {
  free(grid[0]);
  free(grid);
}

/*
  free memory used by a mask

  input:
         mask - array to be freed
*/
void destroyMask(char **mask) {
  free(mask[0]);
  free(mask);
}


/*
  write an array of doubles to a file, tab separator

  input:
         array - array to write
         m     - rows in array
	 n     - columns in array
	 file  - name of file to write to
  output: status code
*/
int array2file(double **array, int m, int n, char *file) {
  int i, j;
  FILE *out = fopen(file, "w");
  if (out == NULL) {
    ERROR("failed to open %s", file);
    return IO_ERR;
  }
  for (i = 0 ; i < m ; i++) {
    for (j = 0 ; j < n ; j++) {
      if (j != 0)
	fprintf(out, "\t");
      fprintf(out, "%f", array[i][j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}

/*
  write an array of ints to a file, tab separator

  input:
         array - array to write
         m     - rows in array
	 n     - columns in array
	 file  - name of file to write to

  output: status code
*/
int intArray2file(int **array, int m, int n, char *file) {
  int i, j;
  FILE *out = fopen(file, "w");
  if (out == NULL) {
    ERROR("failed to open %s", file);
    return IO_ERR;
  }
  for (i = 0 ; i < m ; i++) {
    for (j = 0 ; j < n ; j++) {
      if (j != 0)
	fprintf(out, "\t");
      fprintf(out, "%d", array[i][j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}

/*
  write an array of chars to a file, tab separator

  input:
         array - array to write
         m     - rows in array
	 n     - columns in array
	 file  - name of file to write to

  output: status code
*/
int charArray2file(char **array, int m, int n, char *file) {
  int i, j;
  FILE *out = fopen(file, "w");
  if (out == NULL) {
    ERROR("failed to open %s", file);
    return IO_ERR;
  }
  for (i = 0 ; i < m ; i++) {
    for (j = 0 ; j < n ; j++) {
      if (j != 0)
	fprintf(out, "\t");
      fprintf(out, "%d", (int)array[i][j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}

/*
  apply a mask to a grid and counted by setting all grid to nan and counted tp -1 where mask is 0
  
  inputs:
          grid    - ny x nx array
	  counted - ny x nx array
	  mask    - ny x nx array
	  ny      - number of rows in arrays
          nx      - number of columns in arrays
*/
void applyMask(double **grid, int **counted, char **mask, int ny, int nx) {
  int i, j;
  for (i = 0 ; i < ny ; i++) {
    grid[0][i] = INT_MAX;
    grid[i][0] = INT_MAX;
    counted[0][i] = 0;
    counted[i][0] = 0;
    for (j = 0 ; j < nx ; j++) {
      if (!mask[i][j]) {
	grid[i][j] = INT_MIN;
	counted[i][j] = 0;
      }
    }
  }
}

/*
  calculate wing tip mass of a stadium eigenfunction

  precondition: grid is a quarter stadium with width 2 (qust:2 in verg)

  inputs:
          grid - ny x nx array
	  mask - ny x nx array
	  ny   - number of rows in arrays
	  nx   - number of columns in arrays

 */
double wingTipMass(double **grid, char **mask, int ny, int nx) {
  double wtm = 0.0;
  double dx = 2.0 / nx;

  int i, j;
  for (i = 0 ; i < ny ; i++)
    for (j = 0 ; j < nx ; j++)
      if (j*dx >= 1.1 && mask[i][j])
	wtm += dx*dx * grid[i][j]*grid[i][j];

  return wtm;
}

/*
  allocate a rows x cols matrix of ints
  matrix is initialized to 0
  code adapted from numerical recipes utilites

  inputs:
        rows - number of rows
	cols - number of cols
*/
int **imatrix(int rows, int cols) {
  int **m;
  int i;
  m = (int **)calloc(rows, sizeof(int *));
  if (!m) {
    ERROR("failed to allocate matrix");
    return NULL;
  }
  m[0] = (int *)calloc(rows*cols, sizeof(int));
  if (!m[0]) {
    ERROR("failed to allocate matrix");
    return NULL;
  }
  for(i = 1 ; i < rows ; i++) {
    m[i] = m[i-1] + cols;
  }
  return m;
}

void free_imatrix(int **m) {
  free(m[0]);
  free(m);
}
