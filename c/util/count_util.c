/*
  Utility functions for nodal_domain_count.c and nodal_domain_driver.c

  Kyle Konrad
  4/4/2011
*/

#include "count_util.h"
#include "util.h"
#include "exit_codes.h"
#include <gsl/gsl_matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <wait.h>

extern int verb;

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

  double **grid = dmatrix(ny, nx);
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

  double **grid = dmatrix(ny, nx);
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
  calculate wing tip mass of a stadium eigenfunction

  precondition: grid is a quarter stadium with width 2 (qust:2 in verg)

  inputs:
          grid - ny x nx array
	  ny   - number of rows in arrays
	  nx   - number of columns in arrays

 */
#define IN_QUST(x,y) (((x)-1)*((x)-1)+(y)*(y)<1)
double wingTipMass(double **grid, int ny, int nx) {
  double wtm = 0.0;
  double dx = 2.0 / nx;

  int i, j;
  for (i = 0 ; i < ny ; i++) {
    for (j = 0 ; j < nx ; j++) {
      if (j*dx >= 1.1 && IN_QUST(j*dx, i*dx)) {
	wtm += grid[i][j]*grid[i][j];
      }
    }
  }
  wtm *= dx*dx;
  return wtm;
}
