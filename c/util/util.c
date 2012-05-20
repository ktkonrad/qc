/*
general purpose utility functions
*/

#include <stdlib.h>
#include <stdio.h>
#include "util.h"
#include "exit_codes.h"

extern int verb;

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


/*
malloc an ny x nx array of doubles
input:
       ny   - y dimension of array
       nx   - x dimension of array
*/
double **dmatrix(int ny, int nx) {
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
char **cmatrix(int ny, int nx) {
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
void free_dmatrix(double **grid) {
  free(grid[0]);
  free(grid);
}

/*
  free memory used by a mask

  input:
         mask - array to be freed
*/
void free_cmatrix(char **mask) {
  free(mask[0]);
  free(mask);
}

void free_imatrix(int **m) {
  free(m[0]);
  free(m);
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
