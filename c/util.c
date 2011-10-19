/*
  Utility functions for nodal_domain_count.c and nodal_domain_driver.c

  Kyle Konrad
  4/4/2011
*/

#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

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
      fprintf(stderr, "readOneSta: failed to open %s\n", file);
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
    fprintf(stderr, "readOneSta: incorrect file format in %s\n", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readOneSta: failed to read n_e in %s\n", file);
    return NULL;
  }

  if (fread(&nx, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readOneSta: failed to read nx in %s\n", file);
    return NULL;
  }

  if (fread(&ny, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readOneSta: failed to read ny in %s\n", file);
    return NULL;
  } 

  if (n_e < 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readOneSta: no eigenfunctions in %s\n", file);
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
      fprintf(stderr, "readOneSta: failed to read E_1 in %s\n", file);
      return NULL;
    }
  }

  int j;

  for (i = 0 ; i < nx * ny ; i++) {
    for (j = 0 ; j < n_e ; j++) {
      if (fread(&temp_float, 4, 1, fp) != 1) {
	if (fp != stdin)
	  fclose(fp);
	fprintf(stderr, "readOneSta: failed to read data in %s\n", file);
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
      fprintf(stderr, "readSta: failed to open %s\n", file);
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
    fprintf(stderr, "readSta: incorrect file format in %s\n", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readSta: failed to read n_e in %s\n", file);
    return NULL;
  }

  if (l < 0 || l >= n_e) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readSta: invalid eigenfunction number: %d\n", l);
    return NULL;
  }

  if (fread(&nx, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readSta: failed to read nx in %s\n", file);
    return NULL;
  }

  if (fread(&ny, (size_t)4,1,fp) != 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readSta: failed to read ny in %s\n", file);
    return NULL;
  } 

  if (n_e < 1) {
    if (fp != stdin)
      fclose(fp);
    fprintf(stderr, "readSta: no eigenfunctions in %s\n", file);
    return NULL;
  }

  double **grid = createGrid(ny, nx);
  *ne = n_e;
  *m = ny;
  *n = nx;

  int i;

  for (i = 0 ; i < n_e ; i++) {
    if (fread(&temp_double,sizeof(double),1,fp) != 1) {
      if (fp != stdin)
	fclose(fp); 
      fprintf(stderr, "readSta: failed to read k_%d in %s\n", i, file);
      return NULL;
    }
    if (i == l)
      *k = temp_double;
  }

  int j;

  for (i = 0 ; i < nx * ny ; i++) {
    for (j = 0 ; j < n_e ; j++) {
      if (fread(&temp_float, 4, 1, fp) != 1) {
	if (fp != stdin)
	  fclose(fp);
	fprintf(stderr, "readSta: failed to read data in %s\n", file);
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
    fprintf(stderr, "readSta: failed to open %s\n", file);
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
    fprintf(stderr, "readMask: incorrect file format in %s\n", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    fprintf(stderr, "readMask: failed to read n_e in %s\n", file);
    return NULL;
  }

  if (fread(&nx, (size_t)4,1,fp) != 1) {
    fprintf(stderr, "readMask: failed to read nx in %s\n", file);
    return NULL;
  }

  if (fread(&ny, (size_t)4,1,fp) != 1) {
    fprintf(stderr, "readMask: failed to read ny in %s\n", file);
    return NULL;
  } 

  if (n_e != 1) {
    fprintf(stderr, "readMask: more than one eigenfunction in %s\n", file);
    return NULL;
  }

  char **mask = createMask(ny, nx);
  *m = ny;
  *n = nx;

  if (fread(&temp_double,sizeof(double),1,fp) != 1) { // this is the energy - we don't care about it here
    fprintf(stderr, "readMask: failed to read E_1 in %s\n", file);
    return NULL;
  }

  int i;
  for (i = 0 ; i < nx * ny ; i++) {
    if (fread(&temp_float,4,n_e,fp) != (unsigned int)n_e) {
      fprintf(stderr, "readMask: failed to read data in %s\n", file);
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

  double **grid = (double **)malloc(ny * sizeof(double *));
  int i;
  for (i = 0 ; i < ny ; i++)
    grid[i] = (double *)malloc(nx * sizeof(double));

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

  char **mask = (char **)malloc(ny * sizeof(char *));
  int i;
  for (i = 0 ; i < ny ; i++)
    mask[i] = (char *)malloc(nx * sizeof(char));

  return mask;
}

/*
create a mask for a billiard
input:
       b    - billiard to create mask for
       dx   - grid spacing
output:
       returns - boolean mask array
       ny      - rows in array
       nx      - columns in array
*/
char **createMaskFromBilliard(Billiard b, double dx, int *ny, int *nx) {
  *ny = ceil((b.yh - b.yl) / dx) + 1;
  *nx = ceil((b.xh - b.xl) / dx) + 1;
  char **mask = createMask(*ny, *nx);

  int i, j;
  for (int i = 0 ; i < *ny ; i++)
    for (int j = 0 ; j < *nx ; j++)
      mask[i][j] = inside_billiard(j * dx, i * dx, &b);

  return mask;
}

/*
create a mask for a billiard, scale it
input:
       b     - billiard to create mask for
       dx    - grid spacing
       scale - factor that grid is scaled by (k/k_0)
output:
       returns - boolean mask array
       ny      - rows in array
       nx      - columns in array
*/
char **createScaledMaskFromBilliard(Billiard b, double dx, int *ny, int *nx, double scale) {
  *ny = ceil((b.yh - b.yl) / dx) + 1;
  *nx = ceil((b.xh - b.xl) / dx) + 1;
  char **mask = createMask(*ny, *nx);

  int i, j;
  for (i = 0 ; i < *ny ; i++)
    for (j = 0 ; j < *nx ; j++)
      mask[i][j] = inside_billiard(j * dx / scale, i * dx / scale, &b);

  return mask;
}

/*
  free memory used by a grid

  input:
         grid - array to be freed
         ny - y dimension of grid
*/
void destroyGrid(double **grid, int ny) {
  int i;
  for (i = 0 ; i < ny ; i++)
    free(grid[i]);
  free(grid);
}

/*
  free memory used by a mask

  input:
         mask - array to be freed
         ny - y dimension of mask
*/
void destroyMask(char **mask, int ny) {
  int i;
  for (i = 0 ; i < ny ; i++)
    free(mask[i]);
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
    fprintf(stderr, "array2file: failed to open %s\n", file);
    return 1;
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
    fprintf(stderr, "intArray2file: failed to open %s\n", file);
    return 1;
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
    fprintf(stderr, "charArray2file: failed to open %s\n", file);
    return 1;
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
  for (i = 0 ; i < ny ; i++)
    for (j = 0 ; j < nx ; j++)
      if (!mask[i][j]) {
	grid[i][j] = INFINITY;
	counted[i][j] = -1;
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
