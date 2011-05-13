/* utility to read single eigenfunction from sta file and output to a file
   
   Kyle Konrad
   5/9/2011
*/

#include <stdio.h>
#include <stdlib.h>

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
  write an array of doubles to a file, tab separator

  input:
         array - array to write
         m     - rows in array
	 n     - columns in array
	 file  - name of file to write to
*/
void array2file(double **array, int m, int n, char *file) {
  int i, j;
  FILE *out = fopen(file, "w");
  for (i = 0 ; i < m ; i++) {
    for (j = 0 ; j < n ; j++) {
      if (j != 0)
	fprintf(out, "\t");
      fprintf(out, "%f", array[i][j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: ./convert [infile] [outfile] [i]\n");
    return 1;
  }

  char *infile = argv[1];
  char *outfile = argv[2];
  int i = atoi(argv[3]);

  int ne, ny, nx;
  double k;

  double **func = readSta(infile, &ne, &ny, &nx, &k, i);
  array2file(func, ny, nx, outfile);
 
  return 0;
}
