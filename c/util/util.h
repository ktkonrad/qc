#ifndef _UTIL_H_
#define _UTIL_H_

#include <gsl/gsl_matrix.h> // for gsl_matrix
#include "exit_codes.h"

#define ERROR(fmt, args...) if(verb) { fprintf(stderr, "Error: %s: %s: %d: "fmt"\n", __FILE__, __FUNCTION__,  __LINE__, ## args); }
#define MALLOC_CHECK(ptr) if(!ptr) { ERROR("malloc failed!"); exit(OUT_OF_MEMORY_ERR); }
#define SET(loc, val) do {loc = (char *)malloc((strlen(val)+1)*sizeof(char)); MALLOC_CHECK(loc); strcpy(loc, val);} while (0)

double **readOneSta(char *file, int *m, int *n);
double **readSta(char *file, int *ne, int *m, int *n, double *k, int l);
char **readMask(char *file, int *m, int *n);
int array2file(double **array, int m, int n, char *file);
int intArray2file(int **array, int m, int n, char *file);
int charArray2file(char **array, int m, int n, char *file);
void applyMask(double **grid, int **counted, char **mask, int ny, int nx);
double wingTipMass(double **grid, char **mask, int ny, int nx);
int fillInterpMatrix(double k, double dx, int M, int upsample, gsl_matrix *m);
int **imatrix(int rows, int cols);
void free_imatrix(int **m);
char **cmatrix(int rows, int cols);
void free_cmatrix(char **m);
double **dmatrix(int ny, int nx);
void free_dmatrix(double **grid);

#endif
