#ifndef _COUNT_UTIL_H_
#define _COUNT_UTIL_H_

#include <gsl/gsl_matrix.h> // for gsl_matrix
#include "exit_codes.h"
#include <limits.h>

// special values
#define UNCOUNTED INT_MAX
#define MASKED INT_MIN
// boolean functions to check for counted or masked
#define IS_COUNTED(c) (c < BL_DISCONNECTED)
#define IS_MASKED(c) (c == MASKED)

double **readOneSta(char *file, int *m, int *n);
double **readSta(char *file, int *ne, int *m, int *n, double *k, int l);
int array2file(double **array, int m, int n, char *file);
int intArray2file(int **array, int m, int n, char *file);
int charArray2file(char **array, int m, int n, char *file);
double wingTipMass(double **grid, int ny, int nx);
int fillInterpMatrix(double k, double dx, int M, int upsample, gsl_matrix *m);

#endif
