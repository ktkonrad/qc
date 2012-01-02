#ifndef _UTIL_H_
#define _UTIL_H_

#include "../vergini/billiard.h" // for Billiard
#include <gsl/gsl_matrix.h> // for gsl_matrix

#define ERROR(fmt, args...) fprintf(stderr, "Error: %s: %s: %d: "fmt"\n", __FILE__, __FUNCTION__,  __LINE__, ## args)

double **readOneSta(char *file, int *m, int *n);
double **readSta(char *file, int *ne, int *m, int *n, double *k, int l);
char **readMask(char *file, int *m, int *n);
double **createGrid(int ny, int nx);
char **createMask(int ny, int nx);
char **createMaskFromBilliard(Billiard b, double dx, int *masky, int *maskx);
char **createScaledMaskFromBilliard(Billiard b, double dx, int *masky, int *maskx, double scale);
void destroyGrid(double **grid, int nx);
void destroyMask(char **mask, int nx);
int array2file(double **array, int m, int n, char *file);
int intArray2file(int **array, int m, int n, char *file);
int charArray2file(char **array, int m, int n, char *file);
void applyMask(double **grid, int **counted, char **mask, int ny, int nx);
double wingTipMass(double **grid, char **mask, int ny, int nx);
int fillInterpMatrix(double k, double dx, int M, int upsample, gsl_matrix *m);
int **imatrix(int rows, int cols);
void free_imatrix(int **m);

#endif
