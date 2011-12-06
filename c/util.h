#ifndef _UTIL_H_
#define _UTIL_H_

#include "../vergini/billiard.h" // for Billiard
#include <gsl/gsl_matrix.h> // for gsl_matrix

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
int fillInterpMatrix(int k, double dx, int M, int upsample, gsl_matrix *m);

#endif
