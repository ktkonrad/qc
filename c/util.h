#ifndef _UTIL_H_
#define _UTIL_H_

#include "../vergini/billiard.h"

double **readOneSta(char *file, int *m, int *n);

char **readMask(char *file, int *m, int *n);

double **createGrid(int ny, int nx);

char **createMask(int ny, int nx);

char **createMaskFromBilliard(Billiard b, double dx, int *masky, int *maskx);

void destroyGrid(double **grid, int nx);

void destroyMask(char **mask, int nx);

void array2file(double **array, int m, int n, char *file);

void intArray2file(int **array, int m, int n, char *file);

void charArray2file(char **array, int m, int n, char *file);

void applyMask(double **grid, int **counted, char **mask, int ny, int nx);

#endif
