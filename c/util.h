#ifndef _UTIL_H_
#define _UTIL_H_

double **readSta(char *file, int *m, int *n);

char **readMask(char *file, int *m, int *n);

double **createGrid(int ny, int nx);

char **createMask(int ny, int nx);

void destroyGrid(double **grid, int nx);

void destroyMask(char **mask, int nx);

void array2file(double **array, int m, int n, char *file);

void intArray2file(int **array, int m, int n, char *file);

void applyMask(double **grid, int **counted, char **mask, int ny, int nx);

#endif
