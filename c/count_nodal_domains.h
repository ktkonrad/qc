#ifndef _COUNT_NODAL_DOMAINS_H_
#define _COUNT_NODAL_DOMAINS_H_

#include <gsl/gsl_matrix.h>

int countNodalDomainsInterp(double **grid, char **mask, int ny, int nx, double k, double dx, int M, int upsample);
int findNextUnseen(int **counted, int *i, int *j, int ny, int nx);
void findDomainInterp(double **grid, int **counted, int i, int j, int nd, int ny, int nx, int upsample, gsl_matrix *interp);
void interpolate(double **grid, int **counted, int i, int j, int ny, int nx, int upsample, gsl_matrix *interp);

#endif
