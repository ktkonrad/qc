#ifndef _COUNT_NODAL_DOMAINS_H_
#define _COUNT_NODAL_DOMAINS_H_

#include <gsl/gsl_matrix.h>

#define SMALL_DOMAIN_SIZE 20

typedef struct {
  int small_domain_count;
  int interp_count;
  int boundary_interp_count;
} interp_stats;

int countNodalDomainsInterp(double **grid, char **mask, int ny, int nx, double k, double dx, int M, int upsample, interp_stats *status);
int findNextUnseen(int **counted, int *i, int *j, int ny, int nx);
void findDomainInterp(double **grid, int **counted, int i, int j, int nd, int ny, int nx, int upsample, gsl_matrix *interp, interp_stats *stats);
void interpolate(double **grid, int **counted, int i, int j, int upsample, gsl_matrix *interp, interp_stats *stats);


#endif
