#ifndef _COUNT_NODAL_DOMAINS_H_
#define _COUNT_NODAL_DOMAINS_H_

#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include "util/bit_array.h"

#define SMALL_DOMAIN_SIZE 5

typedef struct {
  int small_domain_count;
  int interp_count;
  int boundary_trouble_count;
  int edge_trouble_count;
} interp_stats;

typedef struct {
  gsl_vector *input;
  gsl_vector *output;
} interp_workspace;

int countNodalDomainsInterp(double **grid, int **counted, int ny, int nx, double alpha, int M, int upsample, interp_stats *stats, FILE *sizefile);
int countNodalDomains(bit_array_t *signs, int **counted, FILE *sizefile);
int countNodalDomainsNoInterp(double **grid, int **counted, int ny, int nx, FILE *sizefile);
int findNextUnseen(int **counted, int *i, int *j, int ny, int nx);
int findDomain(bit_array_t *signs, int **counted, int i, int j, int nd, int ny, int nx);
int findDomainNoInterp(double **grid, int **counted, int i, int j, int nd, int ny, int nx);
void interpolate(double **grid, bit_array_t *counted, int i, int j, int ny, int nx, int upsample, gsl_matrix *interp, interp_stats *stats, interp_workspace *w);
bit_array_t *upsample(double **grid, int **counted, int ny, int nx, double alpha, int M, int upsample, interp_stats *stats);
void fill_interp_input(double **grid, int i, int j, interp_workspace *w, int ny, int nx, int *refl_i, int *refl_j, int *refl_i_sign, int *refl_j_sign);
interp_workspace *new_interp_workspace(int upsample_ratio);
void free_interp_workspace(interp_workspace *w);

#endif
