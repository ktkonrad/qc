#ifndef _INTERP_MATRIX_H_
#define _INTERP_MATRIX_H_

#include <gsl/gsl_matrix.h>

typedef struct {
  double x;
  double y;
} point;

point *stencil(void);
gsl_matrix *bessel_matrix(double alpha, point *points, int npoints, int M, double r_typical);
gsl_matrix *interp_matrix(double alpha, point *points_in, int npoints_in, point *points_out, int npoints_out, int M, double r_typical);
gsl_matrix *create_interp_matrix(double alpha, int M, int upsample);
int pseudoinverse(gsl_matrix *A, gsl_matrix *A_plus);
void dump_matrix(gsl_matrix *m, char *filename);
void dump_vector(gsl_vector *m, char *filename);


#endif
