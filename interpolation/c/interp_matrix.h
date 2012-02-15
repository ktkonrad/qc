#ifndef _INTERP_MATRIX_H_
#define _INTERP_MATRIX_H_

typedef struct {
  double x;
  double y;
} point;

gsl_matrix *bessel_matrix(double k, point *points, int npoints, int M);
gsl_matrix *interp_matrix(double k, point *points_in, int npoints_in, point *points_out, int npoints_out, int M);
int pseudoinverse(gsl_matrix *A, gsl_matrix *A_plus);
void dump_matrix(gsl_matrix *m, char *filename);
void dump_vector(gsl_vector *m, char *filename);


#endif
