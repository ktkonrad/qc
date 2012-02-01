#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_linalg.h>

#include <gsl/gsl_blas.h>

#include <stdio.h>
#include <math.h>

#include "interp_matrix.h"


#define R(x,y) sqrt((x)*(x)+(y)*(y))
#define THETA(x,y) atan2(y,x)
#define EPS (1e-40)

/*
  input: k:       wavenumber
         points:  points to evaluate bessles at
         npoints: number of points in points
	 M:       highest order bessel function to use

  output: return value: matrix A where
          A[i][j] = J_{j'}(k*r_i) * f(j'*\theta_i) where
	    0 <= i < npoints
	    0 <= j <= M
	    j' = j     if j <= M
	         j - M if j > M
	    f = 1   if j = 0
	        sin if 1 <= j <= M
	        cos if j > M
 */
gsl_matrix *bessel_matrix(double k, point *points, int npoints, int M) {
  gsl_matrix *out = gsl_matrix_alloc(npoints, 2*M+1);
  int i,j;
  double x, y, r, theta;
  for (i = 0 ; i < npoints ; i++) {
    x = points[i].x;
    y = points[i].y;
    r = R(x, y);
    theta = THETA(x,y);
    
    gsl_matrix_set(out, i, 0, gsl_sf_bessel_J0(k * r));
    for (j = 1 ; j <= M ; j++) {
      gsl_matrix_set(out, i, j, gsl_sf_bessel_Jn(j, k * r) * sin(j * theta));
    }
    for (j = 1 ; j <= M ; j++) {
      gsl_matrix_set(out, i, j+M, gsl_sf_bessel_Jn(j, k * r) * cos(j * theta));
    }
  }

  return out;
}
/*
  Compute the pseudoinverse of A and store it in A

  A+ = VS+U*
  A m x n
  V m x m
  S n x n
  U m x n
 */
int pseudoinverse(gsl_matrix *A, gsl_matrix *A_plus) {
  int rc;
  float temp;
  int m = (int)(A->size1);
  int n = (int)(A->size2);
  gsl_vector *workspace = gsl_vector_alloc(n);
  gsl_matrix *V = gsl_matrix_alloc(n, n);
  gsl_matrix *S = gsl_matrix_alloc(n, n);
  gsl_matrix *U_times_S_plus_trans = gsl_matrix_alloc(m, n);
  gsl_vector *singular_values = gsl_vector_alloc(n);
  int i;

  rc = gsl_linalg_SV_decomp(A, V, singular_values, workspace); // NOTE: this does a thin SVD
  // Now A = U

  if (rc) {
    // ERROR
    return rc;
  }
  
  // compute S+
  // by taking reciprocal of nonzero entries
  for (i = 0 ; i < n ; i++) {
    temp = gsl_vector_get(singular_values, i);
    if (temp < EPS) {
      break; // singular values are non-negative and form a non-increasing sequence
    }
    gsl_matrix_set(S, i, i, 1 / temp);
  }
  // now S = S+*

  // compute US+*
  rc = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, S, 0.0, U_times_S_plus_trans);
  // now A = US+*

  if (rc) {
    // ERROR
    return rc;
  }

  // compute VS+U
  rc = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, U_times_S_plus_trans, 0.0, A_plus);
  // now A_plus = VS+U
  if (rc) {
    // ERROR
    return rc;
  }


  // cleanup
  gsl_vector_free(workspace);
  gsl_matrix_free(V);
  gsl_matrix_free(S);
  gsl_matrix_free(U_times_S_plus_trans);
  gsl_vector_free(singular_values);

  return 0;
}

/*
    
 */
gsl_matrix *interp_matrix(double k, point *points_in, int npoints_in, point *points_out, int npoints_out, int M) {
  gsl_matrix *A = bessel_matrix(k, points_in, npoints_in, M);
  dump_matrix(A, "A.dat");
  gsl_matrix *B = bessel_matrix(k, points_out, npoints_out, M);
  dump_matrix(B, "B.dat");
  gsl_matrix *A_plus = gsl_matrix_alloc(2*M+1, npoints_in);
  gsl_matrix *interp = gsl_matrix_alloc(npoints_out, npoints_in);
  pseudoinverse(A, A_plus);
  dump_matrix(A_plus, "A_plus.dat");
  
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, B, A_plus, 0.0, interp);
  dump_matrix(interp, "interp.dat");

  gsl_matrix_free(A);
  gsl_matrix_free(A_plus);
  gsl_matrix_free(B);
  return interp;
}

int main(int argc, char **argv) {
#define npoints_in 24
  point points_in[npoints_in] = {{-0.500000,0.500000},{0.500000,0.500000},{-0.500000,-0.500000},{0.500000,-0.500000},{-0.500000,1.500000},{0.500000,1.500000},{-1.500000,0.500000},{-1.500000,-0.500000},{-0.500000,-1.500000},{0.500000,-1.500000},{1.500000,0.500000},{1.500000,-0.500000},{-1.500000,1.500000},{1.500000,1.500000},{-1.500000,-1.500000},{1.500000,-1.500000},{-0.500000,2.500000},{0.500000,2.500000},{-2.500000,0.500000},{-2.500000,-0.500000},{-0.500000,-2.500000},{0.500000,-2.500000},{2.500000,0.500000},{2.500000,-0.500000}};

#define npoints_out 121
  point points_out[npoints_out] = {{-0.500000,0.500000},{-0.500000,0.400000},{-0.500000,0.300000},{-0.500000,0.200000},{-0.500000,0.100000},{-0.500000,0.000000},{-0.500000,-0.100000},{-0.500000,-0.200000},{-0.500000,-0.300000},{-0.500000,-0.400000},{-0.500000,-0.500000},{-0.400000,0.500000},{-0.400000,0.400000},{-0.400000,0.300000},{-0.400000,0.200000},{-0.400000,0.100000},{-0.400000,0.000000},{-0.400000,-0.100000},{-0.400000,-0.200000},{-0.400000,-0.300000},{-0.400000,-0.400000},{-0.400000,-0.500000},{-0.300000,0.500000},{-0.300000,0.400000},{-0.300000,0.300000},{-0.300000,0.200000},{-0.300000,0.100000},{-0.300000,0.000000},{-0.300000,-0.100000},{-0.300000,-0.200000},{-0.300000,-0.300000},{-0.300000,-0.400000},{-0.300000,-0.500000},{-0.200000,0.500000},{-0.200000,0.400000},{-0.200000,0.300000},{-0.200000,0.200000},{-0.200000,0.100000},{-0.200000,0.000000},{-0.200000,-0.100000},{-0.200000,-0.200000},{-0.200000,-0.300000},{-0.200000,-0.400000},{-0.200000,-0.500000},{-0.100000,0.500000},{-0.100000,0.400000},{-0.100000,0.300000},{-0.100000,0.200000},{-0.100000,0.100000},{-0.100000,0.000000},{-0.100000,-0.100000},{-0.100000,-0.200000},{-0.100000,-0.300000},{-0.100000,-0.400000},{-0.100000,-0.500000},{0.000000,0.500000},{0.000000,0.400000},{0.000000,0.300000},{0.000000,0.200000},{0.000000,0.100000},{0.000000,0.000000},{0.000000,-0.100000},{0.000000,-0.200000},{0.000000,-0.300000},{0.000000,-0.400000},{0.000000,-0.500000},{0.100000,0.500000},{0.100000,0.400000},{0.100000,0.300000},{0.100000,0.200000},{0.100000,0.100000},{0.100000,0.000000},{0.100000,-0.100000},{0.100000,-0.200000},{0.100000,-0.300000},{0.100000,-0.400000},{0.100000,-0.500000},{0.200000,0.500000},{0.200000,0.400000},{0.200000,0.300000},{0.200000,0.200000},{0.200000,0.100000},{0.200000,0.000000},{0.200000,-0.100000},{0.200000,-0.200000},{0.200000,-0.300000},{0.200000,-0.400000},{0.200000,-0.500000},{0.300000,0.500000},{0.300000,0.400000},{0.300000,0.300000},{0.300000,0.200000},{0.300000,0.100000},{0.300000,0.000000},{0.300000,-0.100000},{0.300000,-0.200000},{0.300000,-0.300000},{0.300000,-0.400000},{0.300000,-0.500000},{0.400000,0.500000},{0.400000,0.400000},{0.400000,0.300000},{0.400000,0.200000},{0.400000,0.100000},{0.400000,0.000000},{0.400000,-0.100000},{0.400000,-0.200000},{0.400000,-0.300000},{0.400000,-0.400000},{0.400000,-0.500000},{0.500000,0.500000},{0.500000,0.400000},{0.500000,0.300000},{0.500000,0.200000},{0.500000,0.100000},{0.500000,0.000000},{0.500000,-0.100000},{0.500000,-0.200000},{0.500000,-0.300000},{0.500000,-0.400000},{0.500000,-0.500000}};

  int i;
  double dx = 0.001;
  int M = 8;

  for (i = 0 ; i < npoints_in ; i++) {
    points_in[i].x *= dx;
    points_in[i].y *= dx;
  }

  for (i = 0 ; i < npoints_out ; i++) {
    points_out[i].x *= dx;
    points_out[i].y *= dx;
  }

  gsl_matrix *interp = interp_matrix(200, points_in, npoints_in, points_out, npoints_out, M);

  printf("%zd x %zd\n", interp->size1, interp->size2);
  gsl_matrix_fprintf(stdout, interp, "%g");

  gsl_matrix_free(interp);

}

// for debugging
void dump_matrix(gsl_matrix *m, char *filename) {
  FILE *f = fopen(filename, "w");
  gsl_matrix_fprintf(f, m, "%g");
  fclose(f);
}
