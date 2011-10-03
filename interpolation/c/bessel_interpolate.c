#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_linalg.h>

#include <gsl/gsl_blas.h>

#include <stdio.h>

/*
evaluate the solution to the helmholtz equation given by
c1*exp(ikr) + c2*exp(-ikr)
at r
*/
gsl_complex helmholtz_soln_eval(double c1, double c2, double k, double r) {
  return gsl_complex_add(
			 gsl_complex_mul(
					 gsl_complex_rect(c1,0),
					 gsl_complex_exp(
							 gsl_complex_rect(0,k*r)
							 )
					 ),
			 gsl_complex_mul(
					 gsl_complex_rect(c2, 0),
					 gsl_complex_exp(
							 gsl_complex_rect(0,-k*r)
							 )
					 )
			 );
}

gsl_vector_complex *helmholtz_soln_vector(double c1, double c2, double k, gsl_vector *r) {
  gsl_vector_complex *vec = gsl_vector_complex_alloc(r->size);
  int i;
  for (i = 0 ; i < r->size ; i++) {
    gsl_vector_complex_set(vec, i, helmholtz_soln_eval(c1, c2, k, gsl_vector_get(r, i)));
  }
  return vec;
}

/*
return a matrix A where
A[i,j] = J_j(r_i) * exp(imag*j*theta_i)
0 <= i < m (where m == r->size)
0 <= j <= N
*/
gsl_matrix_complex *bessels_matrix(double k, gsl_vector *r, gsl_vector *theta, int N) {
  gsl_matrix_complex *A = gsl_matrix_complex_alloc(r->size, N+1);
  int i, j;
  double *bessel_values = malloc((N+1) * sizeof(double));
  for (i = 0 ; i < r->size ; i++) {
    gsl_sf_bessel_Jn_array(0, N, k * gsl_vector_get(r, i), bessel_values);
    for (j = 0 ; j <= N ; j++) {
      gsl_matrix_complex_set(A, i, j,
			     gsl_complex_mul(
					     gsl_complex_rect(bessel_values[j], 0),
					     gsl_complex_exp(
							     gsl_complex_rect(0, j*gsl_vector_get(theta, i)
									      )
							     )
					     )
			     );
    }
  }
  free(bessel_values);
  return A;
}

/*
compute coefficients using LU decomposition
preconditions:
  bessel_values is square
  f_values has same number of rows as bessel_values
*/
gsl_vector_complex *compute_coeffs(gsl_vector_complex *f_values, gsl_matrix_complex *bessel_values) {
  gsl_vector_complex *coeffs = gsl_vector_complex_alloc(bessel_values->size2);
  int s;
  gsl_permutation *p = gsl_permutation_alloc(bessel_values->size2);
  gsl_linalg_complex_LU_decomp(bessel_values, p, &s);
  gsl_linalg_complex_LU_solve(bessel_values, p, f_values, coeffs);
  return coeffs;
}

/*
compute coefficients using least squares method
*/
gsl_vector_complex *compute_coeffs_least_squares(gsl_vector_complex *f_values, gsl_matrix_complex *bessel_values) {
  gsl_multifit_linear_workspace *workspace = gsl_multifit_linear_alloc(bessel_values->size1, bessel_values->size2)
  gsl_multifit_linear_complex();

}

/*
recompute original function from given bessel coefficients
*/
gsl_vector_complex *recompute_from_coeffs(gsl_vector_complex *coeffs, gsl_vector *r, gsl_vector *theta) {
  gsl_vector_complex *recomputed = gsl_vector_complex_alloc(r->size);
  gsl_vector_complex *phases = gsl_vector_complex_alloc(coeffs->size);
  int i, j;
  double bessel_vals[coeffs->size];
  double current_theta;
  gsl_complex temp;
  for (i = 0 ; i < r->size ; i++) { // loop over points
    temp = gsl_complex_rect(0, 0);
    gsl_sf_bessel_Jn_array(0, coeffs->size, gsl_vector_get(r, i), bessel_vals);
    current_theta = gsl_vector_get(theta, i);
    for (j = 0 ; j < coeffs->size ; j++) { // loop over bessel fns
      temp = gsl_complex_add(temp,
			     gsl_complex_mul(gsl_vector_complex_get(coeffs, j),
					     gsl_complex_mul(
							     gsl_complex_rect(bessel_vals[j], 0),
							     gsl_complex_exp(
									     gsl_complex_rect(0, j * current_theta)
									     )
							     )
					     )
			     );
    }
    gsl_vector_complex_set(recomputed, i, temp);
  }
  return recomputed;
}

/*
print complex vector to given file
*/
my_vector_complex_fprint(FILE *outfile, gsl_vector_complex *vec) {
  int i;
  gsl_complex val;
  for (i = 0 ; i < vec->size ; i++) {
    val = gsl_vector_complex_get(vec, i);
    fprintf(outfile, "%f + %fi\n", GSL_REAL(val), GSL_IMAG(val));
  }
}

/*
generate a grid on [-n*dx, n*dx]^2
with spacing dx
put polar coordinates of grid points into r and theta
*/
void polar_grid(int n, double dx, gsl_vector *r, gsl_vector *theta) {
  double x, y;
  gsl_sf_result r_result, theta_result;
  int i = 0;
  for (x = -n*dx ; x <= n*dx ; x += dx) {
    for (y = -n*dx ; y <= n*dx ; y += dx) {
      if (x == 0 && y == 0) { // gsl_sf_rect_to_polar errors in this case
	gsl_vector_set(r, i, 0);
	gsl_vector_set(theta, i++, 0);	
      }
      else {
	gsl_sf_rect_to_polar(x, y, &r_result, &theta_result); // NOTE: check error values
	gsl_vector_set(r, i, r_result.val);
	gsl_vector_set(theta, i++, theta_result.val);
      }
    }
  }
}

int main() {
  double c1 = 0.5;
  double c2 = 0.3;
  double k = 1.0;
  int n = 2;
  double dx = 0.1;
  int m = (2*n+1)*(2*n+1);
  int N = m - 1;

  gsl_vector *r = gsl_vector_alloc(m);
  gsl_vector *theta = gsl_vector_alloc(m);
 
  polar_grid(n, dx, r, theta);
  
  gsl_matrix_complex *bessel_values = bessels_matrix(k, r, theta, N);
  gsl_vector_complex *f_values = helmholtz_soln_vector(c1, c2, k, r);

  gsl_vector_complex *coeffs = compute_coeffs(f_values, bessel_values);

  gsl_vector_complex *recomputed_f_values = recompute_from_coeffs(coeffs, r, theta);

  bessel_values = bessels_matrix(k, r, theta, N); // recompute this since bessel_values got LU-reduced when finding coeffs
  gsl_vector_complex *recomputed_f_values2 = gsl_vector_complex_alloc(r->size);
  gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1, 0), bessel_values, coeffs, gsl_complex_rect(0, 0), recomputed_f_values2);
  
  my_vector_complex_fprint(stdout, coeffs);
  
  FILE *outfile = fopen("out.dat", "w");
  my_vector_complex_fprint(outfile, f_values);
  fprintf(outfile, "*****\n");
  my_vector_complex_fprint(outfile, recomputed_f_values);
  fprintf(outfile, "*****\n");
  my_vector_complex_fprint(outfile, recomputed_f_values2);
  fprintf(outfile, "*****\n");
  gsl_matrix_complex_fprintf(outfile, bessel_values, "%f");

  fclose(outfile);

  gsl_matrix_complex_free(bessel_values);
  gsl_vector_complex_free(f_values);
  gsl_vector_complex_free(coeffs);
  gsl_vector_complex_free(recomputed_f_values);
  gsl_vector_complex_free(recomputed_f_values2);
  gsl_vector_free(r);
  gsl_vector_free(theta);

}
