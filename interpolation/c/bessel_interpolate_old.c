/*
  testing out some functions to do fourier-bessel expansion
  Kyle Konrad
  8/6/2011
*/


#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#ifndef DEBUG
#define DEBUG 0
#endif
#define PRINT_DEBUG(v, fmt, args...)  if(DEBUG >= v){printf("DEBUG: "fmt, ##args);}

// GLOBALS
double b = 1; // work with functions on [0, b]
double dx = .001; // work with functions sampled with this resolution
int alpha = 2; // use J_alpha basis

/* compute the coefficient of the nth term in Fourier-Bessel series of f
   f_x_weighted is f(z)*z for all z in the input vector x
   n should be positive
 */
double fourier_bessel_coeff(gsl_vector *x, gsl_vector *f_x_weighted, int n) {
  gsl_vector *J_alpha = gsl_vector_alloc(x->size); // alloc J_alpha


  // lambda_n is the nth zero of J_alpha
  double lambda_n = gsl_sf_bessel_zero_Jnu(alpha, n);
  PRINT_DEBUG(2, "lambda: %f\n", lambda_n);

  // J_alpha(lambda_n * x / b) TODO: combine this loop with x loop above
  int i;
  for (i = 0 ; i < J_alpha->size ; i++) {
    gsl_vector_set(J_alpha, i, gsl_sf_bessel_Jn(alpha, lambda_n / b * gsl_vector_get(x, i)));
  }
  
  // <J_alpha, f> (weighted inner product)
  double numerator;
  gsl_blas_ddot(J_alpha, f_x_weighted, &numerator);
  numerator *= dx;

  // ||J_alpha||^2 (norm squared of J_alpha)
  double denominator = b*b/2 * pow(gsl_sf_bessel_Jn(alpha+1, lambda_n), 2);

  double coeff = numerator/denominator;

  PRINT_DEBUG(2, "numerator: %f\n", numerator);
  PRINT_DEBUG(2, "denominator: %f\n", denominator);
  PRINT_DEBUG(1, "c_%d: %f\n", n, coeff);
  double coeff_should = 2 / lambda_n / gsl_sf_bessel_Jn(alpha+1, lambda_n);
  PRINT_DEBUG(1, "s_%d: %f\n", n, coeff_should);

  gsl_vector_free(J_alpha); // free_J_alpha

  return coeff;
}

/* evualate a Fourier-Bessel series
   coeffs are coefficients of terms in series
   x is points to evualate at
   results are stored in interpolated
*/
void evaluate_fourier_bessel_series(gsl_vector *x, gsl_vector *coeffs, gsl_vector *interpolated) {
  if (interpolated->size != x->size) {
    printf("Invalid vector sizes in interpolate_with_coeffs");
    return;
  }

  gsl_vector *temp = gsl_vector_alloc(interpolated->size); // alloc temp

  gsl_vector_set_zero(interpolated);
  
  int n, j;
  for (n = 0 ; n < coeffs->size ; n++) {
    // lambda_n is the nth zero of J_alpha
    if (alpha == 0 && n == 0) {continue;}
    double lambda_n = gsl_sf_bessel_zero_Jnu(alpha, n);
    for (j = 0 ; j < interpolated->size ; j++) {
      gsl_vector_set(temp, j, -gsl_vector_get(coeffs, n) * gsl_sf_bessel_Jn(alpha, lambda_n / b * gsl_vector_get(x, j))); // negative sign should not be necessary!
    }
    gsl_vector_add(interpolated, temp);
  }
	   
  gsl_vector_free(temp); // free temp
}

/* perform function interpolation using Fourier-Bessel series
   use N terms of series
   upsample by factor of s
   store interpolated function values in interpolated
*/
void bessel_interpolate(double (*f)(double), int N, int s, gsl_vector *interpolated) {
  gsl_vector *x = gsl_vector_alloc((int)(b / dx)); // alloc x
  gsl_vector *f_x = gsl_vector_alloc(x->size); // alloc f_x
  gsl_vector *coeffs = gsl_vector_alloc(N); // alloc coeffs
  gsl_vector *x_interpolated = gsl_vector_alloc((x->size)*s+1); // alloc x_interpolated

  // x
  int i;
  for (i = 0 ; i < x->size ; i++) {
    gsl_vector_set(x, i, i * dx + dx / 2); //center points on intervals for more accuracy (TODO: maybe use simposon's rule or something more sophisticated?)
  }

  // f(x)
  for (i = 0 ; i < x->size ; i++) {
    gsl_vector_set(f_x, i, (*f)(gsl_vector_get(x, i)));
  }

  // f(x)*x
  gsl_vector_mul(f_x, x); // result is stored in f_x

  // TODO: make this more effiecient by combining these two loops
  // get coefficients
  int n;
  for (n = 0 ; n < N ; n++) {
    gsl_vector_set(coeffs, n, fourier_bessel_coeff(x, f_x, n+1));
  }
  
  // evaluate
  for (i = 0 ; i < x_interpolated->size ; i++) {
    gsl_vector_set(x_interpolated, i, i*dx/s);
  }    
  evaluate_fourier_bessel_series(x_interpolated, coeffs, interpolated);

  gsl_vector_free(x);  
  gsl_vector_free(f_x);
  gsl_vector_free(coeffs);
  gsl_vector_free(x_interpolated);
}

double square(double x) {
  return x * x;
}

int main(int argc, char **argv) {
  int m = (int)(b / dx + 1);
  int N = 20; // compute N terms of Fourier-Bessel series  
  int s = 10; // interpolate s points between every input point

  gsl_vector *interpolated = gsl_vector_alloc((m-1)*s+1);

  bessel_interpolate(&square, N, s, interpolated);

  int i;
  for (i = 0 ; i < interpolated->size ; i ++) {
    //    printf("%f\n", gsl_vector_get(interpolated, i));
  }
  
}
