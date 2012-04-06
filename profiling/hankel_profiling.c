/*
  code to profile hankel evaluations

  Kyle Konrad
  4/4/2012
*/

#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define GSL_HANKEL(nu, r) (cos((nu)*M_PI)*gsl_sf_bessel_Jnu((nu), (r)) - sin((nu)*M_PI)*gsl_sf_bessel_Ynu((nu), (r)))
#define MAX_NU 12

// globals
int false = 0;

/*
// fortran subroutine
extern void hank103_(double *z, double *h0, double *h1, int *ifexpon);
*/

void gsl_hankel_eval(int n, double *nus, double *rs) {
  // evaluate n hankel functions with given arguments using gsl
  int i;
  double temp;
  int nu;
  double r;

  for (i = 0 ; i < n ; i++) {
    nu = nus[i];
    r = rs[i];
    temp = GSL_HANKEL(nu, r);
  }
}

/*
void fortran_hankel_eval(int n, double *rs, double *h0s, double *h1s) {
  // evaluate n hankel functions with given arguments
  int i;
  int nu;
  double r;

  for (i = 0 ; i < n ; i++) {
    r = rs[i];
    hank103_(&r, &h0s[i], &h1s[i], &false);
  }
}
*/

void math_hankel_eval(int n, double *rs) {
  // evaluate n hankel functions with given arguments using bessel functions in math.h
  int i;
  double temp;
  double r;

  for (i = 0 ; i < n ; i++) {
    r = rs[i];
    temp = j0(r);
    temp = y0(r);
  }
}

void generate_rs(int n, double *rs) {
  int i;
  for (i = 0 ; i < n ; i++) {
    rs[i] = (double)rand()/(double)RAND_MAX; // uniform over [0.0, 1.0]
  }
}

void generate_nus(int n, double *nus) {
  int i;
  for (i = 0 ; i < n ; i++) {
    nus[i] = rand() % MAX_NU; // uniform over [0, MAX_NU)
  }
}

void main(int argc, char **argv) {
  int n;
  double *nus, *rs;
  double *h0s, *h1s;
  clock_t start_time, end_time;
  double elapsed_seconds;

  // argument parsing
  if (argc < 2) {
    perror("usage: ./hankel n");
    exit(-1);
  }
  n = atoi(argv[1]);

  // memory allocation
  nus = malloc(n * sizeof(double));
  rs = malloc(n * sizeof(double));
  h0s = malloc(n * sizeof(double));
  h1s = malloc(n * sizeof(double));
  
  // parameter generation
  generate_rs(n, rs);
  generate_nus(n, nus);

  // gsl timining
  start_time = clock();
  gsl_hankel_eval(n, nus, rs);
  end_time = clock();
  elapsed_seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  printf("gsl completed %d evalutions in %.4f seconds\n", n, elapsed_seconds);
  printf("%g evaluations per second\n", n / elapsed_seconds);
  printf("%.16g seconds per evaluation\n", elapsed_seconds / n);

  // math timining
  start_time = clock();
  math_hankel_eval(n, rs);
  end_time = clock();
  elapsed_seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  printf("math completed %d evalutions in %.4f seconds\n", n, elapsed_seconds);
  printf("%g evaluations per second\n", n / elapsed_seconds);
  printf("%.16g seconds per evaluation\n", elapsed_seconds / n);

  /*
  // fortran timing
  start_time = clock();
  fortran_hankel_eval(n, rs, h0s, h1s);
  end_time = clock();
  elapsed_seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  printf("\nfortran completed %d evalutions in %.4f seconds\n", n, elapsed_seconds);
  printf("%g evaluations per second\n", n / elapsed_seconds);
  printf("%.16g seconds per evaluation\n", elapsed_seconds / n);
  */

  /*
  // memory freeing
  free(rs);
  free(nus);
  free(h0s);
  free(h1s);
  */
}
