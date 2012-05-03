/*
  get accurate timing comparisons for basis evaluation vs coefficient multiplications
 */

#include "../vergini/basis.h"
#include "../vergini/billiard.h"
#include "../vergini/colloc.h"
#include "../vergini/nrutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#define STADIUM 0
#define SINAI 1

int verb = 1;

int main(int argc, char **argv) {
  Basis_Set bas;
  Billiard bil;
  Bdry_Pt_Set bps;
  int i,j,n,N;
  double **coeffs, *x, *y, **psi, *values;
  double k = 500.0, temp;
  time_t start;
  double basis_time, coeff_time;
  char arg;
  double b = 10;
  double k_base = 20;

  if (argc < 7) {
    printf("usage: ./eval -l billiard_args -s basis_args -n number_of_evals\n");
    exit(-1);
  }

  // parse command line args
  while ((arg = getopt(argc, argv, "l:s:n:")) != -1) {
    switch (arg) {
    case 'l':
      if (parse_billiard(optarg, &bil)==-1) {
        printf("failed to parse billiard: %s\n", optarg);    
        exit(-1);
      }
      break;
    case 's':
      if (parse_basis(optarg, &bas)==-1) {
        printf("failed to parse basis: %s\n", optarg);    
        exit(-1);
      }
    case 'n':
      n = atoi(optarg);
      break;
    }
  }
  build_billiard(&bil, k_base);
  build_bdry(b, k_base, &bil, &bps);
  build_basis(k_base, &bil, &bas);

  show_billiard_properties(&bil, k_base);                                                                                      
  show_colloc_properties(&bil, &bps, k_base);                                                                                  
  show_basis_properties(&bil, &bas);

  N = bas.N;

  // initialize
  coeffs = dmatrix(0,n,0,N);
  values = dvector(0,n);
  for (i = 0; i < N ; i++) {
    for (j = 0 ; j < N ; j++) {
      coeffs[i][j] = (double)rand() / RAND_MAX;
    }
  }

  x = dvector(0,n);
  y = dvector(0,n);
  psi = dmatrix(0,n,0,N);
  for (i = 0 ; i < n ; i++) {
      x[i] = (double)rand() / RAND_MAX;
      y[i] = (double)rand() / RAND_MAX;
  }

  // run & time
  start = clock();
  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < bas.N ; j++) {
      psi[i][j] = eval_basis(k*x[i], k*y[i], j, &bas);
    }
  }
  basis_time = (double)(clock() - start) / CLOCKS_PER_SEC;

  start = clock();  
  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < bas.N ; j++) {
      values[i] += psi[i][j] * coeffs[i][j];
    }
  }
  coeff_time = (double)(clock() - start) / CLOCKS_PER_SEC;
  
  // print results
  printf("basis evals took %g seconds\n", basis_time);
  printf("coeff mults took %g seconds\n", coeff_time);
  printf("ratio: %f\n", basis_time / coeff_time);

  free_dvector(x, 0, n);
  free_dvector(y, 0, n);
  free_dmatrix(psi, 0, n, 0, N);
  free_dmatrix(coeffs, 0, n, 0, N);
  free_dvector(values, 0, n);

}
