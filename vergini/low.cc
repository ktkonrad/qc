/* low.cc: Compute low eigenfunctions of the laplacian
 *
 * barnett 10/30/03
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* local low-level utilities... */
#include "cxml.h"
#include "nrutil.h"
#include "useful.h"

/* local higher-level utilities... */
#include "billiard.h"
#include "basis.h"
#include "colloc.h"
#include "grid.h"
#include "matrix.h"



// =============== Build DIRICHLET S, and exact T MATRICES ==================

void build_dirichlet_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			     double k, double *S, double *T)
  /* fills S and T, the area norm and perimeter norm matrices respectively.
     S is given by the Berry/Boasman bdry normalisation,
     correct for Dirichlet BC only.
     However, it is expected to be as good as the exact S for finding
     t-minima in a
     k-sweep, since the gen eigenvalue with smallest t/s is used.
     It is much (ie more than 2 times) faster to build than the exact S,
     only using first derivs.
     Alex Barnett 99/11/16
     10/30/03 modified for objects
     */

{
  int i,j,n, M = p->M, N = s->N;
  double val, ddx, ddy, xi, yi, norx, nory, srni;
  double *A, *B;

  // allocate rectangular matrices...
  A  = dvector(0,M*N-1);
  B  = dvector(0,M*N-1);

  
  // fill 2 matrices used to build S...
  for (i=1;i<=M;++i) {
    xi = k*p->x[i];
    yi = k*p->y[i];
    srni = sqrt(p->rn[i]);
    norx = cos(p->na[i]);
    nory = sin(p->na[i]);
    for(j=1;j<=N;++j) {
      eval_basis_deriv(ddx, ddy, xi, yi, j, s);
      n = i-1 + M*(j-1);
      A[n] = eval_basis(xi, yi, j, s);
      B[n] = srni*k*( norx*ddx + nory*ddy );
    }
  }
  
  // Build S = (1/2k^2).(L/M).B^T.B ...
  dgemm_('T','N',N,N,M, l->Perim/(2*M*k*k), B,M, B,M, 0.0,S,N);

  // T = (L/M) A^T.A ...
  dgemm_('T','N',N,N,M, l->Perim/M,         A,M, A,M, 0.0,T,N);

  free_dvector(A,0,M*N-1);
  free_dvector(B,0,M*N-1);
}




// =============== Build NEUMANN S, and exact T MATRICES ==================

void build_neumann_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			     double k, double *S, double *T)
/* fills S and T, the area norm and neumann perim norm matrices respectively.
   S is given by my Neumann-appropriate area norm formula of 10/29/03.
   Only used values and first derivs.
   
   barnett 10/30/03
 */

{
  int i,j,n, M = p->M, N = s->N;
  double val, ddx, ddy, xi, yi, norx, nory, srni, alphrsi;
  double *A, *B, *C, *D, *E;

  // allocate rectangular matrices...
  A = dvector(0,M*N-1);
  B = dvector(0,M*N-1);
  C = dvector(0,M*N-1);
  D = dvector(0,M*N-1);
  E = dvector(0,M*N-1);

  
  // fill 2 matrices used to build S...
  for (i=1;i<=M;++i) {
    xi = k*p->x[i];
    yi = k*p->y[i];
    srni = sqrt(p->rn[i]);
    alphrsi = p->alpha[i]*p->rs[i];
    norx = cos(p->na[i]);
    nory = sin(p->na[i]);
    for(j=1;j<=N;++j) {
      eval_basis_deriv(ddx, ddy, xi, yi, j, s);
      n = i-1 + M*(j-1);
      val = eval_basis(xi, yi, j, s);
      A[n] = srni * val;
      B[n] = srni * ( -nory*ddx + norx*ddy );
      C[n] = alphrsi * val;
      D[n] = -nory*ddx + norx*ddy;
      E[n] = norx*ddx + nory*ddy;
    }
  }
  
  // Build S = k^2.A^T.A - B^T.B - C^T.D ...
  dgemm_('T','N',N,N,M,  k*k, A,M, A,M, 0.0,S,N);
  dgemm_('T','N',N,N,M, -1.0, B,M, B,M, 1.0,S,N);
  dgemm_('T','N',N,N,M, -1.0, C,M, D,M, 1.0,S,N);
  // scale the sum to be the integral, divide all by 2k^2
  dscal_(N*N, l->Perim/(2*M*k*k), S, 1);

  // Build T = (L/M) E^T.E ...
  dgemm_('T','N',N,N,M, l->Perim/M, E,M, E,M, 0.0,T,N);

  free_dvector(A,0,M*N-1);
  free_dvector(B,0,M*N-1);
  free_dvector(C,0,M*N-1);
  free_dvector(D,0,M*N-1);
  free_dvector(E,0,M*N-1);
}




// ======= BUILD ALEX (arbitrary) S, and T (Dirichlet) MATRICES ===========

void build_alex_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			double k, double *S, double *T)
  /* fills S and T, the perimeter norm and area norm matrices respectively.
     Uses my S formula derived by taking k_2->k_1 limit of unequal k formula,
     involves up to 2nd derivs

     Alex Barnett 00/4/5
     */

{
  int i,j,u,v,n, N=s->N, M=p->M;
  double val, ddx, ddy, dxx, dyy, dxy, xi, yi, norx, nory, srni;
  double *A, *B, *C, *D;

  // allocate rectangular matrices...
  A  = dvector(0,M*N-1);
  C  = dvector(0,M*N-1);
  D  = dvector(0,M*N-1);
  B  = dvector(0,M*N-1);

  // fill 4 matrices used to build S...
  for (i=1;i<=M;++i) {
    xi = k*p->x[i]; // unitless coords.
    yi = k*p->y[i];
    norx = cos(p->na[i]);
    nory = sin(p->na[i]);
    for(j=1;j<=N;++j) {
      eval_basis_everything(val, ddx, ddy, dxx, dyy, dxy, xi, yi, j, s);
      n = i-1 + M*(j-1);
      A[n] = val;
      C[n] = xi*ddx + yi*ddy; // r.Del
      D[n] = k*( norx*ddx + nory*ddy ); // n.Del
      B[n] = k*( norx*(dxx*xi + dxy*yi) + nory*(dxy*xi + dyy*yi) ); // n.DelDel.r
    }
  }
  
  // Building S: S = -A^T.(B + D) ...
  dgemm_('T','N',N,N,M, -1.0, A,M, B,M, 0.0,S,N);
  dgemm_('T','N',N,N,M, -1.0, A,M, D,M, 1.0,S,N);
  // + C^T.D ...
  dgemm_('T','N',N,N,M, 1.0, C,M, D,M, 1.0,S,N);
  
  // scale by (1/2).(L/M)/(2 k^2)...   1/2 is to cancel the adding-transpose.
  dscal_(N*N, l->Perim/(4*M*k*k), S, 1);
 
  // add its transpose to itself...
  for (i=1;i<=N;++i)
    for(j=1;j<=i;++j) {
      S[u=i-1 + N*(j-1)] += S[v=j-1 + N*(i-1)];
      S[v] = S[u];
    }


  // T = (L/M) A^T.A ...
  dgemm_('T','N',N,N,M, l->Perim/M, A,M, A,M, 0.0,T,N);

  free_dvector(A,0,M*N-1);
  free_dvector(C,0,M*N-1);
  free_dvector(D,0,M*N-1);
  free_dvector(B,0,M*N-1);
}



// ------------------------------------------------------------------------

void show_usage()
  /* For low.cc only
   */
{
  fprintf(stderr, "Usage:\n\tlow outfilehead billiard_info b basis_info k_lo k_hi dk mode\n\n\tmode = 0:sweep 1:dir_state -1:neu_state\n");
  show_billiard_usage();
  show_basis_usage();
}



/* ---------------------------- CMD LINE ------------------------------ */
int read_cmd_line(int argc, char **argv, char *head, Billiard *l,\
		  double *b, Basis_Set *s, double *k_lo, double *k_hi,\
		  double *dk, int *mode)
  /* For test.cc only
   * Returns 0 if success
   */
{
  int i=1;          /* pointer to current argument */
  int r=0;

  printf("reading command line...\n\n");

  if (i<argc)
    sscanf(argv[i++], "%s", head);
  else {
    fprintf(stderr, "not enough cmd line args!\n");    return 1;
  }

  r = read_billiard_cmdline(argc-i, argv+i, l);
  if (r==-1)
    return 1;
  i += r;
  
  if (i<argc)
    sscanf(argv[i++], "%lf", b);
  else {
    fprintf(stderr, "not enough cmd line args!\n");
    return 1;
  }

  r = read_basis_cmdline(argc-i, argv+i, s);
  if (r==-1)
    return 1;
  i += r;

  if (i<argc)
    sscanf(argv[i++], "%lf", k_lo);
  else {
    fprintf(stderr, "not enough cmd line args!\n");
    return 1;
  }
  if (i<argc)
    sscanf(argv[i++], "%lf", k_hi);
  else {
    fprintf(stderr, "not enough cmd line args!\n");
    return 1;
  }
  if (i<argc)
    sscanf(argv[i++], "%lf", dk);
  else {
    fprintf(stderr, "not enough cmd line args!\n");
    return 1;
  }
  if (i<argc)
    sscanf(argv[i++], "%d", mode);
  else {
    fprintf(stderr, "not enough cmd line args!\n");
    return 1;
  }

  printf("\toutfilehead=%s, b=%g, k=[%g,%g], dk=%g mode=%d\n\n", \
	 head, *b, *k_lo, *k_hi, *dk, *mode);

  if (i<argc) {
    fprintf(stderr, "too many cmd line args!\n");
    return 1;
  }

  return 0;
}


/* ======================================= MAIN ======================== */
int main(int argc, char **argv)
{
  Billiard bil;
  Bdry_Pt_Set bps;
  Basis_Set bas;
  Grid g;
  FILE *fp;
  int i, j, ne, mode, N, rank;
  double b, k, k_lo, k_hi, dk, **a, *ks, dx, t, *S, *T, *X, *evals;
  char head[LEN], name[LEN];
  
  printf("\nlow: calc low efuncs, 10/30/03\n\n");
  if (read_cmd_line(argc, argv, head, &bil, &b, &bas, &k_lo, &k_hi, &dk,\
		    &mode) != 0) {
    show_usage();
    return 1;
  }

  build_billiard(&bil);
  build_bdry(b, k_hi, &bil, &bps); // use max k for bdry discr
  build_basis(k_hi, &bil, &bas);
  /*show_billiard_properties(&bil, k_hi); 
  show_colloc_properties(&bil, &bps, k_hi); 
  show_basis_properties(&bil, &bas, k_hi); 
   output geom...
  dump_billiard(&bil, head);
  dump_colloc(&bps, head);
  dump_basis(&bas, k_lo, head);
  */

  N = bas.N;
  S = dvector(0,N*N-1);
  T = dvector(0,N*N-1);
  X = dvector(0,N*N-1);
  evals = dvector(1,N);


  if (mode==0) { //........................ sweep ........................

    sprintf(name, "%s.gev", head);
    fp = fopen(name, "w");
    
    for (k=k_lo; k<=k_hi; k+=dk) {
      
      //printf("k=%g\n", k);
      build_dirichlet_S_and_T(&bil, &bps, &bas, k, S, T);
      
      rank = truncated_gen_eig_prob(N, 1e-14, T, S, evals, X);
      
      fprintf(fp, "%.15g %d ", k, rank);
      for (i=1; i<=N; ++i)
	fprintf(fp, "%.15g ", evals[i]);
      fprintf(fp, "\n");
    }
    fclose(fp);
    
  } else if (mode==1 || mode==-1) { //........ efunc at k_lo...............

    ne = 1;
    ks = dvector(1, ne);
    a = dmatrix(1, ne, 1, N);
    
    dx = 0.02;
    make_grid(&bil, &g, ne, dx);

    k = k_lo;

    if (mode==1) // dirichlet
      build_dirichlet_S_and_T(&bil, &bps, &bas, k, S, T);
    else // neumann
      build_neumann_S_and_T(&bil, &bps, &bas, k, S, T);
    
    rank = truncated_gen_eig_prob(N, 1e-14, T, S, evals, X);
    
    printf("ne=%d\n", ne);

    // extract ne states with largest eval...
    for (i=1; i<=ne; ++i) {
      ks[i] = evals[rank+1-i];
      for (j=1; j<=N; ++j)
	a[i][j] = X[j-1 + N*(rank-i)] / sqrt(ks[i]); // norm to S
    }
    eval_grid_vecs(&bil, &bas, k, a, &g, 0); // no clip
    // encode_grid(&bps, &g, ks, k, 1);
    save_grid(&g, ks, head);
  
    free_dvector(ks, 1, ne);
    free_dmatrix(a, 1, ne, 1, N);
    free_grid(&g);
  }

  return 0;
}
