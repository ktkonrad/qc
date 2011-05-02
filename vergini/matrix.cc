/* General matrix / linear algebra for VERGINI package
 *
 *	There should be no billiard- or basis-specific references in this code
 *
 * Oct 2003
 */ 

#include <stdio.h>
#include <math.h>

#include "billiard.h"
#include "basis.h"
#include "colloc.h"

#include "cxml.h"
#include "nrutil.h"
#include "useful.h"
#include "matrix.h"
#include "verb.h"


int truncated_gen_eig_prob(int N, double eps, double *T, double *S, \
			   double *gev, double *X)

/* Gen evals and largest eigenvec of (S - gev.T)x = 0.
 * Truncates basis to rank r, by selecting evals of T which are >= eps.
 * Expects T,S to be zero-offset f77 storage, T to be 
 * Returns r ordered gev's in gev (1-offset NumRec C style array),
 * and corresp r solution vectors in X in columns (f77 storage).
 * Destroys S and T!
 * Norm |x|^2 of any solution vector = 1/(corresponding eval of T).
 * Uses internal workspace O(N(N+10)).
 *
 * Alex Barnett 99/11/6
 */
{
  int i,j,r;
  double inv_sqrt_evT,*work,*temp, *evT, cutoff;

  /* LAPACK options...  */
  int lwork,info;
	
/* allocate... generally vectors are 1-indexed, matrices fortran 0-indexed */
  temp = dvector(0,N*N-1);
  evT = dvector(1,N);
  /* LAPACK workspace... */
  lwork = 80*N;
  work = dvector(0,lwork-1);


/* Diagonalize T... evecs of T returned in T, in columns. */
  dsyev_('V', 'U', N, T, N, evT+1, work, lwork, info);
  if (info!=0)
    fprintf(stderr, " truncated_gen_eig_prob: dsyev(T) info = %d\n",info);
  else if (lwork < (int)work[0]) {
    fprintf(stderr, \
	    "Sub-optimal workspace for dsyev (lwork = %d, optimum = %d)\n",\
	   lwork,(int)work[0]);
  }
  if (verb&0x02)
    printf("\t\tT diagonalised: eigenvalue range [%g,%g]\n", evT[1], evT[N]);
  /*for (j=1;j<=N;++j)
    printf("\t\t\tev[%d] = %g\n", j, evT[j]);*/

  cutoff = eps * evT[N];  /* cutoff now relative to largest eigval - 1/28/05 */

/* decide numerical rank of F ... all ev's above cutoff */
  r = N;
  while (evT[N+1-r] < cutoff) { // evT is 1-indexed
    --r;
  }
  if (verb&0x02)
    printf("\t\trank r=%d\n",r);

/* scale the evecs of T, clipping rank... T now = V.Lambda^{-1/2}, size N*r */
  for(j=1;j<=r;++j) {
    inv_sqrt_evT = 1/sqrt(evT[j+N-r]);
    for (i=1;i<=N;++i)
      T[i-1 + N*(j-1)] = inv_sqrt_evT * T[i-1 + N*(j+N-r-1)];
  }
  if (verb&0x02)
    printf("\t\tevecs of T rescaled\n");

/* compute T^T.S.T, storing it in S... */
  dgemm_('N','N',N,r,N,1.0, S,N, T,N, 0.0,temp,N);
  dgemm_('T','N',r,r,N,1.0, T,N, temp,N, 0.0,S,r);

  if (verb&0x02)
    printf("\t\tS transformed into trunc eigenbasis of T\n");

/* Diagonalize... with evecs in cols out into S. */
  dsyev_('V', 'U', r, S, r, gev+1, work, lwork, info);

  if (info!=0)
    printf("truncated_gen_eig_prob: dsyev(T^T.S.T) info = %d\n",info);
  else if (lwork < (int)work[0]) {
    printf("Sub-optimal workspace for dsyev (lwork = %d, optimum = %d)\n",\
	   lwork,(int)work[0]);
  }
  if (verb&0x02)
    printf("\t\tGeneralized eigenproblem done\n");

/* back-transform solutions to X = T.Ve, where Ve = evecs in cols */
  dgemm_('N','N',N,r,r,1.0, T,N, S,r, 0.0,X,N);


/* de-allocate... */
  free_dvector(temp,0,N*N-1);
  free_dvector(evT,1,N);
  free_dvector(work,0,lwork-1);

  return r;
}




// =============== Build DIRICHLET S, and exact T MATRICES ==================

void build_dirichlet_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			     double k, double *S, double *T)
  /* fills S and T, the area norm and perimeter norm matrices respectively.
     S is given by the BerryWilkinson84 bdry normalisation,
     correct for Dirichlet BC only.
     However, it is expected to be as good as the exact S for finding
     t-minima in a
     k-sweep, since the gen eigenvalue with smallest t/s is used.
     It is much (ie more than 2 times) faster to build than the exact S,
     and only uses up to first derivs.
     Alex Barnett 99/11/16
     10/30/03 modified for objects
     */

{
  int i,j,n, M = p->M, N = s->N;
  double ddx, ddy, xi, yi, norx, nory, srni;
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
      A[n] = srni*eval_basis(xi, yi, j, s);  // srni factor 1/23/05
      B[n] = (1.0/srni)*( norx*ddx + nory*ddy );
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

   Still untested, not sure if believe...
   
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




// --------------------------------------------------------------------------
void diag_scal(double *A, double *D, double *B, int M, int N, double alpha, \
	       double beta)
  // Premultiply M-by-N general matrix A by M-by-M diag matrix D:
  //     B = alpha.D.A + beta.A
  // D is given as an M-length vector. A, B assumed contiguously stored.
  // B may be the same storage space as A.
  // This is identical to scaling the rows of A by the entries in D.
  // No such BLAS routine exists!
  // barnett 8/9/04
{
  int i;

  if (B!=A)                  // compare pointers
    dcopy_(M*N, A, 1, B, 1);  // copy B <- A if need to
  for (i=1; i<=M; ++i)
    dscal_(N, alpha*D[i-1] + beta, B+i-1, M);   // rescale each row
}





// ============ BUILD ALEX (arbitrary-BC) S, and T MATRICES =================

void build_alex_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			double k, double *S, double *T, int weight)
  /* fills S and T, the area norm and perimeter norm matrices respectively.
     Uses my S formula derived by taking k_2->k_1 limit of unequal k formula,
     (thesis H.7), involving up to 2nd derivs.

     Alex Barnett 00/4/5, modified for objects 10/30/03
     7/19/04 added optional bdry weight on T, 8/9/04 added flag for this.
     Now T is correctly weighted by 1/rn if rn<0.
     Neumann option (read from bil object) for T, 1/21/05
     */

{
  int i,j,u,v,n, N=s->N, M=p->M;
  double val, ddx, ddy, dxx, dyy, dxy, xi, yi, norx, nory, sqrtw;
  double *A, *B, *C, *D, *W;

  // allocate rectangular matrices...
  A  = dvector(0,M*N-1);
  C  = dvector(0,M*N-1);
  D  = dvector(0,M*N-1);
  B  = dvector(0,M*N-1);
  // weight vector...
  W = dvector(0,M-1);
  if (verb)
    printf("\tbuild_alex_S_and_T: allocated\n");

  // fill 4 matrices used to build S...
  for (i=1;i<=M;++i) {
    xi = k*p->x[i]; // unitless coords.
    yi = k*p->y[i];
    norx = cos(p->na[i]);
    nory = sin(p->na[i]);
    W[i-1] = weight ? 1.0/p->rn[i] : 1.0;  // choice of weight funcs: 0,1
    for(j=1;j<=N;++j) {
      eval_basis_everything(val, ddx, ddy, dxx, dyy, dxy, xi, yi, j, s);
      n = i-1 + M*(j-1);
      A[n] = val;
      C[n] = xi*ddx + yi*ddy; // r.Del
      D[n] = k*( norx*ddx + nory*ddy ); // n.Del
      B[n] = k*( norx*(dxx*xi + dxy*yi) + nory*(dxy*xi + dyy*yi) ); // n.DelDel.r
    }
  }
  if (verb&0x02)
    printf("\t\tA,B,C filled\n");
  
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
  if (verb&0x02)
    printf("\t\tS done\n");

  if (l->neumann==0) {
  // T = (L/M) A^T.W.A ... (this costs little extra work above building S)
  diag_scal(A, W, B, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, l->Perim/M, A,M, B,M, 0.0,T,N);
  if (verb&0x02)
    printf("\t\tT done\n");
  } else {
  // T = (L/M) D^T.W.D ... since D contains normal derivs.
  diag_scal(D, W, B, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, l->Perim/M, D,M, B,M, 0.0,T,N);
  if (verb&0x02)
    printf("\t\tT (Neumann) done\n");
  }

  free_dvector(A,0,M*N-1);
  free_dvector(C,0,M*N-1);
  free_dvector(D,0,M*N-1);
  free_dvector(B,0,M*N-1);
  free_dvector(W,0,M-1);
}



// ============ BUILD Quasi (arbitrary-BC) S and T MATRICES =================

void build_quasi_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			double k, double *S, double *T, int weight)
  /* fills S and T, the area norm and (weighted) perimeter norm matrices
     respectively. Memory-efficient, reasonably.
     weight = 0: use w=1
              1: use w=1/rn
     Uses my formula from Quasi paper of July 2004, involving 1st derivs only
     barnett 8/9/04
   */

{
  int i,j,u,v,n, N=s->N, M=p->M;
  double val, ddx, ddy, xi, yi, norx, nory;
  double *A, *B, *C, *R, *W;

  // allocate rectangular matrices...
  A  = dvector(0,M*N-1);
  B  = dvector(0,M*N-1);
  C  = dvector(0,M*N-1);
  R = dvector(0,M-1);
  W = dvector(0,M-1);

  if (verb)
    printf("\tbuild_quasi_S_and_T: allocated\n");
  // fill 3 matrices used to build S...
  for (i=1; i<=M; ++i) {
    xi = k*p->x[i]; // unitless coords.
    yi = k*p->y[i];
    norx = cos(p->na[i]);
    nory = sin(p->na[i]);
    R[i-1] = p->rn[i];
    W[i-1] = weight ? 1.0/R[i-1] : 1.0;
    for(j=1; j<=N; ++j) {
      eval_basis_deriv(ddx, ddy, xi, yi, j, s);
      n = i-1 + M*(j-1);
      A[n] = xi*ddx + yi*ddy; // r.Del
      B[n] = k*( nory*ddx - norx*ddy ); // s.Del
      C[n] = k*( norx*ddx + nory*ddy ); // n.Del
    }
  }
  if (verb&0x02)
    printf("\t\tA,B,C filled\n");

  // S = A^T.C ...
  dgemm_('T','N',N,N,M, 1.0, A,M, C,M, 0.0,S,N);
  if (verb&0x02)
    printf("\t\tS = A^T.C done\n");
  // A is now free. Add S transpose to itself... surely must be a BLAS way?
  for (i=1; i<=N; ++i)
    for(j=1; j<=i; ++j) {
      S[u = i-1 + N*(j-1)] += S[v=j-1 + N*(i-1)];
      S[v] = S[u];
    }
  if (verb&0x02)
    printf("\t\tS symmetrized\n");

  // Add -B^T.R.B to S...
  diag_scal(B, R, A, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, -1.0, B,M, A,M, 1.0,S,N);
  if (verb&0x02)
    printf("\t\t-B^T.R.B added to S\n");
  // B is now free. Add -C^T.R.C to S...
  diag_scal(C, R, A, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, -1.0, C,M, A,M, 1.0,S,N);
  if (verb&0x02)
    printf("\t\t-C^T.R.B added to S\n");

  // fill value matrix...
  for (i=1; i<=M; ++i) {
    xi = k*p->x[i]; // unitless coords.
    yi = k*p->y[i];
    for(j=1; j<=N; ++j) {
      val = eval_basis(xi, yi, j, s);
      n = i-1 + M*(j-1);
      A[n] = val;
    }
  }
  if (verb&0x02)
    printf("\t\tvalue matrix A filled\n");
    
  // Add E.A^T.R.A to S...
  diag_scal(A, R, B, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, k*k, A,M, B,M, 1.0,S,N);
  if (verb&0x02)
    printf("\t\tE.A^T.R.A added to S\n");

  // scale S by (L/M)/(2 k^2)...
  dscal_(N*N, l->Perim/(2*M*k*k), S, 1);

  // T = (L/M) A^T.W.A ... (this costs little extra work above building S)
  diag_scal(A, W, B, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, l->Perim/M, A,M, B,M, 0.0,T,N);
  if (verb&0x02)
    printf("\t\tT done\n");

  free_dvector(A,0,M*N-1);
  free_dvector(B,0,M*N-1);
  free_dvector(C,0,M*N-1);
  free_dvector(R,0,M-1);
  free_dvector(W,0,M-1);
}



// ------------------------------------------------ diag -------------------
int diag(int N, double *A, double *evA)
/* Evals of square symm A using lapack.
 * Expects A to be zero-offset f77 storage.
 * Returns N ordered ev's in evA array, (1-offset NumRec C style).
 */
{
  double *work;  /* DSYEV options...  */
  int lwork,info;
	
  // allocate... generally vectors are 1-indexed, matrices fortran 0-indexed
  /* LAPACK workspace... */
  lwork = 80*N;
  work = dvector(0,lwork-1);
  
  // Eigenproblem... 'N' for no evecs, 'U' look in upper tri of A
  dsyev_('V', 'U', N, A, N, evA+1, work, lwork, info);
  
  if (info!=0)
    printf("eigenvals: dsyev info = %d\n",info);
  else if (lwork < (int)work[0]) {
    printf(\
	   "Sub-optimal workspace for dsyev (lwork = %d, optimum = %d)\n",\
	   lwork, (int)work[0]);
  }
  
  // de-allocate...
  free_dvector(work,0,lwork-1);
  return info;
}


void build_vergini_forms(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			double k, double *dFdk, double *F, int weight)
  /* wei=2 is always set to 1/r_n for now, wei=1 gives const 1.
   * Note that absval of r_n is used so still incorrect for nonstar-shaped.
   * barnett 9/9/04
   */
{
  int i, j, a, b, N = s->N, M = p->M;
  double *A, *B, srni, ddx, ddy, beta = 0.1;

  A = dvector(0, M*N-1);
  B = dvector(0, M*N-1);
  if (verb)
    printf("build_vergini_forms at k=%g...\n", k);

  // fill A:   note this could be sped up slightly by swapping ij looping order
  for (i=1;i<=M;++i) {
    srni = (weight==1) ? (beta + (1-beta)*sqrt(p->rn[i])) : sqrt(p->rn[i]);
    for (j=1;j<=N;++j)
      A[i-1 + M*(j-1)] = eval_basis(k*p->x[i], k*p->y[i], j, s) / srni;
  }
  if (verb)
    printf("\tA filled\n");

  // fill B with dA/dk :
  for (i=1;i<=M;++i) {
    srni = (weight==1) ? (beta + (1-beta)*sqrt(p->rn[i])) : sqrt(p->rn[i]);
    for (j=1;j<=N;++j) {
      eval_basis_deriv(ddx, ddy, k*p->x[i], k*p->y[i], j, s);
      B[i-1 + M*(j-1)] = (p->x[i]*ddx + p->y[i]*ddy) / srni;
    }
  }
  if (verb)
    printf("\tdA/dk filled\n");

  // build F = (L/M).A^T.A (does twice the reqd ops since F symm!):
  dgemm_('T','N',N,N,M,l->Perim/M, A,M, A,M, 0.0,F,N);
  if (verb)
    printf("\tF = A^T.A built\n");

  // build dF/dk = (L/M).A^T.dA/dk + Transpose:
  dgemm_('T','N',N,N,M,l->Perim/M, A,M, B,M, 0.0,dFdk,N);
  if (verb)
    printf("\tA^T.dA/dk built\n");
  // add its transpose to itself...  could replace by dgema BLAS someday.
  for (i=1;i<=N;++i)
    for(j=1;j<=i;++j) {
      dFdk[a=i-1 + N*(j-1)] += dFdk[b=j-1 + N*(i-1)];
      dFdk[b] = dFdk[a];
    }
  if (verb)
    printf("\tdF/dk = A^T.dA/dk + transpose built\n");

  free_dvector(A, 0, M*N-1);
  free_dvector(B, 0, M*N-1);
}



// -----------------------------= DIAG GENERAL QUADRATIC FORM =--------------
int diag_quad_form(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, double k, \
		   double *lambdas, double **a, int j_lo, int j_hi, \
		   int weight)
  /* return coeffs of eigenvectors in a, and eigenvalues in lambdas, of
     a quadratic form F. Only writes vectors and values from j_lo to j_hi.
     1/4/04 barnett
  */
{
  int i, j, N = s->N, M = p->M, info;
  double *A, *F, *evF, srni;
  A = dvector(0, M*N-1);
  evF = dvector(1, N);
  F = dvector(0, N*N-1);

  // fill A:
  for (i=1; i<=M; ++i) {
    srni = sqrt(p->rn[i]);
    for (j=1; j<=N; ++j)
      A[i-1 + M*(j-1)] = eval_basis(k*p->x[i], k*p->y[i], j, s) / srni;
  }
  if (verb)
    printf("\tA filled\n");
  // build F = (L/M).A^T.A (does twice the reqd ops since F symm!):
  dgemm_('T', 'N', N, N, M, l->Perim/M, A,M, A,M, 0.0, F, N);
  if (verb)
    printf("\tF= A^T.A built\n");

  info = diag(N, F, evF);
  if (verb)
    printf("\tF diagonalized, info=%d, eigval range=[%g,%g]\n", info, \
	   evF[1], evF[N]);

  // copy eigenvecs back to a array... j labels vector number, i entry.
  for (j=j_lo; j<=j_hi; ++j) {
    lambdas[j+1-j_lo] = evF[j];
    for (i=1; i<=N; ++i)
      a[j+1-j_lo][i] = F[i-1 + (j-1)*N];
  }

  free_dvector(A, 0, M*N-1);
  free_dvector(evF, 1, N);
  free_dvector(F, 0, N*N-1);

  return info;
}



// -----------------------------= DIAG GENERAL QUADRATIC FORM =--------------
int diag_quad_form_FG(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, double k, \
		      double *lambdas, double **a, int j_lo, int j_hi, \
		      int weight, double epsilon)
  /* return coeffs of eigenvectors in a, and eigenvalues in lambdas, of
     a quadratic form F. Only writes vectors and values from j_lo to j_hi.
     5/26/04 barnett
     8/9/04 added weight flag.
     9/9/04 negative wei gives vergini quad form
  */
{
  int i, j, jj, N = s->N, M = p->M, rank;
  double *G, *F, *X, *gev;
  FILE *fp;
  F = dvector(0, N*N-1);
  G = dvector(0, N*N-1);
  X = dvector(0, N*N-1);
  gev = dvector(1, N);

  if (weight>=0) {
  //build_dirichlet_S_and_T(l, p, s, k, G, F); // use for Dirichlet variants
  //build_alex_S_and_T(l, p, s, k, G, F, weight); // is actually bit faster
  build_alex_S_and_T(l, p, s, k, F, G, weight); // use G to truncate
  //build_quasi_S_and_T(l, p, s, k, G, F, weight); // 'better' version
  } else {
    build_vergini_forms(l, p, s, k, G, F, -weight);
  }
  if (verb)
    printf("\tF,G filled (wei=%d)\n", weight);

  // hack to output matrices F, G if needed...
  if (0) {
    fp = fopen("temp.fg", "w");
    for (i=1; i<=N*N; ++i)
      fprintf(fp, "%.16g %.16g\n", F[i-1], G[i-1]);
    fclose(fp);
  }



  rank = truncated_gen_eig_prob(N, epsilon, F, G, gev, X);

  if (verb)
    printf("\tF,G simultaneously diagonalized, eps=%g, rank=%d.\n", epsilon,\
	   rank);

  if (j_hi > rank)
    j_hi = rank;
  if (j_lo > rank)
    j_lo = rank;

  // copy eigenvecs back to a array... j labels vector number, i entry.
  for (j=j_lo; j<=j_hi; ++j) {
    jj = j;
    // jj = rank+1-j; // gets from top of spectrum rather than bottom jj=j;
    lambdas[j+1-j_lo] = gev[jj];
    for (i=1; i<=N; ++i)
      a[j+1-j_lo][i] = X[i-1 + (jj-1)*N];
  }

  free_dvector(gev, 1, N);
  free_dvector(F, 0, N*N-1);
  free_dvector(G, 0, N*N-1);
  free_dvector(X, 0, N*N-1);

  return rank;
}


 
// ------------------------------- VERGINI -------------------------------
int vergini(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, double k, double eps, \
	    double delta_lo, double delta_hi, double *k_mu, double **cf, \
	    int ne, double C4)
  /* Vergini method at wavenumber k.
   Gen evals of (dF/dk - lambda.F)x = 0, with F = A^T.A
   Returns only n_ev relevant states in the delta window, max ne of them.
   Uses O(2N(N+M)) storage.
   barnett 11/17/03
   11/20/03 added asymm delta.
   1/29/05 tried allowing r.n<0 to work correctly.
  */
{
  int i, j, a, b, n, rank, n_ev, n_found, N = s->N, M = p->M;
  double lam, delta, ddx, ddy, srni, norm, *A, *B, *F, *dFdk, *gev, *X, *W;

  A = dvector(0, M*N-1);
  B = dvector(0, M*N-1);
  F = dvector(0, N*N-1);
  dFdk = dvector(0, N*N-1);
  W = dvector(0,M-1);
  if (verb)
    printf("vergini at k=%g...\n", k);

  // fill A:   note this could be sped up slightly by swapping ij looping order
  for (i=1;i<=M;++i) {
    W[i-1] = 1.0/p->rn[i];    // choose weight 1/r.n
    //srni = sqrt(p->rn[i]);
    //srni = (rn[i] + 99*sqrt(rn[i]))/100.0;
    for (j=1;j<=N;++j)
      A[i-1 + M*(j-1)] = eval_basis(k*p->x[i], k*p->y[i], j, s); // / srni;
  }
  if (verb)
    printf("\tA filled\n");

  // F = (L/M) A^T.W.A ... (does twice the reqd ops since F symm!)
  diag_scal(A, W, B, M, N, 1.0, 0.0);
  dgemm_('T','N',N,N,M, l->Perim/M, A,M, B,M, 0.0,F,N);
  if (verb)
    printf("\tF = A^T.W.A built\n");

  // fill B with dA/dk :
  for (i=1;i<=M;++i) {
    //srni = sqrt(p->rn[i]);
    //srni = (rn[i] + 99*sqrt(rn[i]))/100.0;
    for (j=1;j<=N;++j) {
      eval_basis_deriv(ddx, ddy, k*p->x[i], k*p->y[i], j, s);
      B[i-1 + M*(j-1)] = (p->x[i]*ddx + p->y[i]*ddy); // / srni;
    }
  }
  if (verb)
    printf("\tdA/dk filled\n");


  // build dF/dk = (L/M).A^T.W.dA/dk + Transpose:
  diag_scal(A, W, A, M, N, 1.0, 0.0);      // A <- W.A
  dgemm_('T','N',N,N,M,l->Perim/M, A,M, B,M, 0.0,dFdk,N);
  if (verb)
    printf("\tA^T.W.dA/dk built\n");
  // add its transpose to itself...  could replace by dgema BLAS someday.
  for (i=1;i<=N;++i)
    for(j=1;j<=i;++j) {
      dFdk[a=i-1 + N*(j-1)] += dFdk[b=j-1 + N*(i-1)];
      dFdk[b] = dFdk[a];
    }
  if (verb)
    printf("\tdF/dk = A^T.W.dA/dk + transpose built\n");

  free_dvector(A, 0, M*N-1);
  free_dvector(B, 0, M*N-1);
  gev = dvector(1, N);    // generalized lambdas
  X = dvector(0, N*N-1);  // eigenvectors in columns, Fortran storage

  // solve singular gen eig prob...
  rank = truncated_gen_eig_prob(N, eps, F, dFdk, gev, X);
  if (verb)
    printf("\tgen eig prob done, eps = %g, rank = %d\n", eps, rank);

  // convert lambdas to wavenumbers (still in gev array)...
  for (i=1; i<=rank; ++i) {
    lam = (fabs(gev[i])<2.0) ? 2.0 : gev[i]; // set unsafe gev to safe val
    // C4, see ~/bdry/math/polynomial.nb for gen exp in lam: 5.3 is best guess
    delta = 2/lam - (2/k)/(lam*lam) - C4/(lam*lam*lam);
    gev[i] = k - delta;
  }

  // find how many states below... need to allow 0 as well!
  for(a=1;a<=rank && (gev[rank-a+1]>=k-delta_lo && gev[rank-a+1]<k); ++a)
    ;
  a -= 1;
  // find how many states above... need to allow 0 as well!
  for(b=1;b<=rank && (gev[b]<=k+delta_hi && gev[b]>=k); ++b)
    ;
  b -= 1;
  n_found = a+b;
  if (verb)
    printf("\t4th-order correction done (C4=%g): good states below %d, above %d.\n", C4, a, b);

  // copy good coeffs and k_mu into outgoing arrays, norm them according
  // to vergini method...
  n_ev = min(ne, n_found); // don't overrun output arrays
  if (n_ev > 0)
    for (i=1; i<=n_ev; ++i) {
      n = (i<=a) ? rank-a+i : i-a;
      k_mu[i] = gev[n];
      
      // normalise using 4th order expansion of f(k-k0),
      // (C4=-1 typ for qust)...
      delta = k_mu[i] - k;
      norm = sqrt( 2.0*delta*delta*(1 - delta/k_mu[i]) -
                   C4*delta*delta*delta*delta );
      for(j=1; j<=N; ++j)
        cf[i][j] = norm * X[j-1 + N*(n-1)];
    }

  free_dvector(F, 0, N*N-1);
  free_dvector(dFdk, 0, N*N-1);

  return n_found;
}


// ------------------------------ REMOVE SPURIOUS STATES --------------------
#define SPURIOUS_NRM_ERR 0.3
#define SPURIOUS_TEN 0.0003

int spurious(int ne, double *ks, double **cf, double **per, double **ngr, \
	     double *nrm, double *ten, int M, int N, int *ind)
  /* remove spurious states (defined by abnormally low nrm), reshuffling
     all data arrays accordingly, returning the new number of valid states, n.
     Also (1-indexed) indices of nonspurious states returned in ind(1...n)
     11/17/03
  */
{
  int c, i, j, oi;
  
  for (c=1, i=1; i<=ne; ++i)
    if (fabs(nrm[i]-1.0) < SPURIOUS_NRM_ERR && ten[i] < SPURIOUS_TEN)
      ind[c++] = i;
  --c;

  // important to work from bottom up so don't erase data until after used
  for (i=1; i<=c; ++i) {
    if (verb&0x02)
      printf("ind[%d] = %d.\n", i, ind[i]);
    oi = ind[i];
    ks[i] = ks[oi];
    ten[i] = ten[oi];
    nrm[i] = nrm[oi];
    for (j=1; j<=M; ++j) {
      per[i][j] = per[oi][j];
      ngr[i][j] = ngr[oi][j];
    }
    for (j=1; j<=N; ++j)
      cf[i][j] = cf[oi][j];
  }
  
  return c;
}

// ----------------------------- INHOMOGENEOUS SOLVE -----------------------
int inhomog(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, double k, \
		   double eps, double **cf, Basis_Set *s_inh, double **cf_inh)
  /* Solves for optimal (least-squares in L2 sense, with weight w)
     coefficients to represent (with basis s) the
     inhomogeneous field given by coeffs cf_inh (with basis s_inh).
     1/2/04 barnett
  */
{
  int i, j, N = s->N, M = p->M, rank, lwork = 50*max(M,N), info;
  double *work, *A, *A_copy, *b, *x, *sv, sqrtw;
  work = dvector(0, lwork);
  A = dvector(0, M*N-1);
  A_copy = dvector(0, M*N-1);
  b = dvector(0, max(N,M)-1);
  x = dvector(0, max(N,M)-1);
  sv = dvector(0, min(N,M)-1);

  // fill basis effect matrix A
  for (i=1; i<=M; ++i) {
    sqrtw = 1.0; // sqrt w weight factor, eg w=1/(p->rn[i]) for vergini
    for (j=1; j<=N; ++j) 
      A[i-1 + M*(j-1)] = A_copy[i-1 + M*(j-1)] = \
	eval_basis(k*p->x[i], k*p->y[i], j, s) * sqrtw;
  }
  if (verb)
    printf("\tA filled\n");

  // fill RHS vector b (and x)
  for (i=1; i<=M; ++i) {
    sqrtw = 1.0; // MUST match sqrtw used above!
    b[i-1] = eval_vec_spatial(p->x[i], p->y[i], k, cf_inh[1], s_inh) \
      * sqrtw;
    x[i-1] = b[i-1];
  }

  // solve lin sys (one RHS), x and A_copy are modified.
  dgelss_(M, N, 1, A_copy, M, x, max(N,M), sv, eps, rank, work, lwork, info);
  if (verb) {
    printf("\tLeast-squares (dgelss) done: info=%d, rank=%d\n", info, rank);
    printf("\tsing vals in range %g to %g\n", sv[0], sv[min(N,M)-1]);
  }

  for (j=1; j<=N; ++j) // copy answer from x back to coeffs of n=1 state.
    cf[1][j] = x[j-1];

  // Check associated LSQ error: first put b-Ax into b
  dgemv_('N', M, N, -1.0, A, M, x, 1, 1.0, b, 1);
  if (verb)
    printf("\trms (weighted) soln error = %g\n", dnrm2_(M, b, 1) / \
	   sqrt((double)M));
  
  free_dvector(A, 0, M*N-1);
  free_dvector(A_copy, 0, M*N-1);
  free_dvector(b, 0, max(N,M)-1);
  free_dvector(x, 0, max(N,M)-1);
  free_dvector(sv, 0, min(N,M)-1);

  return info;
}
