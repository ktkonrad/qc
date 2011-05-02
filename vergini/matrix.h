/* Linear algebra, solvers, and matrix fillers, for VERGINI package
 *
 * 10/30/03 repackaged, Barnett
 */

#ifndef VERGINI_MATRIX_H
#define VERGINI_MATRIX_H 1

// matrix.cc provides
extern int truncated_gen_eig_prob(int N, double eps, double *T, double *S, \
				  double *gev, double *x);

extern void build_dirichlet_S_and_T(Billiard *l, Bdry_Pt_Set *p, \
				    Basis_Set *s, double k, double *S, \
				    double *T);
extern void build_neumann_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
				  double k, double *S, double *T);
extern void build_alex_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			       double k, double *S, double *T, int weight);
extern void build_quasi_S_and_T(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			       double k, double *S, double *T, int weight);
extern int diag(int N, double *A, double *evA);
extern int diag_quad_form(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			  double k, double *lambdas, double **a, int j_lo, \
			  int j_hi, int weight);
extern int diag_quad_form_FG(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, \
			  double k, double *lambdas, double **a, int j_lo, \
			  int j_hi, int weight, double epsilon);
extern int vergini(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, double k, \
		   double eps, double delta_lo, double delta_hi, \
		   double *k_mu, double **cf, int ne, double C4);
extern int spurious(int ne, double *ks, double **cf, double **per, \
		    double **ngr, double *nrm, double *ten, int M, int N, \
		    int* idx);
extern int inhomog(Billiard *l, Bdry_Pt_Set *p, Basis_Set *s, double k,\
		   double eps, double **cf, Basis_Set *s_inh, double **cf_inh);
#endif // VERGINI_MATRIX_H
