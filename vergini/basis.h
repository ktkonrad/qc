/* Basis set module for VERGINI
 *
 * barnett 9/10/03
 */

#ifndef VERGINI_BASIS_H
#define VERGINI_BASIS_H 1

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>                 // Gnu Sci Lib special funcs

#include "billiard.h"

#define MAX_N_BASIS_PARAMS 20

/* BASIS SET TYPES: NBASES assigned to number of complete basis set types... */
typedef enum
	{RPWS,\
	 RPWS_ODD_ODD,\
	 RPWS_EVEN_EVEN,\
	 RPWS_EPWS_ODD_ODD,\
	 QUST_PISTON1_ODD_ODD,\
	 VERG_STAD_ODD_ODD,\
	 OUTER_BIM_Y0,\
	 OUTER_BIM_Y0_ODD_ODD,\
	 OUTER_BIM_Y0_EVEN_EVEN,\
	 RPWS_REFL,\
	 FOURIER_BESSEL,\
	 NBASES}
	basis_t;

/* GLOBALS: basis set description strings... */
static char *basis_desc[] =
	{"Real plane waves",
	 "Real plane waves, odd-odd symm",
	 "Real plane waves, oven-even symm",
	 "Real PWs + simple evanescent PWs, o-o",
	 "qust Piston1: Real PWs + EPWs along qust straight edge, o-o",
	 "Vergini 1/4 stadium ideal R+EPWs, o-o",
	 "Outer BIM (Y_0 sources)",
	 "Outer BIM (Y_0 sources), o-o",
	 "Outer BIM (Y_0 sources), e-e",
	 "Real plane waves, reflected 2p/pi corner",
	 "Fourier-Bessel reentrant right-angle"};
static char *basis_tags[] =
	{"pw", "pwoo", "pwee", "sepwoo", "pepwoo", "vepwoo", "oyo", "oyooo", \
	 "oyoee", "pwre", "fb"};
static int basis_num_params[] =	{1,1,1,3,4,3,3,3,3,2,3};
static char *basis_param_desc[] =
	{"eta", "eta", "eta", "eta_rpw N_epw alpha", "eta_rpw N_epw A B", \
	 "eta_rpw N_epw beta", "eta kD corner_pack", "eta kD corner_pack", \
	 "eta kD corner_pack", "eta p", "eta start_ang_pi end_ang_pi"};


/* BASIS SCALING FUNCTION TYPES: NSFUNCS assigned to number of func types... */
typedef enum
	{RPW_RE, RPW_IM, RPW_ODD_ODD, RPW_EVEN_EVEN,
	 EPW_RE_ODD_ODD, EPW_IM_ODD_ODD,
	 YO, YO_ODD_ODD, YO_EVEN_EVEN, RPW_REFL, FB, NSFUNCS}
	sfunc_t;


/* BASIS SET OBJECT */
typedef struct Basis_Set {
  basis_t set_type;                  /* what type of basis set is this */
  double param[MAX_N_BASIS_PARAMS];  /* basis cmd line params */
  double k;                          // wavenumber for this basis set
  int N;                             /* size */
  sfunc_t *t;                        /* type array of each basis func */
  double *nx, *ny;                   /* wavevector info for each func */
  double *a, *b, *c;                 /* other arbitrary info, eg location */
  int *d;                            // other arbitrary info of integer type
};


// basis.cc provides functions...
extern double eval_basis(double x, double y, int j, Basis_Set *s);
extern void eval_basis_deriv(double &ddx, double &ddy, double x, double y,\
			     int j, Basis_Set *s);
extern void eval_basis_everything(double &val, double &ddx, double &ddy,\
				  double &dxx, double &dyy, double &dxy,\
				  double xi, double yi, int j, Basis_Set *s);
extern void dump_sfunc_geom(FILE *fp, int j, double k_typ, Basis_Set *s);
extern void alloc_basis(Basis_Set *s, int N);
extern void free_basis(Basis_Set *s);
extern int build_basis(double k_0, Billiard *l, Basis_Set *s);
extern void make_inhomog_basis(double k, Billiard *l, Basis_Set *s, double f,\
			       double kD);
extern double eval_vec_spatial(double x, double y, double k, double *coeff,\
			       Basis_Set *s);
extern void eval_vec_grad_spatial(double &ddx, double &ddy, double x, \
				  double y, double k, double *coeff, \
				  Basis_Set *s);
extern void eval_vec_everything_spatial(double &val, double &ddx, double &ddy,\
					double &dxx, double &dyy, double &dxy,\
					double x, double y, double k, \
					double *coeff, Basis_Set *s);
extern void eval_multi_vecs_spatial(double x, double y, double k, int nv, \
	double **coeff_matrix, double *out, Basis_Set *s);
extern int parse_basis(char *st, Basis_Set *s);
extern void dump_basis(Basis_Set *b, char *head);
extern void show_basis_properties(Billiard *l, Basis_Set *s);
extern void show_basis_usage();


#endif /* VERGINI_BASIS_H */




