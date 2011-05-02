#ifndef VERGINI_COLLOC_H
#define VERGINI_COLLOC_H 1

#include "billiard.h"
#include "basis.h"

/* BILLIARD BDRY POINT SET OBJECT */
typedef struct Bdry_Pt_Set {
  int M;                         /* size */
  double *x,  *y;                /* locations */
  double *nx, *ny;               /* outward normal direction */
  double *na;                    /*                angle */
  double *rn, *rs;               // r dot n, r dot s
  double *alpha;                 // curvature
};

/* Routines provided by colloc: */
extern void free_bdry(Bdry_Pt_Set *p);
extern void alloc_bdry(Bdry_Pt_Set *p);
extern void build_bdry(double b, double k, Billiard *l, Bdry_Pt_Set *p);
extern double eval_bdry(double k, double *cf, Basis_Set *s,\
			double *per, double *ngr, double *tens,\
			Billiard *l, Bdry_Pt_Set *p, int norm_flag, int BCs, \
			int m_c, double &mass);
extern void dump_colloc(Bdry_Pt_Set *p, char *head);
extern void show_colloc_properties(Billiard *l, Bdry_Pt_Set *p, double k); 

#endif /* VERGINI_COLLOC_H */

