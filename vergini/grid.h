/* Grid module for VERGINI
 *
 * barnett 10/30/03
 * 11/16/03 added bdry func save
 */

#ifndef VERGINI_GRID_H
#define VERGINI_GRID_H 1

#include <stdio.h>

#include "billiard.h"
#include "colloc.h"

// GRID OBJECT
typedef struct Grid {
  double ***g;                                // the data (multiple grids)
  int ne;                                     // number of wavefuncs
  int nx, ny;                                 // size (0...nx, 0...ny)
  double dx;                                  // physical grid spacing
  double ox, oy;                              // physical origin location
  int mask;                                   // weight mask subsampling
  double **wm;                                // weight mask data (0<=wm<=1)
};


// grid.cc provides functions:
extern int make_grid(Billiard *l, Grid *g, int ne, double dx, int mask,\
		     double xl, double xh, double yl, double yh);
extern void free_grid(Grid *g);
extern void eval_grid_vecs(Billiard *l, Basis_Set *s, double k, double **a, \
			   Grid *g, int ne, int noff);
extern void eval_grid_vecs_each_k(Billiard *l, Basis_Set *s, double *ks, \
				  double *kos, double **a, Grid *g, int ne, \
				  int noff);
extern void encode_grid(Bdry_Pt_Set *p, Grid *g, double *ks, double k, \
	       int rescale_flag);
extern int save_grid(Grid *g, double *ks, char *head);
extern int save_mask(Billiard *l, Grid *g, char *head, \
		     double xl, double xh, double yl, double yh);
extern int save_oned(Billiard *l, double **a, double *vals, int ne, int M,\
		      double ds, double dsoff, char *head, char *suffix);

#endif // VERGINI_GRID_H
