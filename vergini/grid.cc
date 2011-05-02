/*      Grid module for VERGINI package.
 *      Spatial 2D grid format, writing to viewer format.
 *      Also 1D boundary function grids.
 *
 *      Barnett 10/30/03
 */

#include <stdio.h>

#include "useful.h"
#include "nrutil.h"
#include "billiard.h"
#include "colloc.h"
#include "basis.h"
#include "grid.h"
#include "verb.h"

#define TINY 1e-6

int make_grid(Billiard *l, Grid *g, int ne, double dx, int mask, \
	      double xl, double xh, double yl, double yh)
  /* decide size and alloc a spatial grid array for wavefuncs, using billiard
   * bounding box. (ensure it's all covered).
   * mask>0 makes weight mask be generated with mask*mask grid subsampling,
   * not efficient, but still fast.
   * For mask=1 this has usual crude clip effect (wm=1 iff
   * center of cell is within billiard).
   * Coord offset (for bbox=0): Ctr of leftmost cell is at leftmost edge
   * of billiard - this wastes half a cell of sampling, but is consistent
   * with low.cc
   *
   * Barnett 10/30/03
   * 1/2/04 added bbox
   */
{
  int i, j, a, b;
  double w, cx, cy, off, hx, hy;

  g->ne = ne;
  g->dx = dx;
  if (xl-xh == 0.0) {  // use billiard bounding box
    g->ox = l->xl;
    g->oy = l->yl;
    hx = l->xh;
    hy = l->yh;
  } else {             // use user-specified box
    g->ox = xl;
    g->oy = yl;
    hx = xh;
    hy = yh;
  }
  g->nx = (int)((hx - g->ox)/dx + 1 - TINY);
  g->ny = (int)((hy - g->oy)/dx + 1 - TINY);
  g->g = d3tensor(0, g->nx+1, 0, g->ny+1, 1, g->ne);
  g->mask = mask;

  if (mask>0) {
    if (verb)
      printf("generating weight mask (mask=%d)...\n", mask);
    g->wm = dmatrix(0, g->nx+1, 0, g->ny+1);
    off = dx*(1-mask)/(2*mask); // subsampling offset (<=0)
    for (j=0,cy=g->oy+off; j<=g->ny; ++j,cy+=g->dx)
      for (i=0,cx=g->ox+off; i<=g->nx; ++i,cx+=g->dx) {
	w = 0.0;
	for (b=0; b<mask; ++b)
	  for (a=0; a<mask; ++a)
	    if (inside_billiard(cx+a*dx/mask, cy+b*dx/mask, l))
	      w += 1.0;
	g->wm[i][j] = w/(mask*mask);
      }
    if (verb)
      printf("done.\n");
  }

  return 0;
}  


// ...........................................................................
void free_grid(Grid *g)
{
  free_d3tensor(g->g, 0, g->nx+1, 0, g->ny+1, 1, g->ne);
  if (g->mask)
    free_dmatrix(g->wm, 0, g->nx+1, 0, g->ny+1); 
}


// ...........................................................................
void eval_grid_vecs(Billiard *l, Basis_Set *s, double k, double **a, \
		    Grid *g, int ne, int noff)
/* Fill grid object g with ne 2d grids of wavefuncs whose coeffs given by a.
 * Start at grid number noff+1 (noff is offset).
 * The basis set, and its k, stays fixed for these ne wavefuncs => fast.
 * Use g->mask as in eval_grid_vecs_each_k()
 * Barnett 10/30/03
 * 11/20/03 weight mask usage changed
 * 1/2/04 ne, noff added.
 */
{
  int i, j, n;
  double cx, cy;

  if (verb)
    printf("eval_grid_vecs: ne=%d (nx=%d by ny=%d), mask=%d, noff=%d\n", \
	   ne, g->nx, g->ny, g->mask, noff);

  for (j=0,cy=g->oy; j<=g->ny; ++j,cy+=g->dx)
    for (i=0,cx=g->ox; i<=g->nx; ++i,cx+=g->dx)
      if (g->mask==0 || g->wm[i][j]!=0.0)
	eval_multi_vecs_spatial(cx, cy, k, ne, a, g->g[i][j] + noff, s);
      else
	for (n=1; n<=ne; ++n)
	  g->g[i][j][n+noff] = 0.0;
  
  if (verb)
    printf("done\n");
}



// ...........................................................................
void eval_grid_vecs_each_k(Billiard *l, Basis_Set *s, double *ks, double *kos,\
			   double **a, Grid *g, int ne, int noff)
/* Fill grid object g with ne 2d grids of wavefuncs whose coeffs given by a.
 * Start at grid number noff+1 (noff is offset).
 * A single scaling basis set is used, rescaled to ks wavenumber array.
 * g->mask = 0: evaluate over whole domain
 * g->mask >= 1: eval only where weight mask nonzero
 * Barnett 11/17/03
 * 11/20/03 definition of wm usage changed (no fractional weighting)
 * 1/2/04 ne, noff added.
 */
{
  int i, j, n;
  double cx, cy;

  if (verb)
    printf("eval_grid_vecs_each_k: ne=%d (%d by %d), mask=%d, noff=%d\n", \
	   ne, g->nx, g->ny, g->mask, noff);

  for (j=0,cy=g->oy; j<=g->ny; ++j,cy+=g->dx) {
    if (verb && j%10==0)
      printf("%d%%.\n", (int)(100*j/(double)g->ny + 0.5));
    for (i=0,cx=g->ox; i<=g->nx; ++i,cx+=g->dx)
      if (g->mask==0 || g->wm[i][j]!=0.0) {
	// evaluate wavefuncs there (weight is 1)
	for (n=1; n<=ne; ++n)
	  g->g[i][j][n+noff] = eval_vec_spatial(cx, cy, ks[n], a[n], s);
      } else {
	// put zeros there
	for (n=1; n<=ne; ++n)
	  g->g[i][j][n+noff] = 0.0;
      }
  }
  if (verb)
    printf("done.\n");
}


void encode_grid(Bdry_Pt_Set *p, Grid *g, double *ks, double k, \
	       int rescale_flag)
/* Set pixels in a grid object to NaN at the locations of the collocation
 * points.
 * rescale_flag = 0: don't rescale points as ks array values change
 *                1: do.
 *
 * Barnett 10/30/03
 * 1/3/04 debugged rescaling
 */
{
  int j, i, ix, iy;
  double x, y, sc = 1.0;

  for (j=1; j<=g->ne; ++j) {
    if (rescale_flag)
      sc = ks[j]/k;
    for (i=1; i<=p->M; ++i) {
      // find nearest integer grid location
      ix = (int)((sc*p->x[i] - g->ox)/g->dx + 0.5);
      iy = (int)((sc*p->y[i] - g->oy)/g->dx + 0.5);
      if (ix>=0 && ix<=g->nx && iy>=0 && iy<=g->ny)
	// you're within grid: set data to NaN
	g->g[ix][iy][j] = 1.0/0.0;
    }
  }
}



// ----------------------------------- SAVE ------------------------------
int save_grid(Grid *g, double *ks, char *head)
  /* Output grid in viewer binary format.
   */
{
  FILE *fp;
  int i, j, n, sx = g->nx+1, sy = g->ny+1; // array sizes are (nx+1,ny+1)
  float *v;                                // float vectors for output
  char name[LEN];

  sprintf(name, "%s.sta_bin", head);
  if ((fp = fopen(name, "w")) == NULL) {
    fprintf(stderr, "file: %s, write error...\n", name);
    return 1;
  }

  v = vector(0, g->ne-1);
 
  if (verb)
    fprintf(stderr, "save_grid...\n");
  fprintf(fp,"bin\n"); // tells viewer it's a binary file
  fwrite(&(g->ne), 4,1,fp);
  fwrite(&sx, 4,1,fp);
  fwrite(&sy, 4,1,fp);
  fwrite(ks+1, 8, g->ne, fp); // write k-values (using 1-offset array)
  // write arrays...
  for (j=0; j<sy; ++j)
    for (i=0; i<sx; ++i) {
      // copy double to float vector...
      for (n=1; n<=g->ne; ++n)
	v[n-1] = (float)g->g[i][j][n];
      fwrite(v, 4, g->ne, fp);
    }
  fclose(fp);
  if (verb)
    fprintf(stderr, "done.\n");

  free_vector(v, 0, g->ne-1);
  return 0;
}


// ----------------------------------- SAVE MASK----------------------------
int save_mask(Billiard *l, Grid *g, char *head,\
	      double xl, double xh, double yl, double yh)
  /* Output mask grid in viewer binary format, like a single state.
     Temporary storage for 1 sta array needed.
     11/20/03
   */
{
  int i, j, sx = g->nx+1, sy = g->ny+1, r;
  char maskhead[LEN];
  double fake_k = (double)g->mask; // write int mask as "k" val
  Grid t; // temp grid object

  make_grid(l, &t, 1, g->dx, 0, xl, xh, yl, yh);

  // copy mask from g into data in t... (note ne argument is 1-indexed)
  for (j=0; j<sy; ++j)
    for (i=0; i<sx; ++i)
      if (g->mask>0)
	t.g[i][j][1] = g->wm[i][j];
      else
	t.g[i][j][1] = 1.0;  // g->wm will not have been allocated if mask=0

  sprintf(maskhead, "%s.mask", head); // decide new file head
  r = save_grid(&t, (&fake_k)-1, maskhead); // k array offset by 1 
  free_grid(&t);
  return r;
}


// --------------------------SAVE 1D SAMPLED FUNCS --------------------------
int save_oned(Billiard *l, double **a, double *vals, int ne, int M, \
	       double ds, double dsoff, char *head, char *suffix)
  /* write multiple 1d arrays in viewer format, including param1=Perim,
   * param2=Area.
   * 11/17/03
   */
{
  int i, j;
  char name[LEN];
  FILE *fp;

  sprintf(name, "%s.%s", head, suffix);
  if ((fp = fopen(name, "w")) == NULL) {
    fprintf(stderr, "file: %s, open for write error...\n", name);
    return 1;
  }
   
  /* viewer-file header... set so x-coord gives basis number 1...N */
  fprintf(fp, "1 %d %d %.16g %.16g\n2 %.16g %.16g\n\n", M, ne, dsoff*ds, ds,\
	  l->Perim, l->Area);
  for(j=1; j<=ne; ++j) {
    fprintf(fp, "%d %.16g\n", j, vals[j]);
    for (i=1;i<=M;++i)
      fprintf(fp, "%.16g ", a[j][i]);
    fprintf(fp,"\n");
  }
  if (verb)
    printf("1d dump (%s) written.\n", name);  
  fclose(fp);
  return 0;
}
 
