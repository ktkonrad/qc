/* Point collocation package, for BDRY package.

   Barnett 9/08/03
*/

#include <math.h>

#include "nrutil.h"
#include "useful.h"
#include "colloc.h"
#include "basis.h"
#include "billiard.h"



void alloc_bdry(Bdry_Pt_Set *p)
 /* Allocate a bdry point set
    
    barnett 9/8/03
 */
{
  p->x = dvector(1, p->M);
  p->y = dvector(1, p->M);
  p->nx = dvector(1, p->M);
  p->ny = dvector(1, p->M);
  p->na = dvector(1, p->M);
  p->rn = dvector(1, p->M);
  p->rs = dvector(1, p->M);
  p->alpha = dvector(1, p->M);
}

void free_bdry(Bdry_Pt_Set *p)
 /* Deallocate a bdry point set
    
    barnett 9/8/03
 */
{
  free_dvector(p->x, 1, p->M);
  free_dvector(p->y, 1, p->M);
  free_dvector(p->nx, 1, p->M);
  free_dvector(p->ny, 1, p->M);
  free_dvector(p->na, 1, p->M);
  free_dvector(p->rn, 1, p->M);
  free_dvector(p->rs, 1, p->M);
  free_dvector(p->alpha, 1, p->M);
}


/* ============================ BUILD BOUNDARY ========================== */

void build_bdry(double b, double k, Billiard *l, Bdry_Pt_Set *p)
/*
 * Fills the bdry point struct p, by calculating fractional dist along
 * perim, and sending to desc_by_perim. Same for all billiard types.
 * Having desc_by_perim separate allows perimeter to be measured linearly.
 * Matching point index runs j = 1 to M, where M is set by b, the number of
 * mathing points per wavelength (lambda = 2.pi/k) around perimeter.
 * If first argument < 0, abs val is used to fix M, the number of samples
 *
 * Altered to many billiard TYPES. v1.2
 * Specialised to building matching boundary points only, v1.4
 * Added computation of r_n array for Vergini method, v1.7 99/8/25
 * Readded normal angle array, v1.7 99/10/14
 * Moved to colloc, converted to handle objects, get M from b, 9/10/03.
 */
{
  int j, M;

  /* set M using b... */
  if (b>0)
    M = round_and_clip(b*k*l->Perim / (2*PI), "M in build_bdry");
  else
    M = round_and_clip(-b, "M in build_bdry");
  p->M = M;

  alloc_bdry(p); 

  /* M matching points... r=(x,y) */
  for (j=1; j<=M;++j)
    desc_by_perim((j-0.5)/(double)M, &(p->x[j]), &(p->y[j]), &(p->na[j]), \
		  &(p->alpha[j]), l);
  
  /* use angle of normal to compute n_x, n_y, r.n, r.s ... */
  for (j=1; j<=M;++j) {
    p->nx[j] = cos(p->na[j]);
    p->ny[j] = sin(p->na[j]);
    p->rn[j] = p->x[j] * p->nx[j] + p->y[j] * p->ny[j];
    p->rs[j] = -p->x[j] * p->ny[j] + p->y[j] * p->nx[j];
  }
}

  

// ============ eval and norm one state on bdry, general BCs =================

double eval_bdry(double k, double *cf, Basis_Set *s, double *per, double *ngr,\
		 double *tens, Billiard *l, Bdry_Pt_Set *p, int norm_flag,\
		 int BCs, int m_c, double &mass)
  /*
    Fill per (value) and ngr (normal gradient) for a single state, using
    the coefficient vector cf, with basis set s and colloc point set p.
    per and ngr need to have been allocated previously.
    Returns the gen-BC area norm, returns tension (integral of per) in tens.
    If norm_flag set, normalizes cf, per, ngr, tens, mass corresp to area
    norm of 1.
    BCs = 0 : Dirichlet BCs assumed (BerryWilkinson84 area norm formula).
	  1 : General BCs (uses improved 1st-deriv version).
          2 : General BCs (uses my version of Boasman94 formula, 2nd derivs).
    If m_c!=0, output the area norm accounted for by bdry pts m st m_c<m<=M,
    to the pass-by-reference variable mass.

    Questions:
    * Is there any point in preserving the dnr (2nd deriv) array ?
    
    9/8/03 restarted barnett, 1/16/04 modified for gen BCs
    1/17/04 mass splitting arguments
    3/22/04 new 1st-deriv gen-BC overlap formula
  */
{
  int j;
  double val, ddx, ddy, ddn, dxx, dyy, dxy, dnr, ddr, norm, fac, x, y, nx, ny;
  double sq, sumsq = 0.0;     // sum used to approximate integral.
  double sumsq_c = 0.0;   // the sum part lying in m>m_c
  *tens = 0.0;

  for (j=1; j<=p->M; ++j) {
    x = p->x[j];
    y = p->y[j];
    nx = p->nx[j];
    ny = p->ny[j];
    if (BCs==0 || BCs==1) {
      val = eval_vec_spatial(x, y, k, cf, s);
      eval_vec_grad_spatial(ddx, ddy, x, y, k, cf, s);
    } else
      eval_vec_everything_spatial(val, ddx, ddy, dxx, dyy, dxy, x, y, k, cf, \
				  s);
    per[j] = val;
    *tens += val*val;
    ddn = ddx*nx + ddy*ny;
    ngr[j] = ddn;
    if (BCs==0) // BW84 formula
      sq = p->rn[j] * ddn*ddn;   // sum r_n . d_n(psi)^2 at bdry pts
    else if (BCs==1) { // improved 1st-deriv gen-BC formula
      ddr = ddx*x + ddy*y;
      sq = p->rn[j] * (k*k*val*val - ddx*ddx - ddy*ddy) + 2*ddr*ddn;
    } else { // Barnett version of Boasman 2nd-derivs formula 
      ddr = ddx*x + ddy*y;
      dnr = nx*(dxx*x + dxy*y) + ny*(dxy*x + dyy*y); // n.(DelDelpsi).r
      sq = (ddr - val)*ddn - val*dnr; // Barnett area overlap formula
    }
    if (m_c>0 && j>m_c)
      sumsq_c += sq;  // add to post-m_c part (contrib to mass)
    else
      sumsq += sq; // add to standard part
  }

  /* rescale both to be integral around perimeter... */
  *tens *= l->Perim / p->M;                  // sum->integral conversion factor
  norm = (sumsq+sumsq_c) * l->Perim / (p->M * 2*k*k);  // here prefactor too
  mass = sumsq_c * l->Perim / (p->M * 2*k*k); // contrib to post-m_c part
  
  fac = 1.0/sqrt(norm); // will be faster to multiply rather than divide

  if (norm_flag) {
    /* adjust so norm = 1... */
    for (j=1; j<=p->M; ++j) {
      per[j] *= fac;
      ngr[j] *= fac;
    }
    for (j=1; j<=s->N; ++j)
      cf[j] *= fac;
    *tens /= norm;
    mass /= norm;
  }

  return norm; /* note _original_ norm is returned when norm_flag = 1 */
}




void dump_colloc(Bdry_Pt_Set *p, char *head)
{
  FILE *fp;
  int i;
  char mname[LEN]; 

  sprintf(mname, "%s.m_pts", head);
  if ((fp = fopen(mname, "w")) == NULL)
    printf("file: %s, write error...\n", mname);
  
  fprintf(fp,"# %s from vergini\n\n# gnuplot vec: x y dx dy\n\n", mname);
  
  for(i=1;i<=p->M;++i)
    fprintf(fp,"%g %g %g %g\n", p->x[i], p->y[i], \
	    cos(p->na[i]), sin(p->na[i]));
  
  fprintf(fp,"\n");
  fclose(fp);
}


void show_colloc_properties(Billiard *l, Bdry_Pt_Set *p, double k)
  /* Collocation set properties, including k- and billiard-dependent ones
   */
{
  printf("	Collocation points (M=%d) per wavelength: b = %f\n",\
	 p->M, (2*PI*p->M)/(k*l->Perim));
}
