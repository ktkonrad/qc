/* BNC, Module COLLIDE, does 2d raytracing step
 *
 * barnett 2/20/04
 * included finite arcs and finite line segments, test routine, 6/30/06
 * To Do: * fix up circ_coll to handle from inside collisions, and have sense
 *        * each object return the distance along it so can find q of collision
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "collide.h"
#include "useful.h"


double arc_coll(double x, double y, double vx, double vy, double cx,
		double cy, double R2, double inv_R, double ax, double ay,
		double bx, double by, double *nx, double *ny)
/* Returns d distance to first hit with arc and if hit, writes nx,ny,
 * the `inward' normal at collision point. (cx,cy) is arc center, R radius.
 * Angle range is described by whether the direction vector lies between
 * (in CCW sense)
 * vectors a and b (which are passed in, and must be of length R).
 * Returns d=0 if no hit.
 * Minimal use of sqrt (once) and divide (never)
 */ 
{
  double px, py, dx = cx-x, dy = cy-y;
  double l = dx*vx + dy*vy;
  double xo = dx - l*vx, yo = dy - l*vy;
  double ro2 = xo*xo + yo*yo;

  if (ro2<R2) { /* direction line hits the circle at all */
    double h = sqrt(R2 - ro2);
    //printf("h = %g\n", h); 
    double d = l - h; /* the hit with lowest time value */
    //   printf("d = %g, SMALLDIST=%g\n", d, SMALL_DIST);
    if (d>SMALL_DIST) {
      //printf("d- = %f\n", d); 
      px = d*vx - dx; /* vector from center to collision pt */
      py = d*vy - dy;
      if ((bx-px)*(ay-py)-(by-py)*(ax-px)>=0.0) { /* check angle range */
	*nx = px * inv_R;        /* normalize, outward in circle sense */
	*ny = py * inv_R;
	return d;
      }
    }
    d = l + h; /* the hit with highest time value */
    //printf("d = %g, SMALLDIST=%g\n", d, SMALL_DIST);
    if (d>SMALL_DIST) {
      //printf("d+ = %f\n", d);
      px = d*vx - dx; /* vector from center to collision pt */
      py = d*vy - dy;
      if ((bx-px)*(ay-py)-(by-py)*(ax-px)>=0.0) { /* check angle range */
	*nx = -px * inv_R;        /* normalize, inward in circle sense */
	*ny = -py * inv_R;
	return d;
      }
    }
  }
  /* doesn't hit */
  return 0.0;
}


double lineseg_coll(double x, double y, double vx, double vy,
		    double x1, double y1, double x2, double y2,
		    double inv_l)
/* Given (x,y) start location and unit direction vector (vx,vy),
 * returns distance d to hit with the line segment between (x1,y1) and (x2,y2).
 * Needs inv_l = 1/segment length. Returns d=0 if no hit.
 * Note the line is "1-sided", ie its inward normal is fixed.
 *
 * barnett 6/03/06 based on 9/21/03, version avoiding sqrt. 7/26/06 1-sided.
 */ 
{
  double dx = x1-x, dy = y1-y, px = x2-x1, py = y2-y1;
  double b, c = py*vx - px*vy;

  if (c!=0.0) {
    double inv_c = 1.0/c;
    double d = (py*dx - px*dy) * inv_c;
    if (verb>1)
      fprintf(stderr, "p=(%g,%g) v=(%g,%g) c=%g d=%g b=%g\n", px, py, vx, vy, c, d, (vy*dx - vx*dy)*inv_c);
    if (d>0.0 && (b=(vy*dx - vx*dy)*inv_c)>=0.0 && b<=1.0)       // in segment
      return d;
  }	
  /* doesn't hit */
  return 0.0;
}  


double circ_coll(double x, double y, double vx, double vy, double cx,
		 double cy, double R2, double *nx, double *ny)
/* Returns distance to first hit with circle (center cx,cy, radius sqrt(R2)),
 * and if hit, writes nx,ny, the inward normal at collision point.
 * Returns 0 if no hit. Currently, qugrs only valid for theta>0, ie
 * bouncing off outside of circles, due to sign of normal here.
 * OUT OF DATE.
 */ 
{
  double dx = cx-x, dy = cy-y;
  double l = dx*vx + dy*vy;
  double xo = dx - l*vx, yo = dy - l*vy;
  double ro2 = xo*xo + yo*yo;

  if (ro2<R2) {
    /* hits, but only use if happens in at least small dist in future */
    double d = l - sqrt(R2 - ro2);
    if (d>0.0) {
      /* calc inward normal at impact... (outward for the circle) */
      *nx = d*vx - dx;
      *ny = d*vy - dy;
      double n = sqrt(*nx*(*nx) + *ny*(*ny));
      *nx /= n;
      *ny /= n;
      // printf("\t\tn=%g,%g\n", *nx, *ny);
      return d;
    }
  }
  /* doesn't hit */
  return 0.0;
}


double line_coll(double x, double y, double vx, double vy, double rox,
		 double roy, double px, double py, double *nx, double *ny)
/* Returns distance to hit with line (passing though point ro with
 * unit direction vector p), and writes nx,ny, the
 * inward normal at collision point. Returns 0 if no hit (parallel).
 * Note the line is "2-sided", ie sensitive to hits from both sides,
 * and the normal will be returned 'inward' for the appropriate side.
 *
 * barnett 9/21/03, version avoiding sqrt.   OUT OF DATE.
 */ 
{
  double dx = rox-x, dy = roy-y;
  double c = py*vx - px*vy;

  if (c!=0.0) {
    double d = (py*dx - px*dy) / c;
    if (d>0.0) {
      /* calc inward normal (for coming from 'right' side)... */
      *nx = py;
      *ny = -px;
      if (c>0.0) {
	*nx = -*nx;
	*ny = -*ny;
      }
      return d;
    }    
  }	
  /* doesn't hit */
  return 0.0;
}  


double obj_coll(double x, double y, double vx, double vy, struct Segment *s,
		double *nx, double *ny)
/* Calls above object collision routine appropriate for a Segment type
 */
{
  double d;
  switch (s->type) {
  case 0:
    d = lineseg_coll(x, y, vx, vy, s->x0, s->y0, s->x1, s->y1, s->inv_l);
    *nx = s->nx;
    *ny = s->ny;
    break;
  case 1:
    d = circ_coll(x, y, vx, vy, s->cx, s->cy, s->R2, nx, ny);
    break;
  case 2:
  case -2:
    d = arc_coll(x, y, vx, vy, s->cx, s->cy, s->R2, s->inv_R,  s->ax, s->ay,
		 s->bx, s->by, nx, ny);
    break;
  case 3:
    d = line_coll(x, y, vx, vy, s->x0, s->y0, s->x1-s->x0, s->y1-s->y0,
		  nx, ny);
    break;
  default:
    fprintf(stderr, "Bad segment type (%d) in obj_coll!\n", s->type);
    exit(1);
  }
  return d;
}

double one_bounce(double *x, double *y, double *vx, double *vy, int *obj, \
		  int ns, struct Segment *s)
/* Perform one bounce and update (x,y,vx,vy,obj).
   ns is # of segments, s is segment array.
   Returns distance travelled, or zero if no hit.
   barnett 7/26/06
*/
{
  int i;
  double d[MAX_N_SEGS], nx[MAX_N_SEGS], ny[MAX_N_SEGS];

  // check all segment collisions
  for (i=0; i<ns; ++i) {
    d[i] = obj_coll(*x, *y, *vx, *vy, s + i, nx + i, ny + i);
    if (verb>1)
      fprintf(stderr, "\td[%d] = %g\n", i, d[i]);
  }  

  // remove double-bounces on current obj - means need current obj passed in?
  if (d[*obj] < SMALL_DIST)
    d[*obj] = 0.0;
  *obj = ns;
  double mind = LARGE_DIST;                                // init min dist
  for (i=0; i<ns; ++i) {
    if (d[i] < mind && d[i]!=0.0) {
      mind = d[i];
      *obj = i;
    }
  }
  if (*obj==ns) {
    fprintf(stderr, "one_bounce: no collision! (x,y,vx,vy) = (%g,%g,%g,%g)\n",\
	    *x, *y, *vx, *vy);
    exit(1);
  }
  // advance to next launch coords...
  *x += mind*(*vx);
  *y += mind*(*vy);
  // reflect about (inwards) normal...
  double mx = nx[*obj], my = ny[*obj]; // the inwards normal
  double vt = mx*(*vy) - my*(*vx);
  double vn = -(*vx)*mx - (*vy)*my;
  *vx += 2*vn*mx;
  *vy += 2*vn*my;
  return mind;
}


