/*	Billiard type specific code for VERGINI package: continuous
 *      definitions of billiard types.
 *
 *	Alex Barnett 99/10/14
 *      9/10/03 made Billiard object.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "useful.h"
#include "nrutil.h"
#include "billiard.h"
#include "verb.h"


double perimeter(Billiard *l)
  /* Excludes RAD_FUNC types */
{
  double alpha,R,R1,R2,R3,R4,a,b,t1,t2,t3,t4;

  switch(l->type) {
    
  case STADIUM:
    // stadium: start at bot-R corner, straight side len = alpha*radius
    alpha = l->param[0];
    return 2.0*(PI + alpha);
    
  case QU_STADIUM:
    // quarter stadium: center at origin, straight side len = alpha
    alpha = l->param[0];
    return (PI + alpha)/2;
    
  case QU_STADIUM_ALL:
    // quarter stadium with full perim: center at origin, straight side len = alpha
    alpha = l->param[0];
    return 2.0 + alpha + PI/2;
    
  case DESYM_SINAI:
    // pi/4 Sinai: center at (1,0)
    R = l->param[0];
    return sqrt(2.0) + R*(PI/4 - 1.0);
    
  case DESYM_SINAI_ALL:
    // pi/4 Sinai with full perim: center at (1,0)
    R = l->param[0];
    return sqrt(2.0) + R*(PI/4 - 2.0) + 2.0;
    
  case GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    return 4*(R1*t1 + R2*t2);
    
  case GEN_ASYM_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    t3 = l->param[3];
    t4 = l->param[4];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    R3 = a/sin(t3);
    R4 = 1/sin(t4);
    return 2*(R1*t1 + R2*t2 + R3*t3 + R4*t4);
    
  case QU_GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    return R1*t1 + R2*t2;

  case LINE:
    a = l->param[2] - l->param[0];  // delta x
    b = l->param[3] - l->param[1];  // delta y
    return sqrt(a*a + b*b);
    
  case ELL:
    return 2;

  case HALF_MUSH:
    a = l->param[0];
    b = l->param[1];
    return a + b/2 + (1+b/2)*(1+PI/2);

  default:
    fprintf(stderr,"Bad billiard type in perimeter!\n");
  }
  
  return 0.0;
}


double area(Billiard *l)
  /* Excludes RAD_FUNC types */
{
  double alpha,a,b,R,R1,R2,R3,R4,t1,t2,t3,t4;
  
  switch(l->type) {
    
  case STADIUM:
    // stadium: start at bot-R corner, straight side len = alpha*radius
    alpha = l->param[0];
    return PI + 2*alpha;
    
  case QU_STADIUM:
  case QU_STADIUM_ALL:
    // quarter stadium: center at origin, straight side len = alpha
    alpha = l->param[0];
    return (PI + 2*alpha)/4;
    
  case DESYM_SINAI:
  case DESYM_SINAI_ALL:
    // pi/4 Sinai: center at (1,0)
    R = l->param[0];
    return (4.0 - PI*R*R)/8;
    
  case GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    return 4*a - R1*R1*(2*t1 - sin(2*t1)) - R2*R2*(2*t2 - sin(2*t2));
    
  case GEN_ASYM_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    t3 = l->param[3];
    t4 = l->param[4];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    R3 = a/sin(t3);
    R4 = 1/sin(t4);
    return 4*a - R1*R1*(2*t1 - sin(2*t1))/2 - R2*R2*(2*t2 - sin(2*t2))/2 -\
      R3*R3*(2*t3 - sin(2*t3))/2 - R4*R4*(2*t4 - sin(2*t4))/2;


  case QU_GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    return a - R1*R1*(2*t1 - sin(2*t1))/4 - R2*R2*(2*t2 - sin(2*t2))/4;
    
  case LINE:
    return 0.0; // line has no area since not closed.

  case ELL:
    return 0.75;

  case HALF_MUSH:
    a = l->param[0];
    b = l->param[1];
    return a*b/2 + (PI/4)*SQ(1+b/2);

  default:
    fprintf(stderr,"Bad billiard type in area!\n");
  }
  
  return 0.0;
}


int inside_billiard(double x, double y, Billiard *l)
/*
 * returns true if x,y is inside the current billiard. (ignoring deformation).
 * Included translation 20/8/01
 * 11/17/03 Included desymmetrized testing if in 1st quadrant, for weight-mask
 * 12/13/03 qu_rad_func type
 */
{
  double d,Rsq,a,b,t1,t2,t3,t4,R1,R2,R3,R4,xc,yc,xc2,xc4,yc1,yc3;
  int inside=0, inC1, inC2, inC3, inC4;

  // translation of x,y means default billiard locations used now...
  x -= l->trans_x;
  y -= l->trans_y;

  switch(l->type) {
    
  case STADIUM:
  case QU_STADIUM:
    d = l->param[0]/2;
    if (x>d) {
      inside = (((x-d)*(x-d) + y*y) < 1.0);
    } else if (x>-d) {
      inside = (y<1.0 && y>-1.0);
    } else {
      inside = (((x+d)*(x+d) + y*y) < 1.0);
    }
    break;
    
  case QU_STADIUM_ALL:
    d = l->param[0]/2;
    if (y>=0.0)
      if (x>d) {
	inside = (((x-d)*(x-d) + y*y) < 1.0);
      } else if (x>0.0) {
	inside = y<1.0;
      }
    break;
    
  case DESYM_SINAI:
  case DESYM_SINAI_ALL:
    // pi/4 Sinai: center at (1,0)
    Rsq = l->param[0]*l->param[0];
    if (1.0-y > x)
      inside = (((x-1.0)*(x-1.0) + y*y) > Rsq);
    break;
    
  case GEN_ASYM_SINAI: // need to debug: still cases where gets it wrong!
    // eg, when have some theta<-pi/2.    Only important for wild shapes.
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    t3 = l->param[3];
    t4 = l->param[4];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    R3 = a/sin(t3);
    R4 = 1/sin(t4);
    xc2 = -a - R2*cos(t2);
    yc1 = -1.0 - R1*cos(t1);
    xc4 = a + R4*cos(t4);
    yc3 = 1.0 + R3*cos(t3);
    inside = (fabs(x)<a && fabs(y)<1.0); /* in Rect */
    inC1 = ((y-yc1)*(y-yc1)+x*x < R1*R1); /* in C1 */
    inC2 = ((x-xc2)*(x-xc2)+y*y < R2*R2); /* in C2 */
    inC3 = ((y-yc3)*(y-yc3)+x*x < R3*R3); /* in C3 */
    inC4 = ((x-xc4)*(x-xc4)+y*y < R4*R4); /* in C4 */
    if (t1<0.0)
      inside = inside || (y<-1.0 && inC1);
    else
      inside = inside && !inC1;
    if (t3<0.0)
      inside = inside || (y>1.0 && inC3);
    else
      inside = inside && !inC3;
    if (t2<0.0)
      inside = inside || (x<-a && inC2);
    else
      inside = inside && !inC2;
    if (t4<0.0)
      inside = inside || (x>a && inC4);
    else
      inside = inside && !inC4;
    break;
    
  case GEN_RECT_SINAI:
  case QU_GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    xc = a + R2*cos(t2);
    yc = 1.0 + R1*cos(t1);
    inside = (fabs(x)<=a && fabs(y)<=1.0); /* in Rect */
    inC1 = ((fabs(y)-yc)*(fabs(y)-yc)+x*x < R1*R1); /* in C1 */
    inC2 = ((fabs(x)-xc)*(fabs(x)-xc)+y*y < R2*R2); /* in C2 */
    if (t1<0.0)
      inside = inside || (fabs(y)>1.0 && inC1);
    else
      inside = inside && !inC1;
    if (t2<0.0)
      inside = inside || (fabs(x)>a && inC2);
    else
      inside = inside && !inC2;
    break;
    
  case RAD_FUNC:
  case QU_RAD_FUNC:
    t1 = atan2(y, x);
    if (t1<0.0)
      t1 += 2*PI;
    rad_func(l, t1, &R1, &R2, &R3);
    inside = (x*x + y*y <= R1*R1);
    break;

  case LINE:
    // return true if below the line (if the line has increasing x values)
    a = l->param[2] - l->param[0];  // vector along line
    b = l->param[3] - l->param[1]; 
    xc = x - l->param[0];  // vector from r0 to r
    yc = y - l->param[1];
    inside = (xc*b - yc*a) >= 0.0;
    break;
    
  case ELL:
    inside = (x<1.0 && y<1.0 && x>0.0 && y>0.0);
    if (x>0.5 && y>0.5)
      inside = 0;
    break;

  case HALF_MUSH:
    a = l->param[0];
    b = l->param[1];
    if (y>0)
      inside = (SQ(x+b/2) + y*y < SQ(1+b/2));
    else
      inside = (y>-a) && (x<0) && (x>-b/2);
    break;
 
  default:
    fprintf(stderr,"Bad billiard type in inside_billiard!\n");
    break;
  }

  // handle desymmetrized shapes by testing if in 1st quadrant...
  switch(l->type) {
  case QU_GEN_RECT_SINAI:
  case DESYM_SINAI:
  case DESYM_SINAI_ALL:
  case QU_STADIUM:
  case QU_STADIUM_ALL:
  case QU_RAD_FUNC:
    if (x<0.0 || y<0.0)
      inside = 0;
  }

  return inside;
}



void desc_by_perim(double f, double *cx, double *cy, double *normal,\
		   double *curv, Billiard *l)
/* fills *cx,*cy with location, & *normal with normal angle, of point a given
 * fraction f around the perim of each billiard model type. v1.2
 *	Added computation of r_n array for Vergini method, v1.7 99/8/25
 *	Removed computation of r_n array for Vergini method, v1.7 99/10/14
 *      (normal angle was between -PI and +PI... no longer true).
 * l->Perim is available for use.
 *      9/10/03 converted to Billliard object.
 * 10/30/03 added filling of curvature too.
 */
{
  double d, dh, dc, fs, s1,s2, alpha, R, n_ang;
  double R1,R2,R3,R4, t1,t2,t3,t4, a, b, x,y,x2,x4,y1,y3,s3,s4;

  /* choose the billiard type and calc coords... */
  switch(l->type) {
	
  case STADIUM:
    // stadium: start at bot-R corner, straight side len = alpha*radius
    alpha = l->param[0];
    d = f*l->Perim;
    if (d<PI) {
      *cx = sin(d) + alpha/2;
      *cy = -cos(d);
      n_ang = PI/2.0 + d;
    } else if (f<=0.5) {
      *cx = -(d - PI) + alpha/2;
      *cy = 1.0;
      n_ang = -PI/2.0;
    } else if ((dh=d-(l->Perim/2.0)) < PI) {
      *cx = -alpha/2 - sin(dh);
      *cy = cos(dh);
      n_ang = -PI/2.0 + dh;
    } else {
      *cx = -(l->Perim - d) + alpha/2;
      *cy = -1.0;
      n_ang = PI/2.0;
    }
    n_ang += PI; // flips inwards normals to outwards ones.
    break;
    
  case QU_STADIUM:
    /* quarter stadium: center at origin, size (1+alpha/2) by 1... */
    alpha = l->param[0];
    d = (1-f)*l->Perim;
    if (d < 0.5*PI) {
      *cx = alpha/2 + cos(d);
      *cy = sin(d);
      n_ang = d;
    } else {
      *cx = l->Perim - d;
      *cy = 1.0;
      n_ang = PI/2;
    }
    break;
    
  case QU_STADIUM_ALL:
    /* quarter stadium: center at origin, size (1+alpha/2) by 1 ... */
    alpha = l->param[0];
    d = (alpha + PI)/2 - f*l->Perim;
    if (d < -(1.0 + alpha/2)) {
      *cx = 0.0;
      *cy = -(1.0 + alpha/2 + d);
      n_ang = PI;
    } else if (d < 0.0) {
      *cx = 1.0 + alpha/2 + d;
      *cy = 0.0;
      n_ang = -PI/2;
    } else if (d < 0.5*PI) {
      *cx = alpha/2 + cos(d);
      *cy = sin(d);
      n_ang = d;
    } else {
      *cx = f*l->Perim;
      *cy = 1.0;
      n_ang = PI/2;
    }
    break;
    
  case DESYM_SINAI:
    /* pi/4 Sinai: center at (1,0)
     * n_th is not set up, since only used by j1 basis... */
    R = l->param[0];
    d = (1-f)*l->Perim;
    if (d < R*PI/4 && R!=0.0) {
      *cx = 1.0 - R*cos(d/R);
      *cy = R*sin(d/R);
      n_ang = -d/R;
    } else {
      *cx = (l->Perim - d)/sqrt(2.0);
      *cy = 1.0 - (*cx);
      n_ang = PI/4;
    }
    break;
    
  case DESYM_SINAI_ALL:
    /* pi/4 Sinai, all perim: center at (1,0)
     * n_th is not set up, since only used by j1 basis... */
    R = l->param[0];
    d = (1-f)*l->Perim;
    if (d < 1.0) {
      *cx = 0.0;
      *cy = 1.0 - d;
      n_ang = PI;
    } else if (d < 2.0-R) {
      *cx = d - 1.0;
      *cy = 0.0;
      n_ang = -PI/2;
    } else if (R!=0.0 && d < 2.0+R*(PI/4 - 1.0)) {
      dc = d - (2.0-R);       // length along arc from bottom.
      *cx = 1.0 - R*cos(dc/R);
      *cy = R*sin(dc/R);
      n_ang = -dc/R;
    } else {
      *cx = (l->Perim - d)/sqrt(2.0);
      *cy = 1.0 - (*cx);
      n_ang = PI/4;
    }
    break;
    
  case GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    x = a + R2*cos(t2);
    y = 1.0 + R1*cos(t1);
    fs = (f>0.5) ? f-0.5 : f; // flip into 3pi/8 to 7pi/8 angle range.
    d = fs*l->Perim;
    s2 = 2*R2*t2;
    s1 = 2*R1*t1;
    if (d<s2) {
      // on left wall
      n_ang = PI + t2*(1.0 - 2*d/s2);
      *cx = -x - R2*cos(n_ang);
      *cy = -R2*sin(n_ang);
    } else {
      // on bottom wall
      n_ang = -PI/2 + t1*(1.0 - 2*(d-s2)/s1);
      *cx = -R1*cos(n_ang);
      *cy = -y - R1*sin(n_ang);
    }
    if (f>0.5) {
      // Inversion under C2v if on other 2 walls...
      n_ang = PI + n_ang;
      *cx = -*cx;
      *cy = -*cy;
    }
    break;
    
  case GEN_ASYM_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    t3 = l->param[3];
    t4 = l->param[4];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    R3 = a/sin(t3);
    R4 = 1/sin(t4);
    x2 = -a - R2*cos(t2);
    y1 = -1.0 - R1*cos(t1);
    x4 = a + R4*cos(t4);
    y3 = 1.0 + R3*cos(t3);
    d = f*l->Perim;
    s2 = 2*R2*t2;
    s1 = 2*R1*t1;
    s3 = 2*R3*t3;
    s4 = 2*R4*t4;
    if (d<s2) {
      // on left wall (2)
      n_ang = PI + t2*(1.0 - 2*d/s2);
      *cx = x2 - R2*cos(n_ang);
      *cy = -R2*sin(n_ang);
    } else if (d < s2+s1) {
      // on bottom wall (1)
      n_ang = -PI/2 + t1*(1.0 - 2*(d-s2)/s1);
      *cx = -R1*cos(n_ang);
      *cy = y1 - R1*sin(n_ang);
    } else if (d < s2+s1+s4) {
      // on right wall (4)
      n_ang = t4*(1.0 - 2*(d-(s1+s2))/s4);
      *cx = x4 - R4*cos(n_ang);
      *cy = -R4*sin(n_ang);
    } else {
      // on top wall (3)
      n_ang = +PI/2 + t3*(1.0 - 2*(d-(s1+s2+s4))/s3);
      *cx = -R3*cos(n_ang);
      *cy = y3 - R3*sin(n_ang);
    }
    break;
    
  case QU_GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    x = a + R2*cos(t2);
    y = 1.0 + R1*cos(t1);
    d = f*l->Perim;
    s2 = R2*t2;
    s1 = R1*t1;
    if (d<s1) {
      // on up wall
      n_ang = PI/2 + t1*d/s1;
      *cx = -R1*cos(n_ang);
      *cy = y - R1*sin(n_ang);
      *curv = -1/R1;
    } else {
      // on right wall
      n_ang = t2*(d-l->Perim)/s2;
      *cx = x - R2*cos(n_ang);
      *cy = -R2*sin(n_ang);
      *curv = -1/R2;
    }
    break;

  case RAD_FUNC:
  case QU_RAD_FUNC:
    t1 = inverse_interp(l->nangles, l->angle, l->arclen, f*l->Perim);
    rad_func(l, t1, &R1, &R2, &R3);
    *cx = R1 * cos(t1);
    *cy = R1 * sin(t1);
    n_ang = t1 - atan(R2/R1);
    *curv = pow((2*R2*R2 + R1*(R1-R3)) / (R1*R1 + R2*R2), 1.5);
    break;

  case LINE:
    *cx = f*l->param[2] + (1.0-f)*l->param[0];
    *cy = f*l->param[3] + (1.0-f)*l->param[1];
    n_ang = atan2(l->param[2]-l->param[0],l->param[1]-l->param[3]);
    *curv = 0.0;
    break;

  case ELL:
    d = 2*f  - floor(2*f);
    a = min(d, 1.0 - d)/2;
    *cx = f + a;
    *cy = 1.0 - f + a;
    n_ang = (d<0.5) * PI/2;
    *curv = 0.0;
    break;

  case HALF_MUSH:
    a = l->param[0];
    b = l->param[1];
    d = l->Perim * f;
    t1 = (d - (1+a+b))/(1+b/2);
    if (d < b/2) {
      *cx = -d;
      *cy = -a;
      n_ang = -PI/2;
      *curv = 0.0;
    } else if (d < 1+a+b) {
      *cx = -b/2;
      *cy = -a + (d-b/2);
      n_ang = PI;
      *curv = 0.0;
    } else {
      *cx = -b/2 + (1+b/2)*sin(t1);
      *cy = (1+b/2)*cos(t1);
      n_ang = PI/2 - t1;
      *curv = 1.0/(1+b/2);
    }
    break;

  default:
    fprintf(stderr,"Bad billiard type in desc_by_perim!\n");
    /* dummy values... */
    *cx = 0.0;
    *cy = 0.0;
    *curv = 0.0;
    n_ang = 0.0;
    break;
  }
  
  // achieve translation from billiard's default location...
  *cx += l->trans_x;
  *cy += l->trans_y;
  
  /* return normal angle... */
  *normal = n_ang;
}




// ============================== CORNERS ===================================

int give_corners(double *cf, double *cal, double *cah, Billiard *l)
  /*
   * Fills arrays with info about each corner:
   *   perim frac, start ("low") angle, finish ("high") angle (follows sense)
   * Returns the number of corners.
   * Note that the corner arrays are 0-indexed but index 0 is dummy.
   * The first corner has index 1.
   *
   * Barnett 8/20/01
   * 8/24/03: included C4-symm shapes (with partial corners at start and end
   *          in order to match to lines of symmetry).
   * 9/10/03 adapted for Billiard object.
   */
{
  int n;
  double R, alpha, a,b,t1,t2,t3,t4,R1,R2,R3,R4;
  
  switch(l->type) {
  case STADIUM:
  case QU_STADIUM:
  case RAD_FUNC:
  case QU_RAD_FUNC:
  case LINE:
    n = 0;
    break;

  case QU_STADIUM_ALL:
    n = 3;
    alpha = l->param[0];
    cf[1] = 0.0;
    cal[1] = PI;
    cah[1] = PI/2;
    cf[2] = (alpha + PI)/(2*l->Perim);
    cal[2] = 0.0;
    cah[2] = -PI/2;
    cf[3] = 1.0 - 1.0/l->Perim;
    cal[3] = -PI/2;
    cah[3] = -PI;
    break;
      
    case DESYM_SINAI:
      R = l->param[0];
      n = 2;
      cf[1] = 0.0;  /* not forgetting to connect with y-axis ! */
      cal[1] = PI/2;
      cah[1] = PI/4;
      cf[2] = (sqrt(2.0)-R)/l->Perim;
      cal[2] = PI/4;
      cah[2] = -PI/4;
      break;

    case DESYM_SINAI_ALL:
      R = l->param[0];
      n = 4;
      cf[1] = 0.0;
      cal[1] = PI/4;
      cah[1] = PI;
      cf[2] = (sqrt(2.0)-R)/l->Perim;
      cal[2] = -PI/4;
      cah[2] = PI/4;
      cf[3] = (sqrt(2.0)-R*(1.0 - PI/4))/l->Perim;
      cal[3] = -PI/2;
      cah[3] = 0.0;
      cf[4] = (l->Perim-1.0)/l->Perim;
      cal[4] = -PI;
      cah[4] = -PI/2;
      break;

  case GEN_ASYM_SINAI: // since this goes CCW, needed to flip cah, cal.
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    t3 = l->param[3];
    t4 = l->param[4];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    R3 = a/sin(t3);
    R4 = 1/sin(t4);
    n = 4;
    cf[1] = 0.0;
    cal[1] = PI + t2;
    cah[1] = PI/2 - t3;
    cf[2] = R2*t2/(R1*t1 + R2*t2 + R3*t3 + R4*t4);
    cal[2] = 3*PI/2 + t1;
    cah[2] = PI - t2;
    cf[3] = (R2*t2 + R1*t1)/(R1*t1 + R2*t2 + R3*t3 + R4*t4);
    cal[3] = t4;
    cah[3] = -PI/2 - t1;
    cf[4] = 1.0 - R3*t3/(R1*t1 + R2*t2 + R3*t3 + R4*t4);
    cal[4] = PI/2 + t3;
    cah[4] = -t4;
    break;

    case QU_GEN_RECT_SINAI:
      a = l->param[0];
      t1 = l->param[1];
      t2 = l->param[2];
      R1 = a/sin(t1);
      R2 = 1/sin(t2);
      n = 1;
      cf[1] = R1*t1/(R1*t1 + R2*t2);
      cal[1] = PI/2 + t1;
      cah[1] = -t2;
      break;

  case GEN_RECT_SINAI: // since this goes CCW, needed to flip cah, cal.
      a = l->param[0];
      t1 = l->param[1];
      t2 = l->param[2];
      R1 = a/sin(t1);
      R2 = 1/sin(t2);
      n = 4;
      cf[1] = 0.0;
      cal[1] = PI + t2;
      cah[1] = PI/2 - t1;
      cf[2] = 0.5 * R2*t2/(R1*t1 + R2*t2);
      cal[2] = 3*PI/2 + t1;
      cah[2] = PI - t2;
      cf[3] = 0.5;
      cal[3] = t2;
      cah[3] = -PI/2 - t1;
      cf[4] = 0.5 + cf[2];
      cal[4] = PI/2 + t1;
      cah[4] = -t2;
      break;

  case ELL:
    n = 3;
    cf[1] = 0.25;
    cal[1] = PI/2;
    cah[1] = 0.0;
    cf[2] = 0.5;
    cal[2] = 0.0;
    cah[2] = PI/2;
    cf[3] = 0.75;
    cal[3] = PI/2;
    cah[3] = 0.0;
    break;

  case HALF_MUSH:
    a = l->param[0];
    b = l->param[1];
    n = 2;
    cf[1] = b/2 / l->Perim;
    cal[1] = -PI/2;
    cah[1] = -PI;
    cf[2] = (1+a+b) / l->Perim;
    cal[2] = PI;
    cah[2] = PI/2;
   break;

    default:
      fprintf(stderr,"Bad billiard type in give_corners!\n");
      n = 0;
      break;
  }

  return n;
}




void desc_outer_by_perim(double f, double d, double corner_pack, \
			 double *cx, double *cy, \
			 double *normal, Billiard *l)
  /*
   * Gives location of 'outer bdry point', and its normal direction, given f,
   * the fraction of the way around the 'outer perimeter'. This 'outer
   * perimeter' is a closed curve maintaining a distance d outside the
   * billiard boundary. Includes corner piece arcs (nonreentrant corners),
   * and handling of reentrant corners. Note f increased in clockwise
   * fashion (cah<cal for nonreentrant corner).
   * Needs d, the distance to outer bdry in billiard units,
   * and corner-pack, the ratio of effective to true length in corner arcs.
   *
   * This routine is used for building oyo* basis sets.
   *
   * Currently doesn't handle reentrant corners.
   * 8/24/03: fixed bug which failed if first corner not at s=0.
   * 6/13/04: reentrant corners (cal>cah) added.
   */
{
  int i;
  double p, p_adj, Perim_c = 0.0, cel[MAX_CORNERS+1], x, y, seg_len, dummy;

  l->cf[0] = 0.0;  // used in loop so starts correctly if cf[1]>0
  cel[0] = 0.0;    // again, so zeroth corner is ok
  for (i=1; i<=l->nc; ++i) {
    if (l->cah[i] < l->cal[i]) // nonreentrant (adds length)
      cel[i] = corner_pack * d * (l->cal[i] - l->cah[i]);
    else                       // reentrant (removes length)
      cel[i] = -2*d * tan((l->cah[i] - l->cal[i])/2);
    Perim_c += cel[i];   // effective total corner perim
    //printf("i=%d, cah=%g, cal=%g\n", i, l->cah[i], l->cal[i]);
  }

  p = f * (l->Perim + Perim_c); // desired eff dist around perim.
  //printf("Perim_c=%g, l->Perim=%g, p=%g\n", Perim_c, l->Perim, p);

  // move round perim from s=0 until locate current segment or corner...
  for(i=0; i<=l->nc; ++i) {
    // find current segment length...
    if (i==l->nc) {
      seg_len = (1.0 - l->cf[i])*l->Perim;
      if (cel[i]<0) // remove any length due to reentrant corners
	seg_len += cel[i]/2;
    } else {
      seg_len = (l->cf[i+1] - l->cf[i])*l->Perim;
      if (cel[i]<0) // remove any length due to reentrant corners
	seg_len += cel[i]/2;
      if (cel[i+1]<0) // at other end too
	seg_len += cel[i+1]/2;
    }
    if (p < seg_len) {
      // in segment following corner i...
      // printf("        in segment following corner %d by dist %g\n",i,p);
      p_adj = p; // dist along current seg
      if (cel[i]<0)
	p_adj -= cel[i]/2; // adjust length for this reentrant corner
      desc_by_perim(l->cf[i] + p_adj/l->Perim, &x, &y, normal, &dummy, l);
      *cx = x + cos(*normal)*d;
      *cy = y + sin(*normal)*d;
      break;
    }
    p -= seg_len;
    if (p < cel[i+1]) {
      // in corner i+1...
      // printf("        in corner %d by frac %g\n",i+1,p/cel[i+1]);
      // *normal is undefined here...
      desc_by_perim(l->cf[i+1], &x, &y, normal, &dummy, l);
      *normal = l->cal[i+1] + (l->cah[i+1] - l->cal[i+1])*p/cel[i+1];
      *cx = x + cos(*normal)*d;
      *cy = y + sin(*normal)*d;
      break;
    }
    p -= max(0.0, cel[i+1]);
  }
}



void bounding_box(double *xl,double *xh,double *yl,double *yh, \
	Billiard *l)
  /*
   * Returns bounding box of given billiard type
   * Separated as a routine, 99/10/14.
   * Excludes RAD_FUNC types.
   */
{
  double alpha,R,a,b,R1,R2,R3,R4,t1,t2,t3,t4,xc,yc,y1,x2,y3,x4;

  switch(l->type) {
  case STADIUM:
    alpha = l->param[0];
    *xh = 1.0+alpha/2;
    *xl = -(*xh);
    *yh = 1.0;
    *yl = -(*yh);
    break;
    
  case QU_STADIUM:
  case QU_STADIUM_ALL:
    alpha = l->param[0];
    *xh = 1.0+alpha/2;
    *xl = 0.0;
    *yh = 1.0;
    *yl = 0.0;
    break;
    
  case DESYM_SINAI:
  case DESYM_SINAI_ALL:
    R = l->param[0];
    *xh = 1.0 - R/sqrt(2.0);
    *xl = 0.0;
    *yh = 1.0;
    *yl = 0.0;
    break;
    
  case GEN_RECT_SINAI:
  case QU_GEN_RECT_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    xc = a + R2*cos(t2);
    yc = 1.0 + R1*cos(t1);
    *xh = max(a, xc-R2);
    if (t1<-PI/2)
      *xh = max(*xh, -R1);
    *xl = (l->type==QU_GEN_RECT_SINAI) ? 0.0 : -*xh;
    *yh = max(1.0, yc-R1);
    if (t2<-PI/2)
      *yh = max(*yh, -R2);
    *yl = (l->type==QU_GEN_RECT_SINAI) ? 0.0 : -*yh;
    break;
    
  case GEN_ASYM_SINAI:
    a = l->param[0];
    t1 = l->param[1];
    t2 = l->param[2];
    t3 = l->param[3];
    t4 = l->param[4];
    R1 = a/sin(t1);
    R2 = 1/sin(t2);
    R3 = a/sin(t3);
    R4 = 1/sin(t4);
    x2 = -a - R2*cos(t2);
    y1 = -1.0 - R1*cos(t1);
    x4 = a + R4*cos(t4);
    y3 = 1.0 + R3*cos(t3);
    *xl = min(-a, x2+R2); /* bottom */
    if (t1<-PI/2)
      *xl = min(*xl, -R1);
    if (t3<-PI/2)
      *xl = min(*xl, -R3);
    *xh = max(a, x4-R4); /* top */
    if (t1<-PI/2)
      *xl = max(*xl, R1);
    if (t3<-PI/2)
      *xl = max(*xl, R3);
    *yl = min(-1.0, y1+R1); /* left */
    if (t2<-PI/2)
      *yl = min(*yl, -R2);
    if (t4<-PI/2)
      *yl = min(*yl, -R4);
    *yh = max(1.0, y3-R3); /* right */
    if (t2<-PI/2)
      *yh = max(*yh, R2);
    if (t4<-PI/2)
      *yh = max(*yh, R4);
    break;
    
  case LINE:
    *xl = min(l->param[0], l->param[2]);
    *xh = max(l->param[0], l->param[2]);
    *yl = min(l->param[1], l->param[3]);
    *yh = max(l->param[1], l->param[3]);
    break;

  case ELL:
    *xl = 0.0;
    *yl = 0.0;
    *xh = 1.0;
    *yh = 1.0;
    break;

  case HALF_MUSH:
    a = l->param[0];
    b = l->param[1];
    *xl = -b/2;
    *xh = 1;
    *yl = -a;
    *yh = 1+b/2;
    break;

  default:
    fprintf(stderr,"Bad billiard type in bounding_box!\n");
    break;
  }

  /* achieve translation... */
  *xl += l->trans_x;
  *xh += l->trans_x;
  *yl += l->trans_y;
  *yh += l->trans_y;

}


/* -------------------------- BILLIARD IMPACT ROUTINES -------------------
   (for Vergini grazing-type basis sets)
*/

double reflect_into_second_quadrant(double dir)
  /*
    Flips any direction into the second quandrant, using C2v symmetry.
    dir does NOT have to be [0,2pi] or in any particular range.

    Utility for max_billiard_impact
    */
{
  // put in range [-pi,pi]...
  double mod = fmod(dir, 2*PI);  // changed from drem for Sun math compat

  // flip x-axis
  if (mod<0.0) mod = -mod;
  // flip y-axis
  if (mod>PI/2) mod = PI-mod;

  return mod;
}

double max_billiard_impact(double dir, Billiard *l)
  /*
    Returns the max impact parameter possible which glances the billiard on its
    right-hand side, ie glances in an anti-clockwise sense.
    dir is measured anticlockwise from the +ve x-axis.
    Desymmetrized billiards are treated as the full (unfolded) shape.
    Barnett 99/11/1

    9/10/03 Doesn't work for nonzero trans_{x,y}. Generally obsolete.
    */
{
  double alpha, R, dir2, x, y, a, theta_corner;

  switch(l->type) {
  case STADIUM:
  case QU_STADIUM:
  case QU_STADIUM_ALL:
    alpha = l->param[0];
    dir2 = reflect_into_second_quadrant(dir);
    return 1 + sin(dir2)*alpha/2;
    
  case DESYM_SINAI:
  case DESYM_SINAI_ALL:
    R = l->param[0];
    dir2 = reflect_into_second_quadrant(dir);
    // choose point where ray glances billiard...
    if (dir2 < 3*PI/4) {
      x = 1 - R/sqrt(2.0);
      y = R/sqrt(2.0);
    } else {
      x = 0.0;
      y = 1.0;
    }
    return x*sin(dir2) - y*cos(dir2);
    
  case GEN_RECT_SINAI:
  case GEN_ASYM_SINAI:
  case QU_GEN_RECT_SINAI:
    // point (a,1.0) always determines max impact, if concave.
    a = l->param[0];
    dir2 = reflect_into_second_quadrant(dir);
    theta_corner = atan(1.0/a);
    R = sqrt(1.0 + a*a);
    return R*cos(dir2-theta_corner);

  default:
    fprintf(stderr,"Bad billiard type in max_billiard_impact!\n");
    return 0.0;
  }
}

/* --------------------------------- RADIAL FUNC ----------------------- */

void rad_func(Billiard *l, double t, double *r, double *rt, double *rtt)
  /* Insert radial function types and their 2 derivs here. t = theta.
   *
   * barnett 12/13/03
   */
{
  double e1 = l->param[1], e2 = l->param[2]; // eccentricities
  double beta, a, s, t0, e3, d, c, ct, st;
  int w;
  
  switch ((int)(l->param[0])) {

  case 0: // quarter symmetry (first quadrant)
    *r = 1 + e1*cos(2*t) + e2*cos(4*t);
    *rt = -e1*2*sin(2*t) - e2*4*sin(4*t);
    *rtt = -e1*4*cos(2*t) - e2*16*cos(4*t);
    break;
    
  case 1:  // nonsymmetric radial funcs
    beta = 5*PI/4;  // range of growing angle
    w = 5;          // controls wiggliness, must be odd
    e3 = 0.02;      // gave nice foci for rf param1 = .5
    if (t<beta) {
      s = w*PI/beta;
      t0 = 0.0;
    } else {
      s = PI/(2*PI - beta);
      t0 = 2*PI;
    }
    a = s*(t - t0);
    *r = 1 - e2*cos(a) + e1*cos(2*t) + e3*cos(4*t);
    *rt = e2*s*sin(a) - e1*2*sin(2*t) - e3*4*sin(4*t);
    *rtt = e2*s*s*cos(a) - e1*4*cos(2*t) - e3*16*cos(4*t);
    break;
    
  case 2: // nonsymmetric perturbed ellipse, centered on 1 focus
    beta = 5*PI/4;  // range of growing angle
    w = 1;          // controls wiggliness, must be odd
    if (t<beta) {
      s = w*PI/beta;
      t0 = 0.0;
    } else {
      s = PI/(2*PI - beta);
      t0 = 2*PI;
    }
    a = s*(t - t0);
    d = 1/(1-e1*cos(t));       // do proper ellipse
    c = 1 - e1*e1;        // scale major axis to 2
    *r = c*d - e2*cos(a);    // add to ellipse
    *rt = e2*s*sin(a) - e1*c*sin(t)*d*d;
    *rtt = e2*s*s*cos(a) - e1*c*cos(t)*d*d - 2*e1*e1*c*sin(t)*sin(t)*d*d*d;
    break;
    
  case 3: // nonsymmetric perturbed ellipse, centered in middle
    beta = 5*PI/4;  // range of growing angle
    w = 1;          // controls wiggliness, must be odd
    if (t<beta) {
      s = w*PI/beta;
      t0 = 0.0;
    } else {
      s = PI/(2*PI - beta);
      t0 = 2*PI;
    }
    a = s*(t - t0);
    ct = cos(t);
    st = sin(t);
    d = 1/sqrt(1-e1*e1*ct*ct);       // centered ellipse
    c = sqrt(1 - e1*e1);        // scale x width to 2
    *r = c*d - e2*cos(a);    // add to ellipse
    *rt = e2*s*sin(a) - c*e1*e1*st*ct*d*d*d;
    *rtt = e2*s*s*cos(a) - c*e1*e1*d*d*d*(ct*ct - st*st + \
					  3*st*st*ct*ct*e1*e1*d*d);
    break;
    
  case 4: // Hakan's quadrupolar cavity
    a = 1.0/sqrt(1 + e1*e1/2);
    *r = a*(1 + e1*cos(2*t));
    *rt = -2*a*e1*sin(2*t);
    *rtt = -4*a*e1*cos(2*t);
    break;
    
 default:
    fprintf(stderr, "bad type %g in rad_func!\n", l->param[0]);
  }
  
}



/* --------------------------------- BUILD BILLIARD ----------------------- */

#define B_RAD_FUNC  100 // resolution per wavelength for proto-collocation pts

int build_billiard(Billiard *l, double k_base)
  /*
   * Set up Billiard object, return error status
   *
   * barnett 9/10/03
   * 12/13/03 radial func added, needs access to k_base
   * 5/9/05 trans for rad_funcs debugged
   */
{
  int btp, j, M;
  double t, dt, r, rt, ro, rto, ar, art, rtt, x, y, t_max;

  /* read trans_{x,y} info from param... */
  btp = billiard_trans_param[(int)(l->type)];
  if (btp+1 >= billiard_num_params[(int)(l->type)]) {
    fprintf(stderr, "billiard_trans_param overruns billiard param array!\n");
    return 1;
  } else if (btp == -1) {
    l->trans_x = 0.0;
    l->trans_y = 0.0;
  } else {
    l->trans_x = l->param[btp];
    l->trans_y = l->param[btp+1];
  }

  if (l->type==RAD_FUNC || l->type==QU_RAD_FUNC) {
    // rad funcs need their own set of proto-bdry pts, meas perim etc.

    t_max = 2*PI / ((l->type==QU_RAD_FUNC) ? 4 : 1);  // quarter or not?
    dt = 2*PI/(k_base*B_RAD_FUNC); // for radius approx 1
    M = round_and_clip(t_max/dt, "nangles in build_billiard");
    l->nangles = M;
    l->angle = dvector(0,M);
    l->arclen = dvector(0,M);
    l->angle[0] = 0.0;
    l->arclen[0] = 0.0;
    rad_func(l, 0.0, &ro, &rto, &rtt); // initialize old values of r, rt
    l->Area = 0.0;
    l->xl = 0.0;
    l->yl = 0.0;
    l->xh = ro;
    l->yh = 0.0;
    for (j=1; j<=M; ++j) {
      l->angle[j] = dt * j;
      rad_func(l, dt * j, &r, &rt, &rtt);
      ar = (r + ro)/2; // average value of r for this segment (j-1 to j)
      art = (rt + rto)/2;
      l->arclen[j] = l->arclen[j-1] + dt * sqrt(ar*ar + art*art);
      l->Area += dt * ar*ar/2;
      x = r * cos(l->angle[j]);
      y = r * sin(l->angle[j]);
      if (x > l->xh)
	l->xh = x;
      if (y > l->yh)
	l->yh = y;
      if (x < l->xl) // only needed for full circle case
	l->xl = x;
      if (y < l->yl) // ditto
	l->yl = y;
      ro = r; // keep as old (previous) values
      rto = rt;
    }
    l->Perim = l->arclen[M];
    l->xl += l->trans_x;
    l->xh += l->trans_x;
    l->yl += l->trans_y;
    l->yh += l->trans_y;


  } else {
    l->Perim = perimeter(l);
    l->Area = area(l);
    bounding_box(&(l->xl), &(l->xh), &(l->yl), &(l->yh), l);
  }

  l->nc = give_corners(l->cf, l->cal, l->cah, l);
  if (verb&0x02)
    printf("\tbox: x=[%g,%g] y=[%g,%g]\n", l->xl, l->xh, l->yl, l->yh);
  return 0;
}


/* --------------------- dump bounding box --------------------------------- */
void dump_billiard(Billiard *l, char *head)
{
  FILE *fp;
  char mname[LEN];

  sprintf(mname,"%s.bound_box",head);
  if ((fp = fopen(mname,"w")) == NULL)
    printf("file: %s, write error...\n",mname);
  
  fprintf(fp,"# %s from vergini\n\n# gnuplot lines: x y\n\n", mname);
        
  fprintf(fp,"%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n", \
	  l->xl,l->yl,l->xh,l->yl,l->xh,l->yh,l->xl,l->yh,l->xl,l->yl);
  
  fprintf(fp,"\n");
  fclose(fp);
}




/* -------------------- READ BILLIARD COMMAND LINE ARGS -------------------- */

int parse_billiard(char *s, Billiard *l)
  /*
   * Read string from command line argument into Billiard object type and
   * param array. Returns -1 if failure.
   *
   * Note that b, the collocation point density, is not read here any more.
   *
   * barnett 9/10/03
   * 11/16/03 changed to colon-separated
   */
{
  int n, i;
  char tag[LEN];

  /* Identify what type of billiard we have, call it t... */
  billiard_t t = NBILLIARDS;
  sscanf(s, "%[^:]%n", tag, &i);
  if (verb&0x80)
    printf("tag=%s, i=%d\n", tag, i);
  s += i;
  if (s[0]==':')
    ++s;
  // s now points to 2nd colon-separated item, or \0
  for(i=0; i<NBILLIARDS; ++i)
    if (strcmp(tag, billiard_tags[i]) == 0)
      t = (billiard_t)i;
  if (t==NBILLIARDS) {
    fprintf(stderr, "Invalid billiard type: %s!\n", tag);
    return -1;
  }
  l->type = t;
  if (verb)
    printf("\tBilliard type: %s\n", billiard_desc[(int)t]);

   /* Check for param array overrun... */
  n = billiard_num_params[(int)t];
  if (n > MAX_N_BILLIARD_PARAMS) {
    fprintf(stderr, "billiard_num_params[%d]=%d exceeds max!\n", (int)t, n);
    return -1;
  }
  /* Read in params (all doubles) to param array... */
  if (separated_numbers(s, l->param, n) != n)
    return -1;

  /* feedback... */
  if (verb) {
    printf("\t\t%s\n\t\t", billiard_param_desc[(int)t]);
    for (i=0; i<n; ++i)
      printf("%g ", l->param[i]);
    printf("\n\n");
  }

  return n;
}


void show_billiard_properties(Billiard *l, double k)
  /* includes k-dependent, quantum number stuff
   */
{
  printf("System properties:\n");
  printf("	Physical: area = %.6g, perim = %.6g\n", l->Area, l->Perim);
  printf("	k = %.6g, E = k^2 = %.6g\n", k, k*k);
  printf("	1st-order Weyl level spacing dk = %.6g, dE = %.6g\n",\
	 (2*PI)/(k*l->Area), (4*PI)/l->Area);
  printf("      approx quantum number (area term) %d\n",\
	 (int)(0.5 + k*k*l->Area/(4*PI)));
  printf("      approx quantum number (area + perim terms) %d\n",\
	 (int)(0.5 + k*k*l->Area/(4*PI) - k*l->Perim/(4*PI)));
}  



/* -------------------------- SHOW BILLIARD ARGS USAGE ------------------ */
void show_billiard_usage()
{
  int i;
  
  fprintf(stderr, "\nPossible BILLIARDS, with param input formats.\n");
  for(i=0; i<NBILLIARDS; ++i) {
    fprintf(stderr, "    %s\n", billiard_desc[i]);
    fprintf(stderr, "\t%8s %s\n", billiard_tags[i], billiard_param_desc[i]);
  }
}
