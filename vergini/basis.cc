/*	Basis Set- and Basis Function-specific code for BDRY package.
 *
 *	Alex Barnett 99/10/14
 *	heavily modified for general mix of scaling funcs (EPWs, etc) 99/10/27
 *      show_properties moved across from billiard.cc, 00/4/26
 *      oyooo set and yo_odd_odd type added 8/24/03
 *      9/8/03 Basis_Set object created, converted routines for this.
 *      9/10/03 adapted to Billiard object.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "cxml.h"
#include "nrutil.h"
#include "useful.h"
#include "billiard.h"
#include "basis.h"
#include "verb.h"


// ---------------------------- Evanescent Plane Waves -----------------------

double eval_EPW(double x, double y, double nx, double ny,\
	double cosha, double sinha, double b, int re_or_im)
/*
  Finds val of basic (unsymmetrized) EPW scaling-func at dimless
  coord (x,y).
  PW normal vec is (nx,ny), sinh(EW alpha) is sinha,
  offset along grow direction is b.
  sin_or_cos is flag: =0 for cos (Re part), =1 for sin (Im part).
  
  Used by any EPW_{RE,IM} basis element call.
*/
{
  // b is k_0-scaled impact parameter, here.
  // Don't forget it also goes inside sinha!
  double ampl = exp( -sinha*( -x*ny + y*nx + b ) );
  double phase = cosha*( x*nx +y*ny );

  return re_or_im ? ampl*sin(phase) : ampl*cos(phase);
}


void EPW_grad(double &ddx,double &ddy, double x, double y,\
	      double ct, double st,double ca,double sa,double b, \
	      int re_or_im)
  /*
    Puts grad of Re or Im part of EPW basis func into (ddx,ddy), passed
    in by reference. 
    ct = cos(theta), st = sin(theta), ca = cosh(alpha), sa = sinh(alpha).
    Barnett 99/10/29
    Used by any EPW_{RE,IM} basis_deriv element call.
    */
{
  double ref = eval_EPW(x,y,ct,st,ca,sa,b,0);
  double imf = eval_EPW(x,y,ct,st,ca,sa,b,1);

  if (!re_or_im) {
    // Real part wanted...
    ddx = st*sa*ref - ct*ca*imf;
    ddy = -ct*sa*ref - st*ca*imf;
  } else {
    // Imag part wanted...
    ddx = ct*ca*ref + st*sa*imf;
    ddy = st*ca*ref - ct*sa*imf;
  }

}


void EPW_everything(double &val, double &ddx,double &ddy, \
		    double &dxx,double &dyy, double &dxy, double x,double y,\
		    double ct, double st, double ca,double sa,double b, \
		    int re_or_im)
  /*
    returns val, (ddx,ddy) and dxx,dyy,dxy of Re or Im part of EPW basis func.
    passed in by reference.

    Hessian matrix is ( dxx dxy ), with dxy=dyx (symmetric).
                      ( dyx dyy )

    Barnett 99/11/6
    Used by any EPW_{RE,IM} basis_derivs element call.
    */
{
  double ref = eval_EPW(x,y,ct,st,ca,sa,b,0);
  double imf = eval_EPW(x,y,ct,st,ca,sa,b,1);
  double allfour = ct*st*ca*sa;
  double ca2 = ca*ca;
  double sa2 = sa*sa;
  double ct2 = ct*ct;
  double st2 = st*st;

  if (!re_or_im) {
    // Real part wanted...
    val = ref;
    ddx = st*sa*ref - ct*ca*imf;
    ddy = -ct*sa*ref - st*ca*imf;
    dxx = -( (ct2*ca2 - st2*sa2)*ref + 2*allfour*imf );
    dyy = -( (st2*ca2 - ct2*sa2)*ref - 2*allfour*imf );
    dxy = -( ct*st*(ca2 + sa2)*ref - ca*sa*(ct2 - st2)*imf );
  } else {
    // Imag part wanted...
    val = imf;
    ddx = ct*ca*ref + st*sa*imf;
    ddy = st*ca*ref - ct*sa*imf;
    dxx = -( (ct2*ca2 - st2*sa2)*imf - 2*allfour*ref );
    dyy = -( (st2*ca2 - ct2*sa2)*imf + 2*allfour*ref );
    dxy = -( ct*st*(ca2 + sa2)*imf + ca*sa*(ct2 - st2)*ref );
  }

}


// ----------------------------- Y_0 basis funcs ---------------------------

double eval_yo(double x, double y, double cx, double cy)
{
  double xd = x-cx, yd = y-cy;
  double r = sqrt(xd*xd + yd*yd);    // distance
  // Note prefactor and sign: -(1/4)Y_0.
  return -y0(r)/4;
}

void yo_grad(double &ddx, double &ddy, double x, double y, \
	     double cx, double cy)
{
  double xd = x-cx, yd = y-cy;       // note signs important for grad!
  double r = sqrt(xd*xd + yd*yd);    // distance
  // -(1/4) Y_0 deriv gives (1/4) Y_1...
  double ddr = y1(r)/4;
  ddx = (xd/r) * ddr;
  ddy = (yd/r) * ddr;
}

void yo_everything(double &val, double &ddx, double &ddy, double &dxx, \
		   double &dyy, double &dxy, double x, double y, double cx,\
		   double cy)
  /* Barnett 1/16/04
   */
{
  double xd = x-cx, yd = y-cy;       // note signs important for grad!
  double r = sqrt(xd*xd + yd*yd);    // distance
  double ddr = y1(r)/4;              // use like -Y_1 in formulae
  val = -y0(r)/4;                    // use like Y_0 in formulae
  double ir = 1/r;                   // saves 4 divisions
  double irrbracket = -ir*ir*(val + 2*ir*ddr); // (-Y0+(2/r)Y1)/r^2 in formula
  ddx = xd*ir * ddr;
  ddy = yd*ir * ddr;
  dxx = xd*xd*irrbracket + ir*ddr;
  dyy = yd*yd*irrbracket + ir*ddr;
  dxy = xd*yd*irrbracket;
}



  
// ----------------- General basis funcs & symmetrization --------------------

double eval_basis(double x, double y, int j, Basis_Set *s)
/*
  Generalized to many basis function types allowed within a set, which
  forced passing in of basis element index j.
  Specialized to scaling functions only, to make Vergini method more
  transparent.
  Symmetrization of every basis type except RPWs is handled here.

  Returns phi_j(x,y). (x,y) are dimensionless coords, ie scaled by k.
  Requires Basis_Set to have been built. (Global basis arrays removed).
  
  Barnett 99/10/27
  9/8/03 converted for Basis_Set object.
*/
{
  double phi;
  double ca, sa, a, b, th;
  double ct = s->nx[j], st = s->ny[j];  // plane wave direction unit vector
  int i, re_or_im = 1;

  switch(s->t[j]) { /* Choose functional form based on current basis func */

    // Re and Im parts of a general real plane wave...
  case RPW_RE:
    phi = cos( x*ct + y*st );
    // printf("RPW_RE: %f %f %f %f %f\n",x,y,ct,st,phi);
    break;
  case RPW_IM:
    phi = sin( x*ct + y*st );
    // printf("RPW_IM: %f %f %f %f %f\n",x,y,ct,st,phi);
    break;
	
    // only Re part of a odd-odd symm real plane wave is non-zero...
  case RPW_ODD_ODD:
    phi = sin( x*ct ) * sin( y*st );
    break;

    // only Re part...
  case RPW_EVEN_EVEN:
    phi = cos( x*ct ) * cos( y*st );
    break;

    // Re and Im parts of a odd-odd evan plane wave...
  case EPW_RE_ODD_ODD:
    re_or_im = 0;
  case EPW_IM_ODD_ODD:
    sa = s->a[j];
    ca = s->c[j];
    b = s->b[j];
    phi = eval_EPW(x,y,ct,st,ca,sa,b,re_or_im) \
      - eval_EPW(-x,y,ct,st,ca,sa,b,re_or_im) \
      - eval_EPW(x,-y,ct,st,ca,sa,b,re_or_im) \
      + eval_EPW(-x,-y,ct,st,ca,sa,b,re_or_im);
    break;

  case YO:
    phi = eval_yo(x, y, s->a[j], s->b[j]);
    break;

  case YO_ODD_ODD:
    a = s->a[j];
    b = s->b[j];
    phi = eval_yo(x, y, a, b) - eval_yo(x, y, -a, b) \
      - eval_yo(x, y, a, -b) + eval_yo(x, y, -a, -b);
   break;

  case YO_EVEN_EVEN:
    a = s->a[j];
    b = s->b[j];
    phi = eval_yo(x, y, a, b) + eval_yo(x, y, -a, b) \
      + eval_yo(x, y, a, -b) + eval_yo(x, y, -a, -b);
   break;

  case RPW_REFL:
    phi = sin( x*ct ) * sin( y*st );
    th = PI/s->d[j];
    sa = sin(th);
    ca = cos(th);
    for (i=1; i<s->d[j]; ++i) {
      a = ca*x - sa*y;
      b = sa*x + ca*y;
      x = a; y = b;  // spin (x,y)
      phi += sin( x*ct ) * sin( y*st );
    }
    break;

  case FB:
    th = (atan2(y,x) - s->a[j]); // angle relative to wedge start angle
    if (th<0.0)
      th += 2*PI;
    phi = sin(th*s->b[j]) * gsl_sf_bessel_Jnu(s->b[j], sqrt(x*x + y*y));
    break;

  default:
    fprintf(stderr,"Bad scaling function type in eval_basis!\n");
    phi = 0.0;

  }
  return phi;
}



void eval_basis_deriv(double &ddx, double &ddy, double x, double y, int j, \
		      Basis_Set *s)
  /*
    Calc derivative of scaling basis func wrt dimless scaling coord.
    Returns answer in (ddx,ddy).
    Symmetrization of every basis type except RPWs is handled here.
    Uses analytic formulas which must agree with those in eval_basis.
  */
{
  double ddx_s, ddy_s, ddr, ddth;
  double re,im,ca,sa,a,b, th, r, ir, nu, beslower;
  double ct = s->nx[j], st = s->ny[j];  // plane wave direction unit vector
  double nxa, nxb, nya, nyb;
  int i, re_or_im = 1;

  switch(s->t[j]) {

    // Re and Im parts of a general real plane wave...
  case RPW_RE:
    im = sin( x*ct + y*st );
    ddx = -ct * im;
    ddy = -st * im;
    break;
  case RPW_IM:
    re = cos( x*ct + y*st );
    ddx = ct * re;
    ddy = st * re;
    break;
    
    // only Re part of a odd-odd symm real plane wave is non-zero...
  case RPW_ODD_ODD:
    ddx = ct * cos( x*ct ) * sin( y*st );
    ddy = st * sin( x*ct ) * cos( y*st );
    break;

  case RPW_EVEN_EVEN:
    ddx = -ct * sin( x*ct ) * cos( y*st );
    ddy = -st * cos( x*ct ) * sin( y*st );
    break;

    // Re and Im parts of a odd-odd evan plane wave...
  case EPW_RE_ODD_ODD:
    re_or_im = 0;                 // 0 for Re, 1 for Im.
  case EPW_IM_ODD_ODD:
    sa = s->a[j];
    ca = s->c[j];
    b = s->b[j];
    // add reflections of gradients with correct signs - careful!
    EPW_grad(ddx,ddy, x,y,ct,st,ca,sa,b,re_or_im);
    EPW_grad(ddx_s,ddy_s, -x,y,ct,st,ca,sa,b,re_or_im);
    ddx += ddx_s; ddy -= ddy_s;
    EPW_grad(ddx_s,ddy_s, x,-y,ct,st,ca,sa,b,re_or_im);
    ddx -= ddx_s; ddy += ddy_s;
    EPW_grad(ddx_s,ddy_s, -x,-y,ct,st,ca,sa,b,re_or_im);
    ddx -= ddx_s; ddy -= ddy_s;
    break;

  case YO:
    yo_grad(ddx, ddy, x, y, s->a[j], s->b[j]);
    break;

  case YO_ODD_ODD:
  case YO_EVEN_EVEN:
    a = s->a[j]; // func origin
    b = s->b[j];
    yo_grad(ddx, ddy, x, y, a, b);
    yo_grad(ddx_s, ddy_s, x, y, -a, b); // image of func origin.
    if (s->t[j]==YO_ODD_ODD) {
      ddx -= ddx_s; ddy -= ddy_s;         // sign = reflection parities (o-o)
    } else {
      ddx += ddx_s; ddy += ddy_s;         // (e-e)
    }
    yo_grad(ddx_s, ddy_s, x, y, a, -b);
    if (s->t[j]==YO_ODD_ODD) {
      ddx -= ddx_s; ddy -= ddy_s;
    } else {
      ddx += ddx_s; ddy += ddy_s;
    }
    yo_grad(ddx_s, ddy_s, x, y, -a, -b);
    ddx += ddx_s; ddy += ddy_s;
    break;

  case RPW_REFL:
    ddx = ct * cos( x*ct ) * sin( y*st );
    ddy = st * sin( x*ct ) * cos( y*st );
    th = PI/s->d[j];
    sa = sin(th);
    ca = cos(th);
    nxa = nyb = 1.0; nxb = nya = 0.0; // set up spinning normals nx, ny
    for (i=1; i<s->d[j]; ++i) {
      a = ca*x - sa*y;
      b = sa*x + ca*y;
      x = a; y = b;
      a = ca*nxa - sa*nxb; // spin normal nx in same direc to (x,y)
      b = sa*nxa + ca*nxb;
      nxa = a; nxb = b;
      a = ca*nya - sa*nyb;
      b = sa*nya + ca*nyb;
      nya = a; nyb = b;
      a = ct * cos( x*ct ) * sin( y*st ); // unrotated grad
      b = st * sin( x*ct ) * cos( y*st );
      ddx += a*nxa + b*nxb;
      ddy += a*nya + b*nyb;
    }
    break;

  case FB:
    th = (atan2(y,x) - s->a[j]); // angle relative to wedge start angle
    if (th<0.0)
      th += 2*PI;
    r = sqrt(x*x + y*y);
    ir = 1.0/r;
    b = s->b[j];     // beta param
    // printf("b[%d]=%g\n", j, b);
    a = gsl_sf_bessel_Jnu(b, r);         // temp storage of bessel value
    nu = b-1.0;       // order lowered by 1
    if (nu<0.0) {        // hand-code the nu<0 case (NumRec p.242)
      nu = -nu;
      beslower = cos(nu*PI)*gsl_sf_bessel_Jnu(nu, r) - sin(nu*PI)*gsl_sf_bessel_Ynu(nu, r);
    } else
      beslower = gsl_sf_bessel_Jnu(b-1.0, r);   // J_{nu-1}(r)
    ddr = sin(th*b) * (beslower - b*ir*a); // radial deriv
    ddth = b*cos(th*b) * a;                // angular deriv
    ddx = ddr*x*ir - ddth*y*ir*ir;       // chain rule (r,th)->(x,y)
    ddy = ddr*y*ir + ddth*x*ir*ir;
    break;

  default:
    fprintf(stderr,"Bad scaling function type in eval_basis_deriv!\n");
    
  }
}



void eval_basis_everything(double &val, double &ddx, double &ddy,\
			   double &dxx, double &dyy, double &dxy,\
			   double x, double y, int j, Basis_Set *s)
/*
  Calc value, and 1st and 2nd derivatives, of scaling basis func.
  Value duplicates calc in eval_basis, and 1st deriv duplicates calc in
  eval_basis_deriv. 2nd derviative code is new. The point of bundling
  all the calcs together is speed.
  Returns value (val), gradient (ddx,ddy) and elements of hessian matrix
  dxx, dyy, dxy(=dyx).
  Derivatives are wrt dimensionless scaling coord.
  Analytic formulas must agree with those in eval_basis.
  Symmetrization of every basis type except RPWs is handled here:
    Note signs for Y0 bases match reflection parities since the func origin
    (a,b) was reflected while (x,y) stays fixed. However, for EPW it is only
    possible to change (x,y), so the signs on derivatives are different.

  Ideas for updating code (1/16/04):
  * perform symmetrization independently of basis func type (less errors).
      This would involve a new symmetrization module.
  * make (val,ddx,ddy,dxx,dyy,dxy) a struct.
  * have conditional statements to select highest deriv wanted (will this
      slow down operation much, compared to time to do a bessel eval?)

      NOTE: the second derivs are only needed by Alex-type S,T matrices.
      If you stick to Quasi-type S,T, or vergini only, not needed.

  1/16/04 added Y0 type basis funcs
*/
{
  double val_s,ddx_s,ddy_s,dxx_s,dyy_s,dxy_s;
  double cx,cy,sx,sy;
  double re,im,ca,sa,a,b;
  double ct = s->nx[j], st = s->ny[j];  // plane wave direction unit vector
  int re_or_im = 1;

  switch(s->t[j]) {

    // Re and Im parts of a general real plane wave...
  case RPW_RE:
    im = sin( x*ct + y*st );
    ddx = -ct * im;
    ddy = -st * im;
    re = cos( x*ct + y*st );
    dxx = -ct*ct * re;
    dyy = -st*st * re;
    dxy = -ct*st * re;
    val = re;
    break;
  case RPW_IM:
    re = cos( x*ct + y*st );
    ddx = ct * re;
    ddy = st * re;
    im = sin( x*ct + y*st );
    dxx = -ct*ct * im;
    dyy = -st*st * im;
    dxy = -ct*st * im;
    val = im;
    break;
	
    // only Re part of a odd-odd symm real plane wave is non-zero...
  case RPW_ODD_ODD:
    cx = cos( x*ct );
    cy = cos( y*st );
    sx = sin( x*ct );
    sy = sin( y*st );
    ddx = ct * cx * sy;
    ddy = st * sx * cy;
    re = sx * sy;
    dxx = -ct*ct * re;
    dyy = -st*st * re;
    dxy = ct*st * cx * cy;
    val = re;
    break;

  case RPW_EVEN_EVEN:
    cx = cos( x*ct );
    cy = cos( y*st );
    sx = sin( x*ct );
    sy = sin( y*st );
    ddx = -ct * sx * cy;
    ddy = -st * cx * sy;
    re = cx * cy;
    dxx = -ct*ct * re;
    dyy = -st*st * re;
    dxy = ct*st * sx * sy;
    val = re;
    break;

    // Re and Im parts of a odd-odd evan plane wave...
  case EPW_RE_ODD_ODD:
    re_or_im = 0;                 // 0 for Re, 1 for Im.
  case EPW_IM_ODD_ODD:
    sa = s->a[j];
    ca = s->c[j];
    b = s->b[j];
    EPW_everything(val,ddx,ddy,dxx,dyy,dxy,x,y,ct,st,ca,sa,b,re_or_im);
    // add reflections of gradients with correct signs - careful!
    EPW_everything(val_s,ddx_s,ddy_s,dxx_s,dyy_s,dxy_s,-x,y,ct,st,ca,sa,b,\
		   re_or_im);
    val -= val_s;
    ddx += ddx_s;  ddy -= ddy_s;
    dxx -= dxx_s;  dyy -= dyy_s;  dxy += dxy_s;
    EPW_everything(val_s,ddx_s,ddy_s,dxx_s,dyy_s,dxy_s,x,-y,ct,st,ca,sa,b,\
		   re_or_im);
    val -= val_s;
    ddx -= ddx_s;  ddy += ddy_s;
    dxx -= dxx_s;  dyy -= dyy_s;  dxy += dxy_s;
    EPW_everything(val_s,ddx_s,ddy_s,dxx_s,dyy_s,dxy_s,-x,-y,ct,st,ca,sa,b,\
		   re_or_im);
    val += val_s;
    ddx -= ddx_s;  ddy -= ddy_s;
    dxx += dxx_s;  dyy += dyy_s;  dxy += dxy_s;
    break;

  case YO:
    yo_everything(val, ddx, ddy, dxx, dyy, dxy, x, y, s->a[j], s->b[j]);
    break;

  case YO_ODD_ODD:
  case YO_EVEN_EVEN:
    a = s->a[j]; // func origin
    b = s->b[j];
    yo_everything(val, ddx, ddy, dxx, dyy, dxy, x, y, a, b);
    yo_everything(val_s, ddx_s, ddy_s, dxx_s, dyy_s, dxy_s, x, y, -a, b);
    if (s->t[j]==YO_ODD_ODD) {
      val -= val_s;  ddx -= ddx_s;  ddy -= ddy_s;
      dxx -= dxx_s;  dyy -= dyy_s;  dxy -= dxy_s;
    } else {
      val += val_s;  ddx += ddx_s;  ddy += ddy_s;
      dxx += dxx_s;  dyy += dyy_s;  dxy += dxy_s;
    }
    yo_everything(val_s, ddx_s, ddy_s, dxx_s, dyy_s, dxy_s, x, y, a, -b);
    if (s->t[j]==YO_ODD_ODD) {
      val -= val_s;  ddx -= ddx_s;  ddy -= ddy_s;
      dxx -= dxx_s;  dyy -= dyy_s;  dxy -= dxy_s;
    } else {
      val += val_s;  ddx += ddx_s;  ddy += ddy_s;
      dxx += dxx_s;  dyy += dyy_s;  dxy += dxy_s;
    }
    yo_everything(val_s, ddx_s, ddy_s, dxx_s, dyy_s, dxy_s, x, y, -a, -b);
    val += val_s;  ddx += ddx_s;  ddy += ddy_s;
    dxx += dxx_s;  dyy += dyy_s;  dxy += dxy_s;
    break;

  default:
    fprintf(stderr,"Bad scaling function type in eval_basis_everything!\n");

  }
}


void dump_sfunc_geom(FILE *fp, int j, double k_typ, Basis_Set *s)
/*
  Writes basis sfunc type, followed by vector in gnuplot format x y dx dy,
  illustrating (for plotting) a useful sfunc vector.
  See show_geom.m for how this is can be used.
  
  Barnett 99/10/27
  Sep 03 modified for objects
*/
{
  sfunc_t t = s->t[j];
  double ca, b;
  
  fprintf(fp,"%d ",t);
  
  switch(t) {
    
  case RPW_RE:
  case RPW_IM:
  case RPW_ODD_ODD:
  case RPW_EVEN_EVEN:
  case RPW_REFL:
    fprintf(fp,"%g %g %g %g\n",0.0,0.0, \
	    s->nx[j], s->ny[j]);
    break;
    
  case EPW_RE_ODD_ODD:
  case EPW_IM_ODD_ODD:
    ca = s->c[j];
    b = s->b[j];
    fprintf(fp,"%g %g %g %g\n", b*s->ny[j]/k_typ, -b*s->nx[j]/k_typ, \
	    ca*s->nx[j], ca*s->ny[j]);
    break;
    
  case YO:
  case YO_ODD_ODD:
  case YO_EVEN_EVEN:
    fprintf(fp,"%g %g %g %g\n", s->a[j]/k_typ, s->b[j]/k_typ, s->nx[j], \
	    s->ny[j]);
    break;
    
  case FB:
    fprintf(fp,"%g %g\n", s->a[j], s->b[j]);
    break;
    
  default:
    fprintf(stderr,"Bad scaling function type in dump_sfunc_geom!\n");
    
  }
}



/* ===== NR-style allocation routines for sfunc type ===================== */

#define NR_END 1
#define FREE_ARG char*

sfunc_t *sfunc_tvector(long nl, long nh)
/* allocate an sfunc_t vector with subscript range v[nl..nh] */
{
        sfunc_t *v;

        v=(sfunc_t *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(sfunc_t)));
        if (!v) nrerror("allocation failure in sfunc_tvector()");
        return v-nl+NR_END;
}

void free_sfunc_tvector(sfunc_t *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}


/* ======================== BASIS SET OBJECT CREATE, BUILD DESTROY ======= */

void alloc_basis(Basis_Set *s, int N)
{
  s->N = N;              /* record N in basis set object */
  s->nx = dvector(1,N);
  s->ny = dvector(1,N);
  s->a = dvector(1,N);
  s->b = dvector(1,N);
  s->c = dvector(1,N);
  s->d = ivector(1,N);
  s->t = sfunc_tvector(1,N);
}

void free_basis(Basis_Set *s)
{
  int N = s->N;
  free_dvector(s->nx, 1,N);
  free_dvector(s->ny, 1,N);
  free_dvector(s->a, 1,N);
  free_dvector(s->b, 1,N);
  free_dvector(s->c, 1,N);
  free_ivector(s->d, 1,N);
  free_sfunc_tvector(s->t, 1,N);
}


/* ==================================== BUILD BASIS ======================== */

int build_basis(double k_typ, Billiard *l, Basis_Set *s)
/*
  Allocates and constructs a basis set given input params, billiard and k_typ.
  
  Rewritten for general scaling func type allowed for each basis element.
  The param[0...basis_num_params[basis]-1] array usually from cmd line.
  Barnett 99/10/28
  Basis construction allowed to be billiard- and k_typ-dependent, 99/11/1
  8/24/03: added oyooo set.
  3/5/06: added fourier-bessel re-entrant
*/
{
  int j, ji, N_eta, N_r, N_e, N = s->N, p;
  double N_sc, direction, alpha, beta, A, B, x, y, kd, corner_pack;
  double *param = s->param;                 /* basis set param array */

  s->k = k_typ;
  /* semiclassical basis size */
  N_sc = k_typ*l->Perim/PI;
  N_eta = round_and_clip(1.0 + param[0]*N_sc, "eta, in build_basis");
  // printf("build_basis: N_sc = %g,  N_eta = %d\n", N_sc, N_eta);

  switch (s->set_type) {

    // equally-distributed [0,pi] unit k-vectors, of RE and IM type...
  case RPWS:
    N = N_eta + ((N_eta%2) ? 1 : 0); /* make even */
    alloc_basis(s, N);
    for (j=1;j<=N/2;++j) {
      ji = j + N/2;
      s->t[j] = RPW_RE;				// assign RE sfunc type
      s->t[ji] = RPW_IM;		                // assign IM sfunc type
      direction = j * (2.0*PI/N);			// angle of k-vec
      s->nx[j] = s->nx[ji] = cos(direction);
      s->ny[j] = s->ny[ji] = sin(direction);
    }
    break;
	
    // equally-distributed [0,pi/2] unit k-vectors, odd-odd symm...
  case RPWS_ODD_ODD:
    N = N_eta;
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = RPW_ODD_ODD;
      direction = (j - 0.5) * PI/(2*N);			// angle of k-vec
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    break;

    // equally-distributed [0,pi/2] unit k-vectors, odd-odd symm...
  case RPWS_EVEN_EVEN:
    N = N_eta;
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = RPW_EVEN_EVEN;
      direction = (j - 0.5) * PI/(2*N);			// angle of k-vec
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    break;

    // N_r equally-distributed [0,pi/2] RPWS, N_e similar EPWS (const alpha)...
  case RPWS_EPWS_ODD_ODD:
    N_r = N_eta;
    N_e = round_and_clip(param[1], "N_e (RPWS_EPWS_ODD_ODD)");
    if (N_e % 2 != 0) {
      fprintf(stderr,\
	      "build_basis: evan basis size N_e must be even (N_e = %d)\n",\
	      N_e);
      return 1;
    }
    N = N_r + N_e;
    alloc_basis(s, N);
    // Real PWs...
    for (j=1;j<=N_r;++j) {
      s->t[j] = RPW_ODD_ODD;
      direction = (j - 0.5) * PI/(2*N_r);		// angle of k-vec
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    // Evanescent PWs of RE and IM types, and need alpha, b, arrays...
    for (j=N_r + 1; j <= N_r + N_e/2 ;++j) {
      ji = j + N_e/2;
      // type...
      s->t[j] = EPW_RE_ODD_ODD;
      s->t[ji] = EPW_IM_ODD_ODD;
      // cosh and sinh alpha...
      s->a[j] = s->a[ji] = sinh(param[2]);
      s->c[j] = s->c[ji] = cosh(param[2]);
      // k-vec, ie cos and sin theta, careful: in second quadrant!...
      direction = ( 0.5 + ((j - N_r) - 0.5)/N_e )*PI;
      s->nx[j] = s->nx[ji] = cos(direction);
      s->ny[j] = s->ny[ji] = sin(direction);
      // max impact parameter for billiard type, at given k-vec...
      s->b[j] = s->b[ji] = k_typ * max_billiard_impact(direction, l);
    }
    break;

    // PISTON1 qust basis: N_r biased [0,pi/2] RPWS,
    // N_e EPWS in theta=0 varying alpha
  case QUST_PISTON1_ODD_ODD:
    N_r = N_eta;
    N_e = round_and_clip(param[1], "N_e (QUST_PISTON1_ODD_ODD)");
    if (N_e % 2 != 0) {
      fprintf(stderr,\
	      "build_basis: evan basis size N_e must be even (N_e = %d)\n",\
	      N_e);
      return 1;
    }
    N = N_r + N_e;
    alloc_basis(s, N);
    A = param[2]; // alpha overall size
    B = param[3]; // alpha offset
    // Real PWs...
    for (j=1;j<=N_r;++j) {
      s->t[j] = RPW_ODD_ODD;
      // vergini k-vec dist...
      direction = 0.25*(5 - (j-1.)/N_r)*(j - 0.5) * PI/(2*N_r);
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    // Evanescent PWs of RE and IM types, and need alpha, b, arrays...
    for (j=N_r + 1; j <= N_r + N_e/2 ;++j) {
      ji = j + N_e/2;
      // type...
      s->t[j] = EPW_RE_ODD_ODD;
      s->t[ji] = EPW_IM_ODD_ODD;
      // cosh and sinh alpha...
      alpha = (B + j-N_r)/(A*pow(k_typ, 1./3));
      s->a[j] = s->a[ji] = sinh(alpha);
      s->c[j] = s->c[ji] = cosh(alpha);
      // k-vec, ie cos and sin theta, careful: in second quadrant!...
      direction = PI;
      s->nx[j] = s->nx[ji] = cos(direction);
      s->ny[j] = s->ny[ji] = sin(direction);
      // max impact parameter for billiard type, at given k-vec...
      s->b[j] = s->b[ji] = k_typ * max_billiard_impact(direction, l);
    }
    break;

    // Vergini 1/4 stadium basis: N_r biased [0,pi/2] RPWS,
    // N_e EPWS of varying alpha.
  case VERG_STAD_ODD_ODD:
    N_r = N_eta;
    N_e = round_and_clip(param[1], "N_e (VERG_STAD_ODD_ODD)");
    if (N_e % 2 != 0) {
      fprintf(stderr,\
	      "build_basis: evan basis size N_e must be even (N_e = %d)\n",\
	      N_e);
      return 1;
    }
    N = N_r + N_e;
    alloc_basis(s, N);
    beta = param[2]; // vergini balance param around kink pt
    // Real PWs...
    for (j=1;j<=N_r;++j) {
      s->t[j] = RPW_ODD_ODD;
      // vergini k-vec dist...
      direction = 0.25*(5 - (j-1.)/N_r)*(j - 0.5) * PI/(2*N_r);
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    // Evanescent PWs of RE and IM types, and need alpha, b, arrays...
    for (j=N_r + 1; j <= N_r + N_e/2 ;++j) {
      ji = j + N_e/2;
      // type...
      s->t[j] = EPW_RE_ODD_ODD;
      s->t[ji] = EPW_IM_ODD_ODD;
      // cosh and sinh alpha...
      alpha = (3.0 + j-N_r)/(2*pow(k_typ, 1./3));
      s->a[j] = s->a[ji] = sinh(alpha);
      s->c[j] = s->c[ji] = cosh(alpha);
      // k-vec, ie cos and sin theta, careful: in second quadrant!...
      direction = PI - sqrt(2/(k_typ * s->a[j]))/beta;
      s->nx[j] = s->nx[ji] = cos(direction);
      s->ny[j] = s->ny[ji] = sin(direction);
      // max impact parameter for billiard type, at given k-vec...
      s->b[j] = s->b[ji] = k_typ * max_billiard_impact(direction, l);
    }
    break;

  case OUTER_BIM_Y0:
    N = N_eta;
    kd = param[1];   // how outside it is in terms of typ wavenumbers.
    corner_pack = param[2];   // how many times more dense in corners.
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = YO;
      // get outer bdry location and direction...
      desc_outer_by_perim((j-0.5)/(double)N, kd/k_typ, corner_pack, &x, &y, \
		    &direction, l);
      s->a[j] = k_typ*x;
      s->b[j] = k_typ*y;
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    break;

  case OUTER_BIM_Y0_ODD_ODD:
    N = N_eta;
    kd = param[1];   // how outside it is in terms of typ wavenumbers.
    corner_pack = param[2];   // how many times more dense in corners.
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = YO_ODD_ODD;
      // get outer bdry location and direction...
      desc_outer_by_perim((j-0.5)/(double)N, kd/k_typ, corner_pack, &x, &y, \
		    &direction, l);
      // printf("\t\t\tj=%d: x,y = %g,%g\n", j, x, y);
      s->a[j] = k_typ*x;
      s->b[j] = k_typ*y;
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    break;

  case OUTER_BIM_Y0_EVEN_EVEN:
    N = N_eta;
    kd = param[1];   // how outside it is in terms of typ wavenumbers.
    corner_pack = param[2];   // how many times more dense in corners.
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = YO_EVEN_EVEN;
      // get outer bdry location and direction...
      desc_outer_by_perim((j-0.5)/(double)N, kd/k_typ, corner_pack, &x, &y, \
		    &direction, l);
      s->a[j] = k_typ*x;
      s->b[j] = k_typ*y;
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
    }
    break;

    // equally-distributed [0,pi/(2*p)] unit k-vectors, odd-odd symm...
  case RPWS_REFL:
    N = N_eta;
    p = round_and_clip(param[1], "p in RPWS_REFL");
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = RPW_REFL;
      direction = (j - 0.5) * PI/(2*p*N);	       	// angle of k-vec
      s->nx[j] = cos(direction);
      s->ny[j] = sin(direction);
      s->d[j] = p;
    }
    break;

  case FOURIER_BESSEL:
    A = param[1];    // start angle in units of pi
    B = param[2];    // end angle in units of pi
    // max rad specific to hamush bil, using corner of bounding box...
    // Note the semiclassical bump location of order nu bessel = nu
    N_sc = k_typ * (B-A) * sqrt(SQ(l->xl)+SQ(l->yh));
    N_eta = round_and_clip(1.0 + param[0]*N_sc, "FB eta, in build_basis");
    N = N_eta;
    alloc_basis(s, N);
    for (j=1;j<=N;++j) {
      s->t[j] = FB;              // type
      s->a[j] = PI*A;               // starting angle
      s->b[j] = j/(B-A);      // (ang q# * beta) parameter for each
    }
    break;

  default:
    fprintf(stderr,"Bad basis set type in build_basis!\n");
    return 1;
  }

  return 0;    /* success */
}



// ----------------------------------------------------------------------
void make_inhomog_basis(double k, Billiard *l, Basis_Set *s, double f,
			double kD)
  /* Make a single inhomogeneous source point using a N=1 OYO basis set.
   */
{
  int j = 1; // the one basis function
  double x, y, direction, corner_pack = 0.0;
  s->set_type = OUTER_BIM_Y0;
  alloc_basis(s, 1);
  s->k = k;
  s->t[j] = YO;
  desc_outer_by_perim(f, kD/k, corner_pack, &x, &y, &direction, l);
  s->a[j] = k*x;
  s->b[j] = k*y;
  s->nx[j] = cos(direction);
  s->ny[j] = sin(direction);
}



/* ================================== GENERATING SPATIAL VECs =============
 * NOTE: following code doesn't depend on basis set type or scaling funcs
 */

double eval_vec_spatial(double x, double y, double k, double *coeff, \
			Basis_Set *s)
/*
  Evaluates a single wavefunction at given x,y, given k-value and coeffs.
  The scaling by k is done here and in following routines.
 */
{
  int j, N = s->N;
  double psi;
  
  /* sum the basis funcs mult'd by their coeffs... */
  psi = 0.0;
  for (j=1;j<=N;++j)
    psi += coeff[j] * eval_basis(k*x, k*y, j, s);
  
  return psi;
}


void eval_vec_grad_spatial(double &ddx, double &ddy, double x, double y, \
		      double k, double *coeff, Basis_Set *s)
/*
  evaluates a single wavefunc gradient at given x,y, given k-value, coeffs.
  gradient is returned in pass-by-value (ddx,ddy). Barnett 99/12/2.
*/
{
  int j, N = s->N;
  double ddx_s, ddy_s;
  
  // sum the basis funcs gradients mult'd by their coeffs...
  ddx = ddy = 0.0;
  for (j=1;j<=N;++j) {
    eval_basis_deriv(ddx_s, ddy_s, k*x, k*y, j, s);
    ddx += coeff[j] * ddx_s;
    ddy += coeff[j] * ddy_s;
  }
  // finally scale since grad of a scaling basis...
  ddx *= k;
  ddy *= k;
}


void eval_vec_everything_spatial(double &val, double &ddx, double &ddy, \
				 double &dxx, double &dyy, double &dxy, \
				 double x, double y, double k, double *coeff, \
				 Basis_Set *s)
/*
  Evaluates a single wavefunc at given x,y, given k-value and coeffs.
  Returns value, 1st and 2nd derivs in pass-by-value arguments.
  1/16/04 Barnett
*/
{
  int j, N = s->N;
  double val_s, ddx_s, ddy_s, dxx_s, dyy_s, dxy_s, cj;
  
  // sum the basis func everythings mult'd by their coeffs...
  val = ddx = ddy = dxx = dyy = dxy = 0.0;
  for (j=1; j<=N; ++j) {
    eval_basis_everything(val_s, ddx_s, ddy_s, dxx_s, dyy_s, dxy_s, k*x, k*y, \
			  j, s);
    cj = coeff[j];
    val += cj * val_s;
    ddx += cj * ddx_s;
    ddy += cj * ddy_s;
    dxx += cj * dxx_s;
    dyy += cj * dyy_s;
    dxy += cj * dxy_s;
  }
  // finally scale gradients according to scaling ...
  ddx *= k;
  ddy *= k;
  dxx *= k*k;
  dyy *= k*k;
  dxy *= k*k;
}


void eval_multi_vecs_spatial(double x, double y, double k, int nv, \
	double **coeff_matrix, double *out, Basis_Set *s)
/*
  Efficient routine which calcs all nv spatial wavefuncs at a single point,
  given a SINGLE overall k-value. Therefore the boundary is not fixed,
  rather it scales according to the k_mu of the state mu. The cost does not
  scale with nv; computing many wavefuncs is same as computing one.
  Only values are computed (no derivs).
  
  Input is 1-indexed C matrix, output is 0-indexed float vector.
  Barnett 99/8/28
  10/30/03 output changed to 1-indexed double
 */
{
  int i, j, N = s->N;
  double psi;
  
  /* zero the output... */
  for (i=1;i<=nv;++i)
    out[i] = 0.0;
  
  /* sum the basis funcs mult'd by their coeffs, eval the basis only once... */
  for (j=1;j<=N;++j) {
    psi = eval_basis(k*x, k*y, j, s);
    for (i=1;i<=nv;++i)
      out[i] += coeff_matrix[i][j] * psi;
  }

}



/* --------------------- READ BASIS CMD LINE ARGS ------------------------ */

int parse_basis(char *st, Basis_Set *s)
  /*
   * Read string from command line into Basis_Set object type and
   * param array. Returns -1 if failure.
   *
   * barnett 9/10/03
   * 11/16/03 changed to colon-separated
   */
{
  int i, n;
  char tag[LEN];

  /* Identify what type of basis set we have, call it t... */
  basis_t t = NBASES;
  sscanf(st, "%[^:]%n", tag, &i);
  if (verb&0x80)
    printf("tag=%s, i=%d\n", tag, i);
  st += i;
  if (st[0]==':')
    ++st;
  // st now points to 2nd colon-separated item, or \0
  for(i=0; i<NBASES; ++i)
    if (strcmp(tag, basis_tags[i]) == 0)
      t = (basis_t)i;
  if (t==NBASES) {
    fprintf(stderr, "Invalid basis type: %s!\n", tag);
    return -1;
  }
  s->set_type = t;
  if (verb)
    printf("\tBasis type: %s\n", basis_desc[(int)t]);

  /* Check for param array overrun... */
  n = basis_num_params[(int)t];
  if (n > MAX_N_BASIS_PARAMS) {
    fprintf(stderr, "basis_num_params[%d]=%d exceeds max!\n", (int)t, n);
    return -1;
  }
  /* Read in params (all doubles) to param array... */
  if (separated_numbers(st, s->param, n) != n)
    return -1;

  /* feedback... */
  if (verb) {
    printf("\t\t%s\n\t\t", basis_param_desc[(int)t]);
    for (i=0; i<n; ++i)
      printf("%g ", s->param[i]);
    printf("\n\n");
  }
  return n;
}

void dump_basis(Basis_Set *s, char *head)
{
  FILE *fp;
  int j;
  char mname[LEN];

  sprintf(mname, "%s.b_pts", head);
  if ((fp = fopen(mname, "w")) == NULL)
    printf("file: %s, write error...\n", mname);

  fprintf(fp,"# %s from vergini\n\n", mname);
  fprintf(fp,"# basis sfunc type, followed by gnuplot vec: (x y dx dy)\n\n");
  for (j=1; j<=s->N; ++j)
    dump_sfunc_geom(fp, j, s->k, s);

  fprintf(fp,"\n");
  fclose(fp);
}


void show_basis_properties(Billiard *l, Basis_Set *s)
{
  double N_sc = l->Perim*s->k/PI;
  printf("	Basis: N = %d (%.3g of the semiclassical size %d)\n\n",\
	 s->N, s->N/N_sc, (int)(0.5 + N_sc));
}



/* ------------------------- SHOW BASIS ARGS USAGE ----------------------- */
void show_basis_usage()
{
  int i;

  fprintf(stderr, "\nPossible BASIS SETS, with param input formats.\n");
  for(i=0; i<NBASES; ++i) {
    fprintf(stderr, "    %s\n", basis_desc[i]);
    fprintf(stderr, "\t%8s %s\n", basis_tags[i], basis_param_desc[i]);
  }
}
