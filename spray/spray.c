/* spray.c  -  collect returning orbits by spraying from bdry,
 *             general simply-connected geometry of segments (lines, arcs)
 *
 * compile: gcc spray.c collide.o -o spray -lm
 *
 * current supported segment types are 0, +-2.  isolated circles and lines
 * will bounce ok but not get corret sense or correct q values.
 *
 * barnett 7/3/06, PO search added 8/17/06
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "collide.h"
#include "useful.h"

int verb;  // global verbosity

double parse_closed_curve(FILE *fp, int ns, struct Segment *s, double x0,
			  double y0, int ccw)
/* read geom file (segment entries) into segment object array, compute needed
 * extra lookup info in segment objects, output shape's perimeter.
 * (x0,y0) is initial point from which the curve starts. ccw = +-1, overall
 * sense of curve
 */
{
  int j;
  double theta, L=0.0, x=x0, y=y0;
  struct Segment *sj;

  for (j=0;j<ns;++j) {
    sj = s+j;           // pointer to current segment
    fscanf(fp, "%d", &sj->type);
    sj->x0 = x; sj->y0 = y;
    switch(sj->type) {
    case 0:           // lineseg
      fscanf(fp, "%lf %lf", &(sj->x1),&(sj->y1));
      sj->l = sqrt(SQ(sj->x1-sj->x0) + SQ(sj->y1-sj->y0));
      sj->inv_l = 1.0/sj->l;
      sj->nx = -ccw * sj->inv_l * (sj->y1-sj->y0);  // inwards normal
      sj->ny = ccw * sj->inv_l * (sj->x1-sj->x0);
      x = sj->x1; y = sj->y1;        // use end as next start
      break;
    case 1:           // circ -      not handled properly yet since not single
      fscanf(fp, "%lf %lf %lf", &(sj->cx),&(sj->cy),&(sj->R));
      sj->R2 = SQ(sj->R);
      sj->inv_R = 1.0/sj->R;
      sj->l = 2*PI*sj->R;
      break;
    case 2:           // arc ccw
    case -2:           // arc cw
      fscanf(fp, "%lf %lf %lf %lf", &(sj->x1),&(sj->y1),&(sj->cx),&(sj->cy));
      sj->R = sqrt(SQ(sj->x0-sj->cx) + SQ(sj->y0-sj->cy));
      sj->R2 = SQ(sj->R);
      sj->inv_R = 1.0/sj->R;
      if (sj->type==-2) {                                   // cw
	sj->bx = sj->x0-sj->cx;  sj->by = sj->y0-sj->cy;
	sj->ax = sj->x1-sj->cx;  sj->ay = sj->y1-sj->cy;
      } else {                                              // ccw
	sj->ax = sj->x0-sj->cx;  sj->ay = sj->y0-sj->cy;
	sj->bx = sj->x1-sj->cx;  sj->by = sj->y1-sj->cy;
      }
      sj->th = atan2(sj->ax*sj->by - sj->bx*sj->ay,      // R^2 sin theta
		    sj->ax*sj->bx + sj->ay*sj->by);     // R^2 cos theta
      if (sj->th<0.0)
	sj->th += 2*PI;
      sj->l = sj->th * sj->R;
      x = sj->x1; y = sj->y1;             // use end as next start
      break;
    default:
      fprintf(stderr, "Bad segment type (%d for j=%d) in read_geom!\n",
	      sj->type, j);
      exit(1);
    }
    sj->q0 = L;          // set start q of segment
    L += sj->l;
  }
  if (sqrt(SQ(x-x0)+SQ(y-y0))>1e-15)
    fprintf(stderr, "parse_closed_curve: (%g,%g) not closing to (%g,%g)!\n",
	    x, y, x0, y0);
  return L;
}

int tell_geom(FILE *fp, int ns, struct Segment *s)
// output a segment's info to file pointer
{
  switch (s->type) {
  case 0:
    fprintf(fp, "(0) lineseg: (%g,%g) to (%g,%g), l=%g\n",
	    s->x0, s->y0, s->x1, s->y1, s->l);
    break;
  case 1:
    fprintf(fp, "(1) circ:    c (%g,%g) r %g, l=%g\n",
	    s->cx, s->cy, s->R, s->l);
    break;
  case 2:
    fprintf(fp, "(2) arc ccw: (%g,%g) to (%g,%g), c=(%g,%g) r=%g, th=%g*PI, l=%g\n",
	    s->x0, s->y0, s->x1, s->y1, s->cx, s->cy, s->R, s->th/PI, s->l);
    fprintf(fp, "\t\ta=(%g,%g) b=(%g,%g)\n", s->ax,s->ay, s->bx, s->by);
    break;
  case -2:
    fprintf(fp, "(2) arc cw:  (%g,%g) to (%g,%g), c=(%g,%g) r=%g, th=%g*PI, l=%g\n",
	    s->x0, s->y0, s->x1, s->y1, s->cx, s->cy, s->R, s->th/PI, s->l);
    fprintf(fp, "\t\ta=(%g,%g) b=(%g,%g)\n", s->ax,s->ay, s->bx, s->by);
    break;
  default:
    fprintf(stderr, "Bad segment type (%d) in obj_coll!\n", s->type);
    return 1;
  }
  return 0;
}

int find_bdry_coords(double q, int ccw, int ns, struct Segment *s, double *x, \
		     double *y, double *nx, double *ny)
/* returns (x,y) and (nx,ny) outwards normal given 0<q<L the bdry coord.
 * ns and s give segment objects. ccw  = +-1 gives orientation sense of curve.
 * return value is segment obj # hit: 0,1,...
 * Does not need to be super fast. All normals are inwards-pointing.
 */
{
  int j;
  double d, ct, st, l_tot=0.0, dx, dy;
  struct Segment *sj;
  // find which segment it's in
  for (j=0; l_tot+(s+j)->l<q; l_tot+=(s+j)->l, ++j) ;
  d = q-l_tot; // dist along current segment
  sj = s+j;
  switch (sj->type) {
  case 0:
    *x = (d*sj->x1 + (sj->l-d)*sj->x0) * sj->inv_l; 
    *y = (d*sj->y1 + (sj->l-d)*sj->y0) * sj->inv_l;
    *nx = sj->nx;
    *ny = sj->ny;
    break;
  case 1:
    *x = sj->cx + sj->R*cos(d/sj->R); 
    *y = sj->cy + sj->R*sin(d/sj->R); 
    break;
  case 2:
    ct = cos(d/sj->R);  st = sin(d/sj->R);
    dx = ct*sj->ax - st*sj->ay;                     // rotate a vec ccw
    dy = st*sj->ax + ct*sj->ay;
    *x = sj->cx + dx;
    *y = sj->cy + dy;
    *nx = -ccw * sj->inv_R * dx;
    *ny = -ccw * sj->inv_R * dy;
    break;
  case -2:
    ct = cos(d/sj->R);  st = sin(d/sj->R);
    dx = ct*sj->bx + st*sj->by;                     // rotate b vec cw
    dy = -st*sj->bx + ct*sj->by;
    *x = sj->cx + dx;
    *y = sj->cy + dy;
    *nx = ccw * sj->inv_R * dx;
    *ny = ccw * sj->inv_R * dy;
    break;
  default:
    fprintf(stderr, "Bad segment type (%d) in find_bdry_coords!\n", sj->type);
    return -1;
  }
  return j;
}

int show_usage()
{
  fprintf(stderr, "\nusage: spray [options]\n\nOptions: (output to stdout unless indicated)\n\n");
  fprintf(stderr, "\t-m M\t\tdump geom as M bdry pts\n");
  fprintf(stderr, "\t-S q:N:dq:lmax\tspray N traj's from bdry loc q\n");
  fprintf(stderr, "\t-T q:p:N\toutput n bounces of single traj from (q,p)\n");
  fprintf(stderr, "\t-P q:N:dq:dp:lmax\tspray N traj's from q, collect POs\n");
  exit(1);
}


/* =================================== MAIN =============================== */
int main(int argc, char **argv)
{
  int j, ns, ccw, M=0, N=0, n, verb=1, obj;
  double x0, y0, x, y, nx, ny, vx, vy, q=0.0, p=0.0; // launch location on bdry
  double l, t, l_max, xi, yi, vxi, vyi, dq, dp=0.0, dr2, dv2, cost;
  struct Segment *s;               /* pointer to array of segments */
  FILE *fp;

  /* read cmd line */
  while ((j = getopt(argc, argv, "qm:P:S:T:")) != -1)
    switch (j) {
    case 'm':              // dump geom as bdry pts
      sscanf(optarg, "%d", &M);
      break;
    case 'q':              // quiet
      verb = 0;
      break;
    case 'P':              // periodic orbits q:N:dq:dp:lmax
      sscanf(optarg, "%lf:%d:%lf:%lf:%lf", &q, &N, &dq, &dp, &l_max);
      break;
    case 'S':              // spray q:N:dq:lmax
      sscanf(optarg, "%lf:%d:%lf:%lf", &q, &N, &dq, &l_max);
      break;
    case 'T':              // single trajectory q:p:n 
      sscanf(optarg, "%lf:%lf:%d", &q, &p, &n);
      N=1;
      break;
    case 'h':
    default:
      show_usage();
    }
  if (M==0 && N==0) { // no action was selected
    show_usage();
  }

  // read geom and compute segment params.
  scanf("%d %d %lf %lf", &ccw, &ns, &x0, &y0); // sense, #segs, start location
  if (ns<1) {
    fprintf(stderr, "spray: ns<1 (ns=%d)!\n", ns);
    exit(1);
  }
  if (ccw!=+1 && ccw!=-1) {
    fprintf(stderr, "spray: ccw must be +1 or -1 (ccw=%d)!\n", ccw);
    exit(1);
  }
  s = (struct Segment *)malloc(ns * sizeof(struct Segment));
  double L = parse_closed_curve(stdin, ns, s, x0, y0, ccw);
  if (verb) {
    for (j=0;j<ns;++j) {
      fprintf(stderr, "segment %d:\t", j);
      tell_geom(stderr, ns, s+j);
    }
    fprintf(stderr, "perim = %g\n", L);
  }

  if (M>0) {    // geom:  output regular list of x,y,nx,ny for q in [0,L] ---
    //if ((fp = fopen(mname, "w")) == NULL)
    //  printf("spray: file %s, write error...\n",mname);
    printf("%.15g\n", L);
    for (q=0.;q<L;q+=L/M) {
      obj = find_bdry_coords(q, ccw, ns, s, &x, &y, &nx, &ny);
      printf("%.15g %.15g %.15g %.15g\n", x, y, nx, ny);
    }
  }

  if (N>0 && (q<0 || q>L)) {
    fprintf(stderr, "q=%g is out of bdry range [0,%g]!\n", q, L);
    exit(1);
  }
  if (N==1) {  // single trajectory -----------------------------------------
    if (verb)
      fprintf(stderr, "Task T: single traj q=%g,p=%g,n=%d\n", q, p, n);
    if (p<-1 || p>1) {
      fprintf(stderr, "p=%g is out of range [-1,1]!\n", p);
      exit(1);
    }
    obj = find_bdry_coords(q, ccw, ns, s, &x, &y, &nx, &ny);
    cost = sqrt(1.0-p*p);        // cosine of launch angle
    vx = cost*nx + ccw*p*ny;       // p sense: p>0 is along incr q.
    vy = -ccw*p*nx + cost*ny;
    l = 0.0;                       // total traj length
    printf("%.15g %.15g %.15g %.15g %d %.15g\n", x, y, vx, vy, obj, l); // q,p?
    for (j=1;j<=n;++j) {
      l += one_bounce(&x, &y, &vx, &vy, &obj, ns, s);
      if (j%5==0) {     // expensive: normalize v vector since unstable growth
	double v = 1.0 / sqrt(vx*vx + vy*vy); vx *= v; vy *= v;
      }
      printf("%.15g %.15g %.15g %.15g %d %.15g\n", x, y, vx, vy, obj, l);
    }

  } else if (N>1) { // spray trajectories and collect returning ones ----------
    if (verb)
      if (dp==0.0)
	fprintf(stderr, "Task S: from q=%g spray %d trajs, dq=%g lmax=%g\n", \
		q, N, dq, l_max);
      else
	fprintf(stderr, \
		"Task P: from q=%g spray %d trajs, dq=%g dp=%g lmax=%g\n", \
		q, N, dq, dp, l_max);
    //for (t=-PI/2*(1-1.0/N); t<PI/2; t+=PI/N) {   // uniform spray in angle
    //p = sin(t); cost = cos(t);
    for (p=-1+1.0/N; p<1; p+=2.0/N) {
      obj = find_bdry_coords(q, ccw, ns, s, &xi, &yi, &nx, &ny);
      x = xi; y = yi;
      cost = sqrt(1.0-p*p);        // cosine of launch angle
      t = asin(p);
      vx = cost*nx + ccw*p*ny;       // p sense: p>0 is along incr q.
      vy = -ccw*p*nx + cost*ny;
      vxi = vx; vyi = vy;
      if (verb)
	fprintf(stderr, "\t(%g,%g) launch t=%g p=%g v=(%g,%g)\n",\
		x,y,t,p,vx,vy);
      if (dp==0.0)                     // check only if q close (task S)
	for (l=0.0,j=1; l<l_max; ++j) {       // bounce until l_max reached
	  l += one_bounce(&x, &y, &vx, &vy, &obj, ns, s);  // l increased here
	  if (j%5==0) {           // normalize v
	    double v = 1.0 / sqrt(vx*vx + vy*vy); vx *= v; vy *= v;
	  }
	  if ((dr2=SQ(xi-x) + SQ(yi-y))<dq*dq) {  // hack: close in R^2 not q!
	    dv2 = SQ(vxi-vx)+SQ(vyi-vy);
	    printf("%.15g %d %.15g %.8g %.8g\n", p, j, l, dr2, dv2);
	  }
	}
      else               // check if both q and p close (task P)
	for (l=0.0,j=1; l<l_max; ++j) {       // bounce until l_max reached
	  l += one_bounce(&x, &y, &vx, &vy, &obj, ns, s);  // l increased here
	  if (j%5==0) {           // normalize v
	    double v = 1.0 / sqrt(vx*vx + vy*vy); vx *= v; vy *= v;
	  }
	  if ((dr2=SQ(xi-x)+SQ(yi-y))<dq*dq&&(dv2=SQ(vxi-vx)+SQ(vyi-vy))<dp*dp)
	    printf("%.15g %d %.15g %.8g %.8g\n", p, j, l, dr2, dv2);
	}
    }
  }
  
  return 0;
}
