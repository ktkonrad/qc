
#ifndef BNC_COLLIDE_H
#define BNC_COLLIDE_H 1

#define MAX_N_SEGS 20

/* safety tolerance so that don't get double bounce at same object...*/
#define SMALL_DIST 1e-8
/* distance which is bigger than any distance within billiard... */
#define LARGE_DIST 100



extern int verb;            // globals


/* objects ............................................................... */


struct Segment {             /* segment object */
  int type;                   /* 0 for lineseg, 1 for circ, 2 for arc, etc */
  double l;                   /* its arclength */
  double q0;                  // what q coord on whole perimeter it starts at
  double x0, y0, x1, y1;      /* start, end coords. */
  double cx, cy;            /* center of circ */
  double R, inv_R, R2;        /* radius info */
  double inv_l;              /* length info */ 
  double ax, ay, bx, by, th;     /* arc angle info */
  double nx, ny;               // outwards normal for lines
};


/* Module PROVIDES ...................................  */
extern double circ_coll(double x, double y, double vx, double vy, double cx,\
			double cy, double R2, double *nx, double *ny);
extern double arc_coll(double x, double y, double vx, double vy, double cx,
		double cy, double R2, double inv_R, double ax, double ay,
		       double bx, double by, double *nx, double *ny);
extern double line_coll(double x, double y, double vx, double vy, double rox,
			double roy, double px, double py, double *nx, \
			double *ny);
extern double lineseg_coll(double x, double y, double vx, double vy,
			   double x1, double y1, double x2, double y2,
			   double inv_l);
extern double obj_coll(double x, double y, double vx, double vy,
		       struct Segment *s, double *nx, double *ny);
extern double one_bounce(double *x, double *y, double *vx, double *vy, \
			 int *obj, int ns, struct Segment *s);



#endif /* BNC_COLLIDE_H */
