/* Billiard module for VERGINI
 *
 * barnett 9/10/03
 */

#ifndef VERGINI_BILLIARD_H
#define VERGINI_BILLIARD_H 1

#define MAX_N_BILLIARD_PARAMS 20

// max number of corners...
#define MAX_CORNERS  20

/* BILLIARD TYPES:
 * To add a billiard type, add to the next few lines, + add to the switches
 * in the routines in billiard.cc
 */
typedef enum
        {STADIUM, QU_STADIUM, QU_STADIUM_ALL, DESYM_SINAI,
	 DESYM_SINAI_ALL, GEN_RECT_SINAI, GEN_ASYM_SINAI,
	 QU_GEN_RECT_SINAI, GEN_SEG_FILE, RAD_FUNC, QU_RAD_FUNC,
	 LINE, ELL, HALF_MUSH, NBILLIARDS}
	billiard_t;

/* billiard description strings... */
static char *billiard_desc[] = {
	"Stadium",
	 "1/4 Stadium (for C4 sym basis)",
	 "1/4 Stadium (all-perim)",
	"Desym pi/4 Sinai (for C4 sym basis)",
	"Desym pi/4 Sinai (all-perim)",
	"Gen Rect Sinai",
	"Gen Asym Rect Sinai",
	"1/4 Gen Rect Sinai (for C4 sym basis)",
        "General Segmented from file",
	"Radial Func",
        "Quarter Radial Func",
        "Line",
	"L billiard",
        "Desym Half Mushroom"};
static char *billiard_tags[] =
	{"st","qust","qusta","desin","desina","grs","gas","qugrs","seg",
	 "rf", "qurf", "line", "ell", "hamush"};
static int billiard_num_params[] =
  {3,1,3,1,3,5,7,3,1,5,3,4,0,2};
static char *billiard_param_desc[] =
	{"alpha trans_x trans_y", "alpha", "alpha trans_x trans_y", "radius",
	 "radius trans_x trans_y",
	 "a theta1 theta2 trans_x trans_y",
	 "a theta1 theta2 theta3 theta4 trans_x trans_y",
	 "a theta1 theta2", "seg_file", "type ecc_a ecc_b trans_x trans_y",
	 "type ecc_a ecc_b", "x0 y0 x1 y1", "", "a_stem b_stem"};
/* Which param number (0-indexed) to read trans_x, trans_y from.
   (-1 if no translation)... */
static int billiard_trans_param[] =
  {1,-1,1,-1,1,3,5,-1,-1,3,-1,-1,-1,-1};


/* BILLIARD OBJECT (CONTINUOUS PROPERTIES; see Bdry_Pt_Set for collocation)
 * Prestore frequently-accessed or discrete properties. Others are computed
 * every time, eg the location given bdry coord s is given by desc_by_perim
 */
typedef struct Billiard {
  billiard_t type;
  double Perim;
  double Area;
  double trans_x, trans_y;                        // translation from default.
  double cf[MAX_CORNERS+1], cal[MAX_CORNERS+1];   // corner definition arrays
  double cah[MAX_CORNERS+1], cel[MAX_CORNERS+1];
  int nc;                                         // number of corners
  double xl, xh, yl, yh;                          // bounding box
  double param[MAX_N_BILLIARD_PARAMS];            // incoming params
  int nangles;                                    // # proto colloc pts
  double *angle, *arclen;                         // proto colloc data
  int neumann;                                    // Neumann BCs flag
};


// billiard.cc provides...
extern double perimeter(Billiard *l);
extern double area(Billiard *l);
extern int inside_billiard(double x, double y, Billiard *l);
extern void desc_by_perim(double f, double *cx, double *cy, \
	double *normal, double *curv, Billiard *l);
extern int give_corners(double *cf, double *cal, double *cah,\
			Billiard *l);
extern void desc_outer_by_perim(double f, double d, double corner_pack, \
				double *cx, double *cy, \
				double *normal, Billiard *l);
extern void bounding_box(double *xl, double *xh, double *yl, double *yh, \
			 Billiard *l);
extern double max_billiard_impact(double dir, Billiard *l);
extern void rad_func(Billiard *l, double t, double *r, double *rt,\
		     double *rtt);
extern int build_billiard(Billiard *l, double k_base);
extern int parse_billiard(char *s, Billiard *l);
extern void dump_billiard(Billiard *l, char *head);
extern void show_billiard_properties(Billiard *l, double k);
extern void show_billiard_usage();

#endif /* VERGINI_BILLIARD_H */
