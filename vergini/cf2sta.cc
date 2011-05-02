/* Generate .sta_bin from .sum and .cf file

   Reads cmd line in .sum to get billiard system params, .cf for coeff data

   barnett 8/25/03
*/

#include <stdio.h>
#include <string.h>

#include "./bdry.h"
#include "nrutil.h"

/* window types... */
#define NON_BOX 0
#define BOX 1

/* string lengths... */
#define LEN 256

/* GLOBALS: */

char filehead[100]; /* dummy */

/* basis set "private" data... */
int N;
double *b_nx,*b_ny, *b_a, *b_b, *b_c;
sfunc_t *b_sfunc_t;

/* billiard properties... note these are only stored for ONE param choice. */
double Perim;
double Area;

// overall billiard translation hack (20/8/01):
double trans_x, trans_y;

// This should eventually become local to main...
double billiard_param[MAX_NPARAMS];



/* -------------------------------------------- MAIN --------------------- */
int main(int argc, char **argv)
{
  int n_e, i, j, file_i, window, cmdargc, argp=1, n_params_dummy;
  double k_0, *k, **cfs;
  double gridsize, border, x_lo, x_hi, y_lo, y_hi, dum;
  char head[LEN], cf[LEN], sum[LEN], cmdline[LEN], dummy[LEN];
  char *cmdargv[LEN];
  FILE *fp_cf, *fp_sum;
  char *ptr, spacer[]=" ", delim[]=" \n";

  /* system... */
  basis_t basis;
  billiard_t billiard;
  opmode_t opmode;
  double basis_param[MAX_NPARAMS];
  double opmode_param[MAX_NPARAMS];

  if (argc!=4 && argc!=7) {
      printf("Usage:\n\tcf2sta filehead gridsize border\n");
      printf("   or:\n\tcf2sta filehead gridsize x_lo x_hi y_lo y_hi\n");
      return 1;
  }
  if (argc==7)
    window = BOX;
  else
    window = NON_BOX;
 
  /* open up the two files... */
  sscanf(argv[argp++], "%s", head);
  strcpy(sum, head);
  strcat(sum, ".sum");
  if ((fp_sum = fopen(sum,"r")) == NULL) {
    printf("cf2sta: summary file %s not found!\n", sum);
    return 1;
  }
  strcpy(cf, head);
  strcat(cf, ".cf");
  if ((fp_cf = fopen(cf,"r")) == NULL) {
    printf("cf2sta: coefficient file %s not found!\n", cf);
    return 1;
  }
  
  /* parse cf2sta cmd line... */
  sscanf(argv[argp++], "%lf", &gridsize);
  if (window == BOX) {
    sscanf(argv[argp++], "%lf", &x_lo);
    sscanf(argv[argp++], "%lf", &x_hi);
    sscanf(argv[argp++], "%lf", &y_lo);
    sscanf(argv[argp++], "%lf", &y_hi);
  } else {
    sscanf(argv[argp++], "%lf", &border);
  }
	     
  /* extract original called cmd line argv from .sum... */
  for (i=0;i<4;++i)
    if (fgets(cmdline, LEN, fp_sum) == NULL) {
      printf("cf2sta: error reading summary file header!\n");
      return 1;
    }
  for (i=0;i<3;++i)
    if (fgets(dummy, LEN, fp_sum) == NULL) {
      printf("cf2sta: error reading summary file header!\n");
      return 1;
    }

  /* read .sum calling command line and output feedback to user... */
  /* printf("%s.sum cmdline: %s", head, cmdline); */
  strtok_r(cmdline, delim, &ptr); /* take leading #, cmdline now workspace */
  for (cmdargc=0; (cmdargv[cmdargc]=strtok_r(NULL, delim, &ptr))!=NULL;
       ++cmdargc) ;
  /* reconcatenate cmd line for output... */
  dummy[0] = '\0';
  for(i=0;i<cmdargc;++i) {
    strcat(dummy, cmdargv[i]);
    strcat(dummy, spacer);
  }
  printf("%s.sum reconcat cmdline (cmdargc=%d):\n%s\n\n", head, cmdargc,\
	 dummy);
  if (!cmd_line_read(cmdargc, cmdargv, basis, billiard, opmode,\
		     basis_param, billiard_param, opmode_param, dummy))
    return(1);

  /* set up billiard: */
  // translate corner away from origin is appropriate...
  if (billiard==QU_STADIUM_ALL) {
    trans_x = -0.95;
    trans_y = -0.45;
  } else if (billiard==DESYM_SINAI_ALL) {
    trans_x = -0.25;
    trans_y = -0.30;
  } else {
    trans_x = 0.0;
    trans_y = 0.0;
  }
  Perim = perimeter(billiard);
  Area = area(billiard);
  k_0 = opmode_param[0];
  build_basis(basis_param, k_0, basis, billiard);
  printf("	built N = %d\n\n",N);


  /* read .cf data... */
  fscanf(fp_cf, "1 %d %d %lf %lf %d %lf %lf",&i,&n_e,&dum,&dum,
	     &n_params_dummy,&dum,&dum);
  if (i!=N) {
    printf("cf2sta: N in cf (%d) doesn't match built basis N!", i);
    return 1;
  }
  k = dvector(1, N);  /* alloc coeff and k arrays... */
  cfs = dmatrix(1, n_e, 1, N);

  for (i=1; i<=n_e; ++i) {
    fscanf(fp_cf,"%d %lf",&file_i, &dum);
    if (file_i!=i) {
      printf("cf2sta: file_i (%d) not match i=%d!\n", file_i, i);
      return 1;
    } else {
      for (j=1; j<=N; ++j)
	fscanf(fp_cf, "%lf ", cfs[i] + j);
    }
  }
  printf("%s.cf read (n_e = %d)\n", head, n_e);
  fclose(fp_cf);

  /* read .sum data... (just k values) */
  for (i=1; i<=n_e; ++i)
    fscanf(fp_sum,"%d %lf %lf %lf %lf", &file_i, k+i, &dum, &dum, &dum);
  printf("%s.sum read.\n", head);
  fclose(fp_sum);

  /* generate and output sta_bin data set... */
  if (window==BOX) {
    printf("box mode, writing %s.sta_bin...\n", head);
    eval_fast_save_multi_vecs(n_e, cfs, k, k_0, x_lo, x_hi, y_lo, y_hi, \
			      gridsize, 0, billiard, head);
  } else {
    double bil_xmin, bil_xmax, bil_ymin, bil_ymax;
    bounding_box(&bil_xmin, &bil_xmax, &bil_ymin, &bil_ymax, billiard);
    printf("non-box mode (bounds: %g %g %g %g), writing %s.sta_bin...\n",\
	   bil_xmin, bil_xmax, bil_ymin, bil_ymax, head);
    eval_fast_save_multi_vecs(n_e, cfs, k, k_0, bil_xmin-border,\
			      bil_xmax+border,\
			      bil_ymin-border, bil_ymax+border, \
			      gridsize, (border==0.0), billiard, head);
  } 

  free_dvector(k, 1, N);
  free_dmatrix(cfs, 1, n_e, 1, N);

  return 0;
}

