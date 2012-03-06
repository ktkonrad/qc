/* verg.cc: compute dirichlet eigenproblem via VERGINI scaling method 
 *
 * barnett 11/17/03 based on original code written 1999-2000.
 * 1/16/04 full 2nd-deriv compliance, general BC area norms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>                 // Gnu Sci Lib special funcs
#include <gsl/gsl_errno.h>                 // Gnu Sci Lib error handler

/* local low-level utilities... */
#include "cxml.h"
#include "nrutil.h"
#include "useful.h"

/* local higher-level utilities... */
#include "billiard.h"
#include "basis.h"
#include "colloc.h"
#include "grid.h"
#include "verb.h"
#include "matrix.h"

/* define global verbosity... and set its default
   flags: 0x02 - generally more
          0x80 - cmd line parsing
*/ 
int verb = 1;

// tasks...
#define NO_TASK -1
#define TASK_BASISTEST 0
#define TASK_VERGINI 1
#define TASK_REFVERG 2
#define TASK_INHOMOG 3
#define TASK_EIGENVEC 4
#define TASK_NORMTEST 5
#define TASK_MASSCHOP 6

// ---------------------------------------- USAGE -------------------------
void show_usage()
  /* For verg.cc only
   */
{
  fprintf(stderr, "\nSolve Laplacian eigenproblem via scaling method.\n");
  fprintf(stderr, "Usage:\n\tverg [options] task\nOptions:\n");
  fprintf(stderr, "\t-v\tverbose\n\t-q\tquiet (no stdout), except for ");
  fprintf(stderr, "arguments earlier on cmd line\n");
  fprintf(stderr, "\t-b <B>\tset boundary collocation to B per wavelength\n");
  fprintf(stderr, "\t-k <K>\tset wavenumber to K\n");
  fprintf(stderr, "\t-l <type:a:b:...:z>\tdefine billiard shape (required)\n");
  fprintf(stderr, "\t-s <type:a:b:...:z>\tdefine basis set (required)\n");
  fprintf(stderr, "\t-g <dx:mask>\toutput 2d grids, spacing dx\n\t\t");
  fprintf(stderr, "mask = 0: no clip\n\t\tmask = 1: crude clip\n\t\t");
  fprintf(stderr, "mask > 1: gen weight-mask (output with -m), clip to it\n");
  fprintf(stderr, "\t-f <dx:mask>\tfast output 2d grids (at same k)\n");
  fprintf(stderr, "\t-x <xl:xh:yl:yh>\tdefine 2d output grid box\n");
  fprintf(stderr, "\t-p <Mo>\toutput bdry funcs, Mo samples (or Mo=0: use B)\n");
  fprintf(stderr, "\t-o <head>\theader name for output files\n");
  fprintf(stderr, "\t-m\tdump geometry files\n\t-h\tshow help\n");
  fprintf(stderr, "\t-e <s>\tencode grid output with colloc points\n");
  fprintf(stderr, "\t\ts = 0: don't rescale with k\n\t\ts = 1: do rescale\n");
  fprintf(stderr, "\t-4 <C4>\tset 4th-order vergini coeff (default 4.0)\n");
  fprintf(stderr, "\t-u\tremove spurious states\n");
  fprintf(stderr, "\t-2\tfor mass-chop use old 2nd-deriv overlap formula\n");
  fprintf(stderr, "\t-d\tdon't output coeffs vectors\n");
  fprintf(stderr, "\t-z <eps>\tdefine epsilon for matrix GEP truncation\n");
  fprintf(stderr, "\t-n\tneumann boundary conditions (experimental)\n");
  fprintf(stderr, "Tasks (choose one):\n\t-V <max_delta> or <delta_lo:delta_hi>\tsingle vergini\n");
  fprintf(stderr, "\t-R <delta_lo:delta_hi:delta_ref>\treference, 1 verg per state\n");
  fprintf(stderr, "\t-T <j_lo:j_hi>\ttest output range j of basis funcs\n");
  fprintf(stderr, "\t-N <seed>\ttest random wavefunction area 2-norm\n");
  fprintf(stderr, "\t-E <wei:j_lo:j_hi>\tevecs of operator (don't use with -g)\n");
  fprintf(stderr, "\t\t\tweight: wei=0 gives w=1, wei=1 gives w=1/rn\n");
  fprintf(stderr, "\t\t\tnegative wei uses vergini bdry operator\n");
  fprintf(stderr, "\t-I <f:kD>\tsolve inhomogeneous BVP with single src\n");
  fprintf(stderr, "\t-C <delta_lo:delta_hi:y_c:f_c\tsingle vergini with mass-chop\n");
  fprintf(stderr, "\t\tif bdry out (-p) use closest to same B for Mc as Mo\n");
  fprintf(stderr, "\nExample:\n\t./verg -v -l qugrs:1:0.2:0.4 -b 10 -s oyooo:3:3:1 -k 50 -e 1 -p 10 -V 0.2\n\n");
}


// --------------------------------------- SAVE SUMMARY --------------------
int save_sum(char *cmdline, double *ks, double *kos, double *ten, \
	     double *nrm, int ne, char *head)
  /* dump summary of eigenwavenumbers, quality indicators, to text file.
     Returns 0 if ok.
     barnett 11/21/03
  */
{
  int i;
  char name[LEN];
  FILE *fp;

  sprintf(name, "%s.sum", head);
  if ((fp = fopen(name, "w")) == NULL) {
    fprintf(stderr, "save_sum: %s, open for write error...\n", name);
    return 1;
  }

  fprintf(fp,"# Summary output from verg version 11/21/03\n");
  fprintf(fp,"# calling command line:\n# %s\n#\n", cmdline);	
  fprintf(fp,\
    "# Columns: $1=i, $2=k_i, $3=k0_i, $4=tension_i, $5=perimnorm_i\n#\n");
  
  for(i=1; i<=ne; ++i)
    fprintf(fp, "%d %.16g %.16g %.16g %.16g\n", i, ks[i], kos[i], ten[i],\
	    nrm[i]);
  if (verb)
    printf("summary file (%s) written.\n",name);
  fclose(fp);
  return 0;
}



/* ======================================= MAIN ======================== */
int verg_main(int argc, char **argv)
{
  Billiard bil;
  Bdry_Pt_Set bps, obps;
  Basis_Set bas;
  Grid g;
  int geom = 0, grid = 0, mask = 0, encode = 0, encode_res = 0, bdryfuncs = 0;
  int task = NO_TASK, remove_spurious = 0, sec_deriv_mc = 0, want_coeffs = 1;
  int i, j, j_lo = 1, j_hi = 5, ne, nev, ne_ref, oM, j_ref, seed, *idx;
  int weight = 0;
  double b = 10, k_base = 20, dx = 0.01, delta_lo, delta_hi, delta_ref, gap;
  double k0, xl = 0.0, xh = 0.0, yl = 0.0, yh = 0.0, C4 = 4.0, bo;
  double **a, *ks, *kos, **per, **ngr, *ten, *nrm, *ks_ref, **a_ref;
  double epsilon = 1e-16;   // default
  char cmdline[LEN], *head = "t";
  // For TASK INHOMOG:
  Basis_Set bas_inh;
  double **a_inh, f_inh, kD_inh;
  // For TASK_MASSCHOP:
  Billiard bil_chop;
  Bdry_Pt_Set bps_chop;
  char headchop[LEN];
  int m_c;
  double y_c, f_c, x_fc, y_fc, dummy, *mass, **per_chop, **ngr_chop, alpha;
  
  gsl_set_error_handler_off();       // turn off GNU Sci Lib error handler


  // parse command line....................................................
  grab_cmdline(cmdline, argc, argv);
  opterr = 0;
  optind = 0; // have to do this for getopt to parse an array other than the one from the command line
  bil.type = NBILLIARDS;
  bil.neumann = 0;   // default
  bas.set_type = NBASES;
  while ((i = getopt(argc, argv, \
		     "l:s:b:k:x:o:vhmqg:f:e:p:V:R:T:I:E:N:C:4:u2z:dn")) != -1)
    switch (i) {
    case 'o': // output file head
      head = optarg;
      break;
    case 'g':
    case 'f': // output 2d grid with dx (1=f=fast, 2=g=slow scaled)
      grid = (i=='f') ? 1 : 2;
      sscanf(optarg, "%lf:%d", &dx, &mask);
      break;
    case 'p': // output bdry funcs with oM samples
      bdryfuncs = 1;
      sscanf(optarg, "%d", &oM);
      break;
    case 'b':
      sscanf(optarg, "%lf", &b);
      break;
    case 'k':
      sscanf(optarg, "%lf", &k_base);
      break;
    case 'z':
      sscanf(optarg, "%lf", &epsilon);
      break;
   case 'x':
      sscanf(optarg, "%lf:%lf:%lf:%lf", &xl, &xh, &yl, &yh);
      break;
    case 'l':
      if (parse_billiard(optarg, &bil)==-1) {
	show_billiard_usage();
	show_usage();
	return 1;
      }
      break;
    case 's':
      if (parse_basis(optarg, &bas)==-1) {
	show_basis_usage();
	show_usage();
	return 1;
      }
      break;
    case 'q': // quiet
      verb = 0;
      break;
    case 'v': // very verbose
      verb = 127;
      break;
    case 'h': // help
      show_usage();
      return 0;
      break;
    case 'm': // dump geom
      geom = 1;
      break;
    case 'e': // encode
      encode = 1;
      sscanf(optarg, "%d", &encode_res);
      break;
    case '4': // C4 coeff
      sscanf(optarg, "%lf", &C4);
      break;
    case 'u': // spurious
      remove_spurious = 1;
      break;
    case '2': // 2nd-order derivs for masschop
      sec_deriv_mc = 1;
      break;
    case 'd': // don't want coeffs
      want_coeffs = 0;
      break;
    case 'n':
      bil.neumann = 1;
      if (verb)
	printf("Settting Neumann boundary conditions...\n");
      break;

    case 'T': // TASK basis test ------------------------ TASKS --------------
      task = TASK_BASISTEST;
      sscanf(optarg, "%d:%d", &j_lo, &j_hi);
      break;
    case 'N': // TASK norm test
      task = TASK_NORMTEST;
      sscanf(optarg, "%d", &seed);
      break;
    case 'E': // TASK quadratic form eigenfunctions
      task = TASK_EIGENVEC;
      sscanf(optarg, "%d:%d:%d", &weight, &j_lo, &j_hi);
      break;
    case 'C': // TASK single vergini mass chop
      task = TASK_MASSCHOP;
      sscanf(optarg, "%lf:%lf:%lf:%lf", &delta_lo, &delta_hi, &y_c, &f_c);
      break;
    case 'V': // TASK single vergini, with possible asymm window
      task = TASK_VERGINI;
      if (strchr(optarg, ':')!=NULL)
	sscanf(optarg, "%lf:%lf", &delta_lo, &delta_hi);
      else {
	sscanf(optarg, "%lf", &delta_lo);
	delta_hi = delta_lo;
      }
      if (verb)
	printf("Vergini: range k-%g to k+%g\n", delta_lo, delta_hi);
      break;
    case 'R': // TASK gather vergini
      task = TASK_REFVERG;
      sscanf(optarg, "%lf:%lf:%lf", &delta_lo, &delta_hi, &delta_ref);
      if (verb)
	printf("Reference Vergini: range k-%g to k+%g, delta_ref=%g\n", \
	       delta_lo, delta_hi, delta_ref);
      break;
    case 'I': // TASK Inhomogeneous BVP with single src
      task = TASK_INHOMOG;
      sscanf(optarg, "%lf:%lf", &f_inh, &kD_inh);
      break;
    }
  if (bil.type==NBILLIARDS) {
    fprintf(stderr, "No billiard defined!\n");
    show_usage();
    show_billiard_usage();
    return 1;
  }
  if (bas.set_type==NBASES) {
    fprintf(stderr, "No basis defined!\n");
    show_usage();
    show_basis_usage();
    return 1;
  }
  if (task==NO_TASK) {
    fprintf(stderr, "No task defined!\n");
    show_usage();
    return 1;
  }
  // .............................................................

  build_billiard(&bil, k_base);
  build_bdry(b, k_base, &bil, &bps);
  build_basis(k_base, &bil, &bas);
  if (verb) {
    printf("\nverg: scaling method for Laplacian eigenproblem, 8/9/04\n\n");
    if (verb&0x80)
      printf("cmdline: %s\n\n", cmdline);
    show_billiard_properties(&bil, k_base); 
    show_colloc_properties(&bil, &bps, k_base); 
    show_basis_properties(&bil, &bas);
    printf("Output file head: %s\n\n", head);
    if (grid && (xl-xh)!=0.0)
      printf("\toutput box: x=[%f,%f] y=[%f,%f]\n", xl, xh, yl, yh);
  }
  if (geom) {
  // output geom...
    if (verb)
      printf("dumping geometry files...\n");
    dump_billiard(&bil, head);
    dump_colloc(&bps, head);
    dump_basis(&bas, head);
  }

  // choose ne for allocation, depending on required task...
  switch (task) {
  case TASK_NORMTEST:
    ne = 1;
    break;
  case TASK_BASISTEST:
  case TASK_EIGENVEC:
    ne = j_hi - j_lo + 1;
    break;
  case TASK_VERGINI:
  case TASK_REFVERG:
  case TASK_MASSCHOP:
    // estimate number of states using Weyl, + 50 extra for spurious
    ne = (int)(2*(delta_hi+delta_lo)*k_base*bil.Area/(2*PI) + 50);
  break;
  case TASK_INHOMOG:
    ne = 2; // single function being solved for but 2 sets of per needed.
  }
  if (verb)
    printf("allocating for ne=%d states.\n", ne);
  ks = dvector(1, ne); // wavenumbers k
  kos = dvector(1, ne); // k_0 used for each state
  idx = ivector(1, ne); // index for juggling nonspurious states
  a = dmatrix(1, ne, 1, bas.N); // coeffs

  // set up bdry funcs on an 'obps' (output bdry pt set) with oM samples
  if (!bdryfuncs) {
    // if no output needed, use standard b=6 for norming/spurious-test...
    oM = round_and_clip(6.0 * k_base*bil.Perim / (2*PI), "oM in verg");
    bo = (2*PI*oM)/(k_base*bil.Perim); // bo = b_output needed for masschop
  } else if (oM==0) {
    oM = bps.M;
    bo = b;
  } else
    // Mo is given by cmd line...
    bo = (2*PI*oM)/(k_base*bil.Perim);
  build_bdry(-(double)oM, k_base, &bil, &obps);
  if (verb)
    printf("Output bdry (Mo=%d) has bo = %f\n", oM, \
	   (2*PI*oM)/(k_base*bil.Perim));
  per = dmatrix(1, ne, 1, oM);
  ngr = dmatrix(1, ne, 1, oM);
  ten = dvector(1, ne);
  nrm = dvector(1, ne);

  switch (task) { // ........................................................

  case TASK_BASISTEST: // TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    // make standard basis vectors, and set of same k's... ('calc routine')
    for (j=1; j<=ne; ++j) {
      for (i=1; i<=bas.N; ++i)
	a[j][i] = (i==(j+j_lo-1)) ? 1.0 : 0.0;
      // a[j][i] = 2*random()/MAXINT;
      ks[j] = k_base; // + 0.1*j; // mess around with wavenumbers
      kos[j] = k_base; // + (j>ne/2);
    }
    // evaluate per, ngr values... and general-BC enclosed mass
    for (i=1; i<=ne; ++i) {
      nrm[i] = eval_bdry(k_base*ks[i]/kos[i], a[i], &bas, per[i], ngr[i], \
			 &ten[i], &bil, &obps, 0, 1, 0, dummy);
      if (verb)
	printf("\teval_bdry %d.\n", i);
    }
    break;


  case TASK_NORMTEST: // NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    // a single random wavefunction... j=1
    j = 1;
    srandom((unsigned int)seed);  // allows repeatable random coeffs
    for (i=1; i<=bas.N; ++i)
      a[j][i] = 2*(random()/MAXINT) - 1.0; // bracket order stops overflow! 
    ks[j] = kos[j] = k_base;
    // test area norm formulae using the bdry...
    nrm[j] = eval_bdry(k_base, a[j], &bas, per[j], ngr[j], &ten[j], \
		&bil, &obps, 0, 1, 0, dummy);
    if (verb)
      printf("\teval_bdry %d.\n", j);
    break;
    

  case TASK_EIGENVEC: // EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    // ks carries eigenvalue
    //diag_quad_form(&bil, &bps, &bas, k_base, ks, a, j_lo, j_hi, weight);
    diag_quad_form_FG(&bil, &bps, &bas, k_base, ks, a, j_lo, j_hi, weight,\
		      epsilon);
    // evaluate per, ngr values...
    for (i=1; i<=ne; ++i) {
      //kos[i] = k_base;
      kos[i] = i-1+j_lo; // kos carries state number from diag quad form
      if (bdryfuncs) {
	nrm[i] = eval_bdry(k_base, a[i], &bas, per[i], ngr[i], &ten[i], \
			   &bil, &obps, 0, 1, 0, dummy);
	if (verb)
	  printf("\teval_bdry %d.\n", i);
      } else
	nrm[i] = 0.0;  // dummy
    }
    break;


  case TASK_VERGINI: // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    nev = vergini(&bil, &bps, &bas, k_base, epsilon, delta_lo, delta_hi, \
		  ks, a, ne, C4);
    ne = min(ne, nev); // if too many found, can't use them
    if (verb)
      printf("vergini found %d states in k=[%f,%f]; (%d after clipping)\n", \
	     nev, k_base-delta_lo, k_base+delta_hi, ne);
    // evaluate per, ngr values... nrm is pre-normalization Dirichlet norm.
    for (i=1; i<=ne; ++i) {
      nrm[i] = eval_bdry(ks[i], a[i], &bas, per[i], ngr[i], &ten[i], \
			 &bil, &obps, 1, 0, 0, dummy);
      kos[i] = k_base; // all were done at same base wavenumber...
      if (verb)
	printf("\teval_bdry %d (k=%g, ten=%g, nrm=%g)\n", i, ks[i], \
	       ten[i], nrm[i]);
    }
    if (remove_spurious) {
      // remove spurious states (may reduce ne)...
      ne = spurious(ne, ks, a, per, ngr, nrm, ten, oM, bas.N, idx);
      if (verb)
	printf("kept %d non-spurious states\n", ne);
    }
    for (i=1; i<=ne; ++i)
      kos[i] = k_base; // all were done at same base wavenumber...
    break;


  case TASK_REFVERG: // RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    nev = vergini(&bil, &bps, &bas, k_base, epsilon, delta_lo, delta_hi, \
		  ks, a, ne, C4);
    if (verb)
      printf("vergini found %d states in k=[%f,%f]\n", nev, k_base-delta_lo,\
	     k_base+delta_hi);
    ne = min(ne, nev); // if too many found, can't use them
    // evaluate per, ngr values... nrm is pre-normalization Dirichlet norm.
    for (i=1; i<=ne; ++i) {
      nrm[i] = eval_bdry(ks[i], a[i], &bas, per[i], ngr[i], &ten[i], \
			 &bil, &obps, 1, 0, 0, dummy);
      kos[i] = k_base; // all were done at same base wavenumber...
      if (verb)
	printf("\teval_bdry %d (k=%g, ten=%g, nrm=%g)\n", i, ks[i], \
	       ten[i], nrm[i]);
    }
    if (remove_spurious) {
      // remove spurious states (may reduce ne)...
      ne = spurious(ne, ks, a, per, ngr, nrm, ten, oM, bas.N, idx);
      if (verb)
	printf("kept %d non-spurious states\n", ne);
    }

    // ....... reference part: improve each state with separate vergini call.
    ne_ref = 10; // array size for only a couple states (<ne)
    ks_ref = dvector(1, ne_ref); // reference k array
    a_ref = dmatrix(1, ne_ref, 1, bas.N); // coeffs
    for (i=1; i<=ne; ++i) { // redo fresh vergini for each state
      k0 = ks[i] - delta_ref; // choose k0
      nev = vergini(&bil, &bps, &bas, k0, epsilon, 0.0, 2*delta_ref, \
		  ks_ref, a_ref, ne_ref, C4);
      if (verb)
	printf("ref vergini: found %d (should be 1) states near k=%g\n",\
		nev, ks[i]);
      if (nev>0) { // only replace state if found at least one potential ref...
	// find k nearest original ks[i], call this the reference state...
	gap = 1e6; // larger than delta_ref
	j_ref = -1;
	for (j=1; j<=nev; ++j)
	  if (fabs(ks_ref[j] - ks[i]) < gap) {
	    gap = fabs(ks_ref[j] - ks[i]);
	    j_ref = j;
	  }
	if (j_ref==-1) {
	  fprintf(stderr, "\t\t\tref vergini: problem, no nearest state!\n");
	} else {
	  ks[i] = ks_ref[j_ref];
	  nrm[i] = eval_bdry(ks[i], a_ref[j_ref], &bas, per[i], \
			     ngr[i], &ten[i], &bil, &obps, 1, 0, 0, dummy);
	  kos[i] = k0;
	  if (verb)
	    printf("\tref eval_bdry %d (k=%g, ten=%g, nrm=%g)\n", i, ks[i], \
		   ten[i], nrm[i]);
	  // copy coeffs back into original array...
	  for (j=1; j<=bas.N; ++j)
	    a[i][j] = a_ref[j_ref][j];
	}
      }
    }
    break;


  case TASK_INHOMOG: // IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    /* This task solves for coeffs a which best fit inhomogeneous boundary
       values (generated by a 2nd basis, bas_inh, currently with 1 func).
       1/2/04 barnett
    */
    make_inhomog_basis(k_base, &bil, &bas_inh, f_inh, kD_inh); // set bas_inh
    a_inh = dmatrix(1, 1, 1, 1); // coeff for inhomog src.
    a_inh[1][1] = 1.0; // choose as unity
    // dummy k values...
    ks[1] = ks[2] = kos[1] = kos[2] = k_base;

    // solve for coeffs a given inh wavefunction
    inhomog(&bil, &bps, &bas, k_base, epsilon, a, &bas_inh, a_inh);

    // fill per and ngr arrays.
    i = 1;
    nrm[i] = eval_bdry(ks[i], a[i], &bas, per[i], ngr[i], &ten[i], \
			 &bil, &obps, 0, 1, 0, dummy);
    i = 2;
    nrm[i] = eval_bdry(ks[i], a_inh[1], &bas_inh, per[i], ngr[i], &ten[i], \
			 &bil, &obps, 0, 1, 0, dummy);
    break;


  case TASK_MASSCHOP: // CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    // build line billiard and bdry pt set, depending on bil shape:
    bil_chop.type = LINE;
    if (bil.type==QU_GEN_RECT_SINAI) { // qugrs: uses params f_c and y_c
      // round coord on perim to lie at integer point location on output bps...
      m_c = max(min(int(f_c*obps.M + 0.5), obps.M), 0);
      desc_by_perim(m_c/(double)obps.M, &x_fc, &y_fc, &dummy, &dummy, &bil);
      // make line cut from (0, y_c) to perim location at (x_fc, y_fc)...
      bil_chop.param[0] = 0.0;
      bil_chop.param[1] = y_c;
      bil_chop.param[2] = x_fc;
      bil_chop.param[3] = y_fc;
    } else if (bil.type==QU_STADIUM) { // qust: params f_c and y_c dummy
      alpha = bil.param[0];
      // round coord on perim to lie at integer point location on output bps...
      f_c = alpha / (alpha + PI); // wings start at kink-point
      m_c = max(min(int(f_c*obps.M + 0.5), obps.M), 0);
      // vertical chop line...
      bil_chop.param[0] = alpha/2;
      bil_chop.param[1] = 0;
      bil_chop.param[2] = alpha/2;
      bil_chop.param[3] = 1;
    } else {
      printf("mass_chop only works for qugrs and qust billiard types!\n");
      exit(1);
    }
    if (verb)
      printf("chop intersects at bdry pt index m_c=%d\n", m_c);
    build_billiard(&bil_chop, k_base);
    build_bdry(bo, k_base, &bil_chop, &bps_chop); // closest pt density, use bo
    if (geom) {
      sprintf(headchop, "%s.chop", head);
      dump_billiard(&bil_chop, headchop);
      dump_colloc(&bps_chop, headchop);
    }
    if (verb)
      printf("mass chopping line built (Mc=%d).\n", bps_chop.M);

    // calc eigenstates of original shape...
    nev = vergini(&bil, &bps, &bas, k_base, epsilon, delta_lo, delta_hi, \
		  ks, a, ne, C4);
    ne = min(ne, nev); // if too many found, can't use them
    if (verb)
      printf("vergini found %d states in k=[%f,%f]; (%d after clipping)\n", \
	     nev, k_base-delta_lo, k_base+delta_hi, ne);

    // evaluate per, ngr... norming & collecting only mass beyond split m_c
    mass = dvector(1, ne);
    for (i=1; i<=ne; ++i) {
      nrm[i] = eval_bdry(ks[i], a[i], &bas, per[i], ngr[i], &ten[i], \
			 &bil, &obps, 1, 0, m_c, mass[i]); // Dirichlet BCs
      kos[i] = k_base; // all were done at same base wavenumber...
      if (verb)
	printf("\teval_bdry %d (k=%g, ten=%g, nrm=%g)\n", i, ks[i], \
	       ten[i], nrm[i]);
    }

    if (remove_spurious) {
      // remove spurious states (may reduce ne)...
      ne = spurious(ne, ks, a, per, ngr, nrm, ten, oM, bas.N, idx);
      // juggle mass_chop data accordingly...
      for (i=1; i<=ne; ++i)
	mass[i] = mass[idx[i]];
      if (verb)
	printf("kept %d non-spurious states\n", ne);
    }

    // add on mass contrib from chop line, discarding per, ngr to dummy arrays
    per_chop = dmatrix(1, ne, 1, bps_chop.M);
    ngr_chop = dmatrix(1, ne, 1, bps_chop.M);
    for (i=1; i<=ne; ++i) {
      mass[i] += eval_bdry(ks[i], a[i], &bas, per_chop[i], ngr_chop[i], \
			   &dummy, &bil_chop, &bps_chop, 0, 1+sec_deriv_mc, \
			   0, dummy); // general BCs
      if (verb)
	printf("\teval_chopped_mass %d (mass=%g)\n", i, mass[i]);
    }
    break;


  } // ..................................................................


  if (task==TASK_MASSCHOP)
    save_sum(cmdline, ks, mass, ten, nrm, ne, head); // mass replaces kos
  else
    save_sum(cmdline, ks, kos, ten, nrm, ne, head);
  if (want_coeffs) {
    // output .cf...
    save_oned(&bil, a, kos, ne, bas.N, 1, 1, head, "cf");
  }
  if (bdryfuncs)
    if (task==TASK_MASSCHOP) {
      // only do original bdry part within chop region (m_c < m =< M)...
      // for (i=1; i<=ne; ++i)
      // ngr[i] += m_c; // advance pointers to 1d ngr arrays by m_c
      // save_oned(&bil, ngr, ks, ne, oM-m_c, bil.Perim/oM, 0.5, head, "ngr");
      // chop line data...
      save_oned(&bil_chop, per_chop, ks, ne, bps_chop.M, \
		bil_chop.Perim/bps_chop.M, 0.5, head, "chop.per");
      save_oned(&bil_chop, ngr_chop, ks, ne, bps_chop.M, \
		bil_chop.Perim/bps_chop.M, 0.5, head, "chop.ngr");
    } else {
      save_oned(&bil, per, ks, ne, oM, bil.Perim/oM, 0.5, head, "per");
      save_oned(&bil, ngr, ks, ne, oM, bil.Perim/oM, 0.5, head, "ngr");
      if (geom)
	// write out r.n array at the output bdry pt set... (ks[1] dummy)
	save_oned(&bil, &obps.rn-1, ks, 1, oM, bil.Perim/oM, 0.5, head, "rn");
    }
  if (grid) {
    // 2d spatial grid output...
    make_grid(&bil, &g, ne, dx, mask, xl, xh, yl, yh);
    if (task==TASK_INHOMOG) { // output desired inhomog func in grid n=2
      eval_grid_vecs(&bil, &bas_inh, k_base, a_inh, &g, 1, 1);
      ne = 1; // so that grid calcs only now happen for the 1 state
    }
    if (geom)
      save_mask(&bil, &g, head, xl, xh, yl, yh);
    if (grid==1) // fast, same k
      eval_grid_vecs(&bil, &bas, k_base, a, &g, ne, 0);
    else // slow, correct ks
      eval_grid_vecs_each_k(&bil, &bas, ks, kos, a, &g, ne, 0);
    if (encode)
      encode_grid(&bps, &g, ks, k_base, encode_res);
    save_grid(&g, ks, head);
    free_grid(&g);
  }
  
  return 0;
}
