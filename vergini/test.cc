/* test.cc: Test the modules in VERGINI package, writing out some basis sets
 *
 * barnett 9/10/03
 * 11/16/03 changed cmd line parsing to getopt
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

/* define global verbosity... set its default */
int verb = 1;


void show_usage()
  /* For test.cc only
   */
{
  fprintf(stderr, "\nTest Vergini modules by dumping basis funcs\n");
  fprintf(stderr, "Usage:\n\t.\\test [options] j_lo j_hi\nOptions:\n");
  fprintf(stderr, "\t-v\tverbose\n\t-q\tquiet (no stdout)\n");
  fprintf(stderr, "\t-b <B>\tset boundary collocation to B per wavelength\n");
  fprintf(stderr, "\t-k <K>\tset wavenumber to K\n");
  fprintf(stderr, "\t-l type:a:b:...:z\tdefine billiard shape (required)\n");
  fprintf(stderr, "\t-s type:a:b:...:z\tdefine basis set (required)\n");
  fprintf(stderr, "\t-g <dx>\toutput 2d grids, spacing dx\n\t-h\tshow help\n");
  fprintf(stderr, "\t-p <B>\toutput bdry funcs, B samples per wavelength\n");
  fprintf(stderr, "\t-o <head>\theader name for output files\n");
  fprintf(stderr, "\t-m\tdump geometry files\n");
  fprintf(stderr, "\t-e <s>\tencode grid output with colloc points\n");
  fprintf(stderr, "\t\ts = 0: don't rescale with k\n\t\ts = 1: do rescale\n");
  fprintf(stderr, "\nExample:\n\t./test -v -l qugrs:1:0.2:0.4 -b 10 -s oyooo:3:3:1 -k 50 -e 1 -p 5 -m 1 10\n\n");
}


/* ======================================= MAIN ======================== */
int main(int argc, char **argv)
{
  Billiard bil;
  Bdry_Pt_Set bps, obps;
  Basis_Set bas;
  Grid g;
  int geom = 0, grid = 0, encode = 0, encode_res = 0, bdryfuncs = 0;
  int i, j, j_lo = 1, j_hi = 5, ne, oM;
  double b = 10, bfb = 10, k = 20, dx = 0.01;
  double **a, *ks, *kos, **per, **ngr, *ten;
  char cmdline[LEN], *head = "t";
  
  // parse command line....................................................
  grab_cmdline(cmdline, argc, argv);
  opterr = 0;
  bil.type = NBILLIARDS;
  bas.set_type = NBASES;
  while ((i = getopt (argc, argv, "l:s:b:k:o:vhmqg:e:p:")) != -1)
    switch (i) {
    case 'o': // output file head
      head = optarg;
      break;
    case 'g': // output 2d grid with dx
      grid = 1;
      sscanf(optarg, "%lf", &dx);
      break;
    case 'p': // output bdry funcs with bfb = "boundary func b"
      bdryfuncs = 1;
      sscanf(optarg, "%lf", &bfb);
      break;
    case 'b':
      sscanf(optarg, "%lf", &b);
      break;
    case 'k':
      sscanf(optarg, "%lf", &k);
      break;
    case 'l':
      if (parse_billiard(optarg, &bil)==-1) {
	show_usage();
	show_billiard_usage();
	return 1;
      }
      break;
    case 's':
      if (parse_basis(optarg, &bas)==-1) {
	show_usage();
	show_basis_usage();
	return 1;
      }
      break;
    case 'q': // quiet
      verb = 0;
      break;
    case 'v': // very verbose
      verb = 255;
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
  // read non-option arguments...
  if (optind==argc-2) {
    sscanf(argv[optind], "%d", &j_lo);
    sscanf(argv[optind+1], "%d", &j_hi);
    if (verb)
      printf("j_lo,j_hi = %d,%d\n", j_lo, j_hi);
  } else {
    fprintf(stderr, "incorrect number of arguments!\n");
    return 1;
  }
  // .............................................................

  build_billiard(&bil);
  build_bdry(b, k, &bil, &bps);
  build_basis(k, &bil, &bas);
  if (verb) {
    printf("\ntest: test modules in VERGINI package. 11/16/03\n\n");
    if (verb&0x80)
      printf("cmdline: %s\n\n", cmdline);
    show_billiard_properties(&bil, k); 
    show_colloc_properties(&bil, &bps, k); 
    show_basis_properties(&bil, &bas);
    printf("Output file head: %s\n\n", head);
  }
  if (geom) {
  // output geom...
    if (verb)
      printf("dumping geometry files...\n");
    dump_billiard(&bil, head);
    dump_colloc(&bps, head);
    dump_basis(&bas, head);
  }

  ne = j_hi - j_lo + 1;
  ks = dvector(1, ne); // wavenumbers k
  kos = dvector(1, ne); // k_0 used for each state
  a = dmatrix(1, ne, 1, bas.N); // coeffs

  // set up bdry funcs on an 'obps' (output bdry pt set)...
  build_bdry(bfb, k, &bil, &obps);
  oM = obps.M;
  per = dmatrix(1, ne, 1, oM);
  ngr = dmatrix(1, ne, 1, oM);
  ten = dvector(1, ne);

  // make standard basis vectors, and set of same k's... ('calc routine')
  for (j=1; j<=ne; ++j) {
    for (i=1; i<=bas.N; ++i)
      a[j][i] = (i==(j+j_lo-1)) ? 1.0 : 0.0;
    ks[j] = k;
    kos[j] = k;
  }

  // evaluate per, ngr values...
  for (i=1; i<=ne; ++i) {
    eval_bdry(ks[i], a[i], &bas, per[i], ngr[i], &ten[i], &bil, &obps, 0);
    if (verb)
      printf("eval_bdry %d.\n", i);
  }

  // output .cf...
  save_oned(a, kos, ne, bas.N, 1, 1, head, "cf");
  if (grid) {
    // 2d spatial grid output...
    make_grid(&bil, &g, ne, dx);
    eval_grid_vecs(&bil, &bas, k, a, &g, 0);
    if (encode)
      encode_grid(&bps, &g, ks, k, encode_res);
    save_grid(&g, ks, head);
    free_grid(&g);
  }
  if (bdryfuncs) {
    save_oned(per, ks, ne, oM, bil.Perim/oM, 0.5, head, "per");
    save_oned(ngr, ks, ne, oM, bil.Perim/oM, 0.5, head, "ngr");
  }
  
  return 0;
}
