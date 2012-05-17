/*
  Kyle Konrad
  2/20/2011

  driver program from vergini and count
*/


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "util/util.h"
#include "util/exit_codes.h"
#include "nodal_domain_driver_no_main.h"
#include "../vergini/verg_no_main.h"

#define SET(loc, val) do {loc = (char *)malloc((strlen(val)+1)*sizeof(char)); MALLOC_CHECK(loc); strcpy(loc, val);} while (0)
#define RESET(loc, val) do {loc = (char *)realloc(loc, (strlen(val)+1)*sizeof(char)); MALLOC_CHECK(loc); strcpy(loc, val);} while (0)

#define COUNT_NARGS 14
#define VERG_NARGS 16

// global verbosity
int verb = 1;

// options specified by command line arguments
char *name; // base name of verg output
char *billiard; // billiard shape string
char *basis; // basis set string
char *vc_dx; // grid spacing
char *k; // k_0 value from vergini
double vc_alpha = -1; // k * vc_dx (not passed as arg, just used to compute dx from k)
char *window; // window size on either size of k
char *fourth_order_coeff; // fourth order vergini coefficient
char *bessel_order; // highest order bessel function to use for interpolation
char *vc_upsample; // upsampling ratio to use for interpolation
char *remove_spurious; // -u flag to verg

/*
  print a usage statement
*/
void vc_usage() {
  fprintf(stderr, "VC_USAGE: verg_and_count -n name -l billiardType -s basisSet -d vc_dx -k k -V verginiWidth -M besselOrder -p vc_upsample\n [-q]");
}

/*
  free the memory used for command line args
*/
void freeArgs() {
  free(name);
  free(billiard);
  free(basis);
  free(bessel_order);
  free(vc_upsample);
  free(window);
  free(fourth_order_coeff);
  free(vc_dx);
  free(k);
  free(remove_spurious);
}

/*
  process command line arguments
*/
void processArgs(int argc, char **argv) {
  int i = 0;
  int c;
  opterr = 0;

  if (argc <= 1) {
    vc_usage();
    exit(CMD_LINE_ARG_ERR);
  }

  // set default values
  SET(fourth_order_coeff, "4");
  SET(remove_spurious, "");

  while ((c = getopt(argc, argv, "n:l:s:d:a:k:V:4:M:p:uq")) != -1) {
    switch (c) {
    case 'n':
      SET(name, optarg);
      break;
    case 'l':
      SET(billiard, optarg);
      break;
    case 's':
      SET(basis, optarg);
      break;
    case 'M':
      SET(bessel_order, optarg);
      break;
    case 'p':
      SET(vc_upsample, optarg);
      break;
    case 'V':
      SET(window, optarg);
      break;
    case '4':
      RESET(fourth_order_coeff, optarg);
      break;
    case 'd':
      SET(vc_dx, optarg);
      break;
    case 'a':
      vc_alpha = (double)atof(optarg);
      break;
    case 'k':
      SET(k, optarg);
      break;
    case 'u':
      RESET(remove_spurious, "-u");
      break;
    case 'q':
      verb = 0;
      break;
    case '?':
      switch (optopt) {
      case 'n':
      case 'l':
      case 's':
      case 'd':
      case 'k':
      case 'M':
      case 'p':
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
	break;
      default:
	fprintf(stderr, "Unknown option -%c.\n", optopt);
      }
    default:
      abort();
    }
  }

  // argument validations
  if (name == NULL) {
    ERROR("no name was given");
    exit(CMD_LINE_ARG_ERR);
  }

  if (billiard == NULL) {
    ERROR("no billiard was given");
    exit(CMD_LINE_ARG_ERR);
  }

  if (basis == NULL) {
    ERROR("no basis set was given");
    exit(CMD_LINE_ARG_ERR);
  }

  if (window == NULL) {
    ERROR("no window was given");
    exit(CMD_LINE_ARG_ERR);
  }
  if (k == NULL) {
    ERROR("k not specified or invalid");
    exit(CMD_LINE_ARG_ERR);
  }
  if (vc_dx == NULL) {
    if (vc_alpha <= 0) {
      ERROR("dx and alpha not specified or invalid");
      exit(CMD_LINE_ARG_ERR);
    }
    vc_dx = (char *)malloc(19*sizeof(char)); // 19 for: 0.\d{16}\0
    sprintf(vc_dx, "%.16f", vc_alpha/atof(k));
  }
  if (bessel_order == NULL) {
    ERROR("bessel_order not specified or invalid");
    exit(CMD_LINE_ARG_ERR);
  }
  if (vc_upsample == NULL) {
    ERROR("vc_upsample not specified or invalid");
    exit(CMD_LINE_ARG_ERR);
  }
}

int main(int argc, char **argv) {
  int pid;
  int i;
  int rc = 0;

  processArgs(argc, argv);

  //./verg -o test -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1.000000 -k 200.100000 -V 0.1:0.12 -f 0.001000

  char **verg_args = (char **)malloc((VERG_NARGS+1)*sizeof(char *));
  char *verg_executable = "verg";
  SET(verg_args[0], verg_executable);
  for (i = 1 ; i < VERG_NARGS ; i+=2) {
    verg_args[i] = (char *)malloc(3*sizeof(char));
  }
  strcpy(verg_args[1], "-o");
  SET(verg_args[2], name);
  strcpy(verg_args[3], "-l");
  SET(verg_args[4], billiard);
  strcpy(verg_args[5], "-s");
  SET(verg_args[6], basis);
  strcpy(verg_args[7], "-4");
  SET(verg_args[8], fourth_order_coeff);
  strcpy(verg_args[9], "-k");
  SET(verg_args[10], k);
  strcpy(verg_args[11], "-V");
  SET(verg_args[12], window);
  strcpy(verg_args[13], "-f");
  SET(verg_args[14], vc_dx);
  strcpy(verg_args[15], remove_spurious);
  verg_args[16] = NULL;
  

  verg_main(VERG_NARGS, verg_args);

  
  for (i = 0 ; i < VERG_NARGS ; i++) {
    free(verg_args[i]);
  }  
  free(verg_args);
  

  //   ./count -f test.sta_bin -l qugrs:1.0:0.4:0.7 -d 0.001000 -k 200.100000 -M 9 -u 20 -t
  char **count_args = (char **)malloc((COUNT_NARGS+1)*sizeof(char *));
  char *count_executable = "count";
  SET(count_args[0], count_executable);
  for (i = 1 ; i < COUNT_NARGS ; i+=2) {
    count_args[i] = (char *)malloc(3*sizeof(char));
  }
  strcpy(count_args[1], "-f");
  count_args[2] = (char *)malloc((strlen(name) + 9)*sizeof(char)); // + 8 for .sta_bin
  strcpy(count_args[2], name);
  strcat(count_args[2], ".sta_bin");
  strcpy(count_args[3], "-l");
  SET(count_args[4], billiard);
  strcpy(count_args[5], "-d");
  SET(count_args[6], vc_dx);
  strcpy(count_args[7], "-k");
  SET(count_args[8], k);
  strcpy(count_args[9], "-M");
  SET(count_args[10], bessel_order);
  strcpy(count_args[11], "-u");
  SET(count_args[12], vc_upsample);
  strcpy(count_args[13], "-t");
  count_args[14] = NULL;

  count_main(COUNT_NARGS, count_args);

  
  for (i = 0 ; i < COUNT_NARGS ; i++) {
    free(count_args[i]);
  }  
  free(count_args);

  freeArgs();
}
