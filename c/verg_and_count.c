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


// options specified by command line arguments
char *name; // base name of verg output
char *billiard; // billiard shape string
char *basis; // basis set string
char *vc_dx; // grid spacing
char *k; // k_0 value from vergini
char *window; // window size on either size of k
char *fourth_order_coeff; // fourth order vergini coefficient
char *bessel_order; // highest order bessel function to use for interpolation
char *vc_upsample; // upsampling ratio to use for interpolation
char *remove_spurious; // -u flag to verg

/*
  print a usage statement
*/
void vc_usage() {
  fprintf(stderr, "VC_USAGE: verg_and_count -n name -l billiardType -s basisSet -d vc_dx -k k -V verginiWidth -M besselOrder -p vc_upsample\n");
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

  while ((c = getopt(argc, argv, "n:l:s:d:k:V:4:M:p:u")) != -1) {
    switch (c) {
    case 'n':
      name = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(name, optarg);
      break;
    case 'l':
      billiard = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(billiard, optarg);
      break;
    case 's':
      basis = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(basis, optarg);
      break;
    case 'M':
      bessel_order = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(bessel_order, optarg);
      break;
    case 'p':
      vc_upsample = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(vc_upsample, optarg);
      break;
    case 'V':
      window = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(window, optarg);
      break;
    case '4':
      fourth_order_coeff = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(fourth_order_coeff, optarg);
      break;
    case 'd':
      vc_dx = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(vc_dx, optarg);
      break;
    case 'k':
      k = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(k, optarg);
      break;
    case 'u':
      remove_spurious = (char *)malloc(3*sizeof(char));
      strcpy(remove_spurious, "-u");
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

  if (vc_dx == NULL) {
    ERROR("vc_dx not specified or invalid");
    exit(CMD_LINE_ARG_ERR);
  }
  if (k == NULL) {
    ERROR("k not specified or invalid");
    exit(CMD_LINE_ARG_ERR);
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

#define SET(loc, val) do {loc = (char *)malloc((strlen(val)+1)*sizeof(char)); strcpy(loc, val);} while (0)

#define COUNT_NARGS 13
#define VERG_NARGS 16

int main(int argc, char **argv) {
  int pid;
  int i;
  int rc = 0;

  processArgs(argc, argv);

  //./verg -o test -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1.000000 -k 200.100000 -V 0.1:0.12 -f 0.001000

  char **verg_args = (char **)malloc(VERG_NARGS*sizeof(char *));
  char *verg_executable = "verg";
  verg_args[0] = (char *)malloc(strlen(verg_executable)*sizeof(char));
  strcpy(verg_args[0], verg_executable);
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
  

  verg_main(16, verg_args);

  /*
  for (i = 0 ; i < VERG_NARGS-1 ; i++) {
    free(verg_args[i]);
  }  
  free(verg_args);
  */

  //   ./count -f test.sta_bin -l qugrs:1.0:0.4:0.7 -d 0.001000 -k 200.100000 -M 9 -u 20
  char **count_args = (char **)malloc(COUNT_NARGS * sizeof(char *));
  char *count_executable = "count";
  count_args[0] = (char *)malloc(strlen(count_executable)*sizeof(char));
  strcpy(count_args[0], count_executable);
  for (i = 1 ; i < COUNT_NARGS - 1 ; i+=2) {
    count_args[i] = (char *)malloc(3*sizeof(char));
  }
  strcpy(count_args[1], "-f");
  count_args[2] = (char *)malloc((strlen(name) + 9)*sizeof(char));
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
  count_args[13] = NULL;

  count_main(13, count_args);

  /*
  for (i = 0 ; i < COUNT_NARGS-1 ; i++) {
    free(count_args[i]);
  }  
  free(count_args);
  */

}
