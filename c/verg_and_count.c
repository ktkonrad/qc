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

// options specified by command line arguments
char *name; // base name of verg output
char *billiard; // billiard shape string
char *basis; // basis set string
char *dx; // grid spacing
char *k; // k_0 value from vergini
char *window; // window size on either size of k
char *fourth_order_coeff; // fourth order vergini coefficient
char *bessel_order; // highest order bessel function to use for interpolation
char *upsample; // upsampling ratio to use for interpolation
char *remove_spurious; // -u flag to verg

/*
  print a usage statement
*/
void usage() {
  fprintf(stderr, "USAGE: verg_and_count -n name -l billiardType -s basisSet -d dx -k k -V verginiWidth -M besselOrder -p upsample\n");
}

/*
  process command line arguments
*/
void processArgs(int argc, char **argv) {
  int i = 0;
  int c;
  opterr = 0;

  if (argc <= 1) {
    usage();
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
      upsample = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(upsample, optarg);
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
      dx = (char *)malloc(strlen(optarg)*sizeof(char));
      strcpy(dx, optarg);
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

  if (dx == NULL) {
    ERROR("dx not specified or invalid");
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
  if (upsample == NULL) {
    ERROR("upsample not specified or invalid");
    exit(CMD_LINE_ARG_ERR);
  }
}

int main(int argc, char **argv) {
  int pid;
  int i;
  int rc = 0;

  processArgs(argc, argv);

  //./verg -o test -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1.000000 -k 200.100000 -V 0.1:0.12 -f 0.001000

  char *verg_args[16];
  char *verg_executable = "../vergini/verg";
  verg_args[0] = (char *)malloc(strlen(verg_executable)*sizeof(char));
  for (i = 1 ; i < 16 ; i+=2) {
    verg_args[i] = (char *)malloc(3*sizeof(char));
  }
  strcpy(verg_args[1], "-o");
  verg_args[2] = name;
  strcpy(verg_args[3], "-l");
  verg_args[4] = billiard;
  strcpy(verg_args[5], "-s");
  verg_args[6] = basis;
  strcpy(verg_args[7], "-4");
  verg_args[8] = fourth_order_coeff;
  strcpy(verg_args[9], "-k");
  verg_args[10] = k;
  strcpy(verg_args[11], "-V");
  verg_args[12] = window;
  strcpy(verg_args[13], "-f");
  verg_args[14] = dx;
  verg_args[15] = remove_spurious;
  verg_args[16] = NULL;
  

  pid = fork();
  if (!pid) { // child
    execv(verg_executable, verg_args);
    ERROR("failed to execute verg");
    exit(-1);
  } else { // parent
    wait(&rc);
    if (rc) {
      ERROR("verg failed");
    }
  }

  //   ./count -f test.sta_bin -l qugrs:1.0:0.4:0.7 -d 0.001000 -k 200.100000 -M 9 -u 20
  char *count_args[13];
  char *count_executable = "count";
  count_args[0] = (char *)malloc(strlen(count_executable)*sizeof(char));
  for (i = 1 ; i < 13 ; i+=2) {
    count_args[i] = (char *)malloc(3*sizeof(char));
  }
  strcpy(count_args[1], "-f");
  count_args[2] = (char *)malloc((strlen(name) + 8)*sizeof(char));
  strcpy(count_args[2], name);
  strcat(count_args[2], ".sta_bin");
  strcpy(count_args[3], "-l");
  count_args[4] = billiard;
  strcpy(count_args[5], "-d");
  count_args[6] = dx;
  strcpy(count_args[7], "-k");
  count_args[8] = k;
  strcpy(count_args[9], "-M");
  count_args[10] = bessel_order;
  strcpy(count_args[11], "-u");
  count_args[12] = upsample;
  count_args[13] = NULL;
  pid = fork();
  if (!pid) { // child
    execv(count_executable, count_args); 
  } else { // parent
    wait(&rc);
    if (rc) {
      ERROR("count failed");
    }
  }


}
