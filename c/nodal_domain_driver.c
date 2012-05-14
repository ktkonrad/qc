/*
Driver for count_nodal_domains.c

Kyle Konrad
3/31/2011
*/


#include "count_nodal_domains.h"
#include "random_percolation.h"
#include "util/util.h"
#include "util/util_verg.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h> // for command line parsing with getopt

// vergini code dependencies
#include "../vergini/billiard.h"
int verb = 1;

// options specified by command line arguments
int showTime = 0; // flag: print time it takes for countNodalDomains to run
int outputGrid = 0; // flag: output grid(s) to file(s): grid_[n].dat
char *file; // use grid from this file
int mode = -1; // 0 - generate random_percolation grid; 1 - read grid from sta_bin file
int maskFlag = 0; // flag: use a mask file
char *maskFile; // use mask from this file
int oneFlag = 0; // flag: only count one eigenfunction
int sizeFlag = 0; // flag: write domain sizes to file

Billiard bil; // Billiard shape we are using - defined in vergini code
double dx = -1; // grid spacing
double k_0 = -1; // k_0 value from vergini
double alpha = -1; // k_0 * dx

int besselOrder = -1; // highest order bessel function to use for interpolation
int upsample = -1; // upsampling ratio to use for interpolation
int interp = 1; // boolean whether or not to interpolate

/*
  print a usage statement
*/
void usage() {
  fprintf(stderr, "USAGE: count -f file [-m maskFile | {-l billiardType -k k_0 {-d dx | -a alpha } -M besselOrder -u upsample}] [-t] [-o] [-1] [-n] [-q] [-s]\n");
  fprintf(stderr, "-t: show timing info\n");
  fprintf(stderr, "-o: output grid to file\n");
  fprintf(stderr, "-1: only count first eigenfunction\n");
  fprintf(stderr, "-n: do not use interpolation\n");
  fprintf(stderr, "-q: quiet\n");
  fprintf(stderr, "-s: output nodal domain sizes to file\n");
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

  while ((c = getopt(argc, argv, "f:m:l:d:a:k:M:u:to1nqs")) != -1) {
    switch (c) {
    case 'f':
      SET(file, optarg);
      mode = 1;
      break;
    case 'm':
      SET(maskFile, optarg);
      maskFlag = 1;
      break;
    case 'l':
      if (parse_billiard(optarg, &bil) == -1) {
	fprintf(stderr, "Error: failed to parse billiard args\n");
	usage();
	exit(CMD_LINE_ARG_ERR);
      }
      mode = 2;
      break;
    case 'M':
      besselOrder = atoi(optarg);
      break;
    case 'u':
      upsample = atoi(optarg);
      break;
    case 'd':
      dx = (double)atof(optarg);
      break;
    case 'a':
      alpha = (double)atof(optarg);
      break;
    case 'k':
      k_0 = (double)atof(optarg);
      break;
    case 't':
      showTime = 1;
      break;
    case 'o':
      outputGrid = 1;
      break;
    case '1':
      oneFlag = 1;
      break;
    case 'n':
      interp = 0;
      break;
    case 'q':
      verb = 0;
      break;
    case 's':
      sizeFlag = 1;
      break;
    case '?':
      switch (optopt) {
      case 'f':
      case 'm':
      case 'l':
      case 'd':
      case 'k':
      case 'M':
      case 'u':
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
	break;
      default:
	fprintf(stderr, "Unknown option -%c.\n", optopt);
      }
    default:
      abort();
    }
  }

  if (file == NULL || strcmp(file, "") == 0) {
    ERROR("no sta_bin file was given");
    exit(CMD_LINE_ARG_ERR);
  }

  if (mode == 2) {
    if (k_0 <= 0) {
      ERROR("k_0 not specified or invalid");
      usage();
      exit(CMD_LINE_ARG_ERR);
    }
    if (dx <= 0) {
      if (alpha <= 0) {
        ERROR("dx and alpha not specified or invalid");
        usage();
        exit(CMD_LINE_ARG_ERR);
      }
      dx = alpha / k_0;
    }
    if (besselOrder <= 0) {
      ERROR("besselOrder not specified or invalid");
      usage();
      exit(CMD_LINE_ARG_ERR);
    }
    if (upsample <= 0) {
      ERROR("upsample not specified or invalid");
      usage();
      exit(CMD_LINE_ARG_ERR);
    }
  }
}


/*
output:
        return value: number of nodal domains
*/
int runTest(double **grid, char **mask, int ny, int nx, double k, double dx, int besselOrder, int upsample, interp_stats *stats) {
  char sizefilename[50];
  FILE *sizefile = NULL;
  if (sizeFlag) {
    sprintf(sizefilename, "sizes_%f_%f.dat", k, dx);
    sizefile = fopen(sizefilename, "w");
  }
  int nd;
  if (interp) {
    nd = countNodalDomainsInterp(grid, mask, ny, nx, k, dx, besselOrder, upsample, stats, sizefile);
  }
  else {
    nd = countNodalDomainsNoInterp(grid, mask, ny, nx, sizefile);
  }
  if (sizeFlag) {
    fclose(sizefile);
  }
  return nd;
}



int main(int argc, char **argv) {
  processArgs(argc, argv);
  int ny, nx;
  double **grid;
  interp_stats stats;

  clock_t start = clock();

  if (mode == 1) {
    int count;
    int masky, maskx;
    char **mask = NULL;
    grid = readOneSta(file, &ny, &nx);
    
    if (grid == NULL) {
      ERROR("failed to read grid");
      exit(UTIL_ERR);
    }

    if (maskFlag) {
      mask = readMask(maskFile, &masky, &maskx);
     
      if (maskx != nx || masky != ny) {
	ERROR("mask dimensions do not match grid dimensions");
	exit(DIMENSION_ERR);
      }
    }
    bzero(&stats, sizeof(stats));
    count = runTest(grid, mask, ny, nx, k_0, dx, besselOrder, upsample, &stats);

    destroyGrid(grid);
    free(file);

    if (maskFlag) {
      destroyMask(mask);
      free(maskFile);
    }

    printf("%f,%f,%d,%d,%d,%d,%d\n", k_0, dx, count, stats.small_domain_count, stats.interp_count, stats.boundary_trouble_count, stats.edge_trouble_count);

  }

  if (mode == 2) {
    int rc;
    int count;
    int k_base = 20; // to be passed to build_billiard
    int masky, maskx;
    char **mask = NULL;
    double k, wtm;
    int ne;


    rc = build_billiard(&bil, k_base);
    if (rc != 0) {
      ERROR("failed to build billiard");
      exit(VERGINI_ERR);
    }

    //    printf("%s\t%s\t%s\t%s\t%s\n", "k", "count", "small domains", "interp count", "boundary trouble count", "edge trouble count");

    int i = 0;
    do {
      bzero(&stats, sizeof(stats));
      grid = readSta(file, &ne, &ny, &nx, &k, i); // read eigenfunctions one at atime so we don't have to keep them all in memory at once

      if (ne == 0) {
        break;
      }

      if (grid == NULL) {
	ERROR("failed to read grid");
	exit(IO_ERR);
      }
    
      mask = createScaledMaskFromBilliard(bil, dx, &masky, &maskx, k/k_0);

      if (maskx != nx || masky != ny) {
	ERROR("mask dimensions do not match grid dimensions\n \
	       ny\tnx\tmasky\tmaskx\n \
	       %d\t%d\t%d\t%d",
	      ny,nx,masky,maskx);
	exit(DIMENSION_ERR);
      }

      count = runTest(grid, mask, ny, nx, k, dx, besselOrder, upsample, &stats);
      if (bil.type == QU_STADIUM) {
        wtm = wingTipMass(grid, mask, ny, nx);
      }
      else {
        wtm = 0;
      }

      destroyMask(mask);
      destroyGrid(grid);

      printf("%f, %d, %d, %d, %d, %d, %f\n", k, count, stats.small_domain_count, stats.interp_count, stats.boundary_trouble_count, stats.edge_trouble_count, wtm);

      if (oneFlag)
	break;

    } while (++i < ne);

      free(file);
  }

  if (showTime)
    printf("counting took %f seconds\n", ((double)(clock() - start)) / CLOCKS_PER_SEC);

  exit(0);
}

