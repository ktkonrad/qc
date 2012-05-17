/*
Driver for count_nodal_domains.c using percolation model

Kyle Konrad
3/31/2011
*/


#include "count_nodal_domains.h"
#include "random_percolation.h"
#include "util/util.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h> // for command line parsing with getopt

// global verbosity variable
int verb = 1;

// options specified by command line arguments
int showTime = 0; // flag: print time it takes for countNodalDomains to run
int gridSize = -1; // size of grid in x and y dimensions
int trials = 1; // number of grids to generate and count
int trial = 1; // current trial number
int outputGrid = 0; // flag: output grid(s) to file(s): grid_[n].dat
double k1 = -1; // effective k of grid: k = (gridSize*pi)/(2*sqrt(2))
double k2 = 0; // max of k range
double kstep = 0; // step of k range
FILE *sizefile = NULL; // file to write sizes to
int mode = 0; // 0: single k. 1: k range

/*
  print a usage statement
*/
void usage() {
  fprintf(stderr, "USAGE: count {-n gridSize | -k k | -r k1:kstep:k2} [-N trials] [-t] [-o] [-s]\n");
  fprintf(stderr, "-t: show timing info\n");
  fprintf(stderr, "-o: output grid to file\n");
  fprintf(stderr, "-s: output nodal domain sizes to file\n");
}  

/*
  parse a string "k:kstep:k2"
*/
int parse_k_range(char *str) {
  char *first_colon = strchr(str, ':');
  char *last_colon = strrchr(str, ':');
  if (!first_colon || !last_colon) {
    return -1;
  }

  // replace colons with terminating null characters
  *first_colon = 0;
  *last_colon = 0;

  // parse the pieces
  k1 = atof(str); // str now ends where the first colon was
  kstep = atof(first_colon+1);
  k2 = atof(last_colon+1);

  if (k1 <= 0 || kstep <= 0 || k2 <= 0) {
    return -1;
  }

  return 0;
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

  while ((c = getopt(argc, argv, "n:N:k:r:tos")) != -1) {
    switch (c) {
    case 'n':
      gridSize = atoi(optarg);
      k1 = gridSize * M_PI * sqrt(2) / 2;
      k2 = k1;
      kstep = 1.0;
      break;
    case 'N':
      trials = atoi(optarg);
      break;
    case 'k':
      k1 = atof(optarg);
      k2 = k1;
      kstep = 1.0;
      gridSize = k1/(sqrt(2)*M_PI)*2;
      break;
    case 'r':
      mode = 1;
      if (parse_k_range(optarg)) {
        ERROR("failed to parse k range: %s", optarg);
        exit(CMD_LINE_ARG_ERR);
      }
      break;
    case 't':
      showTime = 1;
      break;
    case 'o':
      outputGrid = 1;
      break;
    case 's':
      sizefile = fopen("perc_sizes.dat", "w");
      break;
    case '?':
      switch (optopt) {
      case 'n':
      case 'N':
      case 'k':
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
	break;
      default:
	fprintf(stderr, "Unknown option -%c.\n", optopt);
      }
    default:
      abort();
    }
  }

  if (k1 <= 0) {
    ERROR("invalid k: %f", k1);
    exit(CMD_LINE_ARG_ERR);
  }
}


/*
output:
        return value: number of nodal domains
*/
int runTest(double **grid, int ny, int nx) {   
  clock_t start = clock();
  int nd = countNodalDomainsNoInterp(grid, NULL, ny, nx, sizefile);
  clock_t end = clock();

  if (showTime) {
    printf("countNodalDomains took %f seconds\n", ((double)(end - start)) / CLOCKS_PER_SEC);
  }
  return nd;
}


int main(int argc, char **argv) {
  processArgs(argc, argv);
  int ny, nx;
  double **grid;
  double k;

  for (k = k1 ; k <= k2 + 1e-15 ; k += kstep) { // + 1e-15 to deal with potential floating point issues
    gridSize = (int)(k/(sqrt(2)*M_PI))*2; // gridSize = 2m (where m is number of sites in y-direction as in Bogomolny 2001)

    nx = gridSize;
    ny = gridSize;
    int counts[trials];

    double mean = 0.0, variance = 0.0;

    int i;
    grid = dmatrix(ny, nx);
    if (!grid) {
      exit(OUT_OF_MEMORY_ERR);
    }
    srand(time(NULL));

    // run the trials and calculate mean count
    for (trial = 0 ; trial < trials ; trial++) {
      randomPercolation(grid, ny, nx);
    
      if (outputGrid) {
        char outfile[50];
        sprintf(outfile, "../data/grid_%d.dat", trial);
	array2file(grid, ny, nx, outfile);
      }
    
      counts[trial] = runTest(grid, ny, nx);
      switch (mode) {
      case 0:
        printf("%d\n", counts[trial]);
        break;
      case 1:
        printf("%f,%d\n", k, counts[trial]);
        break;
      }
      mean += counts[trial];
    }
    free_dmatrix(grid);
    switch (mode) {
    case 0:
      mean /= trials;
      // calculate variance of counts
      for (i = 0 ; i < trials ; i++) {
        variance += (counts[i] - mean) * (counts[i] - mean);
      }
      variance /= trials;
      // printf("mean: %f\n", mean);
      // printf("variance: %f\n", variance);
      printf("scaled mean: %f\n", mean / (ny * nx / 4) * 2 / M_PI);
      printf("scaled variance: %f\n", variance / (ny * nx / 4) * 2 / M_PI);
      break;
    case 1:
      break;
    }
  } 

  if (sizefile) {
    fclose(sizefile);
  }
  exit(0);
}
