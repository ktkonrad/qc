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

// options specified by command line arguments
int showTime = 0; // flag: print time it takes for countNodalDomains to run
int gridSize = -1; // size of grid in x and y dimensions
int trials = 1; // number of grids to generate and count
int trial = 1; // current trial number
int outputGrid = 0; // flag: output grid(s) to file(s): grid_[n].dat

/*
  print a usage statement
*/
void usage() {
  fprintf(stderr, "USAGE: count -n gridSize [-N trials] [-t] [-o]\n");
  fprintf(stderr, "-t: show timing info\n");
  fprintf(stderr, "-o: output grid to file\n");
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

  while ((c = getopt(argc, argv, "n:N:to")) != -1) {
    switch (c) {
    case 'n':
      gridSize = atoi(optarg);
      break;
    case 'N':
      trials = atoi(optarg);
      break;
    case 't':
      showTime = 1;
      break;
    case 'o':
      outputGrid = 1;
      break;
    case '?':
      switch (optopt) {
      case 'n':
      case 'N':
	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
	break;
      default:
	fprintf(stderr, "Unknown option -%c.\n", optopt);
      }
    default:
      abort();
    }
  }

  if (gridSize < 0) {
    ERROR("invalid grid size %d", gridSize);
    exit(CMD_LINE_ARG_ERR);
  }
}


/*

output:
        return value: number of nodal domains
*/
int runTest(double **grid, int ny, int nx) {   
  clock_t start = clock();
  int nd = countNodalDomainsNoInterp(grid, NULL, ny, nx);
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

  nx = gridSize;
  ny = gridSize;
  int counts[trials];

  double mean = 0.0, variance = 0.0;

  int i;
  grid = createGrid(ny, nx);
  if (!grid) {
    exit(OUT_OF_MEMORY_ERR);
  }
  srand(time(NULL));

  // run the trials and calculate mean count
  for (i = 0 ; i < trials ; i++) {
    randomPercolation(grid, ny, nx);
    
    if (outputGrid) {
      char outfile[50];
      sprintf(outfile, "../data/grid_%d.dat", trial++);
	array2file(grid, ny, nx, outfile);
    }
    
    counts[i] = runTest(grid, ny, nx);
    fprintf(stderr, "%d\n", counts[i]);
    mean += counts[i];
  }
  destroyGrid(grid);
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
 
  exit(0);
}
