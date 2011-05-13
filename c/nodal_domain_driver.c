/*
Driver for count_nodal_domains.c

Kyle Konrad
3/31/2011
*/


#include "count_nodal_domains.h"
#include "random_percolation.h"
#include "util.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// vergini code dependencies
#include "../vergini/billiard.h"
int verb;

// options specified by command line arguments
int showTime = 0; // flag: print time it takes for countNodalDomains to run
int gridSize = -1; // size of grid in x and y dimensions
int trials = 1; // number of grids to generate and count
int trial = 1; // current trial number
int outputGrid = 0; // flag: output grid(s) to file(s): grid_[n].dat
char *file; // use grid from this file
int mode = -1; // 0 - generate random_percolation grid; 1 - read grid from sta_bin file
int maskFlag = 0; // flag: use a mask file
char *maskFile; // use mask from this file
int oneFlag = 0; // flag: only count one eigenfunction

Billiard bil;
double dx = -1;
double k_0 = -1;

/*
  print a usage statement
*/
void usage() {
  fprintf(stderr, "USAGE: count {-n gridSize [-N trials] | -f file [-m maskFile | {-l billiardType -d dx [-k k_0]}]} [-t] [-o] [-1]\n");
  fprintf(stderr, "-t: show timing info\n");
  fprintf(stderr, "-o: output grid to file\n");
  fprintf(stderr, "-1: only count first eigenfunction\n");
}

/*
  process command line arguments
*/
void processArgs(int argc, char **argv) {
  int i = 0;

  if (argc <= 1) {
    usage();
    exit(-1);
  }

  while (++i < argc) {
    if (strcmp(argv[i], "-n") == 0) {
      gridSize = atoi(argv[++i]);
      mode = 0;
    }
    else if (strcmp(argv[i], "-N") == 0)
      trials = atoi(argv[++i]);
    else if (strcmp(argv[i], "-f") == 0) {
      file = (char *)malloc(strlen(argv[++i])*sizeof(char));
      strcpy(file, argv[i]);
      mode = 1;
    }
    else if (strcmp(argv[i], "-t") == 0)
      showTime = 1;
    else if (strcmp(argv[i], "-o") == 0)
      outputGrid = 1;
    else if (strcmp(argv[i], "-m") == 0) {
      maskFile = (char *)malloc(strlen(argv[++i])*sizeof(char));
      strcpy(maskFile, argv[i]);
      maskFlag = 1;
    }      
    else if (strcmp(argv[i], "-l") == 0) {
      if (parse_billiard(argv[++i], &bil) == -1) {
	fprintf(stderr, "Error: failed to parse billiard args\n");
	usage();
	exit(5);
      }
      mode = 2;
    }
    else if(strcmp(argv[i], "-d") == 0)
      dx = (double)atof(argv[++i]);
    else if(strcmp(argv[i], "-k") == 0)
      k_0 = (double)atof(argv[++i]);
    else if(strcmp(argv[i++], "-1") == 0)
      oneFlag = 1;
    else
      fprintf(stderr, "WARNING: unknown argument: %s\n", argv[i]);
  }

  if (gridSize < 0 && strcmp(file, "") == 0) {
    usage();
    exit(1);
  }
}


/*
generate a random ny x nx grid and count its nodal domains

output:
        return value: number of nodal domains
*/
int runTest(double **grid, char **mask, int ny, int nx) {   
  clock_t start = clock();
  int nd = countNodalDomains(grid, mask, ny, nx);
  clock_t end = clock();

  if (showTime)
    printf("countNodalDomains took %f seconds\n", ((double)(end - start)) / CLOCKS_PER_SEC);
  return nd;
}



int main(int argc, char **argv) {
  processArgs(argc, argv);
  int ny, nx;
  double **grid;

  if (mode == 0) {
    nx = gridSize;
    ny = gridSize;
    int counts[trials];

    double mean = 0.0, variance = 0.0;

    int i;
    // run the trials and calculate mean count
    for (i = 0 ; i < trials ; i++) {
      grid = createGrid(ny, nx);
      randomPercolation(grid, ny, nx);

      if (outputGrid) {
	char outfile[50];
	sprintf(outfile, "../data/grid_%d.dat", trial++);
	array2file(grid, ny, nx, outfile);
      }

      counts[i] = runTest(grid, NULL, ny, nx);
      destroyGrid(grid, ny);
      mean += counts[i];
    }
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

  }

  if (mode == 1) {
    int count;
    int masky, maskx;
    char **mask = NULL;
    grid = readOneSta(file, &ny, &nx);
    
    if (grid == NULL) {
      fprintf(stderr, "main: FATAL ERROR: failed to read grid\n");
      exit(3);
    }

    if (maskFlag) {
      mask = readMask(maskFile, &masky, &maskx);
     
      if (maskx != nx || masky != ny) {
	fprintf(stderr, "main: FATAL ERROR: mask dimensions do not match grid dimensions\n");
	exit(2);
      }
    }

    count = runTest(grid, mask, ny, nx);

    destroyMask(mask, ny);

    destroyGrid(grid, ny);
    free(file);
    if (maskFlag)
      free(maskFile);

    printf("counted %d nodal domains\n", count);
  }

  if (mode == 2) {
    int count;
    int k_base = 20; // to be passed to build_billiard
    int masky, maskx;
    char **mask = NULL;
    double k, wtm;
    int ne;

    if (dx <= 0) {
      fprintf(stderr, "Error: dx not specified or invalid\n");
      exit(6);
    }

    if (k_0 <= 0) {
      fprintf(stderr, "Error: k_0 not specified or invalid\n");
      exit(7);
    }

    if (build_billiard(&bil, k_base) != 0) {
      fprintf(stderr, "main: ERROR: failed to build billiard\n");
    }
    

    int i = 0;
    do {
      grid = readSta(file, &ne, &ny, &nx, &k, i);

      if (grid == NULL) {
	fprintf(stderr, "main: ERROR: failed to read grid\n");
	exit(3);
      }
    
      mask = createScaledMaskFromBilliard(bil, dx, &masky, &maskx, k/k_0);

      if (maskx != nx || masky != ny) {
	fprintf(stderr, "main: FATAL ERROR: mask dimensions do not match grid dimensions\n");
	fprintf(stderr, "ny\tnx\tmasky\tmaskx\n");
	fprintf(stderr, "%d\t%d\t%d\t%d\n",ny,nx,masky,maskx);
	exit(2);
      }

      count = runTest(grid, mask, ny, nx);

      wtm = wingTipMass(grid, mask, ny, nx);

      destroyMask(mask, ny);
      destroyGrid(grid, ny);

      printf("%d\t%f\t%f\n", count, k, wtm);

      if (oneFlag)
	break;

    } while (++i < ne);  

      free(file);
  }


  return 0;
}
