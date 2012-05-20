/*
  Utility functions with Vergini code dependencies
 */


#include "util.h"
#include "util_verg.h"
#include "count_util.h"
#include "exit_codes.h"
#include "stdio.h"
#include <math.h>

extern int verb;

/*
create a mask for a billiard, scale it
input:
       b     - billiard to create mask for
       dx    - grid spacing
       scale - factor that grid is scaled by (k/k_0)
       ny_p  - where to save number of rows
       nx_p  - where to save number of columns
output:
       returns counted array with values of MASKED or UNCOUNTED
       *ny     - rows in array
       *nx     - columns in array
*/
int **createScaledMaskFromBilliard(Billiard b, double dx, double scale, int *ny_p, int *nx_p) {
  int nx, ny;
  *ny_p = ny = ceil((b.yh - b.yl) / dx) + 1;
  *nx_p = nx = ceil((b.xh - b.xl) / dx) + 1;

  int **counted = imatrix(ny, nx);
  MALLOC_CHECK(counted);

  int i, j;
  for (i = 0 ; i < ny ; i++) {
    for (j = 0 ; j < nx ; j++) {
      if (inside_billiard(j * dx / scale, i * dx / scale, &b)) {
        counted[i][j] = UNCOUNTED;
      } else {
        counted[i][j] = MASKED;
      }
    }
  }
    
  // special case for qugrs billiard: mask out the two points at the tip
  if (b.type == QU_GEN_RECT_SINAI) {
    counted[ny-1][nx-1] = MASKED;
    counted[ny-2][nx-2] = MASKED;
  }

  return counted;
}

/*
  read a mask from sta_bin file
  currently only reads single eigenfunction files

  input:
       file - name of file to read
       ny   - where to save number of rows
       nx   - where to save number of columns
  output:
       returns counted array with values of MASKED or UNCOUNTED
       *ny - rows in array
       *nx - columns in array
*/
int **createMaskFromFile(char *file, int *ny_p, int *nx_p) {
  FILE *fp = fopen(file, "r");
  if (fp == NULL) {
    ERROR("failed to open %s", file);
    return NULL;
  }

  char c;
  int n_e, nx, ny;
  double temp_double;
  float temp_float;

  fscanf(fp,"%c",&c);
  if (c == 'b') {
    fscanf(fp,"%c",&c);
    fscanf(fp,"%c",&c);
    fscanf(fp,"%c",&c);
  }
  else {
    fclose(fp);
    ERROR("incorrect file format in %s", file);
    return NULL;
  }

  if (fread(&n_e, (size_t)4,1,fp) != 1) {
    ERROR("failed to read n_e in %s", file);
    return NULL;
  }

  if (fread(nx_p, (size_t)4,1,fp) != 1) {
    ERROR("failed to read nx in %s", file);
    return NULL;
  }

  if (fread(ny_p, (size_t)4,1,fp) != 1) {
    ERROR("failed to read ny in %s", file);
    return NULL;
  } 

  if (n_e != 1) {
    ERROR("more than one eigenfunction in %s", file);
    return NULL;
  }

  nx = *nx_p;
  ny = *ny_p;
  int **counted = imatrix(ny, nx);
  MALLOC_CHECK(counted)

  if (fread(&temp_double,sizeof(double),1,fp) != 1) { // this is the energy - we don't care about it here
    ERROR("failed to read E_1 in %s", file);
    return NULL;
  }

  int i;
  for (i = 0 ; i < nx * ny ; i++) {
    if (fread(&temp_float,sizeof(float),n_e,fp) != (unsigned int)n_e) {
      ERROR("failed to read data in %s", file);
      return NULL;
    }
    if (temp_float) {
      counted[i/nx][i%nx] = UNCOUNTED;
    } else {
      counted[i/nx][i%nx] = MASKED;
    }
  }

  fclose(fp);
  return counted;
}
