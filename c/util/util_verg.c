/*
  Utility functions with Vergini code dependencies
 */


#include "util.h"
#include "util_verg.h"
#include <math.h>

/*
create a mask for a billiard
input:
       b    - billiard to create mask for
       dx   - grid spacing
output:
       returns - boolean mask array
       ny      - rows in array
       nx      - columns in array
*/
char **createMaskFromBilliard(Billiard b, double dx, int *ny, int *nx) {
  *ny = ceil((b.yh - b.yl) / dx) + 1;
  *nx = ceil((b.xh - b.xl) / dx) + 1;
  char **mask = createMask(*ny, *nx);

  int i, j;
  for (int i = 0 ; i < *ny ; i++)
    for (int j = 0 ; j < *nx ; j++)
      mask[i][j] = inside_billiard(j * dx, i * dx, &b);

  return mask;
}

/*
create a mask for a billiard, scale it
input:
       b     - billiard to create mask for
       dx    - grid spacing
       scale - factor that grid is scaled by (k/k_0)
output:
       returns - boolean mask array
       ny      - rows in array
       nx      - columns in array
*/
char **createScaledMaskFromBilliard(Billiard b, double dx, int *ny, int *nx, double scale) {
  *ny = ceil((b.yh - b.yl) / dx) + 1;
  *nx = ceil((b.xh - b.xl) / dx) + 1;
  char **mask = createMask(*ny, *nx);

  int i, j;
  for (i = 0 ; i < *ny ; i++)
    for (j = 0 ; j < *nx ; j++)
      mask[i][j] = inside_billiard(j * dx / scale, i * dx / scale, &b);

  // special case for qugrs billiard: mask out the two points at the tip
  if (b.type == QU_GEN_RECT_SINAI) {
    mask[*ny-1][*nx-1] = 0;
    mask[*ny-2][*nx-2] = 0;
  }

  return mask;
}
