/* test collide module 6/30/06
 *
 * gcc -o test_collide test_collide.c collide.o -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "collide.h"
int verb;  // global verbosity

/* =================================== MAIN =============================== */
int main(int argc, char **argv)
{
  double y, d, nx, ny, R = 1.0;

  if (0) {
  printf("test lineseg_coll:\n");
  for (y=-1.0;y<2.0;y+=0.1) {
    d = lineseg_coll(0.5, y, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, sqrt(0.5));
    if (d==0.0)
      printf("y=%g: not hit\n", y);
    else    
      printf("y=%g: \td=%g\n", y, d);
  }
  }

  printf("\ntest arc_coll:\n");
  y = 0.5;
  d = arc_coll(0, y, 1.0, 0.0, 1.,0., R*R, 1./R, R,0, 0,R, &nx,&ny);
  printf("y=%g: \td=%g\n", y, d);

  if (0){
    for (y=-1.0;y<2.0;y+=0.1) {
      d = arc_coll(0, y, 1.0, 0.0, 1.,0., R*R, 1./R, R,0, 0,R,&nx,&ny);
      if (d==0.0)
	printf("y=%g: \td=%g\n", y, d);
    else    
      printf("y=%g: \td=%g, \tn=(%g,%g)\n", y, d, nx, ny);
    }
  }
  return 0;
}
