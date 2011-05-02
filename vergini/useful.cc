/* Useful routines */

#include <stdio.h>
#include <string.h>

#include "useful.h"
#include "verb.h"


double lower_limit(double *a, int M)
{
        int i;
        double limit = 1e+308;
        
        for (i=1;i<=M;++i)
                if (a[i]<limit)
                        limit = a[i];

        return limit;
}

double upper_limit(double *a, int M)
{
        int i;
        double limit = -1e+308;
        
        for (i=1;i<=M;++i)
                if (a[i]>limit)
                        limit = a[i];

        return limit;
}

float flower_limit(float *a, int M)
{
        int i;
        float limit = 1e+38;
        
        for (i=1;i<=M;++i)
                if (a[i]<limit)
                        limit = a[i];

        return limit;
}
float fupper_limit(float *a, int M)
{
        int i;
        float limit = -1e+38;
        
        for (i=1;i<=M;++i)
                if (a[i]>limit)
                        limit = a[i];

        return limit;
}


void endian_array_4byte(char *array, int n)
/* flips 4 byte arrays big<->little endian, n values, beginning where
 * ptr points. Processes 4n bytes in total.
 */
{
	char *ptr,temp;
	int i;

	ptr = array;
	for (i=0; i<n; ++i,ptr+=4) {
		temp = ptr[0];
		ptr[0] = ptr[3];
		ptr[3] = temp;
		temp = ptr[1];
		ptr[1] = ptr[2];
		ptr[2] = temp;
	}
}

int endian_4byte(void *in)
/* returns flipped 4 byte value big<->little endian, without destroying orig.
 * Incoming arg is a pointer, returns value. 99/10/15
 */
{
	char *inc,outc[4];

	inc = (char *)in;
	outc[0] = inc[3];
	outc[1] = inc[2];
	outc[2] = inc[1];
	outc[3] = inc[0];

	return *((int *)outc);
}

double endian_8byte(double *in)
/* returns flipped 8 byte double big<->little endian, without destroying orig.
 * Incoming arg is a pointer, returns value. 99/10/15.
 */
{
	char *inc,outc[8];

	inc = (char *)in;
	outc[0] = inc[7];
	outc[1] = inc[6];
	outc[2] = inc[5];
	outc[3] = inc[4];
	outc[4] = inc[3];
	outc[5] = inc[2];
	outc[6] = inc[1];
	outc[7] = inc[0];

	return *((double *)outc);
}


int round_and_clip(double x, char *var_name)
/*
  used to read nearest integer to a double, and clip to >= 1.
  var_name is text label for error feedback.

  9/8/03 changed from make_int_and_clip, moved to useful.
*/
{
  int n = (int)( 0.5 + x );
  
  if (n<1) {
    fprintf(stderr,"%s invalid = %d, make_int_and_clip!\n",var_name,n);
    return 1;
  } else
    return n;
}



// -------------------------------------------------------------------------
int separated_numbers(char *s, double *a, int n)
  /* Read in string as numbers into an array, when there is a separator
   * (colon)
   *
   * barnett 11/16/03
   */
{
  int i, c, off, len = strlen(s);

  if (len<1)
    if (n==0)
      return 0;
    else {
      fprintf(stderr,\
	      "\tNot enough params (0), (needs %d)!\n", n);
      return -1;
    }

  for (i=0, off=0; (i<n && off<len); ++i, off+=c+1) {
    if (verb&0x80)
      printf("i=%d remain=%s\n", i, s+off);
    sscanf(s+off, "%lf%n:", a + i, &c);  // note separator defined here
    if (verb&0x80)
      printf("off=%d c=%d\n", off, c);
  }
  if (verb&0x80)
    printf("off=%d, len=%d\n", off, len);
  /* Check for right number of args... */
  if (i<n) {
    fprintf(stderr,\
      "\tNot enough params (%d), (needs %d)!\n", i, n);
    return -1;
  }
  if (off<=len) {
    fprintf(stderr,\
      "\tToo many params (needs %d)! Remaining chars: %s\n", \
	    i, s+off);
    return -1;
  }

  return n;
}

int grab_cmdline(char *cmdline, int argc, char **argv)
{
  int i;
  char *spacer = " ";

  cmdline[0] = '\0';  // set to empty string
  for (i=0; i<argc; ++i) {
    strcat(cmdline, argv[i]);
    strcat(cmdline,spacer);
  }
  return 0;
}


// ------------------------- interpolation ---------------------------------

double inverse_interp(int M, double *xs, double *ys, double y)
  /* Return x corresponding to y, using piecewise linear approximation to
     monotonically-increasing function y = f(x) given by xs[j], ys[j],
     j = 0...M.
     'Inverse' arises because we are finding x = f^{-1}(y).
     Binary search on ys values is performed.
     barnett 12/13/03
  */
{
  int j, jl = 0, jh = M;
  double frac;

  while (jh-jl > 1) {
    j = (jl+jh)/2;   // integer arithmetic
    if (verb&0x80)
      printf("\t\t\tjl=%d, jh=%d, ys[%d]=%g\n", jl, jh, j, ys[j]);
    if (ys[j] > y)
      jh = j;
    else
      jl = j;
  }
  frac = (y - ys[jl])/(ys[jh] - ys[jl]);
  if (verb&0x80)
      printf("\t\t\tjl=%d, jh=%d, frac=%f\n", jl, jh, frac);
  return (1.0-frac)*xs[jl] + frac*xs[jh];
}
