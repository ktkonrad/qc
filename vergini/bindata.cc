/*
  bindata.c bins (x,y) data by x, writing out 3-column 'yerrorbar' for gnuplot.
  This estimates the mean square of y in each x-bin.
  Output filename has .b appended.

  Use: bindata file bin_width

  If bin_width is < 0, this tells it not to take the absval of x, otherwise it does. 

  Compile: cxx bindata.cc -o bindata -lm

  Barnett 00/1/18
  */

#include <stdio.h>
#include <math.h>

#define MAXBINS 100000

main(int arcg, char *argv[])
{
  int i, num[MAXBINS], max_i, min_i, unsymflag = 0, i_off = MAXBINS/2;
  double dx, x,y, ss[MAXBINS], msy;
  char out[100], header[10][300], *in = argv[1];
  FILE *fp;

  // cmd line read...
  sprintf(out,"%s.b",in);
  sscanf(argv[2],"%lf",&dx);

  // un-symm ?
  if (dx<0) {
    dx = -dx;
    unsymflag = 1;
  }

  // Input...
  if ((fp = fopen(in,"r")) == NULL) {
    printf("file: %s, open for read error...\n",in);
    return 0;
  }

  // read header... (adjust this to corresp to file type read)
  for (i=0;i<6;++i)
    fgets(header[i],300,fp);

  // zero the counting and sum of squares arrays...
  for (i=0;i<MAXBINS;++i) {
    num[i] = 0;
    ss[i] = 0.0;
  }
  max_i = 0;
  min_i = MAXBINS-1;

  // MAIN LOOP: read in and add to a bin...
  while (fscanf(fp,"%lf %lf",&x,&y)==2) {
    if (!unsymflag)
      x = fabs(x);
    if (fabs(x)  > 2.0e-4) { // only count if not on diagonal (O_nm trusty).
      i = (int)floor(x/dx) + i_off;   // bins i goes from
      // x=(i-i_off)*dx to (i-i_off+1)*dx...
      if (i>=MAXBINS) {
	printf("Exceeded highest bin #!\n");
	exit(0);
      }
      if (i<0) {
	printf("Smaller than smallest bin #!\n");
	exit(0);
      }
      ++num[i];
      ss[i] += y*y;
      if (i>max_i)
	max_i = i;
      if (i<min_i)
	min_i = i;
      if (!unsymflag) { // duplicate bins on -ve dk side ....
	x = -x;
      i = (int)floor(x/dx) + i_off;   // bins i goes from
      // x=(i-i_off)*dx to (i-i_off+1)*dx...
      if (i>=MAXBINS) {
	printf("Exceeded highest bin #!\n");
	exit(0);
      }
      if (i<0) {
	printf("Smaller than smallest bin #!\n");
	exit(0);
      }
      ++num[i];
      ss[i] += y*y;
      if (i>max_i)
	max_i = i;
      if (i<min_i)
	min_i = i;
      }
    }
  } ;

  fclose(fp);

  printf("Bindata: done, found %d bins from bin # %d to %d, x_min=%g, x_max=%g\n",\
	 1+max_i-min_i, min_i, max_i, (min_i-i_off)*dx, (max_i-i_off)*dx);

  // Output...
  if ((fp = fopen(out,"w")) == NULL) {
    printf("file: %s, open for write error...\n",out);
    return 0;
  }

  // write header
  for (i=0;i<6;++i)
    fprintf(fp,"%s",header[i]);

  fprintf(fp,"# Rebinned by bindata with bin_width=%e\n",dx);
  if (unsymflag)
    fprintf(fp,"# Unsymmetric binning.\n",dx);
  fprintf(fp,"# Columns : $1=center_x,$2=mean sq O,$3=errorbar on $2\n\n");

  // force a certain number of bins for doron format...
  min_i = i_off - 2000;
  max_i = i_off + 1999;

  // write 3-column format: x_cen, mean sq y, sigma_y (due to sampling)
  for(i=min_i;i<=max_i;++i) {
    msy = (num[i]==0) ? 0.0 : ss[i]/num[i]; // deal with no samples in bin.
    fprintf(fp, "%.12g %.12g %.12g\n",(i-i_off+0.5)*dx, msy, msy/sqrt(1+num[i]));
  }
 
  fclose(fp);

}
