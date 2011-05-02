/*
  binmulti.cc

  bins multiple rectangular overlap matrices stored in a multiple 2d VIEWER file
  and a .sum12 file.

  Usage: binm file bin_width [d]

  Reads: ov/file.trs
         k/file.sum12

	 typical binwidth is 0.02

  Outputs: pnm/file.trs.mb (multiple binned datasets)
           Or, if [d] option present:
	   pnm/file.trs.pnm pnm/file.trs.comega
	   pnm/file.trs.xvals pnm/file.trs.dkvals (Doron format)

  Compile: gcc binmulti.cc -o binm -lm

  Barnett 00/3/29, 00/6/24, 00/6/26
  */

#include <stdio.h>
#include <math.h>

#define MAXBINS 2000
#define MAXS 200
#define MAXI 10000
#define SIZ 200
#define MAXHEAD 10




int output_bin_file(char (*header)[SIZ], char *in, char *ink, double dk,\
		    int ns, int il, int ih, int i_off, double (*ss)[MAXBINS], \
		    int (*num)[MAXBINS], double *x_lam)
{
  int i,s;
  double msy;
  char out[SIZ];
  FILE *fp;

  // decide output filename:
  sprintf(out,"pnm/%s.trs.mb",in);

  if ((fp = fopen(out,"w")) == NULL) {
    printf("file: %s, open for write error...\n",out);
    return 0;
  }

  for (i=0;i<5;++i)
    fputs(header[i],fp);
  fprintf(fp,"# Processed by binm, from files %s and %s\n", in, ink);
  fprintf(fp,"#\n# bin width dk=%f\n", dk);
  //fprintf(fp,"#\n# Columns: $1=x/lambda $2=dk $3=avg P(n|m) $4=err on $3\n#\n");
  fprintf(fp,"#\n# Columns: $1=x/lambda $2=dk $3=avg P(n|m)\n#\n");

  for (s=1; s<=ns; ++s) {
    for (i=il;i<=ih; ++i) {
      // mean sq y: deal with no samples in bin...
      msy = (num[s][i]==0) ? 0.0 : ss[s][i]/num[s][i];
      // fprintf(fp, "%.12g %.12g %.12g %.12g\n",x_lam[s], (i-i_off+0.5)*dk, msy,\
      //      msy/sqrt(1+num[s][i]));
      // for now, omit the errorbar column...
      fprintf(fp, "%.12g %.12g %.12g\n",x_lam[s], (i-i_off+0.5)*dk, msy);
    }
    fprintf(fp, "\n\n"); // increment gnuplot index for 2d multiple plot.
  }
  fclose(fp);
  return 1;
}



// -----------------------------------------------------------------------------
int output_doron_files(char (*header)[SIZ], char *in, char *ink, double dk,\
		       int ns, int il, int ih, int i_off, double (*ss)[MAXBINS], \
		       int (*num)[MAXBINS], double *x_lam)
{
  int i,s;
  double msy;
  char out[SIZ];
  FILE *fp;

  // -------------------------------------pnm-----looping over dk fast, over x slow
  sprintf(out,"pnm/%s.trs.pnm",in);
  if ((fp = fopen(out,"w")) == NULL) {
    printf("file: %s, open for write error...\n",out);
    return 0;
  }
  // omit the first x value since stores the FOPT result
  for (s=2; s<=ns; ++s) {
    for (i=il;i<=ih; ++i) {
      // mean sq y: deal with no samples in bin...
      msy = (num[s][i]==0) ? 0.0 : ss[s][i]/num[s][i];
      fprintf(fp, "%.12g\n", msy);
    }
    fprintf(fp, "\n"); // increment index
  }
  fclose(fp);

  // -------------------------------------rvals-----
  sprintf(out,"pnm/%s.trs.xvals",in);
  if ((fp = fopen(out,"w")) == NULL) {
    printf("file: %s, open for write error...\n",out);
  }
  for (s=2; s<=ns; ++s)
      fprintf(fp, "%.12g\n", x_lam[s]);
  fclose(fp);

  // -------------------------------------dkvals-----
  sprintf(out,"pnm/%s.trs.dkvals",in);
  if ((fp = fopen(out,"w")) == NULL) {
    printf("file: %s, open for write error...\n",out);
  }
    for (i=il;i<=ih; ++i)
      fprintf(fp, "%.12g\n", (i-i_off+0.5)*dk);
  fclose(fp);

  // -------------------------------------comega (FOPT)-----
  sprintf(out,"pnm/%s.trs.comega",in);
  if ((fp = fopen(out,"w")) == NULL) {
    printf("file: %s, open for write error...\n",out);
  }
  s=1; // first x value is FOPT result.
  for (i=il;i<=ih; ++i) {
    // mean sq y: deal with no samples in bin...
    msy = (num[s][i]==0) ? 0.0 : ss[s][i]/num[s][i];
    fprintf(fp, "%.12g\n", msy);
  }
  fclose(fp);

  return 1;
}




// ===============================================================================
main(int argc, char *argv[])
{
  int i, j, s, num[MAXS][MAXBINS], max_i, min_i, unsymflag = 0, i_off = MAXBINS/2;
  int max_max_i, min_min_i;
  int ndims, n1, n2, n1_v, n2_v, n_sets, ijflag, ii, jj, file_s;
  double dk, x, y, x_lam[MAXS], ss[MAXS][MAXBINS], dummy;
  double k1[MAXI], k2[MAXI], k_in;
  char out[SIZ], header[MAXHEAD][SIZ], in[SIZ], ink[SIZ], mode='a';
  FILE *fp;

  // cmd line read...
  if (argc<3 || argc>4) {
    printf("Incorrect usage! See ~/bdry/binmulti.cc.\n");
    return 0;
  }
  sprintf(in,"ov/%s.trs",argv[1]);
  sprintf(ink,"k/%s.sum12",argv[1]);
  sscanf(argv[2],"%lf",&dk);
  if (argc==4)
    sscanf(argv[3],"%c",&mode);
  printf("Multiple Binner: mode=%c\n",mode);

  // Input... read .sum12 file for k-values.
  if ((fp = fopen(ink,"r")) == NULL) {
    printf("binm: %s, open for read error...\n",in);
    return 0;
  }
  printf("reading %s...\n",ink);
  // read header... (adjust this to corresp to file type read)
  for (i=0;i<9;++i)
    fgets(header[i],SIZ,fp);
  // read k data... 
  ijflag = 1;
  while(fscanf(fp,"%d %lf %lf %lf\n",&i, &k_in, &dummy, &dummy)!=EOF) {
    if (i==1)
      ijflag = 1-ijflag;
    if (ijflag) {
      // reading k2... (deformed)
      k2[i] = k_in;
      n2 = i;
      //printf("k2[%d]=%f\n",n2,k2[n2]);
    } else {
      // reading k_i...
      k1[i] = k_in;
      n1 = i;
      //printf("k1[%d]=%f\n",n1,k1[n1]);
    }
  }
  fclose(fp);
  printf("done: n1=%d, n2=%d\n",n1,n2);


  // Input VIEWER file...
  if ((fp = fopen(in,"r")) == NULL) {
    printf("binm: %s, open for read error...\n",in);
    return 0;
  }
  printf("reading %s...\n",in);
  // read viewer file...
  fscanf(fp,"%d %d %d %d 2 %lf %lf",&ndims, &n1_v, &n2_v, &n_sets,
	 &dummy, &dummy);
  if (ndims!=2 || n1_v!=n1 || n2_v!=n2) {
    printf("Illegal file %s: ndims=%d, n1=%d, n2=%d, n_sets=%d\n", in, ndims,\
	   n1, n2, n_sets);
    printf("n1_v=%d, n2_v=%d\n", n1_v, n2_v);
    return(0);
  }



  // ..... LOOP OVER VIEWER 2d ARRAYS...

  printf("    looping over %d overlap matrices, binning...\n",n_sets);
  max_max_i = 0;
  min_min_i = MAXBINS-1;
  for (s=1;s<=n_sets;++s) {

    // read s, val[s] = param x in lambda units...
    fscanf(fp, "%d %lf\n",&file_s, x_lam + s);

    // zero the counting and sum of squares arrays...
    for (i=0;i<MAXBINS;++i) {
      num[s][i] = 0;
      ss[s][i] = 0.0;
    }
    max_i = 0;
    min_i = MAXBINS-1;
    
    // bin a single P_nm = |O_nm|^2 matrix: ii=undeformed, jj=deformed states.
    for (ii=1;ii<=n1;++ii)
      for (jj=1;jj<=n2;++jj) {
	fscanf(fp,"%lf",&y); // read in O[ii][jj] matrix el
	x = k2[jj] - k1[ii];
	if (fabs(x) > 1.0e-9) { // only count if not on diagonal (O_nm trusty).
	  i = (int)(floor(x/dk)) + i_off;
	  // bins i goes from x=(i-i_off)*dk to (i-i_off+1)*dk...
	  if (i>=MAXBINS) {
	    printf("Exceeded highest bin #!\n");
	    exit(0);
	  }
	  if (i<0) {
	    printf("Smaller than smallest bin #!\n");
	    exit(0);
	  }
	  ++num[s][i];
	  ss[s][i] += y*y; // square of matrix el gives Pnm.
	  if (i>max_i)
	    max_i = i;
	  if (i<min_i)
	    min_i = i;
	}
      }

    if (max_i>max_max_i)
      max_max_i = max_i;
    if (min_i<min_min_i)
      min_min_i = min_i;
    
    printf("s=%d: min_i=%d, max_i=%d\n",s,min_i,max_i);
  }
  fclose(fp);

  printf("over all s: min_min_i=%d, max_max_i=%d\n",min_min_i,max_max_i);

  // output...
  if (mode=='a')
    output_bin_file(header, argv[1], ink, dk, \
		    n_sets, min_min_i, max_max_i, i_off, ss, num, x_lam);
  else if (mode=='d')
    output_doron_files(header, argv[1], ink, dk, \
		       n_sets, min_min_i, max_max_i, i_off, ss, num, x_lam);

  return 1;
}
