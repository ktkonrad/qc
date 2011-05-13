/*
 *	program:	viewer.c		executable:	viewer
 *
 *	Multiple 1D and 2D Array Viewer.
 *	Designed with eigenstates and corresponding eigenvalues in mind.
 *
 *	Call with arguments: viewer infile.sta [-e] [GLUT-options]
 *
 *	CURRENT VERSION NOTES:
 *	1d and 2d arrays implemented
 *	expects 2 params exactly.
 *	reads ascii format .sta files AND binary .sta_bin files.
 *
 *
 *	Alex Barnett	7/27/99	(based on soft.c 5/22/1997).
 			99/8/28 changed vec data to float, read .sta_bin files.
			99/11/4 ampl/intensity display options for images.
			99/11/6 antialising, points 1d display options.
			00/3/28 port to alpha linux, allocation bug fix.
			2/28/02 i386, ENDIAN cmdline option.
			4/28/02 1d-mods: spect view cmd line
			10/30/03 colour encoding of Inf values.
			2/25/04 ps/eps output of 2d images.
 */



#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "constants.h"
#include "useful.h"


// largest ratio to 1 of input data allowed
#define MAX_FLOAT 1.0e37

#define PI 3.1415926535897931

/* DATA GLOBALS :----(unbelievably bad global names!)---------------------*/

/* user-available params to read in ... */
double param1,param2;

/* array for filename... */
char filename[100];

/* data-related: #dimensions and array sizes */
int ndims,n,nx,ny,n_e;

/* lattice spacing in real space, aspect, xoffset for 2d graphs */
double xlat, aspectratio, xoffset;

/* the values data */
double *val;

/* the vectors data */
float **vec;

/* 0 for usual endian, 1 for SGI... */
int ENDIAN_FLIP;


/* GRAPHICS GLOBALS -------------------------------------------------------*/

/* animation state... */
int freely_spinning,dragging,scaling,spec_dragging,mouse_still_moving;

/* 3D view globals: */
float th,ph,rh;
/* 2 * tangent of half-width of view in Y-direction... */
float iris,dist;
float zscale,thdot,phdot;

/* 2D image view globals: */
int xof2d,yof2d,xof2ddot,yof2ddot, img_ang;

/* 2D graph view globals: */
double xof,yof,xofdot,yofdot,xview,yview;

/* window size, position */
int width, height, xpos, ypos;
/* spectrum centre energy, and range... */
float eview, erange;

/* old mouse position... */
int omx, omy;
/* viewable state # */
int iv;

/* display toggles */
int contourtoggle;
int disp_spec;

/* dispmode has NUMBER_DISPLAY_MODES settings... */
#define NUMBER_DISPLAY_MODES	3
int dispmode;

/* logtoggle has NUMBER_LOG_MODES settings... */
#define NUMBER_LOG_MODES	4
int logtoggle;






void ratio(void)
/* replaces state i= 2 to n_e by the point-by-point ratio of state i to 1,
	divided by val[i].
	This is computing the ratio of actual O_nm to FOPT estimate.
 */
{
	int i,j;
	
	for (i=n_e;i>1;--i)
		for (j=1;j<=n;++j)
			vec[i][j] /= val[i]*vec[1][j];
}

void unratio(void)
/* undoes whatever ratio() does!
 */
{
	int i,j;
	
	for (i=n_e;i>1;--i)
		for (j=1;j<=n;++j)
			vec[i][j] *= val[i]*vec[1][j];
}

void subtract(void)
/* 
	This is computing the difference of actual O_nm from FOPT estimate,
	divided by val[i]^2.
 */
{
	int i,j;
	
	for (i=n_e;i>1;--i)
		for (j=1;j<=n;++j)
			vec[i][j] = (vec[i][j] - val[i]*vec[1][j]) / (val[i]*val[i]);
}

void unsubtract(void)
/* 
	undoes subtract()!
 */
{
	int i,j;
	
	for (i=n_e;i>1;--i)
		for (j=1;j<=n;++j)
			vec[i][j] = val[i]*val[i]*vec[i][j] + val[i]*vec[1][j];
}






/* --------------------------INFO DUMP------------------------------------ */
void infodump()
/* dumps to stdout... */
{
	int i;

	printf("viewer infodump : %s\n\n",filename);
	if (ndims==1) {
		printf("%d states: %d samples, xoff = %.6g, xlat = %.6g\n\n",\
			n_e,n,xoffset,xlat);
	} else {
		printf("%d states: %d by %d samples, xlat = %.6g\n\n",\
			n_e,nx,ny,xlat);
	}
	printf("param1 = %e, param2 = %e\n\n",param1,param2);

	/* output evals... could output other properties here! */
	printf("index       value\n\n");
	for(i=1;i<=n_e;++i)
		printf("%d  %.15g\n",i,val[i]);
	printf("\n");

}


/* --------------------------FILE HANDLING---------------------------------- */

/* maximum total number of datapoints allowed... 400 Mb */
#define MAX_DATA_SIZE 5e7


void change_size(int nx_old,int ny_old, int n_e_old)
/* call when global sizes nx, ny, n_e may have changed; reallocates space using
 * the NumRec nrutil.c array-handling routines, if necessary.
 */
{
    int n_old;
    
/* the nx+1, ny+1 are arbitrary, purely for display purposes... */
    if (ndims==1) {
	/* 1d array case */
        n = nx;
        n_old = nx_old;
        /* xlat = 1.0/(nx+1); */
    } else {
	/* 2d array case: keeps lattice square */
        n = nx*ny;
        n_old = nx_old*ny_old;
        xlat = 1.0/(nx+1);
        aspectratio = (ny+1)/(float)(nx+1);
    }
    
/* catch maximum sizes allowed... */
    if (n*n_e > MAX_DATA_SIZE) {
        printf("Requested data size of %d exceeds 40Mb, in change_size\n",\
               n*n_e);
        n_e = (int)(MAX_DATA_SIZE/n - 0.5);
        printf("Changing n_e to %d...\n",n_e);
    }
    
    
/* reallocate if either n or n_e has changed... */
    if (n!=n_old || n_e!=n_e_old) {
        if (n_old != 0 && n_e_old!=0) {
	  free_matrix(vec,1,n_e_old,1,n_old);
	  free_dvector(val,1,n_e_old);
	}
        vec = matrix(1,n_e,1,n);
        val = dvector(1,n_e);
    }
}



int loadstates()
/* loads evecs and evals, overwrites previous data.
 *
 * reads evecs and evals from 'filename', resets n, n_e if needed,
 * reallocating all array sizes to fit loaded data size.
 * Also can handle unallocated situation, when ndims is fixed once and for all.
 *
 * returns 0 if loading failed.
 *
 * Added 1d multi-graph read 99/10/19
 */
{
  FILE *fp;
  
  int i,k, file_k, onx, ony, on, on_e, ondims, temp_int;
  int n_params_dummy, size;
  int binary, isok;
  long int badcount=0;
  char first_char;
  float *temp;
  double temp_double, v;


  printf("in loadstates.\n");

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"viewer: file %s not found.\n",filename);
    return 0;
  }
  
  /* old values... */
  onx = nx;
  ony = ny;
  on_e = n_e;
  ondims = ndims;
  
  /* determine if binary (only ndims=2)... and load in appropriate manner. */
  fscanf(fp,"%c",&first_char);
  if (first_char == 'b') {
    binary = 1;
    fscanf(fp,"%c",&first_char);
    fscanf(fp,"%c",&first_char);
    fscanf(fp,"%c",&first_char);
    ndims = 2;
    fprintf(stderr,"%s, %dD binary dataset: ",filename,ndims);
  } else {
    binary = 0;
    /* convert from ASCII char to integer... */
    ndims = (int)first_char - 48;
    fprintf(stderr,"%s, %dD ascii dataset: ",filename,ndims);
  }
  
  /* sets ndims if it's first reading, checks it's unchanged otherwise... */
  if (ondims != 0 && ndims != ondims) {
    fprintf(stderr,"\nError loading %s: ndims cannot change from %d\n",\
	    filename,ondims);
    fclose(fp);
    return 0;
  } else {
    /* continue to read as normal... */
    
    if (ndims==1) {
      
      fscanf(fp,"%d %d %lf %lf %d %lf %lf",&nx,&n_e,&xoffset,&xlat,
	     &n_params_dummy,&param1,&param2);
      ny = 1;
      
      /* check valid sizes... */
      if (nx < 2 || n_e < 1) {
	fprintf(stderr,"n\viewer: (n=%d,n_e=%d) too small in %s.\n",\
		nx, n_e, filename);
	return 0;
      }
      if (nx*n_e > MAX_DATA_SIZE) {
	fprintf(stderr,\
		"\nviewer: too large (n=%d,n_e=%d) in %s.\n",\
		nx, n_e, filename);
	return 0;
      }
      /* n is set to nx*1 here... */
      change_size(onx, ony, on_e);
      printf("n=%d, n_e=%d, xoff=%.4g, xlat=%.4g\n",n,n_e,xoffset,xlat);
      
    } else if (ndims==2) {
      
      if (binary) {
	if (ENDIAN_FLIP) {
	  fread(&temp_int,(size_t)4,1,fp);
	  n_e = endian_4byte(&temp_int);
	  fread(&temp_int,(size_t)4,1,fp);
	  nx = endian_4byte(&temp_int);
	  fread(&temp_int,(size_t)4,1,fp);
	  ny = endian_4byte(&temp_int);
	} else {
	  fread(&n_e,(size_t)4,1,fp);
	  fread(&nx,(size_t)4,1,fp);
	  fread(&ny,(size_t)4,1,fp);
	}
	param1 = param2 = 0.0;			
      } else {
	fscanf(fp,"%d %d %d %d %lf %lf",&nx,&ny,&n_e,
	       &n_params_dummy,&param1,&param2);
      }
      
      /* check valid sizes... */
      if (nx < 2 || ny < 2 || n_e < 1) {
	fprintf(stderr,\
		"\nviewer: one of (nx=%d,ny=%d,n_e=%d) too small in %s.\n",\
		nx, ny, n_e, filename);
	return 0;
      }
      if (nx*ny*n_e > MAX_DATA_SIZE) {
	fprintf(stderr,\
		"\nviewer: too large (nx=%d,ny=%d,n_e=%d) in %s.\n",\
		nx, ny, n_e, filename);
	return 0;
      }
      
      /* n is set to nx*ny here... */
      change_size(onx, ony, on_e);
      printf("nx=%d ny=%d (n=%d), n_e=%d ...\n",nx,ny,n,n_e);
      
    }
    if (!binary)
      printf("param1 = %e, param2 = %e\n",param1,param2);
    
    if (binary) {
      
      if (ENDIAN_FLIP) {
	/* binary read, endian-flipped loop through ordinate,
	   then through states... */
	for (k=1; k<=n_e; ++k) {
	  fread(&temp_double,sizeof(double),1,fp);
	  val[k] = endian_8byte(&temp_double);
	}
	temp = vector(1,n_e);
	for (i=1; i<=n; ++i) {
	  fread(temp+1,4,n_e,fp);
	  endian_array_4byte((char *)(temp+1),n_e);
	  for (k=1; k<=n_e; ++k)
	    vec[k][i] = temp[k];
	}
	free_vector(temp,1,n_e);

      } else {
	/* binary read, loop through ordinate, then through states... */
	fread(val+1,sizeof(double),n_e,fp);
	temp = vector(1,n_e);
	for (i=1; i<=n; ++i) {
	  fread(temp+1,4,n_e,fp);
	  for (k=1; k<=n_e; ++k)
	    vec[k][i] = temp[k];
	}
	free_vector(temp,1,n_e);
      }
      
    } else {
      
      /* looping through states, load in each... (same for 1d & 2d)
	 Removal of illegal or out-of-range float values is done here too.
	 (see MAX_FLOAT definition at top of code).
      */
      for (k=1; k<=n_e; ++k) {
	fscanf(fp,"%d ",&file_k);
	if (k!=file_k) {
	  fprintf(stderr,\
	    "viewer: state # mismatch, expected %d, loaded %d, in file %s.\n",\
	    k,file_k,filename);
	  k=n_e+1; /* ends the loop */
	} else {
	  fscanf(fp,"%lf ",&v);
	  // check if val is finite, set to zero if not...
	  val[file_k] = finite(v) ? v : 0.0;
	  for (i=1; i<=n; ++i) {
	    fscanf(fp,"%lf ",&v); // scan vec in as double, convert to float.
	    isok = (v==0.0) || (finite(v) && fabs(v)>(1.0/MAX_FLOAT) && \
	      fabs(v)<MAX_FLOAT);
	    vec[file_k][i] = isok ? (float)v : 0.0;
	    if (!isok)
	      ++badcount;
	  }
	}
	
      }
    }
    
    if (badcount>0)
      printf("%d data values too large or small, set to zero.\n",badcount);
    printf("loaded.\n");
    fclose(fp);

    return 1;
  }
}



/* -------------------------------- postscript colourmap ------------ */
   /* colour scale... not added wrap-around of out-of-[0,1] values */
float rc(float a) /* RED */
{
  float b = a; /* - floorf(a); */
  return min(1.0,max(0.0,5*b-1.0));
}

float gc(float a) /* GREEN */
{
  float b = a; /* - floorf(a); */
  return min(1.0,max(0.0,5*b-3.0));
}

float bc(float a) /* BLUE */
{
  float b = a; /* - floorf(a);*/
   if (b<0.2)
    return 5*b;
  else if (b<0.8)
    return min(1.0,max(0.0,3.0-5*b));
  else
    return 5*b - 4.0;
}

int byte(float a)
     /* clip & map float from 0-1 to integer in 0-255 range */
{
  return min(0xff, max(0, (int)(a*255.0 + 0.5)));
}


/* ----------------------------------------------- DUMP POSTSCRIPT ------ */
int dump_ps(int i_lo, int i_hi, int n_wid, double gapfrac, int encap, \
	    int colour, int val_label)
 /* Writes postscript output for range [i_lo,i_ho] of states, with n_wid
    images written across, gapfrac fractional gap between them.
    Aspect ratio preserved (square image pixels).
    encap    = true : EPS, otherwise PS.
    colour   = true : use colourmap in rc,gc,gb, otherwise grayscale
    val_label= true : include value labels to 6 digits below each image

    2/26/04 Barnett
    2/27/04 multiple images and value labels
 */
{
  int i, j, l, lx = 0, ly = 0, ne = i_hi - i_lo + 1, nh;
  double xs = 10.0;  /* x-size of one sta image in inches */
  double ys = (xs*ny)/nx;
  double pti = 72.0;   /* pts per inch in PS language */
  double g = gapfrac*xs; /* gap in inches */
  float a, b, *dispvec;
  FILE *fp;
  char ofile[256];

  /* compute nh = # rows in image array... */
  nh = ne/n_wid;
  if (ne > n_wid*nh)
    ++nh;

  if (encap)
    sprintf(ofile, "%s.eps", filename);
  else
    sprintf(ofile, "%s.ps", filename);
  
  if ((fp = fopen(ofile,"w")) == NULL) {
    fprintf(stderr, "Problem writing to file %s ...\n", ofile);
    return 1;
  }

  /* Write header, bound box (including text labels), font size... */
  printf("Writing %d by %d PS image array to %s...\n", n_wid, nh, ofile);
  if (encap)
    fprintf(fp, "%!PS-Adobe-3.0 EPSF-3.0\n%%%%BoundingBox: %f %f %f %f\n", \
	    g*pti, val_label ? 0.0 : g*pti, (xs+g)*n_wid*pti, (ys+g)*nh*pti);
  else
    fprintf(fp, "%!PS-Adobe-3.0\n");
  fprintf(fp, "% Written by viewer v.2/27/04, Alex Barnett\n\n");
  if (colour)
    fprintf(fp, "/DeviceRGB setcolorspace\n");
  else
    fprintf(fp, "/DeviceGray setcolorspace\n");
  fprintf(fp, "/i {%g mul} def\n", pti);
  /* choose text size relative to image size... */
  fprintf(fp, "/Times-Roman findfont 0.03 scalefont setfont\n\n");

  /* Loop over images.................................................... */
  for (l=i_lo; l<=i_hi; ++l) {
    printf("l=%d: image array loc (%d,%d)\n", l, lx, ly);
    /* move and scale for current image... */
    fprintf(fp, "gsave\n%g i %g i translate %g i %g i scale\n", \
	    g + (g+xs)*lx, g + (g+ys)*ly, xs, ys);
    /* write type-1 image dictionary... */
    fprintf(fp, "<<\n /ImageType 1\n /Width %d\n /Height %d\n", nx, ny);
    /* data ordering: left-to right fast, then top to bottom slow... */
    fprintf(fp, " /BitsPerComponent 8\n /ImageMatrix [%d 0 0 %d 0 0]\n", \
	    nx, ny);
    if (colour)
      fprintf(fp, " /Decode [0 1 0 1 0 1]\n");
    else
      fprintf(fp, " /Decode [0 1]\n");
    fprintf(fp, " /Interpolate true\n");
    fprintf(fp, " /DataSource currentfile /ASCIIHexDecode filter\n>>\nimage\n");

    /* Write single image data from l'th vector in vec array... */
    dispvec = vec[l];
    for (j=1; j<=ny; ++j) {
      for (i=1; i<=nx; ++i) {
	/* Plot square of func, use current zscale... */
	a = zscale * (b=dispvec[nx*(j-1) + i])*b;
	/* write ASCII-hex data... (spaces, LF/CR irrelevant) */
	if (colour)
	  fprintf(fp, "%.2x%.2x%.2x", byte(rc(1-a)), byte(gc(1-a)), \
		  byte(bc(1-a)));
	else
	  fprintf(fp, "%.2x", byte(1-a));
      }
      fprintf(fp, "\n");
    }

    /* write a > to close filter, then text label (must fit in BBOX!) */
    fprintf(fp, ">\n");
    if (val_label)
      fprintf(fp, "0.3 -0.05 moveto (k = %.6f) show\n", val[l]);
    fprintf(fp, "grestore\n\n");

    /* update location (lx, ly) in image array... */ 
    ++lx;
    if (lx==n_wid) {
      lx = 0;
      ++ly;
    }
  } /* end image loop ................................................ */
  

  /* Trailer... showpage should be ignored in EPS */
  fprintf(fp, "showpage\n%%EOF\n");
  fclose(fp);
  printf("done.\n");
}




/* --------------------------GRAPHICS ROUTINES------------------------------ */
void initialize_view()
{
	double emin,emax; 


	freely_spinning = dragging = scaling = spec_dragging = 0;
	mouse_still_moving = 0;

	/* dispmode = 2; */
	logtoggle = 1;
	contourtoggle = 0;

/* 3d */
	dist = 2.5;
	/* th is azimuth, ph elevation, rather confusingly... */
	th = 0.3;
	ph = 1.0;
	rh = 0;
	thdot = 0.0;
	phdot = 0.0;
	iris = 0.6;
	zscale = 0.3;

/* 2d images */
	xof2d = 20;
	yof2d = 30;
	xof2ddot = yof2ddot = 0;
	img_ang = 0.0;

/* 2d graphs */
	xof = xoffset + n*xlat/2.0;
	yof = 0.0;
	xofdot = yofdot = 0.0;
	xview = n*xlat;
	yview = 1.0;
	/* tie zscale to yview for 1d graph plots... */
	if (ndims==1)
		zscale = 1.0/yview;

/* use val array to display whole val range: */
	emin = lower_limit(val,n_e);
	emax = upper_limit(val,n_e);

	/* centre of spectrum view */
	eview = (emin+emax)/2.0;
	/* vertical view range of spectrum */
	erange = 1.1*(emax-emin);
	if (erange <= 0.0)
	  erange = 1.0;
       
/* visible state index */
	iv = 1;
}


void texto(GLshort x, GLshort y, char *text)
/* writes bitmap text into graphics window, assuming already in 2d pixel
 * coords projection. Uses GLUT.
 */
{
    char *p;
    
    glRasterPos2s(x,y);
    for (p = text; *p; p++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);
}

void ftexto(float x, float y, char *text)
/* writes bitmap text into graphics window, in floating coords. Uses GLUT.
 * The text is not yet centered at x,y.
 */
{
    char *p;
    
    glRasterPos2f(x,y);
    for (p = text; *p; p++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);
}



#define AXES_LABEL_DIST 1.1

void axes3d()
/* labelled axes draw size 1,1,1 */
{
    glColor3f(0.3,0.3,1.0);
    glBegin(GL_LINE_STRIP);
    glVertex3f(1.0,0.0,0.0);
    glVertex3f(0.0,0.0,0.0);
    glVertex3f(0.0,1.0,0.0);
    glEnd();
    glBegin(GL_LINES);
    glVertex3f(0.0,0.0,0.0);
    glVertex3f(0.0,0.0,1.0);
    glEnd();
    
    glRasterPos3f(AXES_LABEL_DIST,0.0,0.0);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'x');
    glRasterPos3f(0.0,AXES_LABEL_DIST,0.0);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'y');
    glRasterPos3f(0.0,0.0,AXES_LABEL_DIST);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'z');
}


int choose_tics(double lo, double range, double *tics)
/* returns the number of tics, and their locations in zero-indexed tics[].
 * Barnett 99/10/20
 */
{
    int i, n_tics, tic_start;
    double exponent, logr, spacing;
    
    /* adjust the range multiplier here to give good tic density... */
    logr = log10(range * 0.4);
    exponent = floor(logr);
    spacing = pow(10.0, exponent);
    if (logr-exponent > log10(5.0))
        spacing *= 5.0;
    else if (logr-exponent > log10(2.0))
        spacing *= 2.0;
    
    /* (int) and copysign trick is to convert the floor val to an int... */
    tic_start = (int)(copysign(0.5,lo) + 1.0 + floor(lo / spacing));
    
    n_tics = (int)(1.0 + (lo + range - tic_start*spacing)/spacing);
    for (i=0;i<n_tics;++i) {
        tics[i] = spacing * (tic_start + i);
    }
    return n_tics;
}


void axes2d()
/* For 1d graph plots: crosshairs at screen center
 * Barnett 99/10/20
 */
{
    int n_tics, i, len;
    char label[20];
    float hair_squish;
    double x_pix_size, y_pix_size;
    double window_xo, window_yo, box_xl,box_xh, box_yl,box_yh;
    double tics[20], xtic_bot, xtic_top, ytic_bot, ytic_top;
    
    glColor3f(0.3,0.3,1.0);
    
    /* crosshairs */
    hair_squish = sqrt(height/(float)(width));
    glBegin(GL_LINES);
    glVertex2d(xof - hair_squish*xview/20, yof);
    glVertex2d(xof + hair_squish*xview/20, yof);
    glVertex2d(xof, yof - yview/(20*hair_squish));
    glVertex2d(xof, yof + yview/(20*hair_squish));
    glEnd();
    
    /* set up axes & tic locations (box_* coords give axes "box")... */
    x_pix_size = xview/width;
    y_pix_size = yview/height;
    window_xo = xof - xview/2.0;
    window_yo = yof - yview/2.0;
    box_xl = window_xo + 0.1*xview;
    box_yl = window_yo + 0.12*yview;
    box_xh = window_xo + 0.9*xview;
    box_yh = window_yo + 0.9*yview;
    
    /* note xtic_* is a y-coord, and visa versa... */
    /* choosing where axes lines appear */
    xtic_top = box_yl;
    if (box_yl<0.0)
		xtic_top = 0.0;
	if (box_yh<0.0)
		xtic_top = box_yh;
	xtic_bot = xtic_top - 0.02*yview;

	ytic_top = box_xl;
	if (box_xl<0.0)
		ytic_top = 0.0;
	if (box_xh<0.0)
		ytic_top = box_xh;
	ytic_bot = ytic_top - 0.015*xview;

	/* x axis */
	n_tics = choose_tics(box_xl, box_xh - box_xl, tics);
	glBegin(GL_LINES);
	glVertex2d(box_xl,xtic_top);
	glVertex2d(box_xh,xtic_top);
	for (i=0;i<n_tics;++i) {
		glVertex2d(tics[i], xtic_bot);
		glVertex2d(tics[i], xtic_top);
	}
	glEnd();
	for (i=0;i<n_tics;++i) {
		sprintf(label,"%.8g",tics[i]);
		ftexto((float)(tics[i] - 4*strlen(label)*x_pix_size),\
			(float)(xtic_bot - 12*y_pix_size), label);
	}


  if (dispmode!=2) {
	/* y axis */
	n_tics = choose_tics(box_yl, box_yh - box_yl, tics);
	glBegin(GL_LINES);
	glVertex2d(ytic_top,box_yl);
	glVertex2d(ytic_top,box_yh);
	for (i=0;i<n_tics;++i) {
		glVertex2d(ytic_bot, tics[i]);
		glVertex2d(ytic_top, tics[i]);
	}
	glEnd();
	for (i=0;i<n_tics;++i) {
		sprintf(label,"%.8g",tics[i]);
		ftexto((float)(ytic_bot - 8*strlen(label)*x_pix_size),\
			(float)(tics[i] - 4*y_pix_size), label);
	}
  }

}



void pixel_coord_proj()
/* sets up proj and viewing for 2d "pixel coords" in viewport */
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, width, 0.0, height);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	/* predictable rasterization tip from OpenGL Prog Guide 2nd Ed p.601... */
	glTranslatef(0.375,0.375,0.0);
}


void screeninfo(double *evals, int iv, char *bottext)
/* text info to graphics window.
 */
{
	char toptext[100],hairtext[100],*dragstr,*specstr,*spinstr;

	pixel_coord_proj();

	/* green top text */
	glColor3f(0.0,1.0,0.0);

	dragstr = dragging ? "drag " : "";
	specstr = spec_dragging ? "spec " : "";
	spinstr = freely_spinning ? "spin " : "";

	sprintf(toptext,\
		" state %3d of %d: val = %-15.12g zsc = %-9.4g  disp:%d %s%s%s",\
		iv,n_e,evals[iv],(float)zscale,dispmode,dragstr,specstr,spinstr);
	texto(0, height-18, toptext); 

	/* pink bottom text */
	glColor3f(1.0,0.5,1.0);
	texto(0, 9, bottext);
	if (ndims==1) {
		sprintf(hairtext,"+ at (%.8g,%.8g)",xof,yof);
		texto(width-8*35, 9, hairtext);
	}
}



#define SPECT_PIXEL_WIDTH	40

void plotspectrum(double *evals, int lo, int hi, int x)
/* outputs vertical display of energy spectrum: view given by eview, erange.
 * state x has marker. Uses SMOOTH (anti-aliased) lines and polygon.
 */
{
	int i;
	float range,a;

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(1.0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	/* 2d ortho proj view rectangle, origin on window right side midpoint.... */
	gluOrtho2D(-width/(float)SPECT_PIXEL_WIDTH, 0.0, -0.5, 0.5);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1.0, 1/erange, 1.0);
	glTranslatef(-1.5, -eview, 0.0);
	glBegin(GL_LINES);
	glColor3f(0.0, 1.0, 0.0);
	for (i=lo; i<=hi; ++i) {
		glVertex2f(0.0, a = (float)evals[i]);
		glVertex2f(1.0, a);
	}
	glEnd();
	glDisable(GL_LINE_SMOOTH);

	/* display triangle marker... */
	glEnable(GL_POLYGON_SMOOTH);
	glBegin(GL_TRIANGLES);
	glColor3f(0.8, 1.0, 0.8);
	glVertex2f(0.7, a = (float)evals[x]);
	glVertex2f(1.2, a + erange/80);
	glVertex2f(1.2, a - erange/80);
	glEnd();
	glDisable(GL_POLYGON_SMOOTH);

	glDisable(GL_BLEND);
}




void increment(void)
{
/* increment view params which can freely drift.
 * In our case, eye angles (and keep them modulo 2*PI)... */
	th += thdot;
	th = (th > 0.0) ? ((th < 2*PI) ? th : th-2*PI) : th+2*PI;
	ph += phdot;
	ph = (ph > 0.0) ? ((ph < 2*PI) ? ph : ph-2*PI) : ph+2*PI;
	xof2d += xof2ddot;
	yof2d += yof2ddot;
	xof += xofdot;
	yof += yofdot;
}


void animate(void)
{
	if (freely_spinning)
	/* only increment view if mouse_motion isn't already doing it (drag) */
		increment();
	else
	/* don't retrigger animation loop if flag unset... */
		glutIdleFunc(NULL);

/* if no motion, stop 'freely spinning' */
	if (thdot == 0.0 && phdot == 0.0 && xof2ddot == 0 && yof2ddot == 0 \
			&& xofdot == 0.0 && yofdot == 0.0)
		freely_spinning = 0;

	glutPostRedisplay();
}


void vis(int visible)
{
  if (visible == GLUT_VISIBLE) {
    if (freely_spinning)
      glutIdleFunc(animate);
  } else {
    if (freely_spinning)
      glutIdleFunc(NULL);
  }
}





/* ----------------------------- RENDER 1D ARRAY ------------------------------
 */
void render_1darray(int lowest, int highest)

/* 1D graph display. Alex Barnett 99/10/20
 * In multiple modes, displays only state indices from lowest to highest.			
 */
{
  int i,k;
  float *dispvec,color_frac, separate_step;
  
  /* point to the 1->n double array you want displayed... */
  dispvec = vec[iv];
  
  /* use logtoggle to control whether anti-aliased or not... */
  if (logtoggle==2 || logtoggle==4) {
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  } else {
    glDisable(GL_BLEND);
    glClear(GL_COLOR_BUFFER_BIT);
  }
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  /* 2d ortho proj view rectangle.... */
  gluOrtho2D(xof - xview/2.0, xof + xview/2.0,\
	     yof - yview/2.0, yof + yview/2.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  axes2d();
  
  if (dispmode==1) {
    /* use logtoggle to control whether lines or points... */
    if (logtoggle==3 || logtoggle==4) {
      glBegin(GL_POINTS);
      glPointSize(2.0);
    } else
      glBegin(GL_LINE_STRIP);
    glColor3f(1.0, 1.0, 0.0);
    for (i=1;i<=n;++i)
      glVertex2f((float)(xoffset + (i-1)*xlat), dispvec[i]);
    glEnd();
  } else if (dispmode==2) {
    /* multi-plot with y-separation based on current window size */
    separate_step = yview/(2+highest-lowest);
    glTranslatef(0.0, (float)(yof - yview/2.0), 0.0);
    for(k=lowest;k<=highest;++k) {
      glTranslatef(0.0, separate_step, 0.0);
      /* use logtoggle to control whether lines or points... */
      if (logtoggle==3 || logtoggle==4) {
	glBegin(GL_POINTS);
	glPointSize(1.0);
      } else
	glBegin(GL_LINE_STRIP);
      color_frac = (k-lowest)/(float)(highest-lowest);
      if (k==iv)
	glColor3f(1.0, 1.0, 1.0);
      else
	glColor3f(1.0, 1.0 - color_frac, 0.0 + color_frac);
      for (i=1;i<=n;++i)
	glVertex2f((float)(xoffset + (i-1)*xlat), vec[k][i]);
      glEnd();
    }
  }
  
  glDisable(GL_LINE_SMOOTH);
}




/* ----------------------------- RENDER 2D ARRAY ------------------------------
 */
#define TINY	1.0E-30

void plot_image(float *dispvec, int num_pixels)
/* dumps the vector as a 2d image plot, in integer x,y coords.
   Barnett 99/11/4
 */
{
    int i,j;
    float a,b;

    glPointSize(num_pixels*1.0);
    glBegin(GL_POINTS);
    if (logtoggle==1) {
	/* show amplitude - yellow as +ve, red as -ve. */
        for (i=1; i<=nx; ++i)
	    for (j=1; j<=ny; ++j) {
                a=fabsf(zscale*(b=dispvec[nx*(j-1) + i]));
		glColor3f(a,(b>0.0)?a:0.0,0.0);
		if (isinf(b))
		  glColor3f(0.0,1.0,0.0);
                glVertex2i(num_pixels*i,num_pixels*j);
            }
    } else if (logtoggle==2) {
      /* show amplitude: grey = zero, white = +ve, black = -ve. (cyan tint) */
        for (i=1; i<=nx; ++i)
	    for (j=1; j<=ny; ++j) {
                a= 0.5 + zscale*dispvec[nx*(j-1) + i];
                glColor3f(0.0,a,a);
		if (isinf(a))
		  glColor3f(1.0,0.0,0.0);
                glVertex2i(num_pixels*i,num_pixels*j);
            }
    } else if (logtoggle==3) {
	/* show intensity = square of amplitude. */
        for (i=1; i<=nx; ++i)
	    for (j=1; j<=ny; ++j) {
                a = zscale*(b=dispvec[nx*(j-1) + i])*b;
                glColor3f(a,a,0.0);
		if (isinf(b))
		  glColor3f(0.0,1.0,0.0);
                glVertex2i(num_pixels*i,num_pixels*j);
            }
    } else {
	/* show log of absval of amplitude, scaled after logging */
        for (i=1; i<=nx; ++i)
	    for (j=1; j<=ny; ++j) {
                a = 0.5 + zscale*logf(TINY + fabsf(dispvec[nx*(j-1) + i]));
                glColor3f(a,a,a);
		if (isinf(a))
		  glColor3f(0.0,1.0,0.0);
                glVertex2i(num_pixels*i,num_pixels*j);
            }
    }
    glEnd();
}


void render_2darray(int lowest, int highest)

/* 2D Lattice display. NOT VERY EFFICIENT FOR 2D IMAGE PLOTS!
 * In multiple modes, displays only state indices from lowest to highest.			
 */
{
  float a,b;
  float *dispvec;
  int i,j,k,nacross,m,num_pixels,obox,oxbox,oybox;
  
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  /* PROJECTION: 3D if appropriate otherwise 2D proj setup... */
  if (dispmode==1) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(atan(iris/2)*360.0/PI, width/(float)height, 0.05, 20.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0,0.0,-dist);
    glRotatef(-rh*180.0/PI, 0.0,1.0,0.0);
    glRotatef(-ph*180.0/PI, 1.0,0.0,0.0);
    glRotatef(-th*180.0/PI, 0.0,0.0,1.0);
    /* plot axes with true zscale but of size 1 in x,y plane... */
    axes3d();
    /* scale to physical lattice size in x,y plane, physical zscale... */
    glScalef((float)xlat, (float)xlat, zscale);
    /* place origin in centre... */
    glTranslatef(-(nx+1)/2.0,-(ny+1)/2.0,0.0);
  } else
    pixel_coord_proj();
  
  /* point to the 1->n double array you want displayed... */
  dispvec = vec[iv];

  num_pixels = max(1,(int)dist - 1); /* pixel size in 2d image */
  
  /* display current evec, or pot on space lattice 1->nx, 1->ny */
  
  if (dispmode==1) {
    /* 3d points plot...111111111111111111111111111111111111111111111111111111
     */
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_POINTS);
    if (logtoggle==1 || logtoggle==2) {
      for (i=1; i<=nx; ++i)
	for (j=1; j<=ny; ++j)
	  glVertex3f((float)i,(float)j,dispvec[nx*(j-1) + i]);
    } else if (logtoggle==3) {
      for (i=1; i<=nx; ++i)
	for (j=1; j<=ny; ++j)
	  glVertex3f((float)i,(float)j,(a=dispvec[nx*(j-1) + i])*a);
    } else {
      for (i=1; i<=nx; ++i)
	for (j=1; j<=ny; ++j)
	  glVertex3f((float)i,(float)j,\
		     logf(TINY + fabsf(dispvec[nx*(j-1) + i])));
    }
    glEnd();
    
  } else if (dispmode==2) {
    /* 2D pixel plot ...2222222222222222222222222222222222222222222222222222
     */
    glTranslatef((float)xof2d, (float)yof2d, 0.0);
	glRotatef(img_ang, 0.0, 0.0, 1.0);
    plot_image(dispvec,num_pixels);
    /* put box around the one state... */
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINE_LOOP);
    obox = (num_pixels-1)/2;
    oxbox = nx*num_pixels + (num_pixels-1)/2 + 1;
    oybox = ny*num_pixels + (num_pixels-1)/2 + 1;
    glVertex2i(obox,obox);
    glVertex2i(oxbox,obox);
    glVertex2i(oxbox,oybox);
    glVertex2i(obox,oybox);
    glEnd();

    
  } else if (dispmode==3) {
    /* multiple 2D pixel plots ...3333333333333333333333333333333333333333333
     * calcs the # states which fit horizontally across half-screen... */
    nacross = max(1,(int)((float)(width-90)/(nx+8)));
    /* loops through all states... */
    for (k=lowest; k<=highest;++k) {
      m = (int)((float)(k-lowest)/nacross);
      glPushMatrix();
      glTranslatef((float)xof2d + (k-lowest-m*nacross)*(nx+8), \
		   (float)yof2d + m*(ny+16), 0.0);
	glRotatef(img_ang, 0.0, 0.0, 1.0);
      plot_image(vec[k],1);
      if (k==iv) {
	/* put box around current state... */
	glColor3f(1.0,1.0,1.0);
	glBegin(GL_LINE_LOOP);
	glVertex2i(0,0);
        glVertex2i(nx+1,0);
        glVertex2i(nx+1,ny+1);
        glVertex2i(0,ny+1);
        glEnd();
      }
      glPopMatrix();
    }
  }
  
}


int clipstate(int i)
     /* returns a valid state number... */
{
  return max(min(i,n_e),1);
}



void keyfunc( unsigned char key, int x, int y)
/* processes keypresses in window */
{
    switch (key) {
	case 'U':
	case 'u':
		iv = clipstate(iv + 1);
		glutPostRedisplay();
		break;
	case 'D':
	case 'd':
		iv = clipstate(iv - 1);
		glutPostRedisplay();
		break;
	case 'R':
	case 'r':
		if (loadstates() == 0) {
			fprintf(stderr,"viewer: reload %s failed.\n",filename);
			exit(0);
		}
		glutPostRedisplay();
		break;
	case 'I':
	case 'i':
		infodump();
		break;
	case 'S':
	case 's':
	  img_ang += 90.0;
		glutPostRedisplay();
	  break;
	case '1':
		ratio();
		break;
	case '2':
		unratio();
		break;
	case '3':
		subtract();
		break;
	case '4':
		unsubtract();
		break;
	case ' ':
		++dispmode;
		if (dispmode>NUMBER_DISPLAY_MODES) {
			dispmode = 1;
		}
    	/* stop spinning + zero the angular velocity... */
    	freely_spinning = 0;
		thdot = phdot = 0.0;
		xof2ddot = yof2ddot = 0;
		xofdot = yofdot = 0;
		glutPostRedisplay();
		break;
	case 'L':
	case 'l':
		++logtoggle;
		if (logtoggle>NUMBER_LOG_MODES) {
			logtoggle = 1;
		}
		glutPostRedisplay();
		break;
	case 'C':
	case 'c':
		++contourtoggle;
		if (contourtoggle>=2)
			contourtoggle = 0;
		glutPostRedisplay();
		break;
    case 'Q':
    case 'q':
      exit(0);
      break;
    }
}


void mouse_motion(int mx,int my)
     /* deals with mouse motion: mx my are relative to window origin */
{
  
  /* change origin from lower left (OpenGL, IrisGL) to upper left (GLUT) */
  my = height - my;
  
  if (dragging) {
    if (dispmode==1 && ndims==2) {
      /* 3d spin */
      thdot = -(mx-omx)*5/(float)width;
      phdot = (my-omy)*5/(float)width;
    } else if (ndims==2) {
      /* 2d image slide */
      xof2ddot = mx-omx;
      yof2ddot = my-omy;
    } else {
      /* ndims=1, 2d graph slide */
      xofdot = -xview*(mx-omx)/(float)width;
      yofdot = -yview*(my-omy)/(float)height;
    }
    increment();
    /* trigger a single refresh cycle... */
    glutIdleFunc(animate);
  }
  if (scaling) {
    if (ndims==2) {
      dist += (mx-omx)*2/(float)width;
      dist = (dist > -20.0) ? ((dist < 100.0) ? dist : 100.0) : -20.0;
      zscale *= exp((my-omy)*3/(float)width);
    } else {
      xview /= exp((mx-omx)*3/(float)width);
      yview /= exp((my-omy)*3/(float)height);
      /* tie zscale to vertical view... */
      zscale = 1.0/yview;
    }
    /* trigger a single refresh cycle... */
    glutIdleFunc(animate);
  }
  if (spec_dragging) {
    eview -= (my-omy)*erange/(double)height;
    erange *= exp((mx-omx)*3/(float)width);
    /* trigger a single refresh cycle... */
    glutIdleFunc(animate);
  }
  
  /* save test result of whether still 'moving' for mouse_button... */
  mouse_still_moving = (abs(mx-omx)>1 || abs(my-omy)>1);
  
  /* restart mouse relative counting here... */
  omx = mx;
  omy = my;
  
  /* printf("mousemotion: mx = %d,%d\n",mx,my); */
}


void mouse_button(int button, int state, int x, int y)
{
	y = height - y;

	if (button == GLUT_LEFT_BUTTON) {
		if (state == GLUT_DOWN) {
    		freely_spinning = 0;
    		
		/* Shift-Left case... */
			if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
				spec_dragging = 1;
			else
    			dragging = 1;
    		
			omx = x;
			omy = y;
		} else if (state == GLUT_UP) {
   			dragging = 0;
			spec_dragging = 0;
    		if (mouse_still_moving) {
    		/* let freely spin if non-zero release velocity of mouse... */
    			freely_spinning = 1;
    		} else {
    		/* zero the velocities... */
				thdot = phdot = 0.0;
				xof2ddot = yof2ddot = 0;
 				xofdot = yofdot = 0.0;
   		}
		}
	}
	if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {
			scaling = 1;
			omx = x;
			omy = y;
		} else if (state == GLUT_UP) {
			scaling = 0;
		}
	}

	if (button == GLUT_LEFT_BUTTON || button == GLUT_MIDDLE_BUTTON)
		glutPostRedisplay();
		/*glutIdleFunc(animate);*/

/* choose cursor appearance... */
	if (dragging || scaling)
		glutSetCursor(GLUT_CURSOR_INFO);
	else if (spec_dragging)
		glutSetCursor(GLUT_CURSOR_UP_DOWN);
	else
		glutSetCursor(GLUT_CURSOR_CROSSHAIR);
}


void reshape(int w, int h)
{
  glViewport(0, 0, w, h);
  width = w;
  height = h;
}




void right_menu(int value)
/* effects the right button menu options have: --------------------------- */
{
  switch (value) {
  case 1:           /* reset view */
    initialize_view();
    glutPostRedisplay();
    break;
  case 2:
    fprintf(stderr,"\nVIEWER: multiple 1d/2d data array viewer\n");
    fprintf(stderr,"        x86/Linux version, Alex Barnett 2/26/04\n\n");
    fprintf(stderr,"See /usr/local/src/viewer_barnett/README\n\n");
    fprintf(stderr,"    Keystrokes accepted in viewer window:\n\n");
    fprintf(stderr,"U - sweep up states\n");
    fprintf(stderr,"D - sweep down states\n");
    fprintf(stderr,"R - reload data (same filename; size can change)\n");
    fprintf(stderr,"<space> - cycle multiple/single display mode\n");
    fprintf(stderr,"L - cycle 2D display op: value (+yellow/-red)\n");
    fprintf(stderr,"                         value (\"blue\"scale)\n");
    fprintf(stderr,"                         value^2 (\"yellow\"scale)\n");
    fprintf(stderr,"                         log(abs(value)) (greyscale)\n");
    fprintf(stderr,"    OR in 1D, cycle points/lines/aliasing mode\n");
    fprintf(stderr,"S - rotate 2d image by pi/2\n");
    fprintf(stderr,"Left-mouse - move viewpoint (drag/spin)\n");
    fprintf(stderr,"Shift-Left-mouse - adjust spectrum view\n");
    fprintf(stderr,"Middle-mouse - x,y motion adjusts zoom (1D)\n");
    fprintf(stderr,"               x adj distance/zoom, y adj z-scale (2D)\n");
    fprintf(stderr,"Right-mouse - bring up menu\n");
    fprintf(stderr,"I - dump some data properties (currently spectrum)\n\n");
    fprintf(stderr,"1 (2) - do (undo) certain data-changing operations\n");
    fprintf(stderr,"3 (4) - do (undo) certain data-changing operations\n");
    fprintf(stderr,"        See code: keyfunc(), modify for your needs!\n\n");
    break;
  case 3:
    exit(0);
    break;
  case 4:       /* postscript out */
    dump_ps(iv, iv, 1, 0.1, 1, 1, 0);
    break;
  case 5:       /* multiple postscript out */
    dump_ps(1, 16, 4, 0.1, 1, 1, 1);
    break;
  }
}



void display_routine(void)
/* GLUT calls this when redisplay required. */
{
	char userinfo[100];

	/* show data as 1d graphs or 2d lattices ... */
	if (ndims==1)
		render_1darray(1,n_e);
	else
		render_2darray(1,n_e);

	/* could put a calculation routine here. */

	/* show eigenvalue spectrum... */
	if (disp_spec)
	  plotspectrum(val, 1, n_e, iv);

	/* info to appear near bottom of window :  */
	sprintf(userinfo," param1 = %g  param2 = %g",param1,param2);
	/* display params and info... */
	screeninfo(val, iv, userinfo);

	/* frame complete: double buffer swap... */
	glutSwapBuffers();
}




/* ========================================= MAIN =========================== */
int main(int argc, char **argv)
{
   
    
    /* string workspace... */
    char *opts, str[100];

    
    if (argc<2) {
        printf("viewer: multiple 1d/2d data array viewer\n");
        printf("        Alex Barnett 4/28/02\n\n");
        printf("usage:	viewer infile.sta [opts] [-GLUToptions]\n\n");
        printf("if opts includes: e - gives endian flipping of binary reads\n");
        printf("                  s - remove spectrum view\n\n");
        return 0;
    } else {
        printf("\nviewer (v.8/21/03 Alex Barnett)\n");
        printf("Drag right-mouse-button menu to get Help...\n\n");
    }

    /* cmd line defaults */
    ENDIAN_FLIP = 0;    /* 0 for DEC OSF1, alpha, i386... */
    disp_spec = 1;

/* read filename... */
    sscanf(argv[1],"%s",filename);

    if (argc>2) {
      if (argv[2][0]!='-') {
	/* no hyphen means viewer opts are given... */
	opts = argv[2];
	if (strchr(opts,'s')!=NULL)
	  disp_spec = 0;
	/* flip endian (for SGI) and remove -e from arg list... */
	if (strchr(opts,'e')!=NULL)
	  ENDIAN_FLIP = 1;
	--argc;
	++argv;
      }
    }

/* alter calling arguments to remove one of them, for passing to GLUT... */
    --argc;
    ++argv;
    
/* initial loading fixes ndims for ever... */
    ndims = 0;
    

/* load data... */
    if (loadstates() == 0) {
        fprintf(stderr,"viewer: loadstates failed.\n");
        exit(0);
    }
    
    
/* GLUT setup... */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    if (ndims==2)
        glutInitWindowSize(min(1000,nx+100), min(800,ny+60));
    else
        glutInitWindowSize(500,400);
    
    if (n_e==1)
      sprintf(str,"%s (viewer: one %dd array)",filename,ndims);
    else
      sprintf(str,"%s (viewer: %d %dd arrays)",filename,n_e,ndims);


    glutCreateWindow(str);
    
    width  = glutGet(GLUT_WINDOW_WIDTH);
    height = glutGet(GLUT_WINDOW_HEIGHT);
    xpos   = glutGet(GLUT_WINDOW_X);
    ypos   = glutGet(GLUT_WINDOW_Y);
    
/* tell GLUT which interrupt event routines to use: */
    glutMotionFunc(mouse_motion);
    glutMouseFunc(mouse_button);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyfunc);
    glutDisplayFunc(display_routine);

    
/* desired menu options: see right_menu for effects... */
    glutCreateMenu(right_menu);
    glutAddMenuEntry("reset view", 1);
    glutAddMenuEntry("output postscript (current sta)", 4);
    glutAddMenuEntry("output postscript (16 sta, as 4 by 4, val labels)", 5);
    glutAddMenuEntry("help", 2);
    glutAddMenuEntry("quit", 3);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    initialize_view();
    dispmode = 2;
    
/* put lighting model here, in ndims=2 case ? */

    
/* pass all event-handling to GLUT... */
    glutMainLoop();
    return 0;             /* ANSI C requires main to return int. */
}

