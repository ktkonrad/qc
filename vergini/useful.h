#ifndef ALEX_USEFUL_H
#define ALEX_USEFUL_H 1

/* NUMBER CONSTANTS */

/* 2^31 - 1 written as a float */
#define MAXINT		2147483647.0
#define PI		3.14159265358979323846

/* MACROS */

#define max(a,b)	((a) > (b) ? (a) : (b))
#define min(a,b)	((a) < (b) ? (a) : (b))


/* default string lengths */
#define LEN 256


/* max sensible number of doubles to allocate */
#define MAX_N_DOUBLES 64000000

#ifndef SQ
#define SQ(x) ((x)*(x))
#endif

#ifndef max
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#endif

/* Useful routines */

double lower_limit(double *a, int M);
double upper_limit(double *a, int M);
float flower_limit(float *a, int M);
float fupper_limit(float *a, int M);
void endian_array_4byte(char *array, int n);
int endian_4byte(void *in);
double endian_8byte(double *in);
int round_and_clip(double x, char *var_name);
int separated_numbers(char *s, double *a, int n);
int grab_cmdline(char *cmdline, int argc, char **argv);
double inverse_interp(int M, double *xs, double *ys, double y);

#endif /* ALEX_USEFUL_H */
