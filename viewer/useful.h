/* NUMBER CONSTANTS */

/* 2^31 - 1 written as a float */
#define MAXINT		2147483647.0
#define PI		3.14159265358979323846

/* FUNCTION CONSTANTS */

#define max(a,b)	((a) > (b) ? (a) : (b))
#define min(a,b)	((a) < (b) ? (a) : (b))


/* Useful routines */

double lower_limit(double *a, int M);
double upper_limit(double *a, int M);
float flower_limit(float *a, int M);
float fupper_limit(float *a, int M);
void endian_array_4byte(char *array, int n);
int endian_4byte(void *in);
double endian_8byte(double *in);
