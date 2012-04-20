#ifndef ALEX_USEFUL_H
#define ALEX_USEFUL_H 1

/* NUMBER CONSTANTS */

/* 2^31 - 1 written as a float */
#define MAXINT		2147483647.0
#define PI		3.14159265358979323846

/* MACROS */

#define max(a,b)	((a) > (b) ? (a) : (b))
#define min(a,b)	((a) < (b) ? (a) : (b))


#ifndef SQ
#define SQ(x) ((x)*(x))
#endif

#ifndef max
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#endif


#endif /* ALEX_USEFUL_H */
