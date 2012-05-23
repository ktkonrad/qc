#ifndef _UTIL_H_
#define _UTIL_H_

#define ERROR(fmt, args...) if(verb) { fprintf(stderr, "Error: %s: %s: %d: "fmt"\n", __FILE__, __FUNCTION__,  __LINE__, ## args); }
#define MALLOC_CHECK(ptr) if(!ptr) { ERROR("malloc failed!"); exit(OUT_OF_MEMORY_ERR); }
#define SET(loc, val) do {loc = (char *)malloc((strlen(val)+1)*sizeof(char)); MALLOC_CHECK(loc); strcpy(loc, val);} while (0)
#define RESET(loc, val) do {loc = (char *)realloc(loc, (strlen(val)+1)*sizeof(char)); MALLOC_CHECK(loc); strcpy(loc, val);} while (0)

int **imatrix(int rows, int cols);
void free_imatrix(int **m);
char **cmatrix(int rows, int cols);
void free_cmatrix(char **m);
double **dmatrix(int ny, int nx);
void free_dmatrix(double **grid);

#endif
