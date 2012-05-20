/*
  2d bit array

  Kyle Konrad
  4/17/2012
*/

#include "bit_array.h"
#include "util.h"
#include "exit_codes.h"
#include <stdlib.h>
#include <stdio.h>

extern int verb;

bit_array_t *new_bit_array(int ny, int nx) {
  // create an ny x nx bit array
  bit_array_t *arr;

  arr = (bit_array_t *)malloc(sizeof(bit_array_t));
  MALLOC_CHECK(arr);

  arr->data = cmatrix(ny, (nx+7)/8);
  if (!arr->data) {
    ERROR("Failed to allocate bit array data memory");
    free(arr);
    return NULL;
  }

  arr->ny = ny;
  arr->nx = nx;

  return arr;
}

int bit_array_get(bit_array_t *arr, int x, int y) {
  return (arr->data[y][x/8] & (1 << (x & 7))) > 0 ;
}

void bit_array_set(bit_array_t *arr, int x, int y) {
  // no validation of indicies
  arr->data[y][x/8] |= 1 << (x & 7); // x & 7 == x % 8
}

void bit_array_reset(bit_array_t *arr, int x, int y) {
  // no validation indicies
  arr->data[y][x/8] &= ~(1 << (x & 7));
}

void free_bit_array(bit_array_t *arr) {
  free_cmatrix(arr->data);
  free(arr);
}

/*
  write a bit array to a file, tab separator

  input:
         array - array to write
	 file  - name of file to write to
  output: status code
*/
int bit_array2file(bit_array_t *array, char *file) {
  int i, j;
  FILE *out = fopen(file, "w");
  if (out == NULL) {
    ERROR("failed to open %s", file);
    return IO_ERR;
  }
  for (i = 0 ; i < array->ny ; i++) {
    for (j = 0 ; j < array->nx ; j++) {
      if (j != 0)
	fprintf(out, "\t");
      fprintf(out, "%d", bit_array_get(array, j, i));
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}
