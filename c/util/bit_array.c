/*
  2d bit array

  Kyle Konrad
  4/17/2012
*/

#include "util.h"
#include "bit_array.h"

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
  return (arr->data[y][x/8] & (1 << x % 8)) >> x % 8;
}

int bit_array_set(bit_array_t *arr, int x, int y) {
  // no validation of indicies
  arr->data[y][x/8] |= 1 << x % 8;
}

int bit_array_reset(bit_array_t *arr, int x, int y) {
  // no validation indicies
  arr->data[y][x/8] &= ~(1 << x % 8);
}

void free_bit_array(bit_array_t *arr) {
  free_cmatrix(arr->data);
  free(arr);
}
