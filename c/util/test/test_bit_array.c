#include "../bit_array.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NY 9
#define NX 15

int verb = 0;

int main(int argc, char **argv) {
  int i, j, r;
  int int_arr[NY][NX];

  // allocate
  bit_array_t *arr = new_bit_array(NY, NX);
  assert(arr->ny == NY);
  assert(arr->nx == NX);

  // set
  for (i = 0 ; i < NY ; i++) {
    for (j = 0 ; j < NX ; j++) {
      r = rand() % 2;
      if (r) {
        bit_array_set(arr, j, i);
      } else {
        bit_array_reset(arr, j, i);
      }
      int_arr[i][j] = r;
    }
  }

  // get
  for (i = 0 ; i < NY ; i++) {
    for (j = 0 ; j < NX ; j++) {
      assert(bit_array_get(arr, j, i) == int_arr[i][j]);
    }
  }

  // free
  free_bit_array(arr);

  printf("Passed.\n");
}
