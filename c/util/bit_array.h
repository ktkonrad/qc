/*
  2d bit array

  Kyle Konrad
  4/17/2012
 */

typedef struct {
  char **data;
  int ny;
  int nx;
} bit_array_t;

bit_array_t *new_bit_array(int ny, int nx);
void free_bit_array(bit_array_t *arr);
int bit_array_get(bit_array_t *arr, int x, int y);
int bit_array_set(bit_array_t *arr, int x, int y);
int bit_array_reset(bit_array_t *arr, int x, int y);
