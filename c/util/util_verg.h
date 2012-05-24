#ifndef _UTIL_VERG_H_
#define _UTIL_VERG_H_

#include "../../vergini/billiard.h" // for Billiard
#include "bit_array.h"

bit_array_t *createScaledMaskFromBilliard(Billiard *b, double xl, double xh, double yl, double yh, double dx, double upsample_ratio, double scale, int ny, int nx);
bit_array_t *createMaskFromFile(char *file, int *ny, int *nx);

#endif
