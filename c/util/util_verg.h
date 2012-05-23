#ifndef _UTIL_VERG_H_
#define _UTIL_VERG_H_

#include "../../vergini/billiard.h" // for Billiard

int **createScaledMaskFromBilliard(Billiard *b, double xl, double xh, double yl, double yh, double dx, double upsample_ratio, double scale, int ny, int nx);
int **createMaskFromFile(char *file, int *ny, int *nx);

#endif
