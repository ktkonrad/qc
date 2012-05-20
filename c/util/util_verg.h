#ifndef _UTIL_VERG_H_
#define _UTIL_VERG_H_

#include "../../vergini/billiard.h" // for Billiard

int **createScaledMaskFromBilliard(Billiard b, double dx, double scale, int *ny, int *nx);
int **createMaskFromFile(char *file, int *ny, int *nx);

#endif
