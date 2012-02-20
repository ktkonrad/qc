#ifndef _UTIL_VERG_H_
#define _UTIL_VERG_H_

#include "../../vergini/billiard.h" // for Billiard

char **createMaskFromBilliard(Billiard b, double dx, int *masky, int *maskx);
char **createScaledMaskFromBilliard(Billiard b, double dx, int *masky, int *maskx, double scale);

#endif
