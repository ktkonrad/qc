#ifndef _COUNT_NODAL_DOMAINS_H_
#define _COUNT_NODAL_DOMAINS_H_

int countNodalDomains(double **grid, char **mask, int ny, int nx);

int findNextUnseen(int **counted, int *i, int *j, int ny, int nx);

//void findDomainRecursive(int i, int j);

void findDomain(double **grid, int **counted, int i, int j, int nd, int ny, int nx);

#endif
