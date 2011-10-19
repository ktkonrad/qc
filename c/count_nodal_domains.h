#ifndef _COUNT_NODAL_DOMAINS_H_
#define _COUNT_NODAL_DOMAINS_H_

// NOTE: trouble_count is TEMPORARY
int countNodalDomains(double **grid, char **mask, int ny, int nx, int *trouble_count);

int findNextUnseen(int **counted, int *i, int *j, int ny, int nx);

//void findDomainRecursive(int i, int j);

void findDomain(double **grid, int **counted, int i, int j, int nd, int ny, int nx);

#endif
