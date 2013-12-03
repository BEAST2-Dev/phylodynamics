#ifndef __EXPOTREE_H__
#define __EXPOTREE_H__

#include <math.h>

void rExpoTree(double* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, int* parVecLen, double* times, 
    int* ttypes, double* p, double* t0, int* RSImodel, 
    int* Rvflag, int* Rrescale);

void expoTree(int n, double* times, int* ttypes, 
    double* p, int wrklen, double* wrk, int iwrklen, int* iwrk);

#endif
