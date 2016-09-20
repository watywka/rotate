#ifndef SOLVE_H
#define SOLVE_H
#include <stdio.h>

extern int debug;
int solve(int n, double* a, double* b, double* x);
void printm(FILE* fout, int n, double* a, double* b);

#endif
