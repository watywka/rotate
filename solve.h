#ifndef SOLVE_H
#define SOLVE_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

enum FUNC
{
    symm,
    positive_symm,
    hilbert,
    errf
};
extern int debug;
enum FUNC formula(char* str);
void fill(enum FUNC xin,int n, double* a);
int solve(int n, double* a, double* b, double* x);
void printm(FILE* fout, int n, double* a, double* b);

#endif
