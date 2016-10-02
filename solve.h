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
	upper,
	disg,
    errf
};
extern int debug;
enum FUNC formula(char* str);
void fill(enum FUNC xin,int n, long double* a);
int solve(int n, long double* a, long double* b, long double* x);
void printm(FILE* fout, int n, long double* a, long double* b);

#endif
