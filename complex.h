/****************************************************************************
 * Copyright (c) Gregory F. Johnson, Gnu Public Licence v. 2.0.
 * File Name    : complex.h
 ****************************************************************************/
#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdbool.h>
#include <complex.h>

void cprintf(complex *array, int rows, int cols, bool outputPolar);
int readComplexFile(complex **inputvec, FILE *input);
int readDoubleFile(complex **inputvec, FILE *input);
int readPolarFile(complex **inputvec, FILE *input);
int primeFactor(int n);

int FFT(complex *outvec, complex *invec, unsigned int n, bool forward);

#ifndef M_PIl
    #define M_PIl		3.1415926535897932384626433832795029L
#endif

#endif
