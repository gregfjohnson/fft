/****************************************************************************
 * Copyright (c) Gregory F. Johnson, Gnu Public Licence v. 2.0.
 * File Name    : complex.c
 ****************************************************************************/
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include "complex.h"
#include "primes.h"

static bool verbose = false;
static void dbPrint(int chunk, int chunkSize, int i, complex directionalOmega, complex omegaPower,
                    complex invec[], complex outvec[]);

static double cangle(complex n) {
    return atan2(cimag(n), creal(n)) * 180. / M_PIl;
}

void cprintf(complex *array, int rows, int cols, bool outputPolar) {
    char buf[64];
    char **realFormats;
    char **imagFormats;
    int *maxRealLen;
    int *maxImagLen;
    double data;

    realFormats = (char **) malloc(cols * sizeof(char*));
    imagFormats = (char **) malloc(cols * sizeof(char*));
    maxRealLen  = (int *)   calloc(cols, sizeof(int));
    maxImagLen  = (int *)   calloc(cols, sizeof(int));

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int len;

            data = outputPolar ? cabs(array[cols*row + col]) : creal(array[cols*row + col]);
            sprintf(buf, "%.12lf%n", data, &len);

            if (maxRealLen[col] < len) maxRealLen[col] = len;

            data = outputPolar ? cangle(array[cols*row + col]) : cimag(array[cols*row + col]);
            sprintf(buf, "%.12lf%n", data, &len);

            if (maxImagLen[col] < len) maxImagLen[col] = len;
        }
    }

    for (int col = 0; col < cols; ++col) {
        realFormats[col] = (char *) malloc(64);
        imagFormats[col] = (char *) malloc(64);
        sprintf(realFormats[col], "%c%d.12lf", '%', maxRealLen[col]);
        sprintf(imagFormats[col], " %c%d.12lf", '%', maxImagLen[col]);
    }

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {

            data = outputPolar ? cabs(array[cols*row + col]) : creal(array[cols*row + col]);
            printf(realFormats[col], data);

            data = outputPolar ? cangle(array[cols*row + col]) : cimag(array[cols*row + col]);
            printf(imagFormats[col], data);

            if (col < cols - 1) printf("  ");
        }
        printf("\n");
    }

    for (int col = 0; col < cols; ++col) {
        free(realFormats[col]);
        free(imagFormats[col]);
    }
    free(realFormats);
    free(imagFormats);
    free(maxRealLen);
    free(maxImagLen);
}

int primeFactor(int n) {
    for (int i = 0; i < sizeof(primes)/sizeof(primes[0]); ++i) {
        if (n % primes[i] == 0) return primes[i];
    }
    return n;
}

static void recFFT(complex *outvec, complex *invec, unsigned int n, bool forward) {
    if (n == 1) {
        outvec[0] = invec[0];

    } else {
        int chunkCount = primeFactor(n);
        int chunkSize  = n / chunkCount;

        for (int i = 0; i < n; ++i) {
            int indexInChunk = i / chunkCount;
            int chunkIndex   = i % chunkCount;

            outvec[chunkIndex * chunkSize + indexInChunk] = invec[i];
        }

        for (int chunk = 0; chunk < chunkCount; ++chunk) {
            recFFT(&invec[chunk * chunkSize], &outvec[chunk * chunkSize], chunkSize, forward);
        }

        complex directionalOmega =   forward
                                   ? cexp(-2i * M_PIl / n)
                                   : cexp( 2i * M_PIl / n);

        for (int i = 0; i < n; ++i) {
            outvec[i] = invec[i % chunkSize];
        }

        for (int chunk = 1; chunk < chunkCount; ++chunk) {
            for (int i = 0; i < n; ++i) {
                complex omegaPower = cpow(directionalOmega, i * chunk);

                if (verbose) {
                    dbPrint(chunk, chunkSize, i, directionalOmega, omegaPower, invec, outvec);
                }

                outvec[i] += invec[chunk * chunkSize + i%chunkSize] * omegaPower;
            }
        }
    }
}

int FFT(complex *outvec, complex *invec, unsigned int n, bool forward) {
    complex *changeableInvec;

    if (n == 0) return -1;

    changeableInvec = (complex *) malloc(n * sizeof(complex));
    if (changeableInvec == NULL) return -1;

    memcpy(changeableInvec, invec, n * sizeof(complex));

    recFFT(outvec, invec, n, forward);
    free(changeableInvec);

    for (int i = 0; i < n; ++i) {
        outvec[i] /= sqrt(n);
    }

    return 0;
}

int readComplexFile(complex **inputvec, FILE *input) {
    int used = 0;
    int len = 4;
    double real, imag;

    *inputvec = (complex *) malloc(len * sizeof(complex));

    while (2 == fscanf(input, "%lf %lf", &real, &imag)) {
        if (used == len) {
            len *= 2;
            *inputvec = realloc(*inputvec, len * sizeof(complex));
        } 
        (*inputvec)[used++] = real + imag * 1i;
    }

    return used;
}

int readPolarFile(complex **inputvec, FILE *input) {
    int used = 0;
    int len = 4;
    double length, angle;

    *inputvec = (complex *) malloc(len * sizeof(complex));

    while (2 == fscanf(input, "%lf %lf", &length, &angle)) {
        if (used == len) {
            len *= 2;
            *inputvec = realloc(*inputvec, len * sizeof(complex));
        } 
        (*inputvec)[used++] = length * cos(angle * M_PIl / 180.)
                           + (length * sin(angle * M_PIl / 180.)) * 1i;
    }

    return used;
}

int readDoubleFile(complex **inputvec, FILE *input) {
    int used = 0;
    int len = 4;
    double data;

    *inputvec = (complex *) malloc(len * sizeof(complex));

    while (1 == fscanf(input, "%lf", &data)) {
        if (used == len) {
            len *= 2;
            *inputvec = realloc(*inputvec, len * sizeof(complex));
        } 
        (*inputvec)[used++] = data + 0i;
    }

    return used;
}

static void dbPrint(int chunk, int chunkSize, int i, complex directionalOmega, complex omegaPower,
                    complex invec[], complex outvec[])
{
    printf("%d %d <%.2f + %.2fi> <%.2f + %.2fi> "
           "|| <%.2f + %.2fi> + <%.2f + %.2fi> -> <%.2f + %.2fi>\n",
           chunk, i,
           creal(directionalOmega), cimag(directionalOmega),
           creal(omegaPower), cimag(omegaPower),
           creal(invec[chunk * chunkSize + i%chunkSize] * omegaPower),
               cimag(invec[chunk * chunkSize + i%chunkSize] * omegaPower),
           creal(outvec[i]),
               cimag(outvec[i]),
           creal(outvec[i] + invec[chunk * chunkSize + i%chunkSize] * omegaPower),
               cimag(outvec[i] + invec[chunk * chunkSize + i%chunkSize] * omegaPower));
}
#ifdef UNIT_TEST

static void testPrimes() {
    int n = 1000000;
    int pfac;
    do {
        pfac = primeFactor(n);
        printf("pfac %d\n", pfac);
        n /= pfac;
    } while (pfac != n);
    printf("pfac %d\n", n);
}

static void testprint() {
    complex vec[3] = { 1+1i, 2+1i, 3+1i };
    cprintf(vec, 3,1, false);
    printf("\n");

    complex vec6[6] = { 1+1i, 2+1i, 3+1i, -1+1i, -2+1i, 3-1i };
    cprintf(vec6, 3,2, false);
}

static void testFFTWithPrimeFactors() {
    complex  in1[1] = { 1 + 0i, };
    complex out1[1];
    int result = FFT(out1, in1, 1, true);
    printf("\ntest result 1: %d\n", result);
    cprintf(out1, 1,1, false); printf("\n");

    complex  in2[2] = { 1 + 0i, 1 + 0il };
    complex out2[2];
    result = FFT(out2, in2, 2, true);
    printf("\ntest result 2:  %d\n", result);
    cprintf(out2, 2,1, false); printf("\n");
}

static void testHugeFFT() {
    for (int len = 2; len <= 8388608; len = 2*len) {
        complex *bigVec = calloc(len, sizeof(complex)); if (bigVec == NULL) exit(1);
        complex *result = calloc(len, sizeof(complex)); if (result == NULL) exit(1);

        time_t startTime = time(NULL);
        FFT(result, bigVec, len, true);
        printf("len %d; time %ld\n", len, (int) time(NULL) - startTime);

        if (bigVec != NULL) free(bigVec);
        if (result != NULL) free(result);
    }
}

int main(int argc, char **argv) {
    testprint();
    testHugeFFT();
    testPrimes();
    testFFTWithPrimeFactors();

    return 0;
}
#endif
