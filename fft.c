/****************************************************************************
 * Copyright (c) Gregory F. Johnson, Gnu Public Licence v. 2.0.
 * File Name    : fft.c
 ****************************************************************************/
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "complex.h"
#define version_major     0
#define version_minor     1
#define version_increment 0

static void usage() {
    fprintf(stderr, "fft version %d.%d.%d\n", version_major, version_minor, version_increment);
    fprintf(stderr, "usage:  fft\n");
    fprintf(stderr, "            [-h]         # help (this message)\n");
    fprintf(stderr, "            [-check N]   # calculate largest prime factor of N; smaller value => faster fft\n");
    fprintf(stderr, "            [-b]         # inverse (backward) fft\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                         # default input format is real (non-complex) values\n");
    fprintf(stderr, "            [-ip]        # input complex values in polar (radius, radians)\n");
    fprintf(stderr, "            [-ic]        # input complex values in cartesian (real, imaginary)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                         # default output format is cartesian (real, imaginary) values\n");
    fprintf(stderr, "            [-op]        # output complex values in polar (radius, radians)\n");
    fprintf(stderr, "            [-oa]        # output double-precision amplitudes\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            [input_file] # default input is stdin\n");
}

static void printAmplitudes(complex *vec, int length) {
    for (int i = 0; i < length; ++i) {
        printf("%.12lf\n", cabs(vec[i]));
    }
}

int main(int argc, char **argv) {
    FILE *input = NULL;
    complex *inputvec = NULL;
    bool inputPolar = false, inputComplex = false;
    bool outputPolar = false, outputAmplitudes = false;
    int checkSize = -1;
    bool forward = true;
    int result = 0;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-b") == 0) {
            forward = false;

        } else if (strcmp(argv[i], "-check") == 0) {
            if (++i >= argc || 1 != sscanf(argv[i], "%d", &checkSize)
                            || checkSize < 0)
            {
                fprintf(stderr, "-check N\n");
                exit(1);
            }

        } else if (strcmp(argv[i], "-ip") == 0) {
            inputPolar = true;

        } else if (strcmp(argv[i], "-ic") == 0) {
            inputComplex = true;

        } else if (strcmp(argv[i], "-op") == 0) {
            outputPolar = true;

        } else if (strcmp(argv[i], "-oa") == 0) {
            outputAmplitudes = true;

        } else if (strcmp(argv[i], "-h") == 0) {
            usage();
            exit(0);

        } else {
            if (input != NULL) {
                fprintf(stderr, "only one input file allowed.\n");
                exit(1);
            }

            input = fopen(argv[i], "r");
            if (input == NULL) {
                fprintf(stderr, "failed to open input file %s\n", argv[i]);
                exit(1);
            }

        }
    }

    if (input == NULL) {
        input = stdin;
    }

    if (checkSize != -1) {
        int maxFactor = -1;
        do {
            int factor = primeFactor(checkSize);
            if (maxFactor < factor) maxFactor = factor;
            checkSize /= factor;
        } while (checkSize != 1);

        printf("%d\n", maxFactor);
        return(0);
    }

    if (inputPolar && inputComplex) {
        fprintf(stderr, "-ip and -ic are mutually exclusive\n");
        exit(1);
    }

    if (outputPolar && outputAmplitudes) {
        fprintf(stderr, "-op and -oa are mutually exclusive\n");
        exit(1);
    }

    int length;
    complex *outputvec = NULL;

    if (inputPolar) {
        length = readPolarFile(&inputvec, input);

    } else if (inputComplex) {
        length = readComplexFile(&inputvec, input);

    } else {
        length = readDoubleFile(&inputvec, input);
    }

    outputvec = (complex *) malloc(length * sizeof(complex));
    result = FFT(outputvec, inputvec, length, forward);

    if (result < 0) {
        fprintf(stderr, "invalid input\n");

    } else if (outputAmplitudes) {
        printAmplitudes(outputvec, length);

    } else {
        cprintf(outputvec, length, 1, outputPolar);
    }

    if (inputvec)  { free(inputvec); }
    if (outputvec) { free(outputvec); }

    return result;
}
