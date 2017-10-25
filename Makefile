#****************************************************************************
# Copyright (c) Gregory F. Johnson, Gnu Public Licence v. 2.0.
# File Name    : Makefile
#****************************************************************************

all: fft

fft: fft.c complex.o
	gcc -std=c99 -Wall -g -o fft fft.c complex.o -lm

complex.o: complex.c
	gcc -std=c99 -Wall -g -c complex.c

complex_UT: complex.c
	gcc -std=c99 -Wall -g -DUNIT_TEST -o complex_UT complex.c -lm

clean:
	rm -f fft complex_UT *.o core.* *.stackdump
