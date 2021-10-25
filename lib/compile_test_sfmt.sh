#!/bin/bash

#To compile sample1.c with SFMT.c with the period of 2^607, type
#gcc -O3 -DSFMT_MEXP=607 -o test_sfmt SFMT/SFMT.c test_sfmt.c
#If your CPU supports SSE2 and you want to use optimized SFMT for SSE2, type
gcc -O3 -msse2 -DHAVE_SSE2 -DSFMT_MEXP=607 -o test_compile_sfmt SFMT/SFMT.c test_compile_sfmt.c
