#include <stdio.h>
#include <stdlib.h>
#include "SFMT/SFMT.h"

int main(int argc, char* argv[]) {
    int i, cnt, seed, r2;
    double x, y, pi;
    const int NUM = 100;
    static uint32_t r;
    sfmt_t sfmt;
    seed=10;

    // Seed random number generator 
    sfmt_init_gen_rand 	(&sfmt, seed);

    for (i=0; i<NUM; i++){
       r=sfmt_genrand_uint32(&sfmt); //init_gen_rand or init_by_array must be called before this function
        printf("%d\n", (r)%2);

    }
    return 0;
}
