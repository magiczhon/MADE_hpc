#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include "compute_pi.h"

#define nThreads 8  // number of threads to use

unsigned int seeds[nThreads];
int MAX_INT = (1 << 31) - 1;


void seedThreads() {
    int my_thread_id;
    unsigned int seed;
    #pragma omp parallel private (seed, my_thread_id)
    {
        my_thread_id = omp_get_thread_num();
        
        //create seed on thread using current time
        
        unsigned int seed = (unsigned) time(NULL);
        //munge the seed using our thread number so that each thread has its
        //own unique seed, therefore ensuring it will generate a different set of numbers
        seeds[my_thread_id] = (seed & 0xFFFFFFF0) | (my_thread_id + 1);
        
        printf("Thread %d has seed %u\n", my_thread_id, seeds[my_thread_id]);
    }
    
}


bool circle(float x, float y, float r) {
    return x*x + y*y <= r*r;
}


float parallel_compute_pi(int *dots_in_circle, const int all_dots) {

    int dots_in_c = 0;
    int i;
    int tid = 0;
    float x0, y0, seed;

    omp_set_num_threads(nThreads);  
    seedThreads();

    #pragma omp parallel private(x0, y0, i, tid, seed) shared(dots_in_c, all_dots) 
    {
    tid = omp_get_thread_num();   // my thread id
    seed = seeds[tid];            // it is much faster to keep a private copy of our seed
    srand(seed);   

    #pragma omp for reduction(+:dots_in_c)
    for (i = 0; i < all_dots; ++i) {
        x0 = rand_r(&seed) / (float)MAX_INT;
        y0 = rand_r(&seed) / (float)MAX_INT;
        if (circle(x0, y0, 1))
            dots_in_c += 1;
    }
    }
    *dots_in_circle = dots_in_c;
    return *dots_in_circle / (all_dots + 1.0) * 4;
}


float compute_pi(int *dots_in_circle, const int all_dots) {

    int dots_in_c = 0;
    int i;
    int tid = 0;
    float x0, y0, seed;

    srand(time(NULL));   

    for (i = 0; i < all_dots; ++i) {
        x0 = rand_r(&seed) / (float)MAX_INT;
        y0 = rand_r(&seed) / (float)MAX_INT;
        if (circle(x0, y0, 1))
            dots_in_c += 1;
    }
    *dots_in_circle = dots_in_c;
    return *dots_in_circle / (all_dots + 1.0) * 4;
}


