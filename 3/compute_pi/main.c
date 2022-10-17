#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include "compute_pi.h"

int main(int argc, char *argv[]) {
    int const all_dots = 100000000;
    int dots_in_circle = 0;
    float PI = 0.0;
    struct timeval start, end;
    double r_time = 0.0;

    gettimeofday(&start, NULL);
    PI = compute_pi(&dots_in_circle, all_dots);
    gettimeofday(&end, NULL);
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    printf("Dots in circle: %d\nAll dots: %d\nPI: %f\n", dots_in_circle, all_dots, PI);
    printf("Time not parallel: %f\n\n", r_time);
    
    gettimeofday(&start, NULL);
    PI = parallel_compute_pi(&dots_in_circle, all_dots);
    gettimeofday(&end, NULL);
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    printf("Dots in circle: %d\nAll dots: %d\nPI: %f\n", dots_in_circle, all_dots, PI);
    printf("Time parallel    : %f\n\n", r_time);
    return 0;
}