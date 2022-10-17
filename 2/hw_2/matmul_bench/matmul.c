#include <stdio.h>
#include <stdlib.h>
#include "helper.c"


double CalcMatMulTime_ijk(double * A, double * B, double * C, size_t N)
{
    struct timeval start, end;
    double r_time = 0.0;
    size_t i, j, k;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);
    
    for (i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        {
            for(k = 0; k < N; k++)
                C[i * N + j] = C[i * N + j] + A[i * N + k] * B[k * N + j];
        }
    gettimeofday(&end, NULL);
    
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    
    return r_time;
}

double CalcMatMulTime_jik(double * A, double * B, double * C, size_t N)
{
    struct timeval start, end;
    double r_time = 0.0;
    size_t i, j, k;

    ZeroMatrix(&C[0], N);


    gettimeofday(&start, NULL);
    


    for (j = 0; j < N; j++)
        for(i = 0; i < N; i++)
        {
            for(k = 0; k < N; k++)
                C[i * N + j] = C[i * N + j] + A[i * N + k] * B[k * N + j];
        }
    gettimeofday(&end, NULL);
    
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    
    return r_time;
}

double CalcMatMulTime_kij(double * A, double * B, double * C, size_t N)
{
    struct timeval start, end;
    double r_time = 0.0;
    size_t i, j, k;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);
    
    for (k = 0; k < N; k++)
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
                C[i * N + j] = C[i * N + j] + A[i * N + k] * B[k * N + j];
        }
    gettimeofday(&end, NULL);
    
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    
    return r_time;
}

double CalcMatMulTime_kij_opt(double * A, double * B, double * C, size_t N)
{
    struct timeval start, end;
    double r_time = 0.0;
    size_t i, j, k;

    size_t dummy = 0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);
    
    for (k = 0; k < N; k++)
        for(i = 0; i < N; i++)
        {
            dummy = i * N;
            for(j = 0; j < N; j++)
                C[dummy + j] = C[dummy + j] + A[dummy + k] * B[k * N + j];
        }
    gettimeofday(&end, NULL);
    
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    
    return r_time;
}