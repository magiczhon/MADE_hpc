#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void ZeroMatrix(double * A, size_t N)
{
    for(size_t i=0; i < N; i++)
    {
        for(size_t j=0; j<N; j++)
        {
            A[i * N + j] = 0.0;
        }
    }
}

void RandomMatrix(double * A, size_t N)
{
    srand(time(NULL));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i * N + j] = rand() / RAND_MAX;
        }
    }
}