#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

void ZeroMatrix(double * A, size_t N) {
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) {
            A[i * N + j] = 0.0;
        }
    }
}

void IdentityMatrix(double * A, size_t N) {
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) {
            if ((i * N + j) % (N + 1) == 0)
                A[i * N + j] = 1.0;
            else
                A[i * N + j] = 0.0;
        }
    }
}

void CalcMatMulTime_kij_opt(double * A, double * B, double * C, size_t N) {
    size_t i, j, k;

    size_t dummy = 0;

    ZeroMatrix(&C[0], N);


    for (k = 0; k < N; k++)
        for(i = 0; i < N; i++) {
            dummy = i * N;
            for(j = 0; j < N; j++)
                C[dummy + j] = C[dummy + j] + A[dummy + k] * B[k * N + j];
        }

}

void print_matrix(double * A, size_t N) {
    printf("matrix:");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if ((i * N + j) % N == 0) {
                printf("\n");
            }
            printf("%f ", A[i * N + j]);
        }
    }
    printf("\n");
}

double *copy_matrix(double *M, size_t N) {
    double *copy_M = (double *) malloc(N * N * sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            copy_M[i * N + j] = M[i * N + j];
        }
    }
    return copy_M;
}

double *bin_pow_matrix(double *A, size_t N, size_t n_pow) {
    double *tmp_res_A_A, *copy_res, *res, *copy_A_1, *copy_A_2;
    copy_A_1 = copy_matrix(&A[0], N);
    res = (double *) malloc(N * N * sizeof(double));
    IdentityMatrix(&res[0], N);
    while (n_pow)
        if (n_pow & 1) {
            copy_res = copy_matrix(&res[0], N);

            CalcMatMulTime_kij_opt(&copy_A_1[0], &copy_res[0], &res[0], N);

            free(copy_res);
            --n_pow;
        }
        else {
            //a *= a;
            tmp_res_A_A = (double *) malloc(N * N * sizeof(double));
            copy_A_2 = copy_matrix(&copy_A_1[0], N);

            CalcMatMulTime_kij_opt(&copy_A_1[0], &copy_A_2[0], &tmp_res_A_A[0], N);

            free(copy_A_2);
            free(copy_A_1);

            copy_A_1 = tmp_res_A_A;
            n_pow >>= 1;
        }

    return res;
}

void star_graph_4_size(double * A, size_t N) {
    assert(N == 4);
    A[0] = 0; A[1] = 1; A[2] = 1; A[3] = 1;
    A[4] = 1; A[5] = 0; A[6] = 0; A[7] = 0;
    A[8] = 1; A[9] = 0; A[10] = 0; A[11] = 0;
    A[12] = 1; A[13] = 0; A[14] = 0; A[15] = 0;
}

void random_graph_4(double * A, size_t N) {
    assert(N == 4);
    A[0] = 1; A[1] = 1; A[2] = 1; A[3] = 0;
    A[4] = 1; A[5] = 2; A[6] = 1; A[7] = 1;
    A[8] = 1; A[9] = 1; A[10] = 2; A[11] = 1;
    A[12] = 0; A[13] = 1; A[14] = 1; A[15] = 3;
}

int get_number_of_paths_of_length_k(double *A, size_t N, int k) {
    int number_of_paths = 0;

    double *pow_A = bin_pow_matrix(&A[0], N, k);
    // sum matrix elements equal number of path in graph of length k
    for (int i = 0; i < N * N; i++) {
        number_of_paths += pow_A[i];
    }
    return number_of_paths;
}

void read_from_file(double *A, size_t N, char fn[]) {

    FILE *fp;
    int tmp;
    if ((fp = fopen(fn, "r")) == NULL)
    {
        printf("Не удалось открыть файл\n");
        assert(NULL);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(fp, "%d", &tmp);
            A[i * N + j] = tmp;
        }
    }
    fclose(fp);
}

double *prepare_matrix_for_compute_page_rank (double *A, size_t N, double d) {
    double *A_norm;
    double *norm_vec = (double *) malloc(N * sizeof(double));
    A_norm = copy_matrix(&A[0], N);

    // init norm vector
    for (int i = 0; i < N; i++) {
        norm_vec[i] = 0;
    }
    // norming matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            norm_vec[j] += A_norm[i * N + j];
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            if (norm_vec[j] != 0)
                A_norm[i * N + j] /= norm_vec[j];
        }
    }
    free(norm_vec);

    // add dumping
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (A_norm[i * N + j] == 0)
                A_norm[i * N + j] = (1 - d) / N;
            else
                A_norm[i * N + j] = ((1 - d) / N) + d * A_norm[i * N + j];


    return A_norm;
}

void mat_vec_mul(double *A, double *x, double *y, size_t N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            y[i] += A[i * N + j] * x[j];
        }
    }

}

double compute_norm(double *x, size_t N) {
    double norm = 0;
    for (int i = 0; i < N; i++) {
        // norm += x[i] * x[i];
        norm += fabs(x[i]);
    }
    return norm;
}

void devide_vec_on_norm(double *x, size_t N, double norm) {
    for (int i = 0; i < N; i++) {
        x[i] /= norm;
    }
}

double compute_criteria(double *current_v, double *previous_v, size_t N) {
    double *diff_v = (double *) malloc(N * sizeof(double));
    double norm;
    for (int i = 0; i < N; i++) {
        diff_v[i] =  current_v[i] - previous_v[i];
    }
    norm = compute_norm(&diff_v[0], N);
    free(diff_v);
    return norm;
}

void print_vec(double *v, size_t N) {
    printf("vector:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", v[i]);
    }
    printf("\n");
}

double *copy_vector(double *v, size_t N) {
    double *copy_v= (double *) malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        copy_v[i] = v[i];
    }
    return copy_v;
}

double *compute_max_eigenvector (double *A, size_t N, double eps) {
    double *eigenvector = (double *) malloc(N * sizeof(double));
    double *previous_eigenvector;
    double norm;
    double criteria = 1;

    // init eigenvector
    for (int i = 0; i < N; i++) {
        eigenvector[i] = 1;
    }
    previous_eigenvector = copy_vector(eigenvector, N);
    int i = 0;
    while (criteria > eps) {
        //for (int i = 0; i < 100; i++) {
        mat_vec_mul(&A[0], &previous_eigenvector[0], &eigenvector[0], N);
        norm = compute_norm(&eigenvector[0], N);

        if (norm != 0)
            devide_vec_on_norm(&eigenvector[0], N, norm);

        criteria = compute_criteria(&eigenvector[0], &previous_eigenvector[0], N);

        previous_eigenvector = copy_vector(&eigenvector[0], N);
    }

    return eigenvector;
}

int comp(const void *a,const void *b) {
    double *x = (double *) a;
    double *y = (double *) b;

    if (*x < *y) return -1;
    else if (*x > *y) return 1; return 0;
}

double *page_rank(double *A, size_t N) {
    double d = 0.85;
    double eps = 0.0001;
    double *A_norm, *eig_vec;

    A_norm = prepare_matrix_for_compute_page_rank (&A[0], N, d);
    eig_vec = compute_max_eigenvector(A_norm, N, eps);
    qsort(&eig_vec[0], N, sizeof(*eig_vec), comp);

    free(A_norm);
    return eig_vec;
}
