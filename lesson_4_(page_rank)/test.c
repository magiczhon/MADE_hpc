#include <assert.h>
#include "helper.c"


void test_bin_pow() {
    int N = 3;
    size_t pow_n = 7;
    double *A, *C;
    A = (double *) malloc(N * N * sizeof(double));
    A[0] = 1.5;
    A[1] = 2;
    A[2] = 3;
    A[3] = -1;
    A[4] = 0;
    A[5] = 1;
    A[6] = 3;
    A[7] = 1;
    A[8] = 2;

    C = bin_pow_matrix(&A[0], N, pow_n);
//    print_matrix(&C[0], N);
//    printf("c[0] = %f\n", C[0]);
    assert(C[0] == 23449.5859375);
//    printf("c[1] = %f\n", C[1]);
    assert(C[1] == 16304.3125);
//    printf("c[2] = %f\n", C[2]);
    assert(C[2] == 31114.671875);
//    printf("c[3] = %f\n", C[3]);
    assert(C[3] == 408.390625);
//    printf("c[4] = %f\n", C[4]);
    assert(C[4] == 282.875);
//    printf("c[5] = %f\n", C[5]);
    assert(C[5] == 542.78125);
//    printf("c[6] = %f\n", C[6]);
    assert(C[6] == 25407.640625);
//    printf("c[7] = %f\n", C[7]);
    assert(C[7] == 17663.875);
//    printf("c[8] = %f\n", C[8]);
    assert(C[8] == 33708.28125);
    printf("test bin pow matrix complete\n");
}

// task1
void test_path_count_in_graph() {
    size_t N = 4;
    double *A;
    A = (double *) malloc(N * N * sizeof(double));
    // test star graph
    star_graph_4_size(&A[0], N);
    // len 2
    assert(12 == get_number_of_paths_of_length_k(&A[0], N, 2));
    // len 3
    assert(18 == get_number_of_paths_of_length_k(&A[0], N, 3));


    // test random graph
    random_graph_4(A, N);
    // len 4
    assert(1852 == get_number_of_paths_of_length_k(&A[0], N, 4));
    // len 5
    assert(8714 == get_number_of_paths_of_length_k(&A[0], N, 5));
    free(A);
    printf("test path count in graph complete\n");
}

void test_page_rank() {
    double *A, *eig_vec;
    const size_t N = 11;
    A = (double *) malloc(N * N * sizeof(double));

    char fn[] = "matrix_adj.txt";
    read_from_file(&A[0], N, fn);

    eig_vec = page_rank(&A[0], N);

    double s = 0;
    for (int i = 0; i < N; i++)
        s+=eig_vec[i];
    assert(1.0 == s);
    double test_values_eigvec[N] = {0.021991, 0.039404, 0.069534, 0.076080, 0.078978,
                                    0.078978, 0.096391, 0.096391,0.096391, 0.096391, 0.249470};
    for(int i = 0; i < N; i++) {
        // округляем до 2-го знака
        assert(round(test_values_eigvec[i]*100)/100 == round(eig_vec[i]*100)/100);
    }
    free(A);
    free(eig_vec);
    printf("test page rank complete");
}