#include "matrix.h"
#include "eigenvalue.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/**
description : calculate the biggest eigenvalue
    with the pow method.
*/
double pow_method(Matrix* A, bool step)
{
    assert(A->m == A->n);
    Matrix* V = create_matrix(A->m, 1);
    Matrix *U, *last_U = V;
    double lambda1;
    int k;
    const int max_loop = 300;
    fill_matrix(V, 1);
    for(k = 1; k <= max_loop; k++) {
        U = mul(A, last_U);
        if(is_const_similar(U, last_U, &lambda1)) {
            destroy_matrix(&U);
            destroy_matrix(&last_U);
            return lambda1;
        }
        if(step) {
            printf("=== step %d ===\n", k);
            printf("U:\n");
            print_matrix(U);
        }
        destroy_matrix(&last_U);
        last_U = U;
    }
    return 0;
}

