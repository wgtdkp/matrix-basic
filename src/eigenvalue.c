#include "matrix.h"
#include "eigenvalue.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/**
description : calculate the biggest eigenvalue
    with the pow method.
*/
double pow_method(Matrix* A, double delta, bool step)
{
    assert(A->m == A->n);
    Matrix* V = create_matrix(A->m, 1);
    Matrix *U, *last_U = V;
    double lambda1;
    int k, i;
    const int max_loop = 300;
    fill_matrix(V, 1);
    for(k = 1; k <= max_loop; k++) {
        double max = 0;
        U = mul(A, last_U);

        for (i = 0; i < U->m; i++)
            max = MAX(max, fabs(U->mem[i][0]));
        for (i = 0; i < U->m; i++)
            U->mem[i][0] /= max;
        
        if(is_const_similar(U, last_U, delta, &lambda1)) {
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

/**
create a Householder reflect matrix with vector specified.
@n: the size of return matrix, it is bigger than the vector
size, and then the left right sub area of return matrix is an eye.
*/
Matrix* householder(double coeff, int n, Matrix* w)
{
    assert(IS_VECTOR(w));
    int i, j;
    Matrix* ret = create_eye(n);
    for (i = ret->m - 1; i >= 0; i--) {
        for (j = ret->n - 1; j >= 0; j--) {
            ret->mem[n - ret->m + i][n ->ret->n + j] -=
                coeff * w->mem[i][0] * w->mem[j][0];
        }
    }

    return ret;
}

/**
@method: 'H' stands for Householder trans, 'G' stands for Givens trans
*/
void qr_decomp(const Matrix* A, const char method, Matrix* Q, Matrix* R, bool step)
{
    int i, j;
    Matrix *q = NULL, *q_tmp, *r = NULL, *r_tmp;
    Matrix* R_tmp;
    if (NULL == Q || NULL == R)
        return;
    if ('H' == method) {
        double sigma, beta;
        Matrix* x = create_matrix(A->m, 1);
        Matrix* u = x;
        Matrix* h = NULL;
        make_eye(Q);
        copy_inp(R, A);
        for (i = 0; i < R->m - 1; i++) {
            for (j = i; j < R->m; j++)
                x->mem[j-i][0] = R->mem[j][i];
            x->m = R->m - i;

            sigma = SIGN(x[0][0]) * norm(x, 2, false);
            beta = sigma * (sigma + x->mem[0][0]);
            u->mem[0][0] = x->mem[0][0] + sigma;    //it modifies x
            h = householder(1.0/beta, R->m, u);

            mul_inp(Q, h);
            R_tmp = mul(h, R);
            copy_inp(R, R_tmp);
            if (step) {
                printf("H%d:\n", i + 1);
                print_matrix(h);
                printf("Q%d:\n", i + 1);
                print_matrix(Q);
                printf("R%d:\n", i + 1);
                print_matrix(R);
            }
            destroy_matrix(h);
            destroy_matrix(&R_tmp);
        }
        x->m = A->m;    //reset the size of x, otherwise, will cause memory leak
        destroy_matrix(&x);
    } else if('G' == method) {

    }
}