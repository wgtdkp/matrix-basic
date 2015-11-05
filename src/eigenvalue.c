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
static Matrix* householder(double coeff, int n, Matrix* w)
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
 *QR decomposition with Householder transformation
 *@return: return the Q matrix 
 *    and R matrix is stored in A
 *note: this function modifies A
 */
Matrix* qr_decomp_H(Matrix* A, bool step)
{
    int i, j;
    double sigma, beta;
    Matrix* X = create_matrix(A->m, 1);
    Matrix* U = X;
    Matrix* H = NULL;
    Matrix* Q = create_eye(A->m);
    for (i = 0; i < A->m - 1; i++) {
        for (j = i; j < A->m; j++)
            X->mem[j-i][0] = A->mem[j][i];
        X->m = A->m - i;

        sigma = SIGN(X[0][0]) * norm(X, 2, false);
        beta = sigma * (sigma + X->mem[0][0]);
        U->mem[0][0] = X->mem[0][0] + sigma;    //it modifies X
        H = householder(1.0/beta, A->m, U);

        mul_inp_L(Q, H);    //store result in Q
        mul_inp_R(H, A);    //store result in R(A)

        if (step) {
            printf("H%d:\n", i + 1);
            print_matrix(H);
            printf("Q%d:\n", i + 1);
            print_matrix(Q);
            printf("A(R)%d:\n", i + 1);
            print_matrix(A);
        }
        destroy_matrix(&H);
    }
    X->m = A->m;    //reset the size of X, otherwise, will cause memory leak
    destroy_matrix(&X);
    return Q;
}

/**
 *QR decomposition with Givens transformation
 *@return: return the Q matrix
 *    and R matrix is stored in A
 *note: this function modifies A
 */
Matrix* qr_decomp_G(Matrix* A, bool step)
{
    int i, j, k;
    Matrix* P = create_eye(A->m);
    Matrix* Q = create_eye(A->m);

    for (i = 0; i < A->m; i++) {
        for (k = i + 1; k < A->m; k++) {
            double r = A->mem[i][i] * A->mem[i][i] + A->mem[k][i] * A->mem[k][i];
            r = pow(r, 0.5);
            P->mem[i][i] = A->mem[i][i] / r;
            P->mem[i][k] = P->mem[k][i] / r;
            P->mem[k][i] = -P->mem[i][k];
            P->mem[k][k] = P->mem[i][i];

            A->mem[i][i] = r;
            A->mem[k][i] = 0;
            mul_inp_L(Q, P);    //store result in Q
            mul_inp_R(P, A);    //store result in A
            make_eye(P);
            if (step) {
                printf("P%d:\n", i + 1);
                print_matrix(P);
                printf("Q%d:\n", i + 1);
                print_matrix(Q);
                printf("A(R)%d:\n", i + 1);
                print_matrix(A);
            }
        }
    }
    destroy_matrix(&P);
    return Q;
}

void to_hessenberg(Matrix* A, bool step)
{
    int i, j;
    double sigma, beta;
    Matrix* X = create_matrix(A->m - 1, 1);
    Matrix* U = X;
    Matrix* H = NULL;
    Matrix* Q = create_eye(A->m);
    for (i = 0; i < A->m - 1; i++) {
        for (j = i; j < A->m; j++)
            X->mem[j-i][0] = A->mem[j][i];
        X->m = A->m - 1 - i;

        sigma = SIGN(X[0][0]) * norm(X, 2, false);
        beta = sigma * (sigma + X->mem[0][0]);
        U->mem[0][0] = X->mem[0][0] + sigma;    //it modifies X
        H = householder(1.0/beta, A->m, U);

        //Q = QH
        mul_inp_L(Q, H);    //store result in Q

        //A = HAH
        mul_inp_R(H, A);    //store result in R(A)
        mul_inp_L(A, H);

        if (step) {
            printf("H%d:\n", i + 1);
            print_matrix(H);
            printf("Q%d:\n", i + 1);
            print_matrix(Q);
            printf("A(R)%d:\n", i + 1);
            print_matrix(A);
        }
        destroy_matrix(&H);
    }
    X->m = A->m - 1;    //reset the size of X, otherwise, will cause memory leak
    destroy_matrix(&X);
    destroy_matrix(&Q);
    //return Q;    
}

/**
 *QR method to solve the eigenvalue problem
 *@return: the transformed matrix is stored in A
 */
void qr_method(Matrix* A, double delta, bool step)
{
    int i;
    Matirx *Q, *A_last;
    to_hessenberg(A, false);
    if (step) {
        printf("convert A matrix to Hessenberg: \n");
        print_matrix(A);
    }
    for (i = 0; ; i++) {
        A_last = copy(A);
        Q = qr_decomp_H(A, false);
        mul_inp_L(A, Q);
        destroy_matrix(&Q);
        if (is_norm_similar(A_last, A, 2, delta)) {
            destroy_matrix(A_last);
            break;;
        }
        destroy_matrix(&A_last);
        if (step) {
            printf("step %d: \n", i + 1);
            printf("Q matrix: \n");
            print_matrix(Q);
            printf("A matrix: \n");
            print_matrix(A);
        }
    }
}