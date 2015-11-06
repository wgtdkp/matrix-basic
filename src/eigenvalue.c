#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "eigenvalue.h"
#include "norm.h"
#include "nonlinear.h"
#include "poly.h"
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
    for (i = 0; i < w->m; i++) {
        for (j = 0; j < w->m; j++) {
            ret->mem[n - w->m + i][n - w->m + j] -=
                coeff * w->mem[i][0] * w->mem[j][0];
        }
    }

    return ret;
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
    Matrix *R = A, *Q_tmp;
    Matrix* P = create_eye(A->m);
    Matrix* Q = create_eye(A->m);
    for (i = 0; i < A->m; i++) {
        for (k = i + 1; k < A->m; k++) {
            double r = A->mem[i][i] * A->mem[i][i] + A->mem[k][i] * A->mem[k][i];
            r = pow(r, 0.5);
            P->mem[i][i] = A->mem[i][i] / r;
            P->mem[i][k] = A->mem[k][i] / r;
            P->mem[k][i] = -P->mem[i][k];
            P->mem[k][k] = P->mem[i][i];

            mul_inp_R(P, Q);    //store result in Q
            //A->mem[i][i] = r;
            //A->mem[k][i] = 0;
            mul_inp_R(P, A);    //store result in A
            if (step) {
                printf("P%d:\n", i + 1);
                print_matrix(P);
                printf("Q%d:\n", i + 1);
                print_matrix(Q);
                printf("A(R)%d:\n", i + 1);
                print_matrix(A);
                printf("\n");
            }
            make_eye(P);
        }
    }
    Q_tmp = Q;
    Q = inverse(Q);
    destroy_matrix(&Q_tmp);
    for (i = 0; i < Q->m; i++) {
        if (SIGN(R->mem[i][i]) < 0) {
            for (j = 0; j < R->n; j++)
                R->mem[i][j] = 0 - R->mem[i][j];
            for (j = 0; j < Q->m; j++)
                Q->mem[j][i] = 0 - Q->mem[j][i];
        }
    }
    destroy_matrix(&P);
    return Q;
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
    Matrix *U = X, *R = A;
    Matrix* H = NULL;
    Matrix* Q = create_eye(A->m);
    for (i = 0; i < A->m - 1; i++) {
        for (j = i; j < A->m; j++)
            X->mem[j-i][0] = A->mem[j][i];
        X->m = A->m - i;
        sigma = SIGN(X->mem[0][0]) * norm(X, 2, false);
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
            printf("\n");
        }
        destroy_matrix(&H);
    }
    for (i = 0; i < Q->m; i++) {
        if (SIGN(R->mem[i][i]) < 0) {
            for (j = 0; j < R->n; j++)
                R->mem[i][j] = 0 - R->mem[i][j];
            for (j = 0; j < Q->m; j++)
                Q->mem[j][i] = 0 - Q->mem[j][i];
        }
    }

    X->m = A->m;    //reset the size of X, otherwise, will cause memory leak
    destroy_matrix(&X);
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
    for (i = 0; i < A->m - 2; i++) {
        for (j = i; j + 1 < A->m; j++)
            X->mem[j-i][0] = A->mem[j+1][i];
        X->m = A->m - 1 - i;

        sigma = SIGN(X->mem[0][0]) * norm(X, 2, false);
        beta = sigma * (sigma + X->mem[0][0]);
        U->mem[0][0] = X->mem[0][0] + sigma;    //it modifies X
        H = householder(1.0/beta, A->m, U);

        //Q = QH
        mul_inp_R(H, Q);    //store result in Q

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

int eigenvalues(Complex* eigenvals, Matrix* A)
{
    int i, k;
    Poly* equaltion = NULL;
    if (NULL == eigenvals)
        return 0;
    poly_add_inp(&equaltion, create_poly(1, 2));
    poly_add_inp(&equaltion, create_poly(1, 1));
    poly_add_inp(&equaltion, create_poly(1, 0));
    for (i = 0, k = 0; i < A->m; i++, k++) {
        //printf("i = %d\n", i);
        if (i == A->m - 1) {
            eigenvals[k].real = A->mem[i][i];
            eigenvals[k].imag = 0;
        } else {
            equaltion->next->coeff = -(A->mem[i][i] + A->mem[i+1][i+1]);
            equaltion->next->next->coeff = 
                A->mem[i][i] * A->mem[i+1][i+1] - A->mem[i][i+1] * A->mem[i+1][i];
            quadratic(eigenvals + k, equaltion);
            if (DOUBLE_EQUAL_DELTA(eigenvals[k].imag, 0, 1e-4)) { // signle eigenvalue
                eigenvals[k].real = A->mem[i][i];
            } else {
                ++i, ++k;
            }
        }
    }
    destroy_poly(&equaltion);
    return k;
}

/**
 *QR method to solve the eigenvalue problem
 *@return: the transformed matrix is stored in A
 */
void qr_method(Matrix* A, double delta, bool step)
{
    int i, k;
    Matrix* Q;
    int eigen_n, last_eigen_n;
    Complex *eigenvals, *last_eigenvals;
    int size = sizeof(Complex) * 2 * (A->m - 1);
    eigenvals = (Complex*)malloc(size);
    last_eigenvals = (Complex*)malloc(size);
    
    to_hessenberg(A, false);
    if (step) {
        printf("convert A matrix to Hessenberg: \n");
        print_matrix(A);
    }
    for (k = 0; ; k++) {
        Q = qr_decomp_H(A, false);
        mul_inp_L(A, Q);    //A is R now, B = RQ
        
        if (step) {
            printf("step %d: \n", k + 1);
            printf("Q matrix: \n");
            print_matrix(Q);
            printf("A matrix: \n");
            print_matrix(A);
            printf("\n");
        }
        eigen_n = eigenvalues(eigenvals, A);

        //TODO: diff between current eigenvalues and last eigenvalues
        if (eigen_n == last_eigen_n) {
            for (i = 0; i < eigen_n; i++)
                if (!COMPLEX_EQUAL_DELTA(last_eigenvals[i], eigenvals[i], delta))
                    break;
            if (i == eigen_n) break;
        }

        memcpy((void*)last_eigenvals, (void*)eigenvals, sizeof(Complex) * eigen_n);
        last_eigen_n = eigen_n;
        destroy_matrix(&Q);
    }
    destroy_matrix(&Q);
    free(eigenvals);
    free(last_eigenvals);
}

