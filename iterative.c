/**
description: solve linear equation with
    iterative methods.
-- Jacobi
-- Gauss_Seidel
-- SOR
*/


#include "matrix.h"
#include "iterative.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define DELTA   (1e-6)

/**
description: solve linear equation with
    jacobi iteration. a matrix style method
    is straightforward and easy to understand.
    here we give the raw method.
    The ending condition is |norm(X(k+1), INF) - norm(X(k), INF)| < 1e-6
*/
Matrix* jacobi_iter(Matrix* A, Matrix* B, double delta)
{
    assert(A->m == A->n);
    assert(A->m == B->m && B->n == 1);
    Matrix* X[2];
    int k, i, j, cnt = 0;
    const int max_loop = 100;
    double norm_inf, last_norm_inf = (double)INF;
    X[0] = create_matrix(A->m, 1); //X is initied to zero vector
    X[1] = create_matrix(A->m, 1);
    for(k = 1; k <= max_loop; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X[0]->mem[j][0];
            sum -= A->mem[i][i] * X[0]->mem[i][0];
            X[1]->mem[i][0] = (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X[1]->mem[i][0] - X[0]->mem[i][0]));
        }
        swap_matrix(X[0], X[1]);
        if(print_steps) {
            printf("===step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X[0]);
        }
        if(norm_inf >= last_norm_inf) {
            if(++cnt >= 5)
                goto DIVERGENT;
        } 
        else cnt = 0;
        if(norm_inf < delta)
            break;
        last_norm_inf = norm_inf;
    }
    if(k > max_loop)
        goto DIVERGENT;
    destroy_matrix(X[1]);
    return X[0];
DIVERGENT:
    destroy_matrix(X[1]);
    destroy_matrix(X[0]);
    printf("error: the jacobi method is not convergent!\n");
    return NULL;
}

/**
description: solve linear equation with
     gauss seidel iteration. as this method's
     definition, it is easier to implement 
     thought it's definition is more complex 
     than jacobi iteration.
*/
Matrix* gauss_seidel_iter(Matrix* A, Matrix* B, double delta)
{
    assert(A->m == A->n);
    assert(A->m == B->m && B->n == 1);
    int k, i, j, cnt = 0;
    const int max_loop = 100;
    double norm_inf, last_norm_inf = (double)INF;
    Matrix* X = create_matrix(A->m, 1); //X is initied to zero vector
    for(k = 1; k <= max_loop; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double xi, sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X->mem[j][0];
            sum -= A->mem[i][i] * X->mem[i][0];
            xi = X->mem[i][0];
            X->mem[i][0] = (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X->mem[i][0] - xi));
        }
        printf("nrom: %lf\n", norm_inf);
        if(print_steps) {
            printf("===step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X);
        }
        if(norm_inf >= last_norm_inf) {
            if(++cnt >= 5)
                goto DIVERGENT;
        } 
        else cnt = 0;
        if(norm_inf < delta)
            break;
        last_norm_inf = norm_inf;
    }
    if(k > max_loop)
        goto DIVERGENT;
    return X; 
DIVERGENT:    
    destroy_matrix(X);
    printf("error: the Gauss Seidel method is not convergent!\n");
    return NULL; 
}

/**
description: solve linear equation with
    successive over relaxation methond.
*/
Matrix* sor(Matrix* A, Matrix* B, double w, double delta)
{
    assert(A->m == A->n);
    assert(A->m == B->m && B->n == 1);
    int k, i, j, cnt = 0;
    const int max_loop = 100;
    double norm_inf, last_norm_inf = (double)INF;
    Matrix* X = create_matrix(A->m, 1); //X is initied to zero vector
    for(k = 1; k <= max_loop; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double xi, sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X->mem[j][0];
            //sum -= A->mem[i][i] * X->mem[i][0];
            xi = X->mem[i][0];
            X->mem[i][0] += w * (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X->mem[i][0] - xi));
        }
        if(print_steps) {
            printf("===step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X);
        }
        if(norm_inf >= last_norm_inf) {
            if(++cnt >= 5)
                goto DIVERGENT;
        } 
        else cnt = 0;
        if(norm_inf < delta)
            break;
        last_norm_inf = norm_inf;
    }
    if(k > max_loop)
        goto DIVERGENT;
    return X;
DIVERGENT:
    destroy_matrix(X);
    printf("error: the SOR method is not convergent!\n");
    return NULL;
}