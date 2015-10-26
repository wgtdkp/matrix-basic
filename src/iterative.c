/**
description: solve linear equation with
    iterative methods.
-- Jacobi
-- Gauss_Seidel
-- SOR
*/


#include "matrix.h"
#include "iterative.h"
#include "eigenvalue.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define DELTA   (1e-6)

/**
description: straightforward method of checking
    if jacobi iteration is convergent
*/
bool is_jacobi_convergent(Matrix* A, bool step)
{
    assert(IS_SQUARE(A));
    int i, j;
    Matrix *D_inv, *L_U, *J;
    double rho;
    D_inv = create_matrix_n(A->m);
    L_U = create_matrix_n(A->m);
    for(i = 0; i < A->m; i++)
        D_inv->mem[i][i] = 1.0 / A->mem[i][i];
    for(i = 0; i < A->m; i++)
        for(j = 0; j < A->n; j++) {
            if(i !=  j)
                L_U->mem[i][j] = -A->mem[i][j];
        }
    J = mul(D_inv, L_U);
    rho = fabs(pow_method(J, false));
    if(step) {
        printf("J : \n");
        print_matrix(J);
        printf("R(J) : %.8lf\n", rho);
    }
    destroy_matrix(&D_inv);
    destroy_matrix(&L_U);
    destroy_matrix(&J);
    return rho < 1.0;
}

/**
description: solve linear equation with
    jacobi iteration. a matrix style method
    is straightforward and easy to understand.
    here we give the raw method.
    The ending condition is |norm(X(k+1), INF) - norm(X(k), INF)| < 1e-6
*/
Matrix* jacobi_iter(Matrix* A, Matrix* B, double delta, bool step)
{
    assert(IS_EQUATION(A, B));

    Matrix* X[2];
    int k, i, j;
    double norm_inf;
    
    if(!is_jacobi_convergent(A, false)) {
        printf("error: the Jacobi iteration is not convergent!\n");
        return NULL;            
    }
    X[0] = create_matrix(A->m, 1); //X is initied to zero vector
    X[1] = create_matrix(A->m, 1);
    for(k = 1; ; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X[0]->mem[j][0];
            sum -= A->mem[i][i] * X[0]->mem[i][0];
            X[1]->mem[i][0] = (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X[1]->mem[i][0] - X[0]->mem[i][0]));
        }
        swap_matrix(X[0], X[1]);
        if(step) {
            printf("=== step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X[0]);
        }
        if(norm_inf < delta)
            break;
    }
    destroy_matrix(&X[1]);
    return X[0];
}


//it equals to is_sor_convergent() when param 'w' == 1
bool is_gauss_seidel_convergent(Matrix* A, bool step)
{
    assert(IS_SQUARE(A));
    int i, j;
    Matrix *D_L, *D_L_inv, *U, *G;
    double rho;
    D_L = create_matrix_n(A->m);
    U = create_matrix_n(A->m);
    for(i = 0; i < A->m; i++)
        for(j = 0; j < A->n; j++) {
            if(i >= j)
                D_L->mem[i][j] = A->mem[i][j];
            else
                U->mem[i][j] = -A->mem[i][j];
        }

    D_L_inv = inverse(D_L);
    G = mul(D_L_inv, U);
    rho = fabs(pow_method(G, false));
    if(step) {
        printf("inverse of D - L : \n");
        print_matrix(D_L_inv);
        printf("G : \n");
        print_matrix(G);
        printf("R(G) : %.8lf\n", rho);
    }
    destroy_matrix(&U);
    destroy_matrix(&D_L);
    destroy_matrix(&D_L_inv);
    destroy_matrix(&G);
    return rho < 1.0;
}

/**
description: solve linear equation with
     gauss seidel iteration. as this method's
     definition, it is easier to implement 
     thought it's definition is more complex 
     than jacobi iteration.
*/
Matrix* gauss_seidel_iter(Matrix* A, Matrix* B, double delta, bool step)
{
    assert(IS_EQUATION(A, B));

    int k, i, j;
    double norm_inf;
    
    if(!is_gauss_seidel_convergent(A, false)) {
        printf("error: the Gauss-Seidel iteration is not convergent!\n");
        return NULL;            
    }
    Matrix* X = create_matrix(A->m, 1); //X is initied to zero vector
    for(k = 1; ; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double xi, sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X->mem[j][0];
            sum -= A->mem[i][i] * X->mem[i][0];
            xi = X->mem[i][0];
            X->mem[i][0] = (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X->mem[i][0] - xi));
        }
        //printf("nrom: %lf\n", norm_inf);
        if(step) {
            printf("=== step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X);
        }
        if(norm_inf < delta)
            break;
    }
    return X; 
}

bool is_sor_convergent(Matrix* A, double w, bool step)
{
    assert(IS_SQUARE(A));
    int i, j;
    Matrix *D, *L, *U, *wL, *wU;
    Matrix *D_wL, *D_wL_inv, *D_wD, *D_wD_wU;
    Matrix* Lw;
    double rho;
    D = create_matrix_n(A->m);
    L = create_matrix_n(A->m);
    U = create_matrix_n(A->m);
    for(i = 0; i < D->m; i++)
        D->mem[i][i] = A->mem[i][i];
    for(i = 0; i < A->m; i++)
        for(j = 0; j < A->n; j++) {
            if(i > j)
                L->mem[i][j] = -A->mem[i][j];
            else if(i < j)
                U->mem[i][j] = -A->mem[i][j];
        }
    wL = mul_cons(L, w);
    wU = mul_cons(U, w);
    D_wL = sub(D, wL);
    D_wL_inv = inverse(D_wL);
    D_wD = mul_cons(D, 1 - w);
    D_wD_wU = add(D_wD, wU);
    Lw = mul(D_wL_inv, D_wD_wU);
    rho = fabs(pow_method(Lw, false));
    if(step) {
        printf("inverse of D - wL : \n");
        print_matrix(D_wL_inv);
        printf("Lw : \n");
        print_matrix(Lw);
        printf("R(Lw) : %.8lf\n", rho);
    }
    destroy_matrix(&D);
    destroy_matrix(&L);
    destroy_matrix(&U);
    destroy_matrix(&wL);
    destroy_matrix(&wU);
    destroy_matrix(&D_wL);
    destroy_matrix(&D_wL_inv);
    destroy_matrix(&D_wD);
    destroy_matrix(&D_wD_wU);
    destroy_matrix(&Lw);
    return rho < 1.0;
}

/**
description: solve linear equation with
    successive over relaxation methond.
*/
Matrix* sor(Matrix* A, Matrix* B, double w, double delta, bool step)
{
    assert(IS_EQUATION(A, B));

    int k, i, j;
    double norm_inf;
    if(!is_sor_convergent(A, w, false)) {
        printf("error: the SOR iteration is not convergent!\n");
        return NULL;
    }
    Matrix* X = create_matrix(A->m, 1); //X is initied to zero vector
    for(k = 1; ; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double xi, sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X->mem[j][0];
            //sum -= A->mem[i][i] * X->mem[i][0];
            xi = X->mem[i][0];
            X->mem[i][0] += w * (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X->mem[i][0] - xi));
        }
        if(step) {
            printf("=== step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X);
        }
        if(norm_inf < delta)
            break;
    }
    return X;
}

Matrix* sor_homework(Matrix* A, Matrix* B, double w, double delta, bool step)
{
    assert(IS_EQUATION(A, B));

    int k, i, j;
    double norm_inf;
    double xi[] = {0.5, 1.0, -0.5};
    if(!is_sor_convergent(A, w, false)) {
        printf("error: the SOR iteration is not convergent!\n");
        return NULL;
    }
    Matrix* X = create_matrix(A->m, 1); //X is initied to zero vector
    for(k = 1; ; k++) {
        for(i = 0, norm_inf = .0; i < A->m; i++) {
            double sum;//xi, sum;
            for(j = 0, sum = 0; j < A->n; j++)
                sum += A->mem[i][j] * X->mem[j][0];
            //sum -= A->mem[i][i] * X->mem[i][0];
            //xi = X->mem[i][0];

            X->mem[i][0] += w * (B->mem[i][0] - sum) / A->mem[i][i];
            norm_inf = MAX(norm_inf, fabs(X->mem[i][0] - xi[i]));
        }
        if(step) {
            printf("=== step %d ===\n", k);
            printf("delta: %lf\n", norm_inf);
            printf("X: \n");
            print_matrix(X);
        }
        if(norm_inf < delta)
            break;
    }
    return X;
}