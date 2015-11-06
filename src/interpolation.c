#include <stdio.h>
#include "interpolation.h"
#include "direct.h"

Poly* vandermonde(const double* xarr, const double* yarr, int n)
{
    int i, j;
    Poly* res = NULL;   //important
    Matrix* A = create_matrix_n(n); //vandermonde
    Matrix* B = create_matrix(n, 1);
    Matrix* X;
    
    for (i = 0; i < n; i++) {
        A->mem[i][0] = 1;
        for (j = 1; j < n; j++)
            A->mem[i][j] = A->mem[i][j-1] * xarr[i];
        B->mem[i][0] = yarr[i];
    }
    //print_matrix(A);
    //print_matrix(B);
    X = gauss_elim(A, B, false);
    for (i = 0; i < n; i++) {
        poly_add_inp(&res, create_poly(X->mem[i][0], i));
    }
    destroy_matrix(&A);
    destroy_matrix(&B);
    destroy_matrix(&X);
    return res; //we didn't create res, we appended to it
}

Poly* lagrange(const double* xarr, const double* yarr, int n)
{
    int i;
    Poly *omega = create_poly(1, 0);
    Poly *res = NULL, *item = NULL;
    Poly* omega_deriv;
    Poly* tmp = create_poly(1, 1);
    poly_add_inp(&tmp, create_poly(1, 0));
    for (i = 0; i < n; i++) {
        //printf("i = %d\n", i);
        tmp->next->coeff = -xarr[i];
        poly_mul_inp(&omega, tmp);
    }
    //print_poly(omega);
    omega_deriv = poly_copy(omega);
    poly_deriv_inp(&omega_deriv);
    for (i = 0; i < n; i++) {
        //printf("i = %d\n", i);
        tmp->next->coeff = -xarr[i];
        poly_div_inp(&omega, tmp);
        //print_poly(omega);
        item = create_poly(yarr[i] / poly_value(omega_deriv, xarr[i]), 0);
        poly_mul_inp(&item, omega);    
        poly_add_inp(&res, item);    //we don't need to destroy item
        poly_mul_inp(&omega, tmp);   //reset omega
    }
    destroy_poly(&omega);
    destroy_poly(&omega_deriv);
    destroy_poly(&tmp);
    return res;
}

Poly* newton(const double* xarr, const double* yarr, int n)
{
    int i, j;
    Poly *res = NULL, *item = NULL;
    Poly* omega = create_poly(1, 0);
    Poly* omega_deriv = NULL;
    Poly* tmp = create_poly(1, 1);
    poly_add_inp(&tmp, create_poly(1, 0));
    
    for (i = 0; i < n; i++) {
        double coeff;
        tmp->next->coeff = -xarr[i];
        item = poly_copy(omega);
        poly_mul_inp(&omega, tmp);
        omega_deriv = poly_copy(omega);
        poly_deriv_inp(&omega_deriv);
        for (j = 0, coeff = 0; j <= i; j++) {
            coeff += yarr[j] / poly_value(omega_deriv, xarr[j]);
        }
        poly_mul_cons_inp(&item, coeff);
        poly_add_inp(&res, item);
    }
    destroy_poly(&omega);
    destroy_poly(&omega_deriv);
    destroy_poly(&tmp);
    return res;
}

/**
 *difference quotinet
 */
double static inline diff_quot(const double* xarr, const double* yarr, int i, int j)
{
    return (yarr[j] - yarr[i]) / (xarr[j] - xarr[i]);
}

static Poly* construct_spline(const double* xarr, const double* yarr, int i, const Matrix* M)
{
    Poly* res = NULL;
    double h = xarr[i+1] - xarr[i];
    Poly* tmp = create_poly(-1, 1);
    poly_add_inp(&tmp, create_poly(xarr[i+1], 0));
    poly_pow_inp(&tmp, 3);
    poly_mul_cons_inp(&tmp, M->mem[i][0] / (6 * h));
    poly_add_inp(&res, tmp);    // actually, tmp is merged to res

    tmp = create_poly(1, 1);
    poly_add_inp(&tmp, create_poly(-xarr[i], 0));
    poly_pow_inp(&tmp, 3);
    poly_mul_cons_inp(&tmp, M->mem[i+1][0] / (6 * h));
    poly_add_inp(&res, tmp);

    tmp = create_poly(-1, 1);
    poly_add_inp(&tmp, create_poly(xarr[i+1], 0));
    poly_mul_cons_inp(&tmp, (yarr[i] - M->mem[i][0] * h * h / 6) / h);
    poly_add_inp(&res, tmp);

    tmp = create_poly(1, 1);
    poly_add_inp(&tmp, create_poly(-xarr[i], 0));
    poly_mul_cons_inp(&tmp, (yarr[i+1] - M->mem[i+1][0] * h * h / 6) / h);
    poly_add_inp(&res, tmp);
    return res;
}   

/**
 *@type: 0, specify derivative at section boundary;
 *       1, specify second derivative at section boundary
 *       2, periodic function
 */
Poly** spline(const double* xarr, const double* yarr, int n, double a, double b, int type)
{
    int i;
    double lambda, u;
    Poly** res;
    Matrix* M;
    Matrix* d = create_vector(n);
    Matrix* A = create_matrix_n(n);
    
    //init A
    for (i = 0; i < n; i++)
        A->mem[i][i] = 2;
    
    for (i = 1; i < n - 1; i++) {
        lambda = (xarr[i+1] - xarr[i]) / (xarr[i+1] - xarr[i-1]);
        A->mem[i][i+1] = lambda;
        u = 1 - lambda;
        A->mem[i][i-1] = u;
    }

    //init d
    for (i = 1; i < n - 1; i++) {
        double d_tmp = diff_quot(xarr, yarr, i, i+1) - diff_quot(xarr, yarr, i-1, i);
        d_tmp *= 6 / (xarr[i+1] - xarr[i-1]);
        d->mem[i][0]  = d_tmp;
    }

    //init bundary conditions
    if (0 == type) {
        lambda = 1;
        A->mem[0][1] = lambda;
        u = 1;
        A->mem[n-1][n-2] = u;
        d->mem[0][0] = 6 * (diff_quot(xarr, yarr, 0, 1) - a) / (xarr[1] - xarr[0]);
        d->mem[n-1][0] = 6 * (b - diff_quot(xarr, yarr, n-2, n-1)) / (xarr[n-1] - xarr[n-2]);
    } else if (1 == type) {
        lambda = 0, u = 0;
        A->mem[0][1] = lambda;
        A->mem[n-1][n-2] = u;
        d->mem[0][0] = 2 * a;
        d->mem[n-1][0] = 2 * b;
    } else if (2 == type) {

    } else {
        printf("error: bad type = %d!\n", type);
        res = NULL;
        goto RETURN;
    }

    res = (Poly**)malloc(sizeof(Poly*) * (n-1));
    memset((void*)res, 0, sizeof(Poly*) * (n-1));
    //print_matrix(A);
    //print_matrix(d);
    M = chasing_method(A, d, false);
    //print_matrix(M);
    for (i = 0; i < n-1; i++) {
        res[i] = construct_spline(xarr, yarr, i, M);
    }
    destroy_matrix(&M);
RETURN:
    destroy_matrix(&A);
    destroy_matrix(&d);
    return res;
}
