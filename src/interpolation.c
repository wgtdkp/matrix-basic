#include "interpolation.h"
#include "direct.h"


Poly* vandermonde(const double* xarr, const double* yarr, int n)
{
    int i, j;
    Poly* res;
    Matrix* A = create_matrix_n(n + 1); //vandermonde
    Matrix* B = create_matrix(n + 1, 1);
    Matrix* X;
    
    for (i = 0; i < n + 1; i++) {
        A->mem[i][0] = 1;
        for (j = 1; j <= n; j++)
            A->mem[i][j] = A->mem[i][j-1] * xarr[i];
        B->mem[i][0] = yarr[i];
    }

    X = gauss_elim(A, B, false);
    for (i = 0; i < n + 1; i++) {
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
    tmp->next = create_poly(0, 0);
    
    for (i = 0; i < n; i++) {
        tmp->next->coeff = -xarr[i];
        poly_mul_inp(&omega, tmp);
    }
    omega_deriv = copy(omega);
    poly_deriv(&omega_deriv);
    for (i = 0; i < n; i++) {
        tmp->next->coeff = -xarr[i];
        poly_div(&omega, tmp);
        item = create_poly(yarr[i] / poly_value(omega_deriv, xarr[i]), 0);
        poly_mul_inp(&item, omega);    
        poly_add_inp(&res, item);    //we don't need to destroy item
        poly_mul_inp(&omega, tmp);   //reset omega
    }
    destroy_poly(&omega);
    destroy_poly(&omega_deriv);
    destroy_tmp(&tmp);
    return res;
}

Poly* newton(const double* xarr, const double* yarr, int n)
{
    int i, j;
    Poly *res = NULL, *item = NULL;
    Poly *omega = create_poly(1, 0)
    Poly *omega_deriv = NULL;
    Poly *tmp = create_poly(1, 1);
    tmp->next = create_poly(0, 0);
    
    for (i = 0; i < n + 1; i++) {
        double coeff;
        tmp->next->coeff = -xarr[i];
        poly_mul_inp(&omega, tmp);
        omega_deriv = copy(omega);
        poly_deriv(&omega_deriv);
        for (j = 0, coeff = 0; j < i; j++) {
            coeff += yarr[j] / poly_value(omega, xarr[j]);
        }
        item = create_poly(coeff, 0);
        poly_mul_inp(&item, omega);
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

Poly* static construct_spline(const double* xarr, const double* yarr, int i, const Matrix* M)
{
    Poly* res = NULL;
    double h = xarr[i+1] - xarr[i];
    Poly* tmp = create_poly(-1, 1);
    tmp->next = create_poly(xarr[i+1], 0);
    poly_pow_inp(&tmp, 3);
    poly_mul_cons_inp(&tmp, M->mem[i][0] / (6 * h));
    poly_add_inp(&res, tmp);
    destroy_poly(&tmp);

    tmp = create_poly(1, 1);
    tmp->next = (-xarr[i], 0);
    poly_pow_inp(&tmp, 3);
    poly_mul_cons_inp(&tmp, M->mem[i+1][0] / (6 * h));
    poly_add_inp(&res, tmp);
    destroy_poly(&tmp);

    tmp = create_poly(-1, 1);
    tmp->next = create_poly(xarr[i+1], 0);
    poly_mul_cons_inp(&tmp, (yarr[i] - M->mem[i][0] * h * h / 6) / h);
    poly_add_inp(&res, tmp);
    destroy_poly(&tmp);

    tmp = create_poly(1, 1);
    tmp->next = (-xarr[i], 0);
    poly_mul_cons_inp(&tmp, (yarr[i+1] - M->mem[i+1][0] * h * h / 6) / h);
    poly_add_inp(&res, tmp);
    destroy_poly(&tmp);
    return res;
}   

/**
 *@type: 0, specify derivative at section boundary;
         1, specify second derivative at section boundary
         2, periodic function
 */
Poly** spline(const double* xarr, const double* yarr, int n, double a, double b, int type)
{
    int i, j;
    double lambda, u;
    Poly** res;
    Matrix* M;
    Matrix* d = create_vector(n + 1);
    Matrix* A = create_matrix_n(n + 1);
    
    //init A
    for (i = 0; i < n + 1; i++)
        A->mem[i][i] = 2;
    
    for (i = 1; i < n; i++) {
        lambda = (xarr[i+1] - xarr[i]) / (xarr[i+1] - xarr[i-1]);
        A->mem[i][i+1] = lambda;
        u = 1 - lambda;
        A->mem[i][i-1] = u;
    }

    //init d
    for (i = 1; i < n; i++) {
        double d_tmp = diff_quot(xarr, yarr, i, i+1) - diff_quot(xarr, yarr, i-1, i);
        d_tmp *= 6 / (xarr[i+1] - xarr[i-1]);
        d[i]  = d_tmp;
    }

    //init bundary conditions
    if (0 == type) {
        lambda = 1;
        A->mem[0][1] = lambda;
        u = 1;
        A->mem[n][n-1] = u;
        d->mem[0][0] = 6 * (diff_quot(xarr, yarr, 0, 1) - a) / (xarr[1] - xarr[0]);
        d->mem[n][0] = 6 * (b - diff_quot(xarr, yarr, n-1, n)) / (xarr[n] - xarr[n-1]);
    } else if (1 == type) {
        lambda = 0, u = 0;
        A->mem[0][1] = lambda;
        A->mem[n][n-1] = u;
        d->mem[0][0] = 2 * a;
        d->mem[n][0] = 2 * b;
    } else if (2 == type) {

    } else {
        printf("error: bad type = !\n", type);
        res = NULL;
        goto RETURN;
    }

    res = (Poly**)malloc(sizeof(Poly*) * n);
    memset((void*), 0, sizeof(Poly*) * n);
    M = chasing_method(A, d, false);
    for (i = 0; i < n; i++) {
        res[i] = construct_spline(xarr, yarr, i, M);
    }
    destroy_Matrix(&M);
RETURN:
    destroy_Matrix(&A);
    destroy_Matrix(&d);
    return res;
}




