#include "matrix.h"
#include "poly.h"
#include "direct.h"
#include "integration.h"
#include "approximation.h"
#include <math.h>

Poly* legendre(size_t index)
{
    size_t i;
    Poly* res = NULL;
    poly_add_inp(&res, create_poly(1, 2));
    poly_add_inp(&res, create_poly(-1, 0));
    poly_pow_inp(&res, index);

    for (i = 0; i < index; i++)
        ;//poly_diff_inp(&res);

    poly_mul_cons_inp(&res, 1.0 / pow(2, index));
    poly_mul_cons_inp(&res, 1.0 / factorial(index));
    return res;
}

static double wf_power(double x, double* coeff)
{
    return pow(x, coeff[0]);
}

/*
 *description: the direct method of 
 *   best square approximation.
*/
Poly* bsa_direct(func f, double a, double b, size_t n)
{
    size_t i, j;
    Poly* res = NULL;
    Matrix* H = create_matrix_n(n + 1);
    Matrix* D = create_vector(n + 1);
    Matrix* A;
    for (i = 0; i < H->m; i++) {
        for (j = 0; j < H->n; j++) {
            double h = pow(b, i + j + 1) - pow(a, i + j + 1);
            H->mem[i][j] = h / (i + j + 1);
        }
    }
    for (i = 0; i < D->m; i++)
        D->mem[i][0] = integration(f, a, b, 1e-6);
    A = me_gauss_elim(H, D, false);
    
    for (i = 0; i < A->m; i++)
        poly_add_inp(&res, create_poly(A->mem[i][0], i));
    
    destroy_matrix(&H);
    destroy_matrix(&D);
    destroy_matrix(&A);
    return res;
}

Poly* bsa_legendre(func f, double a, double b, size_t n)
{
    return NULL;
}