#ifndef _EIGENVALUE_H_
#define _EIGENVALUE_H_

#include "matrix.h"
#include "nonlinear.h"

double pow_method(Matrix* A, double delta, bool step);
Matrix* qr_decomp_H(Matrix* A, bool step);
Matrix* qr_decomp_G(Matrix* A, bool step);
void to_hessenberg(Matrix* A, bool step);
void qr_method(Matrix* A, double delta, bool step);
int eigenvalues(Complex* eigenvals, Matrix* A);

#endif