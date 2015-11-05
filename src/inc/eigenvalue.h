#ifndef _EIGENVALUE_H_
#define _EIGENVALUE_H_

double pow_method(Matrix* A, double delta, bool step);
Matrix* qr_decomp_H(Matrix* A, bool step);
Matrix* qr_decomp_G(Matrix* A, bool step);
void to_hessenberg(Matrix* A, bool step);
void qr_method(Matrix* A, bool step);
#endif