#ifndef _ITERATIVE_H_
#define _ITERATIVE_H_

Matrix* jacobi_iter(Matrix* A, Matrix* B, double delta, bool step);

Matrix* gauss_seidel_iter(Matrix* A, Matrix* B, double delta, bool step);

Matrix* sor(Matrix* A, Matrix* B, double w, double delta, bool step);


#endif