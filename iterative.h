#ifndef _ITERATIVE_H_
#define _ITERATIVE_H_

Matrix* jacobi_iter(Matrix* A, Matrix* B, double delta);

Matrix* gauss_seidel_iter(Matrix* A, Matrix* B, double delta);

Matrix* sor(Matrix* A, Matrix* B, double w, double delta);


#endif