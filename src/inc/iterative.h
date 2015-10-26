#ifndef _ITERATIVE_H_
#define _ITERATIVE_H_

Matrix* jacobi_iter(Matrix* A, Matrix* B, double delta, bool step);

Matrix* gauss_seidel_iter(Matrix* A, Matrix* B, double delta, bool step);

Matrix* sor(Matrix* A, Matrix* B, double w, double delta, bool step);

bool is_jacobi_convergent(Matrix* A, bool step);

bool is_gauss_seidel_convergent(Matrix* A, bool step);

bool is_sor_convergent(Matrix* A, double w, bool step);

Matrix* sor_homework(Matrix* A, Matrix* B, double w, double delta, bool step);

#endif