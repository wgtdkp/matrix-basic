#ifndef _DIRECT_H_
#define _DIRECT_H_

//高斯消去法解线性方程组 AX=B;
Matrix* gauss_elim(Matrix* A, Matrix* B);

//列主元消去法
Matrix* me_gauss_elim(Matrix* A, Matrix* B);

//矩阵的三角分解法(杜利特尔分解)
Matrix* tri_decomp(Matrix* A, Matrix* B);

//选主原的矩阵三角分解法
Matrix* me_tri_decomp(Matrix* A, Matrix* B);

Matrix* cholesky_decomp(Matrix* A, Matrix* B);

Matrix* en_cholesky_decomp(Matrix* A, Matrix* B);

Matrix* chasing_method(Matrix* A, Matrix* B);

//A为上三角形时，求解线性方程组
Matrix* ru_tri_solv(Matrix* A, Matrix* B);

//A为下三角时， 求解线性方程组
Matrix* lb_tri_solv(Matrix* A, Matrix* B);


#endif