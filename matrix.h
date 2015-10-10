#ifndef _MATRIX_H_
#define _MATRIX_H_

#define DOUBLE_EQUAL(x, y)  (x - y < 1e-6 && x - y > -1e-6)

typedef enum{
    false = 0,
    true = 1
} bool;

extern int print_steps;

typedef struct {
	int m;
	int n;
	double** mem;
} Matrix;

//创建n阶方阵
Matrix* create_matrix_n(int n);

//创建mxn阶方阵
Matrix* create_matrix(int m_, int n_);

//创建n阶单位阵
Matrix* create_eye(int n);

//销毁方阵
void destroy_matrix(Matrix* M);

//矩阵相乘
Matrix* mul(Matrix* A, Matrix* B);

//计算矩阵的行列式
double det(Matrix* M);

//复制矩阵
Matrix* copy(Matrix* M);

//浅复制
Matrix* shallow_copy(Matrix* M);

//打印矩阵
void print_matrix(Matrix* M);

//高斯消去法解线性方程组 AX=B;
Matrix* gauss_elim(Matrix* A, Matrix* B);

//列主元消去法
Matrix* me_gauss_elim(Matrix* A, Matrix* B);

//矩阵的三角分解法(杜利特尔分解)
Matrix* tri_decomp(Matrix* A, Matrix* B);

//选主原的矩阵三角分解法
Matrix* me_tri_decomp(Matrix* A, Matrix* B);

//A为上三角形时，求解线性方程组
Matrix* ru_tri_solv(Matrix* A, Matrix* B);

//A为下三角时， 求解线性方程组
Matrix* lb_tri_solv(Matrix* A, Matrix* B);

Matrix* cholesky_decomp(Matrix* A, Matrix* B);

Matrix* en_cholesky_decomp(Matrix* A, Matrix* B);

Matrix* chasing_method(Matrix* A, Matrix* B);

Matrix* transpose(Matrix* M);


#endif