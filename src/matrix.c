#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>

int step = 0;
static double check_ordered_main_subdet(Matrix* M, bool* res, Comp checker, double x);

/**
创建n阶方阵
*/
Matrix* create_matrix_n(int n) {
    return create_matrix(n, n);
}

/**
创建mxn阶矩阵
m_:矩阵的行数
n_:矩阵的列数
*/
Matrix* create_matrix(int m_, int n_) {
    assert(m_ > 0 && n_ > 0);

    int i;
    Matrix* M = (Matrix*)malloc(sizeof(Matrix));
    if(NULL == M) return NULL;
    M->m = m_;
    M->n = n_;

    //请注意这里实际上是一个大的一维空间
    M->alloc = (double*)malloc(sizeof(double) * M->m * M->n);
    if(NULL == M->alloc) {
        free(M);
        return NULL;
    }
    memset(M->alloc, 0, sizeof(double) * M->m * M->n);

    M->mem = (double**)malloc(sizeof(double*) * M->m);
    if(NULL == M->mem) {
        free(M->alloc);
        free(M);
        return NULL;
    }
    for(i = 0; i < M->m; i++)
        M->mem[i] = M->alloc + M->n * i;

    return M;
}

/**
创建n阶单位阵
*/
Matrix* create_eye(int n) {
    int i;
	Matrix* M = create_matrix_n(n);
	for (i = 0; i < M->m; i++)
		M->mem[i][i] = 1;
	return M;
}

/**
销毁方阵M
*/
void destroy_matrix(Matrix* M) {
    if(NULL == M) return;
    free(M->alloc);
    free(M->mem);
    free(M);
}

/**
复制矩阵, 拷贝当前矩阵M
返回：新生成的矩阵的拷贝
*/
Matrix* copy(Matrix* M) {
    Matrix* replica;
    replica = create_matrix(M->m, M->n);
    if(NULL == replica)
        return NULL;
	memcpy(replica->mem[0], M->mem[0],\
        sizeof(double) * replica->m * replica->n);
    return replica;
}

/**
浅复制, 副本与原矩阵共享同一块memory空间
*/
Matrix* shallow_copy(Matrix* M) {
	return NULL;// return create_matrix_mem(M->m, M->n, M->mem);
}

void fill_matrix(Matrix* M, double x)
{
    int i, j;
    for(i = 0; i < M->m; i++)
        for(j = 0; j < M->n; j++)
            M->mem[i][j] = x;
}

/**
打印矩阵
*/
void print_matrix(Matrix* M) {
    int i, j;
    //printf("%d x %d\n", M->m, M->n);
    for(i = 0; i < M->m; i++) {
        printf("[");
        for(j = 0; j < M->n; j++)
            printf("%lf\t", M->mem[i][j]);
        printf("]\n");
    }
}

/**
矩阵相乘
*/
/*
Matrix* mul(Matrix* A, Matrix* B) {
	assert(A->n == B->m);

	Matrix* ret = create_matrix(A->m, B->n);
	int i, j, k;
	for (i = 0; i < ret->m; i++) {
		for (j = 0; j < ret->n; j++) {
			for (ret->mem[i][j] = 0, k = 0; k < A->n; k++)
				ret->mem[i][j] += A->mem[i][k] * B->mem[k][j];
        }
    }
	return ret;
}
*/

/*
this improved matrix mulplication is better, cause it 
has better cache hit.
*/
Matrix* mul(Matrix* A, Matrix* B) {
    Matrix* ret;
    int i, j, k;
    if(A->n != B->m)
        return NULL;
    ret = create_matrix(A->m, B->n);
    
    for (i = 0; i < ret->m; i++) {
        for (k = 0; k < A->n; k++) {
            int r = A->mem[i][k];
            for (j = 0; j < ret->n; j++)
                ret->mem[i][j] += r * B->mem[k][j];
        }
    }
    return ret;
}

/**
计算矩阵的行列式
1.改进的算法更快
2.当数值较大时，计算存在误差
*/
double det(Matrix* M) {
    assert(M->m == M->n);

    int i, j;
    double ret, cofactor;
	double* memi;
    if(1 == M->n)
        return M->mem[0][0];

    ret = 0, memi = M->mem[M->m - 1];
    for(i = M->m - 1, j = M->n - 1; i >= 0; i--) {
        //int aij = M->mem[i][j];
        //double* memi = M->mem[i];
        //int k;
        //for(k = i; k < M->m - 1; k++)
        //    M->mem[k] = M->mem[k+1];
		//aij = M->mem[i][j];
		swap_dp(&M->mem[i], &memi);
        M->m--;
        M->n--;

        cofactor = det(M);
		cofactor *= ((i + j) & 1) == 1 ? -1 : 1; // Aij = (-1)i+jMij
        ret += memi[j] * cofactor;

        //恢复M矩阵
        M->m++;
        M->n++;
        //for(k = M->m - 1; k > i; k--)
        //    M->mem[k] = M->mem[k-1];
        //M->mem[k] = memi;
    }

	//恢复M矩阵
	for (i = M->m - 1; i > 0; i--)
		M->mem[i] = M->mem[i - 1];
	M->mem[0] = memi;
    return ret;
}




/**
the wrapper of check_ordered_main_subdet()
*/
bool is_ordered_main_subdet(Matrix* M, Comp checker, double x) {
    assert(M->m == M->n);

    bool res = true;
    double det = check_ordered_main_subdet(M, &res, checker, x);
    return res && (*checker)(det, x); 
}

/**
check if all ordered main subdeterminant of M if not zero.
@M: the matrix
@res: the checking result
@comp: the checking method
@x: the parameter for the checking method
return: det M
*/
static double check_ordered_main_subdet(Matrix* M, bool* res, Comp checker, double x) 
{
    int i, j;
    char sign;
    double  ret, *memi;
    if(1 == M->n)
        return M->mem[0][0];
    if(false == *res)
        return 0;   //we can return anything, it is not important now

    ret = 0, memi = M->mem[M->m - 1];
    for(i = M->m - 1, j = M->n - 1; i >= 0; i--) {
        double cofactor;

        swap_dp(&M->mem[i], &memi);
        M->m--;
        M->n--;

        cofactor = check_ordered_main_subdet(M, res, checker, x);
        sign = ((i + j) & 1) == 1 ? -1 : 1; // Aij = (-1)i+jMij
        ret += memi[j] * sign * cofactor;

        //恢复M矩阵
        M->m++;
        M->n++;

        if(i == M->m - 1 && !(*checker)(cofactor, x)) {
            //we can not call return here, it could dirty matrix M
            *res = false; 
        }
    }

    //恢复M矩阵
    for (i = M->m - 1; i > 0; i--)
        M->mem[i] = M->mem[i - 1];
    M->mem[0] = memi;

    return ret;
}

Matrix* transpose(Matrix* M) {
    int i, j;
    Matrix* T = create_matrix(M->n, M->m);
    for(i = 0; i < M->m; i++)
        for(j = 0; j < M->n; j++)
            T->mem[j][i] = M->mem[i][j];
    return T;
}

bool is_symmetrical(Matrix* M) {
    int i, j;
    if(M->m != M->n) return false;
    for(i = 0; i < M->m; i++)
        for(j = i; j < M->n; j++)
            if(!DOUBLE_EQUAL(M->mem[i][j], M->mem[j][i]))
                return false;
    return true;
}

/**
description: check if the two matrixes satisfy:
    M = coeff * N (coeff is an constant)
    If true, the coeff param stores the scale
*/
bool is_const_similar(Matrix* M, Matrix* N, double* coeff)
{
    int i, j;
    if(M->m != N->m || M->n != N->n)
        return false;
    
    for(i = 0; i < M->m; i++)
        for(j = 0; j < M->n; j++) {
            if(!DOUBLE_EQUAL(M->mem[i][j], .0) \
                && !DOUBLE_EQUAL(M->mem[i][j], .0)) {
                *coeff = M->mem[i][j] / N->mem[i][j];
                goto GET;
            }
        }
    return false;
GET:
    for(i = 0; i < M->m; i++)
        for(j = 0; j < M->n; j++) {
            if(!DOUBLE_EQUAL((M->mem[i][j] / N->mem[i][j]) / (*coeff), 1.0))
                return false;
        }
    return true;
}

void swap_matrix(Matrix* M, Matrix* N)
{
    assert(M->m == N->m && M->n == N->n);
    Matrix tmp;
    tmp = *M;
    *M = *N;
    *N = tmp;
}

void swap_dp(double** a, double** b) {
    double* tmp = *a;
    *a = *b;
    *b = tmp;
}

void swap_d(double* a, double* b) {
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

/**
description: get the inverse matrix using gauss elimination. 
return: if inverse matrix does not exist, return NULL.
*/
Matrix* inverse(Matrix* A)
{
    assert(A->m == A->n);
    int i, j, k;
    Matrix* inv = create_eye(A->m);
    for(i = 0; i < A->m - 1; i++) {
        double factor, max_col_value;
        int max_col;
        
        //选择头元素最大的行作为当前行
        max_col_value = A->mem[i][i], max_col = i;
        for(k = i + 1; k < A->m; k++) {
            if(abs(A->mem[k][i]) > max_col_value) {
                max_col_value = abs(A->mem[k][i]);
                max_col = k;
            }
        }
        swap_row(A, i, max_col);
        swap_row(inv, i, max_col);
        if(eq(A->mem[i][i], 0.0)) {
            goto NO_INV;
        }
        
        for(k = i + 1; k < A->m; k++) {
            factor = A->mem[k][i] / A->mem[i][i];
            for(j = i; j < A->n; j++)
                A->mem[k][j] -= A->mem[i][j] * factor;
            for(j = 0; j < A->n; j++)
                inv->mem[k][j] -= inv->mem[i][j] * factor;
        }
    }

    for(i = A->m - 1; i >= 0; i--) {
        for(j = i - 1; j >= 0; j--) {
            double factor = A->mem[j][i] / A->mem[i][i];
            for(k = 0; k < inv->n; k++)
                inv->mem[j][k] -= inv->mem[i][k] * factor;
        }
        for(k = 0; k < inv->n; k++)
            inv->mem[i][k] /= A->mem[i][i];
    }
    return inv;

NO_INV:
    destroy_matrix(inv);
    printf("error: the matrix is a singular "
        "matrix, no inervse matrix!\n");
    return NULL;
}
