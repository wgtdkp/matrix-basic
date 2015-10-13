#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>

int print_steps = 0;

typedef bool (*Comp)(double, double);

static Matrix* gen_p(int ip[], int n);
static void swap_dp(double** a, double** b);
static void swap_d(double* a, double* b);
static bool is_symmetrical(Matrix* M);
static bool is_ordered_main_subdet(Matrix* M, Comp checker, double x);
static double check_ordered_main_subdet(Matrix* M, bool* res, Comp checker, double x);
static bool is_diagnoal_dominance_3(Matrix* M);

static inline bool gt(double x, double y) {
    return x > y;
}

static inline bool eq(double x, double y) {
    return DOUBLE_EQUAL(x, y);
}

static inline bool ne(double x, double y) {
    return !eq(x, y);
}

static inline bool lt(double x, double y) {
    return !gt(x, y) && ne(x, y);
}

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
	Matrix* M = create_matrix_n(n);
	int i;
	for (i = 0; i < M->m; i++)
		M->mem[i][i] = 1;
	return M;
}

/**
销毁方阵M
*/
void destroy_matrix(Matrix* M) {
    free(M->alloc);
    free(M->mem);
    free(M);
}

/**
复制矩阵, 拷贝当前矩阵M
返回：新生成的矩阵的拷贝
*/
Matrix* copy(Matrix* M) {
    Matrix* replica = create_matrix(M->m, M->n);
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

/**
打印矩阵
*/
void print_matrix(Matrix* M) {
    int i, j;
   // printf("%d x %d\n", M->m, M->n);
    for(i = 0; i < M->m; i++) {
        for(j = 0; j < M->n; j++)
            printf("%lf\t", M->mem[i][j]);
        printf("\n");
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
    assert(A->n == B->m);

    Matrix* ret = create_matrix(A->m, B->n);
    int i, j, k;
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
求解A为上三角形的线性方程组
*/
Matrix* ru_tri_solv(Matrix* A, Matrix* B) {
	int i, j;
	Matrix* X = create_matrix(B->m, 1);

	//这里没有将计算的结果存储在B中，虽然这样确实可以节约空间
	for (i = A->m - 1; i >= 0; i--) {
		double s;
		for (s = 0, j = A->n - 1; j > i; j--)
			s += A->mem[i][j] * X->mem[j][0];
		X->mem[i][0] = (B->mem[i][0] - s) / A->mem[i][i];
	}
	return X;
}

/**
求解A为下三角形的线性方程组
*/
Matrix* lb_tri_solv(Matrix* A, Matrix* B) {
	int i, j;
	Matrix* X = create_matrix(B->m, 1);

	//这里没有将计算的结果存储在B中，虽然这样确实可以节约空间
	for (i = 0; i < A->m; i++) {
		double s;
		for (s = 0, j = 0; j < i; j++)
			s += A->mem[i][j] * X->mem[j][0];
		X->mem[i][0] = (B->mem[i][0] - s) / A->mem[i][i];
	}
	return X;
}

/**
高斯消去法解线性方程组 AX=B;
返回：线性方程组的解
*/
Matrix* gauss_elim(Matrix* A, Matrix* B) {
    assert(A->m == A->n);

    int i, j, k;
    for(i = 0; i < A->m - 1; i++) {
        double factor;
        if(DOUBLE_EQUAL(A->mem[i][i], 0.0)) {
            printf("error: gauss elimination failed, \
                because A[%d][%d] is zero\n", i + 1, i + 1);
            return NULL;
        }
        for(k = i + 1; k < A->m; k++) {
            factor = A->mem[k][i] / A->mem[i][i];
            for(j = i; j < A->n; j++)
                A->mem[k][j] -= A->mem[i][j] * factor;
            B->mem[k][0] -= B->mem[i][0] * factor;
        }
        if(0 != print_steps) {
            printf("step:%d\n", i + 1);
            printf("A:\n");
            print_matrix(A);
            printf("B:\n");
            print_matrix(B);
            printf("\n");
        }
    }
	return ru_tri_solv(A, B);
}

/**
列主元消去法
*/
Matrix* me_gauss_elim(Matrix* A, Matrix* B) {
    assert(A->m == A->n);

    int i, j, k;
    for(i = 0; i < A->m - 1; i++) {
        double factor, max_col_value;
        int max_col;
        if(DOUBLE_EQUAL(A->mem[i][i], 0.0)) {
            printf("error: gauss elimination failed, \
                because A[%d][%d] is zero\n", i + 1, i + 1);
            return NULL;
        }
        //选择头元素最大的行作为当前行
        max_col_value = A->mem[i][i], max_col = i;
        for(k = i + 1; k < A->m; k++) {
            if(abs(A->mem[k][i]) > max_col_value) {
                max_col_value = abs(A->mem[k][i]);
                max_col = k;
            }
        }
        swap_dp(&A->mem[i], &A->mem[max_col]); //交换当前行与最大行
		swap_dp(&B->mem[i], &B->mem[max_col]);

        for(k = i + 1; k < A->m; k++) {
            factor = A->mem[k][i] / A->mem[i][i];
            for(j = i; j < A->n; j++)
                A->mem[k][j] -= A->mem[i][j] * factor;
            B->mem[k][0] -= B->mem[i][0] * factor;
        }
        if(0 != print_steps) {
            printf("step:%d\n", i + 1);
            printf("A:\n");
            print_matrix(A);
            printf("B:\n");
            print_matrix(B);
            printf("\n");
        }
    }
	return ru_tri_solv(A, B);
}


/**
矩阵的三角分解法(杜利特尔分解)
*/
Matrix* tri_decomp(Matrix* A, Matrix* B) {

    assert(is_ordered_main_subdet(A, &ne, 0));
    int i, k, r;
	Matrix* L = create_matrix(A->m, A->n);
	Matrix* U = create_matrix(A->m, A->n);
    for(r = 0; r < A->m; r++) {
        for(i = r; i < A->n; i++) {
            double s;
            for(s = 0, k = 0; k < r; k++)
                s += L->mem[r][k] * U->mem[k][i];
            U->mem[r][i] = A->mem[r][i] - s;
            for(s = 0, k = 0; k < r; k++)
                s += L->mem[i][k] * U->mem[k][r];
            if(DOUBLE_EQUAL(U->mem[r][r], 0)) {
				destroy_matrix(L); destroy_matrix(U);
                printf("error: triangular decomposition failed, \
                    because U[%d][%d] is zero\n", r + 1, r + 1);
                return NULL;
            }
            L->mem[i][r] = (A->mem[i][r] - s) / U->mem[r][r];
        }
		if (0 != print_steps) {
			printf("step:%d\n", r + 1);
			printf("L:\n");
			print_matrix(L);
			printf("U:\n");
			print_matrix(U);
			printf("\n");
		}
    }

	printf("the L matrix:\n");
	print_matrix(L);
	printf("the U matrix:\n");
	print_matrix(U);

	//解算方程组
	Matrix* Y = lb_tri_solv(L, B);
	Matrix* X = ru_tri_solv(U, Y);
	destroy_matrix(L); 
    destroy_matrix(U); 
    destroy_matrix(Y);
	return X;
}

/**
选主原的矩阵三角分解法
*/
Matrix* me_tri_decomp(Matrix* A, Matrix* B) {
    assert(is_ordered_main_subdet(A, &ne, 0));

    int i, k, r;
	Matrix* L = create_matrix(A->m, A->n);
	Matrix* U = create_matrix(A->m, A->n);
	Matrix* P = create_eye(A->m);
    for(r = 0; r < A->m; r++) {
        
        //选主元
        double max_s;
        int max_s_index;
        for(i = r; i < A->n; i++) {
            double s;
            for(s = 0, k = 0; k < r; k++)
                s += L->mem[i][k] * U->mem[k][r];
            s = fabs(A->mem[i][r] - s);
            if(i == r || s > max_s) {
                max_s = s;
                max_s_index = i;
            }
        }
        swap_dp(&A->mem[r], &A->mem[max_s_index]);

        /*
        if we store L in matrix A, we don't need to swap L[r] with L[max_s_index]
        and because we calculate U from left to right(not from top to bottom)
        we don't need to swap U[r] with U[max_s_index]
        */
		swap_dp(&L->mem[r], &L->mem[max_s_index]);	
		//swap_dp(&U->mem[r], &U->mem[max_s_index]);
        //ip[r] = max_s_index;	//我们可以省略ip数组
		swap_dp(&P->mem[r], &P->mem[max_s_index]); // this is ip designed for.
		swap_dp(&B->mem[r], &B->mem[max_s_index]);


        U->mem[r][r] = max_s;

        for(i = r; i < A->n; i++) {
            double s;
            for(s = 0, k = 0; k < r; k++)
                s += L->mem[r][k] * U->mem[k][i];
            U->mem[r][i] = A->mem[r][i] - s;
        }

        for(i = r; i < A->n; i++) {
            double s;
            for(s = 0, k = 0; k < r; k++)
                s += L->mem[i][k] * U->mem[k][r];	//这里有重复计算
            if(DOUBLE_EQUAL(U->mem[r][r], 0)) {
				destroy_matrix(L); destroy_matrix(U); destroy_matrix(P);
                printf("error: triangular decomposition failed, \
                    because U[%d][%d] is zero\n", r + 1, r + 1);
                return NULL;
            }
            L->mem[i][r] = (A->mem[i][r] - s) / U->mem[r][r];
        }
		if (0 != print_steps) {
			printf("step:%d\n", r + 1);
			printf("L:\n");
			print_matrix(L);
			printf("U:\n");
			print_matrix(U);
			printf("\n");
		}
    }
    
	printf("the P Matrix:\n");
	print_matrix(P);
	printf("the L matrix:\n");
	print_matrix(L);
	printf("the U matrix:\n");
	print_matrix(U);

	//解算方程组
	Matrix* Y = lb_tri_solv(L, B);
	Matrix* X = ru_tri_solv(U, Y);
	destroy_matrix(L); 
    destroy_matrix(U); 
    destroy_matrix(Y); 
    destroy_matrix(P);
	return X;

}

/**
楚列斯基分解
*/
Matrix* cholesky_decomp(Matrix* A, Matrix* B) {
    //检查系数矩阵A是否为对称正定矩阵
    assert(is_symmetrical(A) && is_ordered_main_subdet(A, &gt, 0));

    int i, j, k;
    Matrix *X, *Y;
    Matrix *U, *L = create_matrix(A->m, A->n);
    for(i = 0; i < A->m; i++) {
        for(j = 0; j <= i; j++) {
            double s;
            for(k = 0, s = 0; k < j; k++)
                s += L->mem[i][k] * L->mem[j][k];
            if(j == i) {
                if(A->mem[i][j] - s < 0) {
                    destroy_matrix(L);
                    printf("error: cholesky dcompsition failed, " 
                        "because A->mem[%d][%d] = sqrt(%lf)\n", 
                        i + 1, i + 1, A->mem[i][j] - s);
                    return NULL;
                }
                L->mem[i][j] = sqrt(A->mem[i][j] - s);
            }
            else
                L->mem[i][j] = (A->mem[i][j] - s) / L->mem[j][j];
        }
        if (0 != print_steps) {
            printf("step:%d\n", i + 1);
            printf("L:\n");
            print_matrix(L);            
            printf("\n");
        }
    }
    
    U = transpose(L);
    printf("the L matrix:\n");
    print_matrix(L);
    printf("the U matrix:\n");
    print_matrix(U);

    Y = lb_tri_solv(L, B);
    X = ru_tri_solv(U, Y);
    destroy_matrix(L); 
    destroy_matrix(U); 
    destroy_matrix(Y);
    return X;
}

/**
enhanced cholesky decomposition, there is no sqrt()
*/
Matrix* en_cholesky_decomp(Matrix* A, Matrix* B) {
    //ordered main subdets should be != 0
    assert(is_symmetrical(A) && is_ordered_main_subdet(A, &ne, 0));

    int i, j, k;
    Matrix *X, *Y, *U;
    Matrix *L, *D, *LD; 

    //we can use an array, but matrix is convenient
    L = create_matrix(A->m, A->n);
    D = create_matrix(A->m, A->n);
    for(i = 0; i < A->m; i++) {
        double s;
        for(j = 0; j < i; j++) {
            for(k = 0, s = 0; k < j; k++)
                s += L->mem[i][k] * D->mem[k][k] * L->mem[j][k];
            L->mem[i][j] = (A->mem[i][j] - s) / D->mem[j][j];
        }
        L->mem[i][i] = 1;
        for(k = 0, s = 0; k < j; k++)
            s += L->mem[i][k] * L->mem[i][k] * D->mem[k][k];
        D->mem[i][i] = A->mem[i][i] - s;
        
        if (0 != print_steps) {
            printf("step:%d\n", i + 1);
            printf("L:\n");
            print_matrix(L);            
            printf("\n");
        }
    }
    
    U = transpose(L);
    printf("the L matrix:\n");
    print_matrix(L);
    printf("the U matrix:\n");
    print_matrix(U);
    printf("the D matrix:\n");
    print_matrix(D);

    LD = mul(L, D);
    Y = lb_tri_solv(LD, B);
    X = ru_tri_solv(U, Y);
    destroy_matrix(L); 
    destroy_matrix(U); 
    destroy_matrix(Y); 
    destroy_matrix(D); 
    destroy_matrix(LD);
    return X;
}

static bool is_diagnoal_dominance_3(Matrix* M) {
    int i, j;
    if(M->m != M->n) return false;
    for(i = 0; i < M->m; i++) {
        if(i - 1 >= 0 && !(abs(M->mem[i][i]) > abs(M->mem[i][i - 1])))
            return false;
        if(i + 1 < M->n && !(abs(M->mem[i][i]) > abs(M->mem[i][i + 1])))
            return false;
    }
    for(i = 0; i < M->m; i++) {
        for(j = 0; j < M->n; j++) {
            if(abs(i - j) > 1 && !DOUBLE_EQUAL(M->mem[i][j], 0))
                return false; 
            if(abs(i - j) <= 1 && DOUBLE_EQUAL(M->mem[i][j], 0))
                return false;
        }
    }
    return true;
}

Matrix* chasing_method(Matrix* A, Matrix* B) {
    assert(is_diagnoal_dominance_3(A));

    int i;
    Matrix *L, *U;
    Matrix *X, *Y;
    L = create_matrix(A->m, A->n);
    U = create_matrix(A->m, A->n);

    /**
    alpha[i]: L[i][i]
    gamma[i]: L[i][i - 1]
    beta[i]:  U[i][i + 1]
    */
    L->mem[0][0] = A->mem[0][0];
    U->mem[0][0] = 1;
    U->mem[0][1] = A->mem[0][1] / A->mem[0][0];
    for(i = 1; i < A->m; i++) {
        L->mem[i][i - 1] = A->mem[i][i - 1];
        L->mem[i][i] = A->mem[i][i] - A->mem[i][i - 1] * U->mem[i - 1][i];
        if(i + 1 < A->m)    // try delete this line!
            U->mem[i][i + 1] = A->mem[i][i + 1] / L->mem[i][i];
        U->mem[i][i] = 1;
        if(0 != print_steps) {
            printf("it is so simple that no steps to show you\n");
        }
    }

    printf("the L matrix:\n");
    print_matrix(L);
    printf("the U matrix:\n");
    print_matrix(U);

    Y = lb_tri_solv(L, B);
    X = ru_tri_solv(U, Y);
    destroy_matrix(L); 
    destroy_matrix(U); 
    destroy_matrix(Y);
    return X;
}



/**
the wrapper of check_ordered_main_subdet()
*/
static bool is_ordered_main_subdet(Matrix* M, Comp checker, double x) {
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
static double check_ordered_main_subdet(
    Matrix* M, 
    bool* res, 
    Comp checker, 
    double x) 
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

static bool is_symmetrical(Matrix* M) {
    int i, j;
    if(M->m != M->n) return false;
    for(i = 0; i < M->m; i++)
        for(j = i; j < M->n; j++)
            if(M->mem[i][j] != M->mem[j][i])
                return false;
    return true;
}

static Matrix* gen_p(int ip[], int n) {
	int i;
    Matrix* P = create_eye(n);
	for (i = 0; i < n; i++)
		swap_dp(&P->mem[i], &P->mem[ip[i]]);
	return P;
}

static void swap_dp(double** a, double** b) {
    double* tmp = *a;
    *a = *b;
    *b = tmp;
}

static void swap_d(double* a, double* b) {
    double tmp = *a;
    *a = *b;
    *b = tmp;
}
