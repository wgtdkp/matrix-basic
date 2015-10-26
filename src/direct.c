/**
description: solve linear equation with
    direct methods.
-- Gauss elimination
-- main element Gauss elimination
-- triangle decomposition
-- main element triangle decomposition
-- cholesky decomposition
-- enhanced cholesky decomposition
-- chasing method
*/

#include "matrix.h"
#include "direct.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

//static Matrix* gen_p(int ip[], int n);

/**
高斯消去法解线性方程组 AX=B;
返回：线性方程组的解
*/
Matrix* gauss_elim(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));

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
            for(j = i; j < A->n; j++) {
                A->mem[k][j] -= A->mem[i][j] * factor;
            }
            B->mem[k][0] -= B->mem[i][0] * factor;
        }
        if(false != step) {
            printf("=== step %d ===\n", i + 1);
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
Matrix* me_gauss_elim(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));

    int i, j, k;
    for(i = 0; i < A->m - 1; i++) {
        double factor, max_col_value;
        int max_col;
        
        //选择头元素最大的行作为当前行
        max_col_value = fabs(A->mem[i][i]), max_col = i;
        for(k = i + 1; k < A->m; k++) {
            if(fabs(A->mem[k][i]) > max_col_value) {
                max_col_value = fabs(A->mem[k][i]);
                max_col = k;
            }
        }
        printf("i: %d, max_col: %d\n", i, max_col);
        swap_dp(&A->mem[i], &A->mem[max_col]); //交换当前行与最大行
        swap_dp(&B->mem[i], &B->mem[max_col]);

        if(DOUBLE_EQUAL(A->mem[i][i], 0.0)) {
            printf("error: me gauss elimination failed, "
                "because A[%d][%d] is zero\n", i + 1, i + 1);
            return NULL;
        }

        for(k = i + 1; k < A->m; k++) {
            factor = A->mem[k][i] / A->mem[i][i];
            for(j = i; j < A->n; j++)
                A->mem[k][j] -= A->mem[i][j] * factor;
            B->mem[k][0] -= B->mem[i][0] * factor;
        }
        if(0 != step) {
            printf("=== step %d ===\n", i + 1);
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
Matrix* tri_decomp(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));
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
                destroy_matrix(&L); destroy_matrix(&U);
                printf("error: triangular decomposition failed, \
                    because U[%d][%d] is zero\n", r + 1, r + 1);
                return NULL;
            }
            L->mem[i][r] = (A->mem[i][r] - s) / U->mem[r][r];
        }
        if (0 != step) {
            printf("=== step %d ===\n", i + 1);
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
    destroy_matrix(&L); 
    destroy_matrix(&U); 
    destroy_matrix(&Y);
    return X;
}

/**
选主原的矩阵三角分解法
*/
Matrix* me_tri_decomp(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));
    assert(is_ordered_main_subdet(A, &ne, 0));

    int i, k, r;
    Matrix* L = create_matrix(A->m, A->n);
    Matrix* U = create_matrix(A->m, A->n);
    Matrix* P = create_eye(A->m);
    for(r = 0; r < A->m; r++) {
        
        //选主元
        double max_s = 0; // make compiler happy
        int max_s_index = 0; // make compiler happy
        for(i = r; i < A->n; i++) {
            double s;
            for(s = 0, k = 0; k < r; k++)
                s += L->mem[i][k] * U->mem[k][r];
            s = fabs(A->mem[i][r] - s);
            if(i == r || s > max_s) {   //i == r, init max_s and max_s_index
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
        //ip[r] = max_s_index;  //我们可以省略ip数组
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
                s += L->mem[i][k] * U->mem[k][r];   //这里有重复计算
            if(DOUBLE_EQUAL(U->mem[r][r], 0)) {
                destroy_matrix(&L); destroy_matrix(&U); destroy_matrix(&P);
                printf("error: triangular decomposition failed, \
                    because U[%d][%d] is zero\n", r + 1, r + 1);
                return NULL;
            }
            L->mem[i][r] = (A->mem[i][r] - s) / U->mem[r][r];
        }
        if (0 != step) {
            printf("=== step %d ===\n", i + 1);
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
    destroy_matrix(&L); 
    destroy_matrix(&U); 
    destroy_matrix(&Y); 
    destroy_matrix(&P);
    return X;

}

/**
楚列斯基分解
*/
Matrix* cholesky_decomp(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));
    assert(is_symmetrical(A) && is_ordered_main_subdet(A, &gt, 0));
    int i, j, k;
    Matrix *X, *Y;
    Matrix *U, *L;
    //检查系数矩阵A是否为对称正定矩阵
    ///if(!is_symmetrical(A) || !is_ordered_main_subdet(A, &gt, 0)) {
    //    printf("error: the matrix is not a symmetric positive definite.\n");
    //    return NULL;
    //}
    L = create_matrix(A->m, A->n);
    for(i = 0; i < A->m; i++) {
        for(j = 0; j <= i; j++) {
            double s;
            for(k = 0, s = 0; k < j; k++)
                s += L->mem[i][k] * L->mem[j][k];
            if(j == i) {
                if(A->mem[i][j] - s < 0) {
                    destroy_matrix(&L);
                    printf("error: cholesky dcompsition failed, " 
                        "because A->mem[%d][%d] = sqrt(%.8lf)\n", 
                        i + 1, i + 1, A->mem[i][j] - s);
                    return NULL;
                }
                L->mem[i][j] = sqrt(A->mem[i][j] - s);
            }
            else
                L->mem[i][j] = (A->mem[i][j] - s) / L->mem[j][j];
        }
        if (0 != step) {
            printf("=== step %d ===\n", i + 1);
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
    destroy_matrix(&L); 
    destroy_matrix(&U); 
    destroy_matrix(&Y);
    return X;
}

/**
enhanced cholesky decomposition, there is no sqrt()
*/
Matrix* en_cholesky_decomp(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));
    assert(is_symmetrical(A) && is_ordered_main_subdet(A, &gt, 0));
    int i, j, k;
    Matrix *X, *Y, *U;
    Matrix *L, *D, *LD; 
    //检查系数矩阵A是否为对称正定矩阵
    //if(!is_symmetrical(A) || !is_ordered_main_subdet(A, &gt, 0)) {
    //    printf("error: the matrix is not a symmetric positive definite.\n");
    //    return NULL;
    //}
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
        
        if (0 != step) {
            printf("=== step %d ===\n", i + 1);
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
    destroy_matrix(&L); 
    destroy_matrix(&U); 
    destroy_matrix(&Y); 
    destroy_matrix(&D); 
    destroy_matrix(&LD);
    return X;
}

Matrix* chasing_method(Matrix* A, Matrix* B, bool step) {
    assert(IS_EQUATION(A, B));
    assert(is_diagnoal_dominance_3(A));
    int i;
    Matrix *L, *U;
    Matrix *X, *Y;
    //if(!is_diagnoal_dominance_3(A)) {
    //    printf("error: the matrix is not a 3 diagnoal dominance.\n");
    //    return NULL;
    //}
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
        if(0 != step) {
            printf("it is so simple that no steps to show you\n");
        }
    }

    printf("the L matrix:\n");
    print_matrix(L);
    printf("the U matrix:\n");
    print_matrix(U);

    Y = lb_tri_solv(L, B);
    X = ru_tri_solv(U, Y);
    destroy_matrix(&L); 
    destroy_matrix(&U); 
    destroy_matrix(&Y);
    return X;
}

/**
求解A为上三角形的线性方程组
*/
Matrix* ru_tri_solv(Matrix* A, Matrix* B) {
    assert(IS_EQUATION(A, B));
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
    assert(IS_EQUATION(A, B));
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

/*
static Matrix* gen_p(int ip[], int n) {
    int i;
    Matrix* P = create_eye(n);
    for (i = 0; i < n; i++)
        swap_dp(&P->mem[i], &P->mem[ip[i]]);
    return P;
}
*/

bool is_diagnoal_dominance_3(Matrix* M) {
    //assert(IS_SQUARE(M));
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