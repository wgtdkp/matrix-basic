/**
description: calculate norm of vector and matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "norm.h"
#include "eigenvalue.h"

double norm(Matrix* M, unsigned int p, bool step)
{
    double ret = .0;
    int i, j;
    if(M->n == 1) {  //M is a vector
        if(INF == p) {
            for(i = 0; i < M->m; i++)
                ret = MAX(ret, fabs(M->mem[i][0]));
        } else {
            for(i = 0, ret = 0; i < M->m; i++)
                ret += pow(fabs(M->mem[i][0]), p);
            ret = pow(ret, 1.0 / p);
        }
    } else {
        assert(M->m == M->n);
        assert(INF == p || 1 == p || 2 == p);
        if(INF == p) {
            for(i = 0; i < M->m; i++) {
                double row_sum;
                for(j = 0, row_sum = 0; j < M->n; j++)
                    row_sum += fabs(M->mem[i][j]);
                ret = MAX(ret, row_sum);
            }
        } else if(1 == p) {
            for(j = 0; j < M->n; j++) {
                double col_sum;
                for(i = 0, col_sum = 0; i < M->m; i++)
                    col_sum += fabs(M->mem[i][j]);
                ret = MAX(ret, col_sum);
            }
        } else if(2 == p) {
            Matrix* MT = transpose(M);
            Matrix* MTM = mul(MT, M);
            ret = pow_method(MTM, false);
            ret = pow(ret, 0.5);
            destroy_matrix(&MT);
            destroy_matrix(&MTM);
        } else if(3 == p) { //Frobenius norm
            for(i = 0; i < M->m; i++)
                for(j = 0; j < M->n; j++)
                    ret += M->mem[i][j] * M->mem[i][j];
            ret = pow(ret, 0.5);
        }
    }
    return ret;
}

bool is_norm_similar(const Matrix* lhs, const Matrix* rhs, int p, double delta)
{
    Matrix* m_sub = sub_inp(lhs, rhs);
    bool res = norm(m_sub, p, false) < delta;
    destroy_matrix(&m_sub);
    return res;
}