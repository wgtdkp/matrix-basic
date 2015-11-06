#ifndef _NORM_H_
#define _NORM_H_

double norm(Matrix* M, unsigned int p, bool step);

/**
 *check if rhs matrix is approaching to lhs matrix,
 *if the norm of the diff matrix is small enough
 *@p: the param of function norm()
 *@delta: the maximum value of the norm
 */
bool static inline is_norm_similar(const Matrix* lhs, const Matrix* rhs, int p, double delta)
{
    Matrix* m_sub = sub(lhs, rhs);
    bool res = norm(m_sub, p, false) < delta;
    destroy_matrix(&m_sub);
    return res;
}
#endif