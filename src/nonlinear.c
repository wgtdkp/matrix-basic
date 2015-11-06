#include "nonlinear.h"
#include "poly.h"

void quadratic(Complex* roots, const Poly* equaltion)
{
    double abc[3], tmp, d;
    const Poly* p;
    if (NULL == roots || 2 != equaltion->index)
        return;
    for (p = equaltion; NULL != p; p = p->next)
        abc[p->index] = p->coeff;
    tmp = abc[1] * abc[1] - 4 * abc[2] * abc[0];
    roots[0].real = -abc[1] / 2;
    roots[1].real = -abc[1] / 2;
    d = sqrt(fabs(tmp));
    if (tmp > 0) {
        roots[0].real += d / 2;
        roots[1].real += -d / 2;
        roots[0].imag = roots[1].imag = 0;
    } else {
        roots[0].imag = d / 2;
        roots[1].imag = -d / 2;
    }
}
