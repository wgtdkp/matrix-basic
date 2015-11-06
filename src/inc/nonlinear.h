#ifndef _NONLINEAR_H_
#define _NONLINEAR_H_

#include <stdio.h>
#include <math.h>
#include "poly.h"


#define COMPLEX_EQUAL_DELTA(x, y, delta)   (DOUBLE_EQUAL_DELTA((x).real, (y).real, (delta)) \
                            && DOUBLE_EQUAL_DELTA((x).imag, (y).imag, (delta)))

typedef struct {
    double real;
    double imag;
} Complex;

static inline Complex* create_complex(double real, double imag)
{
    Complex* comp;
    comp = (Complex*)malloc(sizeof(Complex));
    if (NULL == comp)
        return NULL;
    comp->real = real;
    comp->imag = imag;
    return comp;
}

static inline void print_complex(const Complex* comp)
{
    printf(comp->imag < 0 ? "%.8lf - %.8lfi" : "%.8lf + %.8lfi", 
        comp->real, fabs(comp->imag));
}

void quadratic(Complex* roots, const Poly* equaltion);

#endif