#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include "poly.h"
#include "matrix.h"

typedef struct {
    double x;
    double y;
} Point;

typedef Point Pos;

Poly* vandermonde(const double* xarr, const double* yarr, int n);
Poly* lagrange(const double* xarr, const double* yarr, int n);
Poly* newton(const double* xarr, const double* yarr, int n);
Poly** spline(const double* xarr, const double* yarr, int n, double a, double b, int type);


#endif