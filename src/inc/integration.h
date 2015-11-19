#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_
#include <stddef.h>

typedef double(*func)(double x, void* coeff);

double romberg(func f, double a, double b, double delta);
double simpson(func f, double a, double b, size_t n);
double integration(func f, double a, double b, double delta);
double diff(func f, double x, double delta);

#endif
