#include "integration.h"
#include <memory.h>
#include <stdlib.h>
#include <math.h>


double romberg(func f, double a, double b, double delta)
{
    size_t i, j;
    unsigned long long n;
    double h = b - a, t;
    double T[100] = {0};
    for (t = 0, i = 0, n = 1; i < 64; i++) {
        unsigned long long coeff;
        for (j = 0; j < n; j++)
            t += h * (*f)(a + (j + 0.5) * h, NULL);
        t *= 0.5;
        for (coeff = 1, j = 1; j <= i; j++) {
            double tmp_t = t;
            coeff <<= 2;
            t = (t * coeff - T[j-1]) / (coeff - 1);
            if (fabs(t - T[j-1]) < delta)
                goto END;
            T[j-1] = tmp_t;
        }
        T[i] = t;
        
        n <<= 1;
        h  = (b - a) / n;
    }
END:
    return t;
}

double simpson(func f, double a, double b, size_t n)
{
    size_t k = 0;
    double S;
    double h = (b - a) / n;
    for (S = 0, k = 0; k < n; k++) {
        S += 4 * (*f)(a + (k + 0.5) * h, NULL);
        S += 2 * (*f)(a + (k) * h, NULL);
    }
    S += (*f)(b, NULL) - (*f)(a, NULL);
    return S * h / 6;
}

static double __adaptive(func f, double a, double b, double ss, double delta)
{
    double left, right;
    //printf("ss: %lf, a: %lf, b: %lf\n", ss, a, b);
    left = simpson(f, a, (a + b) / 2, 1);
    right = simpson(f, (a + b) / 2, b, 1);
    //printf("left: %lf, right: %lf\n", left, right);
    //printf("diff: %lf, delta: %lf\n", fabs(ss - (left + right)), 15 * delta);
    //printf("\n");
    if (fabs(ss - (left + right)) < 15 * delta)
        return (16 * (left + right) - ss) / 15;

    left = __adaptive(f, a, (a + b) / 2, left, delta / 2);
    right = __adaptive(f, (a + b) / 2, b, right, delta / 2);

END:
    //printf("*************: %lf\n", (16 * (left + right) - ss) / 15);
    return left + right;
}

double integration(func f, double a, double b, double delta)
{
    double ss = simpson(f, a, b, 1);
    return __adaptive(f, a, b, ss, delta);
}

double diff(func f, double x, double delta)
{
    size_t i, j;
    double h, g;
    double G[100] = {0};
    for (h = 1, i = 0; i < 100; i++) {
        unsigned long long coeff;
        g = ((*f)(x + h, NULL) - (*f)(x - h, NULL)) / (2 * h);
        for (coeff = 1, j = 1; j <= i; j++) {
            double tmp_g = g;
            coeff <<= 2;
            g = (g * coeff - G[j-1]) / (coeff - 1);
            if (fabs(g - G[j-1]) < delta)
                goto END;
            G[j-1] = tmp_g;
        }
        G[i] = g;

        h /= 2;
    }
END:
    //printf("i: %d\n", i);
    return g;
}
