#ifndef _POLY_H_
#define _POLY_H_

#include <stdlib.h>

struct Poly{
    double coeff;
    int index;
    struct Poly* next;
};

typedef struct Poly Poly;

Poly* create_poly(double coeff, int index);
void destroy_poly(Poly** node);


//just destroy the head node, and move to the next
void static inline destroy_poly_head(Poly** head)
{
    Poly* tmp;
    if (NULL == head || NULL == *head)
        return;
    tmp = (*head)->next;
    free(*head);
    *head = tmp;
}

void print_poly(const Poly* poly);

Poly* poly_copy(const Poly* node);
void poly_add_inp(Poly** lhs, Poly* rhs);
void poly_sub_inp(Poly** lhs, Poly* rhs);
void poly_mul_inp(Poly** lhs, const Poly* rhs);
void poly_mul_cons_inp(Poly** lhs, double x);
void poly_div_inp(Poly** plhs, const Poly* rhs);
void poly_diff_inp(Poly** poly);
double poly_value(const Poly* poly, double x);
void poly_pow_inp(Poly** poly, int n);

#endif
