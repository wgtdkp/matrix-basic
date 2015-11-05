#ifndef _POLY_H_
#define _POLY_H_

#include <stdlib.h>

struct Poly{
    double coeff;
    int index;
    struct Poly* next;
} Poly;

typedef struct Poly Poly;

Poly* create_poly(int n);
void destroy_poly(Poly** node);


//just destroy the head node, and move to the next
void static inline destroy_poly_head(Poly** head)
{
    Poly* tmp;
    if (NULL == head || NULL == *head)
        return;
    tmp = *head->next;
    free(*head);
    *head = tmp;
}

Poly* Poly_copy(const Poly* node);
void poly_add_inp(Poly** lhs, Poly* rhs);
void poly_sub_inp(Poly** lhs, Poly* rhs);
void poly_mul_inp(Poly** lhs, const Poly* rhs);
void poly_mul_cons_inp(Poly** lhs, double x);
void poly_div_inp(Poly** plhs, const Poly* rhs);
void poly_deriv(Poly** poly);
double poly_value(Poly* poly, double x);
void poly_pow_inp(Poly** poly, int n);

int static inline poly_max_index(const Poly* poly)
{
    if (0 == poly->n)
        return 0;
    return poly->head->index;
}






#endif