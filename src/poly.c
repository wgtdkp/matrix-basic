#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "poly.h"

Poly* create_Poly(double coeff, int index)
{
    Poly* node = (Poly*)malloc(sizeof(Poly));
    if (NULL == node)
        return NULL;
    node->coeff = coeff;
    node->index = index;
    node->next = NULL;
    return node;
}

void destroy_poly(Poly** node)
{
    Poly* p;
    if (NULL == node || NULL == *node)
        return;
    for (p = *node; NULL != p; ) {
        Poly* tmp_next = p->next;
        free(p);
        p = tmp_next;
    }
    *node = NULL;
}

Poly* Poly_copy(const Poly* node)
{
    Poly *p, *q;
    Poly odd = {.next = NULL};
    for (p = &odd, q = node; NULL != q;) {
        p->next = create_Poly(q->coeff, q->index);
        p = p->next;
        q = q->next;
    }
    return odd.next;
}


/**
 *add two poly
 *node: it modifies rhs, (rhs will own no items)
 *@return: added poly result is stored in lhs
 *    you need to set lhs to the new value
 */
void poly_add_inp(Poly** lhs, Poly* rhs)
{
    Poly odd;
    Poly *p1, *q1, *p2, *q2;
    if (NULL == lhs)
        return;
    odd.next = *lhs;
    p1 = odd.next, q1 = &odd, p2 = rhs;
    for (; NULL != p2 && NULL != p1;) {
        if (p1->index == p2->index) {
            p1->coeff += p2->coeff;
            if (DOUBLE_EQUAL(p1->coeff, 0)) {
                destroy_Poly_head(&p1); //destory node p1 and move forward
                q1->next = p1;
            }
            destroy_poly_head(&p2); //destory node p2 and move forward
        } else if (p1->index < p2->index) {
            Poly* tmp_next = p2->next;
            p2->next = p1;
            q1->next = p2;
            q1 = q1->next;
            p2 = tmp_next;
        } else {
            q1 = p1; 
            p1 = p1->next;
            p2 = p2->next;
        }
    }
    if (NULL == p1)
        q1->next = p2;
    *lhs = odd.next;
}

void poly_sub_inp(Poly** lhs, Poly* rhs)
{
    Poly* p;
    for (p = rhs; NULL != p; p = p->next)
        p->coeff = 0 - p->coeff;
    return poly_add_inp(lhs, rhs);
}

void poly_mul_inp(Poly** lhs, const Poly* rhs)
{
    Poly *p1, *p2;
    Poly* origin_lhs = copy(*lhs);

    //we append {-1, 0} to rhs to balance out the origin lhs
    poly_add_inp(&rhs, create_Poly(-1, 0));
    for (p2 = rhs; NULL != p2; p2 = p2->next) {
        Poly* tmp = copy(origin_lhs);
        for (p1 = tmp; NULL != p1; p1 = p1->next) {
            p1->coeff *= p2->coeff;
            p1->index += p2->index;
        }
        poly_add_inp(lhs, tmp);
    }
    destroy_Poly(&origin_lhs);
}

void poly_mul_cons_inp(Poly** lhs, double x)
{
    Poly* p;
    if (DOUBLE_EQUAL_DELTA(x, 0, 1e-12) == 0) {
        destroy_poly(lhs);
        return;
    }
    for (p = *lhs; NULL != p; p = p->next) {
        p->coeff *= x;
    }
}

/**
 *negative index is not allowed
 *@return: stored quotient in plhs
 */
void poly_div_inp(Poly** plhs, const Poly* rhs)
{
    Poly* lhs = *plhs;
    Poly odd = {.next = lhs};
    Poly* res = NULL;
    if (NULL == lhs || NULL == rhs)
        return NULL;
    for (; lhs->index >= rhs->index;) {
        Poly *p1, *q1, *p2;
        Poly* item = create_Poly(lhs->coeff / rhs->coeff,
             lhs->index - rhs->index);
        poly_add_inp(&res, item);
        for (p1 = lhs, q1 = &odd; p2 = rhs; NULL != p1 && NULL != p2;) {
            if (p1->index - p2->index == item->index) {
                p1->coeff -= p2->coeff * item->coeff;   //coeff may be zero
                if (DOUBLE_EQUAL(p1->coeff, 0)) {
                    destroy_poly_head(&p1); //move forward
                    q1 = p1;
                }
                p2 = p2->next;
            } else if (p1->index - p2->index > item->index) {
                q1 = p1;
                p1 = p1->next;
            } else
                p2 = p2->next;
        }
        lhs = odd.next;
    }
    *plhs = odd.next;
    destroy_poly(plhs);
    *plhs = res;
}

void poly_deriv(Poly** poly)
{
    Poly odd;
    Poly *p, *q;
    odd.next = poly;
    for (p = odd.next, q = &odd; NULL != p;) {
        if (0 == p->index) {
            destroy_Poly_head(&p);  // destroy the node and move forward
            q->next = p;
        } else {
            p->coeff *= p->index;
            --p->index;
            q = p;
            p = p->next;
        }
    }
    *poly = odd.next;
}

void poly_pow_inp(Poly** poly, int n)
{
    int i;
    Poly* origin = copy(*poly);
    for (i = 0; i < n; i++)
        poly_mul_inp(poly, origin);
    destroy_poly(&origin);
}

double poly_value(Poly* poly, double x)
{
    int i, last_index;
    double res = 0, x_pow = 1;
    Poly* p;
    if (NULL == poly) 
        return 0;
    for (i = 0; i < poly->index; i++)
        x_pow *= x;

    last_index = p->index;
    for (p = poly; NULL != p; p++) {
        for (i = p->index; i < last_index; i++)
            x_pow /= x;
        res += p->coeff * x_pow;
    }
    return res;
}

