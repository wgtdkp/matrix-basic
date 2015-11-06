#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "poly.h"
#include "matrix.h"

Poly* create_poly(double coeff, int index)
{
    Poly* poly;
    if (DOUBLE_EQUAL_DELTA(coeff, 0, 1e-12) || index < 0)
        return NULL;

    poly = (Poly*)malloc(sizeof(Poly));
    if (NULL == poly)
        return NULL;
    poly->coeff = coeff;
    poly->index = index;
    poly->next = NULL;
    return poly;
}

void destroy_poly(Poly** poly)
{
    Poly* p;
    if (NULL == poly || NULL == *poly)
        return;
    for (p = *poly; NULL != p; ) {
        Poly* tmp_next = p->next;
        free(p);
        p = tmp_next;
    }
    *poly = NULL;
}

void print_poly(const Poly* poly)
{
    const Poly* p;
    printf("[");
    for (p = poly; NULL != p; p = p->next) {
        printf("(%.6lf, %d), ", p->coeff, p->index);
    }
    printf("]\n");
}

Poly* poly_copy(const Poly* poly)
{
    Poly* p;
    const Poly* q;
    Poly odd = {.next = NULL};
    for (p = &odd, q = poly; NULL != q;) {
        p->next = create_poly(q->coeff, q->index);
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
    Poly *p1, *q1, *p2;
    if (NULL == lhs)
        return;
    if (*lhs == rhs)
        return poly_mul_cons_inp(lhs, 2);
    odd.next = *lhs;
    p1 = odd.next, q1 = &odd, p2 = rhs;
    for (; NULL != p2 && NULL != p1;) {
        if (p1->index == p2->index) {
            p1->coeff += p2->coeff;
            if (DOUBLE_EQUAL(p1->coeff, 0)) {
                destroy_poly_head(&p1); //destory node p1 and move forward
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
        }
    }
    if (NULL == p1)
        q1->next = p2;

    *lhs = odd.next;
}

void poly_sub_inp(Poly** lhs, Poly* rhs)
{
    Poly* p;
    if (NULL == lhs) return;
    if (*lhs == rhs)
        return destroy_poly(lhs);
    for (p = rhs; NULL != p; p = p->next)
        p->coeff = 0 - p->coeff;
    return poly_add_inp(lhs, rhs);
}

void poly_mul_inp(Poly** lhs, const Poly* rhs)
{
    Poly* p1;
    Poly* res = NULL;
    const Poly* p2;
    //think about the situation that: lhs and rhs are the same poly ?
    Poly* origin_lhs = poly_copy(*lhs); 
    
    //we append {-1, 0} to rhs to balance out the origin lhs
    for (p2 = rhs; NULL != p2; p2 = p2->next) {
        Poly* tmp = poly_copy(origin_lhs);
        for (p1 = tmp; NULL != p1; p1 = p1->next) {
            p1->coeff *= p2->coeff;
            p1->index += p2->index;
        }
        poly_add_inp(&res, tmp);
    }
    destroy_poly(lhs);
    *lhs = res;
    destroy_poly(&origin_lhs);
}

void poly_mul_cons_inp(Poly** lhs, double x)
{
    Poly* p;
    if (DOUBLE_EQUAL_DELTA(x, 0, 1e-12)) {
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
        return;
    if (*plhs == rhs) { //lhs and rhs are the same poly
        destroy_poly(plhs);
        *plhs = create_poly(1, 0);
        return;
    }
    for (; NULL != lhs && lhs->index >= rhs->index;) {
        Poly *p1, *q1;
        const Poly* p2;
        Poly* item = create_poly(lhs->coeff / rhs->coeff,
             lhs->index - rhs->index);
        poly_add_inp(&res, poly_copy(item));
        //printf("res: \n");
        //print_poly(res);
        for (p1 = lhs, q1 = &odd, p2 = rhs; NULL != p1 && NULL != p2;) {
            //print_poly(p2);
            if (p1->index - p2->index == item->index) {
                p1->coeff -= p2->coeff * item->coeff;   //coeff may be zero
                if (DOUBLE_EQUAL(p1->coeff, 0)) {
                    destroy_poly_head(&p1); //move forward
                    q1->next = p1;
                }
                p2 = p2->next;
            } else if (p1->index - p2->index > item->index) {
                q1 = p1;
                p1 = p1->next;
            } else
                p2 = p2->next;
        }
        lhs = odd.next;
        //printf("lhs: \n");
        //print_poly(lhs);
    }
    *plhs = odd.next;
    destroy_poly(plhs);
    *plhs = res;
}

void poly_deriv_inp(Poly** poly)
{
    Poly odd;
    Poly *p, *q;
    odd.next = *poly;
    for (p = odd.next, q = &odd; NULL != p;) {
        if (0 == p->index) {
            destroy_poly_head(&p);  // destroy the node and move forward
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
    Poly* origin;
    if (0 == n) {
        *poly = create_poly(1, 0);
        return;
    }
    origin = poly_copy(*poly);
    for (i = 0; i < n - 1; i++)
        poly_mul_inp(poly, origin);
    destroy_poly(&origin);
}

double poly_value(const Poly* poly, double x)
{
    double res = 0;
    const Poly* p;
    if (NULL == poly) 
        return 0;
    for (p = poly; NULL != p; p = p->next)
        res += p->coeff * pow(x, p->index);
    return res;
}
