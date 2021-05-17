#ifndef POLY_H
#define POLY_H

#include <flint/flint.h>
#include <flint/fq.h>
#include <flint/fq_mat.h>

/**
 * The structure for the polynomial over \mathbb{Z}_p, namely
 * f(x_1, \cdots, x_n) = \sum_{i,j=1}^t \xi_{ij} x_i x_j + \sum_{i=1}^t \eta_i x_i
 * 
 */
typedef struct {
    /* The length of the polynomial */
    size_t t; 

    fq_t* xi; 
    fq_t* eta;

} Poly;

/**
 *  The algorithm randomly generates a polynomial
 *  over Z_p.
 * 
 */
void poly_rand_init(Poly* f, size_t t, fq_ctx_t ctx);

/**
 *  The algorithm clears the polynomial.
 */
void poly_clear(Poly* f, fq_ctx_t ctx);

/* evaluate the polynomial and output in res */
void poly_eval(Poly* f, fq_t* x, fq_t res, fq_ctx_t ctx);

#endif