#include "poly.h"


void poly_rand_init(Poly* f, size_t t, fq_ctx_t ctx) {
    size_t i, j;

    f->t = t;
    f->xi = (fq_t*) malloc(sizeof(fq_t) * t * t);
    f->eta = (fq_t*) malloc(sizeof(fq_t) * t);

    flint_rand_t state;
    flint_randinit(state);

    for (i = 0; i < t; i++) {
        for (j = 0; j < t; j++) {
            fq_init(f->xi[i * t + j], ctx);
            fq_randtest(f->xi[i * t + j], state, ctx);
        }
    }

    for (i = 0; i < t; i++) {
        fq_init(f->eta[i], ctx);
        fq_randtest(f->eta[i], state, ctx);
    }

    flint_randclear(state);
}

void poly_clear(Poly* f, fq_ctx_t ctx) {
    size_t i, j;

    for (i = 0; i < f->t; i++) {
        for (j = 0; j < f->t; j++) {
            fq_clear(f->xi[i * f->t + j], ctx);
        }
    }

    for (i = 0; i < f->t; i++) {
        fq_clear(f->eta[i], ctx);
    }

    free(f->eta);
    free(f->xi);
}


void poly_eval(Poly* f, fq_t* x, fq_t res, fq_ctx_t ctx) {
    size_t i, j;

    // Temporary variables
    fq_t t1;
    fq_t t2;
    fq_init(t1, ctx);
    fq_init(t2, ctx);

    // Set to zero
    fq_zero(res, ctx);

    for (i = 0; i < f->t; i++) {
        fq_mul(t1, f->eta[i], x[i], ctx);
        fq_add(res, t1, res, ctx);
    }

    for (i = 0; i < f->t; i++) {
        for (j = 0; j < f->t; j++) {
            fq_mul(t1, x[i], x[j], ctx);
            fq_mul(t2, t1, f->xi[j + i * f->t], ctx);
            fq_add(res, t2, res, ctx);
        }
    }
    
    fq_clear(t1, ctx);
    fq_clear(t2, ctx);
}

void poly_print(Poly* f, fq_ctx_t ctx) {
    size_t i, j;
    printf("Eta: ");
    for (i = 0; i < f->t; i++) {
        printf("%ld: ", i);
        fq_print_pretty(f->eta[i], ctx);
        printf(", ");
    }
        
    printf("Xi: ");
    for (i = 0; i < f->t; i++) {
        for (j = 0; j < f->t; j++) {
            printf("(%ld %ld): ", i, j);
            fq_print_pretty(f->xi[i * f->t + j], ctx);
            printf(", ");
        }
    }
}
