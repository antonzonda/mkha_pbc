#include "tag.h"

void tag_init(Tag* t, size_t n, fq_ctx_t ctx, pairing_t p) {
    size_t r, s, i;

    fq_init(t->y, ctx);

    t->Y = (element_t*) malloc(sizeof(element_t) * n);
    t->Z = (element_t*) malloc(sizeof(element_t) * (n*(n+1)/2));

    for (r = 0; r < n; r++) {
        element_init_G1(t->Y[r], p);
        element_set0(t->Y[r]);
    }

    i = 0;
    for (r = 1; r <= n; r++) {
        // i = ((r-1)*n + (r-1)*(r-2)/2);
        for (s = r; s <= n; s++) {
            element_init_GT(t->Z[i], p);
            element_set1(t->Z[i]);
            i++;
        }
    }
} 

void tag_clear(Tag* t, size_t n, fq_ctx_t ctx) {
    size_t r, s, i;

    fq_clear(t->y, ctx);

    for (r = 0; r < n; r++) {
        element_clear(t->Y[r]);
    }

    i = 0;
    for (r = 1; r <= n; r++) {
        // i = ((r-1)*n + (r-1)*(r-2)/2);
        for (s = r; s <= n; s++) {
            element_clear(t->Z[i]);
            i++;
        }
    }

    free(t->Y);
    free(t->Z);
}

void tag_copy(Tag* res, Tag* src, size_t n, fq_ctx_t ctx) {
    size_t r, s, i;

    fq_set(res->y, src->y, ctx);

    for (r = 0; r < n; r++) {
        element_set(res->Y[r], src->Y[r]);
    }

    i = 0;
    for (r = 1; r <= n; r++) {
        // i = ((r-1)*n + (r-1)*(r-2)/2);
        for (s = r; s <= n; s++) {
            element_set(res->Z[i], src->Z[i]);
            i++;
        }
    }
}
    