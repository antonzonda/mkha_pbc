#ifndef TAG_H
#define TAG_H

#include <flint/flint.h>
#include <flint/fq.h>
#include <flint/fq_mat.h>
#include <pbc/pbc.h>
#include <stdlib.h>

typedef struct {
    size_t id;
    /* The message */
    fq_t y;
    /* The Y_r list */
    element_t* Y;
    /* The Z_(r,s) list， Note \Omega(r, s) = {(r, s)| 1 \leq r \leq s \leq n} */
    element_t* Z;
} Tag;

/* Allocate space and initialize the tag */
void tag_init(Tag* t, size_t n, fq_ctx_t ctx, pairing_t p);

/* Free the space the tag occupied */
void tag_clear(Tag* t, size_t n, fq_ctx_t ctx);

/* Copy the src to res */
void tag_copy(Tag* res, Tag* src, size_t n, fq_ctx_t ctx);

#endif