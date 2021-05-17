#ifndef MKHA_H
#define MKHA_H

#include "prf.h"
#include "poly.h"
#include "tag.h"

/* Public parameters */
typedef struct {
    /* Prime number p in fmpz, represent the message space M = Z_p  */
    fmpz_t pf; 
    /* Represent the identity space ID = [n], user number */
    uint64_t n;
    /* Generator of group g1 */
    element_t g;
    /* Generator of group gt */
    element_t gt;
    /* Info for finite field */
    fq_ctx_t ctx;
    /* pairing */
    pairing_t pairing;
} PublicPara;

/* Verification Key */
typedef struct {

    /* The generated keys */
    Key* K;
    /* The identity of the user who runs the KeyGen algorithm  */
    uint64_t id; 
    /* alpha */
    fq_t alpha;
} VerKey;

typedef struct {
    uint64_t id;
    uint64_t tau;
} Label;

/* Initialize the groups */
int mkha_init();

/* Clean the groups */
int mkha_close();


/** The setup algorithm
 */
void set_up(PublicPara* pp, uint64_t n);


/** The key generation function
 * 
 *  @param[out] K       Key <- KG()
 * 
 */
void key_gen(PublicPara* pp, uint64_t id, VerKey* vk);

/** The authentication algorithm
 * 
 *  @param[in] sk       Secret Key
 *  @param[in] l        Label
 *  @param[in] m        Message
 *  @param[out] sig     Tag
 * 
 */
void auth(VerKey* sk, uint64_t delta, Label* l, fq_t m, Tag* sig, PublicPara* pp);

/**
 *  The evaluation algorithm
 * 
 *  @param[in] f        quadratic function
 *  @param[in] sig      list of tags \sig_i with i \in [t]
 */
void eval(Poly* f, Tag* sig, Tag* sig_out, uint64_t* id_set, PublicPara* pp);

/**
 *  The verification algorithm
 *  outputs 1 if accepted, otherwise 0.
 * 
 */
int ver(Poly* f, Label* l);

/* \cross: Multiplication of two inputs, say x_i x_j */
void GTE_cross(Tag* res, Tag* sig1, Tag* sig2, PublicPara* pp);

/*  Multiplication of an input with a constant, say c \cdot x_i */
void GTE_dot(Tag* res, fq_t c, Tag* sig, PublicPara* pp);

/* The addition of two results from \cross and \dot */
void GTE_add(Tag* res, Tag* sig1, Tag* sig2, PublicPara* pp);

/**
 *  The verification algorithm. 
 *  
 *  @param[in] f        f for the labeled program P
 *  @param[in] l        list of t labels, l_1, \ldots, l_t 
 *  @param[in] Dleta    dataset identifier
 *  @param[in] vk       a set of verification keys
 *  @param[in] m        message
 *  
 */
// int Ver(Poly* f, Label* l, uint64 Delta, VerKey* vk, g1_t m, Tag* sigma, PublicPara* pp);


/* Takes in Key and f outputs omega_f */
void cf_eval_off(Key* K, Label* l, Poly* f, uint64_t* id, Poly* omega_f, PublicPara* pp);

/** Outputs W. 
 *  Note that here, |K| = n, and K_i is the key corresponding to the 
 *  id n.
 **/
void cf_eval_on(Key* K, uint64_t Delta, uint64_t* id_array, Poly* omega_f, element_t W, PublicPara* pp);

#endif 