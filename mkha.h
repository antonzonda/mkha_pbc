#ifndef MKHA_H
#define MKHA_H

#include "prf.h"
#include "poly.h"
#include "tag.h"
#include <pbc/pbc.h>

/* Public parameters */
typedef struct {
    /* Prime number p, represent the message space M = Z_p  */
    mpz_t p;
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


/** The setup algorithm
 */
void set_up(PublicPara* pp, uint64_t n, uint32_t lambda);


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
void eval(Poly* f, Tag* sig, Tag* sig_out, uint64_t* id_set_eval, PublicPara* pp);


/* \cross: Multiplication of two inputs, say x_i x_j */
void GTE_cross(Tag* res, Tag* sig1, Tag* sig2, PublicPara* pp);

/*  Multiplication of an input with a constant, say c \cdot x_i */
void GTE_dot(Tag* res, fq_t c, Tag* sig, PublicPara* pp);

/* The addition of two results from \cross and \dot */
void GTE_add(Tag* res, Tag* sig1, Tag* sig2, PublicPara* pp);

/**
 *  The verification algorithm. 
 *  outputs 1 if accepted, otherwise 0.
 * 
 *
 *  
 *  @param[in] omega    for the simplified labeled program P
 *  @param[in] Dleta    dataset identifier
 *  @param[in] vk       a set of verification keys
 *  @param[in] m        message
 *
 **/  
int eff_ver(Poly* omega, Label* l, uint64_t Delta, VerKey* vk, fq_t m, uint64_t* id_set,
        Tag* sigma, size_t len_t, PublicPara* pp);

/* Takes in Key and f outputs omega_f
 * The function works as VerPrep, 
 * Note that it outputs a polynomial omega_f with length 2n
 */
void cf_eval_off(VerKey* K, Label* l, Poly* f, size_t** id_set, size_t* id_t_list, Poly* omega_f, PublicPara* pp);

/** Outputs W. 
 *  Note that here, |K| = n, and K_i is the key corresponding to the 
 *  id n.
 **/
void cf_eval_on(VerKey* K, uint64_t Delta, Poly* omega_f, element_t W, PublicPara* pp);

/* The naive implementation of GPE
 * Would take longer time,
 * For DEBUG purpose
 **/
void GPE_naive(VerKey* vk, Poly* f, uint64_t Delta, Label* l, size_t* id_t_list, element_t W, PublicPara* pp);

/* Auxiliary function,
 * obtains g^y \prod_{r=1}^n e(Y_r, g)^\alpha_r \prod_{(r, s)\in\Omega_n} Z_{r,s}^{\alpha_r \alpha_s}
 **/
void get_W(Tag* sigma, VerKey* vk, element_t W, PublicPara* pp);

/* Clean the memory of pp and key */
void mkha_clear(PublicPara* pp, VerKey* vk);

#endif 