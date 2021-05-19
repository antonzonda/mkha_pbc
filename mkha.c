#include "mkha.h"

static void fq2mpz(mpz_t out, fq_t in, fq_ctx_t ctx) {
    fmpz_t temp;
    fmpz_init(temp);
    fmpz_set_str(temp, fq_get_str_pretty(in, ctx), 10);
    fmpz_set_mpz(out , temp);
    fmpz_clear(temp);
}

void set_up(PublicPara* pp, uint64_t n, uint32_t lambda) {
    /* Choose a description of pp */
    mpz_init(pp->p);
    pbc_mpz_randomb(pp->p, lambda);

    pbc_param_t par;
    pbc_param_init_a1_gen(par, pp->p);
    // init pairing
    pairing_init_pbc_param(pp->pairing, par);
    
    // init g, gt
    element_init_G1(pp->g, pp->pairing);
    element_random(pp->g);
    element_init_GT(pp->gt, pp->pairing);
    element_random(pp->gt);
    
    pbc_param_clear(par);
}

void key_gen(PublicPara* pp, uint64_t id, VerKey* vk) {
    Key* K;
    KG(K);
    vk->K = K;
    vk->id = id;
    // Set alpha <--Z_p*
    fq_init(vk->alpha, pp->ctx);
    flint_rand_t state;
    flint_randinit(state);
    fq_randtest(vk->alpha, state, pp->ctx);
    flint_randclear(state);
}

void auth(VerKey* sk, uint64_t delta, Label* l, fq_t m, Tag* sig, PublicPara* pp) {
    size_t r, s;
    size_t id = l->id;

    // Set the message y
    sig->id = id;
    fq_init(sig->y, pp->ctx);
    fq_set(sig->y, m, pp->ctx);
    
    // Set the Y_r
    element_t Y;
    element_init_G1(Y, pp->pairing);
    g1_set_infty(Y);

    F(sk->K, (uint8_t*) &delta, (uint8_t*) l, pp->g, Y);

    mpz_t m_mpz; /* Message in bn */
    mpz_init(m_mpz);
    fq2mpz(m_mpz, m, pp->ctx);

    bn_t neg_m_bn; /* -m */
    bn_new(neg_m_bn);
    bn_neg(neg_m_bn, m_bn);

    g1_t x;
    g1_new(x);
    g1_set_infty(x);
    g1_mul(x, x, neg_m_bn);
    g1_norm(x, x);
    g1_add(Y, x, Y);
    g1_norm(Y, Y);
    // Set alpha^{-1}
    fq_t alpha_inv;
    fq_init(alpha_inv, pp->ctx);
    fq_inv(alpha_inv, sk->alpha, pp->ctx);
    bn_t alpha_inv_b;
    bn_new(alpha_inv_b);
    fq2bn(alpha_inv_b, alpha_inv, pp->ctx);

    // Allocate the list
    tag_init(sig, pp->n, pp->ctx);
    g1_mul(sig->Y[id], Y, alpha_inv_b);

    fq_clear(alpha_inv, pp->ctx);
}

void eval(Poly* f, Tag* sig, Tag* sig_out, uint64_t* id_set, PublicPara* pp) {
    size_t i, j, indi;
    indi = 1;
    Tag t1, t2, res[2];
    // Need to initialize t1, t2
    tag_init(&t1, pp->n, pp->ctx, pp->pairing);
    tag_init(&t2, pp->n, pp->ctx, pp->pairing);
    tag_init(&res[0], pp->n, pp->ctx, pp->pairing);
    tag_init(&res[1], pp->n, pp->ctx, pp->pairing);

    // First calculate id_set, that is, let
    // id_set[id] := 1, if id \in some sig
    for (i = 0; i < f->t; i++) {
        id_set[sig[i].id] = 1;
    }

    // Evaluate the x_i x_j part
    for (i = 0; i < f->t; i++) {
        for (j = 0; j < f->t; j++) {
            GTE_cross(&t1, &sig[i], &sig[j], pp);
            if (i == 0 && j == 0) {
                GTE_dot(&res[indi], f->xi[j + f->t * i], &t1, pp);
            } else {
                GTE_dot(&t2, f->xi[j + f->t * i], &t1, pp);
                GTE_add(&res[indi], &res[!indi], &t2, pp);
            }
            indi = !indi;
        }
    }

    // Evaluate the x_i part
    for (i = 0; i < f->t; i++) {
        GTE_dot(&t1, f->eta[i], &sig[i], pp);
        GTE_add(&res[indi], &res[!indi], &t2, pp);
        indi = !indi;
    }

    tag_init(sig_out, pp->n, pp->ctx);
    tag_copy(sig_out, &res[!indi], pp->n, pp->ctx);

    tag_clear(&t1, pp->n, pp->ctx);
    tag_clear(&t2, pp->n, pp->ctx);
    tag_clear(&res[0], pp->n, pp->ctx);
    tag_clear(&res[1], pp->n, pp->ctx);
}


void GTE_cross(Tag* res, Tag* sig1, Tag* sig2, PublicPara* pp) {
    size_t r, s, i;
    size_t n = pp->n;

    fq_mul(res->y, sig1->y, sig2->y, pp->ctx);

    element_t t1, t2; // Some temporary variables
    element_init_G1(t1, pp->pairing);
    g1_set_infty(t1);
    g1_new(t2);
    g1_set_infty(t2);

    element_t t3, t4;
    gt_new(t3);
    gt_set_unity(t3);
    gt_new(t4);
    gt_set_unity(t4);

    // Calculate Y_r
    for (r = 0; r < n; r++) {
        g1_mul(t1, sig1->Y[r], sig1->yb);
        g1_mul(t2, sig2->Y[r], sig2->yb);
        g1_add(res->Y[r], t1, t2);
    }

    // Calculate Z_(r,s)
    for (r = 1; r <= n; r++) {
        i = ((r-1)*n + (r-1)*(r-2)/2);
        s = r;
        pc_map(res->Z[i], sig1->Y[r-1], sig2->Y[s-1]);
        i++;
        for (s = r + 1; s <= n; s++) {
            pc_map(t3, sig1->Y[r-1], sig2->Y[s-1]);
            pc_map(t4, sig1->Y[s-1], sig2->Y[r-1]);
            gt_mul(res->Z[i], t3, t4);
            i++;
        }
    }

    // free the things allocated
    g1_free(t1);
    g1_free(t2);
    gt_free(t3);
    gt_free(t4);
}

void GTE_dot(Tag* res, fq_t c, Tag* sig, PublicPara* pp) {
    size_t r, s, i;
    size_t n = pp->n;

    fq_mul(res->y, c, sig->y, pp->ctx); // Get message y
    mpz_t c_mpz; // c in mpz
    mpz_init(c_mpz);
    fq2mpz(c_mpz, c, pp->ctx);

    // Calculate Y_r
    for (r = 0; r < n; r++) {
        g1_mul(res->Y[r], sig->Y[r], cb);
        g1_norm(res->Y[r], res->Y[r]);
    }

    // Calculate Z_(r,s)
    for (r = 1; r <= n; r++) {
        i = ((r-1)*n + (r-1)*(r-2)/2);
        for (s = r; s <= n; s++) {
            element_pow_mpz(res->Z[i], sig->Z[i], c_mpz);
            i++;
        }
    }

    mpz_clear(c_mpz);
}

void GTE_add(Tag* res, Tag* sig1, Tag* sig2, PublicPara* pp) {
    size_t r, s, i;
    size_t n = pp->n;

    fq_add(res->y, sig1->y, sig2->y, pp->ctx); // Get message y

    for (r = 0; r < n; r++) {
        g1_add(res->Y[r], sig1->Y[r], sig2->Y[r]);
    }

    for (r = 1; r <= n; r++) {
        i = ((r-1)*n + (r-1)*(r-2)/2);
        for (s = r; s <= n; s++) {
            element_mul(res->Z[i], sig1->Z[i], sig2->Z[i]);
            i++;
        }
    }
}

int Ver(Poly* f, Label* l, uint64_t Delta, VerKey* vk, fq_t m, uint64_t* id_set, uint64_t* id_t_list,
        Tag* sigma, PublicPara* pp) {

    int r1 = 0, r2 = 1, r3 = 0;
    size_t i, r, s;
    r1 = fq_equal(m, sigma->y, pp->ctx);

    for (i = 0; i < f->t; i++) {
        if (id_set[l[i].id] == 0) {
            r2 = 0;
            break;
        }
    }

    Poly omega;
    cf_eval_off(vk->K, l, f, id_t_list, &omega, pp);

    gt_t W1;
    gt_new(W1);
    cf_eval_on(vk->K, Delta, id_set, &omega, W1, pp);

    // temporary variables
    gt_t res, t2, t3;
    gt_new(res); gt_new(t2); gt_new(t3);
    bn_t b1, b2, b3;
    bn_new(b1);
    fq_t f1;
    fq_init(f1, pp->ctx);

    // Calculate the right hand side
    gt_t W2;
    gt_new(W2);
    gt_exp(res, pp->gt, sigma->yb);
    for (r = 0; r < pp->n; r++) {
        pc_map(t2, sigma->Y[r], pp->g);
        fq2bn(b1, vk[r].alpha, pp->ctx);
        gt_exp(t3, t2, b1);
        gt_mul(res, t3, res);
    }

    for (r = 1; r <= pp->n; r++) {
        i = ((r-1)*pp->n + (r-1)*(r-2)/2);
        for (s = r; s <= pp->n; s++) {
            fq_mul(f1, vk[r-1].alpha, vk[s-1].alpha, pp->ctx);
            fq2bn(b1, f1, pp->ctx);
            gt_exp(t3, sigma->Z[i], b1);
            gt_mul(res, t3, res);
            i++;
        }
    }

    r3 = gt_cmp(W1, W2);

    gt_free(res); gt_free(t2); gt_free(t3);
    bn_free(b1); fq_clear(f1, pp->ctx);
    gt_free(W1); gt_free(W2);
    return r1 * r2 * r3;
}

void cf_eval_off(Key* K, Label* l, Poly* f, uint64_t* id, Poly* omega_f, PublicPara* pp) {
    size_t i, j;
    omega_f->t = 2 * f->t;
    omega_f->xi = (fq_t *) malloc(sizeof(fq_t) * omega_f->t * omega_f->t);
    omega_f->eta = (fq_t *) malloc(sizeof(fq_t) * omega_f->t);
    
    // Use the pseudo random function
    fq_t u[f->t];
    fq_t v[f->t];
    for (i = 0; i < f->t; i++) {
        fq_init(u[i], pp->ctx);
        fq_init(v[i], pp->ctx);
        PRF_F(u[i], v[i], K[id[i]].k1, (uint8_t*) l, 16, pp->ctx);
    }

    fq_t temp;
    fq_init(temp, pp->ctx);

    // For omega_f eta
    for (i = 0; i < f->t; i++) {
        fq_init(omega_f->eta[i], pp->ctx);
        fq_mul(omega_f->eta[i], f->eta[i], u[i], pp->ctx);

        fq_init(omega_f->eta[i + f->t], pp->ctx);
        fq_mul(omega_f->eta[i + f->t], f->eta[i + f->t], v[i], pp->ctx);
    }

    // For omega_f xi
    for (i = 0; i < f->t; i++) {
        for (j = 0; j < f->t; j++) {
            fq_init(omega_f->xi[j + omega_f->t * i], pp->ctx);
            fq_mul(temp, u[i], u[j], pp->ctx);
            fq_mul(omega_f->xi[j + omega_f->t * i], temp, f->xi[i * f->t + j], pp->ctx);

            fq_init(omega_f->xi[j + f->t + omega_f->t * i], pp->ctx);
            fq_mul(temp, u[i], v[j], pp->ctx);
            fq_mul(omega_f->xi[j + f->t + omega_f->t * i], temp, f->xi[i * f->t + j], pp->ctx);

            fq_init(omega_f->xi[(i + f->t) * omega_f->t + j], pp->ctx);
            fq_mul(temp, v[i], u[j], pp->ctx);
            fq_mul(omega_f->xi[(i + f->t) * omega_f->t + j], temp, f->xi[i * f->t + j], pp->ctx);

            fq_init(omega_f->xi[i + f->t + omega_f->t * (j + f->t)], pp->ctx);
            fq_mul(temp, v[i], v[j], pp->ctx);
            fq_mul(omega_f->xi[i + f->t + omega_f->t * (j + f->t)], temp, f->xi[i * f->t + j], pp->ctx);
        }
    }

    // free
    fq_clear(temp, pp->ctx);
}

void cf_eval_on(Key* K, uint64_t Delta, uint64_t* id_array, Poly* omega_f, gt_t W, PublicPara* pp) {
    size_t i, j, id;
    size_t t = omega_f->t / 2;

    fq_t a[pp->n];
    fq_t b[pp->n];
    
    for (j = 0; j < pp->n; j++) {
        fq_init(a[j], pp->ctx);
        fq_init(b[j], pp->ctx);

        PRF_F(a[j], b[j], K[j].k2, (uint8_t *) &Delta, 8, pp->ctx);
    }

    // Set the input to the \omega_f
    fq_t x_in[omega_f->t];
    for (i = 0; i < t; i++) {
        id = id_array[i];
        fq_init(x_in[i], pp->ctx);
        fq_set(x_in[i], a[id], pp->ctx);

        fq_init(x_in[i + t], pp->ctx);
        fq_set(x_in[i + t], b[id], pp->ctx);
    }
    
    // Evaluation
    fq_t w;
    poly_eval(omega_f, x_in, w, pp->ctx);
}

