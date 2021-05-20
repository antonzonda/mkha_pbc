#include "mkha.h"

static void fq2mpz(mpz_t out, fq_t in, fq_ctx_t ctx) {
    fmpz_t temp;
    fmpz_init(temp);
    fmpz_set_str(temp, fq_get_str_pretty(in, ctx), 10);
    fmpz_get_mpz(out , temp);
    fmpz_clear(temp);
}

static void mpz2fq(fq_t out, mpz_t in, fq_ctx_t ctx) {
    fmpz_t temp;
    fmpz_init(temp);
    fmpz_set_mpz(temp, in);
    fq_set_fmpz(out, temp, ctx);
}

void set_up(PublicPara* pp, uint64_t n, uint32_t lambda) {
    /* Choose a description of pp */
    pp->n = n; // Set n

    // Generate the prime number
    mpz_init(pp->p);
    mp_limb_t p;
    flint_rand_t state;
    flint_randinit(state);
    p = n_randprime(state, lambda, 1);
    mpz_set_ui(pp->p, p);

    // Set pairing
    pbc_param_t par;
    pbc_param_init_a1_gen(par, pp->p);
    // init pairing
    pairing_init_pbc_param(pp->pairing, par);
    
    // init g, gt
    element_init_G1(pp->g, pp->pairing);
    element_random(pp->g);
    element_init_GT(pp->gt, pp->pairing);
    element_random(pp->gt);

    // init ctx
    fmpz_t p1;
    fmpz_init(p1);
    fmpz_set_mpz(p1, pp->p);
    fq_ctx_init(pp->ctx, p1, 1, "ctx");
    
    // free memory
    fmpz_clear(p1);
    pbc_param_clear(par);
    flint_randclear(state);
}

void key_gen(PublicPara* pp, uint64_t id, VerKey* vk) {
    vk->K = (Key*) malloc(sizeof(Key));
    KG(vk->K); // Generate the keys
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
    element_set0(Y);
    element_t R;
    element_init_G1(R, pp->pairing);
    F(sk->K, (uint8_t*) &delta, (uint8_t*) l, pp->g, R, pp->p);

    mpz_t m_mpz; /* Message in bn */
    mpz_init(m_mpz);
    fq2mpz(m_mpz, m, pp->ctx);

    mpz_t neg_m_mpz; /* -m */
    mpz_init(neg_m_mpz);
    mpz_neg(neg_m_mpz, m_mpz);

    element_t x;
    element_init_G1(x, pp->pairing);
    element_set0(x);
    element_mul_mpz(x, pp->g, neg_m_mpz);
    element_add(Y, x, Y);
    
    // Set alpha^{-1}
    fq_t alpha_inv;
    fq_init(alpha_inv, pp->ctx);
    fq_inv(alpha_inv, sk->alpha, pp->ctx);
    mpz_t alpha_inv_m;
    mpz_init(alpha_inv_m);
    fq2mpz(alpha_inv_m, alpha_inv, pp->ctx);

    // Allocate the list
    tag_init(sig, pp->n, pp->ctx, pp->pairing);
    // Set Y[id]
    element_mul_mpz(sig->Y[id], Y, alpha_inv_m);

    // Free memory
    fq_clear(alpha_inv, pp->ctx);
    element_clear(x);
    element_clear(Y);
    element_clear(R);
    mpz_clear(m_mpz); 
    mpz_clear(neg_m_mpz);
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

    tag_init(sig_out, pp->n, pp->ctx, pp->pairing);
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
    element_set0(t1);
    element_init_G1(t2, pp->pairing);
    element_set0(t2);

    element_t t3, t4;
    element_init_GT(t3, pp->pairing);
    element_set1(t3);
    element_init_GT(t4, pp->pairing);
    element_set1(t4);

    // Calculate Y_r
    for (r = 0; r < n; r++) {
        element_mul(t1, sig1->Y[r], sig1->Y[r]);
        element_mul(t2, sig2->Y[r], sig2->Y[r]);
        element_add(res->Y[r], t1, t2);
    }

    // Calculate Z_(r,s)
    for (r = 1; r <= n; r++) {
        i = ((r-1)*n + (r-1)*(r-2)/2);
        s = r;
        element_pairing(res->Z[i], sig1->Y[r-1], sig2->Y[s-1]);
        i++;
        for (s = r + 1; s <= n; s++) {
            element_pairing(t3, sig1->Y[r-1], sig2->Y[s-1]);
            element_pairing(t4, sig1->Y[s-1], sig2->Y[r-1]);
            element_mul(res->Z[i], t3, t4);
            i++;
        }
    }

    // free the things allocated
    element_clear(t1);
    element_clear(t2);
    element_clear(t3);
    element_clear(t4);
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
        element_pow_mpz(res->Y[r], sig->Y[r], c_mpz);
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
        element_add(res->Y[r], sig1->Y[r], sig2->Y[r]);
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

    element_t W1;
    element_init_GT(W1, pp->pairing);
    cf_eval_on(vk->K, Delta, id_set, &omega, W1, pp);

    // temporary variables
    element_t res, t2, t3;
    element_init_GT(res, pp->pairing); element_init_GT(t2, pp->pairing); element_init_GT(t3, pp->pairing);
    mpz_t b1, b2, b3;
    mpz_init(b1);
    fq_t f1;
    fq_init(f1, pp->ctx);

    // Calculate the right hand side
    element_t W2;
    element_init_GT(W2, pp->pairing);
    // convert y to mpz
    mpz_t y_mpz;
    mpz_init(y_mpz);
    fq2mpz(y_mpz, sigma->y, pp->ctx);
    element_pow_mpz(res, pp->gt, y_mpz);  // ??????
    for (r = 0; r < pp->n; r++) {
        element_pairing(t2, sigma->Y[r], pp->g);
        fq2mpz(b1, vk[r].alpha, pp->ctx);
        element_pow_mpz(t3, t2, b1);
        element_mul(res, t3, res);
    }

    for (r = 1; r <= pp->n; r++) {
        i = ((r-1)*pp->n + (r-1)*(r-2)/2);
        for (s = r; s <= pp->n; s++) {
            fq_mul(f1, vk[r-1].alpha, vk[s-1].alpha, pp->ctx);
            fq2mpz(b1, f1, pp->ctx);
            element_pow_mpz(t3, sigma->Z[i], b1);
            element_mul(res, t3, res);
            i++;
        }
    }

    r3 = element_cmp(W1, W2);

    // Free memory
    element_clear(res); element_clear(t2); element_clear(t3);
    mpz_clear(b1); fq_clear(f1, pp->ctx);
    element_clear(W1); element_clear(W2);
    
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
    // Temporary variable
    mpz_t u1, v1;
    mpz_init(u1); mpz_init(v1);

    for (i = 0; i < f->t; i++) {
        fq_init(u[i], pp->ctx);
        fq_init(v[i], pp->ctx);
        PRF_F(u1, v1, K[id[i]].k1, (uint8_t*) l, 16, pp->p);
        mpz2fq(u[i], u1, pp->ctx);
        mpz2fq(v[i], v1, pp->ctx);
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
    mpz_clear(u1); mpz_clear(v1);
    for (i = 0; i < f->t; i++) {
        fq_clear(u[i], pp->ctx);
        fq_clear(v[i], pp->ctx);
    }
}

void cf_eval_on(Key* K, uint64_t Delta, uint64_t* id_array, Poly* omega_f, element_t W, PublicPara* pp) {
    size_t i, j, id;
    size_t t = omega_f->t / 2;

    fq_t a[pp->n];
    fq_t b[pp->n];
    // Temporary variable
    mpz_t a1, b1;
    mpz_init(a1); mpz_init(b1);
    
    for (j = 0; j < pp->n; j++) {
        fq_init(a[j], pp->ctx);
        fq_init(b[j], pp->ctx);

        PRF_F(a1, b1, K[j].k2, (uint8_t *) &Delta, 8, pp->p);
        mpz2fq(a[i], a1, pp->ctx);
        mpz2fq(b[i], b1, pp->ctx);
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

    // free
    mpz_clear(a1); mpz_clear(b1);
    for (j = 0; j < pp->n; j++) {
        fq_clear(a[j], pp->ctx);
        fq_clear(b[j], pp->ctx);
    }
}

