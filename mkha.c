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
    fmpz_t p;
    fmpz_init(p);
    flint_rand_t state;
    flint_randinit(state);
    fmpz_randprime(p, state, lambda, 1);
    fmpz_get_mpz(pp->p, p);
    fmpz_clear(p);

    // Set pairing
    pbc_param_t par;
    pbc_param_init_a1_gen(par, pp->p);
    // init pairing
    pairing_init_pbc_param(pp->pairing, par);
    
    // init g, gt
    element_init_G1(pp->g, pp->pairing);
    element_random(pp->g);
    element_init_GT(pp->gt, pp->pairing);
    element_pairing(pp->gt, pp->g, pp->g);

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
    size_t id = l->id;

    // Allocate the list
    tag_init(sig, pp->n, pp->ctx, pp->pairing);

    // Set the message y
    sig->id = id;
    fq_init(sig->y, pp->ctx);
    fq_set(sig->y, m, pp->ctx);
    
    // Set the Y_r
    element_t Y;
    element_init_G1(Y, pp->pairing);
    element_set1(Y);
    element_t R;
    element_init_G1(R, pp->pairing);
    F(sk->K, (uint8_t*) &delta, (uint8_t*) l, pp->g, R, pp->p);

    mpz_t m_mpz; /* Message in bn */
    mpz_init(m_mpz);
    fq2mpz(m_mpz, m, pp->ctx);

    mpz_t neg_m_mpz; /* -m */
    fq_t m_neg; fq_init(m_neg, pp->ctx);
    fq_neg(m_neg, m, pp->ctx);
    mpz_init(neg_m_mpz);
    fq2mpz(neg_m_mpz, m_neg, pp->ctx);

    element_t x;
    element_init_G1(x, pp->pairing);
    element_pow_mpz(x, pp->g, neg_m_mpz);
    element_mul(Y, x, R);
    
    // Set alpha^{-1}
    fq_t alpha_inv;
    fq_init(alpha_inv, pp->ctx);
    fq_inv(alpha_inv, sk->alpha, pp->ctx);
    mpz_t alpha_inv_m;
    mpz_init(alpha_inv_m);
    fq2mpz(alpha_inv_m, alpha_inv, pp->ctx);

    // Set Y[id]
    element_pow_mpz(sig->Y[id], Y, alpha_inv_m);

    // Free memory
    fq_clear(alpha_inv, pp->ctx);
    element_clear(x);
    element_clear(Y);
    element_clear(R);
    mpz_clear(m_mpz); 
    mpz_clear(neg_m_mpz);
}

void eval(Poly* f, Tag* sig, Tag* sig_out, uint64_t* id_set_eval, PublicPara* pp) {
    size_t i, j;
    Tag t1, t2, res;
    // Need to initialize t1, t2
    tag_init(&t1, pp->n, pp->ctx, pp->pairing);
    tag_init(&t2, pp->n, pp->ctx, pp->pairing);
    tag_init(&res, pp->n, pp->ctx, pp->pairing);

    // First calculate id_set, that is, let
    // id_set[id] := 1, if id \in some sig
    for (i = 0; i < f->t; i++) {
        id_set_eval[sig[i].id] = 1;
    }

    // Evaluate the x_i part
    for (i = 0; i < f->t; i++) {
        GTE_dot(&t1, f->eta[i], &sig[i], pp);
        GTE_add(&res, &res, &t1, pp);
    }
    
    // Evaluate the x_i x_j part
    for (i = 0; i < f->t; i++) {
        for (j = 0; j < f->t; j++) {
            GTE_cross(&t1, &sig[i], &sig[j], pp);
            GTE_dot(&t2, f->xi[j + f->t * i], &t1, pp);
            GTE_add(&res, &res, &t2, pp);
        }
    }

    tag_init(sig_out, pp->n, pp->ctx, pp->pairing);
    tag_copy(sig_out, &res, pp->n, pp->ctx);

    tag_clear(&t1, pp->n, pp->ctx);
    tag_clear(&t2, pp->n, pp->ctx);
    tag_clear(&res, pp->n, pp->ctx);
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

    // For y in mpz
    mpz_t y1, y2;
    mpz_init(y1); mpz_init(y2);
    fq2mpz(y1, sig1->y, pp->ctx); fq2mpz(y2, sig2->y, pp->ctx);

    // Calculate Y_r
    for (r = 0; r < n; r++) {
        element_pow2_mpz(res->Y[r], sig1->Y[r], y2, sig2->Y[r], y1);
    }

    // Calculate Z_(r,s)
    i = 0;
    for (r = 1; r <= n; r++) {
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
    mpz_clear(y1); mpz_clear(y2);
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
    i = 0;
    for (r = 1; r <= n; r++) {
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
        element_mul(res->Y[r], sig1->Y[r], sig2->Y[r]);
    }

    i = 0;
    for (r = 1; r <= n; r++) {
        for (s = r; s <= n; s++) {
            element_mul(res->Z[i], sig1->Z[i], sig2->Z[i]);
            i++;
        }
    }
}

int eff_ver(Poly* omega, Label* l, uint64_t Delta, VerKey* vk, fq_t m, uint64_t* id_set_eval,
        Tag* sigma, size_t len_t, PublicPara* pp) {

    int r1 = 0, r2 = 1, r3 = 0;
    r1 = fq_equal(m, sigma->y, pp->ctx);

    for (size_t i = 0; i < len_t; i++) {
        if (id_set_eval[l[i].id] == 0) {
            r2 = 0;
            break;
        }
    }

    element_t W1;
    element_init_GT(W1, pp->pairing);
    cf_eval_on(vk, Delta, omega, W1, pp);

    // Calculate the right hand side
    element_t W2;
    element_init_GT(W2, pp->pairing);
    get_W(sigma, vk, W2, pp);

    r3 = !(element_cmp(W1, W2));
    // element_printf("%B\t%B", W1, W2);

    // Free memory
    element_clear(W1);
    element_clear(W2);
    
    return r1 * r2 * r3;
}

void cf_eval_off(VerKey* K, Label* l, Poly* f, size_t** id_set, size_t* id_t_list, Poly* omega_f, PublicPara* pp) {
    size_t i, j, r, s;
    size_t n = pp->n;
    size_t t = f->t;
    omega_f->t = 2 * n;
    omega_f->xi = (fq_t *) malloc(sizeof(fq_t) * n * n * 4);
    omega_f->eta = (fq_t *) malloc(sizeof(fq_t) * 2 * n);
    
    // Use the pseudo random function
    fq_t u[f->t];
    fq_t v[f->t];
    // Temporary variable
    mpz_t u1, v1;
    mpz_init(u1); mpz_init(v1);

    // Calculate keys
    for (i = 0; i < f->t; i++) {
        fq_init(u[i], pp->ctx);
        fq_init(v[i], pp->ctx);
        PRF_F(u1, v1, K[id_t_list[i]].K->k1, (uint8_t*) &l[i], 16, pp->p);
        mpz2fq(u[i], u1, pp->ctx);
        mpz2fq(v[i], v1, pp->ctx);
    }

    fq_t temp, temp1;
    fq_init(temp, pp->ctx);
    fq_init(temp1, pp->ctx);
    
    // For omega_f eta
    for (s = 0; s < n; s++) {
        fq_zero(temp1, pp->ctx);

        // For Z_{id_i}, with id_i = s,  
        for (i = 0; i < t / n; i++) {
            fq_mul(temp, u[id_set[s][i]], f->eta[id_set[s][i]], pp->ctx);
            fq_add(temp1, temp1, temp, pp->ctx);
        }
        fq_init(omega_f->eta[s], pp->ctx);
        fq_set(omega_f->eta[s], temp1, pp->ctx);

        fq_zero(temp1, pp->ctx);
        for (i = 0; i < t / n; i++) {
            fq_mul(temp, v[id_set[s][i]], f->eta[id_set[s][i]], pp->ctx);
            fq_add(temp1, temp1, temp, pp->ctx);
        }
        fq_init(omega_f->eta[s + n], pp->ctx);
        fq_set(omega_f->eta[n + s], temp1, pp->ctx);
        // printf("s: %d", s);
    }  

    for (r = 0; r < 2 * n; r++) {
        for (s = 0; s < 2 * n; s++) {
            fq_init(omega_f->xi[2 * r * n + s], pp->ctx);
            fq_zero(omega_f->xi[2 * r * n + s], pp->ctx);
        }
    }

    // For omega_f xi
    for (r = 0; r < n; r++) {
        for (s = 0; s < n; s++) {

            for (i = 0; i < t / n; i++) {
                for (j = 0; j < t / n; j++) {
                    // Consider the set id_i = r and id_j = s;
                    size_t id1, id2;
                    id1 = id_set[r][i];
                    id2 = id_set[s][j];
                    fq_mul(temp, u[id1], u[id2], pp->ctx);
                    fq_mul(temp, temp, f->xi[id1 * t + id2], pp->ctx);
                    fq_add(omega_f->xi[r * 2 * n + s], omega_f->xi[r * 2 * n + s], temp, pp->ctx);
                    
                    fq_mul(temp, u[id1], v[id2], pp->ctx);
                    fq_mul(temp, temp, f->xi[id1 * t + id2], pp->ctx);
                    fq_add(omega_f->xi[r * 2 * n + s + n], omega_f->xi[r * 2 * n + s + n], temp, pp->ctx);

                    fq_mul(temp, v[id1], u[id2], pp->ctx);
                    fq_mul(temp, temp, f->xi[id1 * t + id2], pp->ctx);
                    fq_add(omega_f->xi[(r + n) * 2 * n + s], omega_f->xi[(r + n) * 2 * n + s], temp, pp->ctx);

                    fq_mul(temp, v[id1], v[id2], pp->ctx);
                    fq_mul(temp, temp, f->xi[id1 * t + id2], pp->ctx);
                    fq_add(omega_f->xi[(r + n) * 2 * n + s + n], omega_f->xi[(r + n) * 2 * n + s + n], temp, pp->ctx);
                }
            }
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

void cf_eval_on(VerKey* K, uint64_t Delta, Poly* omega_f, element_t W, PublicPara* pp) {
    size_t j;
    size_t n = pp->n;

    fq_t a[n * 2];
    // fq_t b[pp->n];
    // Temporary variable
    mpz_t a1, b1;
    mpz_init(a1); mpz_init(b1);
    
    for (j = 0; j < n; j++) {
        fq_init(a[j], pp->ctx);
        fq_init(a[j + n], pp->ctx);

        PRF_F(a1, b1, K[j].K->k2, (uint8_t *) &Delta, 8, pp->p);
        mpz2fq(a[j], a1, pp->ctx);
        mpz2fq(a[j + n], b1, pp->ctx);
    }
    
    // Evaluation
    fq_t w;
    fq_init(w, pp->ctx);
    poly_eval(omega_f, a, w, pp->ctx);

    // Result in mpz
    mpz_t w_m;
    mpz_init(w_m);
    fq2mpz(w_m, w, pp->ctx);

    element_pow_mpz(W, pp->gt, w_m);

    // free
    mpz_clear(a1); mpz_clear(b1);
    for (j = 0; j < 2 * n; j++) {
        fq_clear(a[j], pp->ctx);
    }
}

void GPE_naive(VerKey* vk, Poly* f, uint64_t Delta, Label* l, size_t* id_t_list,  element_t W, PublicPara* pp) {
    size_t t = f->t;
    size_t i;

    fq_t x_in[t];
    fq_t res;
    fq_init(res, pp->ctx);
    mpz_t res_mpz;
    mpz_init(res_mpz);

    mpz_t u, v, a, b;
    // allocate the memory
    mpz_init(u); mpz_init(v);
    mpz_init(a); mpz_init(b);
    fq_t u1, v1, a1, b1, t1, t2;
    fq_init(u1, pp->ctx); fq_init(v1, pp->ctx);
    fq_init(a1, pp->ctx); fq_init(b1, pp->ctx);
    fq_init(t1, pp->ctx); fq_init(t2, pp->ctx);


    // run the pseudo random function
    for (i = 0; i < t; i++) {
        fq_init(x_in[i], pp->ctx);
        // F(vk[id_t_list[i]].K, (uint8_t*) &Delta, (uint8_t*) (l + i), pp->g, x_in[i], pp->p);
        PRF_F(u, v, vk[id_t_list[i]].K->k1, (uint8_t*) &l[i], 16, pp->p);
        PRF_F(a, b, vk[id_t_list[i]].K->k2, (uint8_t*) &Delta, 8, pp->p);
        // gmp_printf("t: %d. u, v, a, b: %Zd, %Zd, %Zd, %Zd\n", i, u, v, a, b);

        mpz2fq(u1, u, pp->ctx); mpz2fq(v1, v, pp->ctx);
        mpz2fq(a1, a, pp->ctx); mpz2fq(b1, b, pp->ctx);
        fq_mul(t1, u1, a1, pp->ctx);
        fq_mul(t2, v1, b1, pp->ctx);
        fq_add(x_in[i], t1, t2, pp->ctx);
    }
    
    poly_eval(f, x_in, res, pp->ctx);
    fq2mpz(res_mpz, res, pp->ctx);
    
    element_init_GT(W, pp->pairing);
    element_pow_mpz(W, pp->gt, res_mpz);

    // Free memory
    fq_clear(res, pp->ctx);
    mpz_clear(res_mpz);
    mpz_clear(u); mpz_clear(v);
    mpz_clear(a); mpz_clear(b);
}

void get_W(Tag* sigma, VerKey* vk, element_t W, PublicPara* pp) {
    size_t i, r, s;
    element_t t2, t3;
    element_init_GT(t2, pp->pairing);
    element_init_GT(t3, pp->pairing);

    fq_t f1;
    fq_init(f1, pp->ctx);
    mpz_t b1;
    mpz_init(b1);

    mpz_t y_mpz;
    mpz_init(y_mpz);
    fq2mpz(y_mpz, sigma->y, pp->ctx);
    element_pow_mpz(W, pp->gt, y_mpz);    // gt^y

    for (r = 0; r < pp->n; r++) {
        element_pairing(t2, sigma->Y[r], pp->g);
        fq2mpz(b1, vk[r].alpha, pp->ctx);
        element_pow_mpz(t3, t2, b1);
        element_mul(W, t3, W);
    }

    i = 0;
    for (r = 1; r <= pp->n; r++) {
        for (s = r; s <= pp->n; s++) {
            fq_mul(f1, vk[r-1].alpha, vk[s-1].alpha, pp->ctx);
            fq2mpz(b1, f1, pp->ctx);
            element_pow_mpz(t3, sigma->Z[i], b1);
            element_mul(W, t3, W);
            i++;
        }
    }

    element_clear(t2);
    element_clear(t3);
    fq_clear(f1, pp->ctx);
    mpz_clear(b1);
}


void mkha_clear(PublicPara* pp, VerKey* vk) {
    // Clear vk
    for (size_t i = 0; i < pp->n; i++) {
        key_clear(vk[i].K);
        free(vk[i].K);
        fq_clear(vk[i].alpha, pp->ctx);
    }

    // Clear pp
    element_clear(pp->g);
    element_clear(pp->gt);
    pairing_clear(pp->pairing);
    mpz_clear(pp->p);
    fmpz_clear(pp->pf);
    fq_ctx_clear(pp->ctx);
}
