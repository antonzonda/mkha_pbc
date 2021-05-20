#include "prf.h"

void KG(Key* K) {
    // K = (Key*) malloc(sizeof(Key));
    K->k1 = (uint8_t*) malloc(KEY_LEN);
    K->k2 = (uint8_t*) malloc(KEY_LEN);
    srand(time(NULL));
    // We use rand to obtain a random bytestream
    for (int i = 0; i < KEY_LEN; i++) {
        K->k1[i] = (uint8_t) rand();
        K->k2[i] = (uint8_t) rand();
    }
}

void F(Key* K, uint8_t* delta, uint8_t* l, element_t g, element_t r, mpz_t p) {
    mpz_t u, v, a, b;
    element_t r1, r2;
    // allocate the memory
    mpz_init(u); mpz_init(v);
    mpz_init(a); mpz_init(b);
    element_init_same_as(r1, g);

    // run the pseudo random function
    PRF_F(u, v, K->k1, l, 16, p);
    PRF_F(a, b, K->k2, delta, 8, p);

    element_mul_mpz(r1, g, a);
    element_mul_mpz(r1, r1, u);

    element_mul_mpz(r, g, b);
    element_mul_mpz(r, r, v);
    element_add(r, r, r1);

    // free the memory
    mpz_clear(u); mpz_clear(v);
    mpz_clear(a); mpz_clear(b);
    element_clear(r1);
}

void PRF_F(mpz_t r1, mpz_t r2, uint8_t* k, uint8_t* data, size_t data_size, mpz_t p) {
    uint8_t hash[MD_SIZE];
    gcry_md_hd_t hd;

    gcry_md_open(&hd, HASH_ALG, GCRY_MD_FLAG_HMAC);
    gcry_md_setkey(hd, k, KEY_LEN);
    gcry_md_write(hd, (void*) data, data_size);
    gcry_md_extract(hd, HASH_ALG, hash, MD_SIZE);

    mpz_import(r1, BIN_SIZE, 1, 1, 0, 0, hash);
    mpz_import(r2, BIN_SIZE, 1, 1, 0, 0, hash + BIN_SIZE);

    mpz_mod(r1, r1, p);
    mpz_mod(r2, r2, p);

    gcry_md_close(hd);
}

void clear_key(Key* K) {
    free(K->k1);
    free(K->k2);
    free(K);
}