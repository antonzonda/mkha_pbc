#include "prf.h"

void KG(Key* K) {
    K = (Key*) malloc(sizeof(Key));
    K->k1 = (uint8_t*) malloc(KEY_LEN);
    K->k2 = (uint8_t*) malloc(KEY_LEN);
    rand_bytes(K->k1, KEY_LEN);
    rand_bytes(K->k2, KEY_LEN);
}

void F(Key* K, uint8_t* delta, uint8_t* l, element_t g, element_t r);
    uint8_t mac1[RLC_MD_LEN];
    uint8_t mac2[RLC_MD_LEN];
    bn_t u, v, a, b;
    g1_t r1;

    // Allocate 
    bn_new(u); bn_new(v);
    bn_new(a); bn_new(b);

    md_hmac(mac1, delta, 8, K->k1, KEY_LEN);
    md_hmac(mac2, l, 16, K->k2, KEY_LEN);

    bn_read_bin(u, mac2, BIN_SIZE);
    bn_read_bin(v, mac2 + BIN_SIZE, BIN_SIZE);
    bn_read_bin(a, mac1, BIN_SIZE);
    bn_read_bin(b, mac1 + BIN_SIZE, BIN_SIZE);

    g1_new(r1);
    g1_mul(r1, g, a);
    g1_mul(r1, r1, u);
    g1_norm(r1, r1); /* Not sure what does it do */
    g1_mul(r, g, b);
    g1_mul(r, r, v);
    g1_add(r, r, r1);
    g1_norm(r, r);

    // free
    bn_free(u); bn_free(v);
    bn_free(a); bn_free(b);
}

void PRF_F(fq_t r1, fq_t r2, uint8_t* k, uint8_t* data, size_t data_size, fq_ctx_t ctx);
    uint8_t mac[RLC_MD_LEN];

    md_hmac(mac, data, data_size, k, KEY_LEN);
    bn_read_bin(b1, mac, BIN_SIZE);
    bn_read_bin(b2, mac + BIN_SIZE, BIN_SIZE);

    bn2fq(r1, b1, ctx);
    bn2fq(r2, b2, ctx);

}

void clear_key(Key* K) {
    free(K->k1);
    free(K->k2);
    free(K);
}