#ifndef PRF_H
#define PRF_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <inttypes.h>
#include <pbc/pbc.h>
// #include <gcrypt.h>
#include <time.h>

/* Defines the key length in bytes for key in key space */
#define KEY_LEN 32
/* We use SHA256 */
#define MD_MAP SH256
#define BIN_SIZE 16

/* Represent the two keys chosen in KG */
typedef struct {
    uint8_t* k1;
    uint8_t* k2;
} Key;

/**
 * The key generation algorithm
 * 
 * @param[out] K            Secret key
 * 
 */
void KG(Key* K);

/**
 * The algorithm computes (u, v) <-- F(k, l), 
 * (a, v) <-- F(k', dalta), and out put R = g^(ua+bv). 
 * 
 * @param[in] K             Secret Key
 * @param[out] r            Output
 * 
 */
void F(Key* K, uint8_t* delta, uint8_t* l, element_t g, element_t r, mpz_t p);

/**
 *  The pseudo random function that
 *      K \cross {0, 1}* \to Z_p^2.
 */
void PRF_F(mpz_t r1, mpz_t r2, uint8_t* k, uint8_t* data, size_t data_size, mpz_t p);

/* Clear the memory allocated for K */
void key_clear(Key* K);

#endif