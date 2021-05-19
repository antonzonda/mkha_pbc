#define FP_PRIME 1600
#include "mkha.h"

// t is the #message, while n is #user
// delta is the dataset
void try(uint32_t lambda, size_t t, size_t n, uint64_t delta) {
    // For loops
    size_t i, j, id; 
    size_t* id_t_list = (size_t *) malloc(sizeof(size_t) * t);
    // FLINT state for generating randomness
    flint_rand_t state;
    flint_randinit(state); 

    PublicPara* pp = (PublicPara *) malloc(sizeof(PublicPara));

    // Setup algorithm 
    set_up(pp, n, lambda);

    // Key gen, that is, the i-th key indicates the i-th user
    VerKey* vk = (VerKey *) malloc(sizeof(VerKey) * n);
    for (i = 0; i < n; i++) {
        vk[i].id = i;
        key_gen(pp, i, &(vk[i].K));
    }

    // Now we randomly generate some data for each user
    // Then do the authentication or produce the tag
    // Note that we here supposedly have t = O(n^2)
    Label* l = (Label *) malloc(sizeof(l) * t);
    Tag* s = (Tag *) malloc(sizeof(Tag) * t);
    fq_t* m = (fq_t *) malloc(sizeof(fq_t) * t);
    for (i = 0; i < t; i++) {
        id = i / n;
        id_t_list[i] = id;
        l[i].id = id;
        l[i].tau = i % n;
        // Generate the message
        fq_init(m[i], pp->ctx);
        fq_rand(m[i], state, pp->ctx);
        
        // Generate the tag
        auth(&vk[id], delta, &l[i], m[i], &s[i], pp);
    }

    // Generate a random quadratic function
    Poly *f;
    poly_rand_init(f, t, pp->ctx);
    Tag* s_out = (Tag *) malloc(sizeof(Tag));
    
    // Evaulate
    // eval(f, s, s_out, );
    
    // Verify
    // int result = ver();

    printf("The result of verification is: ");

    // Free all memory
    for (i = 0; i < t; i++) {
        fq_clear(m[i], pp->ctx);
    }
    poly_clear(f, pp->ctx);
    flint_randclear(state);
    free(s); free(s_out);
    free(m); free(l);
    free(vk);
    free(pp);
}

int main() {
    mkha_init();
    try(50, 500, 50, 1);
    mkha_close();
    return 0;
}