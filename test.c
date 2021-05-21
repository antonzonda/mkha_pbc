#define FP_PRIME 1600
#include "mkha.h"

/****** DEBUG SECTION ******/
/* THIS SECTION IS FOR DEBUG PURPOSE */
void print_id_list(uint64_t* ids, size_t len) {
    printf("ID_list: ");
    for (int i = 0; i < len; i++) {
        printf("%d: %u, ", i, ids[i]);
    }
    printf("\n");
}

void print_key(VerKey* vk, size_t len) {
    printf("Keys: ");
    for (int i = 0; i < len; i++) {
        printf("%d: %u, ", i, vk[i].K->k1);
    }
    printf("\n");
}   


/* END OF DEBUG SECTION */

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
        key_gen(pp, i, &(vk[i]));
    }
    // Now we randomly generate some data for each user
    // Then do the authentication or produce the tag
    // Note that we here supposedly have t = O(n^2)
    Label* l = (Label *) malloc(sizeof(Label) * t);
    Tag* s = (Tag *) malloc(sizeof(Tag) * t);
    fq_t* m = (fq_t *) malloc(sizeof(fq_t) * t);
    for (i = 0; i < t; i++) {
        id = i % n;
        id_t_list[i] = id;
        l[i].id = id;
        l[i].tau = i / n;
        // Generate the message
        fq_init(m[i], pp->ctx);
        fq_rand(m[i], state, pp->ctx);
        
        // Generate the tag
        auth(&(vk[id]), delta, &l[i], m[i], &s[i], pp);
    }
    // Generate a random quadratic function
    Poly *f;
    f = (Poly *) malloc(sizeof(Poly)); 
    poly_rand_init(f, t, pp->ctx);
    Tag* s_out = (Tag *) malloc(sizeof(Tag));
    
    // Evaulate
    uint64_t* id_set = (uint64_t *) calloc(n, sizeof(uint64_t));
    eval(f, s, s_out, id_set, pp);
    
    print_id_list(id_t_list, t);
    print_key(vk, n);
    // Verify
    int result = ver(f, l, delta, vk, s_out->y, id_set, id_t_list, s_out, pp);

    printf("The result of verification is: %d.\n", result);

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
    free(f);
}

int main() {

    try(50, 10, 5, 1);
    
    return 0;
}