#define FP_PRIME 1600
#include "mkha.h"
#include <time.h>

/****** DEBUG SECTION ******/
/* THIS SECTION IS FOR DEBUG PURPOSE */
void print_id_list(uint64_t* ids, size_t len) {
    printf("ID_list: ");
    for (size_t i = 0; i < len; i++) {
        printf("%ld: %lu, ", i, ids[i]);
    }
    printf("\n");
}

void print_key(VerKey* vk, size_t len) {
    printf("Keys: ");
    for (size_t i = 0; i < len; i++) {
        printf("%ld: %lu, ", i, (uint64_t) vk[i].K->k1);
    }
    printf("\n");
}   


/***** END OF DEBUG SECTION *****/

// t is #message, n is #user, k is #data_set
// delta is the dataset
// We ENFORCE t be a multiple of n
void try(uint32_t lambda, size_t t, size_t n, size_t k) {
    double eval_time = .0;
    double off_ver_time = .0;
    double on_ver_time = .0;
    double total_ver_time = .0;

    clock_t eval_start, eval_end;
    clock_t off_start, off_end;
    clock_t on_start, on_end;

    // Set t be a multiple of n
    t = t / n * n;
    size_t user_data_len = t / n;

    // For loops
    uint64_t delta;
    size_t i, j; 

    // The id_t_list[i] is the id of the i-th message
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

    // Now allocate the needed components
    Label* l = (Label *) malloc(sizeof(Label) * t);
    Tag* s = (Tag *) malloc(sizeof(Tag) * t);
    fq_t* m = (fq_t *) malloc(sizeof(fq_t) * t);
    Poly *f = (Poly *) malloc(sizeof(Poly)); 
    Tag* s_out = (Tag *) malloc(sizeof(Tag));

    // Memory for the message list
    for (i = 0; i < t; i++) {
        fq_init(m[i], pp->ctx);
    }

    Poly* omega_f = (Poly *) malloc(sizeof(Poly));
    fq_t m_result;
    fq_init(m_result, pp->ctx);
    
    // Generate the label first
    // The id_set[id][tau] is the index in the message list
    size_t** id_set = (size_t **) malloc(sizeof(size_t*) * n);
    for (i = 0; i < n; i++) {
        id_set[i] = (size_t *) malloc(sizeof(size_t) * user_data_len);
        for (j = 0; j < user_data_len; j++) {
            // Set the i
            id_set[i][j] = user_data_len * i + j;
            id_t_list[user_data_len * i + j] = i;

            // Set the label
            l[user_data_len * i + j].id = i;
            l[user_data_len * i + j].tau = j;
        }
    }

    // Generate a random quadratic function
    poly_rand_init(f, t, pp->ctx);
    // poly_print(f, pp->ctx);

    // Now do Offline VerPrep
    off_start = clock();
    cf_eval_off(vk, l, f, id_set, id_t_list, omega_f, pp);
    off_end = clock();
    off_ver_time = (double) (off_end - off_start) / CLOCKS_PER_SEC;

    // print the omega_f
    // poly_print(omega_f, pp->ctx);

    for (size_t k_count = 0; k_count < k; k_count++) {
        delta = k_count * 23 + 315;

        for (i = 0; i < t; i++) {
            // Generate the message
            fq_rand(m[i], state, pp->ctx);
            auth(&(vk[id_t_list[i]]), delta, &l[i], m[i], s + i, pp);
        }
        
        // Evaulate
        uint64_t* id_set_eval = (uint64_t *) calloc(n, sizeof(uint64_t));
        eval(f, s, s_out, id_set_eval, pp);


        // Compute
        eval_start = clock();
        poly_eval(f, m, m_result, pp->ctx);
        eval_end = clock();
        eval_time += (double) (eval_end - eval_start) / CLOCKS_PER_SEC;

        // element_t W3;
        // GPE_naive(vk, f, delta, l, id_t_list, W3, pp);
        // element_printf("W is : %B\n", W3);

        // Verify
        on_start = clock();
        int result = eff_ver(omega_f, l, delta, vk, m_result, id_set_eval, s_out, f->t, pp);    
        on_end = clock();
        on_ver_time += (double) (on_end - on_start) / CLOCKS_PER_SEC;

        if (result != 1) printf("ERROR!!!\n");
        // printf("The result of verification is: %d.\n", result);
        // printf("The message is: "); fq_print_pretty(m_result, pp->ctx); printf("\n");
    }

    total_ver_time = off_ver_time + on_ver_time;
    printf("With n (#user) = %ld, t (#message) = %ld, k (#datasets) = %ld\n", n, t, k);
    printf("The total ver time: %f,\tOn time: %f,\t off time: %f.\n", total_ver_time, on_ver_time, off_ver_time);
    printf("The evaluation time required: %f.\n", eval_time);
    if (eval_time > total_ver_time) 
        printf("Eval time is greater than total ver time, LEGIT.\n");
    else 
        printf("Eval time is less, NOT so good!\n");

    // Free all memory
    // NOTE there is a memory leack caused by pp
    for (i = 0; i < t; i++) {
        fq_clear(m[i], pp->ctx);
    }
    flint_randclear(state);
    fq_clear(m_result, pp->ctx);
    poly_clear(f, pp->ctx);
    poly_clear(omega_f, pp->ctx);
    mkha_clear(pp, vk);
    free(s); free(s_out);
    free(m); free(l);
    free(vk);
    free(pp);
    free(f); free(omega_f);
}

int main() {
    // srand(time(NULL));

    size_t user_size = 2;
    size_t message_size = 55;
    size_t lambda = 64;
    size_t data_set_size = 40;
    try(lambda, message_size, user_size, data_set_size);

    return 0;
}