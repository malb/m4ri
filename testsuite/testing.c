#include "testing.h"

mzd_t *mzd_init_test_matrix_random(rci_t M, rci_t N, rci_t m, rci_t n, rci_t offset, word pattern, mzd_t **A) {
  *A = mzd_init(M, N);

  for(rci_t i=0; i<M; i++) {
    for(rci_t j=0; j<(*A)->width; j++) {
      (*A)->rows[i][j] = pattern;
    }
  }

  mzd_t* a = mzd_init_window(*A, offset, offset, offset + m, offset + n);
  mzd_randomize(a);

  return a;
}

mzd_t *mzd_free_test_matrix_random(mzd_t *A, mzd_t *a) {
  mzd_free(a);
  mzd_free(A);
};

int mzd_check_pattern(mzd_t *A, rci_t m, rci_t n, rci_t offset, word pattern) {

  for(rci_t i=0; i<A->nrows; i++) {
    if (i<offset || i >= m+offset) {

      for(rci_t j=0; j<A->width; j++) 
        if(A->rows[i][j] ^ pattern) {
          return 1;
        }

    } else {

      for(rci_t j=0; j < (offset/RADIX); j++)
        if(A->rows[i][j] ^ pattern) {
          return 1;
        }

      if ( (offset/RADIX) == (offset+n)/RADIX ) {
        word const mask = ~MIDDLE_BITMASK(m, offset%RADIX);
        if( (A->rows[i][offset/RADIX] ^ pattern) & mask ) {
          return 1;
        }

      } else {
        word const mask_begin = ~RIGHT_BITMASK(RADIX - offset%RADIX);
        word const mask_end = ~LEFT_BITMASK((n + offset) % RADIX);

        if( (A->rows[i][offset/RADIX] ^ pattern) & mask_begin ) {
          return 1;
        }

        if( (A->rows[i][(offset+n)/RADIX] ^ pattern) & mask_end ) {
          return 1;       
        }
      }
      
      for(rci_t j=(offset+n)/RADIX+1; j<A->width; j++) 
        if(A->rows[i][j] ^ pattern) {
          return 1;
        }
    } 

  }
  return 0;
}
