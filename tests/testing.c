#include "testing.h"

mzd_t *mzd_init_test_matrix_random(rci_t M, rci_t N, rci_t m, rci_t n, word pattern, mzd_t **A) {
  *A = mzd_init(M, N);

  for (rci_t i = 0; i < M; i++) {
    for (rci_t j = 0; j < (*A)->width; j++) { (*A)->rows[i][j] = pattern; }
  }

  mzd_t *a = mzd_init_window(*A, 0, 0, m, n);
  mzd_randomize(a);

  return a;
}

void mzd_free_test_matrix_random(mzd_t *A, mzd_t *a) {
  mzd_free(a);
  mzd_free(A);
}

int mzd_check_pattern(mzd_t *A, rci_t m, rci_t n, word pattern) {

  for (rci_t i = 0; i < A->nrows; i++) {
    if (i >= m) {
      for (rci_t j = 0; j < A->width; j++)
        if (A->rows[i][j] ^ pattern) { return 1; }
    } else {
      if ((A->rows[i][n / m4ri_radix] ^ pattern) & ~A->high_bitmask) return 1;

      for (rci_t j = n / m4ri_radix + 1; j < A->width; j++)
        if (A->rows[i][j] ^ pattern) { return 1; }
    }
  }
  return 0;
}
