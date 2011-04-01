#ifndef TESTING_H
#define TESTING_H

#include "misc.h"
#include "packedmatrix.h"

mzd_t *mzd_init_test_matrix_random(rci_t M, rci_t N, rci_t m, rci_t n, rci_t offset, word pattern, mzd_t **A);
mzd_t *mzd_free_test_matrix_random(mzd_t *A, mzd_t *a);
int mzd_check_pattern(mzd_t *A, rci_t m, rci_t n, rci_t offset, word pattern);

#endif //TESTING_H
