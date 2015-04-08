#ifndef TESTING_H
#define TESTING_H

#include <m4ri/misc.h>
#include <m4ri/mzd.h>

mzd_t *mzd_init_test_matrix_random(rci_t M, rci_t N, rci_t m, rci_t n, word pattern, mzd_t **A);
void mzd_free_test_matrix_random(mzd_t *A, mzd_t *a);
int mzd_check_pattern(mzd_t *A, rci_t m, rci_t n, word pattern);

#endif //TESTING_H
