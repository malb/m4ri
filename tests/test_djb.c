#include "testing.h"
#include <m4ri/config.h>
#include <m4ri/djb.h>
#include <m4ri/m4ri.h>
#include <stdlib.h>

/**
 * Check that the results of the implementation of Dan Bernstein's "Optimizing
 * linear maps mod 2" matches the matrix multiplication algorithms.
 *
 * \param m Number of rows of A
 * \param l Number of columns of A/number of rows of B
 * \param n Number of columns of B
 */
int mul_test_equality(rci_t m, rci_t l, rci_t n) {
  int ret = 0;
  printf("   mul: m: %4d, l: %4d, n: %4d", m, l, n);

  /* we create two random matrices */
  mzd_t *A = mzd_init(m, l);
  mzd_t *B = mzd_init(l, n);
  mzd_randomize(A);
  mzd_randomize(B);

  /* C = A*B via Strassen */
  mzd_t *C = mzd_mul(NULL, A, B, __M4RI_STRASSEN_MUL_CUTOFF);

  /* C = A*B via DJB */
  djb_t *djb_A = djb_compile(A);
  mzd_t *djb_C = mzd_init(m, n);
  djb_apply_mzd(djb_A, djb_C, B);

  if (mzd_equal(C, djb_C) != TRUE) {
    printf(" Strassen != DJB");
    ret -= 1;
  }

  mzd_free(djb_C);
  djb_free(djb_A);

  mzd_free(C);
  mzd_free(B);
  mzd_free(A);

  if (ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }

  return ret;
}

int main() {
  int status = 0;

  srandom(17);

  status += mul_test_equality(1, 1, 1);
  status += mul_test_equality(1, 128, 128);
  status += mul_test_equality(3, 131, 257);
  status += mul_test_equality(64, 64, 64);
  status += mul_test_equality(128, 128, 128);
  status += mul_test_equality(21, 171, 31);
  status += mul_test_equality(21, 171, 31);
  status += mul_test_equality(193, 65, 65);
  status += mul_test_equality(1025, 1025, 1025);
  status += mul_test_equality(2048, 2048, 4096);
  status += mul_test_equality(4096, 3528, 4096);
  status += mul_test_equality(1024, 1025, 1);
  status += mul_test_equality(1000, 1000, 1000);
  status += mul_test_equality(1000, 10, 20);
  status += mul_test_equality(1710, 1290, 1000);
  status += mul_test_equality(1290, 1710, 200);
  status += mul_test_equality(1290, 1710, 2000);
  status += mul_test_equality(1290, 1290, 2000);
  status += mul_test_equality(1000, 210, 200);

  if (status == 0) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
