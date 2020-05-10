#include "testing.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>
#include <stdlib.h>

/**
 * Check that the results of all implemented multiplication algorithms
 * match up.
 *
 * \param m Number of rows of A
 * \param l Number of columns of A/number of rows of B
 * \param n Number of columns of B
 * \param k Parameter k of M4RM algorithm, may be 0 for automatic choice.
 * \param cutoff Cut off parameter at which dimension to switch from
 * Strassen to M4RM
 */
int mul_test_equality(rci_t m, rci_t l, rci_t n, int k, int cutoff) {
  int ret = 0;
  printf("   mul: m: %4d, l: %4d, n: %4d, k: %2d, cutoff: %4d", m, l, n, k, cutoff);

  /* we create two random matrices */
  mzd_t *A = mzd_init(m, l);
  mzd_t *B = mzd_init(l, n);
  mzd_randomize(A);
  mzd_randomize(B);

  /* C = A*B via Strassen */
  mzd_t *C = mzd_mul(NULL, A, B, cutoff);

  /* D = A*B via M4RM, temporary buffers are managed internally */
  mzd_t *D = mzd_mul_m4rm(NULL, A, B, k);

  if (mzd_equal(C, D) != TRUE) {
    printf(" Strassen != M4RM");
    ret -= 1;
  }

  /* E = A*B via naive cubic multiplication */
  mzd_t *E = mzd_mul_naive(NULL, A, B);

  if (mzd_equal(D, E) != TRUE) {
    printf(" M4RM != Naiv");
    ret -= 1;
  }

  if (mzd_equal(C, E) != TRUE) {
    printf(" Strassen != Naiv");
    ret -= 1;
  }

#if __M4RI_HAVE_OPENMP
  mzd_t *F = mzd_mul_mp(NULL, A, B, cutoff);
  if (mzd_equal(C, F) != TRUE) {
    printf(" MP != Naiv");
    ret -= 1;
  }
  mzd_free(F);
#endif

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  mzd_free(D);
  mzd_free(E);

  if (ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }

  return ret;
}

/**
 * Check that the results of all implemented squaring algorithms match
 * up.
 *
 * \param m Number of rows and columns of A
 * \param k Parameter k of M4RM algorithm, may be 0 for automatic choice.
 * \param cutoff Cut off parameter at which dimension to switch from
 * Strassen to M4RM
 */
int sqr_test_equality(rci_t m, int k, int cutoff) {
  int ret = 0;
  mzd_t *A, *C, *D, *E;

  printf("   sqr: m: %4d, k: %2d, cutoff: %4d", m, k, cutoff);

  /* we create one random matrix */
  A = mzd_init(m, m);
  mzd_randomize(A);

  /* C = A*A via Strassen */
  C = mzd_mul(NULL, A, A, cutoff);

  /* D = A*A via M4RM, temporary buffers are managed internally */
  D = mzd_mul_m4rm(NULL, A, A, k);

  /* E = A*A via naive cubic multiplication */
  E = mzd_mul_naive(NULL, A, A);

  mzd_free(A);

  if (mzd_equal(C, D) != TRUE) {
    printf(" Strassen != M4RM");
    ret -= 1;
  }

  if (mzd_equal(D, E) != TRUE) {
    printf(" M4RM != Naiv");
    ret -= 1;
  }

  if (mzd_equal(C, E) != TRUE) {
    printf(" Strassen != Naiv");
    ret -= 1;
  }

  mzd_free(C);
  mzd_free(D);
  mzd_free(E);

  if (ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }

  return ret;
}

int addmul_test_equality(rci_t m, rci_t l, rci_t n, int k, int cutoff) {
  int ret = 0;
  printf("addmul: m: %4d, l: %4d, n: %4d, k: %2d, cutoff: %4d", m, l, n, k, cutoff);

  /* we create two random matrices */
  mzd_t *A = mzd_init(m, l);
  mzd_t *B = mzd_init(l, n);
  mzd_t *C = mzd_init(m, n);
  mzd_randomize(A);
  mzd_randomize(B);
  mzd_randomize(C);

  /* D = C + A*B via M4RM, temporary buffers are managed internally */
  mzd_t *D = mzd_copy(NULL, C);
  D        = mzd_addmul_m4rm(D, A, B, k);

  /* E = C + A*B via naiv cubic multiplication */
  mzd_t *E = mzd_mul_m4rm(NULL, A, B, k);
  mzd_add(E, E, C);

  if (mzd_equal(D, E) != TRUE) {
    printf(" M4RM != add,mul");
    ret -= 1;
  }

  /* F = C + A*B via naiv cubic multiplication */
  mzd_t *F = mzd_copy(NULL, C);
  F        = mzd_addmul(F, A, B, cutoff);

  if (mzd_equal(E, F) != TRUE) {
    printf(" add,mul = addmul");
    ret -= 1;
  }
  if (mzd_equal(F, D) != TRUE) {
    printf(" M4RM != addmul");
    ret -= 1;
  }

#if __M4RI_HAVE_OPENMP
  mzd_t *G = mzd_copy(NULL, C);
  G        = mzd_addmul_mp(G, A, B, cutoff);
  if (mzd_equal(D, G) != TRUE) {
    printf(" MP != Naiv");
    ret -= 1;
  }
  mzd_free(G);
#endif

  if (ret == 0)
    printf(" ... passed\n");
  else
    printf(" ... FAILED\n");

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  mzd_free(D);
  mzd_free(E);
  mzd_free(F);
  return ret;
}

int addsqr_test_equality(rci_t m, int k, int cutoff) {
  int ret = 0;
  mzd_t *A, *C, *D, *E, *F;

  printf("addsqr: m: %4d, k: %2d, cutoff: %4d", m, k, cutoff);

  /* we create two random matrices */
  A = mzd_init(m, m);
  C = mzd_init(m, m);
  mzd_randomize(A);
  mzd_randomize(C);

  /* D = C + A*B via M4RM, temporary buffers are managed internally */
  D = mzd_copy(NULL, C);
  D = mzd_addmul_m4rm(D, A, A, k);

  /* E = C + A*B via naive cubic multiplication */
  E = mzd_mul_m4rm(NULL, A, A, k);
  mzd_add(E, E, C);

  /* F = C + A*B via naive cubic multiplication */
  F = mzd_copy(NULL, C);
  F = mzd_addmul(F, A, A, cutoff);

  mzd_free(A);
  mzd_free(C);

  if (mzd_equal(D, E) != TRUE) {
    printf(" M4RM != add,mul");
    ret -= 1;
  }
  if (mzd_equal(E, F) != TRUE) {
    printf(" add,mul = addmul");
    ret -= 1;
  }
  if (mzd_equal(F, D) != TRUE) {
    printf(" M4RM != addmul");
    ret -= 1;
  }

  if (ret == 0)
    printf(" ... passed\n");
  else
    printf(" ... FAILED\n");

  mzd_free(D);
  mzd_free(E);
  mzd_free(F);
  return ret;
}

int main() {
  int status = 0;

  srandom(17);

  status += mul_test_equality(1, 1, 1, 0, 1024);
  status += mul_test_equality(1, 128, 128, 0, 0);
  status += mul_test_equality(3, 131, 257, 0, 0);
  status += mul_test_equality(64, 64, 64, 0, 64);
  status += mul_test_equality(128, 128, 128, 0, 64);
  status += mul_test_equality(21, 171, 31, 0, 63);
  status += mul_test_equality(21, 171, 31, 0, 131);
  status += mul_test_equality(193, 65, 65, 8, 64);
  status += mul_test_equality(1025, 1025, 1025, 3, 256);
  status += mul_test_equality(2048, 2048, 4096, 0, 1024);
  status += mul_test_equality(4096, 3528, 4096, 0, 1024);
  status += mul_test_equality(1024, 1025, 1, 0, 1024);
  status += mul_test_equality(1000, 1000, 1000, 0, 256);
  status += mul_test_equality(1000, 10, 20, 0, 64);
  status += mul_test_equality(1710, 1290, 1000, 0, 256);
  status += mul_test_equality(1290, 1710, 200, 0, 64);
  status += mul_test_equality(1290, 1710, 2000, 0, 256);
  status += mul_test_equality(1290, 1290, 2000, 0, 64);
  status += mul_test_equality(1000, 210, 200, 0, 64);

  status += addmul_test_equality(1, 128, 128, 0, 0);
  status += addmul_test_equality(3, 131, 257, 0, 0);
  status += addmul_test_equality(64, 64, 64, 0, 64);
  status += addmul_test_equality(128, 128, 128, 0, 64);
  status += addmul_test_equality(21, 171, 31, 0, 63);
  status += addmul_test_equality(21, 171, 31, 0, 131);
  status += addmul_test_equality(193, 65, 65, 8, 64);
  status += addmul_test_equality(1025, 1025, 1025, 3, 256);
  status += addmul_test_equality(4096, 4096, 4096, 0, 2048);
  status += addmul_test_equality(1000, 1000, 1000, 0, 256);
  status += addmul_test_equality(1000, 10, 20, 0, 64);
  status += addmul_test_equality(1710, 1290, 1000, 0, 256);
  status += addmul_test_equality(1290, 1710, 200, 0, 64);
  status += addmul_test_equality(1290, 1710, 2000, 0, 256);
  status += addmul_test_equality(1290, 1290, 2000, 0, 64);
  status += addmul_test_equality(1000, 210, 200, 0, 64);

  status += sqr_test_equality(1, 0, 1024);
  status += sqr_test_equality(128, 0, 0);
  status += sqr_test_equality(131, 0, 0);
  status += sqr_test_equality(64, 0, 64);
  status += sqr_test_equality(128, 0, 64);
  status += sqr_test_equality(171, 0, 63);
  status += sqr_test_equality(171, 0, 131);
  status += sqr_test_equality(193, 8, 64);
  status += sqr_test_equality(1025, 3, 256);
  status += sqr_test_equality(2048, 0, 1024);
  status += sqr_test_equality(3528, 0, 1024);
  status += sqr_test_equality(1000, 0, 256);
  status += sqr_test_equality(1000, 0, 64);
  status += sqr_test_equality(1710, 0, 256);
  status += sqr_test_equality(1290, 0, 64);
  status += sqr_test_equality(2000, 0, 256);
  status += sqr_test_equality(2000, 0, 64);
  status += sqr_test_equality(210, 0, 64);

  status += addsqr_test_equality(1, 0, 0);
  status += addsqr_test_equality(131, 0, 0);
  status += addsqr_test_equality(64, 0, 64);
  status += addsqr_test_equality(128, 0, 64);
  status += addsqr_test_equality(171, 0, 63);
  status += addsqr_test_equality(171, 0, 131);
  status += addsqr_test_equality(193, 8, 64);
  status += addsqr_test_equality(1025, 3, 256);
  status += addsqr_test_equality(4096, 0, 2048);
  status += addsqr_test_equality(1000, 0, 256);
  status += addsqr_test_equality(1000, 0, 64);
  status += addsqr_test_equality(1710, 0, 256);
  status += addsqr_test_equality(1290, 0, 64);
  status += addsqr_test_equality(2000, 0, 256);
  status += addsqr_test_equality(2000, 0, 64);
  status += addsqr_test_equality(210, 0, 64);

  if (status == 0) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
