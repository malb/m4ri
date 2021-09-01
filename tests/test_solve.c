#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "testing.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>
#include <stdlib.h>

int test_pluq_solve_left(rci_t m, rci_t n, int offsetA, int offsetB) {
  mzd_t *Abase = mzd_init(2048, 2048);
  mzd_t *Bbase = mzd_init(2048, 2048);
  mzd_randomize(Abase);
  mzd_randomize(Bbase);

  mzd_t *A = mzd_init_window(Abase, 0, offsetA, m, m + offsetA);
  mzd_t *B = mzd_init_window(Bbase, 0, offsetB, m, n + offsetB);

  // copy B
  mzd_t *Bcopy = mzd_copy(NULL, B);
  for (rci_t i = 0; i < m; ++i) { mzd_write_bit(A, i, i, 1); }

  mzd_t *Acopy = mzd_copy(NULL, A);
  rci_t r      = mzd_echelonize(Acopy, 1);
  printf("solve_left m: %4d, n: %4d, r: %4d da: %4d db: %4d ", m, n, r, offsetA, offsetB);
  mzd_free(Acopy);
  Acopy = mzd_copy(NULL, A);

  int consistency = mzd_solve_left(A, B, 0, 1);

  // copy B
  mzd_t *X = mzd_copy(NULL, B);
  mzd_t *B1 = mzd_mul(NULL, Acopy, X, 0);
  
  int status = 0;
  if (consistency == 0) {
    status = 1 - mzd_equal(Bcopy, B1);
    if (status == 0) {
      printf("passed\n");
    } else {
      printf("FAILED\n");
    }
  } else {
    printf("skipped (no solution)\n");
  }
  mzd_free(Bcopy);
  mzd_free(B1);
  mzd_free_window(A);
  mzd_free_window(B);
  mzd_free(Acopy);
  mzd_free(Abase);
  mzd_free(Bbase);
  mzd_free(X);
  return status;
}

int test_solve_left_random(rci_t mA, rci_t nA, rci_t nB, int consistent) {
  int mB = MAX(mA, nA);
  mzd_t *A = mzd_init(mA, nA);
  mzd_t *secret = mzd_init(nA, nB);
  mzd_randomize(A);
  mzd_randomize(secret);
  mzd_t *B; 
  if (consistent) {
    B = mzd_mul(NULL, A, secret, 0); 
    assert(B->nrows == A->nrows);
    assert(B->ncols == nB);
  }
  else {
    B = mzd_init(mB, nB);
    mzd_randomize(B);
  }
  // copy A & B
  mzd_t *Acopy = mzd_copy(NULL, A);
  mzd_t *Bcopy = mzd_copy(NULL, B);
  int consistency = !mzd_solve_left(A, B, 0, 1);

  if (consistent && !consistency) {
    printf("failed (solution should have been found)\n");
    return 1;
  }

  if (!consistent && !consistency) {
    printf("skipped (OK, no solution found)\n");
    return 0;
  }
  // copy B
  mzd_t *X = mzd_submatrix(NULL, B, 0, 0, A->ncols, B->ncols);
  mzd_t *B1 = mzd_mul(NULL, Acopy, X, 0);
  mzd_t *Z  = mzd_add(NULL, Bcopy, B1);

  int status = !mzd_is_zero(Z);
  if (status == 0) {
    printf("passed\n");
  } else {
    printf("FAILED (incorrect solution returned)\n");
  }
  mzd_free(Acopy);
  mzd_free(Bcopy);
  mzd_free(B1);
  mzd_free(Z);
  mzd_free(A);
  mzd_free(B);
  mzd_free(X);
  return status;
}

int main() {
  int status = 0;

  srandom(17);

  status += test_solve_left_random(1100, 1000, 1000, 1);
  status += test_solve_left_random(1000, 1000, 1000, 1);
  status += test_solve_left_random(1100, 1100, 1000, 1);
  status += test_solve_left_random(1000, 1000, 1100, 1);
  status += test_solve_left_random(1100, 1000, 1100, 1);

  status += test_solve_left_random(1100, 1000, 1000, 0);
  status += test_solve_left_random(1000, 1000, 1000, 0);
  status += test_solve_left_random(1100, 1100, 1000, 0);
  status += test_solve_left_random(1000, 1000, 1100, 0);
  status += test_solve_left_random(1100, 1000, 1100, 0);


  for (size_t i = 0; i < 100; i++) {
    size_t m = m4ri_random_word() & 511;
    size_t n = m4ri_random_word() & 1023;
    m        = m ? (m) : 1;
    n        = n ? (n) : 1;

    status += test_pluq_solve_left(m, n, 0, 0);
  }

  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
