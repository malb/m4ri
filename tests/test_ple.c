#include "testing.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>
#include <stdlib.h>

int check_ple(mzd_t *A, rci_t *r) {
  mzd_t *Acopy = mzd_copy(NULL, A);

  const rci_t m = A->nrows;
  const rci_t n = A->ncols;

  mzd_t *L = mzd_init(m, m);
  mzd_t *E = mzd_init(m, n);

  mzp_t *P = mzp_init(m);
  mzp_t *Q = mzp_init(n);
  r[0]     = mzd_ple(A, P, Q, 0);
  mzd_apply_p_right_trans_tri(A, Q);
  
  for (rci_t i = 0; i < r[0]; ++i) {
    for (rci_t j = 0; j < i; ++j) mzd_write_bit(L, i, j, mzd_read_bit(A, i, j));
    for (rci_t j = i + 1; j < n; ++j) mzd_write_bit(E, i, j, mzd_read_bit(A, i, j));
  }
  for (rci_t i = r[0]; i < m; ++i)
    for (rci_t j = 0; j < r[0]; ++j) mzd_write_bit(L, i, j, mzd_read_bit(A, i, j));
  for (rci_t i = 0; i < r[0]; ++i) {
    mzd_write_bit(L, i, i, 1);
    mzd_write_bit(E, i, i, 1);
  }

  mzd_apply_p_left(Acopy, P);
  mzd_apply_p_right_trans(Acopy, Q);

  mzd_addmul(Acopy, L, E, 0);
  int status = !mzd_is_zero(Acopy);
  mzd_free(E);
  mzd_free(L);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);

  return status;
}



int test_ple_string(rci_t m, rci_t n, const char *str) {
  printf("ple: testing string m: %5d, n: %5d", m, n);
  mzd_t *A = mzd_from_str(m, n, str);
  rci_t r    = 0;
  int status = check_ple(A, &r);
  printf(", rank: %5d ", r);
  if (status) {
    printf(" ... FAILED\n");
  } else
    printf(" ... passed\n");
  mzd_free(A);
  return status;
}

int test_ple_random(rci_t m, rci_t n) {
  printf("ple: testing random m: %5d, n: %5d", m, n);

  mzd_t *A = mzd_init(m, n);
  mzd_randomize(A);

  rci_t r    = 0;
  int status = check_ple(A, &r);
  printf(", rank: %5d ", r);

  if (status) {
    printf(" ... FAILED\n");
  } else
    printf(" ... passed\n");
  mzd_free(A);
  return status;
}


int test_ple_random_lowrank(rci_t m, rci_t n) {
  printf("ple: testing random m: %5d, n: %5d", m, n);

  rci_t r = MIN(m, n);
  mzd_t *U = mzd_init(m, r / 3);
  mzd_t *V = mzd_init(r / 3, n);
  mzd_randomize(U);
  mzd_randomize(V);
  mzd_t *A = mzd_mul(NULL, U, V, 0);

  r    = 0;
  int status = check_ple(A, &r);
  printf(", rank: %5d ", r);

  if (status) {
    printf(" ... FAILED\n");
  } else
    printf(" ... passed\n");
  mzd_free(A);
  return status;
}


int main() {
  int status = 0;

  srandom(17);

  status += test_ple_string(4, 4, "1000010000100001");
  status += test_ple_string(4, 4, "0001001001001000");
  status += test_ple_string(4, 4, "0000000000000011");
  status += test_ple_string(4, 4, "1111111111111111");
  status += test_ple_string(4, 4, "0001000100011111");
  status += test_ple_string(4, 4, "1111111101110011");
  status += test_ple_string(4, 4, "0110011110101100");

  for (int i = 0; i < 10; i++)
    status += test_ple_random(4, 4);  
  // exit(0);

  status += test_ple_random(63, 63);
  status += test_ple_random(64, 64);
  status += test_ple_random(65, 65);

  status += test_ple_random(128, 128);
  status += test_ple_random(128, 131);
  status += test_ple_random(132, 731);
  status += test_ple_random(150, 150);
  status += test_ple_random(252, 24);
  status += test_ple_random(256, 256);
  status += test_ple_random(1024, 1022);
  status += test_ple_random(1024, 1024);

  status += test_ple_random(128, 1280);
  status += test_ple_random(128, 130);
  status += test_ple_random(132, 132);
  status += test_ple_random(150, 151);
  status += test_ple_random(252, 2);
  status += test_ple_random(256, 251);
  status += test_ple_random(1024, 1025);
  status += test_ple_random(1024, 1021);

  int n = 0;
  while (n*n < __M4RI_PLE_CUTOFF)
    n += 10;
  int rounded = ((n + 1 + m4ri_radix) / m4ri_radix);
  status += test_ple_random(n + 15, m4ri_radix * n + 23);
  status += test_ple_random_lowrank(n + 15, m4ri_radix * n + 23);
  status += test_ple_random_lowrank(rounded, rounded);


  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
