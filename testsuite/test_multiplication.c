#include <stdlib.h>
#include "m4ri.h"

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
int mul_test_equality(int m, int l, int n, int k, int cutoff) {
  int ret  = 0;
  packedmatrix *A, *B, *C, *D, *E;
  
  printf("   mul: m: %4d, l: %4d, n: %4d, k: %2d, cutoff: %4d\n",m,l,n,k,cutoff);

  /* we create two random matrices */
  A = mzd_init(m, l);
  B = mzd_init(l, n);
  mzd_randomize(A);
  mzd_randomize(B);

  /* C = A*B via Strassen */
  C = mzd_mul(NULL, A, B, cutoff);

  /* D = A*B via M4RM, temporary buffers are managed internally */
  D = mzd_mul_m4rm(    NULL, A, B, k);

  /* E = A*B via naiv cubic multiplication */
  E = mzd_mul_naiv(    NULL, A, B);

  mzd_free(A);
  mzd_free(B);

  if (mzd_equal(C, D) != TRUE) {
    printf("FAIL: Strassen != M4RM\n");
    ret -=1;
  }

  if (mzd_equal(D, E) != TRUE) {
    printf("FAIL: M4RM != Naive\n");
    ret -= 1;
  }

  if (mzd_equal(C, E) != TRUE) {
    printf("FAIL: Strassen != Naive\n");
    ret -= 1;
  }

  mzd_free(C);
  mzd_free(D);
  mzd_free(E);
  return ret;

}

int addmul_test_equality(int m, int l, int n, int k, int cutoff) {
  int ret  = 0;
  packedmatrix *A, *B, *C, *D, *E, *F;
  
  printf("addmul: m: %4d, l: %4d, n: %4d, k: %2d, cutoff: %4d\n",m,l,n,k,cutoff);

  /* we create two random matrices */
  A = mzd_init(m, l);
  B = mzd_init(l, n);
  C = mzd_init(m, n);
  mzd_randomize(A);
  mzd_randomize(B);
  mzd_randomize(C);

  /* D = C + A*B via M4RM, temporary buffers are managed internally */
  D = mzd_copy(NULL, C);
  D = mzd_addmul_m4rm(D, A, B, k);

  /* E = C + A*B via naiv cubic multiplication */
  E = mzd_mul_m4rm(NULL, A, B, k);
  mzd_add(E, E, C);

  /* F = C + A*B via naiv cubic multiplication */
  F = mzd_copy(NULL, C);
  F = mzd_addmul(F, A, B, cutoff);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);

  if (mzd_equal(D, E) != TRUE) {
    printf("FAIL: addmul_m4rm != add,mul\n");
    ret -=1;
  }
  if (mzd_equal(E, F) != TRUE) {
    printf("FAIL: add,mul = addmul\n");
    ret -=1;
  }
  if (mzd_equal(F, D) != TRUE) {
    printf("FAILL addmul_m4rm != addmul\n");
    ret -=1;
  }

  mzd_free(D);
  mzd_free(E);
  mzd_free(F);
  return ret;
}

int main(int argc, char **argv) {
  int status = 0;
  m4ri_build_all_codes();
  
  status += mul_test_equality(1, 1, 1, 0, 1024);
  status += mul_test_equality(64, 64, 64, 0, 64);
  status += mul_test_equality(128, 128, 128, 0, 64);
  status += mul_test_equality(21, 171, 31, 0, 63); 
  status += mul_test_equality(21, 171, 31, 0, 131); 
  status += mul_test_equality(193, 65, 65, 10, 64);
  status += mul_test_equality(1025, 1025, 1025, 3, 256);
  status += mul_test_equality(2048, 2048, 4096, 0, 1024);
  status += mul_test_equality(4096, 3528, 4096, 0, 1024);
  status += mul_test_equality(1024, 1025, 1, 0, 1024);
  status += mul_test_equality(1000,1000,1000, 0, 256);
  status += mul_test_equality(1000,10,20, 0, 64);
  status += mul_test_equality(1710,1290,1000, 0, 256);
  status += mul_test_equality(1290, 1710, 200, 0, 64);
  status += mul_test_equality(1290, 1710, 2000, 0, 256);
  status += mul_test_equality(1290, 1290, 2000, 0, 64);
  status += mul_test_equality(1000, 210, 200, 0, 64);

  status += addmul_test_equality(64, 64, 64, 0, 64);
  status += addmul_test_equality(128, 128, 128, 0, 64);
  status += addmul_test_equality(21, 171, 31, 0, 63);
  status += addmul_test_equality(21, 171, 31, 0, 131);
  status += addmul_test_equality(193, 65, 65, 10, 64);
  status += addmul_test_equality(1025, 1025, 1025, 3, 256);
  status += addmul_test_equality(4096, 4096, 4096, 0, 2048);
  status += addmul_test_equality(1000,1000,1000, 0, 256);
  status += addmul_test_equality(1000,10,20, 0, 64);
  status += addmul_test_equality(1710,1290,1000, 0, 256);
  status += addmul_test_equality(1290, 1710, 200, 0, 64);
  status += addmul_test_equality(1290, 1710, 2000, 0, 256);
  status += addmul_test_equality(1290, 1290, 2000, 0, 64);
  status += addmul_test_equality(1000, 210, 200, 0, 64);

  if (status == 0) {
    printf("All tests passed.\n");
  }

  m4ri_destroy_all_codes();
  return 0;
}
