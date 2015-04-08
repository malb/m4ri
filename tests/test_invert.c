#include <m4ri/config.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>

/**
 * Check that inversion works.
 *
 * \param n Number of rows of A
 * \param k Parameter k of M4RM algorithm, may be 0 for automatic choice.
 */
int invert_test(rci_t n, int k) {
  int ret  = 0;
  printf("invert: n: %4d, k: %2d", n, k);

  mzd_t *I2 = mzd_init(n, n);
  mzd_set_ui(I2, 1);

  mzd_t *U = mzd_init(n,n);
  mzd_randomize(U);
  for(rci_t i=0; i<n; i++) {
    mzd_write_bit(U, i, i, 1);
    for (rci_t j=0; j<i; j++)
      mzd_write_bit(U, i, j, 0);
  }

  mzd_t *B = mzd_copy(NULL, U);
  mzd_trtri_upper(B);

  mzd_t *I1 = mzd_mul(NULL, U, B, 0);

  if (mzd_equal(I1, I2) != TRUE) {
    ret += 1;
    printf(" U*~U != 1 ");
  }

  mzd_t *L = mzd_init(n, n);
  mzd_randomize(L);
  for (rci_t i = 0; i < n; ++i) {
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_t *A = mzd_mul(NULL, L, U, 0);

  B = mzd_inv_m4ri(B, A, k);

  I1 = mzd_mul(I1, A, B, 0);

  if (mzd_equal(I1, I2) != TRUE) {
    ret += 1;
    printf(" A*~A != 1 ");
  }

  if(ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }
  mzd_free(U);
  mzd_free(L);
  mzd_free(A);
  mzd_free(B);
  mzd_free(I1);
  mzd_free(I2);

  return ret;

}

int main() {
  int status = 0;
  srandom(17);

  for(int k=0; k<5; k++) {
    status += invert_test(   1,k);
    status += invert_test(   2,k);
    status += invert_test(   3,k);
    status += invert_test(  21,k);
    status += invert_test(  64,k);
    status += invert_test( 128,k);
    status += invert_test( 193,k);
    status += invert_test(1000,k);
    status += invert_test(1024,k);
    status += invert_test(1025,k);
    status += invert_test(1290,k);
    status += invert_test(1710,k);
    status += invert_test(2048,k);
  }
  if (status == 0) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
