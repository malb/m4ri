#include <m4ri/config.h>
#include <stdlib.h>
#include <inttypes.h>
#include "testing.h"
#include <m4ri/m4ri.h>

//#define ABORT_ON_FAIL 1

int smallops_test_add(rci_t M, rci_t N, rci_t m, rci_t n, word pattern) {
  int ret = 0; 

  printf("      mzd_add: M: %4d, N: %4d, m: %4d, n: %4d, pattern: 0x%" PRIx64 " ", M, N, m, n, pattern);

  mzd_t *AA;
  mzd_t *A = mzd_init_test_matrix_random(M, N, m, n, pattern, &AA);

  mzd_t *BB;
  mzd_t *B = mzd_init_test_matrix_random(M, N, m, n, pattern, &BB);

  mzd_t *CC;
  mzd_t *C = mzd_init_test_matrix_random(M, N, m, n, pattern, &CC);

  mzd_t *DD;
  mzd_t *D = mzd_init_test_matrix_random(M, N, m, n, pattern, &DD);

  /* Creation went okay? */

  ret += mzd_check_pattern(AA, m, n, pattern);
  ret += mzd_check_pattern(BB, m, n, pattern);
  ret += mzd_check_pattern(CC, m, n, pattern);
  ret += mzd_check_pattern(DD, m, n, pattern);


  /* Testing equality A+A == 0 */

  mzd_add(C, A, A);

  if(!mzd_is_zero(C)) {
    ret +=1;
  }

  ret += mzd_check_pattern(AA, m, n, pattern);
  ret += mzd_check_pattern(BB, m, n, pattern);
  ret += mzd_check_pattern(CC, m, n, pattern);

  /* Testing equality A+A == 0 but this time C is already zero */

  mzd_add(C, B, B);

  if(!mzd_is_zero(C)) {
    ret +=1;
  }

  ret += mzd_check_pattern(AA, m, n, pattern);
  ret += mzd_check_pattern(BB, m, n, pattern);
  ret += mzd_check_pattern(CC, m, n, pattern);

  /* Testing in place add. C is zero, so afterwards C == A */

  mzd_add(C, C, A);

  if(!mzd_equal(C,A)) {
    ret +=1;
  }

  ret += mzd_check_pattern(AA, m, n, pattern);
  ret += mzd_check_pattern(BB, m, n, pattern);
  ret += mzd_check_pattern(CC, m, n, pattern);

  /* Testing equality C (== A) + A == 0 */

  mzd_add(B, C, A);

  if(!mzd_is_zero(B)) {
    ret +=1;
  }


  if(m == n) {
    /* Testing equality (A + B)^2 == A^2 + BA + AB + B^2 */

    mzd_randomize(A);
    mzd_randomize(B);
    
    mzd_add(C,A,B);
    
    mzd_mul(D,C,C, 0); // (A+B)^2
    
    mzd_mul(C,A,A, 0); 
    mzd_addmul(C, B, A, 0);
    mzd_addmul(C, A, B, 0);
    mzd_addmul(C, B, B, 0);
    
    if(!mzd_equal(C,D)) {
      ret += 1;
    }

    ret += mzd_check_pattern(AA, m, n, pattern);
    ret += mzd_check_pattern(BB, m, n, pattern);
    ret += mzd_check_pattern(CC, m, n, pattern);
    ret += mzd_check_pattern(DD, m, n, pattern);
  }

  mzd_free_test_matrix_random(AA, A);
  mzd_free_test_matrix_random(BB, B);
  mzd_free_test_matrix_random(CC, C);
  mzd_free_test_matrix_random(DD, D);

  if(ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }
#ifdef ABORT_ON_FAIL
  if (ret) abort();
#endif

  return ret;
}


int main() {
  int status = 0;

  srandom(17);

  status += smallops_test_add( 64,  64,  10,  10, 0x03030303030303llu);
  status += smallops_test_add(100, 100,  64,  64, 0x03030303030303llu);

  status += smallops_test_add(1024, 1024, 513,    511, 0x03030303030303llu);
  status += smallops_test_add(1024, 1024, 512, 768+30, 0x03030303030303llu);

  status += smallops_test_add(2048, 2048, 1024, 1024,  0x03030303030303llu);

  if (status == 0) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
