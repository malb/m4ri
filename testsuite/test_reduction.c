#include <stdlib.h>
#include "m4ri.h"

int red_test_equality(int nr, int nc) {
  packedmatrix *A, *B, *C, *D, *E, *F;
  int ret = 0; 

  printf("red: m: %4d, n: %4d\n",nr,nc);

  A = mzd_init(nr, nc);
  mzd_randomize(A);
  B = mzd_copy(NULL, A);
  C = mzd_copy(NULL, A);
  D = mzd_copy(NULL, A);
  E = mzd_copy(NULL, A);
  F = mzd_copy(NULL, A);

  mzd_reduce_m4ri(A, 1, 0, NULL, NULL);

  mzd_reduce_m4ri(B, 1, 8, NULL, NULL);

  mzd_reduce_m4ri(C, 0, 0, NULL, NULL);
  mzd_top_reduce_m4ri(C, 0, NULL, NULL);

  mzd_reduce_m4ri(D, 0, 4, NULL, NULL);
  mzd_top_reduce_m4ri(D, 4, NULL, NULL);

  mzd_reduce_naiv(E, 1);

  mzd_reduce_naiv(F, 0);
  mzd_top_reduce_m4ri(F, 0, NULL, NULL);
  
  if(mzd_equal(A, B) != TRUE) {
    printf("FAIL: A != B\n");
    ret -= 1;
  }
 
  if(mzd_equal(B, C) != TRUE) {
    printf("FAIL: B != C\n");
    ret -= 1;
  }

  if(mzd_equal(D, E) != TRUE) {
    printf("FAIL: D != E\n");
    ret -= 1;
  }

  if(mzd_equal(E, F) != TRUE) {
    printf("FAIL: E != F\n");
    ret -= 1;
  }

  if(mzd_equal(F, A) != TRUE) {
    printf("FAIL: F != A\n");
    ret -= 1;
  }

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  mzd_free(D);
  mzd_free(E);

  return ret;
}

int main(int argc, char **argv) {
  m4ri_build_all_codes();

  int status = 0;
  status += red_test_equality(100, 100);
  status += red_test_equality(21, 171);
  status += red_test_equality(31, 121);
  status += red_test_equality(193, 65);
  status += red_test_equality(1025, 1025);
  status += red_test_equality(2048, 2048);

  if (status == 0) {
    printf("All tests passed.\n");
  }

  m4ri_destroy_all_codes();
  return 0;
}
