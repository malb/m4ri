#include <stdlib.h>
#include "packedmatrix.h"
#include "brilliantrussian.h"
#include "grayflex.h"

int test_equality(nr, nc) {
  packedmatrix *A, *B, *C, *D, *E;
  int ret = 0; 
  A = mzd_init(nr, nc);
  mzd_randomize(A);
  B = mzd_copy(NULL, A);
  C = mzd_copy(NULL, A);
  D = mzd_copy(NULL, A);
  E = mzd_copy(NULL, A);

  mzd_reduce_m4ri(A, 1, 0, NULL, NULL);

  mzd_reduce_m4ri(B, 1, 8, NULL, NULL);

  mzd_reduce_m4ri(C, 0, 0, NULL, NULL);
  mzd_top_reduce_m4ri(C, 0, NULL, NULL);

  mzd_reduce_m4ri(D, 0, 4, NULL, NULL);
  mzd_top_reduce_m4ri(D, 4, NULL, NULL);

  mzd_reduce_naiv(E, 1);
  
  ret = mzd_equal(A, B);
  ret += mzd_equal(A, C);
  ret += mzd_equal(A, D);
  ret += mzd_equal(A, E);

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  mzd_free(D);
  mzd_free(E);

  return ret - 4;
}

int main(int argc, char **argv) {
  m4ri_build_all_codes();

  int eq, failed = 0;
  eq = test_equality(100, 100);
  if (eq != 0) {
    printf("%d, %d, failed\n", 100, 100);
    failed = 1;
  }

  eq = test_equality(100, 120);
  if (eq != 0) {
    printf("%d, %d, failed\n", 100, 120);
    failed = 1;
  }

  eq = test_equality(120, 100);
  if (eq != 0) {
    printf("%d, %d, failed\n", 120, 100);
    failed = 1;
  }

  if (failed == 0) {
    printf("all tests passed.\n");
    failed = 1;
  }

  m4ri_destroy_all_codes();
  return 0;
}
