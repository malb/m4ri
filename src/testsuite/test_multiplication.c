#include <stdlib.h>
#include "packedmatrix.h"
#include "brilliantrussian.h"
#include "grayflex.h"

int main(int argc, char **argv) {
  m4ri_build_all_codes();
  packedmatrix *A, *B, *C, *D, *E;
  int n;
  int eq;

  A = mzd_init(21,210);
  B = mzd_init(210,21);
  mzd_randomize(A);
  mzd_randomize(B);
  C = mzd_mul_strassen(NULL, A,B, 64);
  D = mzd_mul_m4rm(NULL, A,B, 0, NULL, NULL);
  E = mzd_mul_naiv(NULL, A,B);

  eq = mzd_equal(C,D);
  eq += mzd_equal(D,E);
  if (eq==2) {
    printf("all tests passed.\n");
  } else {
    printf("Test failed with eq=%d.\n",eq);
  }

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  mzd_free(D);
  mzd_free(E);

  m4ri_destroy_all_codes();
  return 0;
}
