#include <stdlib.h>
#include "m4ri.h"

int main(int argc, char **argv) {
  m4ri_build_all_codes();

  packedmatrix *A, *B, *C, *D, *E;
  int eq;

  /* we create two random matrices */
  A = mzd_init(21, 210);
  B = mzd_init(210, 17);
  mzd_randomize(A);
  mzd_randomize(B);

  /* C = A*B via Strassen with cutoff = 64 */
  C = mzd_mul_strassen(NULL, A, B, 64);

  /* D = A*B via M4RM, k is automatically chosen, temporary buffers
     are managed internally */
  D = mzd_mul_m4rm(    NULL, A, B, 0, NULL, NULL);

  /* E = A*B via naiv cubic multiplication */
  E = mzd_mul_naiv(    NULL, A, B);

  mzd_free(A);
  mzd_free(B);

  if (mzd_equal(C, D) != TRUE) {
    mzd_free(C);
    mzd_free(D);
    mzd_free(E);
    m4ri_destroy_all_codes();
    m4ri_die("Strassen result differs from M4RM result.");
  }

  if (mzd_equal(D, E) != TRUE) {
    mzd_free(C);
    mzd_free(D);
    mzd_free(E);
    m4ri_destroy_all_codes();
    m4ri_die("M4RM result differs from naiv result.");
  }

  printf("All tests passed.\n");

  mzd_free(C);
  mzd_free(D);
  mzd_free(E);

  m4ri_destroy_all_codes();
  return 0;
}
