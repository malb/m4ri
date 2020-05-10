#include "m4ri/m4ri.h"
#include <stdio.h>

int main(int argc, char const *argv[]) {
  int rows    = 16;
  int columns = 229;
  mzd_t *A    = mzd_init(rows, columns);

  mzd_t *B = mzd_init(rows, columns);
  mzd_t *C;

  mzd_randomize(A);

  mzd_copy(B, A);
  mzd_echelonize(B, 1);

  mzd_copy(B, A);
  C = mzd_init_window(B, 0, 128, 16, 229);
  mzd_echelonize(C, 1);
  mzd_free(C);

  mzd_copy(B, A);
  C = mzd_init_window(B, 0, 64, 16, 229);
  mzd_echelonize(C, 1);
  mzd_free(C);

  mzd_free(B);
  mzd_free(A);

  return 0;
}
