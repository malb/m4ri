#include <stdlib.h>
#include "packedmatrix.h"
#include "brilliantrussian.h"
#include "grayflex.h"

int main(int argc, char **argv) {
  m2t_build_all_codes();
  packedmatrix *A, *B, *C;
  int n;

  if (argc > 1) {
    n = atoi(argv[1]);
  } else {
    n = 10;
  }

  A = m2t_init(1<<n,1<<n);
  B = m2t_init(1<<n,1<<n);
  C = m2t_mul_strassen(NULL, A,B, 1<<11);

  m2t_free(A);
  m2t_free(B);
  m2t_free(C);

  m2t_destroy_all_codes();
  return 0;
}
