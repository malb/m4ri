#include <stdlib.h>
#include "packedmatrix.h"
#include "brilliantrussian.h"
#include "grayflex.h"

int test_equality(nr, nc) {
  packedmatrix *A, *B, *C, *D, *E;
  int ret = 0; 
  A = createMatrix(100,100);
  fillRandomly(A);
  B = clone(A);
  C = clone(A);
  D = clone(A);
  E = clone(A);

  reduceM4RI(A, 1, 0, NULL, NULL);

  reduceM4RI(B, 1, 8, NULL, NULL);

  reduceM4RI(C, 0, 0, NULL, NULL);
  topReduceM4RI(C, 0, NULL, NULL);

  reduceM4RI(D, 0, 4, NULL, NULL);
  topReduceM4RI(D, 4, NULL, NULL);

  reduceGaussian(E, 1);
  
  ret = equalMatrix(A, B);
  ret += equalMatrix(A, C);
  ret += equalMatrix(A, D);
  ret += equalMatrix(A, E);

  destroyMatrix(A);
  destroyMatrix(B);
  destroyMatrix(C);
  destroyMatrix(D);
  destroyMatrix(E);

  return ret - 4;
}

int main(int argc, char **argv) {
  buildAllCodes();

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

  destroyAllCodes();
  return 0;
}
