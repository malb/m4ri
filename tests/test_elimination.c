#include "testing.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>
#include <stdlib.h>

int elim_test_equality(rci_t nr, rci_t nc) {
  int ret = 0;

  printf("elim: m: %4d, n: %4d ", nr, nc);

  mzd_t *A = mzd_init(nr, nc);
  mzd_randomize(A);
  mzd_t *B = mzd_copy(NULL, A);
  mzd_t *C = mzd_copy(NULL, A);
  mzd_t *D = mzd_copy(NULL, A);
  mzd_t *E = mzd_copy(NULL, A);
  mzd_t *F = mzd_copy(NULL, A);
  mzd_t *G = mzd_copy(NULL, A);

  /* M4RI k=auto */
  rci_t ra = mzd_echelonize_m4ri(A, 1, 0);

  /* M4RI k=8 */
  rci_t rb = mzd_echelonize_m4ri(B, 1, 8);

  /* M4RI Upper Triangular k=auto*/
  rci_t rc = mzd_echelonize_m4ri(C, 0, 0);
  mzd_top_echelonize_m4ri(C, 0);

  /* M4RI Upper Triangular k=4*/
  rci_t rd = mzd_echelonize_m4ri(D, 0, 4);
  mzd_top_echelonize_m4ri(D, 4);

  /* Gauss */
  rci_t re = mzd_echelonize_naive(E, 1);

  /* Gauss Upper Triangular */
  rci_t rf = mzd_echelonize_naive(F, 0);
  mzd_top_echelonize_m4ri(F, 0);

  /* PLUQ */
  rci_t rg = mzd_echelonize_pluq(G, 1);

  if (mzd_equal(A, B) != TRUE || ra != rb) {
    printf("A != B ");
    ret -= 1;
  }

  if (mzd_equal(B, C) != TRUE || rb != rc) {
    printf("B != C ");
    ret -= 1;
  }

  if (mzd_equal(C, D) != TRUE || rc != rd) {
    printf("C != D ");
    ret -= 1;
  }

  if (mzd_equal(D, E) != TRUE || rd != re) {
    printf("D != E ");
    ret -= 1;
  }

  if (mzd_equal(E, F) != TRUE || re != rf) {
    printf("E != F ");
    ret -= 1;
  }

  if (mzd_equal(F, G) != TRUE || rf != rg) {
    printf("F != G ");
    ret -= 1;
  }
  if (mzd_equal(G, A) != TRUE || rg != ra) {
    printf("G != A ");
    ret -= 1;
  }

  mzd_free(A);
  mzd_free(B);
  mzd_free(C);
  mzd_free(D);
  mzd_free(E);
  mzd_free(F);
  mzd_free(G);

  if (ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }
  return ret;
}

int main() {
  int status = 0;

  srandom(17);

  status += elim_test_equality(4, 67);
  status += elim_test_equality(17, 121);
  status += elim_test_equality(65, 17);
  status += elim_test_equality(128, 128);
  status += elim_test_equality(1024, 1024);
  status += elim_test_equality(2047, 2047);
  status += elim_test_equality(65, 65);

  status += elim_test_equality(100, 100);
  status += elim_test_equality(21, 171);
  status += elim_test_equality(31, 121);
  status += elim_test_equality(193, 65);
  status += elim_test_equality(1025, 1025);
  status += elim_test_equality(2048, 2048);
  status += elim_test_equality(64, 64);
  status += elim_test_equality(128, 128);
  status += elim_test_equality(4096, 3528);
  status += elim_test_equality(1024, 1025);
  status += elim_test_equality(1000, 1000);
  status += elim_test_equality(1000, 10);
  status += elim_test_equality(1710, 1290);
  status += elim_test_equality(1290, 1710);
  status += elim_test_equality(1290, 1710);
  status += elim_test_equality(1290, 1290);
  status += elim_test_equality(1000, 210);

  if (status == 0) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
