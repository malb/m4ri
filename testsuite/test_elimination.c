#include "config.h"
#include <stdlib.h>
#include "m4ri.h"

int elim_test_equality(rci_t nr, rci_t nc) {
  int ret = 0; 

  printf("elim: m: %4d, n: %4d ", nr.val(), nc.val());

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

  if(mzd_equal(A, B) != TRUE || ra != rb) {
    printf("A != B ");
    ret -= 1;
  }
 
  if(mzd_equal(B, C) != TRUE || rb != rc) {
    printf("B != C ");
    ret -= 1;
  }

  if(mzd_equal(C, D) != TRUE || rc != rd) {
    printf("C != D ");
    ret -= 1;
  }

  if(mzd_equal(D, E) != TRUE || rd != re) {
    printf("D != E ");
    ret -= 1;
  }

  if(mzd_equal(E, F) != TRUE || re != rf) {
    printf("E != F ");
    ret -= 1;
  }

  if(mzd_equal(F, G) != TRUE || rf != rg) {
    printf("F != G ");
    ret -= 1;
  }
  if(mzd_equal(G, A) != TRUE || rg != ra) {
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

  if(ret == 0) {
    printf(" ... passed\n");
  } else {
    printf(" ... FAILED\n");
  }
  return ret;
}

int main(int argc, char **argv) {
  int status = 0;

  status += elim_test_equality(4U, 67U);
  status += elim_test_equality(17U, 121U);
  status += elim_test_equality(65U, 17U);
  status += elim_test_equality(128U, 128U);
  status += elim_test_equality(1024U, 1024U);
  status += elim_test_equality(2047U, 2047U);
  status += elim_test_equality(65U, 65U);

  status += elim_test_equality(100U, 100U);
  status += elim_test_equality(21U, 171U);
  status += elim_test_equality(31U, 121U);
  status += elim_test_equality(193U, 65U);
  status += elim_test_equality(1025U, 1025U);
  status += elim_test_equality(2048U, 2048U);
  status += elim_test_equality(64U, 64U);
  status += elim_test_equality(128U, 128U);
  status += elim_test_equality(4096U, 3528U);
  status += elim_test_equality(1024U, 1025U);
  status += elim_test_equality(1000U,1000U);
  status += elim_test_equality(1000U,10U);
  status += elim_test_equality(1710U,1290U);
  status += elim_test_equality(1290U, 1710U);
  status += elim_test_equality(1290U, 1710U);
  status += elim_test_equality(1290U, 1290U);
  status += elim_test_equality(1000U, 210U);

  if (status == 0) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
