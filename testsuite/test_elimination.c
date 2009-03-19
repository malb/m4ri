#include <stdlib.h>
#include "m4ri/m4ri.h"

int elim_test_equality(int nr, int nc) {
  int ret = 0; 

  printf("elim: m: %4d, n: %4d ",nr,nc);

  mzd_t *A = mzd_init(nr, nc);
  mzd_randomize(A);
  mzd_t *B = mzd_copy(NULL, A);
  mzd_t *C = mzd_copy(NULL, A);
  mzd_t *D = mzd_copy(NULL, A);
  mzd_t *E = mzd_copy(NULL, A);
  mzd_t *F = mzd_copy(NULL, A);
  mzd_t *G = mzd_copy(NULL, A);

  /* M4RI k=auto */
  mzd_echelonize_m4ri(A, 1, 0);

  /* M4RI k=8 */
  mzd_echelonize_m4ri(B, 1, 8);

  /* M4RI Upper Triangular k=auto*/
  mzd_echelonize_m4ri(C, 0, 0);
  mzd_top_echelonize_m4ri(C, 0);

  /* M4RI Upper Triangular k=4*/
  mzd_echelonize_m4ri(D, 0, 4);
  mzd_top_echelonize_m4ri(D, 4);

  /* Gauss */
  mzd_echelonize_naive(E, 1);

  /* Gauss Upper Triangular */
  mzd_echelonize_naive(F, 0);
  mzd_top_echelonize_m4ri(F, 0);

  /* PLUQ */
  mzd_echelonize_pluq(G, 1);
  
  if(mzd_equal(A, B) != TRUE) {
    printf("A != B ");
    ret -= 1;
  }
 
  if(mzd_equal(B, C) != TRUE) {
    printf("B != C ");
    ret -= 1;
  }

  if(mzd_equal(D, E) != TRUE) {
    printf("D != E ");
    ret -= 1;
  }

  if(mzd_equal(E, F) != TRUE) {
    printf("E != F ");
    ret -= 1;
  }

  if(mzd_equal(F, G) != TRUE) {
    printf("F != G ");
    ret -= 1;
  }
  if(mzd_equal(G, A) != TRUE) {
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
  status += elim_test_equality(1000,1000);
  status += elim_test_equality(1000,10);
  status += elim_test_equality(1710,1290);
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
