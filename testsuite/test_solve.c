#include <stdlib.h>
#include "m4ri/m4ri.h"

int test_pluq_solve_left(size_t m, size_t n, size_t offsetA, size_t offsetB){
  packedmatrix* Abase = mzd_init(2048, 2048);
  packedmatrix* Bbase = mzd_init(2048, 2048);
  mzd_randomize(Abase);
  mzd_randomize(Bbase);

  packedmatrix* A = mzd_init_window(Abase, 0, offsetA, m, offsetA + m);
  packedmatrix* B = mzd_init_window(Bbase, 0, offsetB, m, offsetB + n);
  
  size_t i,j;

  // copy B
  packedmatrix* Bcopy = mzd_init(B->nrows, B->ncols);
  for ( i=0; i<B->nrows; ++i)
      for ( j=0; j<B->ncols; ++j)
          mzd_write_bit(Bcopy,i,j, mzd_read_bit (B,i,j));

  for (i=0; i<m; ++i) {
    mzd_write_bit(A,i,i, 1);
  }

  packedmatrix *Acopy = mzd_copy(NULL, A);
  size_t r = mzd_echelonize_m4ri(Acopy,0,0,NULL,NULL);
  printf("solve_left m: %4zu, n: %4zu, r: %4zu da: %2zu db: %2zu ",m, n, r, offsetA, offsetB);
  mzd_free(Acopy);
  Acopy = mzd_copy(NULL, A);
    
  mzd_solve_left(A, B, 0, 1);

  //copy B
  packedmatrix *X = mzd_init(B->nrows,B->ncols);
  for ( i=0; i<B->nrows; ++i)
      for ( j=0; j<B->ncols; ++j)
          mzd_write_bit(X,i,j, mzd_read_bit (B,i,j));

  packedmatrix *B1 = mzd_mul(NULL, Acopy, X, 0);
  packedmatrix *Z = mzd_add(NULL, Bcopy, B1);
  
  int status = 0;
  
  if(r==m) {
    for ( i=0; i<m; ++i)
      for ( j=0; j<n; ++j){
        if (mzd_read_bit (Z,i,j)){
          status = 1;
        }
      }
    if (!status)
      printf("passed\n");
    else
      printf("FAILED\n");

  } else {
    printf("check skipped r!=m\n");
  }
  mzd_free(Bcopy);
  mzd_free(B1);
  mzd_free(Z);

  mzd_free_window(A);
  mzd_free_window(B);
  mzd_free(Acopy);
  mzd_free(Abase);
  mzd_free(Bbase);

  return status;
}

int main(int argc, char **argv) {
  int status = 0;

  status += test_pluq_solve_left(  2,   4,  0,  0);
  status += test_pluq_solve_left(  4,   1,  0,  0);
  status += test_pluq_solve_left( 10,  20,  0,  0);
  status += test_pluq_solve_left( 20,   1,  0,  0);
  status += test_pluq_solve_left( 20,  20,  0,  0);
  status += test_pluq_solve_left( 30,   1,  0,  0);
  status += test_pluq_solve_left( 30,  30,  0,  0);
  status += test_pluq_solve_left( 80,   1,  0,  0);
  status += test_pluq_solve_left( 80,  20,  0,  0);
  status += test_pluq_solve_left( 80,  80,  0,  0);

  status += test_pluq_solve_left(  2,   4,  0,  1);
  status += test_pluq_solve_left(  4,   1,  0,  1);
  status += test_pluq_solve_left( 10,  20,  0,  1);
  status += test_pluq_solve_left( 20,   1,  0,  1);
  status += test_pluq_solve_left( 20,  20,  0,  1);
  status += test_pluq_solve_left( 30,   1,  0,  1);
  status += test_pluq_solve_left( 30,  30,  0,  1);
  status += test_pluq_solve_left( 80,   1,  0,  1);
  status += test_pluq_solve_left( 80,  20,  0,  1);
  status += test_pluq_solve_left( 80,  80,  0,  1);

  status += test_pluq_solve_left(  2,   4,  0, 15);
  status += test_pluq_solve_left(  4,   1,  0, 15);
  status += test_pluq_solve_left( 10,  20,  0, 15);
  status += test_pluq_solve_left( 20,   1,  0, 15);
  status += test_pluq_solve_left( 20,  20,  0, 15);
  status += test_pluq_solve_left( 30,   1,  0, 15);
  status += test_pluq_solve_left( 30,  30,  0, 15);
  status += test_pluq_solve_left( 80,   1,  0, 15);
  status += test_pluq_solve_left( 80,  20,  0, 15);
  status += test_pluq_solve_left( 80,  80,  0, 15);


/*   status += test_pluq_solve_left (10, 20, 15, 0); */
/*   status += test_pluq_solve_left (10, 80, 15, 0); */
/*   status += test_pluq_solve_left (10, 20, 15, 20); */
/*   status += test_pluq_solve_left (10, 80, 15, 20); */
/*   status += test_pluq_solve_left (70, 20, 15, 0); */
/*   status += test_pluq_solve_left (70, 80, 15, 0); */
/*   status += test_pluq_solve_left (70, 20, 15, 20); */
/*   status += test_pluq_solve_left (70, 80, 15, 20); */
/*   status += test_pluq_solve_left (770, 1600, 75, 89); */
/*   status += test_pluq_solve_left (1764, 1345, 198, 123); */

  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
