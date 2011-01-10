#define _GNU_SOURCE
#include <stdlib.h>
#include "m4ri.h"

int test_pluq_solve_left(size_t m, size_t n, size_t offsetA, size_t offsetB){
  mzd_t* Abase = mzd_init(2048, 2048);
  mzd_t* Bbase = mzd_init(2048, 2048);
  mzd_randomize(Abase);
  mzd_randomize(Bbase);

  mzd_t* A = mzd_init_window(Abase, 0, offsetA, m, offsetA + m);
  mzd_t* B = mzd_init_window(Bbase, 0, offsetB, m, offsetB + n);
  
  size_t i,j;

  // copy B
  mzd_t* Bcopy = mzd_init(B->nrows, B->ncols);
  for ( i=0; i<B->nrows; ++i)
      for ( j=0; j<B->ncols; ++j)
          mzd_write_bit(Bcopy,i,j, mzd_read_bit (B,i,j));

  for (i=0; i<m; ++i) {
    mzd_write_bit(A,i,i, 1);
  }

  mzd_t *Acopy = mzd_copy(NULL, A);
  size_t r = mzd_echelonize(Acopy,1);
  printf("solve_left m: %4zu, n: %4zu, r: %4zu da: %4zu db: %4zu ",m, n, r, offsetA, offsetB);
  mzd_free(Acopy);
  Acopy = mzd_copy(NULL, A);
    
  int consistency = mzd_solve_left(A, B, 0, 1);

  //copy B
  mzd_t *X = mzd_init(B->nrows,B->ncols);
  for ( i=0; i<B->nrows; ++i)
      for ( j=0; j<B->ncols; ++j)
          mzd_write_bit(X,i,j, mzd_read_bit (B,i,j));

  mzd_t *B1 = mzd_mul(NULL, Acopy, X, 0);
  mzd_t *Z = mzd_add(NULL, Bcopy, B1);
  
  int status = 0;
  
  if(consistency == 0) {
    status = 1 - mzd_is_zero(Z);
    if (status == 0) {
      printf("passed\n");
    }  else {
      printf("FAILED\n");
    }
  } else {
    printf("skipped (no solution)\n");
  }
  mzd_free(Bcopy);
  mzd_free(B1);
  mzd_free(Z);

  mzd_free_window(A);
  mzd_free_window(B);
  mzd_free(Acopy);
  mzd_free(Abase);
  mzd_free(Bbase);
  mzd_free(X);
  return status;
}

int main(int argc, char **argv) {
  int status = 0;

  for(size_t i=0; i<100; i++) {
    size_t m = random() & 511;
    size_t n = random() & 1023;
    m = m ? (m) : 1;
    n = n ? (n) : 1;

    status += test_pluq_solve_left(  m,   n,  0,  0);
  }

  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
