#include <stdlib.h>
#include "m4ri/m4ri.h"

int test_pluq_solve_left (int m, int n, int offsetA, int offsetB){
  packedmatrix* Abase = mzd_init (2048,2048);
  packedmatrix* Bbase = mzd_init (2048,2048);
  mzd_randomize (Abase);
  mzd_randomize (Bbase);
  packedmatrix* Bbasecopy = mzd_copy (NULL, Bbase);

  packedmatrix* A = mzd_init_window (Abase, 0, offsetA, m, offsetA + m);
  packedmatrix* B = mzd_init_window (Bbase, 0, offsetB, m, offsetB + n);
  
  size_t i,j;

  packedmatrix* W = mzd_init (B->nrows, B->ncols);
  for ( i=0; i<B->nrows; ++i)
      for ( j=0; j<B->ncols; ++j)
          mzd_write_bit(W,i,j, mzd_read_bit (B,i,j));

  for (i=0; i<m; ++i){
    mzd_write_bit(A,i,i, 1);
  }  

  mzd_solve_left(A, B, 2048, 1);

/*  mzd_addmul(W, A, B, 2048); */
  packedmatrix *L = mzd_init(A->nrows,A->nrows);
  for ( i=0; i<m; ++i)
      for ( j=0; j<=i; ++j)
          if (mzd_read_bit (A,i,j)) 
              mzd_write_bit(L,i,j,1);

  packedmatrix *U = mzd_init(A->nrows,A->ncols);
  for ( i=0; i<A->nrows; ++i)
      for ( j=i; j<A->ncols; ++j)
          if (mzd_read_bit (A,i,j)) 
              mzd_write_bit(U,i,j,1);

  packedmatrix *X = mzd_init(B->nrows,B->ncols);
  for ( i=0; i<B->nrows; ++i)
      for ( j=0; j<B->ncols; ++j)
          mzd_write_bit(X,i,j, mzd_read_bit (B,i,j));
  packedmatrix *T = mzd_init(A->nrows,B->ncols);
  mzd_mul(T, U, X, 2048);

  packedmatrix *H = mzd_init(B->nrows,B->ncols);
  mzd_mul(H, L, T, 2048);

  packedmatrix *Z = mzd_init(B->nrows,B->ncols);

  mzd_add(Z, W, H);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (Z,i,j)){
	status = 1;
      }
    }

  mzd_free(L);
  mzd_free(U);
  mzd_free(T);
  mzd_free(H);
  mzd_free(Z);
  mzd_free_window (A);
  mzd_free_window (B);
  mzd_free_window (W);
  mzd_free(Abase);
  mzd_free(Bbase);
  mzd_free(Bbasecopy);

  if (!status)
    printf("passed\n");
  else
    printf("FAILED\n");
  return status;
}

int main(int argc, char **argv) {
  int status = 0;

  printf("SolveLeft: small A even, small B odd  ... ");
  status += test_pluq_solve_left (2, 4, 0, 1);

  printf("SolveLeft: small A even, small B even ... ");
  status += test_pluq_solve_left (10, 20, 0, 0);
  printf("SolveLeft: small A even, large B even ... ");
  status += test_pluq_solve_left (10, 80, 0, 0);

  printf("SolveLeft: small A even, small B odd  ... ");
  status += test_pluq_solve_left (10, 20, 0, 15);
  printf("SolveLeft: small A even, large B odd  ... ");
  status += test_pluq_solve_left (10, 80, 0, 15);

  printf("SolveLeft: small A odd, small B even  ... ");
  status += test_pluq_solve_left (10, 20, 15, 0);
  printf("SolveLeft: small A odd, large B even  ... ");
  status += test_pluq_solve_left (10, 80, 15, 0);

  printf("SolveLeft: small A odd, small B odd   ... ");
  status += test_pluq_solve_left (10, 20, 15, 20);
  printf("SolveLeft: small A odd, large B odd   ... ");
  status += test_pluq_solve_left (10, 80, 15, 20);

  printf("SolveLeft: large A even, small B even ... ");
  status += test_pluq_solve_left (70, 20, 0, 0);
  printf("SolveLeft: large A even, large B even ... ");
  status += test_pluq_solve_left (70, 80, 0, 0);

  printf("SolveLeft: large A even, small B odd  ... ");
  status += test_pluq_solve_left (70, 10, 0, 15);
  printf("SolveLeft: large A even, large B odd  ... ");
  status += test_pluq_solve_left (70, 80, 0, 15);

  printf("SolveLeft: large A odd, small B even  ... ");
  status += test_pluq_solve_left (70, 20, 15, 0);
  printf("SolveLeft: large A odd, large B even  ... ");
  status += test_pluq_solve_left (70, 80, 15, 0);

  printf("SolveLeft: large A odd, small B odd   ... ");
  status += test_pluq_solve_left (70, 20, 15, 20);
  printf("SolveLeft: large A odd, large B odd   ... ");
  status += test_pluq_solve_left (70, 80, 15, 20);

  printf("SolveLeft: larger A odd, larger B odd ... ");
  status += test_pluq_solve_left (770, 1600, 75, 89);
  printf("SolveLeft: larger A odd, larger B odd ... ");
  status += test_pluq_solve_left (1764, 1345, 198, 123);

  
  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
