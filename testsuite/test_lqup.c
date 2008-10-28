#include <stdlib.h>
#include "m4ri.h"


int test_lqup_full_rank (int m, int n){
  packedmatrix* U = mzd_init (m,n);
  packedmatrix* L = mzd_init (m,m);
  packedmatrix* U2 = mzd_init (m,n);
  packedmatrix* L2 = mzd_init (m,m);
  packedmatrix* A = mzd_init (m,n);
  mzd_randomize (U);
  mzd_randomize (L);

  size_t i,j;
  for (i=0; i<m; ++i){
    for (j=0; j<i;++j)
      mzd_write_bit(U,i,j, 0);
    for (j=i+1; j<m;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(U,i,i, 1);
    mzd_write_bit(L,i,i, 1);
  }
  
  mzd_mul(A, L, U, 2048);

  packedmatrix* Acopy = mzd_copy (NULL,A);

  permutation* P = mzp_init(m);
  permutation* Q = mzp_init(n);
  mzd_lqup (A, P, Q, 2048);

  for (i=0; i<m; ++i){
    for (j=0; j<i;++j)
      mzd_write_bit (L2, i, j, mzd_read_bit(A,i,j));
    for (j=i+1; j<n;++j)
      mzd_write_bit (U2, i, j, mzd_read_bit(A,i,j));
  }
  
  for (i=0; i<n; ++i){
    mzd_write_bit(L2,i,i, 1);
    mzd_write_bit(U2,i,i, 1);
  }
  mzd_addmul(Acopy,L2,U2,0);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (Acopy,i,j)){
	status = 1;
      }
    }
  if (status)
    printf(" FAILED.\n");
  else
    printf (" ... passed\n");
  mzd_free(U);
  mzd_free(L);
  mzd_free(U2);
  mzd_free(L2);
  mzd_free(A);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);
  return status;
}

int test_lqup (int m, int n){
  packedmatrix* U = mzd_init (m,n);
  packedmatrix* L = mzd_init (m,m);
  packedmatrix* U2 = mzd_init (m,n);
  packedmatrix* L2 = mzd_init (m,m);
  packedmatrix* A = mzd_init (m,n);
  mzd_randomize (U);
  mzd_randomize (L);

  size_t i,j;
  for (i=0; i<m; ++i){
    mzd_write_bit(U,i,i, 1);
    for (j=0; j<i;++j)
      mzd_write_bit(U,i,j, 0);
    if (i%2)
      for (j=i; j<n;++j)
	mzd_write_bit(U,i,j, 0);
    for (j=i+1; j<m;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  
  mzd_mul(A, L, U, 2048);

  packedmatrix* Acopy = mzd_copy (NULL,A);

  permutation* P = mzp_init(m);
  permutation* Q = mzp_init(n);
  int r;
  r=mzd_lqup (A, P, Q, 2048);

  for (i=0; i<r; ++i){
    for (j=0; j<i;++j)
      mzd_write_bit (L2, i, j, mzd_read_bit(A,i,j));
    for (j=i+1; j<n;++j)
      mzd_write_bit (U2, i, j, mzd_read_bit(A,i,j));
  }
  for (i=r; i<m; i++)
    for (j=0; j<r;++j)
      mzd_write_bit (L2, i, j, mzd_read_bit(A,i,j));
  for (i=0; i<r; ++i){
    mzd_write_bit(L2,i,i, 1);
    mzd_write_bit(U2,i,i, 1);
  }
  printf("\n L2 = \n");
  mzd_print_matrix(L2);
  printf("\n U2 = \n");
  mzd_print_matrix(U2);
  printf("P = \n");
  for (i = 0; i<A->nrows; ++i)
    printf("P[%ld] = %ld\n",i,P->values[i]);
  //  mzd_apply_p_left_trans (Acopy, P);
  //mzd_apply_p_right_trans (Acopy, Q);
  mzd_addmul(Acopy,L2,U2,0);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (Acopy,i,j)){
	status = 1;
      }
    }
  if (status)
    printf(" FAILED.\n");
  else
    printf (" ... passed\n");
  mzd_free(U);
  mzd_free(L);
  mzd_free(U2);
  mzd_free(L2);
  mzd_free(A);
  mzd_free(Acopy);
  mzp_free(P);
  mzp_free(Q);
  return status;
}

int main(int argc, char **argv) {
  int status = 0;

  printf("testing base case full rank n=37");
  status += test_lqup_full_rank(37,37);
  printf("testing base case full rank n=64");
  status += test_lqup_full_rank(64,64);
  printf("testing base case half rank n=64");
  status += test_lqup (64,64);
  printf("testing LU block recursive algorithm n=97");
  status += test_lqup_full_rank(97,97);
  printf("testing LU block recursive algorithm full rank n=128");
  status += test_lqup_full_rank(128,128);
  printf("testing LU block recursive algorithm full rank n=150");
  status += test_lqup_full_rank(150,150);
  printf("testing LU block recursive algorithm full rank n=256");
  status += test_lqup_full_rank(256,256);
  printf("testing LU block recursive algorithm full rank n=1024");
  status += test_lqup_full_rank(1024,1024);

  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
