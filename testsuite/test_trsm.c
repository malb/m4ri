#include <stdlib.h>
#include "m4ri.h"


int testTRSMUpperRight (int m, int n, int offset){
  packedmatrix* Ubase = mzd_init (2048,2048);
  packedmatrix* Bbase = mzd_init (2048,2048);
  mzd_randomize (Ubase);
  mzd_randomize (Bbase);
  packedmatrix* Bbasecopy = mzd_copy (NULL, Bbase);

  packedmatrix* U = mzd_init_window (Ubase, 0, offset, n, offset + n);
  packedmatrix* B = mzd_init_window (Bbase, 0, offset, m, offset + n);
  packedmatrix* W = mzd_copy (NULL, B);

  size_t i,j;
  for (i=0; i<n; ++i){
    for (j=0; j<i;++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }    
  mzd_trsm_upper_right (U, B, 2048);

  mzd_addmul(W, B, U, 2048);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (W,i,j)){
	printf("W [%d,%d] = %d\n", i,j, mzd_read_bit (W,i,j));
	status = 1;
      }
    }
  if (status){
    printf("W = \n");
    mzd_print_matrix(W);
  }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, B, U, 2048);
  mzd_copy (B, W);

  for ( i=0; i<2048; ++i)
    for ( j=0; j<2048/RADIX; ++j){
      if (Bbase->values [Bbase->rowswap[i] + j] != Bbasecopy->values [Bbasecopy->rowswap[i] + j]){
	printf("Bbase [%d,%d] = %llX\n", i,j, Bbase->values[Bbase->rowswap[i]+j]);
	printf("Bbasecopy [%d,%d] = %llX\n", i,j, Bbasecopy->values[Bbasecopy->rowswap[i]+j]);
	status = 1;
      }
    }
  if (status){
    printf("Bbase = \n");
    mzd_print_matrix(Bbase);
    printf("Bbasecopy] = \n");
    mzd_print_matrix(Bbasecopy);
  }
  mzd_free_window (U);
  mzd_free_window (B);
  mzd_free (W);
  mzd_free(Ubase);
  mzd_free(Bbase);
  mzd_free(Bbasecopy);

  if (!status)
    printf("passed\n");
  return status;
}

int testTRSMLowerLeft (int m, int n, int offsetL, int offsetB){
  packedmatrix* Lbase = mzd_init (2048,2048);
  packedmatrix* Bbase = mzd_init (2048,2048);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  packedmatrix* Bbasecopy = mzd_copy (NULL, Bbase);

  packedmatrix* L = mzd_init_window (Lbase, 0, offsetL, m, offsetL + m);
  packedmatrix* B = mzd_init_window (Bbase, 0, offsetB, m, offsetB + n);
  packedmatrix* W = mzd_copy (NULL, B);

  size_t i,j;
  for (i=0; i<m; ++i){
    for (j=i+1; j<m;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }    
  mzd_trsm_lower_left (L, B, 2048);
  
  mzd_addmul(W, L, B, 2048);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (W,i,j)){
	printf("W [%d,%d] = %d\n", i,j, mzd_read_bit (W,i,j));
	status = 1;
      }
    }
  if (status){
    printf("W = \n");
    mzd_print_matrix(W);
    printf("L = \n");
    mzd_print_matrix(L);
    printf("B = \n");
    mzd_print_matrix(B);
  }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, L, B, 2048);

  mzd_copy (B, W);

  for ( i=0; i<2048; ++i)
    for ( j=0; j<2048/RADIX; ++j){
      if (Bbase->values [Bbase->rowswap[i] + j] != Bbasecopy->values [Bbasecopy->rowswap[i] + j]){
	printf("Bbase     [%d,%d] = %llX\n", i,j, Bbase->values[Bbase->rowswap[i]+j]);
	printf("Bbasecopy [%d,%d] = %llX\n", i,j, Bbasecopy->values[Bbasecopy->rowswap[i]+j]);
	status = 1;
      }
    }
  mzd_free_window (L);
  mzd_free_window (B);
  mzd_free_window (W);
  mzd_free(Lbase);
  mzd_free(Bbase);
  mzd_free(Bbasecopy);

  if (!status)
    printf("passed\n");
  return status;
}

int main(int argc, char **argv) {
  m4ri_build_all_codes();

  int status = 0;

  printf("UpperRight: small, even placed ... ");
  status += testTRSMUpperRight (57, 10, 0);
  printf("UpperRight: large, even placed ... ");
  status += testTRSMUpperRight (57, 150, 0);
  printf("UpperRight: small, odd placed ... ");
  status += testTRSMUpperRight (57, 3, 4);
  printf("UpperRight: medium, odd placed ... ");
  status += testTRSMUpperRight (57, 4, 62);
  printf("UpperRight: large, odd placed ... ");
  status += testTRSMUpperRight (57, 80, 60);
  printf("UpperRight: larger, odd placed ... ");
  status += testTRSMUpperRight (1577, 1802, 189);

  printf("LowerLeft: small L even, small B even ... ");
  status += testTRSMLowerLeft (10, 20, 0, 0);
  printf("LowerLeft: small L even, large B even ... ");
  status += testTRSMLowerLeft (10, 80, 0, 0);

  printf("LowerLeft: small L even, small B odd ... ");
  status += testTRSMLowerLeft (10, 20, 0, 15);
  printf("LowerLeft: small L even, large B odd ... ");
  status += testTRSMLowerLeft (10, 80, 0, 15);

  printf("LowerLeft: small L odd, small B even ... ");
  status += testTRSMLowerLeft (10, 20, 15, 0);
  printf("LowerLeft: small L odd, large B even ... ");
  status += testTRSMLowerLeft (10, 80, 15, 0);

  printf("LowerLeft: small L odd, small B odd ... ");
  status += testTRSMLowerLeft (10, 20, 15, 20);
  printf("LowerLeft: small L odd, large B odd ... ");
  status += testTRSMLowerLeft (10, 80, 15, 20);

  printf("LowerLeft: large L even, small B even ... ");
  status += testTRSMLowerLeft (70, 20, 0, 0);
  printf("LowerLeft: large L even, large B even ... ");
  status += testTRSMLowerLeft (70, 80, 0, 0);

  printf("LowerLeft: large L even, small B odd ... ");
  status += testTRSMLowerLeft (70, 10, 0, 15);
  printf("LowerLeft: large L even, large B odd ... ");
  status += testTRSMLowerLeft (70, 80, 0, 15);

  printf("LowerLeft: large L odd, small B even ... ");
  status += testTRSMLowerLeft (70, 20, 15, 0);
  printf("LowerLeft: large L odd, large B even ... ");
  status += testTRSMLowerLeft (70, 80, 15, 0);

  printf("LowerLeft: large L odd, small B odd ... ");
  status += testTRSMLowerLeft (70, 20, 15, 20);
  printf("LowerLeft: large L odd, large B odd ... ");
  status += testTRSMLowerLeft (70, 80, 15, 20);

  printf("LowerLeft: larger L odd, larger B odd ... ");
  status += testTRSMLowerLeft (770, 1600, 75, 89);
  printf("LowerLeft: larger L odd, larger B odd ... ");
  status += testTRSMLowerLeft (1764, 1345, 198, 123);

  
  if (!status) {
    printf("All tests passed.\n");
  }

  m4ri_destroy_all_codes();
  return 0;
}
