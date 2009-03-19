#include <stdlib.h>
#include "m4ri/m4ri.h"


int test_trsm_upper_right (int m, int n, int offset, const char* description){
  printf("upper_right: %s  m: %4d n: %4d offset: %4d ... ",description, m, n, offset);

  mzd_t* Ubase = mzd_init (2048,2048);
  mzd_t* Bbase = mzd_init (2048,2048);
  mzd_randomize (Ubase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* U = mzd_init_window (Ubase, 0, offset, n, offset + n);
  mzd_t* B = mzd_init_window (Bbase, 0, offset, m, offset + n);
  mzd_t* W = mzd_copy (NULL, B);

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
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, B, U, 2048);
  mzd_copy (B, W);

  for ( i=0; i<2048; ++i)
    for ( j=0; j<2048/RADIX; ++j){
      if (Bbase->rows[i][j] != Bbasecopy->rows[i][j]){
	status = 1;
      }
    }
  mzd_free_window (U);
  mzd_free_window (B);
  mzd_free (W);
  mzd_free(Ubase);
  mzd_free(Bbase);
  mzd_free(Bbasecopy);

  if (!status)
    printf("passed\n");
  else
    printf("FAILED\n");
  return status;
}

int test_trsm_lower_right (int m, int n, int offset, const char *description){
  printf("lower_right: %s  m: %4d n: %4d offset: %4d ... ",description, m, n, offset);
  mzd_t* Lbase = mzd_init (2048,2048);
  mzd_t* Bbase = mzd_init (2048,2048);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* L = mzd_init_window (Lbase, 0, offset, n, offset + n);
  mzd_t* B = mzd_init_window (Bbase, 0, offset, m, offset + n);
  mzd_t* W = mzd_copy (NULL, B);

  size_t i,j;
  for (i=0; i<n; ++i){
    for (j=i+1; j<n;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_trsm_lower_right (L, B, 2048);

  mzd_addmul(W, B, L, 2048);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, B, L, 2048);
  mzd_copy (B, W);

  for ( i=0; i<2048; ++i)
    for ( j=0; j<2048/RADIX; ++j){
      if (Bbase->rows[i][j] != Bbasecopy->rows[i][j]){
	status = 1;
      }
    }
  mzd_free_window (L);
  mzd_free_window (B);
  mzd_free (W);
  mzd_free(Lbase);
  mzd_free(Bbase);
  mzd_free(Bbasecopy);

  if (!status)
    printf("passed\n");
  else
    printf("FAILED\n");
  return status;
}


int test_trsm_lower_left (int m, int n, int offsetL, int offsetB){
  mzd_t* Lbase = mzd_init (2048,2048);
  mzd_t* Bbase = mzd_init (2048,2048);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* L = mzd_init_window (Lbase, 0, offsetL, m, offsetL + m);
  mzd_t* B = mzd_init_window (Bbase, 0, offsetB, m, offsetB + n);
  mzd_t* W = mzd_copy (NULL, B);

  size_t i,j;
  for (i=0; i<m; ++i){
    for (j=i+1; j<m;++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_trsm_lower_left(L, B, 2048);
  
  mzd_addmul(W, L, B, 2048);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, L, B, 2048);

  mzd_copy (B, W);

  for ( i=0; i<2048; ++i)
    for ( j=0; j<2048/RADIX; ++j){
      if (Bbase->rows[i][j] != Bbasecopy->rows[i][j]){
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
    printf(" ... passed\n");
  else
    printf(" ... FAILED\n");
  return status;
}



int test_trsm_upper_left (int m, int n, int offsetU, int offsetB, const char *description) {
  printf("upper_left: %s  m: %4d n: %4d offset: %4d ... ",description, m, n, offsetU);
  mzd_t* Ubase = mzd_init (2048,2048);
  mzd_t* Bbase = mzd_init (2048,2048);
  mzd_randomize (Ubase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* U = mzd_init_window (Ubase, 0, offsetU, m, offsetU + m);
  mzd_t* B = mzd_init_window (Bbase, 0, offsetB, m, offsetB + n);
  mzd_t* W = mzd_copy (NULL, B);

  size_t i,j;
  for (i=0; i<m; ++i){
    for (j=0; j<i;++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }    
  mzd_trsm_upper_left(U, B, 2048);
  
  mzd_addmul(W, U, B, 2048);

  int status = 0;
  for ( i=0; i<m; ++i)
    for ( j=0; j<n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, U, B, 2048);

  mzd_copy (B, W);

  for ( i=0; i<2048; ++i)
    for ( j=0; j<2048/RADIX; ++j){
      if (Bbase->rows[i][j] != Bbasecopy->rows[i][j]){
	status = 1;
      }
    }
  mzd_free_window (U);
  mzd_free_window (B);
  mzd_free_window (W);
  mzd_free(Ubase);
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

  status += test_trsm_upper_right(  57,   10,   0, "small, even placed");
  status += test_trsm_upper_right(  57,  150,   0, "large, even placed");
  status += test_trsm_upper_right(  57,    3,   4, " small, odd placed");
  status += test_trsm_upper_right(  57,    4,  62, "medium, odd placed");
  status += test_trsm_upper_right(  57,   80,  60, " large, odd placed");
  status += test_trsm_upper_right(1577, 1802, 189, "larger, odd placed");

  printf("\n");

  status += test_trsm_lower_right(  57,   10,  0,"small, even placed");
  status += test_trsm_lower_right(  57,  150,  0,"large, even placed");
  status += test_trsm_lower_right(  57,    3,  4," small, odd placed");
  status += test_trsm_lower_right(  57,    4, 62,"medium, odd placed");
  status += test_trsm_lower_right(  57,   80, 60," large, odd placed");
  status += test_trsm_lower_right(1577, 1802,189,"larger, odd placed");

  printf("\n");

  printf("LowerLeft: small L even, small B even ");
  status += test_trsm_lower_left (10, 20, 0, 0);
  printf("LowerLeft: small L even, large B even ");
  status += test_trsm_lower_left (10, 80, 0, 0);
  printf("LowerLeft: small L even, small B odd  ");
  status += test_trsm_lower_left (10, 20, 0, 15);
  printf("LowerLeft: small L even, large B odd  ");
  status += test_trsm_lower_left (10, 80, 0, 15);
  printf("LowerLeft: small L odd, small B even  ");
  status += test_trsm_lower_left (10, 20, 15, 0);
  printf("LowerLeft: small L odd, large B even  ");
  status += test_trsm_lower_left (10, 80, 15, 0);
  printf("LowerLeft: small L odd, small B odd   ");
  status += test_trsm_lower_left (10, 20, 15, 20);
  printf("LowerLeft: small L odd, large B odd   ");
  status += test_trsm_lower_left (10, 80, 15, 20);
  printf("LowerLeft: large L even, small B even ");
  status += test_trsm_lower_left (70, 20, 0, 0);
  printf("LowerLeft: large L even, large B even ");
  status += test_trsm_lower_left (70, 80, 0, 0);
  printf("LowerLeft: large L even, small B odd  ");
  status += test_trsm_lower_left (70, 10, 0, 15);
  printf("LowerLeft: large L even, large B odd  ");
  status += test_trsm_lower_left (70, 80, 0, 15);
  printf("LowerLeft: large L odd, small B even  ");
  status += test_trsm_lower_left (70, 20, 15, 0);
  printf("LowerLeft: large L odd, large B even  ");
  status += test_trsm_lower_left (70, 80, 15, 0);
  printf("LowerLeft: large L odd, small B odd   ");
  status += test_trsm_lower_left (70, 20, 15, 20);
  printf("LowerLeft: large L odd, large B odd   ");
  status += test_trsm_lower_left (70, 80, 15, 20);
  printf("LowerLeft: larger L odd, larger B odd ");
  status += test_trsm_lower_left (770, 1600, 75, 89);
  printf("LowerLeft: larger L odd, larger B odd ");
  status += test_trsm_lower_left (1764, 1345, 198, 123);

  printf("\n");

  status += test_trsm_upper_left(  10,  20,  0,  0,"small U even, small B even");
  status += test_trsm_upper_left(  10,  80,  0,  0,"small U even, large B even");
  status += test_trsm_upper_left(  10,  20,  0, 15," small U even, small B odd");
  status += test_trsm_upper_left(  10,  80,  0, 15," small U even, large B odd");
  status += test_trsm_upper_left(  10,  20, 15,  0," small U odd, small B even");
  status += test_trsm_upper_left(  10,  80, 15,  0," small U odd, large B even");
  status += test_trsm_upper_left(  10,  20, 15, 20,"  small U odd, small B odd");
  status += test_trsm_upper_left(  10,  80, 15, 20,"  small U odd, large B odd");
  status += test_trsm_upper_left(  70,  20,  0,  0,"large U even, small B even");
  status += test_trsm_upper_left(  63,   1,  0,  0,"                          ");
  status += test_trsm_upper_left(  70,  80,  0,  0,"large U even, large B even");
  status += test_trsm_upper_left(  70,  10,  0, 15," large U even, small B odd");
  status += test_trsm_upper_left(  70,  80,  0, 15," large U even, large B odd");
  status += test_trsm_upper_left(  70,  20, 15,  0," large U odd, small B even");
  status += test_trsm_upper_left(  70,  80, 15,  0," large U odd, large B even");
  status += test_trsm_upper_left(  70,  20, 15, 20,"  large U odd, small B odd");
  status += test_trsm_upper_left(  70,  80, 15, 20,"  large U odd, large B odd");
  status += test_trsm_upper_left( 770,1600, 75, 89,"larger U odd, larger B odd");
  status += test_trsm_upper_left(1764,1345,198,123,"larger U odd, larger B odd");

  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
