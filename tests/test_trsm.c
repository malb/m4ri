#include <m4ri/config.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>

int test_trsm_upper_right (rci_t m, rci_t n, int offset){
  printf("upper_right:: m: %4d n: %4d offset: %4d ... ", m, n, offset);

  mzd_t* Ubase = mzd_init (2048, 2048);
  mzd_t* Bbase = mzd_init (2048, 2048);
  mzd_randomize(Ubase);
  mzd_randomize(Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* U = mzd_init_window (Ubase, 0, offset, n, n + offset);
  mzd_t* B = mzd_init_window (Bbase, 0, offset, m, n + offset);
  mzd_t* W = mzd_copy (NULL, B);

  for (rci_t i = 0; i < n; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }
  mzd_trsm_upper_right (U, B, 2048);

  mzd_addmul(W, B, U, 2048);

  int status = 0;
  for (rci_t i = 0; i < m; ++i)
    for (rci_t j = 0; j < n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, B, U, 2048);
  mzd_copy (B, W);

  for (rci_t i = 0; i < 2048; ++i)
    for (wi_t j = 0; j < 2048 / m4ri_radix; ++j){
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

int test_trsm_lower_right (rci_t m, rci_t n, int offset){
  printf("lower_right:: m: %4d n: %4d offset: %4d ... ", m, n, offset);
  mzd_t* Lbase = mzd_init (2048, 2048);
  mzd_t* Bbase = mzd_init (2048, 2048);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* L = mzd_init_window(Lbase, 0, offset, n, n + offset);
  mzd_t* B = mzd_init_window (Bbase, 0, offset, m, n + offset);
  mzd_t* W = mzd_copy (NULL, B);

  for (rci_t i = 0; i < n; ++i){
    for (rci_t j = i + 1; j < n; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_trsm_lower_right (L, B, 2048);

  mzd_addmul(W, B, L, 2048);

  int status = 0;
  for (rci_t i = 0; i < m; ++i)
    for (rci_t j = 0; j < n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, B, L, 2048);
  mzd_copy (B, W);

  for (rci_t i = 0; i < 2048; ++i)
    for (wi_t j = 0; j < 2048 / m4ri_radix; ++j){
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


int test_trsm_lower_left (rci_t m, rci_t n, int offsetL, int offsetB){
  printf("lower_left:: m: %4d n: %4d offset L: %4d offset B: %4d ... ", m, n, offsetL, offsetB);
  mzd_t* Lbase = mzd_init (2048, 2048);
  mzd_t* Bbase = mzd_init (2048, 2048);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* L = mzd_init_window (Lbase, 0, offsetL, m, m + offsetL);
  mzd_t* B = mzd_init_window (Bbase, 0, offsetB, m, n + offsetB);
  mzd_t* W = mzd_copy (NULL, B);

  for (rci_t i = 0; i < m; ++i){
    for (rci_t j = i + 1; j < m; ++j)
      mzd_write_bit(L,i,j, 0);
    mzd_write_bit(L,i,i, 1);
  }
  mzd_trsm_lower_left(L, B, 2048);
  
  mzd_addmul(W, L, B, 2048);

  int status = 0;
  for (rci_t i = 0; i < m; ++i)
    for (rci_t j = 0; j < n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }

  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, L, B, 2048);

  mzd_copy (B, W);

  for (rci_t i = 0; i < 2048; ++i)
    for (wi_t j = 0; j < 2048 / m4ri_radix; ++j){
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



int test_trsm_upper_left (rci_t m, rci_t n, int offsetU, int offsetB) {
  printf("upper_left:: m: %4d n: %4d offset U: %4d offset B: %4d ... ", m, n, offsetU, offsetB);
  mzd_t* Ubase = mzd_init (2048, 2048);
  mzd_t* Bbase = mzd_init (2048, 2048);
  mzd_randomize (Ubase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* U = mzd_init_window (Ubase, 0, offsetU, m, m + offsetU);
  mzd_t* B = mzd_init_window (Bbase, 0, offsetB, m, n + offsetB);
  mzd_t* W = mzd_copy (NULL, B);

  for (rci_t i = 0; i < m; ++i){
    for (rci_t j = 0; j < i; ++j)
      mzd_write_bit(U,i,j, 0);
    mzd_write_bit(U,i,i, 1);
  }    
  mzd_trsm_upper_left(U, B, 2048);
  
  mzd_addmul(W, U, B, 2048);

  int status = 0;
  for (rci_t i = 0; i < m; ++i)
    for (rci_t j = 0; j < n; ++j){
      if (mzd_read_bit (W,i,j)){
	status = 1;
      }
    }
  // Verifiying that nothing has been changed around the submatrices
  mzd_addmul(W, U, B, 2048);

  mzd_copy (B, W);

  for (rci_t i = 0; i < 2048; ++i)
    for (wi_t j = 0; j < 2048 / m4ri_radix; ++j){
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

int main() {
  int status = 0;

  srandom(17);

  status += test_trsm_upper_right(  63,   63,   0);
  status += test_trsm_upper_right(  64,   64,   0);
  status += test_trsm_upper_right(  65,   65,   0);
  status += test_trsm_upper_right(  53,   53,   0);
  status += test_trsm_upper_right(  54,   54,   0);
  status += test_trsm_upper_right(  55,   55,   0);
  status += test_trsm_upper_right(  57,   10,   0);
  status += test_trsm_upper_right(  57,  150,   0);
  status += test_trsm_upper_right(  57,    3,   0);
  status += test_trsm_upper_right(  57,    4,  64);
  status += test_trsm_upper_right(  57,   80,  64);
  status += test_trsm_upper_right(1577, 1802, 128);

  printf("\n");

  status += test_trsm_lower_right(  63,   63,   0);
  status += test_trsm_lower_right(  64,   64,   0);
  status += test_trsm_lower_right(  65,   65,   0);
  status += test_trsm_lower_right(  53,   53,   0);
  status += test_trsm_lower_right(  54,   54,   0);
  status += test_trsm_lower_right(  55,   55,   0);
  status += test_trsm_lower_right(  57,   10,   0);
  status += test_trsm_lower_right(  57,  150,   0);
  status += test_trsm_lower_right(  57,    3,   0);
  status += test_trsm_lower_right(  57,    4,  64);
  status += test_trsm_lower_right(  57,   80,  64);
  status += test_trsm_lower_right(1577, 1802, 128);

  printf("\n");

  status += test_trsm_lower_left(  63,   63,   0,   0);
  status += test_trsm_lower_left(  64,   64,   0,   0);
  status += test_trsm_lower_left(  65,   65,   0,   0);
  status += test_trsm_lower_left(  53,   53,   0,   0);
  status += test_trsm_lower_left(  54,   54,   0,   0);
  status += test_trsm_lower_left(  55,   55,   0,   0);
  status += test_trsm_lower_left(  10,   20,   0,   0);
  status += test_trsm_lower_left(  10,   80,   0,   0);
  status += test_trsm_lower_left(  10,   20,   0,   0);
  status += test_trsm_lower_left(  10,   80,   0,   0);
  status += test_trsm_lower_left(  10,   20,   0,   0);
  status += test_trsm_lower_left(  10,   80,   0,   0);
  status += test_trsm_lower_left(  10,   20,   0,   0);
  status += test_trsm_lower_left(  10,   80,   0,   0);
  status += test_trsm_lower_left(  70,   20,   0,   0);
  status += test_trsm_lower_left(  70,   80,   0,   0);
  status += test_trsm_lower_left(  70,   10,   0,   0);
  status += test_trsm_lower_left(  70,   80,   0,   0);
  status += test_trsm_lower_left(  70,   20,   0,   0);
  status += test_trsm_lower_left(  70,   80,   0,   0);
  status += test_trsm_lower_left(  70,   20,   0,   0);
  status += test_trsm_lower_left(  70,   80,  64,  64);
  status += test_trsm_lower_left( 770, 1600,  64, 128);
  status += test_trsm_lower_left(1764, 1345, 256,  64);

  printf("\n");

  status += test_trsm_upper_left(   63,   63,   0,   0);
  status += test_trsm_upper_left(   64,   64,   0,   0);
  status += test_trsm_upper_left(   65,   65,   0,   0);
  status += test_trsm_upper_left(   53,   53,   0,   0);
  status += test_trsm_upper_left(   54,   54,   0,   0);
  status += test_trsm_upper_left(   55,   55,   0,   0);
  status += test_trsm_upper_left(   10,   20,   0,   0);
  status += test_trsm_upper_left(   10,   80,   0,   0);
  status += test_trsm_upper_left(   10,   20,   0,   0);
  status += test_trsm_upper_left(   10,   80,   0,   0);
  status += test_trsm_upper_left(   10,   20,   0,   0);
  status += test_trsm_upper_left(   10,   80,   0,   0);
  status += test_trsm_upper_left(   10,   20,   0,   0);
  status += test_trsm_upper_left(   10,   80,   0,   0);
  status += test_trsm_upper_left(   70,   20,   0,   0);
  status += test_trsm_upper_left(   63,    1,   0,   0);
  status += test_trsm_upper_left(   70,   80,   0,   0);
  status += test_trsm_upper_left(   70,   10,   0,   0);
  status += test_trsm_upper_left(   70,   80,   0,   0);
  status += test_trsm_upper_left(   70,   20,   0,   0);
  status += test_trsm_upper_left(   70,   80,   0,   0);
  status += test_trsm_upper_left(   70,   20,   0,   0);
  status += test_trsm_upper_left(   70,   80,  64,  64);
  status += test_trsm_upper_left(  770, 1600,  64, 128);
  status += test_trsm_upper_left( 1764, 1345, 256,  64);

  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
