#include <m4ri/config.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>

//#define RANDOMIZE

int test_trsm_upper_right (rci_t m, rci_t n, int offset, const char* description){
  printf("upper_right: %s  m: %4d n: %4d offset: %4d ... ", description, m, n, offset);

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

int test_trsm_lower_right (rci_t m, rci_t n, int offset, const char *description){
  printf("lower_right: %s  m: %4d n: %4d offset: %4d ... ", description, m, n, offset);
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


int test_trsm_lower_left (rci_t m, rci_t n, int offsetL, int offsetB, const char *description){
  printf("lower_left: %s  m: %4d n: %4d offset L: %4d offset B: %4d ... ", description, m, n, offsetL, offsetB);
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



int test_trsm_upper_left (rci_t m, rci_t n, int offsetU, int offsetB, const char *description) {
  printf("upper_left: %s  m: %4d n: %4d offset U: %4d offset B: %4d ... ", description, m, n, offsetU, offsetB);
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

  status += test_trsm_upper_right(  63,   63,   0, "  word boundaries, even");
  status += test_trsm_upper_right(  64,   64,   0, "  word boundaries, even");
  status += test_trsm_upper_right(  65,   65,   0, "  word boundaries, even");
  status += test_trsm_upper_right(  53,   53,  10, "word boundaries, offset");
  status += test_trsm_upper_right(  54,   54,  10, "word boundaries, offset");
  status += test_trsm_upper_right(  55,   55,  10, "word boundaries, offset");
  status += test_trsm_upper_right(  57,   10,   0, "     small, even placed");
  status += test_trsm_upper_right(  57,  150,   0, "     large, even placed");
  status += test_trsm_upper_right(  57,    3,   4, "      small, odd placed");
  status += test_trsm_upper_right(  57,    4,  62, "     medium, odd placed");
  status += test_trsm_upper_right(  57,   80,  60, "      large, odd placed");
  status += test_trsm_upper_right(1577, 1802, 189, "     larger, odd placed");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_upper_right(random() & 2047, random() & 2047, random() & 63, "randomized");
  } 
#endif

  printf("\n");

  status += test_trsm_lower_right(  63,   63,  0,"  word boundaries, even");
  status += test_trsm_lower_right(  64,   64,  0,"  word boundaries, even");
  status += test_trsm_lower_right(  65,   65,  0,"  word boundaries, even");
  status += test_trsm_lower_right(  53,   53, 10,"word boundaries, offset");
  status += test_trsm_lower_right(  54,   54, 10,"word boundaries, offset");
  status += test_trsm_lower_right(  55,   55, 10,"word boundaries, offset");
  status += test_trsm_lower_right(  57,   10,  0,"     small, even placed");
  status += test_trsm_lower_right(  57,  150,  0,"     large, even placed");
  status += test_trsm_lower_right(  57,    3,  4,"      small, odd placed");
  status += test_trsm_lower_right(  57,    4, 62,"     medium, odd placed");
  status += test_trsm_lower_right(  57,   80, 60,"      large, odd placed");
  status += test_trsm_lower_right(1577, 1802,189,"     larger, odd placed");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_lower_right(random() & 2047, random() & 2047, random() & 63, "randomized");
  } 
#endif

  printf("\n");

  status += test_trsm_lower_left(  63,   63,   0,   0, "      word boundaries, even");
  status += test_trsm_lower_left(  64,   64,   0,   0, "      word boundaries, even");
  status += test_trsm_lower_left(  65,   65,   0,   0, "      word boundaries, even");
  status += test_trsm_lower_left(  53,   53,  10,  10, "    word boundaries, offset");
  status += test_trsm_lower_left(  54,   54,  10,  10, "    word boundaries, offset");
  status += test_trsm_lower_left(  55,   55,  10,  10, "    word boundaries, offset");
  status += test_trsm_lower_left(  10,   20,   0,   0, " small L even, small B even");
  status += test_trsm_lower_left(  10,   80,   0,   0, " small L even, large B even");
  status += test_trsm_lower_left(  10,   20,   0,  15, "  small L even, small B odd");
  status += test_trsm_lower_left(  10,   80,   0,  15, "  small L even, large B odd");
  status += test_trsm_lower_left(  10,   20,  15,   0, "  small L odd, small B even");
  status += test_trsm_lower_left(  10,   80,  15,   0, "  small L odd, large B even");
  status += test_trsm_lower_left(  10,   20,  15,  20, "   small L odd, small B odd");
  status += test_trsm_lower_left(  10,   80,  15,  20, "   small L odd, large B odd");
  status += test_trsm_lower_left(  70,   20,   0,   0, " large L even, small B even");
  status += test_trsm_lower_left(  70,   80,   0,   0, " large L even, large B even");
  status += test_trsm_lower_left(  70,   10,   0,  15, "  large L even, large B odd");
  status += test_trsm_lower_left(  70,   80,   0,  15, "  large L even, large B odd");
  status += test_trsm_lower_left(  70,   20,  15,   0, "  large L odd, small B even");
  status += test_trsm_lower_left(  70,   80,  15,   0, "  large L odd, large B even");
  status += test_trsm_lower_left(  70,   20,  15,  20, "   large L odd, small B odd");
  status += test_trsm_lower_left(  70,   80,  15,  20, "   large L odd, large B odd");
  status += test_trsm_lower_left( 770, 1600,  75,  89, " larger L odd, larger B odd");
  status += test_trsm_lower_left(1764, 1345, 198, 123, " larger L odd, larger B odd");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_lower_left(random() & 2047, random() & 2047, random() & 63, random() & 63, "randomized");
  } 
#endif


  printf("\n");

  status += test_trsm_upper_left(  63,  63,  0,  0,"    word boundaries, even");
  status += test_trsm_upper_left(  64,  64,  0,  0,"    word boundaries, even");
  status += test_trsm_upper_left(  65,  65,  0,  0,"    word boundaries, even");
  status += test_trsm_upper_left(  53,  53, 10, 10,"  word boundaries, offset");
  status += test_trsm_upper_left(  54,  54, 10, 10,"  word boundaries, offset");
  status += test_trsm_upper_left(  55,  55, 10, 10,"  word boundaries, offset");
  status += test_trsm_upper_left(  10,  20,  0,  0,"small  even, small B even");
  status += test_trsm_upper_left(  10,  80,  0,  0,"small  even, large B even");
  status += test_trsm_upper_left(  10,  20,  0, 15," small  even, small B odd");
  status += test_trsm_upper_left(  10,  80,  0, 15," small  even, large B odd");
  status += test_trsm_upper_left(  10,  20, 15,  0," small  odd, small B even");
  status += test_trsm_upper_left(  10,  80, 15,  0," small  odd, large B even");
  status += test_trsm_upper_left(  10,  20, 15, 20,"  small  odd, small B odd");
  status += test_trsm_upper_left(  10,  80, 15, 20,"  small  odd, large B odd");
  status += test_trsm_upper_left(  70,  20,  0,  0,"large  even, small B even");
  status += test_trsm_upper_left(  63,   1,  0,  0,"                         ");
  status += test_trsm_upper_left(  70,  80,  0,  0,"large  even, large B even");
  status += test_trsm_upper_left(  70,  10,  0, 15," large  even, small B odd");
  status += test_trsm_upper_left(  70,  80,  0, 15," large  even, large B odd");
  status += test_trsm_upper_left(  70,  20, 15,  0," large  odd, small B even");
  status += test_trsm_upper_left(  70,  80, 15,  0," large  odd, large B even");
  status += test_trsm_upper_left(  70,  20, 15, 20,"  large  odd, small B odd");
  status += test_trsm_upper_left(  70,  80, 15, 20,"  large  odd, large B odd");
  status += test_trsm_upper_left( 770,1600, 75, 89,"larger  odd, larger B odd");
  status += test_trsm_upper_left(1764,1345,198,123,"larger  odd, larger B odd");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_upper_left(random() & 2047, random() & 2047, random() & 63, random() & 63, "randomized");
  } 
#endif

  if (!status) {
    printf("All tests passed.\n");
    return 0;
  } else {
    return -1;
  }
}
