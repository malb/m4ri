#include "config.h"
#include <stdlib.h>
#include "m4ri.h"

//#define RANDOMIZE

int test_trsm_upper_right (rci_t m, rci_t n, int offset, const char* description){
  printf("upper_right: %s  m: %4d n: %4d offset: %4d ... ", description, m.val(), n.val(), offset);

  mzd_t* Ubase = mzd_init (2048U, 2048U);
  mzd_t* Bbase = mzd_init (2048U, 2048U);
  mzd_randomize (Ubase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* U = mzd_init_window (Ubase, 0, rci_t(0) + offset, n, n + offset);
  mzd_t* B = mzd_init_window (Bbase, 0, rci_t(0) + offset, m, n + offset);
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
    for (wi_t j = 0U; j < rci_t(2048U) / RADIX; ++j){
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
  printf("lower_right: %s  m: %4d n: %4d offset: %4d ... ", description, m.val(), n.val(), offset);
  mzd_t* Lbase = mzd_init (2048U, 2048U);
  mzd_t* Bbase = mzd_init (2048U, 2048U);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* L = mzd_init_window (Lbase, 0, rci_t(0) + offset, n, n + offset);
  mzd_t* B = mzd_init_window (Bbase, 0, rci_t(0) + offset, m, n + offset);
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
    for (wi_t j = 0U; j < rci_t(2048U) / RADIX; ++j){
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
  printf("lower_left: %s  m: %4d n: %4d offset L: %4d offset B: %4d ... ", description, m.val(), n.val(), offsetL, offsetB);
  mzd_t* Lbase = mzd_init (2048U, 2048U);
  mzd_t* Bbase = mzd_init (2048U, 2048U);
  mzd_randomize (Lbase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* L = mzd_init_window (Lbase, 0, rci_t(0) + offsetL, m, m + offsetL);
  mzd_t* B = mzd_init_window (Bbase, 0, rci_t(0) + offsetB, m, n + offsetB);
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
    for (wi_t j = 0U; j < rci_t(2048U) / RADIX; ++j){
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
  printf("upper_left: %s  m: %4d n: %4d offset U: %4d offset B: %4d ... ", description, m.val(), n.val(), offsetU, offsetB);
  mzd_t* Ubase = mzd_init (2048U, 2048U);
  mzd_t* Bbase = mzd_init (2048U, 2048U);
  mzd_randomize (Ubase);
  mzd_randomize (Bbase);
  mzd_t* Bbasecopy = mzd_copy (NULL, Bbase);

  mzd_t* U = mzd_init_window (Ubase, 0, rci_t(0) + offsetU, m, m + offsetU);
  mzd_t* B = mzd_init_window (Bbase, 0, rci_t(0) + offsetB, m, n + offsetB);
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
    for (wi_t j = 0U; j < rci_t(2048U) / RADIX; ++j){
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

  status += test_trsm_upper_right(  63U,   63U,   0, "  word boundaries, even");
  status += test_trsm_upper_right(  64U,   64U,   0, "  word boundaries, even");
  status += test_trsm_upper_right(  65U,   65U,   0, "  word boundaries, even");
  status += test_trsm_upper_right(  53U,   53U,  10, "word boundaries, offset");
  status += test_trsm_upper_right(  54U,   54U,  10, "word boundaries, offset");
  status += test_trsm_upper_right(  55U,   55U,  10, "word boundaries, offset");
  status += test_trsm_upper_right(  57U,   10U,   0, "     small, even placed");
  status += test_trsm_upper_right(  57U,  150U,   0, "     large, even placed");
  status += test_trsm_upper_right(  57U,    3U,   4, "      small, odd placed");
  status += test_trsm_upper_right(  57U,    4U,  62, "     medium, odd placed");
  status += test_trsm_upper_right(  57U,   80U,  60, "      large, odd placed");
  status += test_trsm_upper_right(1577U, 1802U, 189, "     larger, odd placed");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_upper_right(random() & 2047, random() & 2047, random() & 63, "randomized");
  } 
#endif

  printf("\n");

  status += test_trsm_lower_right(  63U,   63U,  0,"  word boundaries, even");
  status += test_trsm_lower_right(  64U,   64U,  0,"  word boundaries, even");
  status += test_trsm_lower_right(  65U,   65U,  0,"  word boundaries, even");
  status += test_trsm_lower_right(  53U,   53U, 10,"word boundaries, offset");
  status += test_trsm_lower_right(  54U,   54U, 10,"word boundaries, offset");
  status += test_trsm_lower_right(  55U,   55U, 10,"word boundaries, offset");
  status += test_trsm_lower_right(  57U,   10U,  0,"     small, even placed");
  status += test_trsm_lower_right(  57U,  150U,  0,"     large, even placed");
  status += test_trsm_lower_right(  57U,    3U,  4,"      small, odd placed");
  status += test_trsm_lower_right(  57U,    4U, 62,"     medium, odd placed");
  status += test_trsm_lower_right(  57U,   80U, 60,"      large, odd placed");
  status += test_trsm_lower_right(1577U, 1802U,189,"     larger, odd placed");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_lower_right(random() & 2047, random() & 2047, random() & 63, "randomized");
  } 
#endif

  printf("\n");

  status += test_trsm_lower_left(  63U,   63U,   0,   0, "      word boundaries, even");
  status += test_trsm_lower_left(  64U,   64U,   0,   0, "      word boundaries, even");
  status += test_trsm_lower_left(  65U,   65U,   0,   0, "      word boundaries, even");
  status += test_trsm_lower_left(  53U,   53U,  10,  10, "    word boundaries, offset");
  status += test_trsm_lower_left(  54U,   54U,  10,  10, "    word boundaries, offset");
  status += test_trsm_lower_left(  55U,   55U,  10,  10, "    word boundaries, offset");
  status += test_trsm_lower_left(  10U,   20U,   0,   0, " small L even, small B even");
  status += test_trsm_lower_left(  10U,   80U,   0,   0, " small L even, large B even");
  status += test_trsm_lower_left(  10U,   20U,   0,  15, "  small L even, small B odd");
  status += test_trsm_lower_left(  10U,   80U,   0,  15, "  small L even, large B odd");
  status += test_trsm_lower_left(  10U,   20U,  15,   0, "  small L odd, small B even");
  status += test_trsm_lower_left(  10U,   80U,  15,   0, "  small L odd, large B even");
  status += test_trsm_lower_left(  10U,   20U,  15,  20, "   small L odd, small B odd");
  status += test_trsm_lower_left(  10U,   80U,  15,  20, "   small L odd, large B odd");
  status += test_trsm_lower_left(  70U,   20U,   0,   0, " large L even, small B even");
  status += test_trsm_lower_left(  70U,   80U,   0,   0, " large L even, large B even");
  status += test_trsm_lower_left(  70U,   10U,   0,  15, "  large L even, large B odd");
  status += test_trsm_lower_left(  70U,   80U,   0,  15, "  large L even, large B odd");
  status += test_trsm_lower_left(  70U,   20U,  15,   0, "  large L odd, small B even");
  status += test_trsm_lower_left(  70U,   80U,  15,   0, "  large L odd, large B even");
  status += test_trsm_lower_left(  70U,   20U,  15,  20, "   large L odd, small B odd");
  status += test_trsm_lower_left(  70U,   80U,  15,  20, "   large L odd, large B odd");
  status += test_trsm_lower_left( 770U, 1600U,  75,  89, " larger L odd, larger B odd");
  status += test_trsm_lower_left(1764U, 1345U, 198, 123, " larger L odd, larger B odd");

#ifdef RANDOMIZE
  for(size_t i=0; i<256; i++) {
    status += test_trsm_lower_left(random() & 2047, random() & 2047, random() & 63, random() & 63, "randomized");
  } 
#endif


  printf("\n");

  status += test_trsm_upper_left(  63U,  63U,  0,  0,"     word boundaries, even");
  status += test_trsm_upper_left(  64U,  64U,  0,  0,"     word boundaries, even");
  status += test_trsm_upper_left(  65U,  65U,  0,  0,"     word boundaries, even");
  status += test_trsm_upper_left(  53U,  53U, 10, 10,"   word boundaries, offset");
  status += test_trsm_upper_left(  54U,  54U, 10, 10,"   word boundaries, offset");
  status += test_trsm_upper_left(  55U,  55U, 10, 10,"   word boundaries, offset");
  status += test_trsm_upper_left(  10U,  20U,  0,  0,"small U even, small B even");
  status += test_trsm_upper_left(  10U,  80U,  0,  0,"small U even, large B even");
  status += test_trsm_upper_left(  10U,  20U,  0, 15," small U even, small B odd");
  status += test_trsm_upper_left(  10U,  80U,  0, 15," small U even, large B odd");
  status += test_trsm_upper_left(  10U,  20U, 15,  0," small U odd, small B even");
  status += test_trsm_upper_left(  10U,  80U, 15,  0," small U odd, large B even");
  status += test_trsm_upper_left(  10U,  20U, 15, 20,"  small U odd, small B odd");
  status += test_trsm_upper_left(  10U,  80U, 15, 20,"  small U odd, large B odd");
  status += test_trsm_upper_left(  70U,  20U,  0,  0,"large U even, small B even");
  status += test_trsm_upper_left(  63U,   1U,  0,  0,"                          ");
  status += test_trsm_upper_left(  70U,  80U,  0,  0,"large U even, large B even");
  status += test_trsm_upper_left(  70U,  10U,  0, 15," large U even, small B odd");
  status += test_trsm_upper_left(  70U,  80U,  0, 15," large U even, large B odd");
  status += test_trsm_upper_left(  70U,  20U, 15,  0," large U odd, small B even");
  status += test_trsm_upper_left(  70U,  80U, 15,  0," large U odd, large B even");
  status += test_trsm_upper_left(  70U,  20U, 15, 20,"  large U odd, small B odd");
  status += test_trsm_upper_left(  70U,  80U, 15, 20,"  large U odd, large B odd");
  status += test_trsm_upper_left( 770U,1600U, 75, 89,"larger U odd, larger B odd");
  status += test_trsm_upper_left(1764U,1345U,198,123,"larger U odd, larger B odd");

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
