/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008 Clement Pernet <clement.pernet@gmail.com>
*
*  Distributed under the terms of the GNU General Public License (GPL)
*  version 2 or higher.
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
*
********************************************************************/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "trsm.h"
#include "strassen.h"
#include "packedmatrix.h"
#include "parity.h"
#include "stdio.h"

#define TRSM_THRESHOLD RADIX

/*****************
 * UPPER RIGHT
 ****************/

/* 
 * This version assumes that the matrices are at an even position on
 * the RADIX grid and that their dimension is a multiple of RADIX.
 */

void _mzd_trsm_upper_right_even(mzd_t *U, mzd_t *B, const int cutoff);

/*
 * Variant where U and B start at an odd bit position. Assumes that
 * U->ncols < 64
 */

void _mzd_trsm_upper_right_weird(mzd_t *U, mzd_t *B, const int cutoff);

void _mzd_trsm_upper_right_base(mzd_t *U, mzd_t *B, const int cutoff);

void mzd_trsm_upper_right(mzd_t *U, mzd_t *B, const int cutoff) {
  if(U->nrows != B->ncols)
    m4ri_die("mzd_trsm_upper_right: U nrows (%d) need to match B ncols (%d).\n", U->nrows.val(), B->ncols.val());
  if(U->nrows != U->ncols)
    m4ri_die("mzd_trsm_upper_right: U must be square and is found to be (%d) x (%d).\n", U->nrows.val(), U->ncols.val());
  
  _mzd_trsm_upper_right(U, B, cutoff);
}

void _mzd_trsm_upper_right(mzd_t *U, mzd_t *B, const int cutoff) {
  rci_t const nb = B->ncols;
  rci_t const mb = B->nrows;
  int const n1 = RADIX-B->offset;
  if(nb <= n1) {
    _mzd_trsm_upper_right_weird(U, B, cutoff);
    return;
  }
  /**
   \verbatim  
     _________ 
     \U00|   |
      \  |U01|
       \ |   |
        \|___|
         \U11|
          \  |
           \ |
            \|
      _______
     |B0 |B1 |
     |___|___|
   \endverbatim
     * \li U00 and B0 are possibly located at uneven locations.
   * \li Their column dimension is lower than 64.
   * \li The first column of U01, U11, B1 are aligned at words.
   */
  mzd_t *B0  = mzd_init_window (B,  0,  0, mb, (unsigned int)n1);	// FIXME, remove (unsigned int) casts.
  mzd_t *B1  = mzd_init_window (B,  0, (unsigned int)n1, mb, nb);
  mzd_t *U00 = mzd_init_window (U,  0,  0, (unsigned int)n1, (unsigned int)n1);
  mzd_t *U01 = mzd_init_window (U,  0, (unsigned int)n1, (unsigned int)n1, nb);
  mzd_t *U11 = mzd_init_window (U, (unsigned int)n1, (unsigned int)n1, nb, nb);
  
  _mzd_trsm_upper_right_weird (U00, B0, cutoff);
  mzd_addmul (B1, B0, U01, cutoff);
  _mzd_trsm_upper_right_even (U11, B1, cutoff);
  
  mzd_free_window(B0);
  mzd_free_window(B1);
  
  mzd_free_window(U00);
  mzd_free_window(U01);
  mzd_free_window(U11);
}

void _mzd_trsm_upper_right_weird(mzd_t *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  int const offset = B->offset;
  
  for(rci_t i = 1U; i < nb; ++i) {		// FIXME 1U --> 1
    
    /* Computes X_i = B_i + X_{0..i-1} U_{0..i-1,i} */
    
    register word ucol = 0;
    for(rci_t k = 0; k < i; ++k) {
      if(GET_BIT(U->rows[k][0U], i + U->offset))
	SET_BIT(ucol, k + offset);
    }
    /* doing 64 dotproducts at a time, to use the parity64 parallelism */
    rci_t giantstep;
    word tmp[64];
    for(giantstep = 0; giantstep + RADIX < mb; giantstep += RADIX) {
      for(int babystep = 0; babystep < RADIX; ++babystep)
	tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;

      word const dotprod = parity64 (tmp);
      
      for(int babystep = 0; babystep < RADIX; ++babystep)
        if(GET_BIT(dotprod, babystep))
          FLIP_BIT(B->rows[giantstep + babystep][0U], i + offset);
    }      
    for(int babystep = 0; giantstep + babystep < mb; ++babystep){
      tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
      
    }
    for(int babystep = (mb - giantstep).val(); babystep < 64; ++babystep){
      tmp[babystep] = 0;
    }
    
    word const dotprod = parity64 (tmp);
    
    for(int babystep = 0; giantstep + babystep < mb; ++babystep)
      if(GET_BIT(dotprod, babystep))
	FLIP_BIT(B->rows[giantstep + babystep][0U], i + offset);
  }
}

void _mzd_trsm_upper_right_even(mzd_t *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  if(nb <= TRSM_THRESHOLD) {
    /* base case */
    _mzd_trsm_upper_right_base (U, B, cutoff);
    return;
  }
  
  rci_t const nb1 = (((nb - 1) / RADIX + 1U) >> 1) * RADIX;
  
  mzd_t *B0 = mzd_init_window(B,  0,     0,   mb, nb1);
  mzd_t *B1 = mzd_init_window(B,  0,   nb1,   mb, nb);
  mzd_t *U00 = mzd_init_window(U, 0,     0, nb1, nb1);
  mzd_t *U01 = mzd_init_window(U, 0,   nb1, nb1, nb);
  mzd_t *U11 = mzd_init_window(U, nb1, nb1,  nb, nb);
  
  _mzd_trsm_upper_right_even (U00, B0, cutoff);
  mzd_addmul (B1, B0, U01, cutoff);
  _mzd_trsm_upper_right_even (U11, B1, cutoff);
  
  mzd_free_window(B0);
  mzd_free_window(B1);
  
  mzd_free_window(U00);
  mzd_free_window(U01);
  mzd_free_window(U11);
}

void _mzd_trsm_upper_right_base(mzd_t *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  for(rci_t i = 1U; i < nb; ++i) {		// FIXME: 1U --> 1
    /* Computes X_i = B_i + X_{0..i-1} U_{0..i-1,i} */
    register word ucol = 0;
    for(rci_t k = 0; k < i; ++k) {
      if(GET_BIT(U->rows[k][0U], i))
	SET_BIT(ucol, k);
    }
    
    /* doing 64 dotproducts at a time, to use the parity64 parallelism */
    rci_t giantstep;
    word tmp[64];
    for(giantstep = 0; giantstep + RADIX < mb; giantstep += RADIX) {
      for(int babystep = 0; babystep < RADIX; ++babystep)
	tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
      
      word const dotprod = parity64 (tmp);
      
      for(int babystep = 0; babystep < RADIX; ++babystep)
	if(GET_BIT(dotprod, babystep))
	  FLIP_BIT(B->rows[giantstep + babystep][0U], i);
    }
    
    for(int babystep = 0; giantstep + babystep < mb; ++babystep)
      tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
    for(int babystep = (mb - giantstep).val(); babystep < 64; ++babystep)
      tmp[babystep] = 0;
    
    word const dotprod = parity64(tmp);
    for(int babystep = 0; giantstep + babystep < mb; ++babystep)
      if(GET_BIT(dotprod, babystep))
	FLIP_BIT(B->rows[giantstep + babystep][0U], i);
  }
}


/*****************
 * LOWER RIGHT
 ****************/

/*
 * Variant where L and B start at an odd bit position Assumes that
 * L->ncols < 64
 */
void _mzd_trsm_lower_right_weird(mzd_t *L, mzd_t *B, const int cutoff);

/*
 * Variant where L and B start at an even bit position Assumes that
 * L->ncols < 64
 */

void _mzd_trsm_lower_right_even(mzd_t *L, mzd_t *B, const int cutoff);

void _mzd_trsm_lower_right_base(mzd_t *L, mzd_t *B, const int cutoff);

void mzd_trsm_lower_right(mzd_t *L, mzd_t *B, const int cutoff) {
  if(L->nrows != B->ncols)
    m4ri_die("mzd_trsm_lower_right: L nrows (%d) need to match B ncols (%d).\n", L->nrows.val(), B->ncols.val());
  if(L->nrows != L->ncols)
    m4ri_die("mzd_trsm_lower_right: L must be square and is found to be (%d) x (%d).\n", L->nrows.val(), L->ncols.val());
  
  _mzd_trsm_lower_right (L, B, cutoff);
}

void _mzd_trsm_lower_right(mzd_t *L, mzd_t *B, const int cutoff) {
  rci_t const nb = B->ncols;
  rci_t const mb = B->nrows;
  int const n1 = RADIX-B->offset;
  if(nb <= n1)
    _mzd_trsm_lower_right_weird (L, B, cutoff);
  else{
  /**
   \verbatim  
     |\
     | \  
     |  \
     |L00\
     |____\
     |    |\  
     |    | \
     |    |  \
     |L10 |L11\ 
     |____|____\
      _________
     |B0  |B1  |
     |____|____|
   \endverbatim
   * \li L00 and B0 are possibly located at uneven locations.
   * \li Their column dimension is lower than 64.
   * \li The first column of L10, L11, B1 are aligned to words.
   */
    mzd_t *B0  = mzd_init_window (B,  0,  0, mb, (unsigned int)n1);	// FIXME removed casts.
    mzd_t *B1  = mzd_init_window (B,  0, (unsigned int)n1, mb, nb);
    mzd_t *L00 = mzd_init_window (L,  0,  0, (unsigned int)n1, (unsigned int)n1);
    mzd_t *L10 = mzd_init_window (L,  (unsigned int)n1, 0, nb, (unsigned int)n1);
    mzd_t *L11 = mzd_init_window (L, (unsigned int)n1, (unsigned int)n1, nb, nb);
    
    _mzd_trsm_lower_right_even (L11, B1, cutoff);
    mzd_addmul (B0, B1, L10, cutoff);
    _mzd_trsm_lower_right_weird (L00, B0, cutoff);
    
    mzd_free_window(B0);
    mzd_free_window(B1);
    
    mzd_free_window(L00);
    mzd_free_window(L10);
    mzd_free_window(L11);
  }
}

void _mzd_trsm_lower_right_weird(mzd_t *L, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  int const offset = B->offset;
  
  for(rci_t i = nb - 1; i >= 0; --i) {
    
    /* Computes X_i = B_i + X_{i+1,n} L_{i+1..n,i} */
    
    register word ucol = 0;
    for(rci_t k = i + 1; k < nb; ++k) {
      if(GET_BIT(L->rows[k][0U], i + L->offset))
	SET_BIT(ucol, k + offset);
    }
    /* doing 64 dotproducts at a time, to use the parity64 parallelism */
    rci_t giantstep;
    word tmp[64];
    for(giantstep = 0; giantstep + RADIX < mb; giantstep += RADIX) {
      for(int babystep = 0; babystep < RADIX; ++babystep)
	tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
      
      word const dotprod = parity64(tmp);
      
      for(int babystep = 0; babystep < RADIX; ++babystep)
        if(GET_BIT(dotprod, babystep))
          FLIP_BIT(B->rows[giantstep + babystep][0U], i + offset);
    }      
    for(int babystep = 0; giantstep + babystep < mb; ++babystep){
      tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
    }
    for(int babystep = (mb - giantstep).val(); babystep < 64; ++babystep){
      tmp[babystep] = 0;
    }
    
    word const dotprod = parity64(tmp);
    
    for(int babystep = 0; giantstep + babystep < mb; ++babystep)
      if(GET_BIT(dotprod, babystep))
	FLIP_BIT(B->rows[giantstep + babystep ][0U], i + offset);
  }
}

void _mzd_trsm_lower_right_even(mzd_t *L, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  
  if(nb <= TRSM_THRESHOLD){
    /* base case */
    _mzd_trsm_lower_right_base (L, B, cutoff);
  }
  else {
    rci_t const nb1 = (((nb - 1) / RADIX + 1U) >> 1) * RADIX;

    mzd_t *B0 = mzd_init_window(B,  0,     0,   mb, nb1);
    mzd_t *B1 = mzd_init_window(B,  0,   nb1,   mb, nb);
    mzd_t *L00 = mzd_init_window(L, 0,     0, nb1, nb1);
    mzd_t *L10 = mzd_init_window(L, nb1, 0, nb, nb1);
    mzd_t *L11 = mzd_init_window(L, nb1, nb1,  nb, nb);

    _mzd_trsm_lower_right_even (L11, B1, cutoff);
    mzd_addmul (B0, B1, L10, cutoff);
    _mzd_trsm_lower_right_even (L00, B0, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window(L00);
    mzd_free_window(L10);
    mzd_free_window(L11);
  }
}

void _mzd_trsm_lower_right_base(mzd_t *L, mzd_t *B, const int cutoff) { 
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  for(rci_t i = nb - 1; i >= 0; --i) {
    
    /* Computes X_i = B_i + X_{i+1,n} L_{i+1..n,i} */
    
    register word ucol = 0;
    for(rci_t k = i + 1; k < nb; ++k) {
      if(GET_BIT(L->rows[k][0U], i))
	SET_BIT(ucol, k);
    }
    
    /* doing 64 dotproducts at a time, to use the parity64 parallelism */
    rci_t giantstep;
    word tmp[64];
    for(giantstep = 0; giantstep + RADIX < mb; giantstep += RADIX) {
      for(int babystep = 0; babystep < RADIX; ++babystep)
	tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
      
      word const dotprod = parity64 (tmp);
      
      for(int babystep = 0; babystep < RADIX; ++babystep)
	if(GET_BIT(dotprod, babystep))
	  FLIP_BIT(B->rows[giantstep + babystep][0U], i);
    }
    
    for(int babystep = 0; giantstep + babystep < mb; ++babystep)
      tmp[babystep] = B->rows[giantstep + babystep][0U] & ucol;
    for(int babystep = (mb - giantstep).val(); babystep < 64; ++babystep)
      tmp[babystep] = 0;
    
    word const dotprod = parity64 (tmp);
    for(int babystep = 0; giantstep + babystep < mb; ++babystep)
      if(GET_BIT(dotprod, babystep))
	FLIP_BIT(B->rows[giantstep + babystep][0U], i);
  }
}


/*****************
 * LOWER LEFT
 ****************/

/*
 * Variant where U and B start at an odd bit position. Assumes that
 * L->ncols < 64
 */

void _mzd_trsm_lower_left_weird(mzd_t *L, mzd_t *B, const int cutoff);

/* 
 * This version assumes that the matrices are at an even position on
 * the RADIX grid and that their dimension is a multiple of RADIX.
 */

void _mzd_trsm_lower_left_even(mzd_t *L, mzd_t *B, const int cutoff);

void mzd_trsm_lower_left(mzd_t *L, mzd_t *B, const int cutoff) {
  if(L->ncols != B->nrows)
    m4ri_die("mzd_trsm_lower_left: L ncols (%d) need to match B nrows (%d).\n", L->ncols.val(), B->nrows.val());
  if(L->nrows != L->ncols)
    m4ri_die("mzd_trsm_lower_left: L must be square and is found to be (%d) x (%d).\n", L->nrows.val(), L->ncols.val());
  
  _mzd_trsm_lower_left (L, B, cutoff);
}

void _mzd_trsm_lower_left(mzd_t *L, mzd_t *B, const int cutoff) {
  if(!L->offset)
    _mzd_trsm_lower_left_even(L, B, cutoff);
  else{
    rci_t const nb = B->ncols;
    rci_t const mb = B->nrows;
    int const m1 = RADIX - L->offset;
    if(mb <= m1) {
      _mzd_trsm_lower_left_weird (L, B, cutoff);
      return;
    }
    /**
    \verbatim  
    |\           ______
    | \         |      |
    |  \        |  B0  |
    |L00\       |      |
    |____\      |______|
    |    |\     |      |
    |    | \    |      |
    |    |  \   |  B1  |
    |L10 |L11\  |      |
    |____|____\ |______|
    \endverbatim 
    * \li L00 L10 B0 and B1 are possibly located at uneven locations.
    * \li Their column dimension is lower than 64.
    * \li The first column of L01, L11, B1 are aligned to words.
    */
      
    mzd_t *B0  = mzd_init_window (B,  0,  0, (unsigned int)m1, nb);
    mzd_t *B1  = mzd_init_window (B,  (unsigned int)m1, 0, mb, nb);
    mzd_t *L00 = mzd_init_window (L,  0,  0, (unsigned int)m1, (unsigned int)m1);
    mzd_t *L10 = mzd_init_window (L,  (unsigned int)m1, 0, mb, (unsigned int)m1);
    mzd_t *L11 = mzd_init_window (L, (unsigned int)m1, (unsigned int)m1, mb, mb);
    
    _mzd_trsm_lower_left_weird (L00, B0, cutoff);
    mzd_addmul (B1, L10, B0, cutoff);
    _mzd_trsm_lower_left_even (L11, B1, cutoff);
    
    mzd_free_window(B0);
    mzd_free_window(B1);
    
    mzd_free_window(L00);
    mzd_free_window(L10);
    mzd_free_window(L11);
  }
}

void _mzd_trsm_lower_left_weird(mzd_t *L, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  int const Boffset = B->offset;
  int const nbrest = (nb + Boffset) % RADIX;
  if(nb + B->offset > RADIX) {

    // Large B
    word const mask_begin = RIGHT_BITMASK(RADIX-B->offset);
    word const mask_end = LEFT_BITMASK(nbrest);

    // L[0,0] = 1, so no work required for i=0
    for(rci_t i = 1U; i < mb; ++i) {				// FIXME: 1U --> 1
      
      /* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
      /**
       * \todo needs to be optimized!
       **/ 
      wordPtr Lrow = L->rows[i];
      wordPtr Brow = B->rows[i];

      for(rci_t k = 0; k < i; ++k) {
	if(GET_BIT(Lrow[0U], k + L->offset)) {
	  Brow[0U] ^= B->rows[k][0U] & mask_begin;
	  for(wi_t j = 1U; j < B->width - 1U; ++j)
	    Brow[j] ^= B->rows[k][j];
          Brow[B->width - 1U] ^= B->rows[k][B->width - 1U] & mask_end;
	}
      }
    }
  } else { // Small B
    word const mask = MIDDLE_BITMASK(nb, B->offset);
    for(rci_t i = 1U; i < mb; ++i) {			// FIXME 1U --> 1
      /* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
      /**
       * \todo needs to be optimized!
       **/ 
      wordPtr Lrow = L->rows[i];
      wordPtr Brow = B->rows[i];

      for(rci_t k = 0; k < i; ++k) {
	if(GET_BIT(Lrow[0U], k + L->offset)) {
	  Brow[0U] ^= B->rows[k][0U] & mask;
	}
      }
    }
  }
}

void _mzd_trsm_lower_left_even(mzd_t *L, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  int const Boffset = B->offset;
  int const nbrest = (nb + Boffset) % RADIX;

  if(mb <= RADIX){
    /* base case */

    if(nb + B->offset > RADIX) {
      // B is large
      word const mask_begin = RIGHT_BITMASK(RADIX-B->offset);
      word const mask_end = LEFT_BITMASK(nbrest);

      for(rci_t i = 1U; i < mb; ++i) {				// FIXME: 1U --> 1
	/* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
        /**
         * \todo needs to be optimized!
         **/ 
	wordPtr Lrow = L->rows[i];
	wordPtr Brow = B->rows[i];

	for(rci_t k = 0; k < i; ++k) {
	  if(GET_BIT(Lrow[0U], k)){
	    Brow[0U] ^= B->rows[k][0U] & mask_begin;
	    for(wi_t j = 1U; j < B->width - 1U; ++j)
	      Brow[j] ^= B->rows[k][j];
	    Brow[B->width - 1U] ^= B->rows[k][B->width - 1U] & mask_end;
	  }
	}
      }
    } else { // B is small
      word const mask = MIDDLE_BITMASK(nb, B->offset);
      for(rci_t i = 1U; i < mb; ++i) {				// FIXME: 1U --> 1
	/* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
	/** Need to be optimized !!! **/
	wordPtr Lrow = L->rows [i];
	wordPtr Brow = B->rows [i];

	for(rci_t k = 0; k < i; ++k) {
	  if(GET_BIT(Lrow[0U], k)){
	    Brow[0U] ^= B->rows[k][0U] & mask;
	  }
	}
      }
    }
  } else {
    rci_t const mb1 = (((mb - 1) / RADIX + 1U) >> 1) * RADIX;

    mzd_t *B0 = mzd_init_window(B,  0,     0,   mb1, nb);
    mzd_t *B1 = mzd_init_window(B, mb1,    0,   mb,  nb);
    mzd_t *L00 = mzd_init_window(L, 0,     0, mb1, mb1);
    mzd_t *L10 = mzd_init_window(L, mb1,   0, mb, mb1);
    mzd_t *L11 = mzd_init_window(L, mb1, mb1, mb, mb);

    _mzd_trsm_lower_left_even (L00, B0, cutoff);

    mzd_addmul (B1, L10, B0, cutoff);

    _mzd_trsm_lower_left_even (L11, B1, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window(L00);
    mzd_free_window(L10);
    mzd_free_window(L11);
  }
}

/*****************
 * UPPER LEFT
 ****************/

/*
 * Variant where U and B start at an odd bit position
 * Assumes that U->ncols < 64
 */
void _mzd_trsm_upper_left_weird (mzd_t *U, mzd_t *B, const int cutoff);

void _mzd_trsm_upper_left_even(mzd_t *U, mzd_t *B, const int cutoff);

void mzd_trsm_upper_left(mzd_t *U, mzd_t *B, const int cutoff) {
  if(U->ncols != B->nrows)
    m4ri_die("mzd_trsm_upper_left: U ncols (%d) need to match B nrows (%d).\n", U->ncols.val(), B->nrows.val());
  if(U->nrows != U->ncols)
    m4ri_die("mzd_trsm_upper_left: U must be square and is found to be (%d) x (%d).\n", U->nrows.val(), U->ncols.val());
  
  _mzd_trsm_upper_left (U, B, cutoff);
}

void _mzd_trsm_upper_left(mzd_t *U, mzd_t *B, const int cutoff) {
  if(!U->offset)
    _mzd_trsm_upper_left_even (U, B, cutoff);
  else{
    rci_t const nb = B->ncols;
    rci_t const mb = B->nrows;
    int const m1 = RADIX - U->offset;
    if(mb <= m1) {
      _mzd_trsm_upper_left_weird (U, B, cutoff);
      return;
    }
    /**
     \verbatim
     __________   ______
     \ U00|    | |      |
      \   |U01 | |      |
       \  |    | |  B0  |
        \ |    | |      |
         \|____| |______|
          \    | |      |
           \U11| |      |
            \  | |  B1  |
             \ | |      |
              \| |______|
     \endverbatim 
     * \li U00, B0 and B1 are possibly located at uneven locations.
     * \li Their column dimension is greater than 64
     * \li The first column of U01, U11, B0 and B1 are aligned to words.
     */
    
    mzd_t *B0  = mzd_init_window (B,  0,  0, (unsigned int)m1, nb);
    mzd_t *B1  = mzd_init_window (B,  (unsigned int)m1, 0, mb, nb);
    mzd_t *U00 = mzd_init_window (U,  0,  0, (unsigned int)m1, (unsigned int)m1);
    mzd_t *U01 = mzd_init_window (U,  0, (unsigned int)m1, (unsigned int)m1, mb);
    mzd_t *U11 = mzd_init_window (U, (unsigned int)m1, (unsigned int)m1, mb, mb);
    
    _mzd_trsm_upper_left_even (U11, B1, cutoff);
    mzd_addmul (B0, U01, B1, cutoff);
    _mzd_trsm_upper_left_weird (U00, B0, cutoff);
    
    mzd_free_window(B0);
    mzd_free_window(B1);
    
    mzd_free_window(U00);
    mzd_free_window(U01);
    mzd_free_window(U11);
  }
}

void _mzd_trsm_upper_left_weird (mzd_t *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  int const Boffset = B->offset;
  int const nbrest = (nb + Boffset) % RADIX;
  if(nb + Boffset > RADIX) {

    // Large B
    word const mask_begin = RIGHT_BITMASK(RADIX-B->offset);
    word const mask_end = LEFT_BITMASK(nbrest);

    // U[mb-1,mb-1] = 1, so no work required for i=mb-1
    for(rci_t i = mb - 2; i >= 0; --i) {
      
      /* Computes X_i = B_i + U_{i,i+1..mb} X_{i+1..mb}  */
      wordPtr Urow = U->rows[i];
      wordPtr Brow = B->rows[i];

      for(rci_t k = i + 1; k < mb; ++k) {
	if(GET_BIT(Urow[0U], k + U->offset)) {
	  Brow[0U] ^= B->rows[k][0U] & mask_begin;
	  for(wi_t j = 1U; j < B->width - 1U; ++j)
	    Brow[j] ^= B->rows[k][j];
	  Brow[B->width - 1U] ^= B->rows[k][B->width - 1U] & mask_end;
	}
      }
    }
  } else { // Small B
    word const mask = MIDDLE_BITMASK(nb, B->offset);
    // U[mb-1,mb-1] = 1, so no work required for i=mb-1
    for(rci_t i = mb - 2; i >= 0; --i) {
      /* Computes X_i = B_i + U_{i,i+1..mb} X_{i+1..mb}  */
      wordPtr Urow = U->rows[i];
      wordPtr Brow = B->rows[i];

      for(rci_t k = i + 1; k < mb; ++k) {
	if(GET_BIT(Urow[0U], k + U->offset)) {
          Brow[0U] ^= B->rows[k][0U] & mask;
	}
      }
    }
  }
}

void _mzd_trsm_upper_left_even(mzd_t *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;
  int const Boffset = B->offset;
  int const nbrest = (nb + Boffset) % RADIX;

  if(mb <= RADIX) {
    /* base case */
    
    if(nb + B->offset > RADIX) {
      // B is large
      word const mask_begin = RIGHT_BITMASK(RADIX-B->offset);
      word const mask_end = LEFT_BITMASK(nbrest);
      
      // U[mb-1,mb-1] = 1, so no work required for i=mb-1
      for(rci_t i = mb - 2; i >= 0; --i) {
        
        /* Computes X_i = B_i + U_{i,i+1..mb} X_{i+1..mb}  */
        wordPtr Urow = U->rows[i];
        wordPtr Brow = B->rows[i];
        
        for(rci_t k = i + 1; k < mb; ++k) {
          if(GET_BIT(Urow[0U], k)){
            Brow[0U] ^= B->rows[k][0U] & mask_begin;
            for(wi_t j = 1U; j < B->width - 1U; ++j)
              Brow[j] ^= B->rows[k][j];
            Brow[B->width - 1U] ^= B->rows[k][B->width - 1U] & mask_end;
          }
        }
      }
    } else { // B is small
      word const mask = MIDDLE_BITMASK(nb, B->offset);
      // U[mb-1,mb-1] = 1, so no work required for i=mb-1
      for(rci_t i = mb - 2; i >= 0; --i) {
        
	/* Computes X_i = B_i + U_{i,i+1..mb} X_{i+1..mb}  */
	wordPtr Urow = U->rows [i];
	wordPtr Brow = B->rows [i];

	for(rci_t k = i + 1; k < mb; ++k) {
	  if(GET_BIT(Urow[0U], k)){
	    Brow[0U] ^= B->rows[k][0U] & mask;
	  }
	}
      }
    }
  } else if(mb <= MZD_MUL_BLOCKSIZE) {
    _mzd_trsm_upper_left_even_m4r(U, B, 0);
  } else {
    rci_t const mb1 = (((mb-1) / RADIX + 1U) >> 1) * RADIX;

    mzd_t *B0 = mzd_init_window(B,  0,     0,   mb1, nb);
    mzd_t *B1 = mzd_init_window(B, mb1,    0,   mb,  nb);
    mzd_t *U00 = mzd_init_window(U, 0,     0, mb1, mb1);
    mzd_t *U01 = mzd_init_window(U, 0,   mb1, mb1, mb);
    mzd_t *U11 = mzd_init_window(U, mb1, mb1, mb, mb);

    _mzd_trsm_upper_left_even (U11, B1, cutoff);

    _mzd_addmul (B0, U01, B1, cutoff);

    _mzd_trsm_upper_left_even (U00, B0, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window(U00);
    mzd_free_window(U01);
    mzd_free_window(U11);
  }
}
