 /*******************************************************************
 *
 *            M4RI: Method of the Four Russians Inversion
 *
 *       Copyright (C) 2008 Clement Pernet <clement.pernet@gmail.com>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
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

#include "trsm.h"
#include "strassen.h"
#include "packedmatrix.h"
#include "misc.h"
#include "parity.h"
#include "stdio.h"

void mzd_trsm_upper_right (packedmatrix *U, packedmatrix *B, const int cutoff) {
  if(U->nrows != B->ncols)
    m4ri_die("mzd_trsm_upper_right: U nrows (%d) need to match B ncols (%d).\n", U->nrows, B->ncols);
  if(U->nrows != U->ncols)
    m4ri_die("mzd_trsm_upper_right: U must be square and is found to be (%d) x (%d).\n", U->nrows, U->ncols);
  
  if (cutoff <= 0)
    m4ri_die("mzd_trsm_upper_right: cutoff must be > 0.\n");

  _mzd_trsm_upper_right (U, B, cutoff);
}

void _mzd_trsm_upper_right (packedmatrix *U, packedmatrix *B, const int cutoff) {
  size_t nb = B->ncols;
  size_t mb = B->nrows;
  size_t n1 = RADIX-B->offset;
  if (nb <= n1)
    _mzd_trsm_upper_right_weird (U, B, cutoff);
  else{
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
  
   * U00 and B0 are possibly located at uneven locations.
   * Their column dimension is lower than 64
   * The first column of U01, U11, B1 are aligned to words.
   */
    packedmatrix *B0  = mzd_init_window (B,  0,  0, mb, n1);
    packedmatrix *B1  = mzd_init_window (B,  0, n1, mb, nb);
    packedmatrix *U00 = mzd_init_window (U,  0,  0, n1, n1);
    packedmatrix *U01 = mzd_init_window (U,  0, n1, n1, nb);
    packedmatrix *U11 = mzd_init_window (U, n1, n1, nb, nb);
    
    _mzd_trsm_upper_right_weird (U00, B0, cutoff);
    mzd_addmul (B1, B0, U01, cutoff);
    _mzd_trsm_upper_right_even (U11, B1, cutoff);
    
    mzd_free_window(B0);
    mzd_free_window(B1);
    
    mzd_free_window(U00);
    mzd_free_window(U01);
    mzd_free_window(U11);
  }
}

/**
 * Variant where U and B start at an odd bit position
 * Assumes that U->ncols < 64
 */
void _mzd_trsm_upper_right_weird (packedmatrix *U, packedmatrix *B, const int cutoff) {

  size_t mb = B->nrows;
  size_t nb = B->ncols;
  size_t offset = B->offset;
  
  for (size_t i=1; i < nb; ++i) {
    
    /* Computes X_i = B_i + X_{0..i-1} U_{0..i-1,i} */
    
    register word ucol = 0;
    for (size_t k=0; k<i; ++k) {
      if (GET_BIT (U->values[U->rowswap[k]], i + U->offset))
	SET_BIT (ucol, k+offset);
    }
    /* doing 64 dotproducts at a time, to use the parity64 parallelism */
    size_t giantstep;
    word tmp[64];
    for (giantstep = 0; giantstep + RADIX < mb; giantstep += RADIX) {
      for (size_t babystep = 0; babystep < RADIX; ++babystep)
	tmp [babystep] = B->values [B->rowswap [babystep + giantstep]] & ucol;

      word dotprod = parity64 (tmp);
      
      for (size_t babystep = 0; babystep < RADIX; ++babystep)
	  if (GET_BIT (dotprod, babystep))
	    FLIP_BIT (B->values [B->rowswap [giantstep + babystep]], i + offset);
      }      
    for (size_t babystep = 0; babystep < mb - giantstep; ++babystep){
      tmp [babystep] = B->values [B->rowswap [babystep + giantstep]] & ucol;

    }
    for (size_t babystep = mb-giantstep; babystep < 64; ++babystep){
      tmp [babystep] = 0;
    }

    word dotprod = parity64 (tmp);
    
    for (size_t babystep = 0; babystep < mb - giantstep; ++babystep)
      if (GET_BIT (dotprod, babystep))
	FLIP_BIT (B->values [B->rowswap [giantstep + babystep ]], i + offset);
  }
}

/**
 * Variant where U and B start at an even bit position
 * Assumes that U->ncols < 64
 */
void _mzd_trsm_upper_right_even (packedmatrix *U, packedmatrix *B, const int cutoff) {

  size_t mb = B->nrows;
  size_t nb = B->ncols;
    
  if (nb <= RADIX){
    /* base case */
    for (size_t i=1; i < nb; ++i) {

      /* Computes X_i = B_i + X_{0..i-1} U_{0..i-1,i} */

      register word ucol = 0;
      for (size_t k=0; k<i; ++k) {
	if (GET_BIT (U->values[U->rowswap[k]], i))
	  SET_BIT (ucol, k);
      }

      /* doing 64 dotproducts at a time, to use the parity64 parallelism */
      size_t giantstep;
      word tmp[64];
      for (giantstep = 0; giantstep + RADIX < mb; giantstep += RADIX) {
	for (size_t babystep = 0; babystep < RADIX; ++babystep)
	  tmp [babystep] = B->values [B->rowswap [babystep + giantstep]] & ucol;

	word dotprod = parity64 (tmp);

	for (size_t babystep = 0; babystep < RADIX; ++babystep)
	  if (GET_BIT (dotprod, babystep))
	    FLIP_BIT (B->values [B->rowswap [giantstep + babystep]], i);
      }  
      
      for (size_t babystep = 0; babystep < mb - giantstep; ++babystep)
	tmp [babystep] = B->values [B->rowswap [babystep + giantstep]] & ucol;
      for (size_t babystep = mb-giantstep; babystep < 64; ++babystep)
	tmp [babystep] = 0;

      word dotprod = parity64 (tmp);
      for (size_t babystep = 0; babystep < mb - giantstep; ++babystep)
	if (GET_BIT (dotprod, babystep))
	  FLIP_BIT (B->values [B->rowswap [giantstep + babystep ]], i);
    }
  } else {
    size_t nb1 = (((nb-1) / RADIX + 1) >> 1) * RADIX;

    packedmatrix *B0 = mzd_init_window(B,  0,     0,   mb, nb1);
    packedmatrix *B1 = mzd_init_window(B,  0,   nb1,   mb, nb);
    packedmatrix *U00 = mzd_init_window(U, 0,     0, nb1, nb1);
    packedmatrix *U01 = mzd_init_window(U, 0,   nb1, nb1, nb);
    packedmatrix *U11 = mzd_init_window(U, nb1, nb1,  nb, nb);

    _mzd_trsm_upper_right_even (U00, B0, cutoff);
    mzd_addmul (B1, B0, U01, cutoff);
    _mzd_trsm_upper_right_even (U11, B1, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window(U00);
    mzd_free_window(U01);
    mzd_free_window(U11);
  }
}

/* Lower Left Implementations */

void mzd_trsm_lower_left (packedmatrix *L, packedmatrix *B, const int cutoff) {
  if(L->ncols != B->nrows)
    m4ri_die("mzd_trsm_lower_left: L ncols (%d) need to match B nrows (%d).\n", L->ncols, B->nrows);
  if(L->nrows != L->ncols)
    m4ri_die("mzd_trsm_lower_left: L must be square and is found to be (%d) x (%d).\n", L->nrows, L->ncols);
  
  if (cutoff <= 0)
    m4ri_die("mzd_trsm_lower_left: cutoff must be > 0.\n");

  _mzd_trsm_lower_left (L, B, cutoff);
}

void _mzd_trsm_lower_left (packedmatrix *L, packedmatrix *B, const int cutoff) {


  if (!L->offset)
    _mzd_trsm_lower_left_even (L, B, cutoff);
  else{
    size_t nb = B->ncols;
    size_t mb = B->nrows;
    size_t m1 = RADIX - L->offset;
    if (mb <= m1)
      _mzd_trsm_lower_left_weird (L, B, cutoff);
    else{
      /**
       *  
       * |\           ______
       * | \         |      |
       * |  \        |  B0  |
       * |L00\       |      |
       * |____\      |______|
       * |    |\     |      |
       * |    | \    |      |
       * |    |  \   |  B1  |
       * |L10 |L11\  |      |
       * |____|____\ |______|
       * 
       * L00 L10 and B0 are possibly located at uneven locations.
       * Their column dimension is lower than 64
       * The first column of L01, L11, B1 are aligned to words.
       */
      
      packedmatrix *B0  = mzd_init_window (B,  0,  0, m1, nb);
      packedmatrix *B1  = mzd_init_window (B,  m1, 0, mb, nb);
      packedmatrix *L00 = mzd_init_window (L,  0,  0, m1, m1);
      packedmatrix *L10 = mzd_init_window (L,  m1, 0, mb, m1);
      packedmatrix *L11 = mzd_init_window (L, m1, m1, mb, mb);
      
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
}

/**
 * Variant where L and B start at an odd bit position
 * Assumes that L->ncols < 64
 */
void _mzd_trsm_lower_left_weird (packedmatrix *L, packedmatrix *B, const int cutoff) {

  size_t mb = B->nrows;
  size_t nb = B->ncols;
  size_t Boffset = B->offset;
  size_t nbrest = (nb + Boffset) % RADIX;
  if (nb + B->offset >= RADIX) {

    // Large B
    word mask_begin = RIGHT_BITMASK(RADIX-B->offset);
    word mask_end = LEFT_BITMASK(nbrest);

    for (size_t i=1; i < mb; ++i) {
      
      /* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
      // Need to be optimized !!!
      size_t Lidx = L->rowswap[i];
      size_t Bidx = B->rowswap[i];

      for (size_t k=0; k<i; ++k) {
	if (GET_BIT (L->values [Lidx], k + L->offset)){
	  B->values [Bidx] ^= B->values [B->rowswap [k]] & mask_begin;
	  for (size_t j = 1; j < B->width-1; ++j)
	    B->values [Bidx + j] ^= B->values [B->rowswap [k] + j];
	  B->values [Bidx + B->width - 1] ^= B->values [B->rowswap [k] + B->width - 1] & mask_end;
	}
      }
    }
  } else { // Small B

    word mask = ((ONE << nb) - 1) ;
    mask <<= (RADIX-nb-B->offset);

    for (size_t i=1; i < mb; ++i) {
      /* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
      // Need to be optimized !!!
      size_t Lidx = L->rowswap[i];
      size_t Bidx = B->rowswap[i];

      for (size_t k=0; k<i; ++k) {
	if (GET_BIT (L->values [Lidx], k + L->offset)){
	  B->values [Bidx] ^= B->values [B->rowswap [k]] & mask;
	}
      }
    }
  }
}

/**
 * Variant where L and B start at an odd bit position
 * Assumes that L->ncols < 64
 */
void _mzd_trsm_lower_left_even (packedmatrix *L, packedmatrix *B, const int cutoff) {

  size_t mb = B->nrows;
  size_t nb = B->ncols;
  size_t Boffset = B->offset;
  size_t nbrest = (nb + Boffset) % RADIX;

  if (mb <= RADIX){
    /* base case */

    if (nb + B->offset >= RADIX) {
      // B is large
      word mask_begin = RIGHT_BITMASK(RADIX-B->offset);
      word mask_end = LEFT_BITMASK(nbrest);

      for (size_t i=1; i < mb; ++i) {
	/* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
	// Need to be optimized !!!
	size_t Lidx = L->rowswap[i];
	size_t Bidx = B->rowswap[i];

	for (size_t k=0; k<i; ++k) {
	  if (GET_BIT (L->values [Lidx], k)){
	    B->values [Bidx] ^= B->values [B->rowswap [k]] & mask_begin;
	    for (size_t j = 1; j < B->width-1; ++j)
	      B->values [Bidx + j] ^= B->values [B->rowswap [k] + j];
	    B->values [Bidx + B->width - 1] ^= B->values [B->rowswap [k] + B->width - 1] & mask_end;
	  }
	}
      }
    } else { // B is small
      word mask = ((ONE << nb) - 1) ;
      mask <<= (RADIX-nb-B->offset);
      for (size_t i=1; i < mb; ++i) {
	/* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
	// Need to be optimized !!!
	size_t Lidx = L->rowswap [i];
	size_t Bidx = B->rowswap [i];

	for (size_t k=0; k<i; ++k) {
	  if (GET_BIT (L->values [Lidx], k)){
	    B->values [Bidx] ^= B->values [B->rowswap[k]] & mask;
	  }
	}
      }
    }
  } else {
    size_t mb1 = (((mb-1) / RADIX + 1) >> 1) * RADIX;

    packedmatrix *B0 = mzd_init_window(B,  0,     0,   mb1, nb);
    packedmatrix *B1 = mzd_init_window(B, mb1,    0,   mb,  nb);
    packedmatrix *L00 = mzd_init_window(L, 0,     0, mb1, mb1);
    packedmatrix *L10 = mzd_init_window(L, mb1,   0, mb, mb1);
    packedmatrix *L11 = mzd_init_window(L, mb1, mb1, mb, mb);

    _mzd_trsm_lower_left_even (L00, B0, cutoff);

    _mzd_addmul (B1, L10, B0, cutoff);

    _mzd_trsm_lower_left_even (L11, B1, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window(L00);
    mzd_free_window(L10);
    mzd_free_window(L11);
  }
}
