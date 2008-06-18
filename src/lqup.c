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

#include "lqup.h"
#include "strassen.h"
#include "trsm.h"
#include "packedmatrix.h"
#include "misc.h"
#include "parity.h"
#include "stdio.h"

void mzd_lqup (packedmatrix *A, permutation * P, permutation * Q, const int cutoff) {

  if (cutoff <= 0)
    m4ri_die("mzd_trsm_upper_right: cutoff must be > 0.\n");

  _mzd_trsm_upper_right (U, B, cutoff);
}

size_t _mzd_lqup (packedmatrix *A, permutation * P, permutation * Q, const int cutoff) {

  size_t nrows = A->nrows;
  size_t ncols = A->ncols;

  if (ncols < ){
    /**
     * Base case
     */
    
  } else{
    /**
     * Block divide and conquer algorithm
     */
    
    size_t n1 = (((ncols - 1) / RADIX + 1) >> 1) * RADIX;
    packedmatrix *A0  = mzd_init_window_weird (A,  0,  0, nrows,    n1, A->offset);
    packedmatrix *A1  = mzd_init_window_weird (A,  0, n1, nrows, ncols, A->offset);

    size_t r1, r2;

    r1 = _mzd_lqup (A0, P, Q, cutoff);

    if (r1) {
      _mzd_apply_p_left_notrans (A1, P, 0, r1);
      _mzd_trsm_lower_left (U0, A01, cutoff);
      _mzd_addmul (A11, A10, A01, cutoff);
    }

    r2 = _mzd_lqup (A22, P, Q, cutoff);
    
    
  }

  packedmatrix *B0  = mzd_init_window_weird (B,  0,  0, mb, n1, B->offset);
  packedmatrix *B1  = mzd_init_window_weird (B,  0, n1, mb, nb,         0);
  packedmatrix *U00 = mzd_init_window_weird (U,  0,  0, n1, n1, U->offset);
  packedmatrix *U01 = mzd_init_window_weird (U,  0, n1, n1, nb,         0);
  packedmatrix *U11 = mzd_init_window_weird (U, n1, n1, nb, nb,         0);
  
  _mzd_trsm_upper_right_weird (U00, B0, cutoff);
  _mzd_addmul (B1, B0, U01, cutoff);
  _mzd_trsm_upper_right_even (U11, B1, cutoff);
  
  mzd_free_window(B0);
  mzd_free_window(B1);
  
  mzd_free_window(U00);
  mzd_free_window(U01);
  mzd_free_window(U11);
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
      if (GET_BIT (U->values[U->rowswap[k]], i + offset))
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
      
      for (size_t babystep = 0; babystep < mb - giantstep; ++babystep)
	tmp [babystep] = B->values [B->rowswap [babystep + giantstep]] & ucol;
      for (size_t babystep = mb-giantstep; babystep < 64; ++babystep)
	tmp [babystep] = 0;

      word dotprod = parity64 (tmp);
      for (size_t babystep = 0; babystep < mb - giantstep; ++babystep)
	if (GET_BIT (dotprod, babystep))
	  FLIP_BIT (B->values [B->rowswap [giantstep + babystep ]], i + offset);
    }
}

/**
 * Variant where U and B start at an odd bit position
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
    _mzd_addmul (B1, B0, U01, cutoff);
    _mzd_trsm_upper_right_even (U11, B1, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window(U00);
    mzd_free_window(U01);
    mzd_free_window(U11);
  }
}
