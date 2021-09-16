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

#include <stdio.h>

#include "mzd.h"
#include "parity.h"
#include "strassen.h"
#include "triangular.h"
#include "triangular_russian.h"

/*****************
 * UPPER RIGHT
 ****************/

/*
 * Assumes that U->ncols < 64
 */
void _mzd_trsm_upper_right_base(mzd_t const *U, mzd_t *B);

void mzd_trsm_upper_right(mzd_t const *U, mzd_t *B, const int cutoff) {
  if (U->nrows != B->ncols)
    m4ri_die("mzd_trsm_upper_right: U nrows (%d) need to match B ncols (%d).\n", U->nrows,
             B->ncols);
  if (U->nrows != U->ncols)
    m4ri_die("mzd_trsm_upper_right: U must be square and is found to be (%d) x (%d).\n", U->nrows,
             U->ncols);

  _mzd_trsm_upper_right(U, B, cutoff);
}

void _mzd_trsm_upper_right_trtri(mzd_t const *U, mzd_t *B) {
  mzd_t *u = mzd_extract_u(NULL, U);
  mzd_trtri_upper(u);
  mzd_t *C = mzd_mul(NULL, B, u, 0);
  mzd_copy(B, C);
  mzd_free(C);
  mzd_free(u);
}

void _mzd_trsm_upper_right(mzd_t const *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  if (nb <= m4ri_radix) {
    /* base case */
    _mzd_trsm_upper_right_base(U, B);
    return;
  } else if (nb <= __M4RI_MUL_BLOCKSIZE) {
    _mzd_trsm_upper_right_trtri(U, B);
    return;
  }

  rci_t const nb1 = (((nb - 1) / m4ri_radix + 1) >> 1) * m4ri_radix;

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
   */

  mzd_t *B0        = mzd_init_window(B, 0, 0, mb, nb1);
  mzd_t *B1        = mzd_init_window(B, 0, nb1, mb, nb);
  mzd_t const *U00 = mzd_init_window_const(U, 0, 0, nb1, nb1);
  mzd_t const *U01 = mzd_init_window_const(U, 0, nb1, nb1, nb);
  mzd_t const *U11 = mzd_init_window_const(U, nb1, nb1, nb, nb);

  _mzd_trsm_upper_right(U00, B0, cutoff);
  mzd_addmul(B1, B0, U01, cutoff);
  _mzd_trsm_upper_right(U11, B1, cutoff);

  mzd_free_window(B0);
  mzd_free_window(B1);

  mzd_free_window((mzd_t *)U00);
  mzd_free_window((mzd_t *)U01);
  mzd_free_window((mzd_t *)U11);

  __M4RI_DD_MZD(B);
}

void _mzd_trsm_pack(word *tmp, mzd_t const *B, rci_t giantstep, word ucol)
{
#if 0
  for(int babystep = 0; babystep < m4ri_radix; ++babystep) 
    tmp[babystep] = mzd_row(B, giantstep + babystep)[0] & ucol;
#else
  assert(m4ri_radix == 64);
  tmp[ 0] = mzd_row_const(B, giantstep +  0)[0] & ucol; 
  tmp[ 1] = mzd_row_const(B, giantstep +  1)[0] & ucol;
  tmp[ 2] = mzd_row_const(B, giantstep +  2)[0] & ucol;
  tmp[ 3] = mzd_row_const(B, giantstep +  3)[0] & ucol;
  tmp[ 4] = mzd_row_const(B, giantstep +  4)[0] & ucol; 
  tmp[ 5] = mzd_row_const(B, giantstep +  5)[0] & ucol;
  tmp[ 6] = mzd_row_const(B, giantstep +  6)[0] & ucol;
  tmp[ 7] = mzd_row_const(B, giantstep +  7)[0] & ucol;
  tmp[ 8] = mzd_row_const(B, giantstep +  8)[0] & ucol;
  tmp[ 9] = mzd_row_const(B, giantstep +  9)[0] & ucol;
  tmp[10] = mzd_row_const(B, giantstep + 10)[0] & ucol;
  tmp[11] = mzd_row_const(B, giantstep + 11)[0] & ucol;
  tmp[12] = mzd_row_const(B, giantstep + 12)[0] & ucol;
  tmp[13] = mzd_row_const(B, giantstep + 13)[0] & ucol;
  tmp[14] = mzd_row_const(B, giantstep + 14)[0] & ucol;
  tmp[15] = mzd_row_const(B, giantstep + 15)[0] & ucol;
  tmp[16] = mzd_row_const(B, giantstep + 16)[0] & ucol;
  tmp[17] = mzd_row_const(B, giantstep + 17)[0] & ucol;
  tmp[18] = mzd_row_const(B, giantstep + 18)[0] & ucol;
  tmp[19] = mzd_row_const(B, giantstep + 19)[0] & ucol;
  tmp[20] = mzd_row_const(B, giantstep + 20)[0] & ucol;
  tmp[21] = mzd_row_const(B, giantstep + 21)[0] & ucol;
  tmp[22] = mzd_row_const(B, giantstep + 22)[0] & ucol;
  tmp[23] = mzd_row_const(B, giantstep + 23)[0] & ucol;
  tmp[24] = mzd_row_const(B, giantstep + 24)[0] & ucol;
  tmp[25] = mzd_row_const(B, giantstep + 25)[0] & ucol;
  tmp[26] = mzd_row_const(B, giantstep + 26)[0] & ucol;
  tmp[27] = mzd_row_const(B, giantstep + 27)[0] & ucol;
  tmp[28] = mzd_row_const(B, giantstep + 28)[0] & ucol;
  tmp[29] = mzd_row_const(B, giantstep + 29)[0] & ucol;
  tmp[30] = mzd_row_const(B, giantstep + 30)[0] & ucol;
  tmp[31] = mzd_row_const(B, giantstep + 31)[0] & ucol;
  tmp[32] = mzd_row_const(B, giantstep + 32)[0] & ucol;
  tmp[33] = mzd_row_const(B, giantstep + 33)[0] & ucol;
  tmp[34] = mzd_row_const(B, giantstep + 34)[0] & ucol;
  tmp[35] = mzd_row_const(B, giantstep + 35)[0] & ucol;
  tmp[36] = mzd_row_const(B, giantstep + 36)[0] & ucol;
  tmp[37] = mzd_row_const(B, giantstep + 37)[0] & ucol;
  tmp[38] = mzd_row_const(B, giantstep + 38)[0] & ucol;
  tmp[39] = mzd_row_const(B, giantstep + 39)[0] & ucol;
  tmp[40] = mzd_row_const(B, giantstep + 40)[0] & ucol;
  tmp[41] = mzd_row_const(B, giantstep + 41)[0] & ucol;
  tmp[42] = mzd_row_const(B, giantstep + 42)[0] & ucol;
  tmp[43] = mzd_row_const(B, giantstep + 43)[0] & ucol;
  tmp[44] = mzd_row_const(B, giantstep + 44)[0] & ucol;
  tmp[45] = mzd_row_const(B, giantstep + 45)[0] & ucol;
  tmp[46] = mzd_row_const(B, giantstep + 46)[0] & ucol;
  tmp[47] = mzd_row_const(B, giantstep + 47)[0] & ucol;
  tmp[48] = mzd_row_const(B, giantstep + 48)[0] & ucol;
  tmp[49] = mzd_row_const(B, giantstep + 49)[0] & ucol;
  tmp[50] = mzd_row_const(B, giantstep + 50)[0] & ucol;
  tmp[51] = mzd_row_const(B, giantstep + 51)[0] & ucol;
  tmp[52] = mzd_row_const(B, giantstep + 52)[0] & ucol;
  tmp[53] = mzd_row_const(B, giantstep + 53)[0] & ucol;
  tmp[54] = mzd_row_const(B, giantstep + 54)[0] & ucol;
  tmp[55] = mzd_row_const(B, giantstep + 55)[0] & ucol;
  tmp[56] = mzd_row_const(B, giantstep + 56)[0] & ucol;
  tmp[57] = mzd_row_const(B, giantstep + 57)[0] & ucol;
  tmp[58] = mzd_row_const(B, giantstep + 58)[0] & ucol;
  tmp[59] = mzd_row_const(B, giantstep + 59)[0] & ucol;
  tmp[60] = mzd_row_const(B, giantstep + 60)[0] & ucol;
  tmp[61] = mzd_row_const(B, giantstep + 61)[0] & ucol;
  tmp[62] = mzd_row_const(B, giantstep + 62)[0] & ucol;
  tmp[63] = mzd_row_const(B, giantstep + 63)[0] & ucol;
#endif
}

void _mzd_trsm_unpack(mzd_t *B, rci_t giantstep, word dotprod, rci_t i)
{
#if 0
  for(int babystep = 0; babystep < m4ri_radix; ++babystep)
    if(__M4RI_GET_BIT(dotprod, babystep)) 
      __M4RI_FLIP_BIT(mzd_row(B, giantstep + babystep)[0], i);
#else
  assert(m4ri_radix == 64);
  mzd_row(B, giantstep +  0)[0] ^= ((dotprod >>  0) & m4ri_one) << i;
  mzd_row(B, giantstep +  1)[0] ^= ((dotprod >>  1) & m4ri_one) << i;
  mzd_row(B, giantstep +  2)[0] ^= ((dotprod >>  2) & m4ri_one) << i;
  mzd_row(B, giantstep +  3)[0] ^= ((dotprod >>  3) & m4ri_one) << i;
  mzd_row(B, giantstep +  4)[0] ^= ((dotprod >>  4) & m4ri_one) << i;
  mzd_row(B, giantstep +  5)[0] ^= ((dotprod >>  5) & m4ri_one) << i;
  mzd_row(B, giantstep +  6)[0] ^= ((dotprod >>  6) & m4ri_one) << i;
  mzd_row(B, giantstep +  7)[0] ^= ((dotprod >>  7) & m4ri_one) << i;
  mzd_row(B, giantstep +  8)[0] ^= ((dotprod >>  8) & m4ri_one) << i;
  mzd_row(B, giantstep +  9)[0] ^= ((dotprod >>  9) & m4ri_one) << i;
  mzd_row(B, giantstep + 10)[0] ^= ((dotprod >> 10) & m4ri_one) << i;
  mzd_row(B, giantstep + 11)[0] ^= ((dotprod >> 11) & m4ri_one) << i;
  mzd_row(B, giantstep + 12)[0] ^= ((dotprod >> 12) & m4ri_one) << i;
  mzd_row(B, giantstep + 13)[0] ^= ((dotprod >> 13) & m4ri_one) << i;
  mzd_row(B, giantstep + 14)[0] ^= ((dotprod >> 14) & m4ri_one) << i;
  mzd_row(B, giantstep + 15)[0] ^= ((dotprod >> 15) & m4ri_one) << i;
  mzd_row(B, giantstep + 16)[0] ^= ((dotprod >> 16) & m4ri_one) << i;
  mzd_row(B, giantstep + 17)[0] ^= ((dotprod >> 17) & m4ri_one) << i;
  mzd_row(B, giantstep + 18)[0] ^= ((dotprod >> 18) & m4ri_one) << i;
  mzd_row(B, giantstep + 19)[0] ^= ((dotprod >> 19) & m4ri_one) << i;
  mzd_row(B, giantstep + 20)[0] ^= ((dotprod >> 20) & m4ri_one) << i;
  mzd_row(B, giantstep + 21)[0] ^= ((dotprod >> 21) & m4ri_one) << i;
  mzd_row(B, giantstep + 22)[0] ^= ((dotprod >> 22) & m4ri_one) << i;
  mzd_row(B, giantstep + 23)[0] ^= ((dotprod >> 23) & m4ri_one) << i;
  mzd_row(B, giantstep + 24)[0] ^= ((dotprod >> 24) & m4ri_one) << i;
  mzd_row(B, giantstep + 25)[0] ^= ((dotprod >> 25) & m4ri_one) << i;
  mzd_row(B, giantstep + 26)[0] ^= ((dotprod >> 26) & m4ri_one) << i;
  mzd_row(B, giantstep + 27)[0] ^= ((dotprod >> 27) & m4ri_one) << i;
  mzd_row(B, giantstep + 28)[0] ^= ((dotprod >> 28) & m4ri_one) << i;
  mzd_row(B, giantstep + 29)[0] ^= ((dotprod >> 29) & m4ri_one) << i;
  mzd_row(B, giantstep + 30)[0] ^= ((dotprod >> 30) & m4ri_one) << i;
  mzd_row(B, giantstep + 31)[0] ^= ((dotprod >> 31) & m4ri_one) << i;
  mzd_row(B, giantstep + 32)[0] ^= ((dotprod >> 32) & m4ri_one) << i;
  mzd_row(B, giantstep + 33)[0] ^= ((dotprod >> 33) & m4ri_one) << i;
  mzd_row(B, giantstep + 34)[0] ^= ((dotprod >> 34) & m4ri_one) << i;
  mzd_row(B, giantstep + 35)[0] ^= ((dotprod >> 35) & m4ri_one) << i;
  mzd_row(B, giantstep + 36)[0] ^= ((dotprod >> 36) & m4ri_one) << i;
  mzd_row(B, giantstep + 37)[0] ^= ((dotprod >> 37) & m4ri_one) << i;
  mzd_row(B, giantstep + 38)[0] ^= ((dotprod >> 38) & m4ri_one) << i;
  mzd_row(B, giantstep + 39)[0] ^= ((dotprod >> 39) & m4ri_one) << i;
  mzd_row(B, giantstep + 40)[0] ^= ((dotprod >> 40) & m4ri_one) << i;
  mzd_row(B, giantstep + 41)[0] ^= ((dotprod >> 41) & m4ri_one) << i;
  mzd_row(B, giantstep + 42)[0] ^= ((dotprod >> 42) & m4ri_one) << i;
  mzd_row(B, giantstep + 43)[0] ^= ((dotprod >> 43) & m4ri_one) << i;
  mzd_row(B, giantstep + 44)[0] ^= ((dotprod >> 44) & m4ri_one) << i;
  mzd_row(B, giantstep + 45)[0] ^= ((dotprod >> 45) & m4ri_one) << i;
  mzd_row(B, giantstep + 46)[0] ^= ((dotprod >> 46) & m4ri_one) << i;
  mzd_row(B, giantstep + 47)[0] ^= ((dotprod >> 47) & m4ri_one) << i;
  mzd_row(B, giantstep + 48)[0] ^= ((dotprod >> 48) & m4ri_one) << i;
  mzd_row(B, giantstep + 49)[0] ^= ((dotprod >> 49) & m4ri_one) << i;
  mzd_row(B, giantstep + 50)[0] ^= ((dotprod >> 50) & m4ri_one) << i;
  mzd_row(B, giantstep + 51)[0] ^= ((dotprod >> 51) & m4ri_one) << i;
  mzd_row(B, giantstep + 52)[0] ^= ((dotprod >> 52) & m4ri_one) << i;
  mzd_row(B, giantstep + 53)[0] ^= ((dotprod >> 53) & m4ri_one) << i;
  mzd_row(B, giantstep + 54)[0] ^= ((dotprod >> 54) & m4ri_one) << i;
  mzd_row(B, giantstep + 55)[0] ^= ((dotprod >> 55) & m4ri_one) << i;
  mzd_row(B, giantstep + 56)[0] ^= ((dotprod >> 56) & m4ri_one) << i;
  mzd_row(B, giantstep + 57)[0] ^= ((dotprod >> 57) & m4ri_one) << i;
  mzd_row(B, giantstep + 58)[0] ^= ((dotprod >> 58) & m4ri_one) << i;
  mzd_row(B, giantstep + 59)[0] ^= ((dotprod >> 59) & m4ri_one) << i;
  mzd_row(B, giantstep + 60)[0] ^= ((dotprod >> 60) & m4ri_one) << i;
  mzd_row(B, giantstep + 61)[0] ^= ((dotprod >> 61) & m4ri_one) << i;
  mzd_row(B, giantstep + 62)[0] ^= ((dotprod >> 62) & m4ri_one) << i;
  mzd_row(B, giantstep + 63)[0] ^= ((dotprod >> 63) & m4ri_one) << i;
#endif  
}

void _mzd_trsm_upper_right_base(mzd_t const *U, mzd_t *B) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  for (rci_t i = 1; i < nb; ++i) {
    /* Computes X_i = B_i + X_{0..i-1} U_{0..i-1,i} */
    register word ucol = 0;
    for (rci_t k = 0; k < i; ++k) {
      if (__M4RI_GET_BIT(mzd_row_const(U, k)[0], i)) __M4RI_SET_BIT(ucol, k);
    }

    /* doing 64 dotproducts at a time, to use the m4ri_parity64 parallelism */
    rci_t giantstep;
    word tmp[64];

    for (giantstep = 0; giantstep + m4ri_radix < mb; giantstep += m4ri_radix) {
      _mzd_trsm_pack(tmp, B, giantstep, ucol);
      word dotprod = m4ri_parity64(tmp);
      _mzd_trsm_unpack(B, giantstep, dotprod, i);
    }

    for (int babystep = 0; giantstep + babystep < mb; ++babystep)
      tmp[babystep] = mzd_row(B, giantstep + babystep)[0] & ucol;
    for (int babystep = mb - giantstep; babystep < 64; ++babystep) tmp[babystep] = 0;

    word const dotprod = m4ri_parity64(tmp);
    for (int babystep = 0; giantstep + babystep < mb; ++babystep)
      if (__M4RI_GET_BIT(dotprod, babystep)) __M4RI_FLIP_BIT(mzd_row(B, giantstep + babystep)[0], i);
  }

  __M4RI_DD_MZD(B);
}

/*****************
 * LOWER RIGHT
 ****************/

void _mzd_trsm_lower_right_base(mzd_t const *L, mzd_t *B);

void mzd_trsm_lower_right(mzd_t const *L, mzd_t *B, const int cutoff) {
  if (L->nrows != B->ncols)
    m4ri_die("mzd_trsm_lower_right: L nrows (%d) need to match B ncols (%d).\n", L->nrows,
             B->ncols);
  if (L->nrows != L->ncols)
    m4ri_die("mzd_trsm_lower_right: L must be square and is found to be (%d) x (%d).\n", L->nrows,
             L->ncols);

  _mzd_trsm_lower_right(L, B, cutoff);
}

void _mzd_trsm_lower_right(mzd_t const *L, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  if (nb <= m4ri_radix) {
    _mzd_trsm_lower_right_base(L, B);
    return;
  }

  rci_t const nb1 = (((nb - 1) / m4ri_radix + 1) >> 1) * m4ri_radix;

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
   */

  mzd_t *B0        = mzd_init_window(B, 0, 0, mb, nb1);
  mzd_t *B1        = mzd_init_window(B, 0, nb1, mb, nb);
  mzd_t const *L00 = mzd_init_window_const(L, 0, 0, nb1, nb1);
  mzd_t const *L10 = mzd_init_window_const(L, nb1, 0, nb, nb1);
  mzd_t const *L11 = mzd_init_window_const(L, nb1, nb1, nb, nb);

  _mzd_trsm_lower_right(L11, B1, cutoff);
  mzd_addmul(B0, B1, L10, cutoff);
  _mzd_trsm_lower_right(L00, B0, cutoff);

  mzd_free_window(B0);
  mzd_free_window(B1);

  mzd_free_window((mzd_t *)L00);
  mzd_free_window((mzd_t *)L10);
  mzd_free_window((mzd_t *)L11);

  __M4RI_DD_MZD(B);
}

void _mzd_trsm_lower_right_base(mzd_t const *L, mzd_t *B) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  for (rci_t i = nb - 1; i >= 0; --i) {
    /* Computes X_i = B_i + X_{i+1,n} L_{i+1..n,i} */
    register word ucol = 0;
    for (rci_t k = i + 1; k < nb; ++k) {
      if (__M4RI_GET_BIT(mzd_row_const(L, k)[0], i)) __M4RI_SET_BIT(ucol, k);
    }

    /* doing 64 dotproducts at a time, to use the parity64 parallelism */
    rci_t giantstep;
    word tmp[64];
    for (giantstep = 0; giantstep + m4ri_radix < mb; giantstep += m4ri_radix) {
      _mzd_trsm_pack(tmp, B, giantstep, ucol);
      word dotprod = m4ri_parity64(tmp);
      _mzd_trsm_unpack(B, giantstep, dotprod, i);
    }
    for (int babystep = 0; giantstep + babystep < mb; ++babystep)
      tmp[babystep] = mzd_row(B, giantstep + babystep)[0] & ucol;
    for (int babystep = mb - giantstep; babystep < 64; ++babystep) tmp[babystep] = 0;

    word const dotprod = m4ri_parity64(tmp);
    for (int babystep = 0; giantstep + babystep < mb; ++babystep)
      if (__M4RI_GET_BIT(dotprod, babystep)) __M4RI_FLIP_BIT(mzd_row(B, giantstep + babystep)[0], i);
  }

  __M4RI_DD_MZD(B);
}

/*****************
 * LOWER LEFT
 ****************/

void mzd_trsm_lower_left(mzd_t const *L, mzd_t *B, const int cutoff) {
  if (L->ncols != B->nrows)
    m4ri_die("mzd_trsm_lower_left: L ncols (%d) need to match B nrows (%d).\n", L->ncols, B->nrows);
  if (L->nrows != L->ncols)
    m4ri_die("mzd_trsm_lower_left: L must be square and is found to be (%d) x (%d).\n", L->nrows,
             L->ncols);

  _mzd_trsm_lower_left(L, B, cutoff);
}

void _mzd_trsm_lower_left(mzd_t const *L, mzd_t *B, const int cutoff) {
  rci_t const mb   = B->nrows;
  rci_t const nb   = B->ncols;
  int const nbrest = nb % m4ri_radix;

  if (mb <= m4ri_radix) {
    /* base case */
    word const mask_end = __M4RI_LEFT_BITMASK(nbrest);
    for (rci_t i = 1; i < mb; ++i) {
      /* Computes X_i = B_i + L_{i,0..i-1} X_{0..i-1}  */
      word const *Lrow = mzd_row_const(L, i);
      word *Brow = mzd_row(B, i);

      for (rci_t k = 0; k < i; ++k) {
        if (__M4RI_GET_BIT(Lrow[0], k)) {
          for (wi_t j = 0; j < B->width - 1; ++j) Brow[j] ^= mzd_row_const(B, k)[j];
          Brow[B->width - 1] ^= mzd_row_const(B, k)[B->width - 1] & mask_end;
        }
      }
    }
  } else if (mb <= __M4RI_MUL_BLOCKSIZE) {
    _mzd_trsm_lower_left_russian(L, B, 0);
  } else {
    rci_t const mb1 = (((mb - 1) / m4ri_radix + 1) >> 1) * m4ri_radix;

    mzd_t *B0        = mzd_init_window(B, 0, 0, mb1, nb);
    mzd_t *B1        = mzd_init_window(B, mb1, 0, mb, nb);
    mzd_t const *L00 = mzd_init_window_const(L, 0, 0, mb1, mb1);
    mzd_t const *L10 = mzd_init_window_const(L, mb1, 0, mb, mb1);
    mzd_t const *L11 = mzd_init_window_const(L, mb1, mb1, mb, mb);

    _mzd_trsm_lower_left(L00, B0, cutoff);

    mzd_addmul(B1, L10, B0, cutoff);

    _mzd_trsm_lower_left(L11, B1, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window((mzd_t *)L00);
    mzd_free_window((mzd_t *)L10);
    mzd_free_window((mzd_t *)L11);
  }
  __M4RI_DD_MZD(B);
}

/*****************
 * UPPER LEFT
 ****************/

void mzd_trsm_upper_left(mzd_t const *U, mzd_t *B, const int cutoff) {
  if (U->ncols != B->nrows)
    m4ri_die("mzd_trsm_upper_left: U ncols (%d) need to match B nrows (%d).\n", U->ncols, B->nrows);
  if (U->nrows != U->ncols)
    m4ri_die("mzd_trsm_upper_left: U must be square and is found to be (%d) x (%d).\n", U->nrows,
             U->ncols);

  _mzd_trsm_upper_left(U, B, cutoff);
}

void _mzd_trsm_upper_left(mzd_t const *U, mzd_t *B, const int cutoff) {
  rci_t const mb = B->nrows;
  rci_t const nb = B->ncols;

  if (mb <= m4ri_radix) {
    /* base case */

    word const mask_end = B->high_bitmask;

    // U[mb-1,mb-1] = 1, so no work required for i=mb-1
    for (rci_t i = mb - 2; i >= 0; --i) {

      /* Computes X_i = B_i + U_{i,i+1..mb} X_{i+1..mb}  */
      word const *Urow = mzd_row_const(U, i);
      word *Brow = mzd_row(B, i);

      for (rci_t k = i + 1; k < mb; ++k) {
        if (__M4RI_GET_BIT(Urow[0], k)) {
          for (wi_t j = 0; j < B->width - 1; ++j) Brow[j] ^= mzd_row_const(B, k)[j];
          Brow[B->width - 1] ^= mzd_row_const(B, k)[B->width - 1] & mask_end;
        }
      }
    }
  } else if (mb <= __M4RI_MUL_BLOCKSIZE) {
    _mzd_trsm_upper_left_russian(U, B, 0);
  } else {
    rci_t const mb1 = (((mb - 1) / m4ri_radix + 1) >> 1) * m4ri_radix;

    mzd_t *B0        = mzd_init_window(B, 0, 0, mb1, nb);
    mzd_t *B1        = mzd_init_window(B, mb1, 0, mb, nb);
    mzd_t const *U00 = mzd_init_window_const(U, 0, 0, mb1, mb1);
    mzd_t const *U01 = mzd_init_window_const(U, 0, mb1, mb1, mb);
    mzd_t const *U11 = mzd_init_window_const(U, mb1, mb1, mb, mb);

    _mzd_trsm_upper_left(U11, B1, cutoff);

    _mzd_addmul(B0, U01, B1, cutoff);

    _mzd_trsm_upper_left(U00, B0, cutoff);

    mzd_free_window(B0);
    mzd_free_window(B1);

    mzd_free_window((mzd_t *)U00);
    mzd_free_window((mzd_t *)U01);
    mzd_free_window((mzd_t *)U11);
  }

  __M4RI_DD_MZD(B);
}

mzd_t *mzd_trtri_upper(mzd_t *U) {
  if (U->nrows * U->ncols < __M4RI_CPU_L3_CACHE << 1) {
    mzd_trtri_upper_russian(U, 0);
  } else {
    rci_t const n = U->nrows;
    rci_t n2      = (((n - 1) / m4ri_radix + 1) >> 1);

#if __M4RI_HAVE_SSE2
    if (n2 % 2) n2 += 1;
#endif
    n2 *= m4ri_radix;

    assert(n2 < n);

    mzd_t *U00 = mzd_init_window(U, 0, 0, n2, n2);
    mzd_t *U01 = mzd_init_window(U, 0, n2, n2, n);
    mzd_t *U11 = mzd_init_window(U, n2, n2, n, n);

    _mzd_trsm_upper_left(U00, U01, 0);
    _mzd_trsm_upper_right(U11, U01, 0);
    mzd_trtri_upper(U00);
    mzd_trtri_upper(U11);

    mzd_free_window((mzd_t *)U00);
    mzd_free_window((mzd_t *)U01);
    mzd_free_window((mzd_t *)U11);
  }
  return U;
}
