/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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

#include "echelonform.h"
#include "brilliantrussian.h"
#include "ple.h"
#include "triangular.h"

rci_t mzd_echelonize(mzd_t *A, int full) {
  return _mzd_echelonize_m4ri(A, full, 0, 1, __M4RI_ECHELONFORM_CROSSOVER_DENSITY);
}

rci_t mzd_echelonize_m4ri(mzd_t *A, int full, int k) {
  return _mzd_echelonize_m4ri(A, full, k, 0, 1.0);
}

rci_t mzd_echelonize_pluq(mzd_t *A, int full) {
  mzp_t *P = mzp_init(A->nrows);
  mzp_t *Q = mzp_init(A->ncols);

  rci_t r;
  if(full) {
#if 0
    r = mzd_pluq(A, P, Q, 0);
    mzd_t *U = mzd_init_window(A, 0, 0, r, r);
    mzd_t *B = mzd_init_window(A, 0, r, r, A->ncols);
    if(r!=A->ncols) 
      mzd_trsm_upper_left(U, B, 0);
    if(r) 
      mzd_set_ui(U, 0);
    for(rci_t i = 0; i < r; ++i)
      mzd_write_bit(A, i, i, 1);
    mzd_free_window(U);
    mzd_free_window(B);

    if(r) {
      mzd_t *A0 = mzd_init_window(A, 0, 0, r, A->ncols);
      mzd_apply_p_right(A0, Q);
      mzd_free_window(A0);
    } else {
      mzd_apply_p_right(A, Q);
    }
#else
    r = mzd_pluq(A, P, Q, 0);

    mzd_t *U = mzd_init_window(A, 0, 0, r, r);
    const rci_t r_radix = m4ri_radix*(r/m4ri_radix);

    if(r_radix == r && r!=A->ncols) {

      mzd_t *B = mzd_init_window(A, 0, r, r, A->ncols);
      if(r!=A->ncols) 
        mzd_trsm_upper_left(U, B, 0);
      mzd_free_window(B);

    } else if (r_radix != r && r!=A->ncols) {
      assert(r_radix < r);

      if(A->ncols > r_radix+m4ri_radix) {
        mzd_t *B0  = mzd_submatrix(NULL, A, 0, r_radix, r, r_radix+m4ri_radix);
        mzd_t *B0w = mzd_init_window(    A, 0, r_radix, r, r_radix+m4ri_radix);
        mzd_t *B1  = mzd_init_window(A, 0, r_radix+m4ri_radix, r, A->ncols);

        mzd_trsm_upper_left(U, B0, 0);
        mzd_trsm_upper_left(U, B1, 0);

        mzd_copy(B0w, B0);
        mzd_free(B0);
        mzd_free_window(B0w);
        mzd_free_window(B1);

      } else {
        mzd_t *B = mzd_submatrix(NULL, A, 0, r_radix, r, A->ncols);
        mzd_t *Bw = mzd_init_window(A, 0, r_radix, r, A->ncols);

        mzd_trsm_upper_left(U, B, 0);

        mzd_copy(Bw, B);
        mzd_free_window(Bw);
        mzd_free(B);     
      }
    }

    mzd_set_ui(U, 1);
    mzd_free_window(U);
    
    if(r) {
      mzd_t *A0 = mzd_init_window(A, 0, 0, r, A->ncols);
      mzd_apply_p_right(A0, Q);
      mzd_free_window(A0);
    }
#endif

  } else {
    r = mzd_ple(A, P, Q, 0);
    for(rci_t i = 0; i < r; ++i) {
      for(rci_t j = 0; j <= i; j += m4ri_radix) {
        int const length = MIN(m4ri_radix, i - j + 1);
        mzd_clear_bits(A, i, j, length);
      }
      mzd_write_bit(A, i, Q->values[i], 1);
    }
  }

  if(r != A->nrows) {
    mzd_t *R = mzd_init_window(A, r, 0, A->nrows, A->ncols);
    mzd_set_ui(R, 0);
    mzd_free_window(R);
  }

  mzp_free(P);
  mzp_free(Q);

  __M4RI_DD_MZD(A);
  __M4RI_DD_RCI(r);
  return r;
}
