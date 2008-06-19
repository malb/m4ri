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
    m4ri_die("mzd_lqup: cutoff must be > 0.\n");
  if (P->length != A->nrows)
    m4ri_die("mzd_lqup: Permutation P length (%d) must match A nrows (%d)\n",P->length, A->nrows);
  if (Q->length != A->ncols)
    m4ri_die("mzd_lqup: Permutation Q length (%d) must match A ncols (%d)\n",Q->length, A->ncols);
  _mzd_lqup (A, P, Q, cutoff);
}

size_t _mzd_lqup (packedmatrix *A, permutation * P, permutation * Q, const int cutoff) {
  size_t nrows = A->nrows;
  size_t ncols = A->ncols;

  if (ncols < ){
    /* Base case */
    mzd_reduce_m4ri (A, 0, 0, NULL, NULL);
    
  } else{
    /* Block divide and conquer algorithm */
    
    size_t n1 = (((ncols - 1) / RADIX + 1) >> 1) * RADIX;
    packedmatrix *A0  = mzd_init_window_weird (A,  0,  0, nrows,    n1, A->offset);
    packedmatrix *A1  = mzd_init_window_weird (A,  0, n1, nrows, ncols, A->offset);

    size_t r1, r2;
    /* First recursive call */
    r1 = _mzd_lqup (A0, P, Q, cutoff);

    packedmatrix *A10  = mzd_init_window_weird (A,  r1, 0, nrows, r1, A->offset);
    packedmatrix *A01  = mzd_init_window_weird (A,  0, n1, r1, ncols, A->offset);
    if (r1) {
      /* Computation of the Schur complement */
      _mzd_apply_p_left_notrans (A1, P, 0, r1);
      _mzd_trsm_lower_left (U0, A01, cutoff);
      _mzd_addmul (A11, A10, A01, cutoff);
    }

    /* Second recursive call */
    packedmatrix *A11  = mzd_init_window_weird (A,  r1, n1, nrows, ncols, A->offset);
    permutation * P2 = mzd_init_permutation_window (P, r1, nrows);
    permutation * Q2 = mzd_init_permutation_window (Q, n1, ncols);
    r2 = _mzd_lqup (A11, P2, Q2, cutoff);

    /* Update Q */
    for (size_t i = 0; i < ncols - n1; ++i)
      Q2->values += n1;

    /* Update A10 */
    mzd_apply_p_left_trans (A10, P);

    /* Update P */
    for (size_t i = 0; i < nrows - r1; ++i)
      P2->values += r1;
    
    /* Permute the Right side of U to the left */
    permutation * Q2b = mzd_init_permutation_window (P, r1, ncols);
    // TODO : fix A->offset
    packedmatrix* A11b = mzd_init_window_weird (A, r1, r1, nrows, ncols, A->offset);
    packedmatrix* A01b = mzd_init_window_weird (A, 0, r1, r1, ncols, A->offset);
    mzd_col_block_rotate (A1b, 0, n1, n1 + r2, 1, Q2b);
    mzd_col_block_rotate (A0b, 0, n1, n1 + r2, 0, Q2b);
    
    mzd_free_permutation_window (Q2);
    mzd_free_permutation_window (P2);
    mzd_free_window(A1b);
    mzd_free_window(A0);
    mzd_free_window(A1);
    mzd_free_window(A01);
    mzd_free_window(A10);
    mzd_free_window(A11);
  }
  
}
