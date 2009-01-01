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

#include "misc.h"
#include "packedmatrix.h"
#include "trsm.h"
#include "parity.h"
#include <stdio.h>
#include "pluq_mmpf.h"
#include "strassen.h"
#include "lqup.h"

size_t mzd_pluq (packedmatrix *A, permutation * P, permutation * Q, const int cutoff) {
  if (P->length != A->nrows)
    m4ri_die("mzd_pluq: Permutation P length (%d) must match A nrows (%d)\n",P->length, A->nrows);
  if (Q->length != A->ncols)
    m4ri_die("mzd_pluq: Permutation Q length (%d) must match A ncols (%d)\n",Q->length, A->ncols);
  return _mzd_pluq(A, P, Q, cutoff);
}

size_t _mzd_pluq(packedmatrix *A, permutation * P, permutation * Q, const int cutoff) {
  size_t nrows = A->nrows;
  size_t ncols = A->ncols;

  if (ncols <= PLUQ_CUTOFF){
    /* Base case */
    return _mzd_pluq_mmpf(A, P, Q, 0);
  } else{
    /* Block divide and conquer algorithm */
    
    /*                     n1
     *   ------------------------------------------
     *   | A0              | A1                   |
     *   |                 |                      |
     *   |                 |                      |
     *   |                 |                      |
     *   ------------------------------------------
     */

    size_t i, j;
    size_t n1 = (((ncols - 1) / RADIX + 1) >> 1) * RADIX;
    packedmatrix *A0  = mzd_init_window(A,  0,  0, nrows,    n1);
    packedmatrix *A1  = mzd_init_window(A,  0, n1, nrows, ncols);

    size_t r1, r2;
    /* First recursive call */
    r1 = _mzd_pluq(A0, P, Q, cutoff);

    /*           r1           n1
     *   ------------------------------------------
     *   | A00    |           | A01               |
     *   |        |           |                   |
     * r1------------------------------------------ 
     * 
     *   | A01    |           | A11               |
     *   |        |           |                   |
     *   ------------------------------------------
     */

    packedmatrix *A00 = mzd_init_window(A,  0,  0, r1, r1);
    packedmatrix *A10 = mzd_init_window(A, r1,  0, nrows, r1);
    packedmatrix *A01 = mzd_init_window(A,  0, n1, r1, ncols);
    packedmatrix *A11 = mzd_init_window(A, r1, n1, nrows, ncols);

    if (r1) {
      /* Computation of the Schur complement */
      mzd_apply_p_left(A1, P);
      _mzd_trsm_lower_left(A00, A01, cutoff);
      mzd_addmul(A11, A10, A01, cutoff);
    }

    /* Second recursive call */
    permutation *P2 = mzp_init_window(P, r1, nrows);
    permutation *Q2 = mzp_init_window(Q, n1, ncols);

    r2 = _mzd_pluq(A11, P2, Q2, cutoff);

    /* Update A10 */
    mzd_apply_p_left(A10, P2);

    /* Update P */
    for (i = 0; i < nrows - r1; ++i)
      P2->values[i] += r1;

    
    /* Permute the Right side of U to the left */

    /* // fast, buggy solution
     * // Update A01 
     * mzd_apply_p_right(A01, Q2);
     * // Update Q
     * for (i = 0; i < ncols - n1; ++i)
     *  Q2->values[i] += n1;
     * permutation * Q2b = mzp_init_window(Q, r1, ncols);
     * packedmatrix* A01b = mzd_init_window (A, 0, r1, r1, ncols);
     * packedmatrix* A11b = mzd_init_window (A, r1, r1, nrows, ncols);
     * mzd_col_block_rotate (A01b, 0, n1-r1, n1 - r1 + r2, 1);
     * mzd_col_block_rotate (A11b, 0, n1-r1, n1 - r1 + r2, 1);
     * // Update Q
     * mzp_free_window(Q2b);
     * mzd_free_window(A01b);
     * mzd_free_window(A11b);
     * for(i=n1, j=r1; i<n1+r2; i++, j++) {
     *   Q->values[j] = Q->values[i];
     *   mzd_col_swap(A, j, i);
     * }
     */

    // easy working solution
    /* undo permutation */
    mzd_apply_p_right_trans(A11, Q2);


    permutation *tmp = mzp_init(A->ncols);
    for(i=0, j=n1; j<n1+r2; i++, j++) {
      //mzd_col_swap(A, r1 + i, n1 + Q2->values[i]);
      tmp->values[r1+i] = Q2->values[i] + n1;
      Q->values[r1+i] = Q2->values[i] + n1;
    }
    for(i=r1+r2; i<ncols; i++) {
      Q->values[i] = i;
    }
    mzd_apply_p_right(A, tmp);
    mzp_free(tmp);

    mzp_free_window(Q2);
    mzp_free_window(P2);

    mzd_free_window(A0);
    mzd_free_window(A1);
    mzd_free_window(A00);
    mzd_free_window(A01);
    mzd_free_window(A10);
    mzd_free_window(A11);

    return r1 + r2;
  }
}


size_t _mzd_pluq_naive(packedmatrix *A, permutation *P, permutation *Q)  {
  size_t i, j, l, curr_pos;
  int found;

  curr_pos = 0;

  for (curr_pos = 0; curr_pos < A->ncols; ) {
    found = 0;
    /* search for some pivot */
    for (j = curr_pos; j < A->ncols; j++) {
      for (i = curr_pos; i< A->nrows; i++ ) {
	if (mzd_read_bit(A, i, j)) {
          found = 1;
          break;
        }
      }
      if(found)
        break;
    }
    
    if(found) {
      P->values[curr_pos] = i;
      Q->values[curr_pos] = j;
      mzd_row_swap(A, curr_pos, i);
      mzd_col_swap(A, curr_pos, j);
          
      /* clear below but preserve transformation matrix */
      if (curr_pos +1 < A->ncols){
	for(l=curr_pos+1; l<A->nrows; l++) {
	  if (mzd_read_bit(A, l, curr_pos)) {
	    mzd_row_add_offset(A, l, curr_pos, curr_pos+1);
	  }
	}
      }
      curr_pos++;
    } else {
      break;
    }
  }
  for (i=curr_pos; i<A->nrows; ++i)
    P->values[i]=i;
  for (i=curr_pos; i<A->ncols; ++i)
    Q->values[i]=i;
  return curr_pos;
}
 
size_t mzd_echelonize_pluq(packedmatrix *A, int full) {
  permutation *P = mzp_init(A->nrows);
  permutation *Q = mzp_init(A->ncols);

  size_t r = mzd_pluq(A, P, Q, 0);

  if(full) {
    packedmatrix *U = mzd_init_window(A, 0, 0, r, r);
    packedmatrix *B = mzd_init_window(A, 0, r, r, A->ncols);
    if(r!=A->ncols) 
      mzd_trsm_upper_left(U, B, 0);
    if(r!=0) 
      mzd_set_ui(U, 0);
    size_t i;
    for(i=0; i<r; i++) {
      mzd_write_bit(A, i, i, 1);
    }
    mzd_free_window(U);
    mzd_free_window(B);
  } else {
    size_t i, j;
    for(i=0; i<r; i++) {
      for(j=0; j<i; i+=RADIX) {
        const size_t length = MIN(RADIX, i-j);
        mzd_clear_bits(A, i, j, length);
      }
    }
  }
  
  if(r!=A->nrows) {
    packedmatrix *R = mzd_init_window(A, r, 0, A->nrows, A->ncols);
    mzd_set_ui(R, 0);
    mzd_free_window(R);
  }

  mzd_apply_p_right_trans(A, Q);
  mzp_free(P);
  mzp_free(Q);
  return r;
}
