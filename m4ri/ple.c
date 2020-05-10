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

#include "mzd.h"
#include "parity.h"
#include "ple.h"
#include "ple_russian.h"
#include "strassen.h"
#include "triangular.h"
#include <stdio.h>

rci_t mzd_ple(mzd_t *A, mzp_t *P, mzp_t *Q, int const cutoff) {
  if (P->length != A->nrows)
    m4ri_die("mzd_ple: Permutation P length (%d) must match A nrows (%d)\n", P->length, A->nrows);
  if (Q->length != A->ncols)
    m4ri_die("mzd_ple: Permutation Q length (%d) must match A ncols (%d)\n", Q->length, A->ncols);
  return _mzd_ple(A, P, Q, cutoff);
}

rci_t mzd_pluq(mzd_t *A, mzp_t *P, mzp_t *Q, int const cutoff) {
  if (P->length != A->nrows)
    m4ri_die("mzd_pluq: Permutation P length (%d) must match A nrows (%d)\n", P->length, A->nrows);
  if (Q->length != A->ncols)
    m4ri_die("mzd_pluq: Permutation Q length (%d) must match A ncols (%d)\n", Q->length, A->ncols);
  rci_t r = _mzd_pluq(A, P, Q, cutoff);
  return r;
}

rci_t _mzd_pluq(mzd_t *A, mzp_t *P, mzp_t *Q, int const cutoff) {
  rci_t r = _mzd_ple(A, P, Q, cutoff);
  if (r && r < A->nrows) {
    mzd_t *A0 = mzd_init_window(A, 0, 0, r, A->ncols);
    mzd_apply_p_right_trans_tri(A0, Q);
    mzd_free_window(A0);
  } else {
    mzd_apply_p_right_trans_tri(A, Q);
  }
  return r;
}

rci_t _mzd_ple(mzd_t *A, mzp_t *P, mzp_t *Q, int const cutoff) {
  rci_t ncols = A->ncols;

#if 1
  rci_t nrows = mzd_first_zero_row(A);
  for (rci_t i = nrows; i < A->nrows; ++i) P->values[i] = i;
  for (rci_t i = 0; i < A->ncols; ++i) Q->values[i] = i;
  if (!nrows) { return 0; }
#else
  rci_t nrows = A->nrows;
#endif

  if (ncols <= m4ri_radix || A->width * A->nrows <= __M4RI_PLE_CUTOFF) {
    /* this improves data locality and runtime considerably */
    mzd_t *Abar = mzd_copy(NULL, A);
    rci_t r     = _mzd_ple_russian(Abar, P, Q, 0);
    mzd_copy(A, Abar);
    mzd_free(Abar);
    return r;
  }

  {
    /* Block divide and conquer algorithm */

    /*                     n1
     *   ------------------------------------------
     *   | A0              | A1                   |
     *   |                 |                      |
     *   |                 |                      |
     *   |                 |                      |
     *   ------------------------------------------
     */

    rci_t n1 = (((ncols - 1) / m4ri_radix + 1) >> 1) * m4ri_radix;

    mzd_t *A0 = mzd_init_window(A, 0, 0, nrows, n1);
    mzd_t *A1 = mzd_init_window(A, 0, n1, nrows, ncols);

    /* First recursive call */

    mzp_t *P1 = mzp_init_window(P, 0, nrows);
    mzp_t *Q1 = mzp_init_window(Q, 0, A0->ncols);
    rci_t r1  = _mzd_ple(A0, P1, Q1, cutoff);

    /*           r1           n1
     *   ------------------------------------------
     *   | A00    |           | A01               |
     *   |        |           |                   |
     * r1------------------------------------------
     *
     *   | A10    |           | A11               |
     *   |        |           |                   |
     *   ------------------------------------------
     */

    mzd_t *A00 = mzd_init_window(A, 0, 0, r1, r1);
    mzd_t *A10 = mzd_init_window(A, r1, 0, nrows, r1);
    mzd_t *A01 = mzd_init_window(A, 0, n1, r1, ncols);
    mzd_t *A11 = mzd_init_window(A, r1, n1, nrows, ncols);

    if (r1) {
      /* Computation of the Schur complement */
      mzd_apply_p_left(A1, P1);
      _mzd_trsm_lower_left(A00, A01, cutoff);
      mzd_addmul(A11, A10, A01, cutoff);
    }
    mzp_free_window(P1);
    mzp_free_window(Q1);

    /* Second recursive call */
    mzp_t *P2 = mzp_init_window(P, r1, nrows);
    mzp_t *Q2 = mzp_init_window(Q, n1, ncols);

    rci_t r2 = _mzd_ple(A11, P2, Q2, cutoff);

    /*           n
     *   -------------------
     *   |      A0b        |
     *   |                 |
     *   r1-----------------
     *   |      A1b        |
     *   |                 |
     *   -------------------
     */

    /* Update A10 */
    mzd_apply_p_left(A10, P2);

    /* Update P */
    for (rci_t i = 0; i < nrows - r1; ++i) P2->values[i] += r1;

    // Update the A0b block (permutation + rotation)
    for (rci_t i = 0, j = n1; j < ncols; ++i, ++j) Q2->values[i] += n1;

    for (rci_t i = n1, j = r1; i < n1 + r2; ++i, ++j) Q->values[j] = Q->values[i];

    /* Compressing L */

    _mzd_compress_l(A, r1, n1, r2);

    mzp_free_window(Q2);
    mzp_free_window(P2);

    mzd_free_window(A0);
    mzd_free_window(A1);
    mzd_free_window(A00);
    mzd_free_window(A01);
    mzd_free_window(A10);
    mzd_free_window(A11);

    __M4RI_DD_MZD(A);
    __M4RI_DD_MZP(P);
    __M4RI_DD_MZP(Q);
    __M4RI_DD_RCI(r1 + r2);
    return r1 + r2;
  }
}

rci_t _mzd_pluq_naive(mzd_t *A, mzp_t *P, mzp_t *Q) {
  rci_t curr_pos = 0;
  for (curr_pos = 0; curr_pos < A->ncols;) {
    int found = 0;
    /* search for some pivot */
    rci_t i, j;
    for (j = curr_pos; j < A->ncols; ++j) {
      for (i = curr_pos; i < A->nrows; ++i) {
        if (mzd_read_bit(A, i, j)) {
          found = 1;
          break;
        }
      }
      if (found) break;
    }

    if (found) {
      P->values[curr_pos] = i;
      Q->values[curr_pos] = j;
      mzd_row_swap(A, curr_pos, i);
      mzd_col_swap(A, curr_pos, j);

      /* clear below but preserve transformation matrix */
      if (curr_pos + 1 < A->ncols) {
        for (rci_t l = curr_pos + 1; l < A->nrows; ++l) {
          if (mzd_read_bit(A, l, curr_pos)) { mzd_row_add_offset(A, l, curr_pos, curr_pos + 1); }
        }
      }
      ++curr_pos;
    } else {
      break;
    }
  }
  for (rci_t i = curr_pos; i < A->nrows; ++i) P->values[i] = i;
  for (rci_t i = curr_pos; i < A->ncols; ++i) Q->values[i] = i;

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
  __M4RI_DD_MZP(Q);
  __M4RI_DD_RCI(curr_pos);
  return curr_pos;
}

rci_t _mzd_ple_naive(mzd_t *A, mzp_t *P, mzp_t *Q) {
  rci_t col_pos = 0;
  rci_t row_pos = 0;
  /* search for some pivot */
  while (row_pos < A->nrows && col_pos < A->ncols) {
    int found = 0;
    rci_t i, j;
    for (j = col_pos; j < A->ncols; ++j) {
      for (i = row_pos; i < A->nrows; ++i) {
        if (mzd_read_bit(A, i, j)) {
          found = 1;
          break;
        }
      }
      if (found) break;
    }
    if (found) {
      P->values[row_pos] = i;
      Q->values[row_pos] = j;
      mzd_row_swap(A, row_pos, i);
      // mzd_col_swap(A, curr_pos, j);

      /* clear below but preserve transformation matrix */
      if (j + 1 < A->ncols) {
        for (rci_t l = row_pos + 1; l < A->nrows; ++l) {
          if (mzd_read_bit(A, l, j)) { mzd_row_add_offset(A, l, row_pos, j + 1); }
        }
      }
      ++row_pos;
      col_pos = j + 1;
    } else {
      break;
    }
  }
  for (rci_t i = row_pos; i < A->nrows; ++i) P->values[i] = i;
  for (rci_t i = row_pos; i < A->ncols; ++i) Q->values[i] = i;

  /* Now compressing L */
  for (rci_t j = 0; j < row_pos; ++j) {
    if (Q->values[j] > j) {
      // To be optimized by a copy_row function
      mzd_col_swap_in_rows(A, Q->values[j], j, j, A->nrows);
    }
  }

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
  __M4RI_DD_MZP(Q);
  __M4RI_DD_RCI(row_pos);
  return row_pos;
}
