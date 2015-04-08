 /*******************************************************************
 *
 *            M4RI: Linear Algebra over GF(2)
 *
 *       Copyright (C) 2008 Jean-Guillaume.Dumas@imag.fr
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

#include "solve.h"
#include "strassen.h"
#include "ple.h"
#include "triangular.h"
#include "mzp.h"

int mzd_solve_left(mzd_t *A, mzd_t *B, int const cutoff, int const inconsistency_check) {    
  if(A->ncols > B->nrows)
    m4ri_die("mzd_solve_left: A ncols (%d) must be smaller than B nrows (%d).\n", A->ncols, B->nrows);

  return _mzd_solve_left(A, B, cutoff, inconsistency_check);
}
 
int mzd_pluq_solve_left (mzd_t const *A, rci_t rank, 
                         mzp_t const *P, mzp_t const *Q, 
                         mzd_t *B, int const cutoff, int const inconsistency_check) {
  if(A->ncols > B->nrows)
    m4ri_die("mzd_pluq_solve_left: A ncols (%d) need to be lower than B nrows (%d).\n", A->ncols, B->nrows);
  if(P->length != A->nrows)
      m4ri_die("mzd_pluq_solve_left: A nrows (%d) need to match P size (%d).\n", A->nrows, P->length);
  if(Q->length != A->ncols)
      m4ri_die("mzd_pluq_solve_left: A ncols (%d) need to match Q size (%d).\n", A->ncols, P->length);

  return _mzd_pluq_solve_left (A, rank, P, Q, B, cutoff, inconsistency_check);
}

int _mzd_pluq_solve_left(mzd_t const *A, rci_t rank, 
                         mzp_t const *P, mzp_t const *Q, 
                         mzd_t *B, int const cutoff, int const inconsistency_check) {
  /** A is supposed to store L lower triangular and U upper triangular
   *  B is modified in place 
   *  (Bi's in the comments are just modified versions of B)
   *  PLUQ = A
   *  1) P B2 = B1
   *  2) L B3 = B2
   *  3) U B4 = B3
   *  4) Q B5 = B4
   */

  int retval = 0;

  /* P B2 = B1 or B2 = P^T B1 */
  mzd_apply_p_left(B, P);
  
  /* L B3 = B2 */
  
  /* view on the upper part of L */
  mzd_t const *LU = mzd_init_window_const(A, 0, 0, rank, rank);
  mzd_t *Y1 = mzd_init_window(B, 0, 0, rank, B->ncols);
  mzd_trsm_lower_left(LU, Y1, cutoff);

  if (inconsistency_check) { /* Check for inconsistency */    
    /** FASTER without this check; update with the lower part of L
     */
    mzd_t const *H  = mzd_init_window_const(A, rank, 0, A->nrows, rank);
    mzd_t *Y2 = mzd_init_window(B, rank, 0, A->nrows, B->ncols);
    if(A->nrows < B->nrows) {
      mzd_t *Y3 = mzd_init_window(B, A->nrows, 0, B->nrows, B->ncols);
      mzd_set_ui(Y3, 0);
      mzd_free_window(Y3);
    }
    mzd_addmul(Y2, H, Y1, cutoff);
    /*
     * test whether Y2 is the zero matrix
     */
    if(!mzd_is_zero(Y2)) {
      retval = -1;
    }
    mzd_free_window((mzd_t*)H);
    mzd_free_window(Y2);
  }
  /* U B4 = B3 */
  mzd_trsm_upper_left(LU, Y1, cutoff);
  mzd_free_window((mzd_t*)LU);
  mzd_free_window(Y1);
  
  if (!inconsistency_check) {
    /** Default is to set the undefined bits to zero if inconsistency
     * has been checked then Y2 bits are already all zeroes thus this
     * clearing is not needed
     */
    for(rci_t i = rank; i < B->nrows; ++i) {
      for(rci_t j = 0; j < B->ncols; j += m4ri_radix) {
        mzd_clear_bits(B, i, j, MIN(m4ri_radix, B->ncols - j));
      }
    }
  }
  /* Q B5 = B4 or B5 = Q^T B4 */
  mzd_apply_p_left_trans(B, Q);

  /* P L U Q B5 = B1 */
  __M4RI_DD_MZD(B); 
  __M4RI_DD_INT(retval);
  return retval;
}

int _mzd_solve_left(mzd_t *A, mzd_t *B, int const cutoff, int const inconsistency_check) {
  /**
   *  B is modified in place 
   *  (Bi's in the comments are just modified versions of B)
   *  1) PLUQ = A
   *  2) P B2 = B1
   *  3) L B3 = B2
   *  4) U B4 = B3
   *  5) Q B5 = B4
   */
  mzp_t *P = mzp_init(A->nrows);
  mzp_t *Q = mzp_init(A->ncols);
  
  /* PLUQ = A */
  rci_t rank = _mzd_pluq(A, P, Q, cutoff);
  /* 2, 3, 4, 5 */
  int retval = mzd_pluq_solve_left(A, rank, P, Q, B, cutoff, inconsistency_check);
  
  mzp_free(P);
  mzp_free(Q);

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZD(B);
  return retval;
}

mzd_t *mzd_kernel_left_pluq(mzd_t *A, int const cutoff) {
  mzp_t *P = mzp_init(A->nrows);
  mzp_t *Q = mzp_init(A->ncols);

  rci_t r = mzd_pluq(A, P, Q, cutoff);

  if (r == A->ncols) {
    mzp_free(P);
    mzp_free(Q);
    __M4RI_DD_MZD(A);
    return NULL;
  }

  mzd_t *U = mzd_init_window(A, 0, 0, r, r);

  mzd_t *R = mzd_init(A->ncols, A->ncols - r);
  mzd_t *RU = mzd_init_window(R, 0, 0, r, R->ncols);

  for(rci_t i=0; i< r ; i++) {
    for(rci_t j=0; j< RU->ncols; j+=m4ri_radix) {
      const int workload = MIN(m4ri_radix, RU->ncols - j);
      mzd_xor_bits(RU, i, j, workload, mzd_read_bits(A, i, r + j, workload));
    }
  }

  mzd_trsm_upper_left(U, RU, cutoff);

  for(rci_t i = 0; i < R->ncols; ++i) {
    mzd_write_bit(R, r + i, i, 1);
  }
  mzd_apply_p_left_trans(R, Q);
  mzp_free(P);
  mzp_free(Q);
  mzd_free_window(RU);
  mzd_free_window(U);

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZD(R);
  return R;
}
