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

#include "solve.h"
#include "strassen.h"
#include "lqup.h"
#include "trsm.h"
#include "permutation.h"

void mzd_solve_left(packedmatrix *A, packedmatrix *B, const int cutoff, const int inconsistency_check) {    
  if(A->ncols > B->nrows)
    m4ri_die("mzd_solve_left: A ncols (%d) need to be lower than B nrows (%d).\n", A->ncols, B->nrows);

  _mzd_solve_left (A, B, cutoff, inconsistency_check);
}
 
void mzd_pluq_solve_left (packedmatrix *A, size_t rank, 
                          permutation *P, permutation *Q, 
                          packedmatrix *B, const int cutoff, const int inconsistency_check) 
{
  if(A->ncols > B->nrows)
    m4ri_die("mzd_pluq_solve_left: A ncols (%d) need to be lower than B nrows (%d).\n", A->ncols, B->nrows);
  if(P->length != A->nrows)
      m4ri_die("mzd_pluq_solve_left: A nrows (%d) need to match P size (%d).\n", A->nrows, P->length);
  if(Q->length != A->ncols)
      m4ri_die("mzd_pluq_solve_left: A ncols (%d) need to match Q size (%d).\n", A->ncols, P->length);

  _mzd_pluq_solve_left (A, rank, P, Q, B, cutoff, inconsistency_check);
}

void _mzd_pluq_solve_left (packedmatrix *A, size_t rank, 
                           permutation *P, permutation *Q, 
                           packedmatrix *B, const int cutoff, const int inconsistency_check) {
  /** A is supposed to store L lower triangular and U upper triangular
   *  B is modified in place 
   *  (Bi's in the comments are just modified versions of B)
   *  PLUQ = A
   *  1) P B2 = B1
   *  2) L B3 = B2
   *  3) U B4 = B3
   *  4) Q B5 = B4
   */
  
  /* P B2 = B1 or B2 = P^T B1*/
  mzd_apply_p_left_trans(B, P);
  
  /* L B3 = B2 */
  
  /* view on the upper part of L */
  packedmatrix *LU = mzd_init_window(A,0,0,rank,rank);
  packedmatrix *Y1 = mzd_init_window(B,0,0,rank,B->ncols);
  _mzd_trsm_lower_left(LU, Y1, cutoff);
  
  if (inconsistency_check) {
    /* Check for inconsistency */
    /** FASTER without this check
     * 
     * update with the lower part of L 
     */
    packedmatrix *H = mzd_init_window(A,rank,0,A->nrows,rank);
    packedmatrix *Y2 = mzd_init_window(B,rank,0,B->nrows,B->ncols);
    mzd_addmul(Y2, H, Y1, cutoff);
    /*
     * test whether Y2 is the zero matrix
     */
/*     if( !mzd_is_zero(Y2) ) { */
/*       printf("inconsistent system of size %llu x %llu\n", Y2->nrows, Y2->ncols); */
/*       printf("Y2="); */
/*       mzd_print_matrix(Y2); */
/*     } */
    mzd_free_window(H);
    mzd_free_window(Y2);
  }
  /* U B4 = B3 */
  _mzd_trsm_upper_left(LU, Y1, cutoff);
  mzd_free_window(LU);
  mzd_free_window(Y1);
  
  if (! inconsistency_check) {
    /** Default is to set the indefined bits to zero 
     * if inconsistency has been checked then 
     *    Y2 bits are already all zeroes
     * thus this clearing is not needed
     */
    if (rank < B->nrows) mzd_clear_bits(B,rank,0, (B->nrows-rank)*B->ncols);
  }
  
  /* Q B5 = B4 or B5 = Q^T B4*/
  mzd_apply_p_left_trans(B, Q);
  /* P L U Q B5 = B1 */
}


void _mzd_solve_left (packedmatrix *A, packedmatrix *B, const int cutoff, const int inconsistency_check) {
  /**
   *  B is modified in place 
   *  (Bi's in the comments are just modified versions of B)
   *  1) PLUQ = A
   *  2) P B2 = B1
   *  3) L B3 = B2
   *  4) U B4 = B3
   *  5) Q B5 = B4
   */
  permutation * P = mzp_init(A->nrows);
  permutation * Q = mzp_init(A->ncols);
  
  /* PLUQ = A */
  size_t rank = _mzd_pluq(A, P, Q, cutoff);  
  /* 2, 3, 4, 5 */
  mzd_pluq_solve_left(A, rank, P, Q, B, cutoff, inconsistency_check);
  
  mzp_free(P);
  mzp_free(Q);
}

