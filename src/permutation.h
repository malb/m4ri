/**
 * \file permutation.h
 *
 * \brief Permutation matrices.
 * 
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 *
 */
/******************************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de> 
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
******************************************************************************/
#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "misc.h"
#include "packedmatrix.h"

/**
 * \brief Permutations.
 */

typedef struct {
  /**
   * The swap operations in LAPACK format.
   */
  size_t *values;

  /**
   * The length of the swap array.
   */

  size_t length;

} permutation;

/**
 * Construct an identity permutation.
 * 
 * \param length Length of the permutation.
 */

permutation *mzp_init(size_t length);

/**
 * Free a permutation.
 * 
 * \param P Permutation to free.
 */

void mzp_free(permutation *P);


/**
 * \brief Create a window/view into the permutation matrix P.
 *
 * Use mzp_free_permutation_window() to free the window.
 *
 * \param P Permutaiton matrix
 * \param begin Starting index (inclusive)
 * \param end   Ending index   (exclusive)
 *
 */

permutation *mzp_init_window(permutation* P, size_t begin, size_t end);

/**
 * \brief Free a permutation matrix window created with
 * mzp_init_permutation_window().
 * 
 * \param condemned Permutation Matrix
 */

void mzp_free_window(permutation* condemned);

/**
 * Apply the permutation P to A from the left.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_left(packedmatrix *A, permutation *P);

/**
 * Apply the permutation P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_left_trans(packedmatrix *A, permutation *P);

/**
 * Apply the permutation P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_right(packedmatrix *A, permutation *P);

/**
 * Apply the permutation P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_right_trans(packedmatrix *A, permutation *P);

/**
 * Rotate zero columns to the end.
 *
 * Given a matrix M with zero columns from zs up to ze (exclusive) and
 * nonzero columns from ze to de (excluse) with zs < ze < de rotate
 * the zero columns to the end such that the the nonzero block comes
 * before the zero block.
 *
 * \param M Matrix.
 * \param zs Start index of the zero columns.
 * \param ze End index of the zero columns (exclusive).
 * \param de End index of the nonzero columns (exclusive).
 * \param zero_out actually write zero to the end.
 *
 */

void mzd_col_block_rotate(packedmatrix *M, size_t zs, size_t ze, size_t de, int zero_out);

/**
 * Print the permutation P
 *
 * \param P Permutation.
 */

void mzp_print(permutation *P);

#endif //PERMUTATION_H
