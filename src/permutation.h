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

} mzp_t; // note that this is NOT mpz_t

/**
 * Construct an identity permutation.
 * 
 * \param length Length of the permutation.
 */

mzp_t *mzp_init(size_t length);

/**
 * Free a mzp_t.
 * 
 * \param P Mzp_T to free.
 */

void mzp_free(mzp_t *P);


/**
 * \brief Create a window/view into the mzp_t matrix P.
 *
 * Use mzp_free_mzp_t_window() to free the window.
 *
 * \param P Permutaiton matrix
 * \param begin Starting index (inclusive)
 * \param end   Ending index   (exclusive)
 *
 */

mzp_t *mzp_init_window(mzp_t* P, size_t begin, size_t end);

/**
 * \brief Free a mzp_t matrix window created with
 * mzp_init_mzp_t_window().
 * 
 * \param condemned Mzp_T Matrix
 */

void mzp_free_window(mzp_t* condemned);

/**
 * \brief Set the mzp_t P to the identity mzp_t. The only
 * allowed value is 1.
 *
 *
 * \param P Mzp_T
 * \param value 1
 *
 * \note This interface was chosen to be similar to mzd_set_ui().
 */

void mzp_set_ui(mzp_t *P, unsigned int value);


/**
 * Apply the mzp_t P to A from the left.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Mzp_T.
 */

void mzd_apply_p_left(mzd_t *A, mzp_t *P);

/**
 * Apply the mzp_t P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Mzp_T.
 */

void mzd_apply_p_left_trans(mzd_t *A, mzp_t *P);

/**
 * Apply the mzp_t P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Mzp_T.
 */

void mzd_apply_p_right(mzd_t *A, mzp_t *P);

/**
 * Apply the mzp_t P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Mzp_T.
 */

void mzd_apply_p_right_trans(mzd_t *A, mzp_t *P);



/**
 * Apply the mzp_t P to A from the right, but only on the lower triangular
 * part of the matrix A. 
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Mzp_T.
 */
void  mzd_apply_p_right_tri (mzd_t * A, mzp_t * Q);

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
 *
 */

void mzd_col_block_rotate(mzd_t *M, size_t zs, size_t ze, size_t de) ;

/**
 * Print the mzp_t P
 *
 * \param P Mzp_T.
 */

void mzp_print(mzp_t *P);

#endif //PERMUTATION_H
