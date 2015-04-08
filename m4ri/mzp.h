/**
 * \file mzp.h
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

#ifndef M4RI_MZP
#define M4RI_MZP

#include <m4ri/mzd.h>

/**
 * \brief Permutations.
 */

typedef struct mzp_t {
  /**
   * The swap operations in LAPACK format.
   */
  rci_t *values;

  /**
   * The length of the swap array.
   */

  rci_t length;

} mzp_t; // note that this is NOT mpz_t

/**
 * Construct an identity permutation.
 * 
 * \param length Length of the permutation.
 */

mzp_t *mzp_init(rci_t length);

/**
 * Free a mzp_t.
 * 
 * \param P Permutation to free.
 */

void mzp_free(mzp_t *P);

/**
 * \brief Create a window/view into the permutation P.
 *
 * Use mzp_free_mzp_t_window() to free the window.
 *
 * \param P Permutation matrix
 * \param begin Starting index (inclusive)
 * \param end   Ending index   (exclusive)
 *
 */

mzp_t *mzp_init_window(mzp_t *P, rci_t begin, rci_t end);

/**
 * \brief Free a permutation window created with
 * mzp_init_mzp_t_window().
 * 
 * \param condemned Permutation Matrix
 */

void mzp_free_window(mzp_t *condemned);


/**
 * \brief copy permutation Q to P
 *
 * \param P Target permutation matrix (may be NULL)
 * \param Q Source permutation matrix (must not be NULL)
 */

mzp_t *mzp_copy(mzp_t *P, const mzp_t *Q);

/**
 * \brief Set the permutation P to the identity permutation. The only
 * allowed value is 1.
 *
 *
 * \param P Permutation
 * \param value 1
 *
 * \note This interface was chosen to be similar to mzd_set_ui().
 */

void mzp_set_ui(mzp_t *P, unsigned int value);

/**
 * Apply the permutation P to A from the left.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_left(mzd_t *A, mzp_t const *P);

/**
 * Apply the permutation P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_left_trans(mzd_t *A, mzp_t const *P);

/**
 * Apply the permutation P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_right(mzd_t *A, mzp_t const *P);

/**
 * Apply the permutation P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_right_trans(mzd_t *A, mzp_t const *P);

/**
 * Apply the permutation P to A from the right starting at start_row.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 * \param start_row Start swapping at this row.
 * \param start_col Start swapping at this column.
 */

void mzd_apply_p_right_even_capped(mzd_t *A, mzp_t const *P, rci_t start_row, rci_t start_col);

/**
 * Apply the permutation P^T to A from the right starting at start_row.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 * \param start_row Start swapping at this row.
 * \param start_col Start swapping at this column.
 */

void mzd_apply_p_right_trans_even_capped(mzd_t *A, mzp_t const *P, rci_t start_row, rci_t start_col);

/**
 * Apply the mzp_t P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void mzd_apply_p_right_trans(mzd_t *A, mzp_t const *P);

/**
 * Apply the permutation P to A from the right, but only on the upper
 * the matrix A above the main diagonal.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param Q Permutation.
 */
void  mzd_apply_p_right_trans_tri(mzd_t *A, mzp_t const *Q);

/**
 * Print the mzp_t P
 *
 * \param P Permutation.
 */

void mzp_print(mzp_t const *P);

/**
 * Compresses the matrix L in a step in blockwise-recursive PLE
 * decomposition.
 *
 * \param A Matrix.
 * \param r1 Rank of left matrix.
 * \param n1 Column cut which separates left and right matrix.
 * \param r2 Rank of right matrix.
 */

void _mzd_compress_l(mzd_t *A, rci_t r1, rci_t n1, rci_t r2);

#endif // M4RI_MZP
