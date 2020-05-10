/**
 * \file triangular.h
 *
 * \brief Triangular system solving with Matrix routines.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 */

#ifndef M4RI_TRSM_H
#define M4RI_TRSM_H

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

#include <m4ri/mzd.h>

/**
 * \brief Solves X U = B with X and B matrices and U upper triangular.
 *
 * X is stored inplace on B.
 *
 * \attention Note, that the 'right' variants of TRSM are slower than
 * the 'left' variants.
 *
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_upper_right() for implementation details.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void mzd_trsm_upper_right(mzd_t const *U, mzd_t *B, const int cutoff);

/**
 * \brief Solves X U = B with X and B matrices and U upper triangular.
 *
 * X is stored inplace on B.
 *
 * \attention Note, that the 'right' variants of TRSM are slower than
 * the 'left' variants.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */
void _mzd_trsm_upper_right(mzd_t const *U, mzd_t *B, const int cutoff);

/**
 * \brief Solves X L = B with X and B matrices and L lower triangular.
 *
 * X is stored inplace on B.
 *
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_upper_right() for implementation details.
 *
 * \attention Note, that the 'right' variants of TRSM are slower than the 'left'
 * variants.
 *
 * \param L Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void mzd_trsm_lower_right(mzd_t const *L, mzd_t *B, const int cutoff);

/**
 * \brief Solves X L = B with X and B with matrices and L lower
 * triangular.
 *
 * This version assumes that the matrices are at an even position on
 * the m4ri_radix grid and that their dimension is a multiple of m4ri_radix.
 * X is stored inplace on B.
 *
 * \attention Note, that the 'right' variants of TRSM are slower than
 * the 'left' variants.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */
void _mzd_trsm_lower_right(mzd_t const *L, mzd_t *B, const int cutoff);

/**
 * \brief Solves L X = B with X and B matrices and L lower triangular.
 *
 * X is stored inplace on B.
 *
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_lower_left() for implementation details.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void mzd_trsm_lower_left(mzd_t const *L, mzd_t *B, const int cutoff);

/**
 * \brief Solves L X = B with X and B matrices and L lower triangular.
 *
 * X is stored inplace on B.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void _mzd_trsm_lower_left(mzd_t const *L, mzd_t *B, const int cutoff);

/**
 * \brief Solves U X = B with X and B matrices and U upper triangular.
 *
 * X is stored inplace on B.
 *
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_upper_left() for implementation details.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void mzd_trsm_upper_left(mzd_t const *U, mzd_t *B, const int cutoff);

/**
 * \brief Solves U X = B with X and B matrices and U upper triangular.
 *
 * X is stored inplace on B.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */
void _mzd_trsm_upper_left(mzd_t const *U, mzd_t *B, const int cutoff);

/**
 * \brief Invert the upper triangular matrix A by reduction to matrix multiplication.
 *
 * \param A Matrix to be inverted (overwritten).
 *
 * \return Inverse of A or throws an error
 */

mzd_t *mzd_trtri_upper(mzd_t *A);

#endif  // M4RI_TRSM_H
