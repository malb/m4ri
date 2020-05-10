/**
 * \file triangular_russian.h
 * \brief TRSM and TRTRI via Gray code tables.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RI_TRIANGULAR_RUSSIAN
#define M4RI_TRIANGULAR_RUSSIAN

/*******************************************************************
 *
 *                 M4RI:  Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008-2011 Martin Albrecht <martinralbrecht@googlemail.com>
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
 * \brief Solves L X = B with X and B matrices and L lower triangular using Gray code tables.
 *
 * X is stored inplace on B.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param k Size of Gray code tables or zero for automatic choice (recommended).
 */

void _mzd_trsm_lower_left_russian(mzd_t const *L, mzd_t *B, int k);

/**
 * \brief Solves U X = B with X and B matrices and U upper triangular using Gray code tables.
 *
 * X is stored inplace on B.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param k Size of Gray code tables or zero for automatic choice (recommended).
 */

void _mzd_trsm_upper_left_russian(mzd_t const *U, mzd_t *B, int k);

/**
 * \brief Invert the upper triangular matrix A using Kronrod's method.
 *
 * \param A Matrix to be inverted (overwritten).
 * \param k Table size parameter, may be 0 for automatic choice.
 *
 * \return Inverse of A or throws an error
 */

mzd_t *mzd_trtri_upper_russian(mzd_t *A, int k);

#endif  // M4RI_TRIANGULAR_RUSSIAN
