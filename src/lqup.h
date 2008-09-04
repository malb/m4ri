/**
 * \file lqup.h
 *
 * \brief LQUP matrix decomposition routines
 *
 * This is scratch, experimental code.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 *
 * \internal
 */


#ifndef LQUP_H
#define LQUP_H
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

#include "misc.h"
#include "packedmatrix.h"

/**
 * Crossover point for LQUP factorization.
 */

#define LQUP_CUTOFF 1024

/**
 * \brief LQUP matrix decomposition (unfinished).
 *
 * Computes the transposed LQUP matrix decomposition using a block recursive algorithm
 *
 * If (L,Q,U,P) satisfy LQUP = A^T, it returns (L^T, Q^T, U^T, P^T).
 * The Row echelon form (not reduced) can be read from the upper triangular matrix L^T.
 * 
 * This is the wrapper function including bounds checks. See
 * _mzd_lqup for implementation details.
 *
 * The matrix L and U are stored in place over A.
 * L^T is represented by the matrix Q^T L^T Q
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \internal
 */
size_t mzd_lqup(packedmatrix *A, permutation *P, permutation * Q, const int cutoff);

/**
 * \brief LQUP matrix decomposition (unfinished).
 *
 * Computes the transposed LQUP matrix decomposition using a block recursive algorithm
 *
 * If (L,Q,U,P) satisfy LQUP = A^T, it returns (L^T, Q^T, U^T, P^T).
 * The Row echelon form (not reduced) can be read from the upper triangular matrix L^T.
 * 
 * The matrix L and U are stored in place over A.
 * L^T is represented by the matrix Q^T L^T Q
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \internal
 */

size_t _mzd_lqup(packedmatrix *A, permutation * P, permutation * Q, const int cutoff);

/**
 * \brief LQUP matrix decomposition (naiv base case).
 *
 * Computes the LQUP matrix decomposition using a block recursive
 * algorithm
 *
 * If (L,Q,U,P) satisfy LQUP = A, it returns (L, Q, U, P).  The Row
 * echelon form (not reduced) can be read from the upper triangular
 * matrix L.
 * 
 * The matrix L and U are stored in place over A.  L is represented by
 * the matrix Q L Q
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \internal
 */

size_t _mzd_lqup_naiv(packedmatrix *A, permutation * P, permutation * Q);

#endif
