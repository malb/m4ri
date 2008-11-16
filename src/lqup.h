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

/**
 * Crossover point for LQUP factorization.
 */

#define LQUP_CUTOFF 64

/**
 * \brief PLUQ matrix decomposition.
 *
 * Computes the transposed LQUP matrix decomposition using a block
 * recursive algorithm.
 *
 * If (P,L,U,Q) satisfy P^T LU Q^T = A, it returns (P, L, U, Q).
 *
 * The Row echelon form (not reduced) can be read from the upper
 * triangular matrix U.
 * 
 * This is the wrapper function including bounds checks. See _mzd_lqup
 * for implementation details.
 *
 * The matrix L and U are stored in place over A.  U is represented
 * by the matrix Q U Q^T
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */

size_t mzd_lqup(packedmatrix *A, permutation *P, permutation * Q, const int cutoff);

/**
 * \brief LQUP matrix decomposition.
 *
 * Computes the transposed LQUP matrix decomposition using a block
 * recursive algorithm.
 *
 * If (L,Q,U,P) satisfy LQUP = A^T, it returns (L^T, Q^T, U^T, P^T).
 *
 * The Row echelon form (not reduced) can be read from the upper
 * triangular matrix L^T.
 * 
 * The matrix L and U are stored in place over A.  L^T is represented
 * by the matrix Q^T L^T Q
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 */

size_t _mzd_lqup(packedmatrix *A, permutation * P, permutation * Q, const int cutoff);

/**
 * \brief LQUP matrix decomposition (naiv base case).
 *
 * Computes the LQUP matrix decomposition using the naive algorithm.
 *
 * If (L,Q,U,P) satisfy LQUP = A, it returns (L, Q, U, P). 
 * 
 * The Row echelon form (not reduced) can be read from the upper
 * triangular matrix L.
 * 
 * The matrix L and U are stored in place over A.  L is represented by
 * the matrix Q L Q.
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 */

size_t _mzd_lqup_naiv(packedmatrix *A, permutation * P, permutation * Q);

#endif
