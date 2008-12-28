/**
 * \file lqup.h
 *
 * \brief PLUQ matrix decomposition routines.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 */


#ifndef PLUQ_H
#define PLUQ_H
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
 * Crossover point for PLUQ factorization.
 */

#define PLUQ_CUTOFF 512

/**
 * \brief PLUQ matrix decomposition.
 *
 * Computes the transposed PLUQ matrix decomposition using a block
 * recursive algorithm.
 *
 * If (P,L,U,Q) satisfy PLUQ = A, it returns (P^T, L, U, Q^T).
 *
 * The row echelon form (not reduced) can be read from the upper
 * triangular matrix U. See mzd_echelonize_pluq() for details.
 * 
 * This is the wrapper function including bounds checks. See
 * _mzd_pluq() for implementation details.
 *
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */

size_t mzd_pluq(packedmatrix *A, permutation *P, permutation * Q, const int cutoff);

/**
 * \brief PLUQ matrix decomposition.
 *
 * Computes the PLUQ matrix decomposition using a block recursive
 * algorithm.
 *
 * If (P,L,U,Q) satisfy PLUQ = A, this routine returns (P^T,L,U,Q^T).
 *
 * The row echelon form (not reduced) can be read from the upper
 * triangular matrix U*Q. See mzd_echelonize_pluq() for (reduced) row
 * echelon forms using PLUQ factorisation.
 * 
 * The matrix L and U are stored in place over A.
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 */

size_t _mzd_pluq(packedmatrix *A, permutation * P, permutation * Q, const int cutoff);

/**
 * \brief PLUQ matrix decomposition (naive base case).
 *
 * Computes the PLUQ matrix decomposition using the naive algorithm.
 *
 * If (P,L,U,Q) satisfy PLUQ = A, it returns (P^T, L, U, Q^T). 
 * 
 * The matrix L and U are stored in place over A.
 * 
 * \param A Input matrix
 * \param P Output row permutation matrix
 * \param Q Output column permutation matrix
 *
 * \sa mzd_pluq()
 */

size_t _mzd_pluq_naive(packedmatrix *A, permutation * P, permutation * Q);

/**
 * \brief (Reduced) row echelon form using PLUQ factorisation.
 *
 * \param A Matrix.
 * \param full Return the reduced row echelon form, not only upper triangular form.
 *
 * \wordoffset
 * \sa mzd_pluq()
 */


size_t mzd_echelonize_pluq(packedmatrix *A, int full);

#endif
