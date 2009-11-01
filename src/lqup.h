/**
 * \file lqup.h
 *
 * \brief PLUQ matrix decomposition routines.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 * 
 * \note This file should be called pluq.h and will be renamed in the
 * future.
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

#define LQUP_CUTOFF MIN(524288,CPU_L2_CACHE>>3)

/**
 * \brief PLUQ matrix decomposition.
 *
 * If (P,L,U,Q) satisfy PLUQ = A, it returns (P, L, U, Q^T).
 *
 * P and Q must be preallocated but they don't have to be
 * identity permutations. If cutoff is zero a value is chosen
 * automatically. It is recommended to set cutoff to zero for most
 * applications.
 *
 * The row echelon form (not reduced) can be read from the upper
 * triangular matrix U. See mzd_echelonize_pluq() for details.
 * 
 * This is the wrapper function including bounds checks. See
 * _mzd_pluq() for implementation details.
 *
 * \param A Input m x n matrix
 * \param P Output row permutation of length m
 * \param Qt Output column permutation matrix of length n
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa _mzd_pluq() _mzd_pluq_mmpf() mzd_echelonize_pluq()
 *
 * \wordoffset
 *
 * \return Rank of A.
 */

size_t mzd_pluq(mzd_t *A, mzp_t *P, mzp_t * Qt, const int cutoff);


/**
 * \brief LQUP matrix decomposition.
 *
 * Computes the transposed LQUP matrix decomposition using a block
 * recursive algorithm.
 *
 * If (L,Q,U,P) satisfy LQUP = A^T, it returns (L, Q^T, U, P).
 *
 * P and Qt must be preallocated but they don't have to be
 * identity permutations. If cutoff is zero a value is chosen
 * automatically. It is recommended to set cutoff to zero for most
 * applications.
 *
 * This is the wrapper function including bounds checks. See
 * _mzd_lqup() for implementation details.
 *
 * \param A Input m x n matrix
 * \param P Output row permutation of length m
 * \param Qt Output column permutation matrix of length n
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa _mzd_lqup() _mzd_pluq() _mzd_pluq_mmpf() mzd_echelonize_pluq()
 *
 * \wordoffset
 *
 * \return Rank of A.
 */

size_t mzd_lqup(mzd_t *A, mzp_t *P, mzp_t * Q, const int cutoff);

/**
 * \brief PLUQ matrix decomposition.
 *
 * See mzd_pluq() for details.
 *
 * \param A Input matrix
 * \param P Output row mzp_t matrix
 * \param Q Output column mzp_t matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa mzd_pluq()
 *
 * \wordoffset
 * \return Rank of A.
 */

size_t _mzd_pluq(mzd_t *A, mzp_t * P, mzp_t * Q, const int cutoff);

/**
 * \brief LQUP matrix decomposition.
 *
 * See mzd_lqup() for details.
 *
 * \param A Input matrix
 * \param P Output row mzp_t matrix
 * \param Qt Output column mzp_t matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa mzd_lqup()
 *
 * \wordoffset
 * \return Rank of A.
 */

size_t _mzd_lqup(mzd_t *A, mzp_t * P, mzp_t * Qt, const int cutoff);

/**
 * \brief PLUQ matrix decomposition (naive base case).
 *
 * See mzd_pluq() for details.
 * 
 * \param A Input matrix
 * \param P Output row mzp_t matrix
 * \param Q Output column mzp_t matrix
 *
 * \sa mzd_pluq()
 *
 * \wordoffset
 * \return Rank of A.
 */

size_t _mzd_pluq_naive(mzd_t *A, mzp_t * P, mzp_t * Q);

/**
 * \brief LQUP matrix decomposition (naive base case).
 *
 * See mzd_lqup() for details.
 * 
 * \param A Input matrix
 * \param P Output row mzp_t matrix
 * \param Qt Output column mzp_t matrix
 *
 * \sa mzd_lqup()
 *
 * \wordoffset
 * \return Rank of A.
 */

size_t _mzd_lqup_naive(mzd_t *A, mzp_t *P, mzp_t *Qt);

/**
 * \brief (Reduced) row echelon form using PLUQ factorisation.
 *
 * \param A Matrix.
 * \param full Return the reduced row echelon form, not only upper triangular form.
 *
 * \wordoffset
 *
 * \sa mzd_pluq()
 *
 * \return Rank of A.
 */


size_t mzd_echelonize_pluq(mzd_t *A, int full);

#endif
