/**
 * \file ple.h
 *
 * \brief PLE and PLUQ matrix decomposition routines.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 *
 */

#ifndef M4RI_PLUQ_H
#define M4RI_PLUQ_H

/*******************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008, 2009 Clement Pernet <clement.pernet@gmail.com>
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
#include <m4ri/mzp.h>

/**
 * Crossover point for PLUQ factorization.
 */

#define __M4RI_PLE_CUTOFF MIN(524288, __M4RI_CPU_L3_CACHE >> 3)

/**
 * \brief PLUQ matrix decomposition.
 *
 * Returns (P,L,U,Q) satisfying PLUQ = A where P and Q are two
 * permutation matrices, of dimension respectively m x m and n x n, L
 * is m x r unit lower triangular and U is r x n upper triangular.
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
 * \param Q Output column permutation matrix of length n
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa _mzd_pluq() _mzd_pluq_mmpf() mzd_echelonize_pluq()
 *
 * \return Rank of A.
 */

rci_t mzd_pluq(mzd_t *A, mzp_t *P, mzp_t *Q, const int cutoff);

/**
 * \brief PLE matrix decomposition.
 *
 * Computes the PLE matrix decomposition using a block recursive
 * algorithm.
 *
 * Returns (P,L,S,Q) satisfying PLE = A where P is a permutation matrix
 * of dimension m x m, L is m x r unit lower triangular and S is an r
 * x n matrix which is upper triangular except that its columns are
 * permuted, that is S = UQ for U r x n upper triangular and Q is a n
 * x n permutation matrix. The matrix L and S are stored in place over
 * A.
 *
 * P and Q must be preallocated but they don't have to be
 * identity permutations. If cutoff is zero a value is chosen
 * automatically. It is recommended to set cutoff to zero for most
 * applications.
 *
 * This is the wrapper function including bounds checks. See
 * _mzd_ple() for implementation details.
 *
 * \param A Input m x n matrix
 * \param P Output row permutation of length m
 * \param Q Output column permutation matrix of length n
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa _mzd_ple() _mzd_pluq() _mzd_pluq_mmpf() mzd_echelonize_pluq()
 *
 * \return Rank of A.
 */

rci_t mzd_ple(mzd_t *A, mzp_t *P, mzp_t *Q, const int cutoff);

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
 * \return Rank of A.
 */

rci_t _mzd_pluq(mzd_t *A, mzp_t *P, mzp_t *Q, const int cutoff);

/**
 * \brief PLE matrix decomposition.
 *
 * See mzd_ple() for details.
 *
 * \param A Input matrix
 * \param P Output row mzp_t matrix
 * \param Qt Output column mzp_t matrix
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \sa mzd_ple()
 *
 * \return Rank of A.
 */

rci_t _mzd_ple(mzd_t *A, mzp_t *P, mzp_t *Qt, const int cutoff);

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
 * \return Rank of A.
 */

rci_t _mzd_pluq_naive(mzd_t *A, mzp_t *P, mzp_t *Q);

/**
 * \brief PLE matrix decomposition (naive base case).
 *
 * See mzd_ple() for details.
 *
 * \param A Input matrix
 * \param P Output row mzp_t matrix
 * \param Qt Output column mzp_t matrix
 *
 * \sa mzd_ple()
 *
 * \return Rank of A.
 */

rci_t _mzd_ple_naive(mzd_t *A, mzp_t *P, mzp_t *Qt);

#endif  // M4RI_PLUQ_H
