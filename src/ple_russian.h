/**
 * \file ple_mmpf.h
 * \brief PLE and PLUQ factorization using Gray codes.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 * \example testsuite/test_ple.c
 */

#ifndef M4RI_PLE_RUSSIAN
#define M4RI_PLE_RUSSIAN

 /*******************************************************************
 *
 *                 M4RI:  Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008-2011 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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
 * \brief PLE matrix decomposition of A using Gray codes.
 *
 * Returns (P,L,E,Q) satisfying PLE = A where P is a permutation
 * matrix of dimension m x m, L is m x r unit lower triangular and S
 * is an r x n matrix which is upper triangular except that its
 * columns are permuted, that is E = UQ for U r x n upper triangular
 * and Q is a n x n permutation matrix. The matrix L and E are stored
 * in place over A.
 *
 * \param A Matrix.
 * \param P Preallocated row permutation.
 * \param Q Preallocated column permutation.
 * \param k Size of Gray code tables.
 *
 * \wordoffset
 *
 * \return Rank of A.
 */

rci_t _mzd_ple_russian(mzd_t *A, mzp_t *P, mzp_t *Q, int k);

/**
 * \brief PLUQ matrix decomposition of A using Gray codes.
 *
 * Returns (P,L,U,Q) satisfying PLUQ = A where P and Q are two
 * permutation matrices, of dimension respectively m x m and n x n, L
 * is m x r unit lower triangular and U is r x n upper triangular.
 *
 * \param A Matrix.
 * \param P Preallocated row permutation.
 * \param Q Preallocated column permutation.
 * \param k Size of Gray code tables.
 *
 * \wordoffset
 *
 * \return Rank of A.
 */

rci_t _mzd_pluq_russian(mzd_t *A, mzp_t *P, mzp_t *Q, int k);

/**
 * \brief PLE matrix decomposition of a submatrix for up to k columns
 * starting at (r,c).
 *
 * Updates P and Q and modifies A in place. The buffer done afterwards
 * holds how far a particular row was already eliminated.
 *
 * \param A Matrix.
 * \param start_row Row Offset.
 * \param stop_row Up to which row the matrix should be processed (exclusive).
 * \param start_col Column Offset.
 * \param k Size of Gray code tables.
 * \param P Preallocated row permutation.
 * \param Q Preallocated column permutation.
 * \param pivots which column holds the i-th pivot
 * \param done Preallocated temporary buffer.
 * \param done_row Stores the last row which is already reduced processed after function execution.
 * \param splitblock First block which is not considered by this function.
 *
 * \retval knar rank of the considered submatrix
 */

int _mzd_ple_submatrix(mzd_t *A,
                       rci_t const start_row, rci_t const stop_row,
                       rci_t const start_col, int const k,
                       mzp_t *P, mzp_t *Q, rci_t *pivots, rci_t *done, rci_t *done_row,
                       wi_t const splitblock);

/**
 * \brief Helper function for bulk of work in PLE decomposition.
 *
 * This funcion computes M[i] = M[i] + T0[E1[x0]] + T1[E1[x1]] for startrow <=
 * i < stoprow where (1) x0 is read from M[i][startcol:startcol+k0],
 * then (2) T0[E0[x0]] is added, then (3) x1 is read from
 * M[i][startcol:startcol+k1] and then T1[E1[x1]] added.
 *
 * \param M Matrix to work on.
 * \param startrow Start processing in this row.
 * \param stoprow  Stop processing before this row.
 * \param startcol Start processing from this column and read from this column.
 * \param k0       Read k0 bits starting at M[i,startcol]
 * \param T0       Matrix with at least 2^k0 rows
 * \param E0       Permutation on 0..2^k0-1.
 * \param k1       Read k1 bits starting at M[i,startcol]
 * \param T1       Matrix with at least 2^k1 rows
 * \param E1       Permutation on 0..2^k1-1.
 */

void mzd_process_rows2_ple(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol,
                           int const k0, mzd_t const *T0, rci_t const *E0,
                           int const k1, mzd_t const *T1, rci_t const *E1);

void mzd_process_rows3_ple(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol,
                           int const k0, mzd_t const *T0, rci_t const *E0,
                           int const k1, mzd_t const *T1, rci_t const *E1,
			   int const k2, mzd_t const *T2, rci_t const *E2);

void mzd_process_rows4_ple(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol,
                           int const k0, mzd_t const *T0, rci_t const *E0,
                           int const k1, mzd_t const *T1, rci_t const *E1,
                           int const k2, mzd_t const *T2, rci_t const *E2,
                           int const k3, mzd_t const *T3, rci_t const *E3);

/**
 * \brief Extract the k x A::ncols echelon form submatrix of A starting at row r and column c.
 *
 * \param E Storage for k x A::ncols matrix.
 * \param A Source matrix.
 * \param r Row index.
 * \param c Column index.
 * \param k Rank of E.
 * \param k Map from i to column of i-th pivot.
 */
mzd_t *_mzd_ple_to_e(mzd_t *E, mzd_t const *A, rci_t r, rci_t c, int k, rci_t *offsets);

#endif // M4RI_PLE_RUSSIAN
