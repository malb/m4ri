/**
 * \file solve.h
 *
 * \brief System solving with matrix routines.
 *
 * \author Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 */

#ifndef M4RI_SOLVE_H
#define M4RI_SOLVE_H

 /*******************************************************************
 *
 *            M4RI: Linear Algebra over GF(2)
 *
 *       Copyright (C) 2008 Jean-Guillaume.Dumas@imag.fr
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

#include <m4ri/mzp.h>
#include <m4ri/mzd.h>

/**
 * \brief Solves A X = B with A and B matrices. 
 *
 * The solution X is stored inplace on B.
 *
 * \param A Input matrix (overwritten).
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion (default: 0).
 * \param inconsistency_check decide wether or not to perform a check
 *        for incosistency (faster without but output not defined if
 *        system is not consistent).
 * \return 0 if a solution was found, -1 otherwise
 */
int mzd_solve_left(mzd_t *A, mzd_t *B, int const cutoff, int const inconsistency_check);

/**
 * \brief Solves (P L U Q) X = B
 * 
 * A is an input matrix supposed to store both:
 * \li  an upper right triangular matrix U
 * \li  a lower left unitary triangular matrix L.
 *
 * The solution X is stored inplace on B
 *
 * This version assumes that the matrices are at an even position on
 * the m4ri_radix grid and that their dimension is a multiple of m4ri_radix.
 *
 * \param A Input upper/lower triangular matrices.
 * \param rank is rank of A.
 * \param P Input row permutation matrix.
 * \param Q Input column permutation matrix.
 * \param B Input matrix, being overwritten by the solution matrix X.
 * \param cutoff Minimal dimension for Strassen recursion (default: 0).
 * \param inconsistency_check decide whether or not to perform a check
 *        for incosistency (faster without but output not defined if
 *        system is not consistent).  \return 0 if a solution was
 *        found, -1 otherwise
 * \return 0 if a solution was found, -1 otherwise
 */
int mzd_pluq_solve_left (mzd_t const *A, rci_t rank, 
                         mzp_t const *P, mzp_t const *Q, 
                         mzd_t *B, int const cutoff, int const inconsistency_check);

/**
 * \brief  Solves (P L U Q) X = B
 *
 * A is an input matrix supposed to store both:
 * \li an upper right triangular matrix U
 * \li a lower left unitary triangular matrix L.

 * The solution X is stored inplace on B.
 *
 * This version assumes that the matrices are at an even position on
 * the m4ri_radix grid and that their dimension is a multiple of m4ri_radix.
 *
 * \param A Input upper/lower triangular matrices.
 * \param rank is rank of A.
 * \param P Input row permutation matrix.
 * \param Q Input column permutation matrix.
 * \param B Input matrix, being overwritten by the solution matrix X.
 * \param cutoff Minimal dimension for Strassen recursion (default: 0).
 * \param inconsistency_check decide whether or not to perform a check
 *        for incosistency (faster without but output not defined if
 *        system is not consistent).  \return 0 if a solution was
 *        found, -1 otherwise
 * \return 0 if a solution was found, -1 otherwise
 */
int _mzd_pluq_solve_left(mzd_t const *A, rci_t rank, 
                         mzp_t const *P, mzp_t const *Q, 
                         mzd_t *B, int const cutoff, int const inconsistency_check);

/**
 * \brief Solves A X = B with A and B matrices.
 *
 * The solution X is stored inplace on B.
 *
 * This version assumes that the matrices are at an even position on
 * the m4ri_radix grid and that their dimension is a multiple of m4ri_radix.
 *
 * \param A Input matrix (overwritten).
 * \param B Input matrix, being overwritten by the solution matrix X.
 * \param cutoff Minimal dimension for Strassen recursion (default: 0).
 * \param inconsistency_check decide whether or not to perform a check
 *        for incosistency (faster without but output not defined if
 *        system is not consistent).  \return 0 if a solution was
 *        found, -1 otherwise
 * \return 0 if a solution was found, -1 otherwise
 */
int _mzd_solve_left(mzd_t *A, mzd_t *B, int const cutoff, int const inconsistency_check);

/**
 * \brief Solve X for A X = 0.
 *
 * If r is the rank of the nr x nc matrix A, return the nc x (nc-r)
 * matrix X such that A*X == 0 and that the columns of X are linearly
 * independent.
 *
 * \param A Input matrix (overwritten).
 * \param cutoff Minimal dimension for Strassen recursion (default: 0).
 *
 * \sa mzd_pluq()
 *
 * \return X, NULL if kernel is empty
 */

mzd_t *mzd_kernel_left_pluq(mzd_t *A, int const cutoff);

#endif // M4RI_SOLVE_H
