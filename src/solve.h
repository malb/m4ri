/**
 * \file solve.h
 *
 * \brief system solving with matrix routines.
 *
 * \author Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 */
#ifndef SOLVE_H
#define SOLVE_H
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

#include <stdio.h>
#include "misc.h"
#include "permutation.h"
#include "packedmatrix.h"

/**
 * \brief System solving with matrix.
 *
 * Solves A X = B where A, and B are input matrices.  The solution X
 * is stored inplace on B.
 *
 * \param A Input matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 * \param inconsistency_check decide wether or not to check for
 *        incosistency (faster without but output not defined if
 *        system is not consistent).
 *
 */
void mzd_solve_left(packedmatrix *A, packedmatrix *B, const int cutoff, 
                    const int inconsistency_check);

/**
 * \brief System solving with matrix.
 *
 * Solves (P L U Q) X = B where 
 * A is an input matrix supposed to store both:
 *      an upper right triangular matrix U
 *      a lower left unitary triangular matrix L
 * P and Q are permutation matrices
 * B is the input matrix to be solved.
 * The solution X is stored inplace on B
 * This version assumes that the matrices are at an even position on
 * the RADIX grid and that their dimension is a multiple of RADIX.
 *
 * \param A Input upper/lower triangular matrices.
 * \param rank is rank of A.
 * \param P Input row permutation matrix.
 * \param Q Input column permutation matrix.
 * \param B Input matrix, being overwritten by the solution matrix X.
 * \param cutoff Minimal dimension for Strassen recursion.
 * \param inconsistency_check decide whether or not to check for
 *        incosistency (faster without but output not defined if
 *        system is not consistent).
 */
void mzd_pluq_solve_left (packedmatrix *A, size_t rank, 
                          permutation *P, permutation *Q, 
                          packedmatrix *B, const int cutoff, const int inconsistency_check);

/**
 * \brief System solving with matrix.
 *
 * Solves (P L U Q) X = B where 
 * A is an input matrix supposed to store both:
 *      an upper right triangular matrix U
 *      a lower left unitary triangular matrix L
 * P and Q are permutation matrices
 * B is the input matrix to be solved.
 * The solution X is stored inplace on B
 * This version assumes that the matrices are at an even position on
 * the RADIX grid and that their dimension is a multiple of RADIX.
 *
 * \param A Input upper/lower triangular matrices.
 * \param rank is rank of A.
 * \param P Input row permutation matrix.
 * \param Q Input column permutation matrix.
 * \param B Input matrix, being overwritten by the solution matrix X.
 * \param cutoff Minimal dimension for Strassen recursion.
 * \param inconsistency_check decide whether or not to check for
 *        incosistency (faster without but output not defined if
 *        system is not consistent).
 *
 */
void _mzd_pluq_solve_left(packedmatrix *A, size_t rank, 
                          permutation *P, permutation *Q, 
                          packedmatrix *B, const int cutoff, const int inconsistency_check);

/**
 * \brief System solving with matrix.
 *
 * Solves A X = B where A, and B are input matrices.
 * The solution X is stored inplace on B
 * This version assumes that the matrices are at an even position on
 * the RADIX grid and that their dimension is a multiple of RADIX.
 *
 * \param A Input matrix.
 * \param B Input matrix, being overwritten by the solution matrix X.
 * \param cutoff Minimal dimension for Strassen recursion.
 * \param inconsistency_check decide whether or not to check for
 *        incosistency (faster without but output not defined if
 *        system is not consistent).
 *
 */
void _mzd_solve_left(packedmatrix *A, packedmatrix *B, const int cutoff, const int inconsistency_check);


#endif // SOLVE_H
