/**
 * \file trsm.h
 *
 * This is scratch, experimental, unfinished code.
 *
 * \brief Triangular system solving with matrix routines.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 *
 * \internal
 */


#ifndef TRSM_H
#define TRSM_H
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
 * \brief TRiangular System solving with matrix.
 *
 * Solves X U = B where X and B are matrices, and U is upper triangular.
 * X is stored inplace on B
 * * 
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_upper_right for implementation details.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */

void mzd_trsm_upper_right (packedmatrix *U, packedmatrix *B, const int cutoff);

/**
 * \brief TRiangular System solving with matrix.
 *
 * Solves X U = B where X and B are matrices, and U is upper triangular.
 * This version assumes that the matrices are at an even position on
 * the RADIX grid and that their dimension is a multiple of RADIX.
 * X is stored inplace on B
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \internal
 */
void _mzd_trsm_upper_right (packedmatrix *U, packedmatrix *B, const int cutoff);

void _mzd_trsm_upper_right_even (packedmatrix *U, packedmatrix *B, const int cutoff);

void _mzd_trsm_upper_right_weird (packedmatrix *U, packedmatrix *B, const int cutoff);

/**
 * \brief TRiangular System solving with matrix.
 *
 * Solves L X = B where X and B are matrices, and L is lower triangular.
 *  X is stored inplace on B
 *  
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_lower_left for implementation details.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */

void mzd_trsm_lower_left (packedmatrix *L, packedmatrix *B, const int cutoff);

/**
 * \brief TRiangular System solving with matrix.
 *
 * Solves L X = B where X and B are matrices, and L is lower triangular.
 * X is stored inplace on B
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \internal
 */
void _mzd_trsm_lower_left (packedmatrix *L, packedmatrix *B, const int cutoff);

void _mzd_trsm_lower_left_even (packedmatrix *L, packedmatrix *B, const int cutoff);

void _mzd_trsm_lower_left_weird (packedmatrix *L, packedmatrix *B, const int cutoff);

#endif
