/**
 * \file trsm.h
 *
 * \brief Triangular system solving with Matrix routines.
 *
 * \author Clement Pernet <clement.pernet@gmail.com>
 */

#ifndef TRSM_H
#define TRSM_H
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
 * \brief TRiangular System solving with Matrix.
 *
 * Solves X U = B where X and B are matrices, and U is upper
 * triangular.  X is stored inplace on B.
 * 
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_upper_right for implementation details.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void mzd_trsm_upper_right(packedmatrix *U, packedmatrix *B, const int cutoff);

/**
 * \brief TRiangular System solving with Matrix.
 *
 * Solves X U = B where X and B are matrices, and U is upper triangular.
 *
 * \param U Input upper triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \internal
 */
void _mzd_trsm_upper_right(packedmatrix *U, packedmatrix *B, const int cutoff);

/**
 * \brief TRiangular System solving with Matrix.
 *
 * Solves L X = B where X and B are matrices, and L is lower triangular.
 * X is stored inplace on B.
 *  
 * This is the wrapper function including bounds checks. See
 * _mzd_trsm_lower_left for implementation details.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 */

void mzd_trsm_lower_left(packedmatrix *L, packedmatrix *B, const int cutoff);

/**
 * \brief TRiangular System solving with Matrix.
 *
 * Solves L X = B where X and B are matrices, and L is lower
 * triangular.  X is stored inplace on B.
 *
 * \param L Input lower triangular matrix.
 * \param B Input matrix, being overwritten by the solution matrix X
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \internal
 */
void _mzd_trsm_lower_left(packedmatrix *L, packedmatrix *B, const int cutoff);

#endif
