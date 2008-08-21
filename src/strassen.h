/**
 * \file strassen.h
 *
 * \brief Matrix operations using Strassen's formulas including
 * Winograd's improvements.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */


#ifndef STRASSEN_H
#define STRASSEN_H
 /*******************************************************************
 *
 *            M4RI: Method of the Four Russians Inversion
 *
 *       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhu.ac.uk>
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

#include <math.h>
#include "misc.h"
#include "packedmatrix.h"
#include "brilliantrussian.h"

/**
 * \brief Matrix multiplication via the Strassen-Winograd matrix
 * multiplication algorithm, i.e. compute C = AB.
 * 
 * This is the wrapper function including bounds checks. See
 * _mzd_mul_even for implementation details.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 */

packedmatrix *mzd_mul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * \brief Matrix multiplication and in-place addition via the
 * Strassen-Winograd matrix multiplication algorithm, i.e. compute 
 * C = C+ AB.
 * 
 * This is the wrapper function including bounds checks. See
 * _mzd_addmul_even for implementation details.
 *
 * \param C product matrix
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 */

packedmatrix *mzd_addmul(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * \brief Matrix multiplication via the Strassen-Winograd matrix
 * multiplication algorithm, i.e. compute C = AB.
 * 
 * This is the actual implementation. Any matrix where either the
 * number of rows or the number of columns is smaller than cutoff is
 * processed using the M4RM algorithm.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \note This implementation is heavily inspired by the function
 * strassen_window_multiply_c in Sage 3.0; For reference see
 * http://www.sagemath.org
 */

packedmatrix *_mzd_mul_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * \brief Matrix multiplication and in-place addition via the
 * Strassen-Winograd matrix multiplication algorithm, i.e. compute 
 * C = C+ AB.
 * 
 * This is the actual implementation. Any matrix where either the
 * number of rows or the number of columns is smaller than cutoff is
 * processed using the M4RM algorithm.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \note This implementation is heavily inspired by the function
 * strassen_window_multiply_c in Sage 3.0; For reference see
 * http://www.sagemath.org
 */

packedmatrix *_mzd_addmul_even(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * \brief Matrix multiplication and in-place addition via the
 * Strassen-Winograd matrix multiplication algorithm, i.e. compute 
 * C = C + AB.
 * 
 * The matrices A and B are respectively m x k and k x n, and can be not
 * aligned on the RADIX grid.
 * 
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */

packedmatrix *_mzd_addmul (packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * C = A*B + C for matrices with offsets != 0
 *
 * This is scratch code.
 *
 * \internal
 */

packedmatrix *_mzd_addmul_weird_weird (packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * C = A*B + C for A with offset == 0 and B with offset != 0.
 *
 * This is scratch code.
 *
 * \internal
 */

packedmatrix *_mzd_addmul_weird_even (packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * C = A*B + C for A with offset != 0 and B with offset == 0.
 *
 * This is scratch code.
 *
 * \internal
 */

packedmatrix *_mzd_addmul_even_weird (packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * The default cutoff for Strassen-Winograd multiplication. It should
 * hold hold that 2 * (n^2)/8 fits into the L2 cache.
 */

#ifndef STRASSEN_MUL_CUTOFF
#define STRASSEN_MUL_CUTOFF ((int)sqrt((double)(4*CPU_L2_CACHE)))
#endif// STRASSEN_MUL_CUTOFF

#endif //STRASSEN_H
