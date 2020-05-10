/**
 * \file strassen.h
 *
 * \brief Matrix operations using Strassen's formulas including
 * Winograd's improvements.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */

#ifndef M4RI_STRASSEN_H
#define M4RI_STRASSEN_H

/*******************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 *    Copyright (C) 2008 Clement Pernet <pernet@math.washington.edu>
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

#include <m4ri/brilliantrussian.h>
#include <m4ri/mzd.h>
#include <math.h>

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

mzd_t *mzd_mul(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

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

mzd_t *mzd_addmul(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

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

mzd_t *_mzd_mul_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

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

mzd_t *_mzd_addmul_even(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

/**
 * \brief Matrix multiplication and in-place addition via the
 * Strassen-Winograd matrix multiplication algorithm, i.e. compute
 * C = C + AB.
 *
 * The matrices A and B are respectively m x k and k x n, and can be not
 * aligned on the m4ri_radix grid.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 */

mzd_t *_mzd_addmul(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

/**
 * The default cutoff for Strassen-Winograd multiplication. It should
 * hold hold that 2 * (n^2)/8 fits into the L2 cache.
 */

#ifndef __M4RI_STRASSEN_MUL_CUTOFF
#define __M4RI_STRASSEN_MUL_CUTOFF MIN(((int)sqrt((double)(4 * __M4RI_CPU_L3_CACHE))), 4096)
#endif

#endif  // M4RI_STRASSEN_H
