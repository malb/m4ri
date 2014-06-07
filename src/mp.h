/**
 * \file mp.h
 *
 * \brief multicore matrix operations 
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RI_MP_H
#define M4RI_MP_H

/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2014 Martin Albrecht <martinralbrecht@googlemail.com>
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

/**
 * \brief Matrix multiplication via the cubic multiplication algorithm on
 * multiple cores, i.e. compute C = AB.
 * 
 * This is the wrapper function including bounds checks. See _mzd_mul_mp4 for
 * implementation details.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for recursion.
 */

mzd_t *mzd_mul_mp(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

/**
 * \brief Matrix multiplication and in-place addition via the cubic matrix
 * multiplication algorithm on multiple cores, i.e. compute C = C+ AB.
 * 
 * This is the wrapper function including bounds checks. See _mzd_addmul_mp4 for
 * implementation details.
 *
 * \param C product matrix
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for recursion.
 */

mzd_t *mzd_addmul_mp(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);


/**
 * \brief Matrix multiplication and in-place addition via cubic matrix
 *  multiplication algorithm on up to four cores, i.e. compute C = C+ AB.
 * 
 * \param C Product matrix
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for recursion.
 */

mzd_t *_mzd_addmul_mp4(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

/**
 * \brief Matrix multiplication via cubic matrix multiplication algorithm on up
 *  to four cores, i.e. compute C = C+ AB.
 * 
 * \param C Product matrix
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for recursion.
 */

mzd_t *_mzd_mul_mp4(mzd_t *C, mzd_t const *A, mzd_t const *B, int cutoff);

#endif //__M4RI_HAVE_OPENMP
