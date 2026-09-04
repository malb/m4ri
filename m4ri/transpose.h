/**
 * \file transpose.h
 *
 * \brief Matrix transpose using a variant of the cache-oblivious algorithm
 */

#ifndef M4RI_TRANSPOSE_H
#define M4RI_TRANSPOSE_H

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

#include <m4ri/mzd.h>

/**
 * \brief Transpose a matrix.
 *
 * This function uses the fact that:
\verbatim
   [ A B ]T    [AT CT]
   [ C D ]  =  [BT DT]
 \endverbatim
 * and thus rearranges the blocks recursively.
 *
 * \param DST Preallocated return matrix, may be NULL for automatic creation.
 * \param A Matrix
 */
mzd_t *mzd_transpose(mzd_t *DST, mzd_t const *A);

#endif  // M4RI_TRANSPOSE_H
