/**
 * \file pluq_mmpf.h
 * \brief PLUQ factorisation using Gray codes
 *
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */


#ifndef PLUQ_MMPF_H
#define PLUQ_MMPF_H
 /*******************************************************************
 *
 *                 M4RI:  Linear Algebra over GF(2)
 *
 *    Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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

#include "packedmatrix.h"
#include "permutation.h"

/**
 * Perform PLUQ factorization on A using Gray codes.
 *
 * \param A Matrix.
 * \param k Size of Gray code tables.
 * \param P Preallocated row permutation.
 * \param Q Preallocated column permutation.
 *
 * \wordoffset
 */

size_t _mzd_pluq_mmpf(packedmatrix *A, permutation * P, permutation * Q, int k);

/**
 * Perform PLUQ factorization on a submatrix of up to dimension k
 * starting at (r,c).
 *
 * \param A Matrix.
 * \param r Row Offset.
 * \param c Column Offset.
 * \param k Size of Gray code tables.
 * \param P Preallocated row permutation.
 * \param Q Preallocated column permutation.
 *
 * \wordoffset
 */

size_t _mzd_pluq_submatrix(packedmatrix *A, size_t r, size_t c, int k, permutation *P, permutation *Q);

#endif //PLUQ_MMPF_H
