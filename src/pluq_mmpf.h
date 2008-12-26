/**
 * \file pluq_mmpf.h
 * \brief PLUQ factorization using Gray codes.
 *
 *
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 * 
 * \example testsuite/test_lqup.c
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
 * \brief PLUQ matrix decomposition of A using Gray codes.
 *
 * If (P,L,U,Q) satisfy PLUQ = A, this function returns
 * (P^T,L,U,Q^T). The matrix L and U are stored in place over A.
 *
 * \param A Matrix.
 * \param P Preallocated row permutation.
 * \param Q Preallocated column permutation.
 * \param k Size of Gray code tables.
 *
 * \wordoffset
 */

size_t _mzd_pluq_mmpf(packedmatrix *A, permutation * P, permutation * Q, int k);

/**
 * \brief PLUQ matrix decomposition of a submatrix of up to dimension
 * k starting at (r,c).
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
