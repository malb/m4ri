/**
 * \file strassen.h
 * \brief Matrix operations using Strassen's formulas.
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
 *       Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
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

#include "misc.h"
#include "packedmatrix.h"
#include "brilliantrussian.h"

/**
 * \brief Matrix multiplication via Strassen's matrix multiplication formula,
 * i.e. compute C = AB. 
 * 
 * This is the convenient wrapper function.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 */

packedmatrix *mzd_mul_strassen(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

/**
 * \brief Matrix multiplication via Strassen's matrix multiplication
 * formula, i.e. compute C = AB. 
 * 
 * This is the actual implementation.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A
 * \param B Input matrix B
 * \param cutoff Minimal dimension for Strassen recursion.
 *
 * \note This implement is to a very large extend a port of the
 * function strassen_window_multiply_c in Sage 3.0; For reference see
 * http://www.sagemath.org
 */

packedmatrix *_mzd_mul_strassen_impl(packedmatrix *C, packedmatrix *A, packedmatrix *B, int cutoff);

#endif //STRASSEN_H
