/**
 * \file permutation.h
 *
 * This is scratch, experimental, unfinished code.
 *
 * \brief Permutation matrices for the M4RI library.
 * 
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 *
 * \internal
 */
/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de> 
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
******************************************************************************/
#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "misc.h"

/**
 * \brief Permutations.
 */

typedef struct {
  /**
   * The swap operations in LAPACK format.
   */
  size_t *values;

  /**
   * The length of the swap array.
   */

  size_t length;

} permutation;

/**
 * Construct an identity permutation.
 * 
 * \param length Length of the permutation.
 */

permutation *mzp_init(size_t length);

/**
 * Free a permutation.
 * 
 * \param P Permutation to free.
 */

void mzp_free(permutation *P);

#endif //PERMUTATION_H
