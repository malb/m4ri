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

#include "permutation.h"
#include "packedmatrix.h"

permutation *mzp_init(size_t length) {
  size_t i;
  permutation *P = (permutation*)m4ri_mm_malloc(sizeof(permutation));
  P->values = (size_t*)m4ri_mm_malloc(sizeof(size_t)*length);
  P->length = length;
  for (i=0; i<length; i++) {
    P->values[i] = i;
  }
  return P;
}

void mzp_free(permutation *P) {
  m4ri_mm_free(P->values);
  m4ri_mm_free(P);
}
