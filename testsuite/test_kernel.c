/*******************************************************************
*
*                M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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

#include "config.h"
#include <stdlib.h>
#include "m4ri.h"

int test_kernel_left_pluq(rci_t m, rci_t n) {
  mzd_t* A = mzd_init(m, n);
  mzd_randomize(A);
  
  mzd_t *Acopy = mzd_copy(NULL, A);

  rci_t r = mzd_echelonize_m4ri(A, 0, 0);
  printf("kernel_left m: %4zu, n: %4zu, r: %4zu ", m.val(), n.val(), r.val());
  mzd_free(Acopy);
  Acopy = mzd_copy(NULL, A);
    
  mzd_t *X = mzd_kernel_left_pluq(A, 0);
  if (X == NULL) {
    printf("passed\n");
    mzd_free(A);
    mzd_free(Acopy);
    return 0;
  }

  mzd_t *Z = mzd_mul(NULL, Acopy, X, 0);
  
  int status = 1 - mzd_is_zero(Z);
  
  if (!status)
    printf("passed\n");
  else
    printf("FAILED\n");

  mzd_free(A);
  mzd_free(Acopy);
  mzd_free(X);
  mzd_free(Z);
  return status;
}

int main(int argc, char **argv) {
  int status = 0;

  status += test_kernel_left_pluq(  2U,   4U);
  status += test_kernel_left_pluq(  4U,   1U);
  status += test_kernel_left_pluq( 10U,  20U);
  status += test_kernel_left_pluq( 20U,   1U);
  status += test_kernel_left_pluq( 20U,  20U);
  status += test_kernel_left_pluq( 30U,   1U);
  status += test_kernel_left_pluq( 30U,  30U);
  status += test_kernel_left_pluq( 80U,   1U);
  status += test_kernel_left_pluq( 80U,  20U);
  status += test_kernel_left_pluq( 80U,  80U);

  status += test_kernel_left_pluq( 4U,  2U);
  status += test_kernel_left_pluq( 1U,  4U);
  status += test_kernel_left_pluq(20U, 10U);
  status += test_kernel_left_pluq( 1U, 20U);
  status += test_kernel_left_pluq(20U, 20U);
  status += test_kernel_left_pluq( 1U, 30U);
  status += test_kernel_left_pluq(30U, 30U);
  status += test_kernel_left_pluq( 1U, 80U);
  status += test_kernel_left_pluq(20U, 80U);
  status += test_kernel_left_pluq(80U, 80U);

  status += test_kernel_left_pluq(10U, 20U);
  status += test_kernel_left_pluq(10U, 80U);
  status += test_kernel_left_pluq(10U, 20U);
  status += test_kernel_left_pluq(10U, 80U);
  status += test_kernel_left_pluq(70U, 20U);
  status += test_kernel_left_pluq(70U, 80U);
  status += test_kernel_left_pluq(70U, 20U);
  status += test_kernel_left_pluq(70U, 80U);
  status += test_kernel_left_pluq(770U, 1600U);
  status += test_kernel_left_pluq(1764U, 1345U);

  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
