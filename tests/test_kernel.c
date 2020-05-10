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

#include "testing.h"
#include <m4ri/config.h>
#include <m4ri/m4ri.h>
#include <stdlib.h>

int test_kernel_left_pluq(rci_t m, rci_t n) {
  mzd_t *A = mzd_init(m, n);
  mzd_randomize(A);

  mzd_t *Acopy = mzd_copy(NULL, A);

  rci_t r = mzd_echelonize_m4ri(A, 0, 0);
  printf("kernel_left m: %4d, n: %4d, r: %4d ", m, n, r);
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

int main() {
  int status = 0;

  srandom(17);

  status += test_kernel_left_pluq(2, 4);
  status += test_kernel_left_pluq(4, 1);
  status += test_kernel_left_pluq(10, 20);
  status += test_kernel_left_pluq(20, 1);
  status += test_kernel_left_pluq(20, 20);
  status += test_kernel_left_pluq(30, 1);
  status += test_kernel_left_pluq(30, 30);
  status += test_kernel_left_pluq(80, 1);
  status += test_kernel_left_pluq(80, 20);
  status += test_kernel_left_pluq(80, 80);

  status += test_kernel_left_pluq(4, 2);
  status += test_kernel_left_pluq(1, 4);
  status += test_kernel_left_pluq(20, 10);
  status += test_kernel_left_pluq(1, 20);
  status += test_kernel_left_pluq(20, 20);
  status += test_kernel_left_pluq(1, 30);
  status += test_kernel_left_pluq(30, 30);
  status += test_kernel_left_pluq(1, 80);
  status += test_kernel_left_pluq(20, 80);
  status += test_kernel_left_pluq(80, 80);

  status += test_kernel_left_pluq(10, 20);
  status += test_kernel_left_pluq(10, 80);
  status += test_kernel_left_pluq(10, 20);
  status += test_kernel_left_pluq(10, 80);
  status += test_kernel_left_pluq(70, 20);
  status += test_kernel_left_pluq(70, 80);
  status += test_kernel_left_pluq(70, 20);
  status += test_kernel_left_pluq(70, 80);
  status += test_kernel_left_pluq(770, 1600);
  status += test_kernel_left_pluq(1764, 1345);

  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
