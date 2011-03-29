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
#include <stdio.h>
#include <stdlib.h>
#include "m4ri.h"

int test_random(rci_t m, rci_t n)
{
  mzd_t *A = mzd_init(m + 3, n + 64);
  mzd_t *W = mzd_init_window(A, 1, 13, m + 1, n + 13);
  mzd_t *M = mzd_init(m, n);
  printf("randomize m: %4d, n: %4d ", m, n);
  srandom(17);
  mzd_randomize(M);
  srandom(17);
  mzd_randomize(W);
  int failure = !mzd_equal(M, W);
  if (failure)
  {
#if 1
    printf("FAILED\n");
#else
    printf("FAILURE: M != W:\n");
    printf("M %dx%d:\n", m, n);
    mzd_print(M);
    printf("W %dx%d:\n", m, n);
    mzd_print(W);
#endif
  }
  else
    printf("passed\n");
  mzd_free(M);
  mzd_free(A);
  return failure;
}

int main() {
  int status = 0;

  srandom(17);

  for (rci_t n = 0; n < 3 * RADIX; n += RADIX)
  {
    status += test_random(20, n + 1);
    status += test_random(20, n + 2);
    status += test_random(20, n + 32);
    status += test_random(20, n + 50);
    status += test_random(20, n + 51);
    status += test_random(20, n + 52);
    status += test_random(20, n + 63);
    status += test_random(20, n + 64);
    status += test_random(20, n + 65);
  }

  if (!status) {
    printf("All tests passed.\n");
  } else {
    printf("TEST FAILED!\n");
    return 1;
  }

  return 0;
}
