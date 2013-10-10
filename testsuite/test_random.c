/*
 * test_random.c
 *
 * Application to test functionality of mzd_randomize (and mzd_equal).
 * In particular if these function correctly for windowed matrices.
 *
 * Copyright (C) 2011  Carlo Wood  <carlo@alinoe.com>
 * RSA-1024 0x624ACAD5 1997-01-26                    Sign & Encrypt
 * Fingerprint16 = 32 EC A7 B6 AC DB 65 A6  F6 F6 55 DD 1C DC FF 61
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <m4ri/config.h>
#include <stdio.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>

int test_random(rci_t m, rci_t n)
{
  mzd_t *A = mzd_init(m + 3, n + 64);
  mzd_t *W = mzd_init_window(A, 1, 0, m + 1, n);
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

  for (rci_t n = 0; n < 3 * m4ri_radix; n += m4ri_radix)
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
