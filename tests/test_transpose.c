/*
 * test_transpose.c
 *
 * Application to test functionality of mzd_transpose.
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

#include <m4ri/m4ri.h>

#define NSIZE 27
int test_size[NSIZE] = {1, 3, 4, 7, 8, 11, 16, 17, 32, 40, 63, 64, 65, 80, 128, 160, 192, 240, 256, 512, 513, 769, 1000, 2000, 3000, 4000, 5000};

int test_transpose(int i) {
  int failure = 0;
  rci_t m     = test_size[i];
  printf("transpose m: %4d, n: ", m);
  for (int j = 0; j < NSIZE && !failure; ++j) {
    rci_t n = test_size[j];
    printf("%d", n);
    if (j != NSIZE - 1) printf(",");
    int size      = m * n;
    int loop_size = MAX(64 * 64 / size, 2);
    for (int i = 0; i < loop_size && !failure; ++i) {
      mzd_t *A = mzd_init(m, n);
      mzd_t *B = mzd_init(m, n);
      mzd_randomize(A);
      mzd_randomize(B);
      mzd_t *C  = mzd_add(NULL, A, B);
      mzd_t *AT = mzd_init(n, m);
      mzd_randomize(AT);
      mzd_transpose(AT, A);
      mzd_t *BT  = mzd_transpose(NULL, B);
      mzd_t *CT  = mzd_add(NULL, AT, BT);
      mzd_t *CTT = mzd_transpose(NULL, CT);
      if (!mzd_equal(C, CTT)) ++failure;

      // check ATT = A
      mzd_t *ATT = B;
      mzd_randomize(A);
      mzd_randomize(AT);
      mzd_randomize(ATT);
      mzd_transpose(AT, A);
      mzd_transpose(ATT, AT);
      if (!mzd_equal(A, ATT)) ++failure;

      mzd_free(A);
      mzd_free(B);
      mzd_free(C);
      mzd_free(AT);
      mzd_free(BT);
      mzd_free(CT);
      mzd_free(CTT);
    }
  }
  printf("  ");
  if (failure) {
    printf("FAILED\n");
  } else
    printf("passed\n");
  return failure;
}

int main() {
  int status = 0;

  for (int i = 0; i < NSIZE; ++i) {
      status += test_transpose(i);
  }

  if (!status) {
    printf("All tests passed.\n");
  } else {
    printf("TEST FAILED!\n");
    return 1;
  }

  return 0;
}
