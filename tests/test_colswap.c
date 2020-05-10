/*
 * test_colswap.c
 *
 * Application to test functionality of mzd_col_swap.
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

int test_colswap(rci_t c1) {
  int failure = 0;
  printf("col_swap c1: %4d", c1);
  rci_t const rows = 4 * c1 + 1;
  mzd_t *A         = mzd_init(rows, 4 * c1 + 1);
  mzd_randomize(A);
  for (int c2 = 0; c2 < c1 && !failure; c2 += 13) {
    mzd_t *B = mzd_init_window(A, 0, 0, rows, c1 + 1);
    mzd_t *C = mzd_copy(NULL, B);

    mzd_col_swap(B, c1, c2);
    for (int r = 0; r < rows; ++r) {
      if (mzd_read_bit(C, r, c1) != mzd_read_bit(B, r, c2) ||
          mzd_read_bit(C, r, c2) != mzd_read_bit(B, r, c1)) {
        ++failure;
        break;
      }
    }
    mzd_col_swap(B, c2, c1);
    if (!mzd_equal(B, C)) { ++failure; }

    mzd_col_swap_in_rows(B, c1, c2, c2, c1);
    for (int r = c1; r < c2; ++r) {
      if (mzd_read_bit(C, r, c1) != mzd_read_bit(B, r, c2) ||
          mzd_read_bit(C, r, c2) != mzd_read_bit(B, r, c1)) {
        ++failure;
        break;
      }
    }
    mzd_col_swap_in_rows(B, c2, c1, c2, c1);
    if (!mzd_equal(B, C)) { ++failure; }

    mzd_free(C);
    mzd_free(B);
  }
  mzd_free(A);
  printf("  ");
  if (failure) {
    printf("FAILED\n");
  } else
    printf("passed\n");
  return failure;
}

int main() {
  int status = 0;

  for (int c1 = 1; c1 < 400; c1 += 15) { status += test_colswap(c1); }

  if (!status) {
    printf("All tests passed.\n");
  } else {
    printf("TEST FAILED!\n");
    return 1;
  }

  return 0;
}
