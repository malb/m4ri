/******************************************************************************
 *
 *            M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2011 Carlo Wood <carlo@alinoe.com>
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
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mzd.h"
#include "mzp.h"

#if __M4RI_DEBUG_DUMP

static unsigned long dd_sequence_number = 0;

static void entry(char const *function, char const *file, int line) {
#if !__M4RI_DD_QUIET
  printf("Sequence#: %ld; %s @ %s:%d; ", dd_sequence_number, function, file, line);
#endif
  ++dd_sequence_number;
}

static inline void consistency_check_row(mzd_t const *M, rci_t row) {
  assert(row >= 0 && row < M->nrows);
  assert(M->rows[row] == mzd_row(M, row));
  if (mzd_is_windowed(M)) return;
  // Check that the excess bits are zero.
  assert((M->rows[row][M->width - 1] & ~M->high_bitmask) == 0);
  // Check that the padding bits are zero, if any.
  assert(M->width == M->rowstride || M->rows[row][M->width] == 0);
}

static void consistency_check(mzd_t const *M) {
  assert(M->nrows >= 0 && M->ncols >= 0);
  assert(M->width * m4ri_radix >= M->ncols);
  assert((M->width - 1) * m4ri_radix < M->ncols);
  assert(M->width < mzd_paddingwidth || (M->rowstride & 1) == 0);
  // assert((M->blockrows_mask + 1) == (1 << M->blockrows_log));
  assert((1 << M->blockrows_log) * M->rowstride <= __M4RI_MAX_MZD_BLOCKSIZE);
  assert((1 << M->blockrows_log) * M->rowstride > __M4RI_MAX_MZD_BLOCKSIZE / 2);
  assert((M->width > 1 && M->high_bitmask == __M4RI_LEFT_BITMASK((M->ncols) % m4ri_radix)) ||
         (M->width < 2 && M->high_bitmask == __M4RI_MIDDLE_BITMASK(M->ncols, 0)));
  assert(((M->flags & mzd_flag_nonzero_excess) == 0) == ((M->ncols % m4ri_radix == 0)));
  assert((M->flags & mzd_flag_windowed_zeroexcess) == 0 || ((M->ncols) % m4ri_radix == 0));
  assert(
      (((M->flags & mzd_flag_multiple_blocks) == 0) == (mzd_row_to_block(M, M->nrows - 1) == 0)));
  int n         = 0;
  rci_t counted = 0;
  word *ptr     = mzd_first_row(M);
  int row_count = mzd_rows_in_block(M, 0);
  while (1) {
    while (row_count--) {
      assert(ptr == M->rows[counted++]);
      ptr += M->rowstride;
    }
    ++n;
    row_count = mzd_rows_in_block(M, n);
    if (row_count <= 0) break;
    ptr = mzd_first_row_next_block(M, n);
  }
  assert(M->ncols == 0 || counted == M->nrows);
  if (mzd_is_windowed(M)) return;
  assert(M->rowstride == M->width ||
         (M->rowstride == M->width + 1 && M->width >= mzd_paddingwidth));
  for (rci_t r = 0; r < M->nrows; ++r) { consistency_check_row(M, r); }
}

void m4ri_dd_int(char const *function, char const *file, int line, int i) {
  entry(function, file, line);
#if !__M4RI_DD_QUIET
  printf("int: %d\n", i);
#endif
}

void m4ri_dd_rci(char const *function, char const *file, int line, rci_t rci) {
  entry(function, file, line);
#if !__M4RI_DD_QUIET
  printf("rci: %d\n", rci);
#endif
}

void m4ri_dd_rci_array(char const *function, char const *file, int line, rci_t *rciptr, int len) {
  entry(function, file, line);
#if !__M4RI_DD_QUIET
  word hash = 0;
  for (int i = 0; i < len; ++i) hash ^= rotate_word(rciptr[i], i % m4ri_radix);
  printf("rci array (size %d) hash: %llx\n", len, hash);
#endif
}

void m4ri_dd_rawrow(char const *function, char const *file, int line, word const *rowptr,
                    wi_t wide) {
  entry(function, file, line);
#if !__M4RI_DD_QUIET
  word hash = calculate_hash(rowptr, wide);
  printf("raw row (%d words) hash: %llx\n", wide, hash);
#endif
}

void m4ri_dd_row(char const *function, char const *file, int line, mzd_t const *M, rci_t row) {
  entry(function, file, line);
  consistency_check_row(M, row);
#if !__M4RI_DD_QUIET
  word hash = calculate_hash(M->rows[row], M->width);
  printf("row %d hash: %llx\n", row, hash);
#endif
}

void m4ri_dd_mzd(char const *function, char const *file, int line, mzd_t const *M) {
  entry(function, file, line);
  consistency_check(M);
#if !__M4RI_DD_QUIET
  word hash = 0;
  for (rci_t r = 0; r < M->nrows; ++r)
    hash ^= rotate_word(calculate_hash(M->rows[r], M->width), r % m4ri_radix);
  printf("mzd hash: %llx\n", hash);
#endif
}

void m4ri_dd_mzp(char const *function, char const *file, int line, mzp_t const *P) {
  entry(function, file, line);
#if !__M4RI_DD_QUIET
  word hash = 0;
  for (rci_t i = 0; i < P->length; ++i) hash ^= rotate_word(P->values[i], i % m4ri_radix);
  printf("mzp hash: %llx\n", hash);
#endif
}

#endif
