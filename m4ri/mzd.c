/******************************************************************************
 *
 *            M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
 *    Copyright (C) 2009-2013 Martin Albrecht <martinralbrecht+m4ri@googlemail.com>
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

#ifdef __M4RI_HAVE_LIBPNG
#include <png.h>
#endif

#include "mmc.h"
#include "mzd.h"
#include "parity.h"
#include <stdlib.h>
#include <string.h>

/**
 * \brief Cache of mzd_t containers
 */

typedef struct mzd_t_cache {
  mzd_t mzd[64];            /*!< cached matrices */
  struct mzd_t_cache *prev; /*!< previous block */
  struct mzd_t_cache *next; /*!< next block */
  uint64_t used;            /*!< bitmasks which matrices in this block are used */
  unsigned char padding[sizeof(mzd_t) - 2 * sizeof(struct mzd_t_cache *) -
                        sizeof(uint64_t)]; /*!< alignment */
#ifdef __GNUC__
} mzd_t_cache_t __attribute__((__aligned__(64)));
#else
} mzd_t_cache_t;
#endif

#define __M4RI_MZD_T_CACHE_MAX 16
static mzd_t_cache_t mzd_cache;
static mzd_t_cache_t *current_cache = &mzd_cache;

static int log2_floor(uint64_t v) {
  static uint64_t const b[]     = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000, 0xFFFFFFFF00000000};
  static unsigned int const S[] = {1, 2, 4, 8, 16, 32};
  unsigned int r                = 0;
  for (int i = 5; i >= 0; --i) {
    if ((v & b[i])) {
      v >>= S[i];
      r |= S[i];
    }
  }
  return r;
}

/*
 * Return a pointer to a new mzd_t structure.
 * The structure will be 64 byte aligned.
 * Call mzd_t_free to free the structure for next use.
 */

static mzd_t *mzd_t_malloc() {
#if __M4RI_ENABLE_MZD_CACHE == 0
  return (mzd_t *)m4ri_mm_malloc(sizeof(mzd_t));
#else
  mzd_t *ret = NULL;
  int i      = 0;

  if (current_cache->used == (uint64_t)-1) {
    mzd_t_cache_t *cache = &mzd_cache;
    while (cache && cache->used == (uint64_t)-1) {
      current_cache = cache;
      cache         = cache->next;
      i++;
    }
    if (!cache && i < __M4RI_MZD_T_CACHE_MAX) {
      cache = (mzd_t_cache_t *)m4ri_mm_malloc_aligned(sizeof(mzd_t_cache_t), 64);
      memset((char *)cache, 0, sizeof(mzd_t_cache_t));

      cache->prev         = current_cache;
      current_cache->next = cache;
      current_cache       = cache;
    } else if (!cache && i >= __M4RI_MZD_T_CACHE_MAX) {
      /* We have reached the upper limit on the number of caches */
      ret = (mzd_t *)m4ri_mm_malloc(sizeof(mzd_t));
    } else {
      current_cache = cache;
    }
  }
  if (ret == NULL) {
    int free_entry = log2_floor(~current_cache->used);
    current_cache->used |= ((uint64_t)1 << free_entry);
    ret = &current_cache->mzd[free_entry];
  }
  return ret;
#endif  //__M4RI_ENABLE_MZD_CACHE
}

static void mzd_t_free(mzd_t *M) {
#if __M4RI_ENABLE_MZD_CACHE == 0
  m4ri_mm_free(M);
#else
  int foundit          = 0;
  mzd_t_cache_t *cache = &mzd_cache;
  while (cache) {
    size_t entry = M - cache->mzd;
    if (entry < 64) {
      cache->used &= ~((uint64_t)1 << entry);
      if (cache->used == 0) {
        if (cache == &mzd_cache) {
          current_cache = cache;
        } else {
          if (cache == current_cache) { current_cache = cache->prev; }
          cache->prev->next = cache->next;
          if (cache->next) cache->next->prev = cache->prev;
          m4ri_mm_free(cache);
        }
      }
      foundit = 1;
      break;
    }
    cache = cache->next;
  }
  if (!foundit) { m4ri_mm_free(M); }
#endif  //__M4RI_ENABLE_MZD_CACHE
}

mzd_t *mzd_init(rci_t r, rci_t c) {
  assert(sizeof(mzd_t) == 64);
  mzd_t *A = mzd_t_malloc();
  A->nrows         = r;
  A->ncols         = c;
  A->width         = (c + m4ri_radix - 1) / m4ri_radix;
  A->rowstride     = ((A->width & 1) == 0) ? A->width : A->width + 1;
  A->high_bitmask  = __M4RI_LEFT_BITMASK(c % m4ri_radix);
  A->flags         = (A->high_bitmask != m4ri_ffff) ? mzd_flag_nonzero_excess : 0;
  if (r && c) {
    size_t block_words = r * A->rowstride;
    A->data = m4ri_mmc_calloc(block_words, sizeof(word));
  } else {
    A->data = NULL;
  }
  return A;
}

mzd_t *mzd_init_window(mzd_t *M, const rci_t lowr, const rci_t lowc, const rci_t highr,
                       const rci_t highc) {
  assert(lowc % m4ri_radix == 0);

  mzd_t *W = mzd_t_malloc();

  rci_t nrows = MIN(highr - lowr, M->nrows - lowr);
  rci_t ncols = highc - lowc;
  W->nrows = nrows;
  W->ncols = ncols;
  W->rowstride = M->rowstride;
  W->width = (ncols + m4ri_radix - 1) / m4ri_radix;
  W->high_bitmask = __M4RI_LEFT_BITMASK(ncols % m4ri_radix);
  W->flags = mzd_flag_windowed;
  if (ncols % m4ri_radix != 0) W->flags |= mzd_flag_nonzero_excess;
  W->data = M->data + lowr * M->rowstride + (lowc / m4ri_radix);
  __M4RI_DD_MZD(W);
  return W;
}

void mzd_free(mzd_t *A) {
  if (!mzd_is_windowed(A)) {
    size_t block_words = A->nrows * A->rowstride;
    m4ri_mmc_free(A->data, block_words * sizeof(word));
  }
  mzd_t_free(A);
}

void mzd_row_add(mzd_t *M, rci_t sourcerow, rci_t destrow) {
  mzd_row_add_offset(M, destrow, sourcerow, 0);
}

void mzd_row_clear_offset(mzd_t *M, rci_t row, rci_t coloffset) {
  wi_t const startblock = coloffset / m4ri_radix;
  word temp;
  word *truerow = mzd_row(M, row);
  /* make sure to start clearing at coloffset */
  if (coloffset % m4ri_radix) {
    temp = truerow[startblock];
    temp &= __M4RI_RIGHT_BITMASK(m4ri_radix - coloffset);
  } else {
    temp = 0;
  }
  truerow[startblock] = temp;
  for (wi_t i = startblock + 1; i < M->width; ++i) { truerow[i] = 0; }

  __M4RI_DD_ROW(M, row);
}

rci_t mzd_gauss_delayed(mzd_t *M, rci_t startcol, int full) {
  rci_t startrow = startcol;
  rci_t pivots   = 0;
  for (rci_t i = startcol; i < M->ncols; ++i) {
    for (rci_t j = startrow; j < M->nrows; ++j) {
      if (mzd_read_bit(M, j, i)) {
        mzd_row_swap(M, startrow, j);
        ++pivots;

        for (rci_t ii = full ? 0 : startrow + 1; ii < M->nrows; ++ii) {
          if (ii != startrow) {
            if (mzd_read_bit(M, ii, i)) { mzd_row_add_offset(M, ii, startrow, i); }
          }
        }
        startrow = startrow + 1;
        break;
      }
    }
  }

  __M4RI_DD_MZD(M);
  __M4RI_DD_RCI(pivots);
  return pivots;
}

rci_t mzd_echelonize_naive(mzd_t *M, int full) { return mzd_gauss_delayed(M, 0, full); }

/**
 * Transpose a 64 x 64 matrix with width 1.
 *
 * \param dst First word of destination matrix.
 * \param src First word of source matrix.
 * \param rowstride_dst Rowstride of matrix dst.
 * \param rowstride_src Rowstride of matrix src.
 *
 * Rows of both matrices are expected to fit exactly in a word (offset == 0)
 * and lay entirely inside a single block.
 *
 * \note This function also works when dst == src.
 */

static inline void _mzd_copy_transpose_64x64(word *dst, word const *src, wi_t rowstride_dst,
                                             wi_t rowstride_src) {
  /*
   * m runs over the values:
   *   0x00000000FFFFFFFF
   *   0x0000FFFF0000FFFF
   *   0x00FF00FF00FF00FF
   *   0x0F0F0F0F0F0F0F0F
   *   0x3333333333333333
   *   0x5555555555555555,
   * alternating j zeroes with j ones.
   *
   * Assume we have a matrix existing of four jxj matrices ((0,0) is in the top-right corner,
   * this is the memory-model view, see the layout on
   * http://m4ri.sagemath.org/doxygen/structmzd__t.html):
   * ...[A1][B1][A0][B0]
   * ...[C1][D1][C0][D0]
   *          . [A2][B2]
   *        .   [C2][B2]
   *      .         .
   *                .
   * The following calulates the XOR between A and D,
   * and subsequently applies that to A and D respectively,
   * swapping A and D as a result.
   * Therefore wk starts at the first row and then has rowstride
   * added j times, running over the rows of A, then skips C
   * by adding j * rowstride to continue with the next A below C.
   */

  word m               = __M4RI_CONVERT_TO_WORD(0xFFFFFFFF);
  wi_t j_rowstride_dst = rowstride_dst * 64;
  wi_t j_rowstride_src = rowstride_src * 32;
  word *const end      = dst + j_rowstride_dst;
  // We start with j = 32, and a one-time unrolled loop, where
  // we copy from src and write the result to dst, swapping
  // the two 32x32 corner matrices.
  int j = 32;
  j_rowstride_dst >>= 1;
  word *RESTRICT wk = dst;
  for (word const *RESTRICT wks = src; wk < end; wk += j_rowstride_dst, wks += j_rowstride_src) {
    for (int k = 0; k < j; ++k, wk += rowstride_dst, wks += rowstride_src) {
      word xor                = ((*wks >> j) ^ *(wks + j_rowstride_src)) & m;
      *wk                     = *wks ^ (xor << j);
      *(wk + j_rowstride_dst) = *(wks + j_rowstride_src) ^ xor;
    }
  }
  // Next we work in-place in dst and swap the corners of
  // each of the last matrices, all in parallel, for all
  // remaining values of j.
  m ^= m << 16;
  for (j = 16; j != 0; j = j >> 1, m ^= m << j) {
    j_rowstride_dst >>= 1;
    for (wk = dst; wk < end; wk += j_rowstride_dst) {
      for (int k = 0; k < j; ++k, wk += rowstride_dst) {
        word xor = ((*wk >> j) ^ *(wk + j_rowstride_dst)) & m;
        *wk ^= xor << j;
        *(wk + j_rowstride_dst) ^= xor;
      }
    }
  }
}

/**
 * Transpose two 64 x 64 matrix with width 1.
 *
 * \param dst1 First word of destination matrix 1.
 * \param dst2 First word of destination matrix 2.
 * \param src1 First word of source matrix 1.
 * \param src2 First word of source matrix 2.
 * \param rowstride_dst Rowstride of destination matrices.
 * \param rowstride_src Rowstride of source matrices.
 *
 * Rows of all matrices are expected to fit exactly in a word (offset == 0)
 * and lay entirely inside a single block.
 *
 * \note This function also works to transpose in-place.
 */

static inline void _mzd_copy_transpose_64x64_2(word *RESTRICT dst1, word *RESTRICT dst2,
                                               word const *RESTRICT src1, word const *RESTRICT src2,
                                               wi_t rowstride_dst, wi_t rowstride_src) {
  word m               = __M4RI_CONVERT_TO_WORD(0xFFFFFFFF);
  wi_t j_rowstride_dst = rowstride_dst * 64;
  wi_t j_rowstride_src = rowstride_src * 32;
  word *const end      = dst1 + j_rowstride_dst;
  int j                = 32;
  word *RESTRICT wk[2];
  word const *RESTRICT wks[2];
  word xor [2];

  j_rowstride_dst >>= 1;
  wk[0]  = dst1;
  wk[1]  = dst2;
  wks[0] = src1;
  wks[1] = src2;

  do {

    for (int k = 0; k < j; ++k) {
      xor[0]                     = ((*wks[0] >> j) ^ *(wks[0] + j_rowstride_src)) & m;
      xor[1]                     = ((*wks[1] >> j) ^ *(wks[1] + j_rowstride_src)) & m;
      *wk[0]                     = *wks[0] ^ (xor[0] << j);
      *wk[1]                     = *wks[1] ^ (xor[1] << j);
      *(wk[0] + j_rowstride_dst) = *(wks[0] + j_rowstride_src) ^ xor[0];
      *(wk[1] + j_rowstride_dst) = *(wks[1] + j_rowstride_src) ^ xor[1];
      wk[0] += rowstride_dst;
      wk[1] += rowstride_dst;
      wks[0] += rowstride_src;
      wks[1] += rowstride_src;
    }

    wk[0] += j_rowstride_dst;
    wk[1] += j_rowstride_dst;
    wks[0] += j_rowstride_src;
    wks[1] += j_rowstride_src;

  } while (wk[0] < end);

  m ^= m << 16;
  for (j = 16; j != 0; j = j >> 1, m ^= m << j) {

    j_rowstride_dst >>= 1;
    wk[0] = dst1;
    wk[1] = dst2;

    do {

      for (int k = 0; k < j; ++k) {
        xor[0] = ((*wk[0] >> j) ^ *(wk[0] + j_rowstride_dst)) & m;
        xor[1] = ((*wk[1] >> j) ^ *(wk[1] + j_rowstride_dst)) & m;
        *wk[0] ^= xor[0] << j;
        *wk[1] ^= xor[1] << j;
        *(wk[0] + j_rowstride_dst) ^= xor[0];
        *(wk[1] + j_rowstride_dst) ^= xor[1];
        wk[0] += rowstride_dst;
        wk[1] += rowstride_dst;
      }

      wk[0] += j_rowstride_dst;
      wk[1] += j_rowstride_dst;

    } while (wk[0] < end);
  }
}

static unsigned char log2_ceil_table[64] = {
    0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};

static inline int log2_ceil(int n) { return log2_ceil_table[n - 1]; }

static word const transpose_mask[6] = {
    0x5555555555555555ULL, 0x3333333333333333ULL, 0x0F0F0F0F0F0F0F0FULL,
    0x00FF00FF00FF00FFULL, 0x0000FFFF0000FFFFULL, 0x00000000FFFFFFFFULL,
};

/**
 * Transpose 64/j matrices of size jxj in parallel.
 *
 * Where j equals n rounded up to the nearest power of 2.
 * The input array t must be of size j (containing the rows i of all matrices in t[i]).
 *
 * t[0..{j-1}]  = [Al]...[A1][A0]
 *
 * \param t An array of j words.
 * \param n The number of rows in each matrix.
 *
 * \return log2(j)
 */

static inline int _mzd_transpose_Nxjx64(word *RESTRICT t, int n) {
  int j  = 1;
  int mi = 0;  // Index into the transpose_mask array.

  while (j < n)  // Don't swap with entirely undefined data (where [D] exists entirely of
                 // non-existant rows).
  {
    // Swap 64/j matrices of size jxj in 2j rows. Thus,
    // <---- one word --->
    // [Al][Bl]...[A0][B0]
    // [Cl][Dl]...[C0][D0], where l = 64/j - 1 and each matrix [A], [B] etc is jxj.
    // Then swap [A] and [D] in-place.

    // m runs over the values in transpose_mask, so that at all
    // times m exists of j zeroes followed by j ones, repeated.
    word const m = transpose_mask[mi];
    int k        = 0;  // Index into t[].
    do {
      // Run over all rows of [A] and [D].
      for (int i = 0; i < j; ++i, ++k) {
        // t[k] contains row i of all [A], and t[k + j] contains row i of all [D]. Swap them.
        word xor = ((t[k] >> j) ^ t[k + j]) & m;
        t[k] ^= xor << j;
        t[k + j] ^= xor;
      }
      k += j;         // Skip [C].
    } while (k < n);  // Stop if we passed all valid input.

    // Double the size of j and repeat this for the next 2j rows until all
    // n rows have been swapped (possibly with non-existant rows).
    j <<= 1;
    ++mi;
  }

  return mi;
}

/**
 * Transpose a n x 64 matrix with width 1.
 *
 * \param dst First word of destination matrix.
 * \param src First word of source matrix.
 * \param rowstride_dst Rowstride of destination matrix.
 * \param rowstride_src Rowstride of source matrix.
 * \param n Number of rows in source matrix, must be less than 64.
 *
 * Rows of all matrices are expected have offset zero
 * and lay entirely inside a single block.
 *
 * \note This function also works to transpose in-place.
 */

static inline void _mzd_copy_transpose_lt64x64(word *RESTRICT dst, word const *RESTRICT src,
                                               wi_t rowstride_dst, wi_t rowstride_src, int n) {
  // Preload the n input rows into level 1, using a minimum of cache lines (compact storage).
  word t[64];
  word const *RESTRICT wks = src;
  int k;
  for (k = 0; k < n; ++k) {
    t[k] = *wks;
    wks += rowstride_src;
  }
  // see https://bitbucket.org/malb/m4ri/issues/53
  for (; k < 64; ++k) { t[k] = 0; }
  if (n > 32) {
    while (k < 64) t[k++] = 0;
    _mzd_copy_transpose_64x64(dst, t, rowstride_dst, 1);
    return;
  }
  int log2j = _mzd_transpose_Nxjx64(t, n);
  // All output bits are now transposed, but still might need to be shifted in place.
  // What we have now is 64/j matrices of size jxj. Thus,
  // [Al]...[A1][A0], where l = 64/j - 1.
  // while the actual output is:
  // [A0]
  // [A1]
  // ...
  // [Al]
  word const m      = __M4RI_LEFT_BITMASK(n);
  word *RESTRICT wk = dst;
  switch (log2j) {
  case 5: {
    wi_t const j_rowstride_dst = 32 * rowstride_dst;
    for (int k = 0; k < 32; ++k) {
      wk[0]               = t[k] & m;
      wk[j_rowstride_dst] = (t[k] >> 32) & m;
      wk += rowstride_dst;
    }
    break;
  }
  case 4: {
    wi_t const j_rowstride_dst = 16 * rowstride_dst;
    for (int k = 0; k < 16; ++k) {
      wk[0]                   = t[k] & m;
      wk[j_rowstride_dst]     = (t[k] >> 16) & m;
      wk[2 * j_rowstride_dst] = (t[k] >> 32) & m;
      wk[3 * j_rowstride_dst] = (t[k] >> 48) & m;
      wk += rowstride_dst;
    }
    break;
  }
  case 3: {
    wi_t const j_rowstride_dst = 8 * rowstride_dst;
    for (int k = 0; k < 8; ++k) {
      wk[0]                   = t[k] & m;
      wk[j_rowstride_dst]     = (t[k] >> 8) & m;
      wk[2 * j_rowstride_dst] = (t[k] >> 16) & m;
      wk[3 * j_rowstride_dst] = (t[k] >> 24) & m;
      wk[4 * j_rowstride_dst] = (t[k] >> 32) & m;
      wk[5 * j_rowstride_dst] = (t[k] >> 40) & m;
      wk[6 * j_rowstride_dst] = (t[k] >> 48) & m;
      wk[7 * j_rowstride_dst] = (t[k] >> 56) & m;
      wk += rowstride_dst;
    }
    break;
  }
  case 2: {
    wi_t const j_rowstride_dst = 4 * rowstride_dst;
    for (int k = 0; k < 4; ++k) {
      word *RESTRICT wk2 = wk;
      word tk            = t[k];
      for (int i = 0; i < 2; ++i) {
        wk2[0]                   = tk & m;
        wk2[j_rowstride_dst]     = (tk >> 4) & m;
        wk2[2 * j_rowstride_dst] = (tk >> 8) & m;
        wk2[3 * j_rowstride_dst] = (tk >> 12) & m;
        wk2[4 * j_rowstride_dst] = (tk >> 16) & m;
        wk2[5 * j_rowstride_dst] = (tk >> 20) & m;
        wk2[6 * j_rowstride_dst] = (tk >> 24) & m;
        wk2[7 * j_rowstride_dst] = (tk >> 28) & m;
        wk2 += 8 * j_rowstride_dst;
        tk >>= 32;
      }
      wk += rowstride_dst;
    }
    break;
  }
  case 1: {
    wi_t const j_rowstride_dst = 2 * rowstride_dst;
    for (int k = 0; k < 2; ++k) {
      word *RESTRICT wk2 = wk;
      word tk            = t[k];
      for (int i = 0; i < 8; ++i) {
        wk2[0]                   = tk & m;
        wk2[j_rowstride_dst]     = (tk >> 2) & m;
        wk2[2 * j_rowstride_dst] = (tk >> 4) & m;
        wk2[3 * j_rowstride_dst] = (tk >> 6) & m;
        wk2 += 4 * j_rowstride_dst;
        tk >>= 8;
      }
      wk += rowstride_dst;
    }
    break;
  }
  case 0: {
    word *RESTRICT wk2 = wk;
    word tk            = t[0];
    for (int i = 0; i < 16; ++i) {
      wk2[0]                 = tk & m;
      wk2[rowstride_dst]     = (tk >> 1) & m;
      wk2[2 * rowstride_dst] = (tk >> 2) & m;
      wk2[3 * rowstride_dst] = (tk >> 3) & m;
      wk2 += 4 * rowstride_dst;
      tk >>= 4;
    }
    break;
  }
  }
}

/**
 * Transpose a 64 x n matrix with width 1.
 *
 * \param dst First word of destination matrix.
 * \param src First word of source matrix.
 * \param rowstride_dst Rowstride of destination matrix.
 * \param rowstride_src Rowstride of source matrix.
 * \param n Number of columns in source matrix, must be less than 64.
 *
 * Rows of all matrices are expected have offset zero
 * and lay entirely inside a single block.
 *
 * \note This function also works to transpose in-place.
 */

static inline void _mzd_copy_transpose_64xlt64(word *RESTRICT dst, word const *RESTRICT src,
                                               wi_t rowstride_dst, wi_t rowstride_src, int n) {
  word t[64];
  int log2j                = log2_ceil(n);
  word const *RESTRICT wks = src;
  switch (log2j) {
  case 6: {
    _mzd_copy_transpose_64x64(t, src, 1, rowstride_src);
    word *RESTRICT wk = dst;
    for (int k = 0; k < n; ++k) {
      *wk = t[k];
      wk += rowstride_dst;
    }
    return;
  }
  case 5: {
    wi_t const j_rowstride_src = 32 * rowstride_src;
    for (int k = 0; k < 32; ++k) {
      t[k] = wks[0] | (wks[j_rowstride_src] << 32);
      wks += rowstride_src;
    }
    break;
  }
  case 4: {
    wi_t const j_rowstride_src = 16 * rowstride_src;
    for (int k = 0; k < 16; ++k) {
      t[k] = wks[0] | (wks[j_rowstride_src] << 16);
      t[k] |= (wks[2 * j_rowstride_src] << 32) | (wks[3 * j_rowstride_src] << 48);
      wks += rowstride_src;
    }
    break;
  }
  case 3: {
    wi_t const j_rowstride_src = 8 * rowstride_src;
    word tt;
    for (int k = 0; k < 8; ++k) {
      tt   = wks[0] | (wks[j_rowstride_src] << 8);
      t[k] = (wks[2 * j_rowstride_src] << 16) | (wks[3 * j_rowstride_src] << 24);
      tt |= (wks[4 * j_rowstride_src] << 32) | (wks[5 * j_rowstride_src] << 40);
      t[k] |= (wks[6 * j_rowstride_src] << 48) | (wks[7 * j_rowstride_src] << 56);
      wks += rowstride_src;
      t[k] |= tt;
    }
    break;
  }
  case 2: {
    word const *RESTRICT wks2 = wks + 60 * rowstride_src;
    t[0]                      = wks2[0];
    t[1]                      = wks2[rowstride_src];
    t[2]                      = wks2[2 * rowstride_src];
    t[3]                      = wks2[3 * rowstride_src];
    for (int i = 0; i < 15; ++i) {
      wks2 -= 4 * rowstride_src;
      t[0] <<= 4;
      t[1] <<= 4;
      t[2] <<= 4;
      t[3] <<= 4;
      t[0] |= wks2[0];
      t[1] |= wks2[rowstride_src];
      t[2] |= wks2[2 * rowstride_src];
      t[3] |= wks2[3 * rowstride_src];
    }
    break;
  }
  case 1: {
    wks += 62 * rowstride_src;
    t[0] = wks[0];
    t[1] = wks[rowstride_src];
    for (int i = 0; i < 31; ++i) {
      wks -= 2 * rowstride_src;
      t[0] <<= 2;
      t[1] <<= 2;
      t[0] |= wks[0];
      t[1] |= wks[rowstride_src];
    }
    break;
  }
  case 0: {
    word tt[2];
    tt[0] = wks[0];
    tt[1] = wks[rowstride_src];
    for (int i = 2; i < 64; i += 2) {
      wks += 2 * rowstride_src;
      tt[0] |= wks[0] << i;
      tt[1] |= wks[rowstride_src] << i;
    }
    *dst = tt[0] | (tt[1] << 1);
    return;
  }
  }
  int j = 1 << log2j;
  _mzd_transpose_Nxjx64(t, j);
  word *RESTRICT wk = dst;
  for (int k = 0; k < n; ++k) {
    *wk = t[k];
    wk += rowstride_dst;
  }
}

/**
 * Transpose a n x m matrix with width 1, offset 0 and m and n less than or equal 8.
 *
 * \param dst First word of destination matrix.
 * \param src First word of source matrix.
 * \param rowstride_dst Rowstride of destination matrix.
 * \param rowstride_src Rowstride of source matrix.
 * \param n Number of rows in source matrix, must be less than or equal 8.
 * \param m Number of columns in source matrix, must be less than or equal 8.
 *
 * Rows of all matrices are expected to have offset zero
 * and lay entirely inside a single block.
 *
 * \note This function also works to transpose in-place.
 */

static inline void _mzd_copy_transpose_le8xle8(word *RESTRICT dst, word const *RESTRICT src,
                                               wi_t rowstride_dst, wi_t rowstride_src, int n, int m,
                                               int maxsize) {
  int end                  = maxsize * 7;
  word const *RESTRICT wks = src;
  word w                   = *wks;
  int shift                = 0;
  for (int i = 1; i < n; ++i) {
    wks += rowstride_src;
    shift += 8;
    w |= (*wks << shift);
  }
  word mask = 0x80402010080402ULL;
  word w7   = w >> 7;
  shift     = 7;
  --m;
  do {
    word xor = (w ^ w7) & mask;
    mask >>= 8;
    w ^= (xor << shift);
    shift += 7;
    w7 >>= 7;
    w ^= xor;
  } while (shift < end);
  word *RESTRICT wk = dst + m * rowstride_dst;
  for (int shift = 8 * m; shift > 0; shift -= 8) {
    *wk = (unsigned char)(w >> shift);
    wk -= rowstride_dst;
  }
  *wk = (unsigned char)w;
}

/**
 * Transpose a n x m matrix with width 1, offset 0 and m and n less than or equal 16.
 *
 * \param dst First word of destination matrix.
 * \param src First word of source matrix.
 * \param rowstride_dst Rowstride of destination matrix.
 * \param rowstride_src Rowstride of source matrix.
 * \param n Number of rows in source matrix, must be less than or equal 16.
 * \param m Number of columns in source matrix, must be less than or equal 16.
 *
 * Rows of all matrices are expected to have offset zero
 * and lay entirely inside a single block.
 *
 * \note This function also works to transpose in-place.
 */

static inline void _mzd_copy_transpose_le16xle16(word *RESTRICT dst, word const *RESTRICT src,
                                                 wi_t rowstride_dst, wi_t rowstride_src, int n,
                                                 int m, int maxsize) {
  int end                  = maxsize * 3;
  word const *RESTRICT wks = src;
  word t[4];
  int i = n;
  do {
    t[0] = wks[0];
    if (--i == 0) {
      t[1] = 0;
      t[2] = 0;
      t[3] = 0;
      break;
    }
    t[1] = wks[rowstride_src];
    if (--i == 0) {
      t[2] = 0;
      t[3] = 0;
      break;
    }
    t[2] = wks[2 * rowstride_src];
    if (--i == 0) {
      t[3] = 0;
      break;
    }
    t[3] = wks[3 * rowstride_src];
    if (--i == 0) break;
    wks += 4 * rowstride_src;
    for (int shift = 16;; shift += 16) {
      t[0] |= (*wks << shift);
      if (--i == 0) break;
      t[1] |= (wks[rowstride_src] << shift);
      if (--i == 0) break;
      t[2] |= (wks[2 * rowstride_src] << shift);
      if (--i == 0) break;
      t[3] |= (wks[3 * rowstride_src] << shift);
      if (--i == 0) break;
      wks += 4 * rowstride_src;
    }
  } while (0);
  word mask = 0xF0000F0000F0ULL;
  int shift = 12;
  word xor [4];
  do {
    xor[0] = (t[0] ^ (t[0] >> shift)) & mask;
    xor[1] = (t[1] ^ (t[1] >> shift)) & mask;
    xor[2] = (t[2] ^ (t[2] >> shift)) & mask;
    xor[3] = (t[3] ^ (t[3] >> shift)) & mask;
    mask >>= 16;
    t[0] ^= (xor[0] << shift);
    t[1] ^= (xor[1] << shift);
    t[2] ^= (xor[2] << shift);
    t[3] ^= (xor[3] << shift);
    shift += 12;
    t[0] ^= xor[0];
    t[1] ^= xor[1];
    t[2] ^= xor[2];
    t[3] ^= xor[3];
  } while (shift < end);
  _mzd_transpose_Nxjx64(t, 4);
  i                 = m;
  word *RESTRICT wk = dst;
  do {
    wk[0] = (uint16_t)t[0];
    if (--i == 0) break;
    wk[rowstride_dst] = (uint16_t)t[1];
    if (--i == 0) break;
    wk[2 * rowstride_dst] = (uint16_t)t[2];
    if (--i == 0) break;
    wk[3 * rowstride_dst] = (uint16_t)t[3];
    if (--i == 0) break;
    wk += 4 * rowstride_dst;
    for (int shift = 16;; shift += 16) {
      wk[0] = (uint16_t)(t[0] >> shift);
      if (--i == 0) break;
      wk[rowstride_dst] = (uint16_t)(t[1] >> shift);
      if (--i == 0) break;
      wk[2 * rowstride_dst] = (uint16_t)(t[2] >> shift);
      if (--i == 0) break;
      wk[3 * rowstride_dst] = (uint16_t)(t[3] >> shift);
      if (--i == 0) break;
      wk += 4 * rowstride_dst;
    }
  } while (0);
}

/**
 * Transpose a n x m matrix with width 1, offset 0 and m and n less than or equal 32.
 *
 * \param dst First word of destination matrix.
 * \param src First word of source matrix.
 * \param rowstride_dst Rowstride of destination matrix.
 * \param rowstride_src Rowstride of source matrix.
 * \param n Number of rows in source matrix, must be less than or equal 32.
 * \param m Number of columns in source matrix, must be less than or equal 32.
 *
 * Rows of all matrices are expected to have offset zero
 * and lay entirely inside a single block.
 *
 * \note This function also works to transpose in-place.
 */

static inline void _mzd_copy_transpose_le32xle32(word *RESTRICT dst, word const *RESTRICT src,
                                                 wi_t rowstride_dst, wi_t rowstride_src, int n,
                                                 int m) {
  word const *RESTRICT wks = src;
  word t[16];
  int i = n;
  if (n > 16) {
    i -= 16;
    for (int j = 0; j < 16; ++j) {
      t[j] = *wks;
      wks += rowstride_src;
    }
    int j = 0;
    do {
      t[j++] |= (*wks << 32);
      wks += rowstride_src;
    } while (--i);
  } else {
    int j;
    for (j = 0; j < n; ++j) {
      t[j] = *wks;
      wks += rowstride_src;
    }
    for (; j < 16; ++j) t[j] = 0;
  }
  _mzd_transpose_Nxjx64(t, 16);
  int one_more      = (m & 1);
  word *RESTRICT wk = dst;
  if (m > 16) {
    m -= 16;
    for (int j = 0; j < 16; j += 2) {
      *wk               = (t[j] & 0xFFFF) | ((t[j] >> 16) & 0xFFFF0000);
      wk[rowstride_dst] = (t[j + 1] & 0xFFFF) | ((t[j + 1] >> 16) & 0xFFFF0000);
      wk += 2 * rowstride_dst;
    }
    for (int j = 1; j < m; j += 2) {
      *wk               = ((t[j - 1] >> 16) & 0xFFFF) | ((t[j - 1] >> 32) & 0xFFFF0000);
      wk[rowstride_dst] = ((t[j] >> 16) & 0xFFFF) | ((t[j] >> 32) & 0xFFFF0000);
      wk += 2 * rowstride_dst;
    }
    if (one_more) { *wk = ((t[m - 1] >> 16) & 0xFFFF) | ((t[m - 1] >> 32) & 0xFFFF0000); }
  } else {
    for (int j = 1; j < m; j += 2) {
      *wk               = (t[j - 1] & 0xFFFF) | ((t[j - 1] >> 16) & 0xFFFF0000);
      wk[rowstride_dst] = (t[j] & 0xFFFF) | ((t[j] >> 16) & 0xFFFF0000);
      wk += 2 * rowstride_dst;
    }
    if (one_more) { *wk = (t[m - 1] & 0xFFFF) | ((t[m - 1] >> 16) & 0xFFFF0000); }
  }
}

static inline void _mzd_copy_transpose_le64xle64(word *RESTRICT dst, word const *RESTRICT src,
                                                 wi_t rowstride_dst, wi_t rowstride_src, int n,
                                                 int m) {
  word const *RESTRICT wks = src;
  word t[64];
  int k;
  for (k = 0; k < n; ++k) {
    t[k] = *wks;
    wks += rowstride_src;
  }
  while (k < 64) t[k++] = 0;
  _mzd_copy_transpose_64x64(t, t, 1, 1);
  word *RESTRICT wk = dst;
  for (int k = 0; k < m; ++k) {
    *wk = t[k];
    wk += rowstride_dst;
  }
  return;
}

static inline void _mzd_copy_transpose_small(word *RESTRICT fwd, word const *RESTRICT fws,
                                                 wi_t rowstride_dst, wi_t rowstride_src, rci_t nrows,
                                                 rci_t ncols, rci_t maxsize) {
  assert(maxsize < 64);
  if (maxsize <= 8) {
    _mzd_copy_transpose_le8xle8(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols, maxsize);
  } else if (maxsize <= 16) {
    _mzd_copy_transpose_le16xle16(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols, maxsize);
  } else if (maxsize <= 32) {
    _mzd_copy_transpose_le32xle32(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols);
  } else {
    _mzd_copy_transpose_le64xle64(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols);
  }
}


void _mzd_transpose_base(word *RESTRICT fwd, word const *RESTRICT fws, wi_t rowstride_dst, 
                            wi_t rowstride_src, rci_t nrows, rci_t ncols, rci_t maxsize) {
  assert(maxsize >= 64);
    // Note that this code is VERY sensitive. ANY change to _mzd_transpose can easily
    // reduce the speed for small matrices (up to 64x64) by 5 to 10%.   
    if (nrows >= 64) {
      /*
       * This is an interesting #if ...
       * I recommend to investigate the number of instructions, and the clocks per instruction,
       * as function of various sizes of the matrix (most likely especially the number of columns
       * (the size of a row) will have influence; also always use multiples of 64 or even 128),
       * for both cases below.
       *
       * To measure this run for example:
       *
       * ./bench_mzd -m 10 -x 10 -p PAPI_TOT_INS,PAPI_L1_TCM,PAPI_L2_TCM mzd_transpose 32000 32000
       * ./bench_mzd -m 10 -x 100 -p PAPI_TOT_INS,PAPI_L1_TCM,PAPI_L2_TCM mzd_transpose 128 10240
       * etc (increase -x for smaller sizes to get better accuracy).
       *
       * --Carlo Wood
       */
#if 1
      int js = ncols & nrows & 64;  // True if the total number of whole 64x64 matrices is odd.
      wi_t const rowstride_64_dst      = 64 * rowstride_dst;
      word *RESTRICT fwd_current       = fwd;
      word const *RESTRICT fws_current = fws;
      if (js) {
        js = 1;
        _mzd_copy_transpose_64x64(fwd, fws, rowstride_dst, rowstride_src);
        if ((nrows | ncols) == 64) {
          return;
        }
        fwd_current += rowstride_64_dst;
        ++fws_current;
      }
      rci_t const whole_64cols = ncols / 64;
      // The use of delayed and even, is to avoid calling _mzd_copy_transpose_64x64_2 twice.
      // This way it can be inlined without duplicating the amount of code that has to be loaded.
      word *RESTRICT fwd_delayed       = NULL;
      word const *RESTRICT fws_delayed = NULL;
      int even                         = 0;
      while (1) {
        for (int j = js; j < whole_64cols; ++j) {
          if (!even) {
            fwd_delayed = fwd_current;
            fws_delayed = fws_current;
          } else {
            _mzd_copy_transpose_64x64_2(fwd_delayed, fwd_current, fws_delayed, fws_current,
                                        rowstride_dst, rowstride_src);
          }
          fwd_current += rowstride_64_dst;
          ++fws_current;
          even = !even;
        }
        nrows -= 64;
        if (ncols % 64) {
          _mzd_copy_transpose_64xlt64(fwd + whole_64cols * rowstride_64_dst, fws + whole_64cols,
                                      rowstride_dst, rowstride_src, ncols % 64);
        }
        fwd += 1;
        fws += 64 * rowstride_src;
        if (nrows < 64) break;
        js          = 0;
        fws_current = fws;
        fwd_current = fwd;
      }
#else
      // The same as the above, but without using _mzd_copy_transpose_64x64_2.
      wi_t const rowstride_64_dst = 64 * DST->rowstride;
      rci_t const whole_64cols    = ncols / 64;
      assert(nrows >= 64);
      do {
        for (int j = 0; j < whole_64cols; ++j) {
          _mzd_copy_transpose_64x64(fwd + j * rowstride_64_dst, fws + j, DST->rowstride,
                                    A->rowstride);
        }
        nrows -= 64;
        if (ncols % 64) {
          _mzd_copy_transpose_64xlt64(fwd + whole_64cols * rowstride_64_dst, fws + whole_64cols,
                                      DST->rowstride, A->rowstride, ncols % 64);
        }
        fwd += 1;
        fws += 64 * A->rowstride;
      } while (nrows >= 64);
#endif
    }

    if (nrows == 0) {
      return;
    }

    // Transpose the remaining top rows. Now 0 < nrows < 64.

    while (ncols >= 64) {
      _mzd_copy_transpose_lt64x64(fwd, fws, rowstride_dst, rowstride_src, nrows);
      ncols -= 64;
      fwd += 64 * rowstride_dst;
      fws += 1;
    }

    if (ncols == 0) {
      return ;
    }

  maxsize = MAX(nrows, ncols);
  
  // Transpose the remaining corner. Now both 0 < nrows < 64 and 0 < ncols < 64.
  _mzd_copy_transpose_small(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols, maxsize);
}

/* return the smallest multiple of k larger than n/2 */
static inline rci_t split_round(rci_t n, rci_t k) {
  rci_t half = n / 2;
  return ((half + (k - 1)) / k) * k;
}

static void _mzd_transpose_notsmall(word *RESTRICT fwd, word const *RESTRICT fws, wi_t rowstride_dst, 
                            wi_t rowstride_src, rci_t nrows, rci_t ncols, rci_t maxsize) {
  assert(maxsize >= 64);

  if (maxsize <= 512) {  // just one big block
    _mzd_transpose_base(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols, maxsize);
  } else {
    rci_t large_size = split_round(maxsize, (maxsize <= 768) ? 64 : 512);
    wi_t offset = large_size / m4ri_radix;
      if (nrows >= ncols) {
        word const *RESTRICT fws_up = fws; 
        word const *RESTRICT fws_down = fws + large_size * rowstride_src;
        word *RESTRICT fwd_left = fwd;
        word *RESTRICT fwd_right = fwd + offset; 
        rci_t maxsize_up = MAX(large_size, ncols);
        rci_t maxsize_down = MAX(nrows - large_size, ncols);
        _mzd_transpose_notsmall(fwd_left, fws_up, rowstride_dst, rowstride_src, large_size, ncols, maxsize_up);
        _mzd_transpose_notsmall(fwd_right, fws_down, rowstride_dst, rowstride_src, nrows - large_size, ncols, maxsize_down);
      } else {
        word const *RESTRICT fws_left = fws; 
        word const *RESTRICT fws_right = fws + offset;
        word *RESTRICT fwd_up = fwd;
        word *RESTRICT fwd_down = fwd + large_size * rowstride_dst; 
        rci_t maxsize_left = MAX(nrows, large_size);
        rci_t maxsize_right = MAX(nrows, ncols - large_size);
        _mzd_transpose_notsmall(fwd_up, fws_left, rowstride_dst, rowstride_src, nrows, large_size, maxsize_left);
        _mzd_transpose_notsmall(fwd_down, fws_right, rowstride_dst, rowstride_src, nrows, ncols - large_size, maxsize_right);
    }
  }
}

static void _mzd_transpose(word *RESTRICT fwd, word const *RESTRICT fws, wi_t rowstride_dst, 
                            wi_t rowstride_src, rci_t nrows, rci_t ncols, rci_t maxsize) {
  // rationale: small blocks corresponds to the word size
  //            two big blocks fit in L1 cache (512 --> 8KB).
  
  if (maxsize < 64) {  // super-fast path for very small matrices
    _mzd_copy_transpose_small(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols, maxsize);
  } else {
    _mzd_transpose_notsmall(fwd, fws, rowstride_dst, rowstride_src, nrows, ncols, maxsize);
  }
}



mzd_t *mzd_transpose(mzd_t *DST, mzd_t const *A) {
  if (DST == NULL) {
    DST = mzd_init(A->ncols, A->nrows);
  } else if (__M4RI_UNLIKELY(DST->nrows != A->ncols || DST->ncols != A->nrows)) {
    m4ri_die("mzd_transpose: Wrong size for return matrix.\n");
  }

  if (A->nrows == 0 || A->ncols == 0)
    return mzd_copy(DST, A);

  rci_t maxsize = MAX(A->nrows, A->ncols);
  if (__M4RI_LIKELY(!mzd_is_dangerous_window(DST))) {
    _mzd_transpose(DST->data, A->data, DST->rowstride, A->rowstride, A->nrows, A->ncols, maxsize);
    return DST;
  }
  
  mzd_t *D = mzd_init(DST->nrows, DST->ncols);
  _mzd_transpose(D->data, A->data, D->rowstride, A->rowstride, A->nrows, A->ncols, maxsize);
  mzd_copy(DST, D);
  mzd_free(D);
  return DST;
}

mzd_t *mzd_mul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions.\n");
    }
  }
  if (B->ncols < m4ri_radix - 10) { /* this cutoff is rather arbitrary */
    mzd_t *BT = mzd_transpose(NULL, B);
    _mzd_mul_naive(C, A, BT, 1);
    mzd_free(BT);
  } else {
    _mzd_mul_va(C, A, B, 1);
  }
  return C;
}

mzd_t *mzd_addmul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  if (C->nrows != A->nrows || C->ncols != B->ncols) {
    m4ri_die("mzd_addmul_naive: Provided return matrix has wrong dimensions.\n");
  }

  if (B->ncols < m4ri_radix - 10) { /* this cutoff is rather arbitrary */
    mzd_t *BT = mzd_transpose(NULL, B);
    _mzd_mul_naive(C, A, BT, 0);
    mzd_free(BT);
  } else {
    _mzd_mul_va(C, A, B, 0);
  }
  return C;
}

mzd_t *_mzd_mul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B, const int clear) {
  wi_t eol;

  if (clear) {
    word const mask_end = C->high_bitmask;
    /* improves performance on x86_64 but is not cross plattform */
    /* asm __volatile__ (".p2align 4\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop"); */
    for (rci_t i = 0; i < C->nrows; ++i) {
      wi_t j = 0;
      word *row = mzd_row(C, i);
      for (; j < C->width - 1; ++j) { row[j] = 0; }
      row[j] &= ~mask_end;
    }
  }

  if (C->ncols % m4ri_radix) {
    eol = (C->width - 1);
  } else {
    eol = (C->width);
  }

  word parity[64];
  for (int i = 0; i < 64; ++i) { parity[i] = 0; }
  wi_t const wide     = A->width;
  int const blocksize = __M4RI_MUL_BLOCKSIZE;
  for (rci_t start = 0; start + blocksize <= C->nrows; start += blocksize) {
    for (rci_t i = start; i < start + blocksize; ++i) {
      word const *a = mzd_row_const(A, i);
      word *c = mzd_row(C, i);
      for (rci_t j = 0; j < m4ri_radix * eol; j += m4ri_radix) {
        for (int k = 0; k < m4ri_radix; ++k) {
          word const *b = mzd_row_const(B, j + k);
          parity[k] = a[0] & b[0];
          for (wi_t ii = wide - 1; ii >= 1; --ii) parity[k] ^= a[ii] & b[ii];
        }
        c[j / m4ri_radix] ^= m4ri_parity64(parity);
      }

      if (eol != C->width) {
        word const mask_end = C->high_bitmask;
        /* improves performance on x86_64 but is not cross plattform */
        /* asm __volatile__ (".p2align 4\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop");
         */
        for (int k = 0; k < (C->ncols % m4ri_radix); ++k) {
          word const *b = mzd_row_const(B, m4ri_radix * eol + k);
          parity[k] = a[0] & b[0];
          for (wi_t ii = 1; ii < A->width; ++ii) parity[k] ^= a[ii] & b[ii];
        }
        c[eol] ^= m4ri_parity64(parity) & mask_end;
      }
    }
  }

  for (rci_t i = C->nrows - (C->nrows % blocksize); i < C->nrows; ++i) {
    word const *a = mzd_row_const(A, i);
    word *c = mzd_row(C, i);
    for (rci_t j = 0; j < m4ri_radix * eol; j += m4ri_radix) {
      for (int k = 0; k < m4ri_radix; ++k) {
        word const *b = mzd_row_const(B, j + k);
        parity[k] = a[0] & b[0];
        for (wi_t ii = wide - 1; ii >= 1; --ii) parity[k] ^= a[ii] & b[ii];
      }
      c[j / m4ri_radix] ^= m4ri_parity64(parity);
    }

    if (eol != C->width) {
      word const mask_end = C->high_bitmask;
      /* improves performance on x86_64 but is not cross plattform */
      /* asm __volatile__ (".p2align 4\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop"); */
      for (int k = 0; k < (C->ncols % m4ri_radix); ++k) {
        word const *b = mzd_row_const(B, m4ri_radix * eol + k);
        parity[k] = a[0] & b[0];
        for (wi_t ii = 1; ii < A->width; ++ii) parity[k] ^= a[ii] & b[ii];
      }
      c[eol] ^= m4ri_parity64(parity) & mask_end;
    }
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *_mzd_mul_va(mzd_t *C, mzd_t const *v, mzd_t const *A, int const clear) {
  if (clear) mzd_set_ui(C, 0);

  rci_t const m = v->nrows;
  rci_t const n = v->ncols;

  for (rci_t i = 0; i < m; ++i)
    for (rci_t j = 0; j < n; ++j)
      if (mzd_read_bit(v, i, j)) mzd_combine(C, i, 0, C, i, 0, A, j, 0);

  __M4RI_DD_MZD(C);
  return C;
}

void mzd_randomize(mzd_t *A) {
  wi_t const width    = A->width - 1;
  word const mask_end = A->high_bitmask;
  for (rci_t i = 0; i < A->nrows; ++i) {
    word *row = mzd_row(A, i);
    for (wi_t j = 0; j < width; ++j) row[j] = m4ri_random_word();
    row[width] ^= (row[width] ^ m4ri_random_word()) & mask_end;
  }

  __M4RI_DD_MZD(A);
}

void mzd_randomize_custom(mzd_t *A, m4ri_random_callback rc, void *data) {
  wi_t const width    = A->width - 1;
  word const mask_end = A->high_bitmask;
  for (rci_t i = 0; i < A->nrows; ++i) {
    word *row = mzd_row(A, i);
    for (wi_t j = 0; j < width; ++j) row[j] = rc(data);
    row[width] ^= (row[width] ^ rc(data)) & mask_end;
  }

  __M4RI_DD_MZD(A);
}

void mzd_set_ui(mzd_t *A, unsigned int value) {
  word const mask_end = A->high_bitmask;

  for (rci_t i = 0; i < A->nrows; ++i) {
    word *row = mzd_row(A, i);
    for (wi_t j = 0; j < A->width - 1; ++j) row[j] = 0;
    row[A->width - 1] &= ~mask_end;
  }

  if (value % 2 == 0) {
    __M4RI_DD_MZD(A);
    return;
  }

  rci_t const stop = MIN(A->nrows, A->ncols);
  for (rci_t i = 0; i < stop; ++i) { mzd_write_bit(A, i, i, 1); }

  __M4RI_DD_MZD(A);
}

int mzd_equal(mzd_t const *A, mzd_t const *B) {
  if (A->nrows != B->nrows) return FALSE;
  if (A->ncols != B->ncols) return FALSE;
  if (A == B) return TRUE;

  wi_t Awidth = A->width - 1;
  word const mask_end = A->high_bitmask;

  for (rci_t i = 0; i < A->nrows; ++i) {
    word const *rowa = mzd_row_const(A, i);
    word const *rowb = mzd_row_const(B, i);
    for (wi_t j = 0; j < Awidth; ++j) {
      if (rowa[j] != rowb[j]) return FALSE;
    }
    if (((rowa[Awidth] ^ rowb[Awidth]) & mask_end)) return FALSE;
  }
  return TRUE;
}

int mzd_cmp(mzd_t const *A, mzd_t const *B) {
  if (A->nrows < B->nrows) return -1;
  if (B->nrows < A->nrows) return 1;
  if (A->ncols < B->ncols) return -1;
  if (B->ncols < A->ncols) return 1;

  const word mask_end = A->high_bitmask;
  const wi_t n        = A->width - 1;

  /* Columns with large index are "larger", but rows with small index
     are more important than with large index. */

  for (rci_t i = 0; i < A->nrows; i++) {
    word const *rowa = mzd_row_const(A, i);
    word const *rowb = mzd_row_const(B, i);
    if ((rowa[n] & mask_end) < (rowb[n] & mask_end))
      return -1;
    else if ((rowa[n] & mask_end) > (rowb[n] & mask_end))
      return 1;

    for (wi_t j = n - 1; j >= 0; j--) {
      if (rowa[j] < rowb[j])
        return -1;
      else if (rowa[j] > rowb[j])
        return 1;
    }
  }
  return 0;
}

mzd_t *mzd_copy(mzd_t *N, mzd_t const *P) {
  if (N == P) return N;

  if (N == NULL) {
    N = mzd_init(P->nrows, P->ncols);
  } else {
    if (N->nrows < P->nrows || N->ncols < P->ncols)
      m4ri_die("mzd_copy: Target matrix is too small.");
  }
  wi_t const wide = P->width - 1;
  word mask_end   = P->high_bitmask;
  for (rci_t i = 0; i < P->nrows; ++i) {
    word const *p_truerow = mzd_row_const(P, i);
    word *n_truerow = mzd_row(N, i);
    for (wi_t j = 0; j < wide; ++j) n_truerow[j] = p_truerow[j];
    n_truerow[wide] = (n_truerow[wide] & ~mask_end) | (p_truerow[wide] & mask_end);
  }
  __M4RI_DD_MZD(N);
  return N;
}

/* This is sometimes called augment */
mzd_t *mzd_concat(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  if (A->nrows != B->nrows) { m4ri_die("mzd_concat: Bad arguments to concat!\n"); }

  if (C == NULL) {
    C = mzd_init(A->nrows, A->ncols + B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != (A->ncols + B->ncols)) {
    m4ri_die("mzd_concat: C has wrong dimension!\n");
  }

  for (rci_t i = 0; i < A->nrows; ++i) {
    word *dst_truerow = mzd_row(C, i);
    word const *src_truerow = mzd_row_const(A, i);
    for (wi_t j = 0; j < A->width; ++j) { dst_truerow[j] = src_truerow[j]; }
  }

  for (rci_t i = 0; i < B->nrows; ++i) {
    for (rci_t j = 0; j < B->ncols; ++j) {
      mzd_write_bit(C, i, j + A->ncols, mzd_read_bit(B, i, j));
    }
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *mzd_stack(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  if (A->ncols != B->ncols) {
    m4ri_die("mzd_stack: A->ncols (%d) != B->ncols (%d)!\n", A->ncols, B->ncols);
  }

  if (C == NULL) {
    C = mzd_init(A->nrows + B->nrows, A->ncols);
  } else if (C->nrows != (A->nrows + B->nrows) || C->ncols != A->ncols) {
    m4ri_die("mzd_stack: C has wrong dimension!\n");
  }

  for (rci_t i = 0; i < A->nrows; ++i) {
    word const *src_truerow = mzd_row_const(A, i);
    word *dst_truerow = mzd_row(C, i);
    for (wi_t j = 0; j < A->width; ++j) { dst_truerow[j] = src_truerow[j]; }
  }

  for (rci_t i = 0; i < B->nrows; ++i) {
    word *dst_truerow = mzd_row(C, A->nrows + i);
    word const *src_truerow = mzd_row_const(B, i);
    for (wi_t j = 0; j < B->width; ++j) { dst_truerow[j] = src_truerow[j]; }
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *mzd_invert_naive(mzd_t *INV, mzd_t const *A, mzd_t const *I) {
  mzd_t *H;

  H = mzd_concat(NULL, A, I);

  rci_t x = mzd_echelonize_naive(H, TRUE);

  if (x == 0) {
    mzd_free(H);
    return NULL;
  }

  INV = mzd_submatrix(INV, H, 0, A->ncols, A->nrows, 2 * A->ncols);

  mzd_free(H);

  __M4RI_DD_MZD(INV);
  return INV;
}

mzd_t *mzd_add(mzd_t *ret, mzd_t const *left, mzd_t const *right) {
  if (left->nrows != right->nrows || left->ncols != right->ncols) {
    m4ri_die("mzd_add: rows and columns must match.\n");
  }
  if (ret == NULL) {
    ret = mzd_init(left->nrows, left->ncols);
  } else if (ret != left) {
    if (ret->nrows != left->nrows || ret->ncols != left->ncols) {
      m4ri_die("mzd_add: rows and columns of returned matrix must match.\n");
    }
  }
  return _mzd_add(ret, left, right);
}

mzd_t *_mzd_add(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  rci_t const nrows = MIN(MIN(A->nrows, B->nrows), C->nrows);

  if (C == B) {  // swap
    mzd_t const *tmp = A;
    A                = B;
    B                = tmp;
  }

  word const mask_end = C->high_bitmask;

  switch (A->width) {
  case 0: return C;
  case 1:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] ^= ((rowa[0] ^ rowb[0] ^ rowc[0]) & mask_end);
    }
    break;
  case 2:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] ^= ((rowa[1] ^ rowb[1] ^ rowc[1]) & mask_end);
    }
    break;
  case 3:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] = rowa[1] ^ rowb[1];
      rowc[2] ^= ((rowa[2] ^ rowb[2] ^ rowc[2]) & mask_end);
    }
    break;
  case 4:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] = rowa[1] ^ rowb[1];
      rowc[2] = rowa[2] ^ rowb[2];
      rowc[3] ^= ((rowa[3] ^ rowb[3] ^ rowc[3]) & mask_end);
    }
    break;
  case 5:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] = rowa[1] ^ rowb[1];
      rowc[2] = rowa[2] ^ rowb[2];
      rowc[3] = rowa[3] ^ rowb[3];
      rowc[4] ^= ((rowa[4] ^ rowb[4] ^ rowc[4]) & mask_end);
    }
    break;
  case 6:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] = rowa[1] ^ rowb[1];
      rowc[2] = rowa[2] ^ rowb[2];
      rowc[3] = rowa[3] ^ rowb[3];
      rowc[4] = rowa[4] ^ rowb[4];
      rowc[5] ^= ((rowa[5] ^ rowb[5] ^ rowc[5]) & mask_end);
    }
    break;
  case 7:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] = rowa[1] ^ rowb[1];
      rowc[2] = rowa[2] ^ rowb[2];
      rowc[3] = rowa[3] ^ rowb[3];
      rowc[4] = rowa[4] ^ rowb[4];
      rowc[5] = rowa[5] ^ rowb[5];
      rowc[6] ^= ((rowa[6] ^ rowb[6] ^ rowc[6]) & mask_end);
    }
    break;
  case 8:
    for (rci_t i = 0; i < nrows; ++i) {
      word const * rowa = mzd_row_const(A, i);
      word const * rowb = mzd_row_const(B, i);
      word * rowc = mzd_row(C, i);
      rowc[0] = rowa[0] ^ rowb[0];
      rowc[1] = rowa[1] ^ rowb[1];
      rowc[2] = rowa[2] ^ rowb[2];
      rowc[3] = rowa[3] ^ rowb[3];
      rowc[4] = rowa[4] ^ rowb[4];
      rowc[5] = rowa[5] ^ rowb[5];
      rowc[6] = rowa[6] ^ rowb[6];
      rowc[7] ^= ((rowa[7] ^ rowb[7] ^ rowc[7]) & mask_end);
    }
    break;

  default:
    for (rci_t i = 0; i < nrows; ++i) { mzd_combine_even(C, i, 0, A, i, 0, B, i, 0); }
  }

  __M4RI_DD_MZD(C);
  return C;
}

mzd_t *mzd_submatrix(mzd_t *S, mzd_t const *M, rci_t const startrow, rci_t const startcol,
                     rci_t const endrow, rci_t const endcol) {
  rci_t const nrows = endrow - startrow;
  rci_t const ncols = endcol - startcol;

  if (S == NULL) {
    S = mzd_init(nrows, ncols);
  } else if ((S->nrows < nrows) | (S->ncols < ncols)) {
    m4ri_die("mzd_submatrix: got S with dimension %d x %d but expected %d x %d\n", S->nrows,
             S->ncols, nrows, ncols);
  }

  if (startcol % m4ri_radix == 0) {

    wi_t const startword = startcol / m4ri_radix;
    /* we start at the beginning of a word */
    if (ncols / m4ri_radix != 0) {
      for (rci_t x = startrow, i = 0; i < nrows; ++i, ++x) {
        memcpy(mzd_row(S, i), mzd_row_const(M, x) + startword, sizeof(word) * (ncols / m4ri_radix));
      }
    }
    if (ncols % m4ri_radix) {
      word const mask_end = __M4RI_LEFT_BITMASK(ncols % m4ri_radix);
      for (rci_t x = startrow, i = 0; i < nrows; ++i, ++x) {
        /* process remaining bits */
        word temp                      = mzd_row_const(M, x)[startword + ncols / m4ri_radix] & mask_end;
        mzd_row(S, i)[ncols / m4ri_radix] = temp;
      }
    }
  } else {
    wi_t j;
    for (rci_t i = 0; i < nrows; i++) {
      word *srow = mzd_row(S, i);
      for (j = 0; j + m4ri_radix < ncols; j += m4ri_radix)
        srow[j / m4ri_radix] = mzd_read_bits(M, startrow + i, startcol + j, m4ri_radix);
      srow[j / m4ri_radix] &= ~S->high_bitmask;
      srow[j / m4ri_radix] |=
          mzd_read_bits(M, startrow + i, startcol + j, ncols - j) & S->high_bitmask;
    }
  }
  __M4RI_DD_MZD(S);
  return S;
}

int mzd_is_zero(mzd_t const *A) {
  word status   = 0;
  word mask_end = A->high_bitmask;
  for (rci_t i = 0; i < A->nrows; ++i) {
    word const *row = mzd_row_const(A, i);
    for (wi_t j = 0; j < A->width - 1; ++j) status |= row[j];
    status |= row[A->width - 1] & mask_end;
    if (status) return 0;
  }
  return !status;
}

void mzd_copy_row(mzd_t *B, rci_t i, mzd_t const *A, rci_t j) {
  assert(B->ncols >= A->ncols);
  wi_t const width = MIN(B->width, A->width) - 1;

  word const *a = mzd_row_const(A, j);
  word *b       = mzd_row(B, i);

  word const mask_end = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);

  if (width != 0) {
    for (wi_t k = 0; k < width; ++k) b[k] = a[k];
    b[width] = (b[width] & ~mask_end) | (a[width] & mask_end);

  } else {
    b[0] = (a[0] & mask_end) | (b[0] & ~mask_end);
  }

  __M4RI_DD_ROW(B, i);
}

int mzd_find_pivot(mzd_t const *A, rci_t start_row, rci_t start_col, rci_t *r, rci_t *c) {
  rci_t const nrows   = A->nrows;
  rci_t const ncols   = A->ncols;
  word data           = 0;
  rci_t row_candidate = 0;
  if (A->ncols - start_col < m4ri_radix) {
    for (rci_t j = start_col; j < A->ncols; j += m4ri_radix) {
      int const length = MIN(m4ri_radix, ncols - j);
      for (rci_t i = start_row; i < nrows; ++i) {
        word const curr_data = mzd_read_bits(A, i, j, length);
        if (m4ri_lesser_LSB(curr_data, data)) {
          row_candidate = i;
          data          = curr_data;
        }
      }
      if (data) {
        *r = row_candidate;
        for (int l = 0; l < length; ++l) {
          if (__M4RI_GET_BIT(data, l)) {
            *c = j + l;
            break;
          }
        }
        __M4RI_DD_RCI(*r);
        __M4RI_DD_RCI(*c);
        __M4RI_DD_INT(1);
        return 1;
      }
    }
  } else {
    /* we definitely have more than one word */
    /* handle first word */
    int const bit_offset   = (start_col % m4ri_radix);
    wi_t const word_offset = start_col / m4ri_radix;
    word const mask_begin  = __M4RI_RIGHT_BITMASK(m4ri_radix - bit_offset);
    for (rci_t i = start_row; i < nrows; ++i) {
      word const *row = mzd_row_const(A, i);
      word const curr_data = row[word_offset] & mask_begin;
      if (m4ri_lesser_LSB(curr_data, data)) {
        row_candidate = i;
        data          = curr_data;
        if (__M4RI_GET_BIT(data, bit_offset)) { break; }
      }
    }
    if (data) {
      *r = row_candidate;
      data >>= bit_offset;
      assert(data);
      for (int l = 0; l < (m4ri_radix - bit_offset); ++l) {
        if (__M4RI_GET_BIT(data, l)) {
          *c = start_col + l;
          break;
        }
      }
      __M4RI_DD_RCI(*r);
      __M4RI_DD_RCI(*c);
      __M4RI_DD_INT(1);
      return 1;
    }
    /* handle complete words */
    for (wi_t wi = word_offset + 1; wi < A->width - 1; ++wi) {
      for (rci_t i = start_row; i < nrows; ++i) {
        word const *row = mzd_row_const(A, i);
        word const curr_data = row[wi];
        if (m4ri_lesser_LSB(curr_data, data)) {
          row_candidate = i;
          data          = curr_data;
          if (__M4RI_GET_BIT(data, 0)) break;
        }
      }
      if (data) {
        *r = row_candidate;
        for (int l = 0; l < m4ri_radix; ++l) {
          if (__M4RI_GET_BIT(data, l)) {
            *c = wi * m4ri_radix + l;
            break;
          }
        }
        __M4RI_DD_RCI(*r);
        __M4RI_DD_RCI(*c);
        __M4RI_DD_INT(1);
        return 1;
      }
    }
    /* handle last word */
    int const end_offset = (A->ncols % m4ri_radix) ? (A->ncols % m4ri_radix) : m4ri_radix;
    word const mask_end  = __M4RI_LEFT_BITMASK(end_offset % m4ri_radix);
    wi_t wi              = A->width - 1;
    for (rci_t i = start_row; i < nrows; ++i) {
      word const *row = mzd_row_const(A, i);
      word const curr_data = row[wi] & mask_end;
      if (m4ri_lesser_LSB(curr_data, data)) {
        row_candidate = i;
        data          = curr_data;
        if (__M4RI_GET_BIT(data, 0)) break;
      }
    }
    if (data) {
      *r = row_candidate;
      for (int l = 0; l < end_offset; ++l) {
        if (__M4RI_GET_BIT(data, l)) {
          *c = wi * m4ri_radix + l;
          break;
        }
      }
      __M4RI_DD_RCI(*r);
      __M4RI_DD_RCI(*c);
      __M4RI_DD_INT(1);
      return 1;
    }
  }
  __M4RI_DD_RCI(*r);
  __M4RI_DD_RCI(*c);
  __M4RI_DD_INT(0);
  return 0;
}

#define MASK(c) (((uint64_t)(-1)) / (__M4RI_TWOPOW(__M4RI_TWOPOW(c)) + 1))
#define COUNT(x, c) ((x)&MASK(c)) + (((x) >> (__M4RI_TWOPOW(c))) & MASK(c))

static inline int m4ri_bitcount(word w) {
  uint64_t n = __M4RI_CONVERT_TO_UINT64_T(w);
  n          = COUNT(n, 0);
  n          = COUNT(n, 1);
  n          = COUNT(n, 2);
  n          = COUNT(n, 3);
  n          = COUNT(n, 4);
  n          = COUNT(n, 5);
  return (int)n;
}

double _mzd_density(mzd_t const *A, wi_t res, rci_t r, rci_t c) {
  size_t count = 0;
  size_t total = 0;

  if (A->width == 1) {
    for (rci_t i = r; i < A->nrows; ++i)
      for (rci_t j = c; j < A->ncols; ++j)
        if (mzd_read_bit(A, i, j)) ++count;
    return ((double)count) / (1.0 * A->ncols * A->nrows);
  }

  if (res == 0) res = A->width / 100;
  if (res < 1) res = 1;

  for (rci_t i = r; i < A->nrows; ++i) {
    word const *truerow = mzd_row_const(A, i);
    for (rci_t j = c; j < m4ri_radix; ++j)
      if (mzd_read_bit(A, i, j)) ++count;
    total += m4ri_radix;

    for (wi_t j = MAX(1, c / m4ri_radix); j < A->width - 1; j += res) {
      count += m4ri_bitcount(truerow[j]);
      total += m4ri_radix;
    }
    for (int j = 0; j < A->ncols % m4ri_radix; ++j)
      if (mzd_read_bit(A, i, m4ri_radix * (A->ncols / m4ri_radix) + j)) ++count;
    total += A->ncols % m4ri_radix;
  }

  return (double)count / total;
}

double mzd_density(mzd_t const *A, wi_t res) { return _mzd_density(A, res, 0, 0); }

rci_t mzd_first_zero_row(mzd_t const *A) {
  word const mask_end = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
  wi_t const end      = A->width - 1;
  for (rci_t i = A->nrows - 1; i >= 0; --i) {
    word const *row = mzd_row_const(A, i);
    word tmp = row[0];
    for (wi_t j = 1; j < end; ++j) tmp |= row[j];
    tmp |= row[end] & mask_end;
    if (tmp) {
      __M4RI_DD_INT(i + 1);
      return i + 1;
    }
  }
  __M4RI_DD_INT(0);
  return 0;
}

mzd_t *mzd_extract_u(mzd_t *U, mzd_t const *A) {
  rci_t k = MIN(A->nrows, A->ncols);
  if (U != NULL) { assert(U->nrows == k && U->ncols == k); }
  U = mzd_submatrix(U, A, 0, 0, k, k);
  for (rci_t i = 1; i < U->nrows; i++) {
    word *row = mzd_row(U, i);
    for (wi_t j = 0; j < i / m4ri_radix; j++) { row[j] = 0; }
    if (i % m4ri_radix) mzd_clear_bits(U, i, (i / m4ri_radix) * m4ri_radix, i % m4ri_radix);
  }
  return U;
}

mzd_t *mzd_extract_l(mzd_t *L, mzd_t const *A) {
  rci_t k = MIN(A->nrows, A->ncols);
  if (L != NULL) { assert(L->nrows == k && L->ncols == k); }
  L = mzd_submatrix(L, A, 0, 0, k, k);
  for (rci_t i = 0; i < L->nrows - 1; i++) {
    word *row = mzd_row(L, i);
    if (m4ri_radix - (i + 1) % m4ri_radix)
      mzd_clear_bits(L, i, i + 1, m4ri_radix - (i + 1) % m4ri_radix);
    for (wi_t j = (i / m4ri_radix + 1); j < L->width; j++) { row[j] = 0; }
  }
  return L;
}
