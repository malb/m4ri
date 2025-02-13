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

#include "mzd.h"

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
