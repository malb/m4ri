/**
 * \file mzd.h
 * \brief Dense matrices over GF(2) represented as a bit field.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <martinralbrecht+m4ri@googlemail.com>
 * \author Carlo Wood <carlo@alinoe.com>
 */

#ifndef M4RI_MZD
#define M4RI_MZD

/*******************************************************************
 *
 *                M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
 *    Copyright (C) 2008-2013 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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
 *
 ********************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <m4ri/m4ri_config.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>

#if __M4RI_HAVE_SSE2
#include <emmintrin.h>
#endif

#include <m4ri/debug_dump.h>

/**
 * \brief Matrix multiplication block-ing dimension.
 *
 * Defines the number of rows of the matrix A that are
 * processed as one block during the execution of a multiplication
 * algorithm.
 */

#define __M4RI_MUL_BLOCKSIZE MIN(((int)sqrt((double)(4 * __M4RI_CPU_L3_CACHE))) / 2, 2048)


/**
 * \brief Dense matrices over GF(2).
 *
 * The most fundamental data type in this library.
 */

typedef struct mzd_t {

  rci_t nrows; /*!< Number of rows. */
  rci_t ncols; /*!< Number of columns. */
  wi_t width;  /*!< Number of words with valid bits: width = ceil(ncols / m4ri_radix) */

  /**
   * Offset in words between rows.
   */

  wi_t rowstride;

  /**
   * Booleans to speed up things.
   *
   * The bits have the following meaning:
   *
   * 1: Has non-zero excess.
   * 2: Is windowed, but has zero offset.
   * 3: Is windowed, but has zero excess.
   * 4: Is windowed, but owns the blocks allocations.
   * 5: Spans more than 1 block.
   */

  uint8_t flags;

  /* ensures sizeof(mzd_t) == 64 */
  uint8_t padding[63 - 2 * sizeof(rci_t) - 2 * sizeof(wi_t) - sizeof(word) - sizeof(void *)];

  word high_bitmask;   /*!< Mask for valid bits in the word with the highest index (width - 1). */
  word *data;
} mzd_t;

/* clang-format off */
/*
  The rows are stacked consecutively in memory.  The first word of the first row
  is given by M->data.

  Each row start M->rowstride words after the previous row (see the mzd_row function).

  Windows are "views" into an underlying matrix, with potentially a smaller number
  of rows and columns.  The window has the same rowstride as the underlying matrix,
  but they have a different data pointer.  Because a window may start anywhere,
  M->data cannot be expected to be aligned, at least not in windows.

  There may or may not be extra (full) words between rows (always when the matrix
  is a window, not necessarily when the matrix is not a window, depending on what
  we do about alignment).

  If (M->ncols % m4ri_radix != 0), then there are "excess bits" at the end of the row.
  M->high_bitmask tells which bits are meaningful.  M4RI policy is the following:
  * if the matrix is not a window, the excess bits MUST be zero. In this case, they
    MAY be overwritten (with fresh zeroes).
  * if the matrix is a window with non-zero excess, then the excess bits MUST be
    preserved (because they may be significant in the underlying matrix).


                  <------------------------ M->rowstride (in words) -------------------->

                  .-----------------------------------------------------------------------.
    M->data ----> |                                                        [excess 0 bits]|
                  |                                                        [excess 0 bits]|
                  |      .-------------------------------------.           [excess 0 bits]|
    W->data ------+----> |window                  [excess bits]|           [excess 0 bits]|
                  |      |                        [excess bits]|           [excess 0 bits]|
                  |      |                        [excess bits]|           [excess 0 bits]|
                  |      |                        [excess bits]|           [excess 0 bits]|
                  |      .-------------------------------------'           [excess 0 bits]|
                  |                                                        [excess 0 bits]|
                  `-----------------------------------------------------------------------'
*/
/* clang-format on */

/**
 * \brief flag when ncols%64 != 0
 */
static uint8_t const mzd_flag_nonzero_excess = 0x2;

/**
 * \brief flag for windowed matrix
 */

static uint8_t const mzd_flag_windowed = 0x4;

/**
 * \brief Test if a matrix is windowed.
 *
 * \param M Matrix
 *
 * \return a non-zero value if the matrix is windowed, otherwise return zero.
 */
static inline int mzd_is_windowed(mzd_t const *M) {
  return M->flags & mzd_flag_windowed;
}

/**
 * \brief Test if a matrix is windowed with non-zero excess.
 *
 * \param M Matrix
 *
 * \return a non-zero value if the matrix is windowed with non-zero excess, otherwise return zero.
 */
static inline int mzd_is_dangerous_window(mzd_t const *M) {
  uint8_t const danger = mzd_flag_windowed | mzd_flag_nonzero_excess;
  return (M->flags & danger) == danger;
}


/**
 * \brief Get pointer to first word of row.
 *
 * \param M Matrix
 * \param row The row index.
 *
 * \return pointer to first word of the row.
 */

static inline word *mzd_row(mzd_t *M, rci_t row) {
  return M->data + M->rowstride * row;
}

static inline word const * mzd_row_const(mzd_t const *M, rci_t row) {
  return mzd_row((mzd_t *)M, row);
}

/**
 * \brief Create a new matrix of dimension r x c.
 *
 * Use mzd_free to kill it.
 *
 * \param r Number of rows
 * \param c Number of columns
 *
 */

mzd_t *mzd_init(rci_t const r, rci_t const c);

/**
 * \brief Free a matrix created with mzd_init.
 *
 * \param A Matrix
 */

void mzd_free(mzd_t *A);

/**
 * \brief Create a window/view into the matrix M.
 *
 * A matrix window for M is a meta structure on the matrix M. It is
 * setup to point into the matrix so M \em must \em not be freed while the
 * matrix window is used.
 *
 * This function puts the restriction on the provided parameters that
 * all parameters must be within range for M which is not enforced
 * currently .
 *
 * Use mzd_free_window to free the window.
 *
 * \param M Matrix
 * \param lowr Starting row (inclusive)
 * \param lowc Starting column (inclusive, must be multiple of m4ri_radix)
 * \param highr End row (exclusive)
 * \param highc End column (exclusive)
 *
 */

mzd_t *mzd_init_window(mzd_t *M, rci_t const lowr, rci_t const lowc, rci_t const highr,
                       rci_t const highc);

/**
 * \brief Create a const window/view into a const matrix M.
 *
 * See mzd_init_window, but for constant M.
 */

static inline mzd_t const *mzd_init_window_const(mzd_t const *M, rci_t const lowr, rci_t const lowc,
                                                 rci_t const highr, rci_t const highc) {
  return mzd_init_window((mzd_t *)M, lowr, lowc, highr, highc);
}

/**
 * \brief Free a matrix window created with mzd_init_window.
 *
 * \param A Matrix
 */

#define mzd_free_window mzd_free

/**
 * \brief Swap the two rows rowa and rowb starting at startblock.
 *
 * \param M Matrix with a zero offset.
 * \param rowa Row index.
 * \param rowb Row index.
 * \param startblock Start swapping only in this block.
 */

static inline void _mzd_row_swap(mzd_t *M, rci_t const rowa, rci_t const rowb,
                                 wi_t const startblock) {
  if ((rowa == rowb) || (startblock >= M->width)) { return; }

  wi_t width = M->width - startblock - 1;
  word *a    = mzd_row(M, rowa) + startblock;
  word *b    = mzd_row(M, rowb) + startblock;
  word tmp;
  word const mask_end = M->high_bitmask;

  for (wi_t i = 0; i < width; ++i) {
    tmp  = a[i];
    a[i] = b[i];
    b[i] = tmp;
  }
  tmp = (a[width] ^ b[width]) & mask_end;
  a[width] ^= tmp;
  b[width] ^= tmp;

  __M4RI_DD_ROW(M, rowa);
  __M4RI_DD_ROW(M, rowb);
}

/**
 * \brief Swap the two rows rowa and rowb.
 *
 * \param M Matrix
 * \param rowa Row index.
 * \param rowb Row index.
 */

static inline void mzd_row_swap(mzd_t *M, rci_t const rowa, rci_t const rowb) {
  _mzd_row_swap(M, rowa, rowb, 0);
}

/**
 * \brief copy row j from A to row i from B.
 *
 * The offsets of A and B must match and the number of columns of A
 * must be less than or equal to the number of columns of B.
 *
 * \param B Target matrix.
 * \param i Target row index.
 * \param A Source matrix.
 * \param j Source row index.
 */

void mzd_copy_row(mzd_t *B, rci_t i, mzd_t const *A, rci_t j);


/**
 * \brief Swap the two columns cola and colb but only between start_row and stop_row.
 *
 * \param M Matrix.
 * \param cola Column index.
 * \param colb Column index.
 * \param start_row Row index.
 * \param stop_row Row index (exclusive).
 */

static inline void mzd_col_swap_in_rows(mzd_t *M, rci_t const cola, rci_t const colb,
                                        rci_t const start_row, rci_t const stop_row) {
  if (cola == colb) { return; }

  rci_t const _cola = cola;
  rci_t const _colb = colb;

  wi_t const a_word = _cola / m4ri_radix;
  wi_t const b_word = _colb / m4ri_radix;

  int const a_bit = _cola % m4ri_radix;
  int const b_bit = _colb % m4ri_radix;

  word *RESTRICT ptr  = mzd_row(M, start_row);
  int max_bit         = MAX(a_bit, b_bit);
  int count_remaining = stop_row - start_row;
  int min_bit         = a_bit + b_bit - max_bit;
  int offset          = max_bit - min_bit;
  word mask           = m4ri_one << min_bit;
  int count           = count_remaining;

  // Apparently we're calling with start_row == stop_row sometimes (seems a bug to me).
  if (count <= 0) { return; }

  if (a_word == b_word) {
    while (1) {
      count_remaining -= count;
      assert(count_remaining == 0);
      ptr += a_word;
      int fast_count = count / 4;
      int rest_count = count - 4 * fast_count;
      word xor_v[4];
      wi_t const rowstride = M->rowstride;
      while (fast_count--) {
        xor_v[0] = ptr[0];
        xor_v[1] = ptr[rowstride];
        xor_v[2] = ptr[2 * rowstride];
        xor_v[3] = ptr[3 * rowstride];
        xor_v[0] ^= xor_v[0] >> offset;
        xor_v[1] ^= xor_v[1] >> offset;
        xor_v[2] ^= xor_v[2] >> offset;
        xor_v[3] ^= xor_v[3] >> offset;
        xor_v[0] &= mask;
        xor_v[1] &= mask;
        xor_v[2] &= mask;
        xor_v[3] &= mask;
        xor_v[0] |= xor_v[0] << offset;
        xor_v[1] |= xor_v[1] << offset;
        xor_v[2] |= xor_v[2] << offset;
        xor_v[3] |= xor_v[3] << offset;
        ptr[0] ^= xor_v[0];
        ptr[rowstride] ^= xor_v[1];
        ptr[2 * rowstride] ^= xor_v[2];
        ptr[3 * rowstride] ^= xor_v[3];
        ptr += 4 * rowstride;
      }
      while (rest_count--) {
        word xor_v = *ptr;
        xor_v ^= xor_v >> offset;
        xor_v &= mask;
        *ptr ^= xor_v | (xor_v << offset);
        ptr += rowstride;
      }
      break;
    }
  } else {
    word *RESTRICT min_ptr;
    wi_t max_offset;
    if (min_bit == a_bit) {
      min_ptr    = ptr + a_word;
      max_offset = b_word - a_word;
    } else {
      min_ptr    = ptr + b_word;
      max_offset = a_word - b_word;
    }
    while (1) {
      count_remaining -= count;
      assert(count_remaining == 0);
      wi_t const rowstride = M->rowstride;
      while (count--) {
        word xor_v = (min_ptr[0] ^ (min_ptr[max_offset] >> offset)) & mask;
        min_ptr[0] ^= xor_v;
        min_ptr[max_offset] ^= xor_v << offset;
        min_ptr += rowstride;
      }
      break;
    }
  }

  __M4RI_DD_MZD(M);
}

/**
 * \brief Swap the two columns cola and colb.
 *
 * \param M Matrix.
 * \param cola Column index.
 * \param colb Column index.
 */
static inline void mzd_col_swap(mzd_t *M, rci_t const cola, rci_t const colb) {
  mzd_col_swap_in_rows(M, cola, colb, 0, M->nrows);
}


/**
 * \brief Read the bit at position M[row,col].
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 *
 * \note No bounds checks whatsoever are performed.
 *
 */

static inline BIT mzd_read_bit(mzd_t const *M, rci_t const row, rci_t const col) {
  word const * truerow = mzd_row_const(M, row);
  return __M4RI_GET_BIT(truerow[col / m4ri_radix], col % m4ri_radix);
}

/**
 * \brief Write the bit value to position M[row,col]
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 * \param value Either 0 or 1
 *
 * \note No bounds checks whatsoever are performed.
 *
 */

static inline void mzd_write_bit(mzd_t *M, rci_t const row, rci_t const col, BIT const value) {
    word * truerow = mzd_row(M, row);
  __M4RI_WRITE_BIT(truerow[col / m4ri_radix], col % m4ri_radix, value);
}

/**
 * \brief XOR n bits from values to M starting a position (x,y).
 *
 * \param M Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param n Number of bits (<= m4ri_radix);
 * \param values Word with values;
 */

static inline void mzd_xor_bits(mzd_t *M, rci_t const x, rci_t const y, int const n,
                                word values) {
  int const spot   = y % m4ri_radix;
  wi_t const block = y / m4ri_radix;
  word *row = mzd_row(M, x);
  row[block] ^= values << spot;
  int const space = m4ri_radix - spot;
  if (n > space) { row[block + 1] ^= values >> space; }
}

/**
 * \brief AND n bits from values to M starting a position (x,y).
 *
 * \param M Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param n Number of bits (<= m4ri_radix);
 * \param values Word with values;
 */

static inline void mzd_and_bits(mzd_t *M, rci_t const x, rci_t const y, int const n,
                                word values) {
  /* This is the best way, since this will drop out once we inverse the bits in values: */
  values >>= (m4ri_radix - n); /* Move the bits to the lowest columns */

  int const spot   = y % m4ri_radix;
  wi_t const block = y / m4ri_radix;
  word *row = mzd_row(M, x);
  row[block] &= values << spot;
  int const space = m4ri_radix - spot;
  if (n > space) { row[block + 1] &= values >> space; }
}

/**
 * \brief Clear n bits in M starting a position (x,y).
 *
 * \param M Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param n Number of bits (0 < n <= m4ri_radix);
 */

static inline void mzd_clear_bits(mzd_t *M, rci_t const x, rci_t const y, int const n) {
  assert(n > 0 && n <= m4ri_radix);
  word values      = m4ri_ffff >> (m4ri_radix - n);
  int const spot   = y % m4ri_radix;
  wi_t const block = y / m4ri_radix;
  word *row = mzd_row(M, x);
  row[block] &= ~(values << spot);
  int const space = m4ri_radix - spot;
  if (n > space) { row[block + 1] &= ~(values >> space); }
}

/**
 * \brief Add the rows sourcerow and destrow and stores the total in the row
 * destrow, but only begins at the column coloffset.
 *
 * \param M Matrix
 * \param dstrow Index of target row
 * \param srcrow Index of source row
 * \param coloffset Start column (0 <= coloffset < M->ncols)
 *
 * \warning This function expects that there is at least one word worth of work.
 */

static inline void mzd_row_add_offset(mzd_t *M, rci_t dstrow, rci_t srcrow, rci_t coloffset) {
  assert(dstrow < M->nrows && srcrow < M->nrows && coloffset < M->ncols);
  wi_t const startblock = coloffset / m4ri_radix;
  wi_t wide             = M->width - startblock;
  word *src             = mzd_row(M, srcrow) + startblock;
  word *dst             = mzd_row(M, dstrow) + startblock;
  word const mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - coloffset % m4ri_radix);
  word const mask_end   = M->high_bitmask;

  *dst++ ^= *src++ & mask_begin;
  --wide;

#if __M4RI_HAVE_SSE2
  int not_aligned = __M4RI_ALIGNMENT(src, 16) != 0; /* 0: Aligned, 1: Not aligned */
  if (wide > not_aligned + 1)                       /* Speed up for small matrices */
  {
    if (not_aligned) {
      *dst++ ^= *src++;
      --wide;
    }
    /* Now wide > 1 */
    __m128i *__src     = (__m128i *)src;
    __m128i *__dst     = (__m128i *)dst;
    __m128i *const eof = (__m128i *)((uintptr_t)(src + wide) & ~(uintptr_t)0xF);
    do {
      __m128i xmm1 = _mm_xor_si128(*__dst, *__src);
      *__dst++     = xmm1;
    } while (++__src < eof);
    src  = (word *)__src;
    dst  = (word *)__dst;
    wide = ((sizeof(word) * wide) % 16) / sizeof(word);
  }
#endif
  wi_t i = -1;
  while (++i < wide) { dst[i] ^= src[i]; }
  /*
   * Revert possibly non-zero excess bits.
   * Note that i == wide here, and wide can be 0.
   * But really, src[wide - 1] is M->rows[srcrow][M->width - 1] ;)
   * We use i - 1 here to let the compiler know these are the same addresses
   * that we last accessed, in the previous loop.
   */
  dst[i - 1] ^= src[i - 1] & ~mask_end;

  __M4RI_DD_ROW(M, dstrow);
}

/**
 * \brief Add the rows sourcerow and destrow and stores the total in
 * the row destrow.
 *
 * \param M Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 *
 * \note this can be done much faster with mzd_combine.
 */

void mzd_row_add(mzd_t *M, rci_t const sourcerow, rci_t const destrow);

/**
 * \brief Transpose a matrix.
 *
 * This function uses the fact that:
\verbatim
   [ A B ]T    [AT CT]
   [ C D ]  =  [BT DT]
 \endverbatim
 * and thus rearranges the blocks recursively.
 *
 * \param DST Preallocated return matrix, may be NULL for automatic creation.
 * \param A Matrix
 */

mzd_t *mzd_transpose(mzd_t *DST, mzd_t const *A);

/**
 * \brief Naive cubic matrix multiplication.
 *
 * That is, compute C such that C == AB.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \note Normally, if you will multiply several times by b, it is
 * smarter to calculate bT yourself, and keep it, and then use the
 * function called _mzd_mul_naive
 *
 */
mzd_t *mzd_mul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B);

/**
 * \brief Naive cubic matrix multiplication and addition
 *
 * That is, compute C such that C == C + AB.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \note Normally, if you will multiply several times by b, it is
 * smarter to calculate bT yourself, and keep it, and then use the
 * function called _mzd_mul_naive
 */

mzd_t *mzd_addmul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B);

/**
 * \brief Naive cubic matrix multiplication with the pre-transposed B.
 *
 * That is, compute C such that C == AB^t.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Pre-transposed input matrix B.
 * \param clear Whether to clear C before accumulating AB
 */

mzd_t *_mzd_mul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B, int const clear);

/**
 * \brief Matrix multiplication optimized for v*A where v is a vector.
 *
 * \param C Preallocated product matrix.
 * \param v Input matrix v.
 * \param A Input matrix A.
 * \param clear If set clear C first, otherwise add result to C.
 *
 */
mzd_t *_mzd_mul_va(mzd_t *C, mzd_t const *v, mzd_t const *A, int const clear);

/**
 * \brief Fill matrix M with uniformly distributed bits.
 *
 * \param M Matrix
 */

void mzd_randomize(mzd_t *M);

/**
 * \brief Random callback that produces uniformly distributed random
 * words on every call.
 *
 * \param data callback data
 *
 * \return uniformly distributed random word
 */
typedef word (*m4ri_random_callback)(void *data);

/**
 * \brief Fill matrix M with uniformly distributed bits.
 *
 * \param M Matrix
 * \param rc callback
 * \param data callback data passed to every call to rc
 */
void mzd_randomize_custom(mzd_t *M, m4ri_random_callback rc, void *data);

/**
 * \brief Set the matrix M to the value equivalent to the integer
 * value provided.
 *
 * Specifically, this function does nothing if value%2 == 0 and
 * returns the identity matrix if value%2 == 1.
 *
 * If the matrix is not square then the largest possible square
 * submatrix is set to the identity matrix.
 *
 * \param M Matrix
 * \param value Either 0 or 1
 */

void mzd_set_ui(mzd_t *M, unsigned int const value);

/**
 * \brief Gaussian elimination.
 *
 * This will do Gaussian elimination on the matrix m but will start
 * not at column 0 necc but at column startcol. If full=FALSE, then it
 * will do triangular style elimination, and if full=TRUE, it will do
 * Gauss-Jordan style, or full elimination.
 *
 * \param M Matrix
 * \param startcol First column to consider for reduction.
 * \param full Gauss-Jordan style or upper triangular form only.
 */

rci_t mzd_gauss_delayed(mzd_t *M, rci_t const startcol, int const full);

/**
 * \brief Gaussian elimination.
 *
 * This will do Gaussian elimination on the matrix m.  If full=FALSE,
 *  then it will do triangular style elimination, and if full=TRUE,
 *  it will do Gauss-Jordan style, or full elimination.
 *
 * \param M Matrix
 * \param full Gauss-Jordan style or upper triangular form only.
 *
 * \sa mzd_echelonize_m4ri(), mzd_echelonize_pluq()
 */

rci_t mzd_echelonize_naive(mzd_t *M, int const full);

/**
 * \brief Return TRUE if A == B.
 *
 * \param A Matrix
 * \param B Matrix
 */

int mzd_equal(mzd_t const *A, mzd_t const *B);

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined mathematically and
 * relatively arbitrary since elements of GF(2) don't have an
 * ordering.
 */

int mzd_cmp(mzd_t const *A, mzd_t const *B);

/**
 * \brief Copy matrix  A to DST.
 *
 * \param DST May be NULL for automatic creation.
 * \param A Source matrix.
 */

mzd_t *mzd_copy(mzd_t *DST, mzd_t const *A);

/**
 * \brief Concatenate B to A and write the result to C.
 *
 * That is,
 *
 \verbatim
 [ A ], [ B ] -> [ A  B ] = C
 \endverbatim
 *
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 *
 * \note This is sometimes called augment.
 */

mzd_t *mzd_concat(mzd_t *C, mzd_t const *A, mzd_t const *B);

/**
 * \brief Stack A on top of B and write the result to C.
 *
 * That is,
 *
 \verbatim
 [ A ], [ B ] -> [ A ] = C
                 [ B ]
 \endverbatim
 *
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 */

mzd_t *mzd_stack(mzd_t *C, mzd_t const *A, mzd_t const *B);

/**
 * \brief Copy a submatrix.
 *
 * Note that the upper bounds are not included.
 *
 * \param S Preallocated space for submatrix, may be NULL for automatic creation.
 * \param M Matrix
 * \param lowr start rows
 * \param lowc start column
 * \param highr stop row (this row is \em not included)
 * \param highc stop column (this column is \em not included)
 */
mzd_t *mzd_submatrix(mzd_t *S, mzd_t const *M, rci_t const lowr, rci_t const lowc,
                     rci_t const highr, rci_t const highc);

/**
 * \brief Invert the matrix target using Gaussian elimination.
 *
 * To avoid recomputing the identity matrix over and over again, I may
 * be passed in as identity parameter.
 *
 * \param INV Preallocated space for inversion matrix, may be NULL for automatic creation.
 * \param A Matrix to be reduced.
 * \param I Identity matrix.
 */

mzd_t *mzd_invert_naive(mzd_t *INV, mzd_t const *A, mzd_t const *I);

/**
 * \brief Set C = A+B.
 *
 * C is also returned. If C is NULL then a new matrix is created which
 * must be freed by mzd_free.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

mzd_t *mzd_add(mzd_t *C, mzd_t const *A, mzd_t const *B);

/**
 * \brief Same as mzd_add but without any checks on the input.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

mzd_t *_mzd_add(mzd_t *C, mzd_t const *A, mzd_t const *B);

/**
 * \brief Same as mzd_add.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

#define mzd_sub mzd_add

/**
 * \brief Same as mzd_sub but without any checks on the input.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

#define _mzd_sub _mzd_add

/**
 * Get n bits starting a position (x,y) from the matrix M.
 *
 * \param M Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param n Number of bits (<= m4ri_radix);
 */

static inline word mzd_read_bits(mzd_t const *M, rci_t const x, rci_t const y, int const n) {
  int const spot   = y % m4ri_radix;
  wi_t const block = y / m4ri_radix;
  int const spill  = spot + n - m4ri_radix;
  word const *row = mzd_row_const(M, x);
  word temp        = (spill <= 0)
                  ? row[block] << -spill
                  : (row[block + 1] << (m4ri_radix - spill)) | (row[block] >> spill);
  return temp >> (m4ri_radix - n);
}

/**
 * \brief a_row[a_startblock:] += b_row[b_startblock:] for offset 0
 *
 * Adds a_row of A, starting with a_startblock to the end, to
 * b_row of B, starting with b_startblock to the end. This gets stored
 * in A, in a_row, starting with a_startblock.
 *
 * \param A destination matrix
 * \param a_row destination row for matrix C
 * \param a_startblock starting block to work on in matrix C
 * \param B source matrix
 * \param b_row source row for matrix B
 * \param b_startblock starting block to work on in matrix B
 *
 */

static inline void mzd_combine_even_in_place(mzd_t *A, rci_t const a_row, wi_t const a_startblock,
                                             mzd_t const *B, rci_t const b_row,
                                             wi_t const b_startblock) {

  wi_t wide = A->width - a_startblock - 1;

  word *a = mzd_row(A, a_row) + a_startblock;
  word const *b = mzd_row_const(B, b_row) + b_startblock;

#if __M4RI_HAVE_SSE2
  if (wide > 2) {
    /** check alignments **/
    if (__M4RI_ALIGNMENT(a, 16)) {
      *a++ ^= *b++;
      wide--;
    }

    if (__M4RI_ALIGNMENT(a, 16) == 0 && __M4RI_ALIGNMENT(b, 16) == 0) {
      __m128i *a128      = (__m128i *)a;
      __m128i *b128      = (__m128i *)b;
      const __m128i *eof = (__m128i *)((uintptr_t)(a + wide) & ~(uintptr_t)0xF);

      do {
        *a128 = _mm_xor_si128(*a128, *b128);
        ++b128;
        ++a128;
      } while (a128 < eof);

      a    = (word *)a128;
      b    = (word *)b128;
      wide = ((sizeof(word) * wide) % 16) / sizeof(word);
    }
  }
#endif  // __M4RI_HAVE_SSE2

  if (wide > 0) {
    wi_t n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *(a++) ^= *(b++);
      case 7: *(a++) ^= *(b++);
      case 6: *(a++) ^= *(b++);
      case 5: *(a++) ^= *(b++);
      case 4: *(a++) ^= *(b++);
      case 3: *(a++) ^= *(b++);
      case 2: *(a++) ^= *(b++);
      case 1: *(a++) ^= *(b++);
      } while (--n > 0);
    }
  }

  *a ^= *b & A->high_bitmask;

  __M4RI_DD_MZD(A);
}

/**
 * \brief c_row[c_startblock:] = a_row[a_startblock:] + b_row[b_startblock:] for offset 0
 *
 * Adds a_row of A, starting with a_startblock to the end, to
 * b_row of B, starting with b_startblock to the end. This gets stored
 * in C, in c_row, starting with c_startblock.
 *
 * \param C destination matrix
 * \param c_row destination row for matrix C
 * \param c_startblock starting block to work on in matrix C
 * \param A source matrix
 * \param a_row source row for matrix A
 * \param a_startblock starting block to work on in matrix A
 * \param B source matrix
 * \param b_row source row for matrix B
 * \param b_startblock starting block to work on in matrix B
 *
 */

static inline void mzd_combine_even(mzd_t *C, rci_t const c_row, wi_t const c_startblock,
                                    mzd_t const *A, rci_t const a_row, wi_t const a_startblock,
                                    mzd_t const *B, rci_t const b_row, wi_t const b_startblock) {

  wi_t wide = A->width - a_startblock - 1;
  word const *a = mzd_row_const(A, a_row) + a_startblock;
  word const *b = mzd_row_const(B, b_row) + b_startblock;
  word *c   = mzd_row(C, c_row) + c_startblock;

#if __M4RI_HAVE_SSE2
  if (wide > 2) {
    /** check alignments **/
    if (__M4RI_ALIGNMENT(a, 16)) {
      *c++ = *b++ ^ *a++;
      wide--;
    }

    if ((__M4RI_ALIGNMENT(b, 16) | __M4RI_ALIGNMENT(c, 16)) == 0) {
      __m128i *a128      = (__m128i *)a;
      __m128i *b128      = (__m128i *)b;
      __m128i *c128      = (__m128i *)c;
      const __m128i *eof = (__m128i *)((uintptr_t)(a + wide) & ~(uintptr_t)0xF);

      do {
        *c128 = _mm_xor_si128(*a128, *b128);
        ++c128;
        ++b128;
        ++a128;
      } while (a128 < eof);

      a    = (word *)a128;
      b    = (word *)b128;
      c    = (word *)c128;
      wide = ((sizeof(word) * wide) % 16) / sizeof(word);
    }
  }
#endif  // __M4RI_HAVE_SSE2

  if (wide > 0) {
    wi_t n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *(c++) = *(a++) ^ *(b++);
      case 7: *(c++) = *(a++) ^ *(b++);
      case 6: *(c++) = *(a++) ^ *(b++);
      case 5: *(c++) = *(a++) ^ *(b++);
      case 4: *(c++) = *(a++) ^ *(b++);
      case 3: *(c++) = *(a++) ^ *(b++);
      case 2: *(c++) = *(a++) ^ *(b++);
      case 1: *(c++) = *(a++) ^ *(b++);
      } while (--n > 0);
    }
  }
  *c ^= ((*a ^ *b ^ *c) & C->high_bitmask);

  __M4RI_DD_MZD(C);
}

/**
 * \brief row3[col3:] = row1[col1:] + row2[col2:]
 *
 * Adds row1 of SC1, starting with startblock1 to the end, to
 * row2 of SC2, starting with startblock2 to the end. This gets stored
 * in DST, in row3, starting with startblock3.
 *
 * \param C destination matrix
 * \param c_row destination row for matrix dst
 * \param c_startblock starting block to work on in matrix dst
 * \param A source matrix
 * \param a_row source row for matrix sc1
 * \param a_startblock starting block to work on in matrix sc1
 * \param B source matrix
 * \param b_row source row for matrix sc2
 * \param b_startblock starting block to work on in matrix sc2
 *
 */
static inline void mzd_combine(mzd_t *C, rci_t const c_row, wi_t const c_startblock, mzd_t const *A,
                               rci_t const a_row, wi_t const a_startblock, mzd_t const *B,
                               rci_t const b_row, wi_t const b_startblock) {

  if ((C == A) & (a_row == c_row) & (a_startblock == c_startblock)) {
    mzd_combine_even_in_place(C, c_row, c_startblock, B, b_row, b_startblock);
  } else {
    mzd_combine_even(C, c_row, c_startblock, A, a_row, a_startblock, B, b_row, b_startblock);
  }
  return;
}

/**
 * \brief Get n bits starting a position (x,y) from the matrix M.
 *
 * This function is in principle the same as mzd_read_bits,
 * but it explicitely returns an 'int' and is used as
 * index into an array (Gray code).
 */

static inline int mzd_read_bits_int(mzd_t const *M, rci_t const x, rci_t const y, int const n) {
  return __M4RI_CONVERT_TO_INT(mzd_read_bits(M, x, y, n));
}

/**
 * \brief Zero test for matrix.
 *
 * \param A Input matrix.
 *
 */
int mzd_is_zero(mzd_t const *A);

/**
 * \brief Clear the given row, but only begins at the column coloffset.
 *
 * \param M Matrix
 * \param row Index of row
 * \param coloffset Column offset
 */

void mzd_row_clear_offset(mzd_t *M, rci_t const row, rci_t const coloffset);

/**
 * \brief Find the next nonzero entry in M starting at start_row and start_col.
 *
 * This function walks down rows in the inner loop and columns in the
 * outer loop. If a nonzero entry is found this function returns 1 and
 * zero otherwise.
 *
 * If and only if a nonzero entry is found r and c are updated.
 *
 * \param M Matrix
 * \param start_row Index of row where to start search
 * \param start_col Index of column where to start search
 * \param r Row index updated if pivot is found
 * \param c Column index updated if pivot is found
 */

int mzd_find_pivot(mzd_t const *M, rci_t start_row, rci_t start_col, rci_t *r, rci_t *c);

/**
 * \brief Return the number of nonzero entries divided by nrows *
 * ncols
 *
 * If res = 0 then 100 samples per row are made, if res > 0 the
 * function takes res sized steps within each row (res = 1 uses every
 * word).
 *
 * \param A Matrix
 * \param res Resolution of sampling (in words)
 */

double mzd_density(mzd_t const *A, wi_t res);

/**
 * \brief Return the number of nonzero entries divided by nrows *
 * ncols considering only the submatrix starting at (r,c).
 *
 * If res = 0 then 100 samples per row are made, if res > 0 the
 * function takes res sized steps within each row (res = 1 uses every
 * word).
 *
 * \param A Matrix
 * \param res Resolution of sampling (in words)
 * \param r Row to start counting
 * \param c Column to start counting
 */

double _mzd_density(mzd_t const *A, wi_t res, rci_t r, rci_t c);

/**
 * \brief Return the first row with all zero entries.
 *
 * If no such row can be found returns nrows.
 *
 * \param A Matrix
 */

rci_t mzd_first_zero_row(mzd_t const *A);

/**
 * \brief Return hash value for matrix.
 *
 * \param A Matrix
 */

static inline word mzd_hash(mzd_t const *A) {
  word hash = 0;
  for (rci_t r = 0; r < A->nrows; ++r) {
    hash ^= rotate_word(calculate_hash(mzd_row_const(A, r), A->width), r % m4ri_radix);
  }
  return hash;
}

/**
 * Return upper triangular submatrix of A
 *
 * \param U Output matrix, if NULL a new matrix will be returned
 * \param A Source matrix
 *
 * \return U
 */

mzd_t *mzd_extract_u(mzd_t *U, mzd_t const *A);

/**
 * Return lower triangular submatrix of A
 *
 * \param L Output matrix, if NULL a new matrix will be returned
 * \param A Source matrix
 *
 * \return L
 */

mzd_t *mzd_extract_l(mzd_t *L, mzd_t const *A);

#endif  // M4RI_MZD
