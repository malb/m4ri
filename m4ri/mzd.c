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
#include "transpose.h"

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
