/******************************************************************************
*
*            M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*    Copyright (C) 2009,2010 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <stdlib.h>
#include <string.h>
#include "packedmatrix.h"
#include "parity.h"

#define SAFECHAR (m4ri_radix + m4ri_radix / 4 + 1)

/*
 * Return r such that x elements fit into r blocks of length y.
 */
#define DIV_CEIL(x,y) (((x) % (y)) ? (x) / (y) + 1 : (x) / (y))

/*
 * Maximum number of bytes allocated in one malloc() call.
 * This value must fit in an int, even though it's type is size_t.
 */
#define __M4RI_MM_MAX_MALLOC (((size_t)1)<<30)

/*
 * Enable memory block cache (default: enabled)
 */
#define __M4RI_ENABLE_MMC

/*
 * Number of blocks that are cached.
 */
#define __M4RI_MMC_NBLOCKS 16

/*
 * Maximal size of blocks stored in cache.
 */
#define __M4RI_MMC_THRESHOLD CPU_L2_CACHE

/*
 * The mmc memory management functions check a cache for re-usable
 * unused memory before asking the system for it.
 */
typedef struct _mm_block {
  /**
   * Size in bytes of the data.
   */
  size_t size;

  /**
   * Pointer to buffer of data.
   */
  void *data;

} mmb_t;

#ifdef __M4RI_ENABLE_MMC
/*
 * The actual memory block cache.
 */
static mmb_t m4ri_mmc_cache[__M4RI_MMC_NBLOCKS];
#endif // __M4RI_ENABLE_MMC

/*
 * Allocate size bytes.
 *
 * size: Number of bytes.
 *
 * Returns pointer to allocated memory block.
 */
static void *m4ri_mmc_malloc(size_t size) {

#ifdef __M4RI_ENABLE_MMC
  void *ret = NULL;
#endif

#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif

#ifdef __M4RI_ENABLE_MMC
  mmb_t *mm = m4ri_mmc_cache;
  if (size <= __M4RI_MMC_THRESHOLD) {
    for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
      if(mm[i].size == size) {
        ret = mm[i].data;
        mm[i].data = NULL;
        mm[i].size = 0;
        break;
      }
    }
  }
#endif // __M4RI_ENABLE_MMC

#ifdef HAVE_OPENMP
 }
#endif

#ifdef __M4RI_ENABLE_MMC
 if (ret)
   return ret;
 else
   return m4ri_mm_malloc(size);
#else 
 return m4ri_mm_malloc(size);
#endif
}

/*
 * Allocate an array of count times size zeroed bytes.
 *
 * count: Number of elements.
 * size: Number of bytes per element.
 *
 * Returns pointer to allocated memory block.
 */
static inline void *m4ri_mmc_calloc(size_t count, size_t size) {
  void *ret = m4ri_mmc_malloc(count * size);
  memset((char*)ret, 0, count * size);
  return ret;
}

/*
 * Free the data pointed to by condemned of the given size.
 *
 * condemned: Pointer to memory.
 * size: Number of bytes.
 */
static void m4ri_mmc_free(void *condemned, size_t size) {
#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
#ifdef __M4RI_ENABLE_MMC
  static int j = 0;
  mmb_t *mm = m4ri_mmc_cache;
  if (size < __M4RI_MMC_THRESHOLD) {
    for(int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
      if(mm[i].size == 0) {
        mm[i].size = size;
        mm[i].data = condemned;
        return;
      }
    }
    m4ri_mm_free(mm[j].data);
    mm[j].size = size;
    mm[j].data = condemned;
    j = (j+1) % __M4RI_MMC_NBLOCKS;
  } else {
    m4ri_mm_free(condemned);
  }
#else
  m4ri_mm_free(condemned);
#endif // __M4RI_ENABLE_MMC
#ifdef HAVE_OPENMP
 }
#endif
}

/*
 * Cleans up memory block cache.
 *
 * This function is called automatically when the shared library is loaded.
 */
void m4ri_mmc_cleanup(void) {
#ifdef HAVE_OPENMP
#pragma omp critical
{
#endif
#ifdef __M4RI_ENABLE_MMC
  mmb_t *mm = m4ri_mmc_cache;
  for(int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (mm[i].size)
      m4ri_mm_free(mm[i].data);
    mm[i].size = 0;
  }
#endif // __M4RI_ENABLE_MMC
#ifdef HAVE_OPENMP
 }
#endif
}

mzd_t *mzd_init(rci_t r, rci_t c) {

  mzd_t *A = (mzd_t *)m4ri_mmc_malloc(sizeof(mzd_t));

  A->width = DIV_CEIL(c, m4ri_radix);

#ifdef HAVE_SSE2
  int incw = 0;
  /* make sure each row is 16-byte aligned */
  if ((A->width & 1)) {
    A->width++;
    incw = 1;
  }
#endif

  A->ncols = c;
  A->nrows = r;
  A->offset = 0;

  A->rows = (word**)m4ri_mmc_calloc(r + 1, sizeof(word*)); // We're overcomitting here.

  if(r && c) {
    /* we allow more than one malloc call so he have to be a bit clever
       here */
    
    int const bytes_per_row = A->width * sizeof(word);
    int const max_rows_per_block = __M4RI_MM_MAX_MALLOC / bytes_per_row;
    assert(max_rows_per_block);
    int rest = r % max_rows_per_block;
    
    int const nblocks = (rest == 0) ? r / max_rows_per_block : r / max_rows_per_block + 1;
    A->blocks = (mmb_t*)m4ri_mmc_calloc(nblocks + 1, sizeof(mmb_t));
    for(int i = 0; i < nblocks - 1; ++i) {
      A->blocks[i].size = __M4RI_MM_MAX_MALLOC;
      A->blocks[i].data = m4ri_mmc_calloc(1, __M4RI_MM_MAX_MALLOC);
      for(rci_t j = 0; j < max_rows_per_block; ++j)
      {
	int offset = A->width * j;				// Offset to start of row j within block.
	word *block_start = (word*)A->blocks[i].data;		// Start of block i.
	rci_t r = j + i * max_rows_per_block;
        A->rows[r] = block_start + offset;			// FIXME: rowstride candidate
      }
    }
    if(rest == 0)
      rest = max_rows_per_block;

    A->blocks[nblocks-1].size = rest * bytes_per_row;
    A->blocks[nblocks-1].data = m4ri_mmc_calloc(rest, bytes_per_row);
    for(rci_t j = 0; j < rest; ++j) {
      A->rows[j + max_rows_per_block * (nblocks - 1)] = (word*)(A->blocks[nblocks - 1].data) + A->width * j;
    }
#ifdef M4RI_WRAPWORD
    for (rci_t i = 0; i < A->nrows; ++i)
      word::init_array(A->rows[i], A->width);
#endif
  } else {
    A->blocks = NULL;
  }

#ifdef HAVE_SSE2
  if (incw) {
    A->width--;
  }
#endif

  return A;
}

mzd_t *mzd_init_window (mzd_t const *m, rci_t lowr, rci_t lowc, rci_t highr, rci_t highc) {
  rci_t nrows, ncols;
  mzd_t *window;
  window = (mzd_t*)m4ri_mmc_malloc(sizeof(mzd_t));

  nrows = MIN(highr - lowr, m->nrows - lowr);
  ncols = highc - lowc;
  
  window->ncols = ncols;
  window->nrows = nrows;

  window->offset = (lowc + m->offset) % m4ri_radix;
  wi_t const offset = (lowc + m->offset) / m4ri_radix;
  
  window->width = (ncols + window->offset) / m4ri_radix;
  if ((ncols + window->offset) % m4ri_radix)
    window->width++;
  window->blocks = NULL;

  if(nrows)
    window->rows = (word**)m4ri_mmc_calloc(nrows + 1, sizeof(word*));
  else
    window->rows = NULL;

  for(rci_t i = 0; i < nrows; ++i) {
    window->rows[i] = m->rows[lowr + i] + offset;
  }
  
  return window;
}


void mzd_free(mzd_t *A) {
  if(A->rows)
    m4ri_mmc_free(A->rows, (A->nrows + 1) * sizeof(word*));
  if(A->blocks) {
    int i;
    for(i = 0; A->blocks[i].size; ++i) {
      m4ri_mmc_free(A->blocks[i].data, A->blocks[i].size);
    }
    m4ri_mmc_free(A->blocks, (i + 1) * sizeof(mmb_t));
  }
  m4ri_mmc_free(A, sizeof(mzd_t));
}

void mzd_print( mzd_t const *M ) {
  char temp[SAFECHAR];
  for (rci_t i = 0; i < M->nrows; ++i) {
    printf("[");
    word *row = M->rows[i];
    if(M->offset == 0) {
      for (wi_t j = 0; j < M->width - 1; ++j) {
        m4ri_word_to_str(temp, row[j], 1);
        printf("%s|", temp);
      }
      row = row + M->width - 1;
      int const wide = (M->ncols % m4ri_radix) ? M->ncols % m4ri_radix : m4ri_radix;
      for (int j = 0; j < wide; ++j) {
        if(j != 0 && (j % 4) == 0)
          printf(":");
        if (__M4RI_GET_BIT(*row, j)) 
          printf("1");
        else
          printf(" ");
      }
    } else {
      for (rci_t j = 0; j < M->ncols; ++j) {
        if(j != 0 && (j % 4) == 0)
          printf(((j % m4ri_radix) == 0) ? "|" : ":");
        if(mzd_read_bit(M, i, j))
          printf("1");
        else
          printf(" ");
      }
    }
    printf("]\n");
  }
}

void mzd_print_tight(mzd_t const *M) {
  assert(M->offset == 0);

  char temp[SAFECHAR];
  word *row;

  for (rci_t i = 0; i < M->nrows; ++i) {
    printf("[");
    row = M->rows[i];
    for (wi_t j = 0; j < M->ncols / m4ri_radix; ++j) {
      m4ri_word_to_str(temp, row[j], 0);
      printf("%s", temp);
    }
    row = row + M->width - 1;
    for (int j = 0; j < (M->ncols % m4ri_radix); ++j) {
      printf("%d", __M4RI_GET_BIT(*row, j));
    }
    printf("]\n");
  }
}

void mzd_row_add(mzd_t *M, rci_t sourcerow, rci_t destrow) {
  mzd_row_add_offset(M, destrow, sourcerow, 0);
}

rci_t mzd_gauss_delayed(mzd_t *M, rci_t startcol, int full) {
  assert(M->offset == 0);

  rci_t startrow = startcol;
  rci_t pivots = 0;
  for (rci_t i = startcol; i < M->ncols ; ++i) {
    for(rci_t j = startrow ; j < M->nrows; ++j) {
      if (mzd_read_bit(M, j, i)) {
	mzd_row_swap(M, startrow, j);
	++pivots;

	for(rci_t ii = full ? 0 : startrow + 1;  ii < M->nrows; ++ii) {
	  if (ii != startrow) {
	    if (mzd_read_bit(M, ii, i)) {
	      mzd_row_add_offset(M, ii, startrow, i);
	    }
	  }
	}
	startrow = startrow + 1;
	break;
      }
    }
  }

  return pivots;
}

rci_t mzd_echelonize_naive(mzd_t *M, int full) { 
  return mzd_gauss_delayed(M, 0, full); 
}

/**
 * Transpose the 128 x 128-bit matrix SRC and write the result in DST.
 */
static inline mzd_t *_mzd_transpose_direct_128(mzd_t *DST, mzd_t const *SRC) {
  assert(DST->offset == 0);
  assert(SRC->offset == 0);

  /* we do one recursion level 
   * [AB] -> [AC]
   * [CD]    [BD]
   */
  for(int k = 0; k < 64; ++k)  {
    DST->rows[   k][0] = SRC->rows[   k][0]; //A
    DST->rows[64+k][0] = SRC->rows[   k][1]; //B
    DST->rows[   k][1] = SRC->rows[64+k][0]; //C
    DST->rows[64+k][1] = SRC->rows[64+k][1]; //D
  }

  /* now transpose each block A,B,C,D separately, cf. Hacker's Delight */
  word t[4];
  word m = __M4RI_CONVERT_TO_WORD(0xFFFFFFFF);
  for (int j = 32; j != 0; j = j >> 1, m = m ^ (m << j)) {
    for (int k = 0; k < 64; k = (k + j + 1) & ~j) {
      t[0] = ((DST->rows[k][0] >> j) ^ DST->rows[k+j][0]) & m;
      t[1] = ((DST->rows[k][1] >> j) ^ DST->rows[k+j][1]) & m;
      t[2] = ((DST->rows[64+k][0] >> j) ^ DST->rows[64+k+j][0]) & m;
      t[3] = ((DST->rows[64+k][1] >> j) ^ DST->rows[64+k+j][1]) & m;

      DST->rows[k][0] ^= t[0] << j;			// A
      DST->rows[k][1] ^= t[1] << j;			// C

      DST->rows[k+j][0] ^= t[0];			// A
      DST->rows[k+j][1] ^= t[1];			// C

      DST->rows[64+k][0] ^= t[2] << j;			// B
      DST->rows[64+k][1] ^= t[3] << j;			// D

      DST->rows[64+k+j][0] ^= t[2];			// B
      DST->rows[64+k+j][1] ^= t[3];			// D
    }
  }
  return DST;
}


static inline mzd_t *_mzd_transpose_direct(mzd_t *DST, mzd_t const *A) {
  if(A->offset || DST->offset) {
    for(rci_t i = 0; i < A->nrows; ++i) {
      for(rci_t j = 0; j < A->ncols; ++j) {
        mzd_write_bit(DST, j, i, mzd_read_bit(A, i, j));
      }
    }
    return DST;
  }

  if (A->nrows == 128 && A->ncols == 128 && m4ri_radix == 64) {
    _mzd_transpose_direct_128(DST, A);
    return DST;
  }

  int const spill = DST->ncols % m4ri_radix;
  int const have_incomplete_word = (spill != 0);			/* 0: all words are full; 1: last word is incomplete */
  wi_t const complete_words = DST->width - have_incomplete_word;
  //wi_t const rowdiff = A->rows[1] - A->rows[0];			/* Assume that the distance between every row is the same */
  for (rci_t i = 0; i < DST->nrows; ++i)
  {
    int const shift = i % m4ri_radix;
    wi_t const wordi = i / m4ri_radix;
    word *dstp = &DST->rows[i][complete_words + 1];			/* If there is no incomplete word, then the first k loop will be empty (k = -1) */
    int k = spill - 1;
    rci_t j = DST->ncols - spill;
    word *ap = &A->rows[j + k][wordi];
    word collect = 0;
    if (spill == 0)
    {
      k = m4ri_radix - 1;
      --dstp;
      j -= m4ri_radix;
    }
    /* Make k even... */
    else if ((spill & 1))
    {
      collect = ((*ap >> shift) & m4ri_one) << k;
      //ap -= rowdiff;
      --k;
      ap = &A->rows[j + k][wordi];
    }
    for (; j >= 0; j -= m4ri_radix)
    {
      /* ...so that we can unroll this loop a factor of two */
      for (; k > 0; k -=2)
      {
        collect |= ((*ap >> shift) & m4ri_one) << k;
	//ap -= rowdiff;
	ap = &A->rows[j + k - 1][wordi];
        collect |= ((*ap >> shift) & m4ri_one) << (k - 1);
	//ap -= rowdiff;
	//FIXME (this test is too slow, use rowstride)
	if (j > 0 || k > 1)
	  ap = &A->rows[j + k - 2][wordi];
      }
      k = m4ri_radix - 1;
      *--dstp = collect;
      collect = 0;
    }
  }
  return DST;
}

static inline mzd_t *_mzd_transpose(mzd_t *DST, mzd_t const *X) {
  assert(X->offset == 0);

  rci_t const nr = X->nrows;
  rci_t const nc = X->ncols;
  int const cutoff = 128; // must be >= 128.

  if(nr <= cutoff || nc <= cutoff) {
    mzd_t *x = mzd_copy(NULL, X);
    _mzd_transpose_direct(DST, x);
    mzd_free(x);
    return DST;
  }

  /* we cut at multiples of 128 if possible, otherwise at multiples of 64 */
  rci_t nr2 = (X->nrows > 256) ? 2 * m4ri_radix * (X->nrows / (4 * m4ri_radix)) : m4ri_radix * (X->nrows / (2 * m4ri_radix));
  rci_t nc2 = (X->ncols > 256) ? 2 * m4ri_radix * (X->ncols / (4 * m4ri_radix)) : m4ri_radix * (X->ncols / (2 * m4ri_radix));

  mzd_t *A = mzd_init_window(X,    0,   0, nr2, nc2);
  mzd_t *B = mzd_init_window(X,    0, nc2, nr2,  nc);
  mzd_t *C = mzd_init_window(X,  nr2,   0,  nr, nc2);
  mzd_t *D = mzd_init_window(X,  nr2, nc2,  nr,  nc);

  mzd_t *AT = mzd_init_window(DST,   0,   0, nc2, nr2);
  mzd_t *CT = mzd_init_window(DST,   0, nr2, nc2,  nr);
  mzd_t *BT = mzd_init_window(DST, nc2,   0,  nc, nr2);
  mzd_t *DT = mzd_init_window(DST, nc2, nr2,  nc,  nr);

  _mzd_transpose(AT, A);
  _mzd_transpose(BT, B);
  _mzd_transpose(CT, C);
  _mzd_transpose(DT, D);

  mzd_free_window(A); mzd_free_window(B);
  mzd_free_window(C); mzd_free_window(D);

  mzd_free_window(AT); mzd_free_window(CT);
  mzd_free_window(BT); mzd_free_window(DT);
  
  return DST;
}

mzd_t *mzd_transpose(mzd_t *DST, mzd_t const *A) {
  if (DST == NULL) {
    DST = mzd_init( A->ncols, A->nrows );
  } else {
    if (DST->nrows != A->ncols || DST->ncols != A->nrows) {
      m4ri_die("mzd_transpose: Wrong size for return matrix.\n");
    }
  }
  if(A->offset || DST->offset)
    return _mzd_transpose_direct(DST, A);
  else
    return _mzd_transpose(DST, A);
}

mzd_t *mzd_mul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  if (C == NULL) {
    C = mzd_init(A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols) {
      m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions.\n");
    }
  }
  if(B->ncols < m4ri_radix-10) { /* this cutoff is rather arbitrary */
    mzd_t *BT = mzd_transpose(NULL, B);
    _mzd_mul_naive(C, A, BT, 1);
    mzd_free (BT);
  } else {
    _mzd_mul_va(C, A, B, 1);
  }
  return C;
}

mzd_t *mzd_addmul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  if (C->nrows != A->nrows || C->ncols != B->ncols) {
    m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions.\n");
  }

  if(B->ncols < m4ri_radix-10) { /* this cutoff is rather arbitrary */
    mzd_t *BT = mzd_transpose(NULL, B);
    _mzd_mul_naive(C, A, BT, 0);
    mzd_free (BT);
  } else {
    _mzd_mul_va(C, A, B, 0);
  }
  return C;
}

mzd_t *_mzd_mul_naive(mzd_t *C, mzd_t const *A, mzd_t const *B, const int clear) {
  assert(A->offset == 0);
  assert(B->offset == 0);
  assert(C->offset == 0);
  wi_t eol;
  word *a, *b, *c;

  if (clear) {
    word const mask_end = __M4RI_LEFT_BITMASK(C->ncols % m4ri_radix);
    /* improves performance on x86_64 but is not cross plattform */
    /* asm __volatile__ (".p2align 4\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop"); */
    for (rci_t i = 0; i < C->nrows; ++i) {
      wi_t j = 0;
      for (; j < C->width - 1; ++j) {
  	C->rows[i][j] = 0;
      }
      C->rows[i][j] &= ~mask_end;
    }
  }

  if(C->ncols % m4ri_radix) {
    eol = (C->width - 1);
  } else {
    eol = (C->width);
  }

  word parity[64];
  for (int i = 0; i < 64; ++i) {
    parity[i] = 0;
  }
  wi_t const wide = A->width;
  int const blocksize = __M4RI_MZD_MUL_BLOCKSIZE;
  for (rci_t start = 0; start + blocksize <= C->nrows; start += blocksize) {
    for (rci_t i = start; i < start + blocksize; ++i) {
      a = A->rows[i];
      c = C->rows[i];
      for (rci_t j = 0; j < m4ri_radix * eol; j += m4ri_radix) {
	for (int k = 0; k < m4ri_radix; ++k) {
          b = B->rows[j + k];
          parity[k] = a[0] & b[0];
          for (wi_t ii = wide - 1; ii >= 1; --ii)
	    parity[k] ^= a[ii] & b[ii];
        }
        c[j / m4ri_radix] ^= m4ri_parity64(parity);
      }
      
      if (eol != C->width) {
	word const mask_end = __M4RI_LEFT_BITMASK(C->ncols % m4ri_radix);
        /* improves performance on x86_64 but is not cross plattform */
	/* asm __volatile__ (".p2align 4\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop"); */
        for (int k = 0; k < (C->ncols % m4ri_radix); ++k) {
          b = B->rows[m4ri_radix * eol + k];
          parity[k] = a[0] & b[0];
          for (wi_t ii = 1; ii < A->width; ++ii)
            parity[k] ^= a[ii] & b[ii];
        }
        c[eol] ^= m4ri_parity64(parity) & mask_end;
      }
    }
  }

  for (rci_t i = C->nrows - (C->nrows % blocksize); i < C->nrows; ++i) {
    a = A->rows[i];
    c = C->rows[i];
    for (rci_t j = 0; j < m4ri_radix * eol; j += m4ri_radix) {
      for (int k = 0; k < m4ri_radix; ++k) {
        b = B->rows[j+k];
        parity[k] = a[0] & b[0];
        for (wi_t ii = wide - 1; ii >= 1; --ii)
          parity[k] ^= a[ii] & b[ii];
      }
      c[j/m4ri_radix] ^= m4ri_parity64(parity);
    }
    
    if (eol != C->width) {
      word const mask_end = __M4RI_LEFT_BITMASK(C->ncols % m4ri_radix);
      /* improves performance on x86_64 but is not cross plattform */
      /* asm __volatile__ (".p2align 4\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop\n\tnop"); */
      for (int k = 0; k < (C->ncols % m4ri_radix); ++k) {
        b = B->rows[m4ri_radix * eol + k];
        parity[k] = a[0] & b[0];
        for (wi_t ii = 1; ii < A->width; ++ii)
          parity[k] ^= a[ii] & b[ii];
      }
      c[eol] ^= m4ri_parity64(parity) & mask_end;
    }
  }

  return C;
}

mzd_t *_mzd_mul_va(mzd_t *C, mzd_t const *v, mzd_t const *A, int const clear) {
  assert(C->offset == 0);
  assert(A->offset == 0);
  assert(v->offset == 0);

  if(clear)
    mzd_set_ui(C, 0);

  rci_t const m = v->nrows;
  rci_t const n = v->ncols;
  
  for(rci_t i = 0; i < m; ++i)
    for(rci_t j = 0; j < n; ++j)
      if (mzd_read_bit(v,i,j))
        mzd_combine(C,i,0, C,i,0, A,j,0);
  return C;
}

void mzd_randomize(mzd_t *A) {
  wi_t const width = A->width - 1;
  int const offset = A->offset;
  if(offset) {
    if(width == 0) {
      word const mask = __M4RI_MIDDLE_BITMASK(A->ncols, offset);
      for(rci_t i = 0; i < A->nrows; ++i)
	A->rows[i][0] ^= (A->rows[i][0] ^ (m4ri_random_word() << offset)) & mask;
    } else {
      word const mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - offset);
      word const mask_end = __M4RI_LEFT_BITMASK((A->ncols + offset) % m4ri_radix);
      int const need_last_bits = ((m4ri_one << offset) & mask_end) != 0;
      for(rci_t i = 0; i < A->nrows; ++i) {
	word prev_random_word;
	word random_word = m4ri_random_word();
	A->rows[i][0] ^= (A->rows[i][0] ^ (random_word << offset)) & mask_begin;
	for(wi_t j = 1; j < width; ++j) {
	  prev_random_word = random_word;
	  random_word = m4ri_random_word();
	  A->rows[i][j] = (random_word << offset) | (prev_random_word >> (m4ri_radix - offset));
	}
	prev_random_word = random_word;
	random_word = 0;
	if (need_last_bits)
	  random_word = m4ri_random_word();
	A->rows[i][width] ^= (A->rows[i][width] ^ ((random_word << offset) | (prev_random_word >> (m4ri_radix - offset)))) & mask_end;
      }
    }
  } else {
    word const mask_end = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
    for(rci_t i = 0; i < A->nrows; ++i) {
      for(wi_t j = 0; j < width; ++j)
	A->rows[i][j] = m4ri_random_word();
      A->rows[i][width] ^= (A->rows[i][width] ^ m4ri_random_word()) & mask_end;
    }
  }
}

void mzd_set_ui( mzd_t *A, unsigned int value) {
  word const mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->offset);
  word const mask_end = __M4RI_LEFT_BITMASK((A->ncols + A->offset) % m4ri_radix);
  
  if(A->width == 1) {
    for(rci_t i = 0; i < A->nrows; ++i) {
      for(rci_t j = 0 ; j < A->ncols; ++j)
        mzd_write_bit(A,i,j, 0);
    }
  } else {
    for (rci_t i = 0; i < A->nrows; ++i) {
      word *row = A->rows[i];
      row[0] &= ~mask_begin;
      for(wi_t j = 1; j < A->width - 1; ++j)
        row[j] = 0;
      row[A->width - 1] &= ~mask_end;
    }
  }

  if(value % 2 == 0)
    return;

  rci_t const stop = MIN(A->nrows, A->ncols);
  for (rci_t i = 0; i < stop; ++i) {
    mzd_write_bit(A, i, i, 1);
  }
}

int mzd_equal(mzd_t const *A, mzd_t const *B) {
  if (A->nrows != B->nrows) return FALSE;
  if (A->ncols != B->ncols) return FALSE;
  if (A == B) return TRUE;

  wi_t const width = A->width - 1;

  if (A->offset == B->offset) {
    int const non_zero_offset = (A->offset != 0);
    if (non_zero_offset < width) {
      for (rci_t i = 0; i < A->nrows; ++i) {
	for (wi_t j = non_zero_offset; j < width; ++j) {
	  if (A->rows[i][j] != B->rows[i][j])
	    return FALSE;
	}
      }
    }
    word const mask_end = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
    if (non_zero_offset) {
      word mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->offset);
      if (!width)
	mask_begin &= mask_end;
      for (rci_t i = 0; i < A->nrows; ++i) {
	if (((A->rows[i][0] ^ B->rows[i][0]) & mask_begin))
	  return FALSE;
      }
      if (!width)
	return TRUE;
    }
    for (rci_t i = 0; i < A->nrows; ++i) {
      if (((A->rows[i][width] ^ B->rows[i][width]) & mask_end))
	return FALSE;
    }
  } else {
    int shift = B->offset - A->offset;
    if (shift < 0) {
      mzd_t const *tmp = A;
      A = B;
      B = tmp;
      shift = -shift;
    }
    int const non_zero_offset = (A->offset != 0);
    if (non_zero_offset < width) {
      for (rci_t i = 0; i < A->nrows; ++i) {
	for (wi_t j = non_zero_offset; j < width; ++j) {
	  word Bval = (B->rows[i][j] >> shift) | (B->rows[i][j + 1] << (m4ri_radix - shift));
	  if (A->rows[i][j] != Bval)
	    return FALSE;
	}
      }
    }
    word const mask_end = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
    if (non_zero_offset) {
      word mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->offset);
      if (!width)
	mask_begin &= mask_end;
      if (1 < B->width) {
	for (rci_t i = 0; i < A->nrows; ++i) {
	  word Bval = (B->rows[i][0] >> shift) | (B->rows[i][1] << (m4ri_radix - shift));
	  if (((A->rows[i][0] ^ Bval) & mask_begin))
	    return FALSE;
	}
      } else {
	for (rci_t i = 0; i < A->nrows; ++i) {
	  word Bval = B->rows[i][0] >> shift;
	  if (((A->rows[i][0] ^ Bval) & mask_begin))
	    return FALSE;
	}
      }
      if (!width)
	return TRUE;
    }
    if (width + 1 < B->width)
    {
      for (rci_t i = 0; i < A->nrows; ++i) {
	word Bval = (B->rows[i][width] >> shift) | (B->rows[i][width + 1] << (m4ri_radix - shift));
	if (((A->rows[i][width] ^ Bval) & mask_end))
	  return FALSE;
      }
    } else {
      for (rci_t i = 0; i < A->nrows; ++i) {
	word Bval = B->rows[i][width] >> shift;
	if (((A->rows[i][width] ^ Bval) & mask_end))
	  return FALSE;
      }
    }
  }
  return TRUE;
}

int mzd_cmp(mzd_t const *A, mzd_t const *B) {
  assert(A->offset == 0);
  assert(B->offset == 0);

  if(A->nrows < B->nrows) return -1;
  if(B->nrows < A->nrows) return 1;
  if(A->ncols < B->ncols) return -1;
  if(B->ncols < A->ncols) return 1;

  /* Columns with large index are "larger", but rows with small
     index are more important than with large index. */
  for(rci_t i = 0; i < A->nrows; ++i) {
    for(wi_t j = A->width - 1; j >= 0; --j) {
      if (__M4RI_CONVERT_TO_UINT64_T(A->rows[i][j]) < __M4RI_CONVERT_TO_UINT64_T(B->rows[i][j]))
	return -1;
      else if (__M4RI_CONVERT_TO_UINT64_T(A->rows[i][j]) > __M4RI_CONVERT_TO_UINT64_T(B->rows[i][j]))
	return 1;
    }
  }
  return 0;
}

void mzd_copy_row_weird_to_even(mzd_t *B, rci_t i, mzd_t const *A, rci_t j);

mzd_t *mzd_copy(mzd_t *N, mzd_t const *P) {
  if (N == P)
    return N;

  if (!P->offset){
    if (N == NULL) {
      N = mzd_init(P->nrows, P->ncols);
    } else {
      if (N->nrows < P->nrows || N->ncols < P->ncols)
	m4ri_die("mzd_copy: Target matrix is too small.");
    }
    word *p_truerow, *n_truerow;
    wi_t const wide = P->width - 1;
    word mask = __M4RI_LEFT_BITMASK(P->ncols % m4ri_radix);
    for (rci_t i = 0; i < P->nrows; ++i) {
      p_truerow = P->rows[i];
      n_truerow = N->rows[i];
      for (wi_t j = 0; j < wide; ++j)
        n_truerow[j] = p_truerow[j];
      n_truerow[wide] = (n_truerow[wide] & ~mask) | (p_truerow[wide] & mask);
    }
  } else { // P->offset > 0
    if (N == NULL) {
      N = mzd_init(P->nrows, P->ncols+ P->offset);
      N->ncols -= P->offset;
      N->offset = P->offset;
      N->width=P->width;
    } else {
      if (N->nrows < P->nrows || N->ncols < P->ncols)
	m4ri_die("mzd_copy: Target matrix is too small.");
    }
    if(N->offset == P->offset) {
      for(rci_t i = 0; i < P->nrows; ++i) {
        mzd_copy_row(N, i, P, i);
      }
    } else if(N->offset == 0) {
      for(rci_t i = 0; i < P->nrows; ++i) {
        mzd_copy_row_weird_to_even(N, i, P, i);
      }
    } else {
      m4ri_die("mzd_copy: completely unaligned copy not implemented yet.");
    }
  }
  return N;
}

/* This is sometimes called augment */
mzd_t *mzd_concat(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  assert(A->offset == 0);
  assert(B->offset == 0);
  
  if (A->nrows != B->nrows) {
    m4ri_die("mzd_concat: Bad arguments to concat!\n");
  }

  if (C == NULL) {
    C = mzd_init(A->nrows, A->ncols + B->ncols);
  } else if (C->nrows != A->nrows || C->ncols != (A->ncols + B->ncols)) {
    m4ri_die("mzd_concat: C has wrong dimension!\n");
  }

  for (rci_t i = 0; i < A->nrows; ++i) {
    word *dst_truerow = C->rows[i];
    word *src_truerow = A->rows[i];
    for (wi_t j = 0; j < A->width; ++j) {
      dst_truerow[j] = src_truerow[j];
    }
  }

  for (rci_t i = 0; i < B->nrows; ++i) {
    for (rci_t j = 0; j < B->ncols; ++j) {
      mzd_write_bit(C, i, j + A->ncols, mzd_read_bit(B, i, j));
    }
  }

  return C;
}

mzd_t *mzd_stack(mzd_t *C, mzd_t const *A, mzd_t const *B) {
  assert(A->offset == 0);
  assert(B->offset == 0);

  if (A->ncols != B->ncols) {
    m4ri_die("mzd_stack: A->ncols (%d) != B->ncols (%d)!\n", A->ncols, B->ncols);
  }

  if (C == NULL) {
    C = mzd_init(A->nrows + B->nrows, A->ncols);
  } else if (C->nrows != (A->nrows + B->nrows) || C->ncols != A->ncols) {
    m4ri_die("mzd_stack: C has wrong dimension!\n");
  }
  
  for(rci_t i = 0; i < A->nrows; ++i) {
    word *src_truerow = A->rows[i];
    word *dst_truerow = C->rows[i];
    for (wi_t j = 0; j < A->width; ++j) {
      dst_truerow[j] = src_truerow[j]; 
    }
  }

  for(rci_t i = 0; i < B->nrows; ++i) {
    word *dst_truerow = C->rows[A->nrows + i];
    word *src_truerow = B->rows[i];
    for (wi_t j = 0; j < B->width; ++j) {
      dst_truerow[j] = src_truerow[j]; 
    }
  }
  return C;
}

mzd_t *mzd_invert_naive(mzd_t *INV, mzd_t const *A, mzd_t const *I) {
  assert(A->offset == 0);
  mzd_t *H;

  H = mzd_concat(NULL, A, I);

  rci_t x = mzd_echelonize_naive(H, TRUE);

  if (x == 0) { 
    mzd_free(H); 
    return NULL; 
  }
  
  INV = mzd_submatrix(INV, H, 0, A->ncols, A->nrows, 2 * A->ncols);

  mzd_free(H);
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

  if (C == B) { //swap
    mzd_t const *tmp = A;
    A = B;
    B = tmp;
  }

  if (C->offset | A->offset | B->offset) {
    for(rci_t i = 0; i < nrows; ++i) {
      mzd_combine_weird(C,i,0, A,i,0, B,i,0);
    }
    return C;
  }

  word const mask_end = __M4RI_LEFT_BITMASK(C->ncols%m4ri_radix);

  switch(A->width) {
  case 0:
    return C;
  case 1:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] ^= ((A->rows[i][0] ^ B->rows[i][0] ^ C->rows[i][0]) & mask_end);
    }
    break;
  case 2:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] ^= ((A->rows[i][1] ^ B->rows[i][1] ^ C->rows[i][1]) & mask_end);
    }
    break;
  case 3:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] = A->rows[i][1] ^ B->rows[i][1];
      C->rows[i][2] ^= ((A->rows[i][2] ^ B->rows[i][2] ^ C->rows[i][2]) & mask_end);
    }
    break;
  case 4:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] = A->rows[i][1] ^ B->rows[i][1];
      C->rows[i][2] = A->rows[i][2] ^ B->rows[i][2];
      C->rows[i][3] ^= ((A->rows[i][3] ^ B->rows[i][3] ^ C->rows[i][3]) & mask_end);
    }
    break;
  case 5:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] = A->rows[i][1] ^ B->rows[i][1];
      C->rows[i][2] = A->rows[i][2] ^ B->rows[i][2];
      C->rows[i][3] = A->rows[i][3] ^ B->rows[i][3];
      C->rows[i][4] ^= ((A->rows[i][4] ^ B->rows[i][4] ^ C->rows[i][4]) & mask_end);
    }
    break;
  case 6:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] = A->rows[i][1] ^ B->rows[i][1];
      C->rows[i][2] = A->rows[i][2] ^ B->rows[i][2];
      C->rows[i][3] = A->rows[i][3] ^ B->rows[i][3];
      C->rows[i][4] = A->rows[i][4] ^ B->rows[i][4];
      C->rows[i][5] ^= ((A->rows[i][5] ^ B->rows[i][5] ^ C->rows[i][5]) & mask_end);
    }
    break;
  case 7:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] = A->rows[i][1] ^ B->rows[i][1];
      C->rows[i][2] = A->rows[i][2] ^ B->rows[i][2];
      C->rows[i][3] = A->rows[i][3] ^ B->rows[i][3];
      C->rows[i][4] = A->rows[i][4] ^ B->rows[i][4];
      C->rows[i][5] = A->rows[i][5] ^ B->rows[i][5];
      C->rows[i][6] ^= ((A->rows[i][6] ^ B->rows[i][6] ^ C->rows[i][6]) & mask_end);
    }
    break;
  case 8:
    for(rci_t i = 0; i < nrows; ++i) {
      C->rows[i][0] = A->rows[i][0] ^ B->rows[i][0];
      C->rows[i][1] = A->rows[i][1] ^ B->rows[i][1];
      C->rows[i][2] = A->rows[i][2] ^ B->rows[i][2];
      C->rows[i][3] = A->rows[i][3] ^ B->rows[i][3];
      C->rows[i][4] = A->rows[i][4] ^ B->rows[i][4];
      C->rows[i][5] = A->rows[i][5] ^ B->rows[i][5];
      C->rows[i][6] = A->rows[i][6] ^ B->rows[i][6];
      C->rows[i][7] ^= ((A->rows[i][7] ^ B->rows[i][7] ^ C->rows[i][7]) & mask_end);
    }
    break;

  default:
    for(rci_t i = 0; i < nrows; ++i) {
      mzd_combine_even(C,i,0, A,i,0, B,i,0);
    }
  }
  return C;
}

mzd_t *mzd_submatrix(mzd_t *S, mzd_t const *M, rci_t const startrow, rci_t const startcol, rci_t const endrow, rci_t const endcol) {
  rci_t const nrows = endrow - startrow;
  rci_t const ncols = endcol - startcol;

  if (S == NULL) {
    S = mzd_init(nrows, ncols);
  } else if(S->nrows < nrows || S->ncols < ncols) {
    m4ri_die("mzd_submatrix: got S with dimension %d x %d but expected %d x %d\n", S->nrows, S->ncols, nrows, ncols);
  }
  assert(M->offset == S->offset);

  wi_t const startword = (startcol + M->offset) / m4ri_radix;

  /* we start at the beginning of a word */
  if ((startcol + M->offset) % m4ri_radix == 0) {
    if(ncols / m4ri_radix != 0) {
      for(rci_t x = startrow, i = 0; i < nrows; ++i, ++x) {
        memcpy(S->rows[i], M->rows[x] + startword, sizeof(word) * (ncols / m4ri_radix));
      }
    }
    if (ncols % m4ri_radix) {
      word const mask_end = __M4RI_LEFT_BITMASK(ncols % m4ri_radix);
      for(rci_t x = startrow, i = 0; i < nrows; ++i, ++x) {
        /* process remaining bits */
	word temp = M->rows[x][startword + ncols / m4ri_radix] & mask_end;
	S->rows[i][ncols / m4ri_radix] = temp;
      } 
    }
    /* startcol is not the beginning of a word */
  } else { 
    int const spot = (startcol + M->offset) % m4ri_radix;
    for(rci_t x = startrow, i = 0; i < nrows; ++i, ++x) {
      word *truerow = M->rows[x];

      /* process full words first */
      for(wi_t colword = 0; colword < (ncols / m4ri_radix); ++colword) {
	wi_t block = colword + startword;
	word temp = (truerow[block] >> spot) | (truerow[block + 1] << (m4ri_radix - spot)); 
	S->rows[i][colword] = temp;
      }
      /* process remaining bits (lazy) */
      wi_t colword = ncols / m4ri_radix;
      for (int y = 0; y < ncols % m4ri_radix; ++y) {
	BIT bit = mzd_read_bit(M, x, startcol + colword * m4ri_radix + y);
	mzd_write_bit(S, i, colword * m4ri_radix + y, bit);
      }
    }
  }
  return S;
}

void mzd_combine(mzd_t *C,       rci_t const c_row, wi_t const c_startblock,
		 mzd_t const *A, rci_t const a_row, wi_t const a_startblock, 
		 mzd_t const *B, rci_t const b_row, wi_t const b_startblock) {

  /** 
   * \todo respect ncols at the end 
   */

  if(C->offset | A->offset | B->offset) {
    mzd_combine_weird(C, c_row, c_startblock, A, a_row, a_startblock, B, b_row, b_startblock);
    return;
  }


  if( C == A && a_row == c_row && a_startblock == c_startblock) {
    mzd_combine_even_in_place(C, c_row, c_startblock, B, b_row, b_startblock);
    return;
  }

  mzd_combine_even(C, c_row, c_startblock, A, a_row, a_startblock, B, b_row, b_startblock);
  return;
}

void mzd_col_swap(mzd_t *M, rci_t const cola, rci_t const colb) {
  if (cola == colb)
    return;

  rci_t const _cola = cola + M->offset;
  rci_t const _colb = colb + M->offset;

  wi_t const a_word = _cola / m4ri_radix;
  wi_t const b_word = _colb / m4ri_radix;
  int const a_bit = _cola % m4ri_radix;
  int const b_bit = _colb % m4ri_radix;

  if(a_word == b_word) {
    for (rci_t i = 0; i < M->nrows; ++i) {
      word *base = (M->rows[i] + a_word);
      register word b = *base;
      register word x = ((b >> a_bit) ^ (b >> b_bit)) & m4ri_one; // XOR temporary
      *base = b ^ ((x << a_bit) | (x << b_bit));
    }
    return;
  }

  word const a_bm = m4ri_one << a_bit;
  word const b_bm = m4ri_one << b_bit;

  if(a_bit > b_bit) {
    int const offset = a_bit - b_bit;

    for (rci_t i = 0; i < M->nrows; ++i) {
      word *base = M->rows[i];
      word a = *(base + a_word);
      word b = *(base + b_word);

      a ^= (b & b_bm) << offset;
      b ^= (a & a_bm) >> offset;
      a ^= (b & b_bm) << offset;

      *(base + a_word) = a;
      *(base + b_word) = b;
    }
  } else {
    int const offset = b_bit - a_bit;
    for (rci_t i = 0; i < M->nrows; ++i) {
      word *base = M->rows[i];
      word a = *(base + a_word);
      word b = *(base + b_word);

      a ^= (b & b_bm) >> offset;
      b ^= (a & a_bm) << offset;
      a ^= (b & b_bm) >> offset;
      *(base + a_word) = a;
      *(base + b_word) = b;
    }
  }

}


int mzd_is_zero(mzd_t const *A) {
  /* Could be improved: stopping as the first non zero value is found (status!=0) */
  rci_t const mb = A->nrows;
  rci_t const nb = A->ncols;
  int const Aoffset = A->offset;
  int const nbrest = (nb + Aoffset) % m4ri_radix;
  word status = 0;
  if (nb + Aoffset >= m4ri_radix) {
    // Large A
    word mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix-Aoffset);
    if (Aoffset == 0)
      mask_begin = ~mask_begin;
    word mask_end = __M4RI_LEFT_BITMASK(nbrest);
    for (rci_t i = 0; i < mb; ++i) {
        status |= A->rows[i][0] & mask_begin;
        for (wi_t j = 1; j < A->width - 1; ++j)
            status |= A->rows[i][j];
        status |= A->rows[i][A->width - 1] & mask_end;
    }
  } else {
    // Small A
    word mask = __M4RI_MIDDLE_BITMASK(nb,Aoffset);
    for (rci_t i = 0; i < mb; ++i) {
      status |= A->rows[i][0] & mask;
    }
  }
  
  return !status;
}

void mzd_copy_row_weird_to_even(mzd_t *B, rci_t i, mzd_t const *A, rci_t j) {
  assert(B->offset == 0);
  assert(B->ncols >= A->ncols);

  word *b = B->rows[j];

  int const rest = A->ncols % m4ri_radix;

  rci_t c;
  for(c = 0; c + m4ri_radix <= A->ncols; c += m4ri_radix) {
    b[c / m4ri_radix] = mzd_read_bits(A, i, c, m4ri_radix);
  }
  if (rest) {
    word const temp = mzd_read_bits(A, i, c, rest);
    b[c / m4ri_radix] &= __M4RI_LEFT_BITMASK(m4ri_radix - rest);
    b[c / m4ri_radix] |= temp;
  }
}

void mzd_copy_row(mzd_t *B, rci_t i, mzd_t const *A, rci_t j) {
  assert(B->offset == A->offset);
  assert(B->ncols >= A->ncols);
  wi_t const width = MIN(B->width, A->width) - 1;

  word *a = A->rows[j];
  word *b = B->rows[i];
 
  word const mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->offset);
  word const mask_end = __M4RI_LEFT_BITMASK((A->ncols + A->offset) % m4ri_radix);

  if (width != 0) {
    b[0] = (b[0] & ~mask_begin) | (a[0] & mask_begin);
    for(wi_t k = 1; k < width; ++k)
      b[k] = a[k];
    b[width] = (b[width] & ~mask_end) | (a[width] & mask_end);
    
  } else {
    b[0] = (b[0] & ~mask_begin) | (a[0] & mask_begin & mask_end) | (b[0] & ~mask_end);
  }
}


void mzd_row_clear_offset(mzd_t *M, rci_t row, rci_t coloffset) {
  coloffset += M->offset;
  wi_t const startblock = coloffset / m4ri_radix;
  word temp;

  /* make sure to start clearing at coloffset */
  if (coloffset%m4ri_radix) {
    temp = M->rows[row][startblock];
    temp &= __M4RI_RIGHT_BITMASK(m4ri_radix - coloffset);
  } else {
    temp = 0;
  }
  M->rows[row][startblock] = temp;
  for (wi_t i = startblock + 1; i < M->width; ++i) {
    M->rows[row][i] = 0;
  }
}


int mzd_find_pivot(mzd_t const *A, rci_t start_row, rci_t start_col, rci_t *r, rci_t *c) { 
  assert(A->offset == 0);
  rci_t const nrows = A->nrows;
  rci_t const ncols = A->ncols;
  word data = 0;
  rci_t row_candidate = 0;
  if(A->ncols - start_col < m4ri_radix) {
    for(rci_t j = start_col; j < A->ncols; j += m4ri_radix) {
      int const length = MIN(m4ri_radix, ncols - j);
      for(rci_t i = start_row; i < nrows; ++i) {
        word const curr_data = mzd_read_bits(A, i, j, length);
        if (m4ri_lesser_LSB(curr_data, data)) {
          row_candidate = i;
          data = curr_data;
        }
      }
      if(data) {
        *r = row_candidate;
        for(int l = 0; l < length; ++l) {
          if(__M4RI_GET_BIT(data, l)) {
            *c = j + l;
            break;
          }
        }
        return 1;
      }
    }
  } else {
    /* we definitely have more than one word */
    /* handle first word */
    int const bit_offset = (start_col % m4ri_radix);
    wi_t const word_offset = start_col / m4ri_radix;
    word const mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix-bit_offset);
    for(rci_t i = start_row; i < nrows; ++i) {
      word const curr_data = A->rows[i][word_offset] & mask_begin;
      if (m4ri_lesser_LSB(curr_data, data)) {
        row_candidate = i;
        data = curr_data;
        if(__M4RI_GET_BIT(data,bit_offset)) {
          break;
        }
      }
    }
    if(data) {
      *r = row_candidate;
      data >>= bit_offset;
      assert(data);
      for(int l = 0; l < (m4ri_radix - bit_offset); ++l) {
        if(__M4RI_GET_BIT(data, l)) {
          *c = start_col + l;
          break;
        }
      }
      return 1;
    }
    /* handle complete words */
    for(wi_t wi = word_offset + 1; wi < A->width - 1; ++wi) {
      for(rci_t i = start_row; i < nrows; ++i) {
        word const curr_data = A->rows[i][wi];
        if (m4ri_lesser_LSB(curr_data, data)) {
          row_candidate = i;
          data = curr_data;
          if(__M4RI_GET_BIT(data, 0))
            break;
        }
      }
      if(data) {
        *r = row_candidate;
        for(int l = 0; l < m4ri_radix; ++l) {
          if(__M4RI_GET_BIT(data, l)) {
            *c = wi * m4ri_radix + l;
            break;
          }
        }
        return 1;
      }
    }
    /* handle last word */
    int const end_offset = (A->ncols % m4ri_radix) ? (A->ncols % m4ri_radix) : m4ri_radix;
    word const mask_end = __M4RI_LEFT_BITMASK(end_offset % m4ri_radix);
    wi_t wi = A->width - 1;
    for(rci_t i = start_row; i < nrows; ++i) {
      word const curr_data = A->rows[i][wi] & mask_end;
      if (m4ri_lesser_LSB(curr_data, data)) {
        row_candidate = i;
        data = curr_data;
        if(__M4RI_GET_BIT(data,0))
          break;
      }
    }
    if(data) {
      *r = row_candidate;
      for(int l = 0; l < end_offset; ++l) {
        if(__M4RI_GET_BIT(data, l)) {
          *c = wi * m4ri_radix + l;
          break;
        }
      }
      return 1;
    }
  }
  return 0;
}


#define MASK(c)    (((uint64_t)(-1)) / (__M4RI_TWOPOW(__M4RI_TWOPOW(c)) + 1))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (__M4RI_TWOPOW(c))) & MASK(c))

static inline int m4ri_bitcount(word w)  {
   uint64_t n = __M4RI_CONVERT_TO_UINT64_T(w);
   n = COUNT(n, 0);
   n = COUNT(n, 1);
   n = COUNT(n, 2);
   n = COUNT(n, 3);
   n = COUNT(n, 4);
   n = COUNT(n, 5);
   return (int)n;
}


double _mzd_density(mzd_t const *A, wi_t res, rci_t r, rci_t c) {
  size_t count = 0;
  size_t total = 0;
  
  if(A->width == 1) {
    for(rci_t i = r; i < A->nrows; ++i)
      for(rci_t j = c; j < A->ncols; ++j)
        if(mzd_read_bit(A, i, j))
          ++count;
    return ((double)count)/(1.0 * A->ncols * A->nrows);
  }

  if(res == 0)
    res = A->width / 100;
  if (res < 1)
    res = 1;

  for(rci_t i = r; i < A->nrows; ++i) {
    word *truerow = A->rows[i];
    for(rci_t j = c; j < m4ri_radix-A->offset; ++j)
      if(mzd_read_bit(A, i, j))
        ++count;
    total += m4ri_radix - A->offset;

    for(wi_t j = MAX(1, (c + A->offset) / m4ri_radix); j < A->width - 1; j += res) {
      count += m4ri_bitcount(truerow[j]);
      total += m4ri_radix;
    }
    for(int j = 0; j < (A->ncols + A->offset) % m4ri_radix; ++j)
      if(mzd_read_bit(A, i, m4ri_radix * ((A->ncols + A->offset) / m4ri_radix) + j))
        ++count;
    total += (A->ncols + A->offset) % m4ri_radix;
  }

  return (double)count / total;
}

double mzd_density(mzd_t const *A, wi_t res) {
  return _mzd_density(A, res, 0, 0);
}

rci_t mzd_first_zero_row(mzd_t const *A) {
  word const mask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->offset);
  word const mask_end = __M4RI_LEFT_BITMASK((A->ncols + A->offset) % m4ri_radix);
  wi_t const end = A->width - 1;
  word *row;

  for(rci_t i = A->nrows - 1; i >= 0; --i) {
    row = A->rows[i];
    word tmp = row[0] & mask_begin;
    for (wi_t j = 1; j < end; ++j)
      tmp |= row[j];
    tmp |= row[end] & mask_end;
    if(tmp)
      return i + 1;
  }
  return 0;
}
