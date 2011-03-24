/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008-2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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

#include <assert.h>

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include "pls_mmpf.h"
#include "brilliantrussian.h"
#include "grayflex.h"
#include "xor.h"

static inline rci_t _max_value(rci_t *data, int length) {
  rci_t max = 0;
  for(int i = 0; i < length; ++i) {
    max = MAX(max, data[i]);
  }
  return max;
}

int _mzd_pls_submatrix(mzd_t *A, 
                       rci_t const start_row, rci_t const stop_row, 
                       rci_t const start_col, int const k, 
                       mzp_t *P, mzp_t *Q,
                       rci_t *done, rci_t *done_row, wi_t const splitblock)  {
  word bm[4 * MAXKAY];
  wi_t os[4 * MAXKAY];

  /* we're essentially constructing a submatrix but cheaply */
  wi_t const width = A->width;
  rci_t const ncols = A->ncols;

  if (A->width > splitblock) {
    A->width = splitblock;
    A->ncols = splitblock * RADIX;
  }

  int curr_pos;
  for(curr_pos = 0; curr_pos < k; ++curr_pos) {
    os[curr_pos] = (start_col + curr_pos) / RADIX;
    bm[curr_pos] = ONE << ((start_col + curr_pos) % RADIX);
    int found = 0;
    /* search for some pivot */
    rci_t i;
    for(i = start_row + curr_pos; i < stop_row; ++i) {
      word const tmp = mzd_read_bits(A, i, start_col, curr_pos + 1);
      if(tmp) {
        wordPtr Arow = A->rows[i];
        /* clear before but preserve transformation matrix */
        for (rci_t l = 0; l < curr_pos; ++l)
          if(done[l.val()] < i) {
            if((Arow[os[l.val()]] & bm[l.val()]))
              mzd_row_add_offset(A, i, start_row + l, start_col + l + 1);
            done[l.val()] = i; /* encode up to which row we added for l already */
          }
        if(mzd_read_bit(A, i, start_col + curr_pos)) {
          found = 1;
          break;
        }
      }
    }
    if(!found) {
      break;
    }

    P->values[(start_row + curr_pos).val()] = i;
    mzd_row_swap(A, i, start_row + curr_pos);

    Q->values[(start_row + curr_pos).val()] = start_col + curr_pos;
    done[curr_pos] = i;
  }
  
  /* finish submatrix */
  *done_row = _max_value(done, curr_pos);
  for(rci_t c2 = 0; c2 < curr_pos && start_col + c2 < A->ncols -1; ++c2)
    for(rci_t r2 = done[c2.val()] + 1; r2 <= *done_row; ++r2)
      if(mzd_read_bit(A, r2, start_col + c2))
        mzd_row_add_offset(A, r2, start_row + c2, start_col + c2 + 1);

  /* reset to original size */
  A->ncols = ncols;
  A->width = width;

  return curr_pos;
}

/* create a table of all 2^k linear combinations */
void mzd_make_table_pls(mzd_t *M, rci_t r, rci_t c, int k, mzd_t *T, rci_t *Le, rci_t *Lm) {
  assert(T->blocks[1].size == 0);
  wi_t const blockoffset= c / RADIX;
  int const twokay= TWOPOW(k);
  wi_t const wide = T->width - blockoffset;
  wi_t const count = (wide + 7U) / 8U;
  int const entry_point = wide % 8U;

  wordPtr ti, ti1, m;

  ti1 = T->rows[0] + blockoffset;
  ti = ti1 + T->width;
#ifdef HAVE_SSE2
  unsigned long incw = 0;
  if (T->width & 1) incw = 1;
  ti += incw;
#endif

  Le[0] = 0;
  Lm[0] = 0;
  for (unsigned int i = 1; i < twokay; ++i) {		// FIXME: make int again
    rci_t rowneeded = r + codebook[k]->inc[i - 1];
    m = M->rows[rowneeded] + blockoffset;

    /* Duff's device loop unrolling */
    wi_t n = count;
    switch (entry_point) {
    case 0: do { *(ti++) = *(m++) ^ *(ti1++);
    case 7:      *(ti++) = *(m++) ^ *(ti1++);
    case 6:      *(ti++) = *(m++) ^ *(ti1++);
    case 5:      *(ti++) = *(m++) ^ *(ti1++);
    case 4:      *(ti++) = *(m++) ^ *(ti1++);
    case 3:      *(ti++) = *(m++) ^ *(ti1++);
    case 2:      *(ti++) = *(m++) ^ *(ti1++);
    case 1:      *(ti++) = *(m++) ^ *(ti1++);
      } while (--n > 0);
    }
#ifdef HAVE_SSE2
    ti += incw; ti1 += incw;
#endif
    ti += blockoffset;
    ti1 += blockoffset;

    /* U is a basis but not the canonical basis, so we need to read what
       element we just created from T */
    Le[mzd_read_bits_int(T,i,c,k)] = i;
    Lm[codebook[k]->ord[i]] = i;
    
  }
  /* We need fix the table to update the transformation matrix
     correctly; e.g. if the first row has [1 0 1] and we clear a row
     below with [1 0 1] we need to encode that this row is cleared by
     adding the first row only ([1 0 0]). */
  for(unsigned int i = 1; i < twokay; ++i) {	// FIXME: make int again
    word const correction = CONVERT_TO_WORD(codebook[k]->ord[i]);
    mzd_xor_bits(T, i,c, k, correction);
  }
}

void mzd_process_rows2_pls(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t *T0, rci_t *E0, mzd_t *T1, rci_t *E1) {
  int const ka = k / 2;
  int const kb = k - k / 2;
  wi_t const blocknuma = startcol / RADIX;
  wi_t const blocknumb = (startcol + ka) / RADIX;
  wi_t const blockoffset = blocknumb - blocknuma;
  wi_t wide = M->width - blocknuma;

  if(wide < 3) {
    mzd_process_rows(M, startrow, stoprow, startcol, ka, T0, E0);
    mzd_process_rows(M, startrow, stoprow, startcol + ka, kb, T1, E1);
    return;
  }

  wide -= 2U;
#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = E0[ mzd_read_bits_int(M, r, startcol, ka) ];
    wordPtr t0 = T0->rows[x0] + blocknuma;
    wordPtr m0 = M->rows[r+0] + blocknuma;
    m0[0U] ^= t0[0U];
    m0[1U] ^= t0[1U];
    rci_t const x1 = E1[ mzd_read_bits_int(M, r, startcol+ka, kb) ];
    wordPtr t1 = T1->rows[x1] + blocknumb;
    for(wi_t i = blockoffset; i < 2; ++i) {
      m0[i] ^= t1[i - blockoffset];
    }

    t0 += 2U;
    t1 += 2U - blockoffset;
    m0 += 2U;

    _mzd_combine2(m0,t0,t1,wide);
  }
}

void mzd_process_rows3_pls(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t *T0, rci_t *E0, mzd_t *T1, rci_t *E1, mzd_t *T2, rci_t *E2) {
  int const rem = k % 3;
  int const ka = k / 3 + ((rem >= 2) ? 1 : 0);
  int const kb = k / 3 + ((rem >= 1) ? 1 : 0);
  int const kc = k / 3;
  wi_t const blocknuma = startcol / RADIX;
  wi_t const blocknumb = (startcol + ka) / RADIX;
  wi_t const blocknumc = (startcol + ka + kb) / RADIX;
  wi_t const blockoffsetb = blocknumb - blocknuma;
  wi_t const blockoffsetc = blocknumc - blocknuma;
  wi_t wide = M->width - blocknuma;

  if(wide < 4) {
    mzd_process_rows(M, startrow, stoprow, startcol, ka, T0, E0);
    mzd_process_rows(M, startrow, stoprow, startcol + ka, kb, T1, E1);
    mzd_process_rows(M, startrow, stoprow, startcol + ka + kb, kc, T2, E2);
    return;
  }

  wide -= 3U;
#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = E0[ mzd_read_bits_int(M, r, startcol, ka) ];
    wordPtr t0 = T0->rows[x0] + blocknuma;
    wordPtr m0 = M->rows[r] + blocknuma;
    m0[0U] ^= t0[0U];
    m0[1U] ^= t0[1U];
    m0[2U] ^= t0[2U];

    t0 += 3U;

    rci_t const x1 = E1[ mzd_read_bits_int(M, r, startcol+ka, kb) ];
    wordPtr t1 = T1->rows[x1] + blocknumb;
    for(wi_t i = blockoffsetb; i < 3; ++i) {
      m0[i] ^= t1[i-blockoffsetb];
    }
    t1 += 3U - blockoffsetb;

    rci_t const x2 = E2[ mzd_read_bits_int(M, r, startcol+ka+kb, kc) ];
    wordPtr t2 = T2->rows[x2] + blocknumc;
    for(wi_t i = blockoffsetc; i < 3; ++i) {
      m0[i] ^= t2[i-blockoffsetc];
    }
    t2 += 3U - blockoffsetc;

    m0 += 3U;

    _mzd_combine3(m0,t0,t1,t2,wide);
  }
}

void mzd_process_rows4_pls(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, int k, mzd_t *T0, rci_t *E0, mzd_t *T1, rci_t *E1, mzd_t *T2, rci_t *E2, mzd_t *T3, rci_t *E3) {
  int const rem = k % 4;
  int const ka = k / 4 + ((rem >= 3) ? 1 : 0);
  int const kb = k / 4 + ((rem >= 2) ? 1 : 0);
  int const kc = k / 4 + ((rem >= 1) ? 1 : 0);
  int const kd = k / 4;
  wi_t const blocknuma = startcol / RADIX;
  wi_t const blocknumb = (startcol + ka) / RADIX;
  wi_t const blocknumc = (startcol + ka + kb) / RADIX;
  wi_t const blocknumd = (startcol + ka + kb + kc) / RADIX;
  wi_t const blockoffsetb = blocknumb - blocknuma;
  wi_t const blockoffsetc = blocknumc - blocknuma;
  wi_t const blockoffsetd = blocknumd - blocknuma;
  wi_t wide = M->width - blocknuma;

  if(wide < 5) {
    mzd_process_rows(M, startrow, stoprow, startcol, ka, T0, E0);
    mzd_process_rows(M, startrow, stoprow, startcol + ka, kb, T1, E1);
    mzd_process_rows(M, startrow, stoprow, startcol + ka + kb, kc, T2, E2);
    mzd_process_rows(M, startrow, stoprow, startcol + ka + kb + kc, kd, T3, E3);
    return;
  }
  wide -= 4U;
#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(rci_t r = startrow; r < stoprow; ++r) {
    rci_t const x0 = E0[mzd_read_bits_int(M, r, startcol, ka)];
    wordPtr t0 = T0->rows[x0] + blocknuma;
    wordPtr m0 = M->rows[r] + blocknuma;
    m0[0U] ^= t0[0U];
    m0[1U] ^= t0[1U];
    m0[2U] ^= t0[2U];
    m0[3U] ^= t0[3U];

    t0 += 4U;

    rci_t const x1 = E1[ mzd_read_bits_int(M, r, startcol+ka, kb) ];
    wordPtr t1 = T1->rows[x1] + blocknumb;
    for(wi_t i = blockoffsetb; i < 4; ++i) {
      m0[i] ^= t1[i - blockoffsetb];
    }
    t1 += 4U - blockoffsetb;

    rci_t const x2 = E2[ mzd_read_bits_int(M, r, startcol+ka+kb, kc) ];
    wordPtr t2 = T2->rows[x2] + blocknumc;
    for(wi_t i = blockoffsetc; i < 4; ++i) {
      m0[i] ^= t2[i - blockoffsetc];
    }
    t2 += 4U - blockoffsetc;

    rci_t const x3 = E3[ mzd_read_bits_int(M, r, startcol+ka+kb+kc, kd) ];
    wordPtr t3 = T3->rows[x3] + blocknumd;
    for(wi_t i = blockoffsetd; i < 4; ++i) {
      m0[i] ^= t3[i - blockoffsetd];
    }
    t3 += 4U - blockoffsetd;

    m0 += 4U;

    _mzd_combine4(m0, t0, t1, t2, t3, wide);
  }
}

void _mzd_finish_pls_done_pivots(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const start_col, wi_t const addblock, int const k) {
  /* perform needed row swaps */
  for(rci_t i = start_row; i < start_row + k; ++i) {
    _mzd_row_swap(A, i, P->values[i.val()], addblock);
  }

  for(int i = 1; i < k; ++i) {
    word const tmp = mzd_read_bits(A, start_row + i, start_col, i);
    wordPtr target = A->rows[start_row + i];
    for(int j = 0; j < i; ++j) {
      if((tmp & ONE << j)) {
        wordConstPtr source = A->rows[start_row + j];
        for(wi_t w = addblock; w < A->width; ++w) {
          target[w] ^= source[w];
        }
      }
    }
  }
}

void _mzd_finish_pls_done_rest1(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                                int k0, mzd_t *T0, rci_t const *M0) {

  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;
  
  for(rci_t i = start_row + k0; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col,k0)];
    wordConstPtr s0 = T0->rows[x0] + addblock;
    wordPtr t = A->rows[i] + addblock;
    _mzd_combine(t, s0, wide);
  }
}


void _mzd_finish_pls_done_rest2(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                                int k0, mzd_t *T0, rci_t const *M0,
                                int k1, mzd_t *T1, rci_t const *M1) {

  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;
  
  for(rci_t i = start_row + k0 + k1; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col,k0)];
    rci_t x1 = M1[mzd_read_bits_int(A,i,start_col+k0,k1)];
    wordConstPtr s0 = T0->rows[x0] + addblock;
    wordConstPtr s1 = T1->rows[x1] + addblock;
    wordPtr t = A->rows[i] + addblock;
    _mzd_combine2(t, s0, s1, wide);
  }
}


void _mzd_finish_pls_done_rest3(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                                int k0, mzd_t *T0, rci_t const *M0,
                                int k1, mzd_t *T1, rci_t const *M1,
                                int k2, mzd_t *T2, rci_t const *M2) {

  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;

  for(rci_t i = start_row + k0 + k1 + k2; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col, k0)];
    rci_t x1 = M1[mzd_read_bits_int(A,i,start_col+k0, k1)];
    rci_t x2 = M2[mzd_read_bits_int(A,i,start_col+k0+k1, k2)];
    wordConstPtr s0 = T0->rows[x0] + addblock;
    wordConstPtr s1 = T1->rows[x1] + addblock;
    wordConstPtr s2 = T2->rows[x2] + addblock;
    wordPtr t = A->rows[i] + addblock;
    _mzd_combine3(t, s0, s1, s2, wide);
  }
}


void _mzd_finish_pls_done_rest4(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                                int k0, mzd_t *T0, rci_t const *M0,
                                int k1, mzd_t *T1, rci_t const *M1,
                                int k2, mzd_t *T2, rci_t const *M2,
                                int k3, mzd_t *T3, rci_t const *M3) {

  wi_t const wide = A->width - addblock;
  if(wide <= 0)
    return;

  for(rci_t i = start_row + k0 + k1 + k2 + k3; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col, k0)];
    rci_t x1 = M1[mzd_read_bits_int(A,i,start_col+k0, k1)];
    rci_t x2 = M2[mzd_read_bits_int(A,i,start_col+k0+k1, k2)];
    rci_t x3 = M3[mzd_read_bits_int(A,i,start_col+k0+k1+k2, k3)];
    wordConstPtr s0 = T0->rows[x0] + addblock;
    wordConstPtr s1 = T1->rows[x1] + addblock;
    wordConstPtr s2 = T2->rows[x2] + addblock;
    wordConstPtr s3 = T3->rows[x3] + addblock;
    wordPtr t = A->rows[i] + addblock;
    _mzd_combine4(t, s0, s1, s2, s3, wide);
  }
}



/* extract U from A for table creation */
mzd_t *_mzd_pls_to_u(mzd_t *U, mzd_t *A, rci_t r, rci_t c, int k) {
  /* this function call is now rather cheap, but it could be avoided
     completetly if needed */
  assert(U->offset == 0);
  assert(A->offset == 0);
  rci_t startcol = (c / RADIX) * RADIX;
  mzd_submatrix(U, A, r, 0, r+k, A->ncols);

  for(rci_t i = 0; i < k; ++i)
    for(rci_t j = startcol; j < c + i; ++j) 
      mzd_write_bit(U, i, j,  0);
  return U;
}

/* method of many people factorisation */
rci_t _mzd_pls_mmpf(mzd_t *A, mzp_t *P, mzp_t *Q, int k) {
  assert(A->offset == 0);
  rci_t const nrows = A->nrows;
  rci_t const ncols = A->ncols; 
  rci_t curr_row = 0;
  rci_t curr_col = 0;
  int kbar = 0;
  rci_t done_row = 0;

  if(k == 0) {
    k = m4ri_opt_k(nrows, ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.5 * TWOPOW(k) * A->ncols > CPU_L2_CACHE / 2.0)
      k -= 1;
  }

  int kk = 4 * k;

  for(rci_t i = 0; i < ncols; ++i) 
    Q->values[i.val()] = i;

  for(rci_t i = 0; i < A->nrows; ++i)
    P->values[i.val()] = i;

  mzd_t *T0 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T1 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T2 = mzd_init(TWOPOW(k), ncols);
  mzd_t *T3 = mzd_init(TWOPOW(k), ncols);
  mzd_t *U = mzd_init((unsigned int)kk, ncols);	// FIXME

  /* these are the elimination lookups */
  rci_t *E0 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *E1 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *E2 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *E3 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));

  /* these are the multiplication lookups */
  rci_t *M0 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *M1 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *M2 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));
  rci_t *M3 = (rci_t*)m4ri_mm_calloc(TWOPOW(k), sizeof(rci_t));

  rci_t *done = (rci_t*)m4ri_mm_malloc(kk * sizeof(rci_t));

  /**
   * The algorithm proceeds as follows 
   */

  while(curr_col < ncols && curr_row < nrows) {
    if(curr_col + kk > ncols)
      kk = (ncols - curr_col).val();

    /**
     * 1. compute PLS factorisation for the kbar x kbar submatrix A00
\verbatim
       RADIX * splitblock
--------------------------------------
| A00  |  A10                        |
|      |                             |
-------------------------------------- kbar
| A01  |  A11                        |
|      |                             | 
-------------------------------------- done_row
| A02  | A21                         |
|      |                             |
|      |                             |
|      |                             |
|      |                             |
|      |                             |
--------------------------------------
\endverbatim
     */
    wi_t splitblock = (curr_col + kk) / RADIX + 1U;

    kbar = _mzd_pls_submatrix(A, curr_row, nrows, curr_col, kk, P, Q, done, &done_row, splitblock);

    /**
     * 2. update A10 
     */

    _mzd_finish_pls_done_pivots(A, P, curr_row, curr_col, splitblock, kbar);

    /**
     * 3. extract U from A0 = (A00 | A10) 
     */

    _mzd_pls_to_u(U, A, curr_row, curr_col, kbar);

    if(kbar > 3 * k) {
      int const rem = kbar % 4;
  
      int const ka = kbar / 4 + ((rem >= 3) ? 1 : 0);
      int const kb = kbar / 4 + ((rem >= 2) ? 1 : 0);
      int const kc = kbar / 4 + ((rem >= 1) ? 1 : 0);
      int const kd = kbar / 4;

      rci_t const first_col = 0;
      mzd_make_table_pls(U, first_col,          curr_col,          ka, T0, E0, M0);
      mzd_make_table_pls(U, first_col+ka,       curr_col+ka,       kb, T1, E1, M1);
      mzd_make_table_pls(U, first_col+ka+kb,    curr_col+ka+kb,    kc, T2, E2, M2);
      mzd_make_table_pls(U, first_col+ka+kb+kc, curr_col+ka+kb+kc, kd, T3, E3, M3);

      _mzd_finish_pls_done_rest4(A, P, curr_row, done_row+1, curr_col, splitblock,
                                 ka, T0, M0,
                                 kb, T1, M1,
                                 kc, T2, M2,
                                 kd, T3, M3);

      if (kbar == kk) {
        mzd_process_rows4_pls(A, done_row + 1, nrows, curr_col, kbar, T0, E0, T1, E1, T2, E2, T3, E3);
      } else {
        curr_col += 1; 
      }

    } else if(kbar > 2 * k) {
      int const rem = kbar % 3;

      int const ka = kbar / 3 + ((rem >= 2) ? 1 : 0);
      int const kb = kbar / 3 + ((rem >= 1) ? 1 : 0);
      int const kc = kbar / 3;

      rci_t const first_col = 0;
      mzd_make_table_pls(U, first_col,       curr_col,       ka, T0, E0, M0);
      mzd_make_table_pls(U, first_col+ka,    curr_col+ka,    kb, T1, E1, M1);
      mzd_make_table_pls(U, first_col+ka+kb, curr_col+ka+kb, kc, T2, E2, M2);

      _mzd_finish_pls_done_rest3(A, P, curr_row, done_row+1, curr_col, splitblock,
                                 ka, T0, M0,
                                 kb, T1, M1,
                                 kc, T2, M2);


      if (kbar == kk) {
        mzd_process_rows3_pls(A, done_row + 1, nrows, curr_col, kbar, T0, E0, T1, E1, T2, E2);
      } else {
        curr_col += 1; 
      }

    } else if(kbar > k) {
      int const ka = kbar / 2;
      int const kb = kbar - ka;

      rci_t const first_col = 0;
      mzd_make_table_pls(U, first_col,    curr_col,    ka, T0, E0, M0);
      mzd_make_table_pls(U, first_col+ka, curr_col+ka, kb, T1, E1, M1);

      _mzd_finish_pls_done_rest2(A, P, curr_row, done_row+1, curr_col, splitblock,
                                 ka, T0, M0,
                                 kb, T1, M1);

      if(kbar == kk) {
        mzd_process_rows2_pls(A, done_row + 1, nrows, curr_col, kbar, T0, E0, T1, E1);
      } else {
        curr_col += 1; 
      }

    } else if(kbar > 0) {

      /**
       * 4. generate multiplication and inversion tables T amd Tm from A0 
       */
      rci_t const first_col = 0;
      mzd_make_table_pls(U, first_col, curr_col, kbar, T0, E0, M0);

      /**
       * 5. update A11 using A10 and the multiplication table M
       */
      _mzd_finish_pls_done_rest1(A, P, curr_row, done_row+1, curr_col, splitblock, kbar, T0, M0);


      if(done_row < nrows) {
        /**
         * 6. update A2 = (A20 | A21) using the elimination table E
         */        
        mzd_process_rows(A, done_row + 1, nrows, curr_col, kbar, T0, E0);
      } else {
        curr_col += 1; 
      }

    } else {
      
      curr_col += 1;
      rci_t i = curr_row;
      rci_t j = curr_col;
      int found = mzd_find_pivot(A, curr_row, curr_col, &i, &j);
      if(found) {
        P->values[curr_row.val()] = i;
        Q->values[curr_row.val()] = j;
        mzd_row_swap(A, curr_row, i);
        wi_t const wrd = j / RADIX;
        word const bm = ONE << (j % RADIX);
	if (j + 1 < A->ncols)
	  for(rci_t l = curr_row + 1; l < nrows; ++l)
	    if(A->rows[l][wrd] & bm)
	      mzd_row_add_offset(A, l, curr_row, j + 1);
        curr_col = j + 1;
        ++curr_row;
      } else {
        break;
      }
    }
    curr_col += kbar;
    curr_row += kbar;
    if (kbar > 0)
      if (kbar == kk && kk < 4 * k)
        kk = kbar + 1;
      else
        kk = kbar;
    else if(kk > 2)
      kk = kk / 2;
  }

  /* Now compressing L */
  for (rci_t j = 0; j < curr_row; ++j){
    if (Q->values[j.val()]>j) {
      mzd_col_swap_in_rows(A, Q->values[j.val()], j, j, curr_row);
    }
  }
  mzp_t *Qbar = mzp_init_window(Q, 0, curr_row);
  mzd_apply_p_right_trans_even_capped(A, Qbar, curr_row, 0);
  mzp_free_window(Qbar);

  mzd_free(U);
  mzd_free(T0);
  mzd_free(T1);
  mzd_free(T2);
  mzd_free(T3);
  m4ri_mm_free(E0);  m4ri_mm_free(M0);
  m4ri_mm_free(E1);  m4ri_mm_free(M1);
  m4ri_mm_free(E2);  m4ri_mm_free(M2);
  m4ri_mm_free(E3);  m4ri_mm_free(M3);
  m4ri_mm_free(done);
  return curr_row;
}

rci_t _mzd_pluq_mmpf(mzd_t *A, mzp_t *P, mzp_t *Q, int const k) {
  rci_t r  = _mzd_pls_mmpf(A, P, Q, k);
  mzd_apply_p_right_trans_tri(A, Q);
  return r;
}
