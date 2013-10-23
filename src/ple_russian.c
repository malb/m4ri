/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008-2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include "m4ri_config.h"

#if __M4RI_HAVE_SSE2
#include <emmintrin.h>
#endif

#include "ple_russian.h"
#include "brilliantrussian.h"
#include "graycode.h"
#include "xor.h"
#ifndef NDEBUG
#include "mmc.h"
#endif


/** the number of tables used in PLE decomposition **/
#define __M4RI_PLE_NTABLES 8

ple_table_t *ple_table_init(int k, rci_t ncols) {
  ple_table_t *T = (ple_table_t*)m4ri_mm_malloc(sizeof(ple_table_t));
  T->T =  mzd_init(__M4RI_TWOPOW(k), ncols);
  T->M = (rci_t*)m4ri_mm_malloc(__M4RI_TWOPOW(k)*sizeof(rci_t));
  T->E = (rci_t*)m4ri_mm_malloc(__M4RI_TWOPOW(k)*sizeof(rci_t));
  T->B =  (word*)m4ri_mm_malloc(__M4RI_TWOPOW(k)*sizeof(word) );
  return T;
}

void ple_table_free(ple_table_t *T) {
  m4ri_mm_free(T->B);
  m4ri_mm_free(T->M);
  m4ri_mm_free(T->E);
  mzd_free(T->T);
  m4ri_mm_free(T);
}

static inline rci_t _max_value(rci_t *data, int length) {
  rci_t max = 0;
  for(int i = 0; i < length; ++i) {
    max = MAX(max, data[i]);
  }
  return max;
}

static inline void _kk_setup(int const kk, int const knar, int *k_, int *knar_, int const *pivots, int const ntables) {
  assert((ntables <= __M4RI_PLE_NTABLES) & (ntables > 0));

  int lb[__M4RI_PLE_NTABLES], ub[__M4RI_PLE_NTABLES];
  int rem = kk % ntables;

  int s = 1;

  switch(ntables) {
  case 8: k_[ntables - 8] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[7] = 0; assert((k_[ntables - 8] > 0));
  case 7: k_[ntables - 7] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[6] = 0; assert((k_[ntables - 7] > 0));
  case 6: k_[ntables - 6] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[5] = 0; assert((k_[ntables - 6] > 0));
  case 5: k_[ntables - 5] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[4] = 0; assert((k_[ntables - 5] > 0));
  case 4: k_[ntables - 4] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[3] = 0; assert((k_[ntables - 4] > 0));
  case 3: k_[ntables - 3] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[2] = 0; assert((k_[ntables - 3] > 0));
  case 2: k_[ntables - 2] = kk / ntables + ((rem >= (ntables - s++)) ? 1 : 0); knar_[1] = 0; assert((k_[ntables - 2] > 0));
  case 1: k_[ntables - 1] = kk / ntables; knar_[0] = 0;  assert((k_[ntables - 1] > 0));
  }

  lb[0] = 0;
  ub[0] = k_[0];
  for(int i = 1; i < ntables; i++) {
    lb[i] = lb[i-1] + k_[i-1];
    ub[i] = ub[i-1] + k_[i];
  }

  for(int i=0; i<knar; i++) {
    for(int j=0;j<ntables;j++)
      if (pivots[i] >= lb[j] && pivots[i] < ub[j]) {
        knar_[j]++;
    }
  }
}

int _mzd_ple_submatrix(mzd_t *A,
                       rci_t const start_row, rci_t const stop_row,
                       rci_t const start_col, int const k,
                       mzp_t *P, mzp_t *Q, rci_t *pivots,
                       rci_t *done, rci_t *done_row, wi_t const splitblock) {
  word bm[__M4RI_PLE_NTABLES * __M4RI_MAXKAY];
  wi_t os[__M4RI_PLE_NTABLES * __M4RI_MAXKAY];

  /* we're essentially constructing a submatrix but cheaply */
  wi_t const width = A->width;
  rci_t const ncols = A->ncols;
  int const flags = A->flags;
  word high_bitmask = A->high_bitmask;

  if (A->width > splitblock) {
    A->width = splitblock;
    A->ncols = splitblock * m4ri_radix;
    A->flags &= mzd_flag_multiple_blocks;
    A->flags |= (mzd_flag_windowed_zerooffset | mzd_flag_windowed_zeroexcess);
    A->high_bitmask = m4ri_ffff;
    /* No need to set mzd_flag_windowed_ownsblocks, because we won't free A until it's elements are restored below. */
  }

  int curr_pos;
  int rank = 0;
  for(curr_pos = 0; curr_pos < k; ++curr_pos) {
    os[curr_pos] = (start_col + curr_pos) / m4ri_radix;
    bm[curr_pos] = m4ri_one << ((start_col + curr_pos) % m4ri_radix);
    int found = 0;
    /* search for some pivot */
    rci_t i;
    for(i = start_row + rank; i < stop_row; ++i) {
      word const tmp = mzd_read_bits(A, i, start_col, curr_pos + 1);
      if(tmp) {
        word *Arow = A->rows[i];
        /* clear before but preserve transformation matrix */
        for (rci_t l = 0; l < rank; ++l)
          if(done[l] < i) {
            if((Arow[os[pivots[l]]] & bm[pivots[l]]))
              mzd_row_add_offset(A, i, start_row + l, start_col + pivots[l] + 1);
            done[l] = i; /* encode up to which row we added for l already */
          }
        if(mzd_read_bit(A, i, start_col + curr_pos)) {
          found = 1;
          break;
        }
      }
    }
    if (found) {
      P->values[start_row + rank] = i;
      mzd_row_swap(A, i, start_row + rank);

      Q->values[start_row + rank] = start_col + curr_pos;
      pivots[rank] = curr_pos;
      done[rank] = i;
      rank++;
    }
  }

  /* finish submatrix */
  *done_row = _max_value(done, rank);
  for(rci_t c2 = 0; c2 < rank && start_col + pivots[c2] < A->ncols -1; ++c2)
    for(rci_t r2 = done[c2] + 1; r2 <= *done_row; ++r2)
      if(mzd_read_bit(A, r2, start_col + pivots[c2]))
        mzd_row_add_offset(A, r2, start_row + c2, start_col + pivots[c2] + 1);

  /* reset to original size */
  A->ncols = ncols;
  A->width = width;
  A->flags = flags;
  A->high_bitmask = high_bitmask;

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
  __M4RI_DD_MZP(Q);
  __M4RI_DD_INT(curr_pos);
  return rank;
}

/* create a table of all 2^k linear combinations */
void mzd_make_table_ple(mzd_t const *A, rci_t r, rci_t writecol, int k, int knar, ple_table_t *table, rci_t *offsets, int base, rci_t readcol) {

  mzd_t *T = table->T;
  rci_t *E = table->E;
  rci_t *M = table->M;
  word  *B = table->B;

  // Note that this restricts the number of columns of any matrix to __M4RI_MAX_MZD_BLOCKSIZE *
  // radix / twokay = 268 million.
  assert(!(T->flags & mzd_flag_multiple_blocks));

  wi_t const writeblock  = writecol / m4ri_radix;
  wi_t const readblock   = readcol  / m4ri_radix;

  assert(writeblock - readblock <= 1);

  int const twokay = __M4RI_TWOPOW(knar);
  wi_t const wide  = T->width - writeblock;
  wi_t const count = (wide + 7) / 8;
  int const entry_point = wide % 8;
  wi_t const next_row_offset = writeblock + T->rowstride - T->width;

  word *a;
  word *ti1 = T->rows[0] + writeblock;
  word *ti = ti1         + T->rowstride;

  M[0] = 0; E[0] = 0; B[0] = 0;
  for (int i = 1; i < twokay; ++i) {
    T->rows[i][readblock] = 0; /* we make sure that we can safely add from readblock */

    rci_t rowneeded = r + m4ri_codebook[knar]->inc[i - 1];
    a = A->rows[rowneeded] + writeblock;

    /* Duff's device loop unrolling */
    wi_t n = count;
    switch (entry_point) {
    case 0: do { *(ti++) = *(a++) ^ *(ti1++);
    case 7:      *(ti++) = *(a++) ^ *(ti1++);
    case 6:      *(ti++) = *(a++) ^ *(ti1++);
    case 5:      *(ti++) = *(a++) ^ *(ti1++);
    case 4:      *(ti++) = *(a++) ^ *(ti1++);
    case 3:      *(ti++) = *(a++) ^ *(ti1++);
    case 2:      *(ti++) = *(a++) ^ *(ti1++);
    case 1:      *(ti++) = *(a++) ^ *(ti1++);
      } while (--n > 0);
    }
    ti  += next_row_offset;
    ti1 += next_row_offset;

    /* U is a basis but not the canonical basis, so we need to read what element we just created from T */
    E[ mzd_read_bits_int(T, i, writecol, k) ] = i;
    M[ m4ri_spread_bits(m4ri_codebook[k]->ord[i], offsets, knar, base) ] = i;
  }

  /* We need fix the table to update the transformation matrix correctly; e.g. if the first row has
     [1 0 1] and we clear a row below with [1 0 1] we need to encode that this row is cleared by
     adding the first row only ([1 0 0]). */

  int const bits_to_read = MIN(m4ri_radix, T->ncols - readcol);
  for(int i = 1; i < twokay; ++i) {
    word const fix = m4ri_spread_bits(__M4RI_CONVERT_TO_WORD(m4ri_codebook[k]->ord[i]), offsets, knar, base);
    mzd_xor_bits(T, i, writecol, k, fix);
    B[i] = mzd_read_bits(T, i, readcol, bits_to_read);
  }

  __M4RI_DD_MZD(T);
  __M4RI_DD_RCI_ARRAY(E, twokay);
  __M4RI_DD_RCI_ARRAY(M, twokay);
}

#define N 2
#include "ple_russian_template.h"
#undef N

#define N 3
#include "ple_russian_template.h"
#undef N

#define N 4
#include "ple_russian_template.h"
#undef N

#define N 5
#include "ple_russian_template.h"
#undef N

#define N 6
#include "ple_russian_template.h"
#undef N

#define N 7
#include "ple_russian_template.h"
#undef N

#define N 8
#include "ple_russian_template.h"
#undef N

void _mzd_ple_a10(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const start_col,
                  wi_t const addblock, int const k, rci_t *pivots) {
  if(addblock == A->width)
    return;

  /* perform needed row swaps */
  for(rci_t i = start_row; i < start_row + k; ++i) {
    _mzd_row_swap(A, i, P->values[i], addblock);
  }

  for(int i = 1; i < k; ++i) {
    word const tmp = mzd_read_bits(A, start_row + i, start_col, pivots[i]);
    word *target = A->rows[start_row + i];
    for(int j = 0; j < i; ++j) {
      if((tmp & m4ri_one << pivots[j])) {
        word const *source = A->rows[start_row + j];
        for(wi_t w = addblock; w < A->width; ++w) {
          target[w] ^= source[w];
        }
      }
    }
  }

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
}

void _mzd_ple_a11_1(mzd_t *A,
                    rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock,
                    int const k, ple_table_t const *T0) {

  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;

  for(rci_t i = start_row; i < stop_row; ++i) {
    rci_t x0 = T0->M[mzd_read_bits_int(A,i,start_col, k)];
    word const *s0 = T0->T->rows[x0] + addblock;
    word *t = A->rows[i] + addblock;
    _mzd_combine(t, s0, wide);
  }

  __M4RI_DD_MZD(A);
}


/* extract E from A for table creation */
mzd_t *_mzd_ple_to_e(mzd_t *E, mzd_t const *A, rci_t r, rci_t c, int k, rci_t *offsets) {
  /* this function call is now rather cheap, but it could be avoided
     completetly if needed */
  rci_t startcol = (c / m4ri_radix) * m4ri_radix;
  mzd_submatrix(E, A, r, 0, r+k, A->ncols);

  for(rci_t i = 0; i < k; ++i) {
    for(rci_t j = startcol; j < c + offsets[i]; j+=m4ri_radix)
      mzd_clear_bits(E, i, j, MIN(c + offsets[i] - j, m4ri_radix));
  }

  __M4RI_DD_MZD(E);
  return E;
}

/* method of many people factorisation */
rci_t _mzd_ple_russian(mzd_t *A, mzp_t *P, mzp_t *Q, int k) {
  rci_t const nrows = A->nrows;
  rci_t const ncols = A->ncols;
  rci_t curr_row = 0;
  rci_t curr_col = 0;
  rci_t done_row = 0;

  int knar = 0;

  /** compute good k **/

  if(k == 0) {
    /* __M4RI_CPU_L2_CACHE == __M4RI_PLE_NTABLES * 2^k * B->width * 8 */
    k = (int)log2((__M4RI_CPU_L2_CACHE/8)/(double)A->width/(double)__M4RI_PLE_NTABLES);

    rci_t const klog = round(0.75 * log2_floor(MIN(nrows, ncols)));

    if(klog < k) k = klog;
    if (k<2)     k = 2;
    else if(k>8) k = 8;
  }
  int kk = __M4RI_PLE_NTABLES * k;
  assert(kk <= m4ri_radix);

  /** initialise permutations as identity **/

  for(rci_t i = 0; i < ncols; ++i)
    Q->values[i] = i;

  for(rci_t i = 0; i < nrows; ++i)
    P->values[i] = i;

  ple_table_t *T[__M4RI_PLE_NTABLES];

  for(int i=0; i<__M4RI_PLE_NTABLES; i++)
    T[i] = ple_table_init(k, ncols);

  mzd_t *U = mzd_init(kk, ncols);

  /* these are the elimination lookups */

  rci_t *done   = (rci_t*)m4ri_mm_malloc(kk * sizeof(rci_t));
  rci_t *pivots = (rci_t*)m4ri_mm_malloc(kk * sizeof(rci_t));

  /**
   * The algorithm proceeds as follows
   */

  while(curr_col < ncols && curr_row < nrows) {
    if(curr_col + kk > ncols)
      kk = ncols - curr_col;

    /**
     * 1. compute PLE factorisation for the knar x knar submatrix A00
\verbatim
       m4ri_radix * splitblock
--------------------------------------
| A00  |  A10                        |
|      |                             |
-------------------------------------- knar
| A01  |  A11                        |
|      |                             |
-------------------------------------- done_row
| A02  |  A12                        |
|      |                             |
|      |                             |
|      |                             |
|      |                             |
|      |                             |
--------------------------------------
\endverbatim
     */
    /* + 8 blocks comes from 64 byte cache line, rest is minimum required and maximum possible */
    wi_t splitblock = MIN( MAX( (curr_col + kk) / m4ri_radix + 1, (curr_col / m4ri_radix) + 8 ), A->width );

    knar = _mzd_ple_submatrix(A, curr_row, nrows, curr_col, kk, P, Q, pivots, done, &done_row, splitblock);

    /**
     * 2. update A10
     */

    _mzd_ple_a10(A, P, curr_row, curr_col, splitblock, knar, pivots);

    /**
     * 3. extract U from A0 = (A00 | A10)
     */

    _mzd_ple_to_e(U, A, curr_row, curr_col, knar, pivots);


    // treat no-pivot-was-found case
    if (knar == 0) {
      curr_col += kk;
      curr_row += knar;

      rci_t i = curr_row;
      rci_t j = curr_col;
      int found = mzd_find_pivot(A, curr_row, curr_col, &i, &j);
      if(found) {
        P->values[curr_row] = i;
        Q->values[curr_row] = j;
        mzd_row_swap(A, curr_row, i);
        wi_t const wrd = j / m4ri_radix;
        word const bm = m4ri_one << (j % m4ri_radix);
        if (j + 1 < A->ncols)
          for(rci_t l = curr_row + 1; l < nrows; ++l)
            if(A->rows[l][wrd] & bm)
              mzd_row_add_offset(A, l, curr_row, j + 1);
        curr_col = j + 1;
        ++curr_row;
      } else {
        break;
      }
      continue;
    }

    int k_[__M4RI_PLE_NTABLES], knar_[__M4RI_PLE_NTABLES], ntables = 0;

    if (__M4RI_PLE_NTABLES >= 8 && kk >= 7*k && kk >= 8)
      ntables = 8;
    else if (__M4RI_PLE_NTABLES >= 7 && kk >= 6*k && kk >= 7)
      ntables = 7;
    else if (__M4RI_PLE_NTABLES >= 6 && kk >= 5*k && kk >= 6)
      ntables = 6;
    else if (__M4RI_PLE_NTABLES >= 5 && kk >= 4*k && kk >= 5)
      ntables = 5;
    else if (__M4RI_PLE_NTABLES >= 4 && kk >= 3*k && kk >= 4)
      ntables = 4;
    else if (__M4RI_PLE_NTABLES >= 3 && kk >= 2*k && kk >= 3)
      ntables = 3;
    else if (__M4RI_PLE_NTABLES >= 2 && kk >=   k && kk >= 2)
      ntables = 2;
    else
      ntables = 1;

    _kk_setup(kk, knar, k_, knar_, pivots, ntables);

    /**
     * 4. generate multiplication and inversion tables T amd E from U
     */

    rci_t i_knar = 0;
    rci_t i_curr_col = curr_col;
    rci_t *i_pivots = pivots;
    int i_base = 0;
    for(int i=0; i<ntables; i++) {
      mzd_make_table_ple(U, i_knar, i_curr_col, k_[i], knar_[i], T[i], i_pivots,  i_base, curr_col);
      i_knar += knar_[i];
      i_curr_col += k_[i];
      i_pivots += knar_[i];
      i_base += k_[i];
    }

    switch(ntables) {
    case 8:
      /**
       * 5. update A1 = (A01 | A11) */
      _mzd_ple_a11_8(A, curr_row + knar , done_row + 1, curr_col, splitblock, k_, (const ple_table_t**)T);
      /**
       * 6. update A2 = (A02 | A12) */
      if (done_row < nrows) _mzd_process_rows_ple_8(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 7:
      _mzd_ple_a11_7(A, curr_row + knar , done_row + 1, curr_col, splitblock, k_, (const ple_table_t**)T);

      if (done_row < nrows) _mzd_process_rows_ple_7(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 6:
      _mzd_ple_a11_6(A, curr_row + knar , done_row + 1, curr_col, splitblock, k_, (const ple_table_t**)T);
      if (done_row < nrows) _mzd_process_rows_ple_6(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 5:
      _mzd_ple_a11_5(A, curr_row + knar, done_row + 1, curr_col, splitblock, k_, (const ple_table_t**)T);

      if (done_row < nrows) _mzd_process_rows_ple_5(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 4:
      _mzd_ple_a11_4(A, curr_row + knar, done_row + 1, curr_col, splitblock, k_, (const ple_table_t**)T);

      if (done_row < nrows) _mzd_process_rows_ple_4(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 3:
      _mzd_ple_a11_3(A, curr_row + knar, done_row+1, curr_col, splitblock, k_, (const ple_table_t**)T);

      if (done_row < nrows) _mzd_process_rows_ple_3(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 2:
      _mzd_ple_a11_2(A, curr_row + knar, done_row+1, curr_col, splitblock, k_, (const ple_table_t**)T);

      if(done_row < nrows) _mzd_process_rows_ple_2(A, done_row + 1, nrows, curr_col, k_, (const ple_table_t**)T);
      break;

    case 1:
      _mzd_ple_a11_1(A, curr_row + knar, done_row+1, curr_col, splitblock, kk, T[0]);

      if(done_row < nrows) mzd_process_rows(A, done_row + 1, nrows, curr_col, kk, T[0]->T, T[0]->E);
      break;
    default:
      m4ri_die("ntables = %d not supported.\n",ntables);
    }

    curr_col += kk;
    curr_row += knar;
  }

  /* Now compressing L */
  for (rci_t j = 0; j < curr_row; ++j){
    if (Q->values[j] > j) {
      mzd_col_swap_in_rows(A, Q->values[j], j, j, curr_row);
    }
  }
  mzp_t *Qbar = mzp_init_window(Q, 0, curr_row);
  mzd_apply_p_right_trans_even_capped(A, Qbar, curr_row, 0);
  mzp_free_window(Qbar);

  mzd_free(U);
  for(int i=0; i<__M4RI_PLE_NTABLES; i++) {
    ple_table_free(T[i]);
  }
  m4ri_mm_free(done);   m4ri_mm_free(pivots);

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
  __M4RI_DD_MZP(Q);
  __M4RI_DD_RCI(curr_row);

  return curr_row;
}

rci_t _mzd_pluq_russian(mzd_t *A, mzp_t *P, mzp_t *Q, int const k) {
  rci_t r = _mzd_ple_russian(A, P, Q, k);
  mzd_apply_p_right_trans_tri(A, Q);
  return r;
}

