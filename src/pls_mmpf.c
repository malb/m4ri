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

#include "m4ri_config.h"
#include <assert.h>

#if __M4RI_HAVE_SSE2
#include <emmintrin.h>
#endif

#include "pls_mmpf.h"
#include "brilliantrussian.h"
#include "grayflex.h"
#include "xor.h"
#ifndef NDEBUG
#include "mmc.h"
#endif

/** the number of tables used in PLE decomposition **/
#define __M4RI_PLS_NTABLES 4

static inline rci_t _max_value(rci_t *data, int length) {
  rci_t max = 0;
  for(int i = 0; i < length; ++i) {
    max = MAX(max, data[i]);
  }
  return max;
}

static inline void _kk_setup(int const kk, int const knar, int *k_, int *knar_, int const *pivots, int const ntables) {
  int i,j, rem;
  int lb[4], ub[4];

  switch(ntables) {
  case 4:
    rem = kk % 4;
    k_[0] = kk / 4 + ((rem >= 3) ? 1 : 0);
    k_[1] = kk / 4 + ((rem >= 2) ? 1 : 0);
    k_[2] = kk / 4 + ((rem >= 1) ? 1 : 0);
    k_[3] = kk / 4;

    knar_[0] = 0;
    knar_[1] = 0;
    knar_[2] = 0;
    knar_[3] = 0;

    lb[0] = 0;
    lb[1] = lb[0]+k_[0];
    lb[2] = lb[1]+k_[1];
    lb[3] = lb[2]+k_[2];
    
    ub[0] =     0+k_[0];
    ub[1] = ub[0]+k_[1];
    ub[2] = ub[1]+k_[2];
    ub[3] = ub[2]+k_[3];

    assert((k_[0] > 0) && (k_[1] > 0) && (k_[2] > 0) && (k_[3] > 0));
    break;

  case 3:
    rem = kk % 3;
    k_[0] = kk / 3 + ((rem >= 2) ? 1 : 0);
    k_[1] = kk / 3 + ((rem >= 1) ? 1 : 0);
    k_[2] = kk / 3;

    knar_[0] = 0;
    knar_[1] = 0;
    knar_[2] = 0;

    lb[0] = 0;
    lb[1] = lb[0]+k_[0];
    lb[2] = lb[1]+k_[1];
    
    ub[0] =     0+k_[0];
    ub[1] = ub[0]+k_[1];
    ub[2] = ub[1]+k_[2];

    assert((k_[0] > 0) && (k_[1] > 0) && (k_[2] > 0));
    break;

  case 2: 
    k_[0] = kk / 2;
    k_[1] = kk - k_[0];

    knar_[0] = 0;
    knar_[1] = 0;

    lb[0] = 0;
    lb[1] = lb[0]+k_[0];
    
    ub[0] =     0+k_[0];
    ub[1] = ub[0]+k_[1];

    assert((k_[0] > 0) && (k_[1] > 0));
    break;

  case 1:
    k_[0] = kk;
    knar_[0] = knar;
    lb[0] = 0;
    ub[0] =     0+k_[0];

    break;

  default:
    break;
  }

  for(i=0; i<knar; i++) {
    for(j=0;j<ntables;j++)
      if (pivots[i] >= lb[j] && pivots[i] < ub[j]) {
        knar_[j]++;
    }
  }
}

int _mzd_pls_submatrix(mzd_t *A, 
                       rci_t const start_row, rci_t const stop_row, 
                       rci_t const start_col, int const k, 
                       mzp_t *P, mzp_t *Q, rci_t *pivots, 
                       rci_t *done, rci_t *done_row, wi_t const splitblock)  {
  word bm[__M4RI_PLS_NTABLES * __M4RI_MAXKAY];
  wi_t os[__M4RI_PLS_NTABLES * __M4RI_MAXKAY];

  /* we're essentially constructing a submatrix but cheaply */
  wi_t const width = A->width;
  rci_t const ncols = A->ncols;
  int const flags = A->flags;
  word low_bitmask = A->low_bitmask;
  word high_bitmask = A->high_bitmask;

  if (A->width > splitblock) {
    A->width = splitblock;
    A->ncols = splitblock * m4ri_radix;
    assert(A->offset == 0);
    A->flags &= mzd_flag_multiple_blocks;
    A->flags |= (mzd_flag_windowed_zerooffset | mzd_flag_windowed_zeroexcess);
    A->high_bitmask = A->low_bitmask = m4ri_ffff;
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
  A->low_bitmask = low_bitmask;
  A->high_bitmask = high_bitmask;

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
  __M4RI_DD_MZP(Q);
  __M4RI_DD_INT(curr_pos);
  return rank;
}

/* create a table of all 2^k linear combinations */
void mzd_make_table_pls(mzd_t const *M, rci_t r, rci_t c, int k, int knar, mzd_t *T, rci_t *Le, rci_t *Lm, rci_t *offsets, int base) {

  // Note that this restricts the number of columns of any matrix to
  // __M4RI_MAX_MZD_BLOCKSIZE * radix / twokay = 268 million.

  assert(!(T->flags & mzd_flag_multiple_blocks));
  wi_t const blockoffset= c / m4ri_radix;
  int const twokay= __M4RI_TWOPOW(knar);
  wi_t const wide = T->width - blockoffset;
  wi_t const count = (wide + 7) / 8;
  int const entry_point = wide % 8;
  wi_t const next_row_offset = blockoffset + T->rowstride - T->width;

  word *ti, *ti1, *m;

  ti1 = T->rows[0] + blockoffset;
  ti = ti1 + T->rowstride;

  Le[0] = 0;
  Lm[0] = 0;
  for (int i = 1; i < twokay; ++i) {
    rci_t rowneeded = r + m4ri_codebook[knar]->inc[i - 1];
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
    ti += next_row_offset;
    ti1 += next_row_offset;

    /* U is a basis but not the canonical basis, so we need to read what
       element we just created from T */
    Le[mzd_read_bits_int(T,i,c,k)] = i;
    Lm[m4ri_spread_bits(m4ri_codebook[k]->ord[i],offsets,knar,base)] = i;
    
  }
  /* We need fix the table to update the transformation matrix
     correctly; e.g. if the first row has [1 0 1] and we clear a row
     below with [1 0 1] we need to encode that this row is cleared by
     adding the first row only ([1 0 0]). */
  for(int i = 1; i < twokay; ++i) {
    word const correction = m4ri_spread_bits(__M4RI_CONVERT_TO_WORD(m4ri_codebook[k]->ord[i]), offsets, knar,base);
    mzd_xor_bits(T, i,c, k, correction);
  }

  __M4RI_DD_MZD(T);
  __M4RI_DD_RCI_ARRAY(Le, twokay);
  __M4RI_DD_RCI_ARRAY(Lm, twokay);
}

void mzd_process_rows2_pls(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, 
                           int const k0, mzd_t const *T0, rci_t const *E0, 
                           int const k1, mzd_t const *T1, rci_t const *E1) {

  wi_t const blocknuma = startcol / m4ri_radix;
  wi_t const blocknumb = (startcol + k0) / m4ri_radix;
  wi_t const blockoffset = blocknumb - blocknuma;
  wi_t wide = M->width - blocknuma;

  if(wide < 3) {
    mzd_process_rows(M, startrow, stoprow, startcol, k0, T0, E0);
    mzd_process_rows(M, startrow, stoprow, startcol + k0, k1, T1, E1);
    return;
  }

  rci_t r;

  wide -= 2;
#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r = startrow; r < stoprow; ++r) {
    rci_t const x0 = E0[ mzd_read_bits_int(M, r, startcol, k0) ];
    word const *t0 = T0->rows[x0] + blocknuma;
    word *m0 = M->rows[r+0] + blocknuma;
    m0[0] ^= t0[0];
    m0[1] ^= t0[1];
    rci_t const x1 = E1[ mzd_read_bits_int(M, r, startcol+k0, k1) ];
    word const *t1 = T1->rows[x1] + blocknumb;
    for(wi_t i = blockoffset; i < 2; ++i) {
      m0[i] ^= t1[i - blockoffset];
    }

    t0 += 2;
    t1 += 2 - blockoffset;
    m0 += 2;

    _mzd_combine2(m0,t0,t1,wide);
  }

  __M4RI_DD_MZD(M);
}

void mzd_process_rows3_pls(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, 
                           int const k0, mzd_t const *T0, rci_t const *E0, 
                           int const k1, mzd_t const *T1, rci_t const *E1,
			   int const k2, mzd_t const *T2, rci_t const *E2) {

  wi_t const blocknuma = startcol / m4ri_radix;
  wi_t const blocknumb = (startcol + k0) / m4ri_radix;
  wi_t const blocknumc = (startcol + k0 + k1) / m4ri_radix;
  wi_t const blockoffsetb = blocknumb - blocknuma;
  wi_t const blockoffsetc = blocknumc - blocknuma;
  wi_t wide = M->width - blocknuma;

  if(wide < 4) {
    mzd_process_rows(M, startrow, stoprow, startcol, k0, T0, E0);
    mzd_process_rows(M, startrow, stoprow, startcol + k0, k1, T1, E1);
    mzd_process_rows(M, startrow, stoprow, startcol + k0 + k1, k2, T2, E2);
    return;
  }

  rci_t r;

  wide -= 3;
#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r = startrow; r < stoprow; ++r) {
    rci_t const x0 = E0[ mzd_read_bits_int(M, r, startcol, k0) ];
    word const *t0 = T0->rows[x0] + blocknuma;
    word *m0 = M->rows[r] + blocknuma;
    m0[0] ^= t0[0];
    m0[1] ^= t0[1];
    m0[2] ^= t0[2];

    t0 += 3;

    rci_t const x1 = E1[ mzd_read_bits_int(M, r, startcol+k0, k1) ];
    word const *t1 = T1->rows[x1] + blocknumb;
    for(wi_t i = blockoffsetb; i < 3; ++i) {
      m0[i] ^= t1[i-blockoffsetb];
    }
    t1 += 3 - blockoffsetb;

    rci_t const x2 = E2[ mzd_read_bits_int(M, r, startcol+k0+k1, k2) ];
    word const *t2 = T2->rows[x2] + blocknumc;
    for(wi_t i = blockoffsetc; i < 3; ++i) {
      m0[i] ^= t2[i-blockoffsetc];
    }
    t2 += 3 - blockoffsetc;

    m0 += 3;

    _mzd_combine3(m0,t0,t1,t2,wide);
  }

  __M4RI_DD_MZD(M);
}

void mzd_process_rows4_pls(mzd_t *M, rci_t startrow, rci_t stoprow, rci_t startcol, 
                           int const k0, mzd_t const *T0, rci_t const *E0, 
                           int const k1, mzd_t const *T1, rci_t const *E1, 
                           int const k2, mzd_t const *T2, rci_t const *E2, 
                           int const k3, mzd_t const *T3, rci_t const *E3) {

  wi_t const blocknuma = startcol / m4ri_radix;
  wi_t const blocknumb = (startcol + k0) / m4ri_radix;
  wi_t const blocknumc = (startcol + k0 + k1) / m4ri_radix;
  wi_t const blocknumd = (startcol + k0 + k1 + k2) / m4ri_radix;
  wi_t const blockoffsetb = blocknumb - blocknuma;
  wi_t const blockoffsetc = blocknumc - blocknuma;
  wi_t const blockoffsetd = blocknumd - blocknuma;
  wi_t wide = M->width - blocknuma;

  if(wide < 5) {
    mzd_process_rows(M, startrow, stoprow, startcol,  k0, T0, E0);
    mzd_process_rows(M, startrow, stoprow, startcol + k0,  k1, T1, E1);
    mzd_process_rows(M, startrow, stoprow, startcol + k0 + k1,  k2, T2, E2);
    mzd_process_rows(M, startrow, stoprow, startcol + k0 + k1 + k2, k3, T3, E3);
    return;
  }
  wide -= 4;

  rci_t r;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r = startrow; r < stoprow; ++r) {
    rci_t const x0 = E0[mzd_read_bits_int(M, r, startcol, k0)];
    word *t0 = T0->rows[x0] + blocknuma;
    word *m0 = M->rows[r] + blocknuma;
    m0[0] ^= t0[0];
    m0[1] ^= t0[1];
    m0[2] ^= t0[2];
    m0[3] ^= t0[3];

    t0 += 4;

    rci_t const x1 = E1[ mzd_read_bits_int(M, r, startcol+k0, k1) ];
    word *t1 = T1->rows[x1] + blocknumb;
    for(wi_t i = blockoffsetb; i < 4; ++i) {
      m0[i] ^= t1[i - blockoffsetb];
    }
    t1 += 4 - blockoffsetb;

    rci_t const x2 = E2[ mzd_read_bits_int(M, r, startcol+k0+k1, k2) ];
    word *t2 = T2->rows[x2] + blocknumc;
    for(wi_t i = blockoffsetc; i < 4; ++i) {
      m0[i] ^= t2[i - blockoffsetc];
    }
    t2 += 4 - blockoffsetc;

    rci_t const x3 = E3[ mzd_read_bits_int(M, r, startcol+k0+k1+k2, k3) ];
    word *t3 = T3->rows[x3] + blocknumd;
    for(wi_t i = blockoffsetd; i < 4; ++i) {
      m0[i] ^= t3[i - blockoffsetd];
    }
    t3 += 4 - blockoffsetd;

    m0 += 4;

    _mzd_combine4(m0, t0, t1, t2, t3, wide);
  }

  __M4RI_DD_MZD(M);
}

void _mzd_pls_a10(mzd_t *A, mzp_t const *P, rci_t const start_row, rci_t const start_col, 
                  wi_t const addblock, int const k, rci_t *pivots) {
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

void _mzd_pls_a11_1(mzd_t *A, 
                    rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                    int const k, int const knar, mzd_t const *T0, rci_t const *M0) {

  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;
  
  for(rci_t i = start_row + knar; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col, k)];
    word const *s0 = T0->rows[x0] + addblock;
    word *t = A->rows[i] + addblock;
    _mzd_combine(t, s0, wide);
  }

  __M4RI_DD_MZD(A);
}


void _mzd_pls_a11_2(mzd_t *A,
                    rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                    int const k0, int const knar0, mzd_t const *T0, rci_t const *M0, 
                    int const k1, int const knar1, mzd_t const *T1, rci_t const *M1) {

  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;
  
  for(rci_t i = start_row + knar0 + knar1; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col,k0)];
    rci_t x1 = M1[mzd_read_bits_int(A,i,start_col+k0,k1)];
    word const *s0 = T0->rows[x0] + addblock;
    word const *s1 = T1->rows[x1] + addblock;
    word *t = A->rows[i] + addblock;
    _mzd_combine2(t, s0, s1, wide);
  }

  __M4RI_DD_MZD(A);
}


void _mzd_pls_a11_3(mzd_t *A,
                    rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                    int const k0, int const knar0, mzd_t const *T0, rci_t const *M0, 
                    int const k1, int const knar1, mzd_t const *T1, rci_t const *M1, 
                    int const k2, int const knar2, mzd_t const *T2, rci_t const *M2) { 
  wi_t const wide = A->width - addblock;
  if (wide <= 0)
    return;

  for(rci_t i = start_row + knar0 + knar1 + knar2; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col, k0)];
    rci_t x1 = M1[mzd_read_bits_int(A,i,start_col+k0, k1)];
    rci_t x2 = M2[mzd_read_bits_int(A,i,start_col+k0+k1, k2)];
    word const *s0 = T0->rows[x0] + addblock;
    word const *s1 = T1->rows[x1] + addblock;
    word const *s2 = T2->rows[x2] + addblock;
    word *t = A->rows[i] + addblock;
    _mzd_combine3(t, s0, s1, s2, wide);
  }

  __M4RI_DD_MZD(A);
}


void _mzd_pls_a11_4(mzd_t *A,
                    rci_t const start_row, rci_t const stop_row, rci_t const start_col, wi_t const addblock, 
                    int const k0, int const knar0, mzd_t const *T0, rci_t const *M0, 
                    int const k1, int const knar1, mzd_t const *T1, rci_t const *M1, 
                    int const k2, int const knar2, mzd_t const *T2, rci_t const *M2,
                    int const k3, int const knar3, mzd_t const *T3, rci_t const *M3) { 

  wi_t const wide = A->width - addblock;
  if(wide <= 0)
    return;

  for(rci_t i = start_row + knar0 + knar1 + knar2 + knar3; i < stop_row; ++i) {
    rci_t x0 = M0[mzd_read_bits_int(A,i,start_col, k0)];
    rci_t x1 = M1[mzd_read_bits_int(A,i,start_col+k0, k1)];
    rci_t x2 = M2[mzd_read_bits_int(A,i,start_col+k0+k1, k2)];
    rci_t x3 = M3[mzd_read_bits_int(A,i,start_col+k0+k1+k2, k3)];
    word const *s0 = T0->rows[x0] + addblock;
    word const *s1 = T1->rows[x1] + addblock;
    word const *s2 = T2->rows[x2] + addblock;
    word const *s3 = T3->rows[x3] + addblock;
    word *t = A->rows[i] + addblock;
    _mzd_combine4(t, s0, s1, s2, s3, wide);
  }

  __M4RI_DD_MZD(A);
}

/* extract U from A for table creation */
mzd_t *_mzd_pls_to_u(mzd_t *U, mzd_t const *A, rci_t r, rci_t c, int k, rci_t *offsets) {
  /* this function call is now rather cheap, but it could be avoided
     completetly if needed */
  assert(U->offset == 0);
  assert(A->offset == 0);
  rci_t startcol = (c / m4ri_radix) * m4ri_radix;
  mzd_submatrix(U, A, r, 0, r+k, A->ncols);

  for(rci_t i = 0; i < k; ++i)
    for(rci_t j = startcol; j < c + offsets[i]; ++j) 
      mzd_write_bit(U, i, j,  0);

  __M4RI_DD_MZD(U);
  return U;
}

/* method of many people factorisation */
rci_t _mzd_pls_mmpf(mzd_t *A, mzp_t *P, mzp_t *Q, int k) {
  assert(A->offset == 0);

  rci_t const nrows = A->nrows;
  rci_t const ncols = A->ncols; 
  rci_t curr_row = 0;
  rci_t curr_col = 0;
  rci_t done_row = 0;

  int knar = 0;

  /** compute good k **/

  if(k == 0) {
    k = m4ri_opt_k(nrows, ncols, 0);
    if (k >= 7)
      k = 7;
    if (0.5 * __M4RI_TWOPOW(k) * A->ncols > __M4RI_CPU_L2_CACHE / 2.0)
      k -= 1;
  }

  int kk = __M4RI_PLS_NTABLES * k;

  /** initialise permutations as identity **/

  for(rci_t i = 0; i < ncols; ++i) 
    Q->values[i] = i;

  for(rci_t i = 0; i < nrows; ++i)
    P->values[i] = i;

  mzd_t *T0 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T1 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T2 = mzd_init(__M4RI_TWOPOW(k), ncols);
  mzd_t *T3 = mzd_init(__M4RI_TWOPOW(k), ncols);

  mzd_t *U = mzd_init(kk, ncols);

  /* these are the elimination lookups */

  rci_t *E0 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *E1 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *E2 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *E3 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));

  /* these are the multiplication lookups */

  rci_t *M0 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *M1 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *M2 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));
  rci_t *M3 = (rci_t*)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(rci_t));

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
    wi_t splitblock = (curr_col + kk) / m4ri_radix + 1;

    knar = _mzd_pls_submatrix(A, curr_row, nrows, curr_col, kk, P, Q, pivots, done, &done_row, splitblock);

    /**
     * 2. update A10 
     */

    _mzd_pls_a10(A, P, curr_row, curr_col, splitblock, knar, pivots);

    /**
     * 3. extract U from A0 = (A00 | A10) 
     */

    _mzd_pls_to_u(U, A, curr_row, curr_col, knar, pivots);


    // treat no pivot was found case
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

    if (kk >= 3*k && kk >= 4) {

      int k_[4], knar_[4];
      _kk_setup(kk, knar, k_, knar_, pivots, 4);

      /**
       * 4. generate multiplication and inversion tables T amd E from U 
       */
    
      mzd_make_table_pls(U, 0,                            curr_col,                   k_[0], knar_[0], T0, E0, M0, pivots, 0);
      mzd_make_table_pls(U, 0+knar_[0],                   curr_col+k_[0],             k_[1], knar_[1], T1, E1, M1, pivots+knar_[0], k_[0]);
      mzd_make_table_pls(U, 0+knar_[0]+knar_[1],          curr_col+k_[0]+k_[1],       k_[2], knar_[2], T2, E2, M2, pivots+knar_[0]+knar_[1], k_[0]+k_[1]);
      mzd_make_table_pls(U, 0+knar_[0]+knar_[1]+knar_[2], curr_col+k_[0]+k_[1]+k_[2], k_[3], knar_[3], T3, E3, M3, pivots+knar_[0]+knar_[1]+knar_[2], k_[0]+k_[1]+k_[2]);
      
      /**
       * 5. update A1 = (A01 | A11)
       */
    
      _mzd_pls_a11_4(A, curr_row, done_row+1, curr_col, splitblock,
                     k_[0], knar_[0], T0, M0, 
                     k_[1], knar_[1], T1, M1, 
                     k_[2], knar_[2], T2, M2,
                     k_[3], knar_[3], T3, M3);
      
      /**
       * 6. update A2 = (A02 | A12)
       */

      if (done_row < nrows) {
        mzd_process_rows4_pls(A, done_row + 1, nrows, curr_col, 
                              k_[0], T0, E0, 
                              k_[1], T1, E1, 
                              k_[2], T2, E2,
                              k_[3], T3, E3);
      }
          
    } else if (kk >= 2*k && kk >= 3) {

      int k_[3], knar_[3];
      _kk_setup(kk, knar, k_, knar_, pivots, 3);

      mzd_make_table_pls(U, 0,                   curr_col,             k_[0], knar_[0], T0, E0, M0, pivots, 0);
      mzd_make_table_pls(U, 0+knar_[0],          curr_col+k_[0],       k_[1], knar_[1], T1, E1, M1, pivots+knar_[0], k_[0]);
      mzd_make_table_pls(U, 0+knar_[0]+knar_[1], curr_col+k_[0]+k_[1], k_[2], knar_[2], T2, E2, M2, pivots+knar_[0]+knar_[1], k_[0]+k_[1]);
 
      _mzd_pls_a11_3(A, curr_row, done_row+1, curr_col, splitblock,
                     k_[0], knar_[0], T0, M0, 
                     k_[1], knar_[1], T1, M1, 
                     k_[2], knar_[2], T2, M2);
      
    
      if (done_row < nrows) {
        mzd_process_rows3_pls(A, done_row + 1, nrows, curr_col, 
                              k_[0], T0, E0, 
                              k_[1], T1, E1, 
                              k_[2], T2, E2);
      }

    } else if (kk >= k && kk >= 2) {

      int k_[2], knar_[2];
      _kk_setup(kk, knar, k_, knar_, pivots, 2);

      mzd_make_table_pls(U, 0,          curr_col,       k_[0], knar_[0], T0, E0, M0, pivots, 0);
      mzd_make_table_pls(U, 0+knar_[0], curr_col+k_[0], k_[1], knar_[1], T1, E1, M1, pivots+knar_[0], k_[0]);
      
      _mzd_pls_a11_2(A, curr_row, done_row+1, curr_col, splitblock,
                     k_[0], knar_[0], T0, M0, k_[1], knar_[1], T1, M1);

      if(done_row < nrows) {
        mzd_process_rows2_pls(A, done_row + 1, nrows, curr_col,
                              k_[0], T0, E0,
                              k_[1], T1, E1);
      }
    } else {
      mzd_make_table_pls(U, 0, curr_col, kk, knar, T0, E0, M0, pivots, 0);
      _mzd_pls_a11_1(A, curr_row, done_row+1, curr_col, splitblock, kk, knar, T0, M0);
      
      if(done_row < nrows) {
        mzd_process_rows(A, done_row + 1, nrows, curr_col, kk, T0, E0);
      }
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
  mzd_free(T0);
  mzd_free(T1);
  mzd_free(T2);
  mzd_free(T3);
  m4ri_mm_free(E0);  m4ri_mm_free(M0);
  m4ri_mm_free(E1);  m4ri_mm_free(M1);
  m4ri_mm_free(E2);  m4ri_mm_free(M2);
  m4ri_mm_free(E3);  m4ri_mm_free(M3);
  m4ri_mm_free(done);
  m4ri_mm_free(pivots);

  __M4RI_DD_MZD(A);
  __M4RI_DD_MZP(P);
  __M4RI_DD_MZP(Q);
  __M4RI_DD_RCI(curr_row);

  return curr_row;
}

rci_t _mzd_pluq_mmpf(mzd_t *A, mzp_t *P, mzp_t *Q, int const k) {
  rci_t r = _mzd_pls_mmpf(A, P, Q, k);
  mzd_apply_p_right_trans_tri(A, Q);
  return r;
}
