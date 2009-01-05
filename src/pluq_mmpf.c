/*******************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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

#include <assert.h>

#include "misc.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include "pluq_mmpf.h"
#include "brilliantrussian.h"
#include "grayflex.h"

size_t _mzd_pluq_submatrix(packedmatrix *A, size_t start_row, size_t start_col, int k, permutation *P, permutation *Q, size_t *done)  {
  size_t i, j, l, curr_pos;
  int found;

  for(curr_pos=0; curr_pos < (size_t)k; curr_pos++) {
    found = 0;
    /* search for some pivot */
    for(j = start_col + curr_pos; j < start_col + k; j++) {
      for(i = start_row + curr_pos; i < A->nrows; i++) {
        /* clear before but preserve transformation matrix */
        for(l = 0; l < curr_pos; l++)
          if(done[l] < i) {
            if(mzd_read_bit(A, i, start_col + l))
              mzd_row_add_offset(A, i, start_row + l, start_col + l + 1);
            done[l] = i; /* encode up to which row we added for l already */
          }
        
	if(mzd_read_bit(A, i, j)) {
          found = 1;
          break;
        }
      }
      if(found)
        break;
    }
    if(!found) {
      return curr_pos;
    }

    if (i > P->values[start_row + curr_pos])
      P->values[start_row + curr_pos] = i;
    mzd_row_swap(A, i, start_row + curr_pos);

    if (j > Q->values[start_col + curr_pos])
      Q->values[start_col + curr_pos] = j;
    mzd_col_swap(A, start_col + curr_pos, j);

    done[curr_pos] = i;
  }
  return curr_pos;
}

/* create a table of all 2^k linear combinations */
void mzd_make_table_pluq( packedmatrix *M, size_t r, size_t c, int k, packedmatrix *T, size_t *L) {
  const size_t blockoffset= c/RADIX;
  size_t i, rowneeded;
  size_t twokay= TWOPOW(k);
  size_t wide = T->width - blockoffset;

  word *ti, *ti1, *m;

  ti1 = T->values + blockoffset;
  ti = ti1 + T->width;
#ifdef HAVE_SSE2
  unsigned long incw = 0;
  if (T->width & 1) incw = 1;
  ti += incw;
#endif

  L[0]=0;
  for (i=1; i<twokay; i++) {
    rowneeded = r + codebook[k]->inc[i-1];
    m = M->values + M->rowswap[rowneeded] + blockoffset;

    /* Duff's device loop unrolling */
    register int n = (wide + 7) / 8;
    switch (wide % 8) {
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
    ti+=incw; ti1+=incw;
#endif
    ti += blockoffset;
    ti1 += blockoffset;

    /* U is a basis but not the canonical basis, so we need to read what
       element we just created from T*/
    L[(int)mzd_read_bits(T,i,c,k)] = i;
    
  }
  /* We need fix the table to update the transformation matrix
     correctly; e.g. if the first row has [1 0 1] and we clear a row
     below with [1 0 1] we need to encode that this row is cleared by
     adding the first row only ([1 0 0]).*/
  for(i=1; i < twokay; i++) {
    const word correction = (word)codebook[k]->ord[i];
    mzd_xor_bits(T, i,c, k, correction);
  }
}

void mzd_process_rows2_pluq(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, packedmatrix *T0, size_t *L0, packedmatrix *T1, size_t *L1) {
  size_t r;
  const int ka = k/2;
  const int kb = k-k/2;
  const size_t blocknuma=startcol/RADIX;
  const size_t blocknumb=(startcol+ka)/RADIX;
  const size_t blockoffset = blocknumb - blocknuma;
  size_t wide = M->width - blocknuma;

  if(wide < 4) {
    mzd_process_rows(M, startrow, stoprow, startcol, ka, T0, L0);
    mzd_process_rows(M, startrow, stoprow, startcol + ka, kb, T1, L1);
    return;
  }

  wide -= 2;
  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ (int)mzd_read_bits(M, r, startcol, ka) ];
    word *t0 = T0->values + T0->rowswap[x0] + blocknuma;
    word *m0 = M->values + M->rowswap[r+0] + blocknuma;
    m0[0] ^= t0[0];
    m0[1] ^= t0[1];
/*     m0[2] ^= t0[2]; */
/*     m0[3] ^= t0[3]; */
    const int x1 = L1[ (int)mzd_read_bits(M, r, startcol+ka, kb) ];
    word *t1 = T1->values + T1->rowswap[x1] + blocknumb;
    for(size_t i=blockoffset; i<2; i++) {
      m0[i] ^= t1[i-blockoffset];
    }

    t0+=2;
    t1+=2-blockoffset;
    m0+=2;

    register int n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++;
      case 7:    *m0++ ^= *t0++ ^ *t1++;
      case 6:    *m0++ ^= *t0++ ^ *t1++;
      case 5:    *m0++ ^= *t0++ ^ *t1++;
      case 4:    *m0++ ^= *t0++ ^ *t1++;
      case 3:    *m0++ ^= *t0++ ^ *t1++;
      case 2:    *m0++ ^= *t0++ ^ *t1++;
      case 1:    *m0++ ^= *t0++ ^ *t1++;
      } while (--n > 0);
    }
  }
}

/* extract U from A for table creation */
packedmatrix *_mzd_pluq_to_u(packedmatrix *U, packedmatrix *A, size_t r, size_t c, int k) {
  /* this function call is now rather cheap, but it could be avoided
     completetly if needed */
  assert(U->offset == 0);
  assert(A->offset == 0);
  size_t i, j;
  size_t startcol = (c/RADIX)*RADIX;
  mzd_submatrix(U, A, r, 0, r+k, A->ncols);

  for(i=0; i<(size_t)k; i++)
    for(j=startcol; j<c+i; j++) 
      mzd_write_bit(U, i, j,  0);
  return U;
}

static inline size_t _max_value(size_t *data, size_t length) {
  size_t max = 0;
  for(size_t i=0; i<length; i++) {
    max = MAX(max, data[i]);
  }
  return max;
}

/* method of many people factorisation */
size_t _mzd_pluq_mmpf(packedmatrix *A, permutation * P, permutation * Q, int k) {
  assert(A->offset == 0);
  const size_t nrows = A->nrows; 
  const size_t ncols = A->ncols; 
  size_t curr_pos = 0;
  size_t kbar = 0;
  size_t done_row = 0;

  if(k == 0) {
    k = m4ri_opt_k(nrows, ncols, 0);
    if(k>3)
      k-=2;
  }
  int kk = 2*k;

  for(size_t i = 0; i<ncols; i++)
    Q->values[i] = i;
  for(size_t i = 0; i<nrows; i++)
    P->values[i] = i;

  packedmatrix *T0 = mzd_init(TWOPOW(k), ncols);
  packedmatrix *T1 = mzd_init(TWOPOW(k), ncols);
  packedmatrix *U = mzd_init(kk, ncols);

  size_t *L0 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L1 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *done = (size_t *)m4ri_mm_malloc(kk * sizeof(size_t));

  while(curr_pos < MIN(ncols,nrows)) {
    if(curr_pos + kk > ncols)
      kk = ncols - curr_pos;

    /* 1. compute PLUQ factorisation for a kxk submatrix */
    kbar = _mzd_pluq_submatrix(A, curr_pos, curr_pos, kk, P, Q, done);
    /* 1.5. finish submatrix*/
    done_row = _max_value(done, kbar);
    for(size_t c2=0; c2<kbar && curr_pos + c2 < A->ncols -1; c2++)
      for(size_t r2=done[c2]+1; r2<=done_row; r2++)
        if(mzd_read_bit(A, r2, curr_pos + c2))
          mzd_row_add_offset(A, r2, curr_pos + c2, curr_pos + c2 + 1);

    /* 2. extract U */
    _mzd_pluq_to_u(U, A, curr_pos, curr_pos, kbar);
    
    if(kbar > (size_t)k) {
      const int ka = kbar/2;
      const int kb = kbar - ka;
      /* 2. generate table T */
      mzd_make_table_pluq(U, 0, curr_pos, ka, T0, L0);
      mzd_make_table_pluq(U, 0+ka, curr_pos + ka, kb, T1, L1);
      /* 3. use that table to process remaining rows below */
      mzd_process_rows2_pluq(A, done_row + 1, nrows, curr_pos, kbar, T0, L0, T1, L1);
      curr_pos += kbar;
    } else if(kbar > 0) {
      /* 2. generate table T */
      mzd_make_table_pluq(U, 0, curr_pos, kbar, T0, L0);
      /* 3. use that table to process remaining rows below */
      mzd_process_rows(A, done_row + 1, nrows, curr_pos, kbar, T0, L0);
      curr_pos += kbar;
    } else {
      /* we have at least kk columns which are all zero! */
      int found = 0;
      size_t stop_col = curr_pos + kk;
      while(curr_pos < stop_col) {
        size_t i = curr_pos;
        size_t j = curr_pos;
        found = mzd_find_pivot(A, curr_pos, curr_pos, &i, &j);
        if(found) {
          P->values[curr_pos] = i;
          Q->values[curr_pos] = j;
          mzd_row_swap(A, curr_pos, i);
          mzd_col_swap(A, curr_pos, j);
          for(size_t l = i+1; l<A->nrows; l++)
            if(mzd_read_bit(A, l, curr_pos))
              mzd_row_add_offset(A, l, curr_pos, curr_pos + 1);
          
          curr_pos++;
        } else {
          break;
        }
      }
      if(curr_pos != stop_col)
        break;
    }
  }

  mzd_free(U);
  mzd_free(T0);
  mzd_free(T1);
  m4ri_mm_free(L0);
  m4ri_mm_free(L1);
  m4ri_mm_free(done);
  return curr_pos;
}
