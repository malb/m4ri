 /*******************************************************************
 *
 *            M4RI: Method of the Four Russians Inversion
 *
 *       Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
 *       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhu.ac.uk>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
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
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include <assert.h>

#include "brilliantrussian.h"
#include "grayflex.h"
#include "misc.h"

/**
 * Get k bits starting a position (x,y) from the matrix m.
 *
 * \param m Source matrix.
 * \param x Starting row.
 * \param y Starting column.
 * \param k Number of bits.
 */ 

static inline int _mzd_get_bits(const packedmatrix *m, const int x, const int y, const int k) {
  int truerow = m->rowswap[ x ];
  word temp;

  /* there are two possible situations. Either all bits are in one
   * word or they are spread across two words. */

  if ( (y%RADIX + k -1 ) < RADIX ) {
    /* everything happens in one word here */
    temp =  m->values[ y / RADIX + truerow ]; /* get the value */
    temp <<= y%RADIX; /* clear upper bits */
    temp >>= RADIX - k; /* clear lower bits and move to correct position.*/
    return (int)temp;

  } else { 
    /* two words are affected */
    const int block = y / RADIX + truerow; /* correct block */
    const int spot = (y + k ) % RADIX; /* correct offset */
    /* make room by shifting spot times to the right, and add stuff from the second word */
    temp = (m->values[block] << spot) | ( m->values[block + 1] >> (RADIX - spot) ); 
    return ((int)temp & ((1<<k)-1)); /* clear upper bits and return */
   }
}

/**
 * \brief Perform Gaussian reduction on a submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row endrow (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param endrow Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static int _mzd_gauss_submatrix(packedmatrix *A, int r, int c, int k, int endrow) {
  int i,j,l;
  int start_row = r;
  int found;
  for (j=c; j<c+k; j++) {
    found = 0;
    for (i=start_row; i< endrow; i++) {
      /* first we need to clear the first columns first */
      for (l=0; l<j-c; l++)
        if (mzd_read_bit(A, i, c+l))
          mzd_row_add_offset(A, r+l, i, c+l);
      
      /* pivot? */
      if (mzd_read_bit(A, i, j)) {
        mzd_row_swap(A, i, start_row);
        /* clear above */
        for (l=r; l<start_row; l++) {
          if (mzd_read_bit(A, l, j)) {
            mzd_row_add_offset(A, start_row, l, j);
          }
        }
        start_row++;
        found = 1;
        break;
      }
    }
    if (found==0) {
      return j - c;
    }
  }
  return j - c;
}

void mzd_make_table( packedmatrix *M, int r, int c, int k, packedmatrix *T, int *L) {
  const int homeblock= c/RADIX;
  int i, j, rowneeded, id;
  int twokay= TWOPOW(k);
  unsigned int wide = T->width - homeblock;

  word *ti, *ti1, *m;

  ti1 = T->values + homeblock;
  ti = ti1 + T->width;
#ifdef HAVE_SSE2
  unsigned long incw = 0;
  if (T->width & 1) incw = 1;
  ti += incw;
#endif

  L[0]=0;
  for (i=1; i<twokay; i++) {
    rowneeded = r + codebook[k]->inc[i-1];
    id = codebook[k]->ord[i];
    L[id] = i;

    if (rowneeded >= M->nrows) {
      for (j = wide-1; j >= 0; j--) {
        *ti++ = *ti1++;
      }
#ifdef HAVE_SSE2
      ti+=incw; ti1+=incw;
#endif
    } else {
      m = M->values + M->rowswap[rowneeded] + homeblock;

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
      ti += homeblock;
      ti1 += homeblock;
    }
  }
}

void mzd_process_rows(packedmatrix *M, int startrow, int stoprow, int startcol, int k, packedmatrix *T, int *L) {
  int r;
  const int blocknum=startcol/RADIX;
  int wide = M->width - blocknum;

  for (r=startrow; r+2<=stoprow; r+=2) {
    const int x0 = L[ _mzd_get_bits(M, r+0, startcol, k) ];
    const int x1 = L[ _mzd_get_bits(M, r+1, startcol, k) ];
    
    word *m0 = M->values + M->rowswap[r+0] + blocknum;
    word *t0 = T->values + T->rowswap[x0] + blocknum;

    word *m1 = M->values + M->rowswap[r+1] + blocknum;
    word *t1 = T->values + T->rowswap[x1] + blocknum;

    register int n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 7:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 6:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 5:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 4:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 3:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 2:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      case 1:    *m0++ ^= *t0++; *m1++ ^= *t1++;
      } while (--n > 0);
    }
  }

  for( ; r<stoprow; r++) {
    const int x0 = L[ _mzd_get_bits(M, r, startcol, k) ];
    word *m0 = M->values + M->rowswap[r] + blocknum;
    word *t0 = T->values + T->rowswap[x0] + blocknum;

    register int n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *m0++ ^= *t0++;
      case 7:    *m0++ ^= *t0++;
      case 6:    *m0++ ^= *t0++;
      case 5:    *m0++ ^= *t0++;
      case 4:    *m0++ ^= *t0++;
      case 3:    *m0++ ^= *t0++;
      case 2:    *m0++ ^= *t0++;
      case 1:    *m0++ ^= *t0++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows2(packedmatrix *M, int startrow, int stoprow, int startcol, int k, packedmatrix *T0, int *L0, packedmatrix *T1, int *L1) {
  int r;
  const int blocknum=startcol/RADIX;
  const int wide = M->width - blocknum;

  const int ka = k/2;
  const int kb = k-k/2;

  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ _mzd_get_bits(M, r, startcol, ka)];
    const int x1 = L1[ _mzd_get_bits(M, r, startcol+ka, kb)];
    if(x0 == 0 && x1 == 0)
      continue;
    word * m0 = M->values + M->rowswap[r] + blocknum;
    const word *t0 = T0->values + T0->rowswap[x0] + blocknum;
    const word *t1 = T1->values + T1->rowswap[x1] + blocknum;

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

void mzd_process_rows3(packedmatrix *M, int startrow, int stoprow, int startcol, int k, packedmatrix *T0, int *L0, packedmatrix *T1, int *L1, packedmatrix *T2, int *L2) {
  int r;
  const int blocknum=startcol/RADIX;
  const int wide = M->width - blocknum;

  int rem = k%3;
  
  const int ka = k/3 + ((rem>=2) ? 1 : 0);
  const int kb = k/3 + ((rem>=1) ? 1 : 0);
  const int kc = k/3;

  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ _mzd_get_bits(M, r, startcol, ka)];
    const int x1 = L1[ _mzd_get_bits(M, r, startcol+ka, kb)];
    const int x2 = L2[ _mzd_get_bits(M, r, startcol+ka+kb, kc)];
    if(x0 == 0 && x1 == 0 && x2 == 0) 
      continue;

    word * m0 = M->values + M->rowswap[r] + blocknum;
    const word *t0 = T0->values + T0->rowswap[x0] + blocknum;
    const word *t1 = T1->values + T1->rowswap[x1] + blocknum;
    const word *t2 = T2->values + T2->rowswap[x2] + blocknum;

    register int n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows4(packedmatrix *M, int startrow, int stoprow, int startcol, int k, packedmatrix *T0, int *L0, packedmatrix *T1, int *L1, packedmatrix *T2, int *L2, packedmatrix *T3, int *L3) {
  int r;
  const int blocknum=startcol/RADIX;
  const int wide = M->width - blocknum;

  int rem = k%4;
  
  const int ka = k/4 + ((rem>=3) ? 1 : 0);
  const int kb = k/4 + ((rem>=2) ? 1 : 0);
  const int kc = k/4 + ((rem>=1) ? 1 : 0);
  const int kd = k/4;

  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ _mzd_get_bits(M, r, startcol, ka)];
    const int x1 = L1[ _mzd_get_bits(M, r, startcol+ka, kb)];
    const int x2 = L2[ _mzd_get_bits(M, r, startcol+ka+kb, kc)];
    const int x3 = L3[ _mzd_get_bits(M, r, startcol+ka+kb+kc, kd)];
    if(x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0) 
      continue;

    word * m0 = M->values + M->rowswap[r] + blocknum;
    const word *t0 = T0->values + T0->rowswap[x0] + blocknum;
    const word *t1 = T1->values + T1->rowswap[x1] + blocknum;
    const word *t2 = T2->values + T2->rowswap[x2] + blocknum;
    const word *t3 = T3->values + T3->rowswap[x3] + blocknum;

    register int n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++;
      } while (--n > 0);
    }
  }
}

int mzd_reduce_m4ri(packedmatrix *A, int full, int k, packedmatrix *T, int *L) {
  /**
   * The algorithm works as follows:
   *
   * Step 1.Denote the first column to be processed in a given
   * iteration as \f$a_i\f$. Then, perform Gaussian elimination on the
   * first \f$3k\f$ rows after and including the \f$i\f$-th row to
   * produce an identity matrix in \f$a_{i,i} ... a_{i+k-1,i+k-1},\f$
   * and zeroes in \f$a_{i+k,i} ... a_{i+3k-1,i+k-1}\f$.
   *
   * Step 2. Construct a table consisting of the \f$2^k\f$ binary strings of
   * length k in a Gray code.  Thus with only \f$2^k\f$ vector
   * additions, all possible linear combinations of these k rows
   * have been precomputed.
   *
   *
   * Step 3. One can rapidly process the remaining rows from \f$i +
   * 3k\f$ until row \f$m\f$ (the last row) by using the table. For
   * example, suppose the \f$j\f$-th row has entries \f$a_{j,i}
   * ... a_{j,i+k-1}\f$ in the columns being processed. Selecting the
   * row of the table associated with this k-bit string, and adding it
   * to row j will force the k columns to zero, and adjust the
   * remaining columns from \f$ i + k\f$ to n in the appropriate way,
   * as if Gaussian elimination had been performed.
   *
   * Step 4. While the above form of the algorithm will reduce a
   * system of boolean linear equations to unit upper triangular form,
   * and thus permit a system to be solved with back substitution, the
   * M4RI algorithm can also be used to invert a matrix, or put the
   * system into reduced row echelon form (RREF). Simply run Step 3
   * on rows \f$0 ... i-1\f$ as well as on rows \f$i + 3k
   * ... m\f$. This only affects the complexity slightly, changing the
   * 2.5 coeffcient to 3
   */

  const int ncols = A->ncols; 
  int r = 0;
  int c = 0;
  int kbar = 0;

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, A->ncols, 0);
    if (k>5) {
      k -= 4;
    }
  }
  int kk = 4*k;

  packedmatrix *T0 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T1 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T2 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T3 = mzd_init(TWOPOW(k), A->ncols);
  int *L0 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  int *L1 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  int *L2 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  int *L3 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));

  while(c<ncols) {
    if(c+kk > A->ncols) {
      kk = ncols - c;
    }
    kbar = _mzd_gauss_submatrix(A, r, c, kk, A->nrows);

    if (kbar>3*k) {
      const int rem = kbar%4;
      const int ka = kbar/4 + ((rem>=3) ? 1 : 0);
      const int kb = kbar/4 + ((rem>=2) ? 1 : 0);
      const int kc = kbar/4 + ((rem>=1) ? 1 : 0);
      const int kd = kbar/4;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_process_rows4(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      if(full)
        mzd_process_rows4(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      
    } else if (kbar>2*k) {
      int rem = kbar%3;
      int ka = kbar/3 + ((rem>=2) ? 1 : 0);
      int kb = kbar/3 + ((rem>=1) ? 1 : 0);
      int kc = kbar/3;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_process_rows3(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2);
      if(full)
        mzd_process_rows3(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2);
      
    } else if (kbar>k) {
      const int ka = kbar/2;
      const int kb = kbar - ka;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_process_rows2(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1);
      if(full)
        mzd_process_rows2(A, 0, r, c, kbar, T0, L0, T1, L1);
      
    } else if(kbar > 0) {
      mzd_make_table(A, r, c, kbar, T0, L0);
      mzd_process_rows(A, r+kbar, A->nrows, c, kbar, T0, L0);
      if(full)
        mzd_process_rows(A, 0, r, c, kbar, T0, L0);
    }
    r += kbar;
    c += kbar;
    if(kk!=kbar) {
      c++;
    }
  }

  mzd_free(T0);
  m4ri_mm_free(L0);
  mzd_free(T1);
  m4ri_mm_free(L1);
  mzd_free(T2);
  m4ri_mm_free(L2);
  mzd_free(T3);
  m4ri_mm_free(L3);
  return r;
}

void mzd_top_reduce_m4ri(packedmatrix *A, int k, packedmatrix *T, int *L) {
  const int ncols = A->ncols; 
  int r = 0;
  int c = 0;
  int kbar = 0;

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, A->ncols, 0);
    if (k>5) {
      k -= 4;
    }
  }
  int kk = 4*k;

  packedmatrix *T0 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T1 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T2 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T3 = mzd_init(TWOPOW(k), A->ncols);
  int *L0 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  int *L1 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  int *L2 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));
  int *L3 = (int *)m4ri_mm_calloc(TWOPOW(k), sizeof(int));

  while(c<ncols) {
    if(c+kk > A->ncols) {
      kk = ncols - c;
    }
    kbar = _mzd_gauss_submatrix(A, r, c, kk, A->nrows);

    if (kbar>3*k) {
      const int rem = kbar%4;
      const int ka = kbar/4 + ((rem>=3) ? 1 : 0);
      const int kb = kbar/4 + ((rem>=2) ? 1 : 0);
      const int kc = kbar/4 + ((rem>=1) ? 1 : 0);
      const int kd = kbar/4;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_process_rows4(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3);
      
    } else if (kbar>2*k) {
      int rem = kbar%3;
      int ka = kbar/3 + ((rem>=2) ? 1 : 0);
      int kb = kbar/3 + ((rem>=1) ? 1 : 0);
      int kc = kbar/3;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_process_rows3(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2);
      
    } else if (kbar>k) {
      const int ka = kbar/2;
      const int kb = kbar - ka;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_process_rows2(A, 0, r, c, kbar, T0, L0, T1, L1);
      
    } else if(kbar > 0) {
      mzd_make_table(A, r, c, kbar, T0, L0);
      mzd_process_rows(A, 0, r, c, kbar, T0, L0);
    }
    r += kbar;
    c += kbar;
    if(kk!=kbar) {
      c++;
    }
  }

  mzd_free(T0);
  m4ri_mm_free(L0);
  mzd_free(T1);
  m4ri_mm_free(L1);
  mzd_free(T2);
  m4ri_mm_free(L2);
  mzd_free(T3);
  m4ri_mm_free(L3);
}

packedmatrix *mzd_invert_m4ri(packedmatrix *m, packedmatrix *I, int k) {
  packedmatrix *big = mzd_concat(NULL, m, I);
  int size=m->ncols;
  if (k == 0) {
    k = m4ri_opt_k(m->nrows, m->ncols, 0);
  }
  int twokay=TWOPOW(k);
  int i;
  packedmatrix *T=mzd_init(twokay, size*2);
  int *L=(int *)m4ri_mm_malloc(twokay * sizeof(int));
  packedmatrix *answer;
  
  mzd_reduce_m4ri(big, TRUE, k, T, L);
  
  for(i=0; i < size; i++) {
    if (!mzd_read_bit(big, i,i )) {
      answer = NULL;
      break;
    }
  }
  if (i == size)
    answer=mzd_submatrix(NULL, big, 0, size, size, size*2);
  
  m4ri_mm_free(L);
  mzd_free(T);
  mzd_free(big);
  
  return answer;
}

packedmatrix *mzd_mul_m4rm_t(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  packedmatrix *AT, *BT, *CT;
  
  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm_t: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  
  AT = mzd_transpose(NULL, A);
  BT = mzd_transpose(NULL, B);
  
  CT = mzd_init(B->ncols, A->nrows);
  CT = _mzd_mul_m4rm_impl(CT, BT, AT, k, 0);
  
  mzd_free(AT);
  mzd_free(BT);

  C = mzd_transpose(C, CT);
  mzd_free(CT);
  return C;
}

packedmatrix *mzd_mul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  int a = A->nrows;
  int c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C (%d x %d) has wrong dimensions.\n", C->nrows, C->ncols);
  }
  return _mzd_mul_m4rm_impl(C, A, B, k, TRUE);
}

packedmatrix *mzd_addmul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  int a = A->nrows;
  int c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm A ncols (%d) need to match B nrows (%d) .\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C has wrong dimensions.\n");
  }
  return _mzd_mul_m4rm_impl(C, A, B, k, FALSE);
}

packedmatrix *_mzd_mul_m4rm_impl_old(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, int clear) {
  int i,j, a_nr, a_nc, b_nc;
  int truerow;
  unsigned int x;


  a_nr = A->nrows;
  a_nc = A->ncols;
  b_nc = B->ncols;

  if (b_nc < RADIX-10) {
    return mzd_mul_naiv(C, A, B);
  }

  int wide = C->width;

  /* clear first */
  if (clear) {
    for (i=0; i<C->nrows; i++) {
      truerow = C->rowswap[i];
      for (j=0; j<C->width; j++) {
  	C->values[truerow + j] = 0;
     }
    }
  }

  const int blocksize = MZD_MUL_BLOCKSIZE;

  if (k == 0) {
    k = m4ri_opt_k(blocksize, a_nc, b_nc);
  }

  packedmatrix *T = mzd_init(TWOPOW(k), b_nc);
  int *L = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));

  long s, start;
  const long end = a_nc/k;

    for (start=0; start + blocksize <= a_nr; start += blocksize) {
      for(i=0; i < end; i++) {
        mzd_make_table( B, i*k, 0, k, T, L);
        for(s = 0; s < blocksize; s++) {
          j = start + s;
          x = L[ _mzd_get_bits(A, j, i*k, k) ];
          word * const c = C->values + C->rowswap[j];
          const word * const t = T->values + T->rowswap[x];
          for(int ii=0; ii<wide ; ii++)
            c[ii] ^= t[ii];
        }
      }
    }
  
    for(i=0; i < a_nc/k; i++) {
      mzd_make_table( B, i*k, 0, k, T, L);
      for(s = 0; s < a_nr-start; s++) {
        j = start + s;
        x = L[ _mzd_get_bits(A, j, i*k, k) ];
        /* mzd_combine( C,j,0, C,j,0,  T,x,0); */
        word *c = C->values + C->rowswap[j];
      const word *t = T->values + T->rowswap[x];
      for(int ii=0; ii<wide ; ii++)
        c[ii] ^= t[ii];
      }
    }
    
    /* handle rest */
    if (a_nc%k) {
      mzd_make_table( B, a_nc/k * k , 0, a_nc%k, T, L);
    
    for(j = 0; j<a_nr; j++) {
      x = L[ _mzd_get_bits(A, j, i*k, a_nc%k) ];
      mzd_combine(C,j,0, C,j,0, T,x,0);
    }
    }
    
    mzd_free(T);
    m4ri_mm_free(L);
    return C;
}

#ifdef HAVE_SSE2
static inline void _mzd_combine8_sse2(word *c, word *t1, word *t2, word *t3, word *t4, word *t5, word *t6, word *t7, word *t8, int wide) {
  int i;
  /* assuming t1 ... t8 are aligned, but c might not be */
  if (ALIGNMENT(c,16)==0) {
    __m128i *__c = (__m128i*)c;
    __m128i *__t1 = (__m128i*)t1;
    __m128i *__t2 = (__m128i*)t2;
    __m128i *__t3 = (__m128i*)t3;
    __m128i *__t4 = (__m128i*)t4;
    __m128i *__t5 = (__m128i*)t5;
    __m128i *__t6 = (__m128i*)t6;
    __m128i *__t7 = (__m128i*)t7;
    __m128i *__t8 = (__m128i*)t8;
    const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xF);
    __m128i xmm1;
    
    while(__c < eof) {
      xmm1 = _mm_xor_si128(*__c, *__t1++);
      xmm1 = _mm_xor_si128(xmm1, *__t2++);
      xmm1 = _mm_xor_si128(xmm1, *__t3++);
      xmm1 = _mm_xor_si128(xmm1, *__t4++);
      xmm1 = _mm_xor_si128(xmm1, *__t5++);
      xmm1 = _mm_xor_si128(xmm1, *__t6++);
      xmm1 = _mm_xor_si128(xmm1, *__t7++);
      xmm1 = _mm_xor_si128(xmm1, *__t8++);
      *__c++ = xmm1;
    }
    c  = (word*)__c;
    t1 = (word*)__t1;
    t2 = (word*)__t2;
    t3 = (word*)__t3;
    t4 = (word*)__t4;
    t5 = (word*)__t5;
    t6 = (word*)__t6;
    t7 = (word*)__t7;
    t8 = (word*)__t8;
    wide = ((sizeof(word)*wide)%16)/sizeof(word);
  }
  for(i=0; i<wide; i++) {
    c[i] ^= t1[i] ^ t2[i] ^ t3[i] ^ t4[i] ^ t5[i] ^ t6[i] ^ t7[i] ^ t8[i];
  }
}
#endif

#ifdef HAVE_SSE2
static inline void _mzd_combine4_sse2(word *c, word *t1, word *t2, word *t3, word *t4, int wide) {
  int i;
  /* assuming t1 ... t8 are aligned, but c might not be */
  if (ALIGNMENT(c,16)==0) {
    __m128i *__c = (__m128i*)c;
    __m128i *__t1 = (__m128i*)t1;
    __m128i *__t2 = (__m128i*)t2;
    __m128i *__t3 = (__m128i*)t3;
    __m128i *__t4 = (__m128i*)t4;
    const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xF);
    __m128i xmm1;
    
    while(__c < eof) {
      xmm1 = _mm_xor_si128(*__c, *__t1++);
      xmm1 = _mm_xor_si128(xmm1, *__t2++);
      xmm1 = _mm_xor_si128(xmm1, *__t3++);
      xmm1 = _mm_xor_si128(xmm1, *__t4++);
      *__c++ = xmm1;
    }
    c  = (word*)__c;
    t1 = (word*)__t1;
    t2 = (word*)__t2;
    t3 = (word*)__t3;
    t4 = (word*)__t4;
    wide = ((sizeof(word)*wide)%16)/sizeof(word);
  }
  for(i=0; i<wide; i++) {
    c[i] ^= t1[i] ^ t2[i] ^ t3[i] ^ t4[i];
  }
}
#endif //HAVE_SSE2

#ifdef HAVE_SSE2

#ifdef M4RM_GRAY8
#define _MZD_COMBINE _mzd_combine8_sse2(c, t1, t2, t3, t4, t5, t6, t7, t8, wide)
#else //M4RM_GRAY8
#define _MZD_COMBINE _mzd_combine4_sse2(c, t1, t2, t3, t4, wide)
#endif //M4RM_GRAY8

#else //HAVE_SSE2

#ifdef M4RM_GRAY8
#define _MZD_COMBINE for(ii=0; ii<wide ; ii++) c[ii] ^= t1[ii] ^ t2[ii] ^ t3[ii] ^ t4[ii] ^ t5[ii] ^ t6[ii] ^ t7[ii] ^ t8[ii]
#else //M4RM_GRAY8
#define _MZD_COMBINE for(ii=0; ii<wide ; ii++) c[ii] ^= t1[ii] ^ t2[ii] ^ t3[ii] ^ t4[ii]
#endif //M4RM_GRAY8

#endif //HAVE_SSE2

packedmatrix *_mzd_mul_m4rm_impl(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, int clear) {
  /**
   * The algorithm proceeds as follows:
   * 
   * Step 1. Make a Gray code table of all the \f$2^k\f$ linear combinations
   * of the \f$k\f$ rows of \f$B_i\f$.  Call the \f$x\f$-th row
   * \f$T_x\f$.
   *
   * Step 2. Read the entries 
   *    \f$a_{j,(i-1)k+1}, a_{j,(i-1)k+2} , ... , a_{j,(i-1)k+k}.\f$
   *
   * Let \f$x\f$ be the \f$k\f$ bit binary number formed by the
   * concatenation of \f$a_{j,(i-1)k+1}, ... , a_{j,ik}\f$.
   *
   * Step 3. for \f$h = 1,2, ... , c\f$ do
   *   calculate \f$C_{jh} = C_{jh} + T_{xh}\f$.
   */

  int i,j;
  int ii;
  unsigned int x1, x2, x3, x4;
  word *t1, *t2, *t3, *t4;

#ifdef M4RM_GRAY8
  unsigned int x5, x6, x7, x8;
  word *t5, *t6, *t7, *t8;
#endif

  word *c;

  int a_nr = A->nrows;
  int a_nc = A->ncols;
  int b_nc = B->ncols;

  if (b_nc < RADIX-10) {
    return mzd_mul_naiv(C, A, B);
  }

  int wide = C->width;

  /* clear first */
  int truerow;
  if (clear) {
    for (i=0; i<C->nrows; i++) {
      truerow = C->rowswap[i];
      for (j=0; j<C->width; j++) {
  	C->values[truerow + j] = 0;
     }
    }
  }

  const int blocksize = MZD_MUL_BLOCKSIZE;

  if (k == 0) {
    k = m4ri_opt_k(blocksize, a_nc, b_nc);
#ifdef M4RM_GRAY8
    if (k>3)
      k -= 2;
#else
    if (k>2)
      k -= 1;
#endif
  }
  packedmatrix *T1 = mzd_init(TWOPOW(k), b_nc);
  int *L1 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
  packedmatrix *T2 = mzd_init(TWOPOW(k), b_nc);
  int *L2 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
  packedmatrix *T3 = mzd_init(TWOPOW(k), b_nc);
  int *L3 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
  packedmatrix *T4 = mzd_init(TWOPOW(k), b_nc);
  int *L4 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));

#ifdef M4RM_GRAY8
  packedmatrix *T5 = mzd_init(TWOPOW(k), b_nc);
  int *L5 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
  packedmatrix *T6 = mzd_init(TWOPOW(k), b_nc);
  int *L6 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
  packedmatrix *T7 = mzd_init(TWOPOW(k), b_nc);
  int *L7 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
  packedmatrix *T8 = mzd_init(TWOPOW(k), b_nc);
  int *L8 = (int *)m4ri_mm_malloc(TWOPOW(k) * sizeof(int));
#endif

  /* process stuff that fits into multiple of k first, but blockwise (babystep-giantstep)*/
  long babystep, giantstep;
#ifdef M4RM_GRAY8
  const int kk = 8*k;
#else
  const int kk = 4*k;
#endif
  const long end = a_nc/kk;

  for (giantstep=0; giantstep + blocksize <= a_nr; giantstep += blocksize) {
    for(i=0; i < end; i++) {
      mzd_make_table( B, i*kk, 0, k, T1, L1);
      mzd_make_table( B, i*kk+k, 0, k, T2, L2);
      mzd_make_table( B, i*kk+k+k, 0, k, T3, L3);
      mzd_make_table( B, i*kk+k+k+k, 0, k, T4, L4);
#ifdef M4RM_GRAY8
      mzd_make_table( B, i*kk+k+k+k+k, 0, k, T5, L5);
      mzd_make_table( B, i*kk+k+k+k+k+k, 0, k, T6, L6);
      mzd_make_table( B, i*kk+k+k+k+k+k+k, 0, k, T7, L7);
      mzd_make_table( B, i*kk+k+k+k+k+k+k+k, 0, k, T8, L8);
#endif   
      for(babystep = 0; babystep < blocksize; babystep++) {
        j = giantstep + babystep;
        x1 = L1[ _mzd_get_bits(A, j, i*kk, k) ];
        x2 = L2[ _mzd_get_bits(A, j, i*kk+k, k) ];
        x3 = L3[ _mzd_get_bits(A, j, i*kk+k+k, k) ];
        x4 = L4[ _mzd_get_bits(A, j, i*kk+k+k+k, k) ];
#ifdef M4RM_GRAY8 
        x5 = L5[ _mzd_get_bits(A, j, i*kk+k+k+k+k, k) ];
        x6 = L6[ _mzd_get_bits(A, j, i*kk+k+k+k+k+k, k) ];
        x7 = L7[ _mzd_get_bits(A, j, i*kk+k+k+k+k+k+k, k) ];
        x8 = L8[ _mzd_get_bits(A, j, i*kk+k+k+k+k+k+k+k, k) ];
#endif
        c = C->values + C->rowswap[j];
        t1 = T1->values + T1->rowswap[x1];
        t2 = T2->values + T2->rowswap[x2];
        t3 = T3->values + T3->rowswap[x3];
        t4 = T4->values + T4->rowswap[x4];
#ifdef M4RM_GRAY8
        t5 = T5->values + T5->rowswap[x5];
        t6 = T6->values + T6->rowswap[x6];
        t7 = T7->values + T7->rowswap[x7];
        t8 = T8->values + T8->rowswap[x8];
#endif
        _MZD_COMBINE;
      }
    }
  }
  
  for(i=0; i < end; i++) {
    mzd_make_table( B, i*kk, 0, k, T1, L1);
    mzd_make_table( B, i*kk+k, 0, k, T2, L2);
    mzd_make_table( B, i*kk+k+k, 0, k, T3, L3);
    mzd_make_table( B, i*kk+k+k+k, 0, k, T4, L4);
#ifdef M4RM_GRAY8
    mzd_make_table( B, i*kk+k+k+k+k, 0, k, T5, L5);
    mzd_make_table( B, i*kk+k+k+k+k+k, 0, k, T6, L6);
    mzd_make_table( B, i*kk+k+k+k+k+k+k, 0, k, T7, L7);
    mzd_make_table( B, i*kk+k+k+k+k+k+k+k, 0, k, T8, L8);
#endif
    for(babystep = 0; babystep < a_nr - giantstep; babystep++) {
      j = giantstep + babystep;
      x1 = L1[ _mzd_get_bits(A, j, i*kk, k) ];
      x2 = L2[ _mzd_get_bits(A, j, i*kk+k, k) ];
      x3 = L3[ _mzd_get_bits(A, j, i*kk+k+k, k) ];
      x4 = L4[ _mzd_get_bits(A, j, i*kk+k+k+k, k) ];
#ifdef M4RM_GRAY8
      x5 = L5[ _mzd_get_bits(A, j, i*kk+k+k+k+k, k) ];
      x6 = L6[ _mzd_get_bits(A, j, i*kk+k+k+k+k+k, k) ];
      x7 = L7[ _mzd_get_bits(A, j, i*kk+k+k+k+k+k+k, k) ];
      x8 = L8[ _mzd_get_bits(A, j, i*kk+k+k+k+k+k+k+k, k) ];
#endif
      c = C->values + C->rowswap[j];
      t1 = T1->values + T1->rowswap[x1];
      t2 = T2->values + T2->rowswap[x2];
      t3 = T3->values + T3->rowswap[x3];
      t4 = T4->values + T4->rowswap[x4];
#ifdef M4RM_GRAY8
      t5 = T5->values + T5->rowswap[x5];
      t6 = T6->values + T6->rowswap[x6];
      t7 = T7->values + T7->rowswap[x7];
      t8 = T8->values + T8->rowswap[x8];
#endif
      _MZD_COMBINE;
    }
  }

  /* handle stuff that doesn't fit into multiple of kk */
  if (a_nc%kk) {
    for (i=end*kk/k; i < (a_nc)/k; i++) {
      mzd_make_table( B, i*k, 0, k, T1, L1);
      for(j = 0; j<a_nr; j++) {
        x1 = L1[ _mzd_get_bits(A, j, i*k, k) ];
        c = C->values + C->rowswap[j];
        t1 = T1->values + T1->rowswap[x1];
        for(ii=0; ii<wide; ii++) {
          c[ii] ^= t1[ii];
        }
      }
    }
    /* handle stuff that doesn't fit into multiple of k */
    if (a_nc%k) {
      mzd_make_table( B, a_nc/k * k , 0, a_nc%k, T1, L1);
      for(j = 0; j<a_nr; j++) {
        x1 = L1[ _mzd_get_bits(A, j, i*k, a_nc%k) ];
        c = C->values + C->rowswap[j];
        t1 = T1->values + T1->rowswap[x1];
        for(ii=0; ii<wide; ii++) {
          c[ii] ^= t1[ii];
        }
      }
    }
  }

  mzd_free(T1);
  m4ri_mm_free(L1);
  mzd_free(T2);
  m4ri_mm_free(L2);
  mzd_free(T3);
  m4ri_mm_free(L3);
  mzd_free(T4);
  m4ri_mm_free(L4);
#ifdef M4RM_GRAY8
  mzd_free(T5);
  m4ri_mm_free(L5);
  mzd_free(T6);
  m4ri_mm_free(L6);
  mzd_free(T7);
  m4ri_mm_free(L7);
  mzd_free(T8);
  m4ri_mm_free(L8);
#endif
  return C;
}

