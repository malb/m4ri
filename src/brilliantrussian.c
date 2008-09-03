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

#include "misc.h"

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include <assert.h>

#include "brilliantrussian.h"
#include "grayflex.h"


/**
 * \brief Perform Gaussian reduction to reduced row echelon form on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row endrow (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzd_gauss_submatrix_full(packedmatrix *A, size_t r, size_t c, size_t end_row, int k) {
  size_t i,j,l;
  size_t start_row = r;
  int found;
  for (j=c; j<c+k; j++) {
    found = 0;
    for (i=start_row; i< end_row; i++) {
      /* first we need to clear the first columns */
      for (l=0; l<j-c; l++)
        if (mzd_read_bit(A, i, c+l))
          mzd_row_add_offset(A, i, r+l, c+l);
      
      /* pivot? */
      if (mzd_read_bit(A, i, j)) {
        mzd_row_swap(A, i, start_row);
        /* clear above */
        for (l=r; l<start_row; l++) {
          if (mzd_read_bit(A, l, j)) {
            mzd_row_add_offset(A, l, start_row, j);
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

/**
 * \brief Perform Gaussian reduction to upper triangular matrix on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row end_row (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzd_gauss_submatrix(packedmatrix *A, size_t r, size_t c, size_t end_row, int k) {
  size_t i,j,l;
  size_t start_row = r;
  int found;
  for (j=c; j<c+k; j++) {
    found = 0;
    for (i=start_row; i< end_row; i++) {
      /* first we need to clear the first columns */
      for (l=0; l<j-c; l++)
        if (mzd_read_bit(A, i, c+l))
          mzd_row_add_offset(A, i, r+l, c+l);
      
      /* pivot? */
      if (mzd_read_bit(A, i, j)) {
        mzd_row_swap(A, i, start_row);
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

/**
 * \brief Given a submatrix in upper triangular form compute the
 * reduced row echelon form.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row end_row (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzd_gauss_submatrix_top(packedmatrix *A, size_t r, size_t c, int k) {
  size_t j,l;
  size_t start_row = r;
  for (j=c; j<c+k; j++) {
    for (l=r; l<start_row; l++) {
      if (mzd_read_bit(A, l, j)) {
        mzd_row_add_offset(A, l, start_row, j);
      }
    }
    start_row++;
  }
  return k;
}

static inline void _mzd_copy_back_rows(packedmatrix *A, packedmatrix *U, size_t r, size_t c, size_t k) {
  size_t startblock = c/RADIX;
  size_t width = A->width - startblock;
  size_t i, j;
  for (i=0 ; i < k ; i++) {
    const word * const src = U->values + U->rowswap[i] + startblock;
    word *const dst = A->values + A->rowswap[r+i] + startblock;
    for (j=0; j< width; j++) {
      dst[j] = src[j];
    }
  }
}

void mzd_make_table( packedmatrix *M, size_t r, size_t c, int k, packedmatrix *T, size_t *L) {
  const size_t homeblock= c/RADIX;
  size_t i, j, rowneeded, id;
  size_t twokay= TWOPOW(k);
  size_t wide = T->width - homeblock;

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
      for (j = 0; j < wide; j++) {
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

void mzd_process_rows(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, packedmatrix *T, size_t *L) {
  size_t r;
  const size_t blocknum=startcol/RADIX;
  size_t wide = M->width - blocknum;

  for (r=startrow; r+2<=stoprow; r+=2) {
    const int x0 = L[ (int)mzd_read_bits(M, r+0, startcol, k) ];
    const int x1 = L[ (int)mzd_read_bits(M, r+1, startcol, k) ];
    
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
    const int x0 = L[ (int)mzd_read_bits(M, r, startcol, k) ];
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

void mzd_process_rows2(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, packedmatrix *T0, size_t *L0, packedmatrix *T1, size_t *L1) {
  size_t r;
  const size_t blocknum=startcol/RADIX;
  const size_t wide = M->width - blocknum;

  const int ka = k/2;
  const int kb = k-k/2;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ (int)mzd_read_bits(M, r, startcol, ka)];
    const int x1 = L1[ (int)mzd_read_bits(M, r, startcol+ka, kb)];
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

void mzd_process_rows3(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, packedmatrix *T0, size_t *L0, packedmatrix *T1, size_t *L1, packedmatrix *T2, size_t *L2) {
  size_t r;
  const size_t blocknum=startcol/RADIX;
  const size_t wide = M->width - blocknum;

  int rem = k%3;
  
  const int ka = k/3 + ((rem>=2) ? 1 : 0);
  const int kb = k/3 + ((rem>=1) ? 1 : 0);
  const int kc = k/3;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ (int)mzd_read_bits(M, r, startcol, ka)];
    const int x1 = L1[ (int)mzd_read_bits(M, r, startcol+ka, kb)];
    const int x2 = L2[ (int)mzd_read_bits(M, r, startcol+ka+kb, kc)];
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

void mzd_process_rows4(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, 
                       packedmatrix *T0, size_t *L0, packedmatrix *T1, size_t *L1, packedmatrix *T2, size_t *L2, packedmatrix *T3, size_t *L3) {
  size_t r;
  const size_t blocknum=startcol/RADIX;
  const size_t wide = M->width - blocknum;

  int rem = k%4;
  
  const int ka = k/4 + ((rem>=3) ? 1 : 0);
  const int kb = k/4 + ((rem>=2) ? 1 : 0);
  const int kc = k/4 + ((rem>=1) ? 1 : 0);
  const int kd = k/4;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ (int)mzd_read_bits(M, r, startcol, ka)];
    const int x1 = L1[ (int)mzd_read_bits(M, r, startcol+ka, kb)];
    const int x2 = L2[ (int)mzd_read_bits(M, r, startcol+ka+kb, kc)];
    const int x3 = L3[ (int)mzd_read_bits(M, r, startcol+ka+kb+kc, kd)];
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

void mzd_process_rows5(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, 
                       packedmatrix *T0, size_t *L0, packedmatrix *T1, size_t *L1, packedmatrix *T2, size_t *L2, packedmatrix *T3, size_t *L3,
                       packedmatrix *T4, size_t *L4) {
  size_t r;
  const size_t blocknum=startcol/RADIX;
  const size_t wide = M->width - blocknum;
  int rem = k%5;
  
  const int ka = k/5 + ((rem>=4) ? 1 : 0);
  const int kb = k/5 + ((rem>=3) ? 1 : 0);
  const int kc = k/5 + ((rem>=2) ? 1 : 0);
  const int kd = k/5 + ((rem>=1) ? 1 : 0);
  const int ke = k/5;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r=startrow; r<stoprow; r++) {
    
    const int x0 = L0[ (int)mzd_read_bits(M, r, startcol, ka)];
    const int x1 = L1[ (int)mzd_read_bits(M, r, startcol+ka, kb)];
    const int x2 = L2[ (int)mzd_read_bits(M, r, startcol+ka+kb, kc)];
    const int x3 = L3[ (int)mzd_read_bits(M, r, startcol+ka+kb+kc, kd)];
    const int x4 = L4[ (int)mzd_read_bits(M, r, startcol+ka+kb+kc+kd, ke)];

    if(x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0) 
      continue;

    word * m0 = M->values + M->rowswap[r] + blocknum;
    const word *t0 = T0->values + T0->rowswap[x0] + blocknum;
    const word *t1 = T1->values + T1->rowswap[x1] + blocknum;
    const word *t2 = T2->values + T2->rowswap[x2] + blocknum;
    const word *t3 = T3->values + T3->rowswap[x3] + blocknum;
    const word *t4 = T4->values + T4->rowswap[x4] + blocknum;
    
    register int n = (wide + 7) / 8;
    switch (wide % 8) {
    case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++;
      } while (--n > 0);
    }
  }
}

void mzd_process_rows6(packedmatrix *M, size_t startrow, size_t stoprow, size_t startcol, int k, 
                       packedmatrix *T0, size_t *L0, packedmatrix *T1, size_t *L1, packedmatrix *T2, size_t *L2, packedmatrix *T3, size_t *L3,
                       packedmatrix *T4, size_t *L4, packedmatrix *T5, size_t *L5) {
  size_t r;
  const size_t blocknum=startcol/RADIX;
  const size_t wide = M->width - blocknum;

  int rem = k%6;
  
  const int ka = k/6 + ((rem>=5) ? 1 : 0);
  const int kb = k/6 + ((rem>=4) ? 1 : 0);
  const int kc = k/6 + ((rem>=3) ? 1 : 0);
  const int kd = k/6 + ((rem>=2) ? 1 : 0);
  const int ke = k/6 + ((rem>=1) ? 1 : 0);;
  const int kf = k/6;

#ifdef HAVE_OPENMP
#pragma omp parallel for private(r) shared(startrow, stoprow) schedule(dynamic,32) if(stoprow-startrow > 128)
#endif
  for(r=startrow; r<stoprow; r++) {
    const int x0 = L0[ (int)mzd_read_bits(M, r, startcol, ka)];
    const int x1 = L1[ (int)mzd_read_bits(M, r, startcol+ka, kb)];
    const int x2 = L2[ (int)mzd_read_bits(M, r, startcol+ka+kb, kc)];
    const int x3 = L3[ (int)mzd_read_bits(M, r, startcol+ka+kb+kc, kd)];
    const int x4 = L4[ (int)mzd_read_bits(M, r, startcol+ka+kb+kc+kd, ke)];
    const int x5 = L5[ (int)mzd_read_bits(M, r, startcol+ka+kb+kc+kd+ke, kf)];
    
    if(x0 == 0 && x1 == 0 && x2 == 0 && x3 == 0 && x4 == 0 && x5 == 0) 
      continue;

    word * m0 = M->values + M->rowswap[r] + blocknum;
    const word *t0 = T0->values + T0->rowswap[x0] + blocknum;
    const word *t1 = T1->values + T1->rowswap[x1] + blocknum;
    const word *t2 = T2->values + T2->rowswap[x2] + blocknum;
    const word *t3 = T3->values + T3->rowswap[x3] + blocknum;
    const word *t4 = T4->values + T4->rowswap[x4] + blocknum;
    const word *t5 = T5->values + T5->rowswap[x5] + blocknum;

    register int n = (wide + 7) / 8;
    switch (wide % 8) {
      case 0: do { *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 7:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 6:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 5:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 4:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 3:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 2:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      case 1:    *m0++ ^= *t0++ ^ *t1++ ^ *t2++ ^ *t3++ ^ *t4++ ^ *t5++;
      } while (--n > 0);
    }
  }
}

int mzd_reduce_m4ri(packedmatrix *A, int full, int k, packedmatrix *T, size_t *L) {
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

  const size_t ncols = A->ncols; 
  size_t r = 0;
  size_t c = 0;
  int kbar = 0;

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, A->ncols, 0);
    if (k>=7)
      k = 7;
    if ( (6*(1<<k)*A->ncols / 8.0) > CPU_L2_CACHE / 2.0 )
      k -= 1;
  }
  /*printf("k: %d\n",k);*/
  int kk = 6*k;

  packedmatrix *U  = mzd_init(kk, A->ncols);
  packedmatrix *T0 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T1 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T2 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T3 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T4 = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *T5 = mzd_init(TWOPOW(k), A->ncols);
  size_t *L0 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L1 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L2 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L3 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L4 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L5 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));

  while(c<ncols) {
    if(c+kk > A->ncols) {
      kk = ncols - c;
    }
    if (full) {
      kbar = _mzd_gauss_submatrix_full(A, r, c, A->nrows, kk);
    } else {
      kbar = _mzd_gauss_submatrix(A, r, c, A->nrows, kk);
      U = mzd_submatrix(U, A, r, 0, r+kbar, A->ncols);
      _mzd_gauss_submatrix_top(A, r, c, kbar);
    }

    if (kbar>5*k) {
      const int rem = kbar%6;
      const int ka = kbar/6 + ((rem>=5) ? 1 : 0);
      const int kb = kbar/6 + ((rem>=4) ? 1 : 0);
      const int kc = kbar/6 + ((rem>=3) ? 1 : 0);
      const int kd = kbar/6 + ((rem>=2) ? 1 : 0);
      const int ke = kbar/6 + ((rem>=1) ? 1 : 0);;
      const int kf = kbar/6;

      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
      mzd_make_table(A, r+ka+kb+kc+kd+ke, c, kf, T5, L5);
      mzd_process_rows6(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);
      if(full)
        mzd_process_rows6(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);

  } else if (kbar>4*k) {
      const int rem = kbar%5;
      const int ka = kbar/5 + ((rem>=4) ? 1 : 0);
      const int kb = kbar/5 + ((rem>=3) ? 1 : 0);
      const int kc = kbar/5 + ((rem>=2) ? 1 : 0);
      const int kd = kbar/5 + ((rem>=1) ? 1 : 0);
      const int ke = kbar/5;
      mzd_make_table(A, r, c, ka, T0, L0);
      mzd_make_table(A, r+ka, c, kb, T1, L1);
      mzd_make_table(A, r+ka+kb, c, kc, T2, L2);
      mzd_make_table(A, r+ka+kb+kc, c, kd, T3, L3);
      mzd_make_table(A, r+ka+kb+kc+kd, c, ke, T4, L4);
      mzd_process_rows5(A, r+kbar, A->nrows, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      if(full)
        mzd_process_rows5(A, 0, r, c, kbar, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      
    } else if (kbar>3*k) {
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

    if (!full) {
      _mzd_copy_back_rows(A, U, r, c, kbar);
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
  mzd_free(T4);
  m4ri_mm_free(L4);
  mzd_free(T5);
  m4ri_mm_free(L5);
  mzd_free(U);
  return r;
}

void mzd_top_reduce_m4ri(packedmatrix *A, int k, packedmatrix *T, size_t *L) {
  const size_t ncols = A->ncols; 
  size_t r = 0;
  size_t c = 0;
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
  size_t *L0 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L1 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L2 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L3 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));

  while(c<ncols) {
    if(c+kk > A->ncols) {
      kk = ncols - c;
    }
    kbar = _mzd_gauss_submatrix_full(A, r, c, A->nrows, kk);

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
  size_t size=m->ncols;
  if (k == 0) {
    k = m4ri_opt_k(m->nrows, m->ncols, 0);
  }
  size_t twokay=TWOPOW(k);
  size_t i;
  packedmatrix *T=mzd_init(twokay, size*2);
  size_t *L=(size_t *)m4ri_mm_malloc(twokay * sizeof(size_t));
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
  CT = _mzd_mul_m4rm(CT, BT, AT, k, 0);
  
  mzd_free(AT);
  mzd_free(BT);

  C = mzd_transpose(C, CT);
  mzd_free(CT);
  return C;
}

packedmatrix *mzd_mul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  size_t a = A->nrows;
  size_t c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm: A ncols (%d) need to match B nrows (%d).\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C (%d x %d) has wrong dimensions.\n", C->nrows, C->ncols);
  }
  return _mzd_mul_m4rm(C, A, B, k, TRUE);
}

packedmatrix *mzd_addmul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k) {
  size_t a = A->nrows;
  size_t c = B->ncols;

  if(A->ncols != B->nrows) 
    m4ri_die("mzd_mul_m4rm A ncols (%d) need to match B nrows (%d) .\n", A->ncols, B->nrows);
  if (C == NULL) {
    C = mzd_init(a, c);
  } else {
    if (C->nrows != a || C->ncols != c)
      m4ri_die("mzd_mul_m4rm: C has wrong dimensions.\n");
  }
  return _mzd_mul_m4rm(C, A, B, k, FALSE);
}

#ifdef HAVE_SSE2
static inline void _mzd_combine8(word *c, word *t1, word *t2, word *t3, word *t4, word *t5, word *t6, word *t7, word *t8, int wide) {
  size_t i;
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
#else

#define _mzd_combine8(c,t1,t2,t3,t4,t5,t6,t7,t8,wide) for(ii=0; ii<wide ; ii++) c[ii] ^= t1[ii] ^ t2[ii] ^ t3[ii] ^ t4[ii] ^ t5[ii] ^ t6[ii] ^ t7[ii] ^ t8[ii]

#endif

#ifdef HAVE_SSE2
static inline void _mzd_combine4(word *c, word *t1, word *t2, word *t3, word *t4, size_t wide) {
  size_t i;
  /* assuming t1 ... t4 are aligned, but c might not be */
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
#else

#define _mzd_combine4(c, t1, t2, t3, t4, wide) for(ii=0; ii<wide ; ii++) c[ii] ^= t1[ii] ^ t2[ii] ^ t3[ii] ^ t4[ii]

#endif //HAVE_SSE2

#ifdef HAVE_SSE2
static inline void _mzd_combine2(word *c, word *t1, word *t2, size_t wide) {
  size_t i;
  /* assuming t1 ... t2 are aligned, but c might not be */
  if (ALIGNMENT(c,16)==0) {
    __m128i *__c = (__m128i*)c;
    __m128i *__t1 = (__m128i*)t1;
    __m128i *__t2 = (__m128i*)t2;
    const __m128i *eof = (__m128i*)((unsigned long)(c + wide) & ~0xF);
    __m128i xmm1;
    
    while(__c < eof) {
      xmm1 = _mm_xor_si128(*__c, *__t1++);
      xmm1 = _mm_xor_si128(xmm1, *__t2++);
      *__c++ = xmm1;
    }
    c  = (word*)__c;
    t1 = (word*)__t1;
    t2 = (word*)__t2;
    wide = ((sizeof(word)*wide)%16)/sizeof(word);
  }
  for(i=0; i<wide; i++) {
    c[i] ^= t1[i] ^ t2[i];
  }
}
#else

#define _mzd_combine2(c, t1, t2, wide) for(ii=0; ii<wide ; ii++) c[ii] ^= t1[ii] ^ t2[ii]

#endif //HAVE_SSE2


#ifdef M4RM_GRAY8
#define _MZD_COMBINE _mzd_combine8(c, t1, t2, t3, t4, t5, t6, t7, t8, wide)
#else //M4RM_GRAY8
#define _MZD_COMBINE _mzd_combine4(c, t1, t2, t3, t4, wide)
#endif //M4RM_GRAY8

packedmatrix *_mzd_mul_m4rm(packedmatrix *C, packedmatrix *A, packedmatrix *B, int k, int clear) {
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

  size_t i,j;
  size_t ii;
  unsigned int x1, x2, x3, x4;
  word *t1, *t2, *t3, *t4;

#ifdef M4RM_GRAY8
  unsigned int x5, x6, x7, x8;
  word *t5, *t6, *t7, *t8;
#endif

  word *c;

  size_t a_nr = A->nrows;
  size_t a_nc = A->ncols;
  size_t b_nc = B->ncols;

  if (b_nc < RADIX-10) {
    if(clear)
      return mzd_mul_naiv(C, A, B);
    else
      return mzd_addmul_naiv(C, A, B);
  }

  size_t wide = C->width;

  /* clear first */
  size_t truerow;
  if (clear) {
    for (i=0; i<C->nrows; i++) {
      truerow = C->rowswap[i];
      for (j=0; j<C->width-1; j++) {
  	C->values[truerow + j] = 0;
      }
      C->values[truerow + j] &= ~LEFT_BITMASK(C->ncols);
    }
  }

  const size_t blocksize = MZD_MUL_BLOCKSIZE;

  if (k == 0) {
    k = m4ri_opt_k(blocksize, a_nc, b_nc);
#ifdef M4RM_GRAY8
    if (k>3)
      k -= 2;
    /* reduce k further if that has a chance of hitting L1 */
    const size_t tsize = (int)(0.8*(TWOPOW(k) * b_nc));
    if(tsize > CPU_L1_CACHE && tsize/2 <= CPU_L1_CACHE)
      k -= 1;
#else
    if (k>2)
      k -= 1;
#endif
  }

#ifndef M4RM_GRAY8
  size_t *buffer = (size_t*)m4ri_mm_malloc(4 * TWOPOW(k) * sizeof(size_t));
#else
  size_t *buffer = (size_t*)m4ri_mm_malloc(8 * TWOPOW(k) * sizeof(size_t));
#endif

  packedmatrix *T1 = mzd_init(TWOPOW(k), b_nc);
  size_t *L1 = buffer;
  packedmatrix *T2 = mzd_init(TWOPOW(k), b_nc);
  size_t *L2 = buffer + 1*TWOPOW(k);
  packedmatrix *T3 = mzd_init(TWOPOW(k), b_nc);
  size_t *L3 = buffer + 2*TWOPOW(k);
  packedmatrix *T4 = mzd_init(TWOPOW(k), b_nc);
  size_t *L4 = buffer + 3*TWOPOW(k);

#ifdef M4RM_GRAY8
  packedmatrix *T5 = mzd_init(TWOPOW(k), b_nc);
  size_t *L5 = buffer + 4*TWOPOW(k);
  packedmatrix *T6 = mzd_init(TWOPOW(k), b_nc);
  size_t *L6 = buffer + 5*TWOPOW(k);
  packedmatrix *T7 = mzd_init(TWOPOW(k), b_nc);
  size_t *L7 = buffer + 6*TWOPOW(k);
  packedmatrix *T8 = mzd_init(TWOPOW(k), b_nc);
  size_t *L8 = buffer + 7*TWOPOW(k);
#endif

  /* process stuff that fits into multiple of k first, but blockwise (babystep-giantstep)*/
  size_t babystep, giantstep;
#ifdef M4RM_GRAY8
  const int kk = 8*k;
#else
  const int kk = 4*k;
#endif
  const size_t end = a_nc/kk;

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
        x1 = L1[ (int)mzd_read_bits(A, j, i*kk, k) ];
        x2 = L2[ (int)mzd_read_bits(A, j, i*kk+k, k) ];
        x3 = L3[ (int)mzd_read_bits(A, j, i*kk+k+k, k) ];
        x4 = L4[ (int)mzd_read_bits(A, j, i*kk+k+k+k, k) ];
#ifdef M4RM_GRAY8
        x5 = L5[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k, k) ];
        x6 = L6[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k+k, k) ];
        x7 = L7[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k+k+k, k) ];
        x8 = L8[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k+k+k+k, k) ];
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
      x1 = L1[ (int)mzd_read_bits(A, j, i*kk, k) ];
      x2 = L2[ (int)mzd_read_bits(A, j, i*kk+k, k) ];
      x3 = L3[ (int)mzd_read_bits(A, j, i*kk+k+k, k) ];
      x4 = L4[ (int)mzd_read_bits(A, j, i*kk+k+k+k, k) ];
#ifdef M4RM_GRAY8
      x5 = L5[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k, k) ];
      x6 = L6[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k+k, k) ];
      x7 = L7[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k+k+k, k) ];
      x8 = L8[ (int)mzd_read_bits(A, j, i*kk+k+k+k+k+k+k+k, k) ];
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
        x1 = L1[ (int)mzd_read_bits(A, j, i*k, k) ];
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
        x1 = L1[ (int)mzd_read_bits(A, j, i*k, a_nc%k) ];
        c = C->values + C->rowswap[j];
        t1 = T1->values + T1->rowswap[x1];
        for(ii=0; ii<wide; ii++) {
          c[ii] ^= t1[ii];
        }
      }
    }
  }

  mzd_free(T1);
  mzd_free(T2);
  mzd_free(T3);
  mzd_free(T4);
#ifdef M4RM_GRAY8
  mzd_free(T5);
  mzd_free(T6);
  mzd_free(T7);
  mzd_free(T8);
#endif
  m4ri_mm_free(buffer);
  return C;
}

/*
 * Experimental scratch code, do not call.
 */

void _mzd_lqup_submatrix_finish(packedmatrix *A, size_t start_col, int k) {
  size_t c,r2,r;
  for(r=0 ; r < (size_t)k; r++) {

    // clear up to submatrix up to start_col
    if (start_col >= RADIX)
      for(c=0; c<start_col/RADIX-1; c++)
        A->values[A->rowswap[r] + c] = 0;
    mzd_clear_bits(A, r, RADIX*(start_col/RADIX), start_col%RADIX);
    
    // clear L
    for (c=0; c<r; c++)
      mzd_write_bit(A, r, start_col + c, 0);

    // clear U
    for(r2=0; r2<r ; r2++) {
      if(mzd_read_bit(A, r2, start_col + r))
        mzd_row_add_offset(A, r2, r, start_col + r);
    }
    // clear the pivot bit
    mzd_write_bit(A, r, start_col+r, 0);
  }
}

size_t _mzd_lqup_submatrix(packedmatrix *A, size_t r, size_t c, size_t end_row, int k, permutation *P, permutation *Q)  {
  size_t i,j,l;
  size_t start_row = r;
  int found;
  for (j=c; j<c+k; j++) {
    found = 0;
    for (i=start_row; i< end_row; i++) {
      if (mzd_read_bit(A, i, j)) {
        P->values[start_row] = i;
        mzd_row_swap_offset(A, i, start_row, j);
        /* clear below but preserve transformation matrix */
        for(l=start_row+1; l<end_row; l++) {
          if (mzd_read_bit(A, l, j))
            mzd_row_add_offset(A, l, start_row, j+1);
        }
        start_row++;
        found = 1;
        break;
      }
    }
    if(!found) {
      return j-c;
    }
  }
  return j - c;
}

size_t _mzd_lqup_m4rf(packedmatrix *A, int k, permutation * P, permutation * Q) {
  const size_t ncols = A->ncols; 
  size_t r = 0;
  size_t c = 0;
  int kbar = 0;

  if (k == 0) {
    k = m4ri_opt_k(A->nrows, A->ncols, 0);
  }

  if (Q == NULL)
    Q = mzp_init(A->ncols);

  packedmatrix *T = mzd_init(TWOPOW(k), A->ncols);
  packedmatrix *I = mzd_init(k, A->ncols);
  size_t *L = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));

  while(c<ncols) {
    if(c+k > A->ncols)
      k = ncols - c;

    /* 1. compute LQUP factorisation for a kxk submatrix */
    kbar = _mzd_lqup_submatrix(A, r, c, MIN(A->nrows,r+k), k, P, Q);
    printf("kbar: %d c: %d\n",kbar, (int)c);

    if(kbar > 0) {
      /* 2. compute RREF for LQUP submatrix to generate the table T */
      mzd_set_ui(I, 0);
      I = mzd_submatrix(I, A, r, 0, r+kbar, A->ncols);
      _mzd_lqup_submatrix_finish(I, c, kbar);
      mzd_print_matrix(I);

      /* 3. generate table T */
      mzd_make_table(I, 0, c, kbar, T, L);

      /* 4. use that table to process remaining rows below */
      mzd_process_rows(A, r+kbar, A->nrows, c, kbar, T, L);
    }
    
    r += kbar;
    c += kbar;
    if(kbar==0) {
      // we would need to do something about Q[i]
      c++;
    }
    printf("A\n");
    mzd_print_matrix(A);
  }

  mzd_free(I);
  mzd_free(T);
  m4ri_mm_free(L);
  return r;
}
